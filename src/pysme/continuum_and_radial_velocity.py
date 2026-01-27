# -*- coding: utf-8 -*-
"""
Determine continuum based on continuum mask
and fit best radial velocity to observation
"""

import logging
import warnings

import emcee
import numpy as np
from scipy.constants import speed_of_light
from scipy.interpolate import splev, splrep
from scipy.ndimage import binary_dilation
from scipy.optimize import least_squares
from scipy.signal import correlate
from tqdm import tqdm

from pysme.util import disable_progress_bars

from .iliffe_vector import Iliffe_vector
from .sme import MASK_VALUES
from .sme_synth import SME_DLL
from .util import show_progress_bars

logger = logging.getLogger(__name__)

c_light = speed_of_light * 1e-3  # speed of light in km/s


class ContinuumNormalizationAbstract:
    def __init__(self):
        self.cscale_type = "none"
        pass

    def __call__(self, sme, x_syn, y_syn, segments, rvel=0, least_squares_kwargs=None):
        if least_squares_kwargs is None:
            least_squares_kwargs = {}
        raise NotImplementedError

    def apply(self, wave, smod, cwave, cscale, segments):
        return apply_continuum(wave, smod, cwave, cscale, self.cscale_type, segments)


class ContinuumNormalizationMask(ContinuumNormalizationAbstract):
    def __init__(self):
        super().__init__()
        self.cscale_type = "mask"

    def __call__(
        self,
        sme,
        x_syn,
        y_syn,
        segments,
        rvel=0,
        least_squares_kwargs=None,
    ):
        """
        Fit a polynomial to the spectrum points marked as continuum
        The degree of the polynomial fit is determined by sme.cscale_flag

        Parameters
        ----------
        sme : SME_Struct
            input sme structure with sme.sob, sme.wave, and sme.mask
        segment : int
            index of the wavelength segment to use, or -1 when dealing with the whole spectrum

        Returns
        -------
        cscale : array of size (ndeg + 1,)
            polynomial coefficients of the continuum fit, in numpy order, i.e. largest exponent first
        """

        if segments < 0:
            return sme.cscale

        if least_squares_kwargs is None:
            least_squares_kwargs = {}

        if "spec" not in sme or "wave" not in sme:
            # If there is no observation, we have no continuum scale
            warnings.warn("Missing data for continuum fit")
            cscale = [1]
        elif sme.cscale_flag in ["none", -3]:
            cscale = [1]
        elif sme.cscale_flag in ["fix", -1, -2]:
            # Continuum flag is set to no continuum
            cscale = sme.cscale[segments]
        else:
            # fit a line to the continuum points
            ndeg = sme.cscale_degree

            # Extract points in this segment
            x, y = sme.wave, sme.spec
            if "mask" in sme:
                m = sme.mask
            else:
                m = sme.spec.copy()
                m[:] = sme.mask_value["line"]

            if "uncs" in sme:
                u = sme.uncs
            else:
                u = sme.spec.copy()
                u[:] = 1
            x, y, m, u = x[segments], y[segments], m[segments], u[segments]

            # Set continuum mask
            if np.all((m & MASK_VALUES.CONT) == 0):
                # If no continuum mask has been set
                # Use the effective wavelength ranges of the lines to determine continuum points
                logger.info(
                    "No Continuum mask was set in segment %s, "
                    "Using effective wavelength range of lines to find continuum instead",
                    segments,
                )
                cont = self.get_continuum_mask(x, y, sme.linelist, mask=m)
                # Save mask for next iteration
                m[cont == MASK_VALUES.CONT] = MASK_VALUES.CONT
                logger.debug("Continuum mask points: %i", np.count_nonzero(cont == 2))

            cont = (m & MASK_VALUES.CONT) != 0

            x = x - x[0]
            x, y, u = x[cont], y[cont], u[cont]

            # Fit polynomial
            try:
                func = lambda coef: (np.polyval(coef, x) - y) / u
                c0 = np.polyfit(x, y, deg=ndeg)
                res = least_squares(func, x0=c0, **least_squares_kwargs)
                cscale = res.x
            except TypeError:
                warnings.warn("Could not fit continuum, set continuum mask?")
                cscale = [1]

        return cscale

    def get_continuum_mask(self, wave, synth, linelist, threshold=0.1, mask=None):
        """
        Use the effective wavelength range of the lines,
        to find wavelength points that should be unaffected by lines
        However one usually has to ignore the weak lines, as most points are affected by one line or another
        Therefore keep increasing the threshold until enough lines have been found (>10%)

        Parameters
        ----------
        wave : array of size (n,)
            wavelength points
        linelist : LineList
            LineList object that was input into the Radiative Transfer
        threshold : float, optional
            starting threshold, lines with depth below this value are ignored
            the actual threshold is increased until enough points are found (default: 0.1)

        Returns
        -------
        mask : array(bool) of size (n,)
            True for points between lines and False for points within lines
        """

        if "depth" not in linelist.columns:
            raise ValueError(
                "No depth specified in the linelist, can't auto compute the mask"
            )

        if threshold <= 0:
            threshold = 0.01

        if mask is None:
            mask = np.full(len(wave), MASK_VALUES.LINE)

        # TODO make this better
        dll = SME_DLL()
        dll.linelist = linelist

        width = dll.GetLineRange()
        # TODO: optimize this
        temp = False
        while np.count_nonzero(temp) < len(wave) * 0.1:
            temp = np.full(len(wave), True)
            for i, line in enumerate(width):
                if linelist["depth"][i] > threshold:
                    w = (wave >= line[0]) & (wave <= line[1])
                    temp[w] = False

            # TODO: Good value to increase threshold by?
            temp[mask == MASK_VALUES.BAD] = False
            threshold *= 1.1

        mask[temp] = MASK_VALUES.CONT

        logger.debug("Ignoring lines with depth < %f", threshold)
        return mask


class ContinuumNormalizationMCMC(ContinuumNormalizationAbstract):
    def __init__(self):
        super().__init__()
        self.cscale_tyoe = "mcmc"

    def __call__(
        self,
        sme,
        x_syn,
        y_syn,
        segments,
        rvel=0,
        least_squares_kwargs=None,
    ):
        """
        Fits both radial velocity and continuum level simultaneously
        by comparing the synthetic spectrum to the observation

        The best fit is determined using robust least squares between
        a shifted and scaled synthetic spectrum and the observation

        Parameters
        ----------
        sme : SME_Struct
            contains the observation
        segment : int
            wavelength segment to fit
        x_syn : array of size (ngrid,)
            wavelength of the synthetic spectrum
        y_syn : array of size (ngrid,)
            intensity of the synthetic spectrum

        Returns
        -------
        vrad : float
            radial velocity in km/s
        vrad_unc : float
            radial velocity uncertainty in km/s
        cscale : array of size (ndeg+1,)
            polynomial coefficients of the continuum
        cscale_unc : array if size (ndeg + 1,)
            uncertainties of the continuum coefficients
        """

        if np.isscalar(segments):
            segments = [segments]
        nseg = len(segments)

        if least_squares_kwargs is None:
            least_squares_kwargs = {}

        if sme.cscale_flag in ["none", "fix"] and sme.vrad_flag in [
            "none",
            "fix",
        ]:
            vrad, vunc, cscale, cunc = null_result(nseg, sme.cscale_degree)
            if sme.vrad_flag == "fix":
                vrad = sme.vrad[segments]
            if sme.cscale_flag == "fix":
                cscale = sme.cscale[segments]
            return vrad, vunc, cscale, cunc

        if "spec" not in sme or "wave" not in sme:
            # No observation no radial velocity
            logger.warning("Missing data for radial velocity/continuum determination")
            return null_result(nseg, sme.cscale_degree)

        if "mask" not in sme:
            sme.mask = np.full(sme.spec.size, MASK_VALUES.LINE)
        if "uncs" not in sme:
            sme.uncs = np.full(sme.spec.size, 1.0)

        if np.all(sme.mask_bad[segments].ravel()):
            warnings.warn(
                "Only bad pixels in this segments, can't determine radial velocity/continuum",
                UserWarning,
            )
            return null_result(nseg, sme.cscale_degree)

        if x_syn.ndim == 1:
            x_syn = x_syn[None, :]
        if y_syn.ndim == 1:
            y_syn = y_syn[None, :]

        if x_syn.shape[0] != nseg or y_syn.shape[0] != nseg:
            raise ValueError(
                "Size of synthetic spectrum, does not match the number of requested segments"
            )

        mask = sme.mask_good[segments]
        x_obs = Iliffe_vector(
            [w[m] for w, m in zip(sme.wave[segments].segments, mask.segments)]
        )
        y_obs = sme.spec[segments][mask]
        x_num = x_obs - sme.wave[segments][:, 0]

        if x_obs.size <= sme.cscale_degree:
            warnings.warn(
                "Not enough good pixels to determine radial velocity/continuum"
            )
            return null_result(nseg)

        if sme.cscale_flag in [-3, "none"]:
            cflag = False
            cscale = np.ones((nseg, 1))
            ndeg = 0
        elif sme.cscale_flag in [-1, -2, "fix"]:
            cflag = False
            cscale = sme.cscale[segments]
            ndeg = cscale.shape[1] - 1
        elif sme.cscale_flag in [0, "constant"]:
            ndeg = 0
            cflag = True
        elif sme.cscale_flag in [1, "linear"]:
            ndeg = 1
            cflag = True
        elif sme.cscale_flag in [2, "quadratic"]:
            ndeg = 2
            cflag = True
        else:
            raise ValueError("cscale_flag not recognized")

        if cflag:
            if sme.cscale is not None:
                cscale = sme.cscale[segments]
            else:
                cscale = np.zeros(nseg, ndeg + 1)
                for i, seg in enumerate(segments):
                    cscale[i, -1] = np.nanpercentile(y_obs[seg], 95)

        # Even when the vrad_flag is set to whole
        # you still want to fit the rv of the individual segments
        # just for the continuum fit
        if sme.vrad_flag == "none":
            vrad = np.zeros(len(segments))
            vflag = False
        elif sme.vrad_flag == "whole":
            vrad = sme.vrad[:1]
            vflag = True
        elif sme.vrad_flag == "each":
            vrad = sme.vrad[segments]
            vflag = True
        elif sme.vrad_flag == "fix":
            vrad = sme.vrad[segments]
            vflag = False
        else:
            raise ValueError(f"Radial velocity Flag not understood {sme.vrad_flag}")

        # Limit shift to half an order
        x1, x2 = (
            x_obs[:, 0],
            x_obs[:, [s // 4 for s in x_obs.shape[1]]].ravel(),
        )
        rv_limit = np.abs(c_light * (1 - x2 / x1))
        if sme.vrad_flag == "whole":
            rv_limit = np.min(rv_limit)

        # Use Cross corellatiom as a first guess for the radial velocity
        # This uses the median as the continuum
        if vrad is None:
            y_tmp = np.interp(x_obs, x_syn, y_syn, left=1, right=1)
            corr = correlate(y_obs - np.median(y_obs), y_tmp - 1, mode="same")
            offset = np.argmax(corr)
            x1, x2 = x_obs[offset], x_obs[len(x_obs) // 2]
            vrad = c_light * (1 - x2 / x1)
            if np.abs(vrad) >= rv_limit:
                logger.warning(
                    "Radial Velocity could not be estimated from cross correlation, using initial guess of 0 km/h. Please check results!"
                )
                vrad = 0

        def log_prior(rv, cscale, nwalkers):
            prior = np.zeros(nwalkers)
            # Add the prior here
            # TODO reject too large/small rv values in a prior
            where = np.full(nwalkers, False)
            if vflag:
                where |= np.any(np.abs(rv) > rv_limit, axis=1)
            if cflag:
                where |= np.any(cscale[:, :, -1] < 0, axis=1)
                if ndeg == 1:
                    where |= np.any(
                        (cscale[:, :, -1] + cscale[:, :, -2] * x_num[:, -1]) < 0,
                        axis=1,
                    )
                elif ndeg == 2:
                    for i in range(nseg):
                        where |= np.any(
                            cscale[:, i, None, -1]
                            + cscale[:, i, None, -2] * x_num[i]
                            + cscale[:, i, None, -3] * x_num[i] ** 2
                            < 0,
                            axis=1,
                        )
            prior[where] = -np.inf
            return prior

        def log_prob(par, sep, nseg, ndeg):
            """
            par : array of shape (nwalkers, ndim)
                ndim = 1 for radial velocity + continuum polynomial coeficients
            """
            nwalkers = par.shape[0]
            rv = par[:, :sep] if vflag else vrad[None, :]
            if rv.shape[0] == 1 and nwalkers > 1:
                rv = np.tile(rv, [nwalkers, 1])
            if rv.shape[1] == 1 and nseg > 1:
                rv = np.tile(rv, [1, nseg])
            if cflag:
                cs = par[:, sep:]
                cs.shape = nwalkers, nseg, ndeg + 1
            else:
                cs = cscale[None, ...]

            prior = log_prior(rv, cs, nwalkers)

            # Apply RV shift
            rv_factor = np.sqrt((1 - rv / c_light) / (1 + rv / c_light))
            total = np.zeros(nwalkers)
            for i in range(nseg):
                x = x_obs[i][None, :] * rv_factor[:, i, None]
                model = np.interp(x, x_syn[i], y_syn[i], left=0, right=0)

                # Apply continuum
                y = np.zeros_like(x_num[i])[None, :]
                for j in range(ndeg + 1):
                    y = y * x_num[i] + cs[:, i, j, None]
                model *= y

                # Ignore the non-overlapping parts of the spectrum
                mask = model == 0
                npoints = mask.shape[1] - np.count_nonzero(mask, axis=1)
                resid = (model - y_obs[i]) ** 2
                resid[mask] = 0
                prob = -0.5 * np.sum(resid, axis=-1)
                # Need to rescale here, to account for the ignored points before
                with np.errstate(divide='ignore', invalid='ignore'):
                    prob *= mask.shape[1] / npoints
                prob[np.isnan(prob)] = -np.inf
                total += prob
            return prior + total

        sep = len(vrad) if vflag else 0
        ndim, p0, scale = 0, [], []
        if vflag:
            ndim += len(vrad)
            p0 += list(vrad)
            scale += [1] * len(vrad)
        if cflag:
            ndim += cscale.size
            p0 += list(cscale.ravel())
            scale += [0.001] * cscale.size
        p0 = np.array(p0)[None, :]
        scale = np.array(scale)[None, :]

        max_n = 10000
        ncheck = 100
        nburn = 300
        nwalkers = max(2 * ndim + 1, 10)
        p0 = p0 + np.random.randn(nwalkers, ndim) * scale
        # If the original guess is good then DEMove is much faster, and sometimes just as good
        # However StretchMove is much more robust to the initial starting value
        moves = [
            (emcee.moves.DEMove(), 0.8),
            (emcee.moves.DESnookerMove(), 0.2),
        ]

        sampler = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            log_prob,
            vectorize=True,
            moves=moves,
            args=(sep, nseg, ndeg),
        )
        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty((max_n // ncheck + 1, ndim))
        # This will be useful to testing convergence
        # old_tau = 0

        # Now we'll sample for up to max_n steps
        with tqdm(
            leave=False, desc="RV", total=max_n, disable=not show_progress_bars
        ) as t:
            for _ in sampler.sample(p0, iterations=max_n):
                t.update()
                # Only check convergence every 100 steps
                if sampler.iteration < 2 * nburn or sampler.iteration % ncheck != 0:
                    continue

                # Compute the autocorrelation time so far
                # Using tol=0 means that we'll always get an estimate even
                # if it isn't trustworthy
                tau = sampler.get_autocorr_time(
                    tol=0, discard=sampler.iteration - ncheck
                )
                autocorr[index] = tau
                index += 1

                # Check convergence
                converged = np.all(tau * 100 < sampler.iteration - nburn)
                # converged &= np.all(np.abs(old_tau - tau) < 0.01 * tau)
                # old_tau = tau
                if converged:
                    break

        if sampler.iteration == max_n:
            logger.warning(
                "The radial velocity did not converge within the limit. Check the results!"
            )

        samples = sampler.get_chain(flat=True, discard=nburn)
        _, vrad_unc, _, cscale_unc = null_result(nseg, ndeg)
        if vflag:
            vmin, vrad, vmax = np.percentile(samples[:, :sep], (32, 50, 68), axis=0)
            vrad_unc[:, 0] = vrad - vmin
            vrad_unc[:, 1] = vmax - vrad

        if cflag:
            vmin, cscale, vmax = np.percentile(samples[:, sep:], (32, 50, 68), axis=0)
            vmin.shape = cscale.shape = vmax.shape = nseg, ndeg + 1

            cscale_unc[..., 0] = cscale - vmin
            cscale_unc[..., 1] = vmax - cscale

        if sme.vrad_flag == "whole":
            vrad = np.tile(vrad, [nseg])

        # Squeeze single-segment results to match caller expectations
        if nseg == 1:
            vrad = vrad[0]
            vrad_unc = vrad_unc[0]
            cscale = cscale[0]
            cscale_unc = cscale_unc[0]

        return vrad, vrad_unc, cscale, cscale_unc


class ContinuumNormalizationMatch(ContinuumNormalizationAbstract):
    def __init__(self):
        super().__init__()
        self.cscale_type = "match"
        self.mask = False
        # We over exaggerate the weights on the top of the spectrum
        # this works well to determine the continuum
        # assuming that there is something there
        self.top_factor = 1000
        self.bottom_factor = 1

    def __call__(
        self,
        sme,
        x_syn,
        y_syn,
        segments,
        rvel=0,
        least_squares_kwargs=None,
    ):
        """
        Fit a continuum when no continuum points exist

        Parameters
        ----------
        sme : SME_Struct
            sme structure with observation data
        segment : int
            index of the wavelength segment to fit
        x_syn : array of size (n,)
            wavelengths of the synthetic spectrum
        y_syn : array of size (n,)
            intensity of the synthetic spectrum
        rvel : float, optional
            radial velocity in km/s to apply to the wavelength (default: 0)

        Returns
        -------
        continuum : array of size (ndeg,)
            continuum fit polynomial coefficients
        """

        if sme.cscale_flag == "none":
            return [1]
        elif sme.cscale_flag == "fix":
            return sme.cscale[segments]

        if least_squares_kwargs is None:
            least_squares_kwargs = {}

        # else
        x = sme.wave[segments]
        y = sme.spec[segments]
        u = sme.uncs[segments]

        if self.mask:
            m = sme.mask_cont[segments]
        else:
            m = sme.mask_good[segments]

        # Interpolate the synthetic spectrum to the observed wavelength grid
        yp = apply_radial_velocity([x], [x_syn], [y_syn], [rvel], [0], copy=True)[0]

        # Apply the mask to all vectors
        xs = x - x[0]
        xs, y, u, yp = xs[m], y[m], u[m], yp[m]

        if sme.telluric is not None:
            # by definition sme.telluric is on the same wavelength grid as sme.spec
            # also it should be in the correct restframe
            tell = sme.telluric[segments][m]
            yp *= tell

        # TODO: what should this be?
        if self.top_factor != 0:
            from scipy.ndimage import median_filter

            diff = median_filter(y, 5) - y
            std = 1.5 * np.median(np.abs(diff[diff != 0] - np.median(diff[diff != 0])))
            snr = 1 / std
            factor = 1000 + (snr - 50) * 90
            factor = np.clip(factor, 100, 100_000)
            mod = np.nanmedian(u) / np.nanmedian(y) * (np.nanpercentile(y, 95) - y) ** 2
            u = u + factor * mod

        deg = sme.cscale_degree
        p0 = sme.cscale[segments]

        def func(p):
            model = yp * np.polyval(p, xs)
            resid = (model - y) / u
            # resid[resid > 0] *= self.top_factor
            resid[resid < 0] *= self.bottom_factor
            return resid

        try:
            res = least_squares(
                func,
                x0=p0,
                **least_squares_kwargs,
            )
            popt = res.x
        except (RuntimeError, ValueError) as ex:
            # logger.warning("Failed to determine the continuum: " + ex.msg)
            logger.warning("Failed to determine the continuum: %s", str(ex))
            popt = p0

        return popt[None, :]


class ContinuumNormalizationMatchLines(ContinuumNormalizationMatch):
    def __init__(self):
        super().__init__()
        self.cscale_type = "matchlines"
        self.top_factor = 0
        self.bottom_factor = 10


class ContinuumNormalizationMatchMask(ContinuumNormalizationMatch):
    def __init__(self):
        super().__init__()
        self.cscale_type = "match+mask"
        self.mask = True
        self.top_factor = 0


class ContinuumNormalizationMatchLinesMask(ContinuumNormalizationMatchLines):
    def __init__(self):
        super().__init__()
        self.cscale_type = "matchlines+mask"
        self.mask = True


class ContinuumNormalizationSpline(ContinuumNormalizationAbstract):
    def __init__(self):
        super().__init__()
        self.cscale_type = "spline"
        self.mask = False

    def __call__(
        self,
        sme,
        x_syn,
        y_syn,
        segments,
        rvel,
        least_squares_kwargs=None,
    ):
        if sme.cscale_flag in ["none"]:
            return np.ones(len(sme.spec[segments]))
        elif sme.cscale_flag in ["fix"]:
            if sme.cscale is not None:
                return sme.cscale[segments]
            else:
                return np.ones(len(sme.spec[segments]))

        if least_squares_kwargs is None:
            least_squares_kwargs = {}
        least_squares_kwargs.setdefault("f_scale", 0.01)

        w = sme.wave[segments]
        s = sme.spec[segments]
        u = sme.uncs[segments]

        # Apply RV correction to the synthetic spectrum
        y = apply_radial_velocity([w], [x_syn], [y_syn], [rvel], [0], copy=True)[0]
        # and don't forget the telluric spectrum if available
        if sme.telluric is not None:
            tell = sme.telluric[segments]
            y *= tell

        # Apply the bpm to all arrays
        if self.mask:
            m = sme.mask_cont[segments]
        else:
            m = sme.mask_good[segments]

        wm, sm, ym, um = w[m], s[m], y[m], u[m]

        # Fit the spline like in the polynomial
        # so that synth * spline = obs
        func = lambda p: (ym * splev(wm, (t, p, 3)) - sm) / um
        # We use splrep to find the intial guess for the number of knots and their
        # positions
        # TODO: how do we know how many points to use?
        if sme.cscale_flag == "constant":
            tlen = int(np.round(w.max() - w.min()))
        elif sme.cscale_flag == "linear":
            tlen = 3
        elif sme.cscale_flag == "quadratic":
            tlen = 4
        elif sme.cscale_flag == "cubic":
            tlen = 5
        elif sme.cscale_flag == "quintic":
            tlen = 6
        elif sme.cscale_flag == "quantic":
            tlen = 7
        else:
            tlen = int(sme.cscale_flag)
        # We need to avoid the 0 points
        m2 = ym != 0
        # Get a first guess by dividing by ym
        t = np.linspace(wm[m2].min(), wm[m2].max(), tlen)[1:-1]
        t, c, k = splrep(wm[m2], sm[m2] / ym[m2], w=1 / um[m2], k=3, t=t)
        # Then get a real fit using the function
        res = least_squares(
            func,
            x0=c,
            **least_squares_kwargs,
        )
        # And finally evaluate the continuum
        c = res.x
        coef = splev(w, (t, c, k))
        return coef


class ContinuumNormalizationSplineMask(ContinuumNormalizationSpline):
    def __init__(self):
        super().__init__()
        self.cscale_type = "spline+mask"
        self.mask = True


def apply_radial_velocity(wave, wmod, smod, vrad, segments, copy=False):
    if copy:
        wmod = [np.copy(w) for w in wmod]
        smod = [np.copy(s) for s in smod]

    if vrad is None:
        return smod
    for il in segments:
        if vrad[il] is not None:
            rv_factor = np.sqrt((1 + vrad[il] / c_light) / (1 - vrad[il] / c_light))
            wmod[il] *= rv_factor
        smod[il] = np.interp(wave[il], wmod[il], smod[il])
    return smod


def apply_continuum(wave, smod, cwave, cscale, cscale_type, segments, copy=False):
    '''
    Apply continuum to the spectra according to cscale.
    '''
    if copy:
        smod = np.copy(smod)

    if cscale is None:
        return smod
    for il in segments:
        if cscale[il] is not None and not np.all(cscale[il] == 0):
            if cscale_type in ["spline", "spline+mask"]:
                if len(cscale[il]) != len(smod[il]):
                    cs = np.interp(wave[il], cwave[il], cscale[il])
                else:
                    cs = cscale[il]
                smod[il] *= cs
            else:
                x = wave[il] - wave[il][0]
                smod[il] *= np.polyval(cscale[il], x)
    return smod


def apply_radial_velocity_and_continuum(
    wave, wmod, smod, vrad, cscale, cscale_type, segments, copy=False
):
    """
    Apply the radial velocity and continuum corrections
    to a syntheic spectrum to match the observation

    Parameters
    ----------
    wave : array
        final wavelength array of the observation
    wmod : array
        wavelength array of the synthethic spectrum
    smod : array
        flux array of the synthetic spectrum
    vrad : array, None
        radial velocities in km/s for each segment
    cscale : array, None
        continnum scales for each segment, exact meaning depends on cscale_type
    cscale_type : str
        defines the continuum correction behaviour
    segments : array
        the segments to apply the correction to

    Returns
    -------
    smod : array
        the corrected synthetic spectrum
    """
    smod = apply_radial_velocity(wave, wmod, smod, vrad, segments, copy=copy)
    # The radial velocity shift also interpolates onto the wavelength grid
    smod = apply_continuum(wave, smod, wave, cscale, cscale_type, segments, copy=copy)
    return smod


def null_result(nseg, ndeg=0, ctype=None):
    vrad, vrad_unc = np.zeros(nseg), np.zeros((nseg, 2))
    if ctype in ["spline", "spline+mask"]:
        cscale = [np.ones(ndeg[i]) for i in range(nseg)]
        cscale = Iliffe_vector(values=cscale)
        cscale_unc = [np.zeros(ndeg[i]) for i in range(nseg)]
        cscale_unc = Iliffe_vector(values=cscale_unc)
    else:
        cscale, cscale_unc = np.zeros((nseg, ndeg + 1)), np.zeros((nseg, ndeg + 1, 2))
        cscale[:, -1] = 1
    return vrad, vrad_unc, cscale, cscale_unc


def cross_correlate_segment(x_obs, y_obs, x_syn, y_syn, mask, mask_wider, rv_bounds):
    # Interpolate synthetic observation onto the observation wavelength grid
    y_tmp = np.interp(x_obs, x_syn, y_syn)

    # Normalize both spectra
    y_obs_tmp = np.copy(y_obs)
    y_obs_tmp /= np.nanpercentile(y_obs_tmp, 95)
    y_obs_tmp -= 1
    y_obs_tmp[~mask] = 0

    y_tmp_tmp = np.copy(y_tmp)
    y_tmp_tmp /= np.nanpercentile(y_tmp_tmp, 95)
    y_tmp_tmp -= 1
    mask_wider = np.interp(x_obs, x_syn, mask_wider) > 0.5
    y_tmp_tmp[~mask_wider] = 0

    # Perform cross correaltion between normalized spectra
    corr = correlate(y_obs_tmp, y_tmp_tmp, mode="same")

    # Determine the radial velocity offset
    # and only retain the area within the bounds
    x_mid = x_obs[len(x_obs) // 2]
    x_shift = c_light * (1 - x_mid / x_obs)

    idx = (x_shift >= rv_bounds[0]) & (x_shift <= rv_bounds[1])
    x_shift = x_shift[idx]
    corr = corr[idx]

    return x_shift, corr


def determine_radial_velocity(
    sme,
    x_syn,
    y_syn,
    segment,
    cscale=None,
    whole=False,
    least_squares_kwargs=None,
):
    """
    Calculate radial velocity by using cross correlation and
    least-squares between observation and synthetic spectrum

    Parameters
    ----------
    sme : SME_Struct
        sme structure with observed spectrum and flags
    segment : int
        which wavelength segment to handle, -1 if its using the whole spectrum
    cscale : array of size (ndeg,)
        continuum coefficients, as determined by e.g. determine_continuum
    x_syn : array of size (n,)
        wavelength of the synthetic spectrum
    y_syn : array of size (n,)
        intensity of the synthetic spectrum

    Raises
    ------
    ValueError
        if sme.vrad_flag is not recognized

    Returns
    -------
    rvel : float
        best fit radial velocity for this segment/whole spectrum
        or None if no observation is present
    """

    if least_squares_kwargs is None:
        least_squares_kwargs = {}
    least_squares_kwargs.setdefault("bounds", (-100, 100))
    least_squares_kwargs.setdefault("jac", "3-point")
    least_squares_kwargs.setdefault("loss", "soft_l1")
    least_squares_kwargs.setdefault("method", "dogbox")
    least_squares_kwargs.setdefault("x_scale", "jac")
    least_squares_kwargs.setdefault("ftol", 1e-8)
    least_squares_kwargs.setdefault("xtol", 1e-8)
    least_squares_kwargs.setdefault("gtol", 1e-8)

    if "spec" not in sme or "wave" not in sme:
        # No observation no radial velocity
        warnings.warn("Missing data for radial velocity determination")
        rvel = 0
    elif sme.vrad_flag == "none":
        # vrad_flag says don't determine radial velocity
        rvel = sme.vrad[segment]
    elif sme.vrad_flag == "whole" and not whole:
        # We are inside a segment, but only want to determine rv at the end
        rvel = 0
    elif sme.vrad_flag == "fix":
        rvel = sme.vrad[segment]
    else:
        # Fit radial velocity
        # Extract data
        x, y = sme.wave, sme.spec
        if "mask" in sme:
            m = sme.mask
        else:
            m = Iliffe_vector(
                [np.full(s, MASK_VALUES.LINE, dtype="i4") for s in sme.spec.shape[1]]
            )

        if "uncs" in sme:
            u = sme.uncs
        else:
            u = sme.spec.copy()
            u[:] = 1

        # Only this one segment
        x_obs = x[segment]
        y_obs = y[segment].copy()
        u_obs = u[segment]
        mask = m[segment]
        if sme.telluric is not None:
            tell = sme.telluric[segment]
        else:
            tell = 1

        if sme.vrad_flag == "each":
            # apply continuum
            y_syn = apply_continuum(
                {segment: x_syn},
                {segment: y_syn},
                sme.wave,
                {segment: cscale},
                sme.cscale_type,
                [segment],
            )[segment]
        elif sme.vrad_flag == "whole":
            # All segments
            y_syn = apply_continuum(
                x_syn,
                y_syn,
                sme.wave,
                cscale,
                sme.cscale_type,
                range(len(y_obs)),
            )
        else:
            raise ValueError(
                f"Radial velocity flag {sme.vrad_flag} not recognised, expected one of 'each', 'whole', 'none'"
            )

        # Filter the mask for the correct sections
        # Need to use temporary m instead of mask, so that we can check
        m = (mask & MASK_VALUES.VRAD) != 0
        if not np.any(m):
            # No VRAD specified, use LINE instead
            m = (mask & MASK_VALUES.LINE) != 0
        mask = m
        mask &= u_obs != 0

        # Widen the mask by roughly the amount expected from the rv_bounds
        bounds = least_squares_kwargs["bounds"]
        rv = max(bounds)
        rv_factor = np.sqrt((1 + rv / c_light) / (1 - rv / c_light))
        # mask_wider = mask.copy()
        if sme.vrad_flag == "each":
            iterations = int(
                np.ceil(np.nanmedian((x_obs * rv_factor - x_obs) / np.gradient(x_obs)))
            )
            iterations = max(1, iterations)
            mask_wider = binary_dilation(mask, iterations=iterations)
            mask_wider = np.interp(x_syn, x_obs, mask_wider) > 0.5
        elif sme.vrad_flag == "whole":
            mask_wider = [None] * len(x_obs)
            for i in range(len(x_obs)):
                iterations = int(
                    np.ceil(
                        np.nanmedian(
                            (x_obs[i] * rv_factor - x_obs[i]) / np.gradient(x_obs[i])
                        )
                    )
                )
                iterations = max(1, iterations)
                mask_wider[i] = binary_dilation(mask[i], iterations=iterations)
                mask_wider[i] = np.interp(x_syn[i], x_obs[i], mask_wider[i]) > 0.5
            mask_wider = Iliffe_vector(mask_wider)
        else:
            raise ValueError

        # Get a first rough estimate from cross correlation
        if sme.vrad_flag == "each":
            x_shift, corr = cross_correlate_segment(
                x_obs, y_obs, x_syn, y_syn, mask, mask_wider, bounds
            )
        else:
            # If using several segments we run the cross correlation for each
            # segment seperately and then sum the results
            nseg = len(segment)
            shift = [None for _ in range(nseg)]
            corr = [None for _ in range(nseg)]
            for i in range(len(x_obs)):
                shift[i], corr[i] = cross_correlate_segment(
                    x_obs[i],
                    y_obs[i],
                    x_syn[i],
                    y_syn[i],
                    mask[i],
                    mask_wider[i],
                    bounds,
                )

            n_min = min(len(s) for s in shift)
            x_shift = np.linspace(bounds[0], bounds[1], n_min)
            corrs_interp = [np.interp(x_shift, s, c) for s, c in zip(shift, corr)]
            corrs_interp = np.array(corrs_interp)
            corr = np.sum(corrs_interp, axis=0)

            # Concatenate the segments for the last step
            mask = np.concatenate(mask)
            mask_wider = np.concatenate(mask_wider)
            x_obs = np.concatenate(x_obs)
            y_obs = np.concatenate(y_obs)
            u_obs = np.concatenate(u_obs)
            x_syn = np.concatenate(x_syn)
            y_syn = np.concatenate(y_syn)
            if sme.telluric is not None:
                tell = np.concatenate(tell)

        # Retrieve the initial cross correlation guess
        offset = np.argmax(corr)
        rvel = x_shift[offset]

        # Normalize the spectra to 1
        y_obs = y_obs / np.nanpercentile(y_obs, 95)
        y_obs -= 1
        y_syn = y_syn / np.nanpercentile(y_syn, 95)
        y_syn -= 1

        # Apply mask
        y_obs[~mask] = 0
        y_syn[~mask_wider] = 0
        u_obs = u_obs.copy()
        u_obs[~mask] = 1

        # Then minimize the least squares for a better fit
        # as cross correlation can only find
        def func(rv):
            rv_factor = np.sqrt((1 - rv / c_light) / (1 + rv / c_light))
            shifted = interpolator(x_obs * rv_factor)
            resid = (y_obs - shifted * tell) / u_obs
            resid[~mask] = 0
            resid = np.nan_to_num(resid, copy=False, nan=0, posinf=0, neginf=0)
            return resid

        interpolator = lambda x: np.interp(x, x_syn, y_syn)
        try:
            res = least_squares(func, x0=rvel, **least_squares_kwargs)
            rvel = res.x[0]
        except ValueError:
            logger.warning(f"Could not determine radial velocity for segment {segment}")
    return rvel


def match_rv_continuum(sme, segments, x_syn, y_syn):
    """
    Match both the continuum and the radial velocity of observed/synthetic spectrum

    Note that the parameterization of the continuum is different to old SME !!!

    Parameters
    ----------
    sme : SME_Struct
        input sme structure with all the parameters
    segment : int
        index of the wavelength segment to match, or -1 when dealing with the whole spectrum
    x_syn : array of size (n,)
        wavelength of the synthetic spectrum
    y_syn : array of size (n,)
        intensitz of the synthetic spectrum

    Returns
    -------
    rvel : float
        new radial velocity
    cscale : array of size (ndeg + 1,)
        new continuum coefficients
    """

    cont_func = {
        "mask": ContinuumNormalizationMask,
        "match": ContinuumNormalizationMatch,
        "match+mask": ContinuumNormalizationMatchMask,
        "matchlines": ContinuumNormalizationMatchLines,
        "matchlines+mask": ContinuumNormalizationMatchLinesMask,
        "spline": ContinuumNormalizationSpline,
        "spline+mask": ContinuumNormalizationSplineMask,
        "mcmc": ContinuumNormalizationMCMC,
    }

    vrad, vrad_unc, cscale, cscale_unc = null_result(
        sme.nseg, sme.cscale_degree, sme.cscale_type
    )
    if sme.cscale_flag == "none" and sme.vrad_flag == "none":
        return cscale, cscale_unc, vrad, vrad_unc

    if np.isscalar(segments):
        segments = [segments]

    if not callable(sme.vrad_flag):
        radial_velocity = determine_radial_velocity
    else:
        radial_velocity = sme.vrad_flag

    if not callable(sme.cscale_type):
        continuum_normalization = cont_func[sme.cscale_type]()
    else:
        continuum_normalization = sme.cscale_type

    least_squares_kwargs = dict(
        bounds=sme.vrad_bounds,
        loss=sme.vrad_loss,
        method=sme.vrad_method,
        jac=sme.vrad_jac,
        x_scale=sme.vrad_xscale,
        ftol=sme.vrad_ftol,
        xtol=sme.vrad_xtol,
        gtol=sme.vrad_gtol,
    )

    if sme.vrad_flag == "none":
        pass
    elif sme.vrad_flag == "fix":
        vrad[segments] = sme.vrad[segments]
    elif sme.vrad_flag == "each":
        for s in segments:
            # We only use the continuum mask for the continuum fit,
            # we need the lines for the radial velocity
            vrad[s] = radial_velocity(
                sme,
                x_syn[s],
                y_syn[s],
                s,
                cscale[s],
                whole=False,
                least_squares_kwargs=least_squares_kwargs,
            )
    elif sme.vrad_flag == "whole":
        s = segments
        vrad[s] = radial_velocity(
            sme,
            [x_syn[s] for s in s],
            [y_syn[s] for s in s],
            s,
            cscale[s],
            whole=True,
            least_squares_kwargs=least_squares_kwargs,
        )
    else:
        raise ValueError

    if sme.cscale_flag == "none":
        pass
    elif sme.cscale_flag == "fix":
        if sme.cscale is not None:
            cscale[segments] = sme.cscale[segments]
        else:
            pass
    else:
        if sme.cscale_type == "mcmc":
            for s in segments:
                (
                    vrad[s],
                    vrad_unc[s],
                    cscale[s],
                    cscale_unc[s],
                ) = continuum_normalization(
                    sme,
                    x_syn[s],
                    y_syn[s],
                    s,
                    rvel=vrad[s],
                )
        else:
            for s in segments:
                cscale[s] = continuum_normalization(
                    sme,
                    x_syn[s],
                    y_syn[s],
                    s,
                    rvel=vrad[s],
                    least_squares_kwargs=dict(
                        bounds=sme.cscale_bounds,
                        loss=sme.cscale_loss,
                        method=sme.cscale_method,
                        jac=sme.cscale_jac,
                        x_scale=sme.cscale_xscale,
                        ftol=sme.cscale_ftol,
                        xtol=sme.cscale_xtol,
                        gtol=sme.cscale_gtol,
                    ),
                )

    # Keep values from unused segments
    select = np.arange(sme.nseg)
    mask = np.full(select.shape, True)
    for seg in segments:
        mask &= select != seg
    vrad[mask] = sme.vrad[mask]
    if sme.cscale_type in ["spline", "spline+mask"]:
        for i in range(len(mask)):
            if (
                mask[i]
                and sme.cscale is not None
                and len(cscale[i]) == len(sme.cscale[i])
            ):
                cscale[i] = sme.cscale[i]
    else:
        cscale[mask] = sme.cscale[mask]

    return cscale, cscale_unc, vrad, vrad_unc
