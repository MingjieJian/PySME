# -*- coding: utf-8 -*-
"""
Calculates the spectrum, based on a set of stellar parameters
And also determines the best fit parameters
"""

import json
import logging
import warnings
from numbers import Number, Real
from os.path import splitext

import numpy as np
from scipy.constants import speed_of_light
from scipy.optimize import OptimizeWarning, least_squares
from scipy.optimize._numdiff import approx_derivative
from scipy.stats import norm
from tqdm import tqdm

from . import __file_ending__
from .abund import Abund
from .atmosphere.atmosphere import AtmosphereError
from .atmosphere.krzfile import KrzFile
from .atmosphere.savfile import SavFile
from .large_file_storage import setup_lfs
from .nlte import DirectAccessFile
from .sme import MASK_VALUES
from .synthesize import Synthesizer
from .util import print_to_log, show_progress_bars

# # Debug usage
# from memory_profiler import profile

logger = logging.getLogger(__name__)

clight = speed_of_light * 1e-3  # km/s
warnings.filterwarnings("ignore", category=OptimizeWarning)


class VariableNumber(np.lib.mixins.NDArrayOperatorsMixin):
    """
    This is essentially a pointer to the value, but will interact with operators
    as though it was the value. This is useful to pass one f_scale parameter
    to the least squares fit, that can then be modified during the iterations.
    """

    def __init__(self, value):
        self.value = value

    _HANDLED_TYPES = (np.ndarray, Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use ArrayLike instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle ArrayLike objects.
            if not isinstance(x, self._HANDLED_TYPES + (VariableNumber,)):
                return NotImplemented

        # Defer to the implementation of the ufunc on unwrapped values.
        inputs = tuple(x.value if isinstance(x, VariableNumber) else x for x in inputs)
        if out:
            kwargs["out"] = tuple(
                x.value if isinstance(x, VariableNumber) else x for x in out
            )
        result = getattr(ufunc, method)(*inputs, **kwargs)

        if type(result) is tuple:
            # multiple return values
            return result
        elif method == "at":
            # no return value
            return None
        else:
            # one return value
            return result


class SME_Solver:
    def __init__(self, filename=None, restore=False):
        self.config, self.lfs_atmo, self.lfs_nlte = setup_lfs()
        self.synthesizer = Synthesizer(
            config=self.config,
            lfs_atmo=self.lfs_atmo,
            lfs_nlte=self.lfs_nlte,
        )

        # Various parameters to keep track of during solving
        self.filename = filename
        self.iteration = 0
        self.parameter_names = []
        self.update_linelist = False
        self.derived_param = None
        self._latest_residual = None
        self._latest_jacobian = None
        self.restore = restore
        self.f_scale = 0.2
        # For displaying the progressbars
        self.progressbar = None
        self.progressbar_jacobian = None

    @property
    def nparam(self):
        return len(self.parameter_names)

    def restore_func(self, sme):
        fname = self.filename.rsplit(".", 1)[0]
        fname = f"{fname}_iter.json"
        try:
            with open(fname) as f:
                data = json.load(f)
            # The keys are string, but we want the max in int, so we need to convert back and forth
            iteration = str(max(int(i) for i in data.keys()))
            for fp in self.parameter_names:
                sme[fp] = data[iteration].get(fp, sme[fp])
            logger.warning(f"Restoring existing backup data from {fname}")
        except:
            pass
        return sme

    def backup(self, sme):
        fname = self.filename.rsplit(".", 1)[0]
        fname = f"{fname}_iter.json"
        try:
            with open(fname) as f:
                data = json.load(f)
        except:
            data = {}
        data[self.iteration] = {fp: sme[fp] for fp in self.parameter_names}
        try:
            with open(fname, "w") as f:
                json.dump(data, f)
        except:
            pass

    def _residuals(
        self, param, sme, spec, uncs, mask, segments="all", isJacobian=False, linelist_mode='all', smelib_lineinfo_mode=0, **_
    ):
        """
        Calculates the synthetic spectrum with sme_func and
        returns the residuals between observation and synthetic spectrum

        residual = (obs - synth) / uncs

        Parameters
        ----------
        param : list(float) of size (n,)
            parameter values to use for synthetic spectrum, order is the same as names
        names : list(str) of size (n,)
            names of the parameters to set, as defined by SME_Struct
        sme : SME_Struct
            sme structure holding all relevant information for the synthetic spectrum generation
        spec : array(float) of size (m,)
            observed spectrum
        uncs : array(float) of size (m,)
            uncertainties of the observed spectrum
        mask : array(bool) of size (k,)
            mask to apply to the synthetic spectrum to select the same points as spec
            The size of the synthetic spectrum is given by sme.wave
            then mask must have the same size, with m True values
        isJacobian : bool, optional
            Flag to use when within the calculation of the Jacobian (default: False)
        fname : str, optional
            filename of the intermediary product (default: "sme.npy")
        fig : Figure, optional
            plotting interface, fig.add(x, y, title) will be called each non jacobian iteration

        Returns
        -------
        resid : array(float) of size (m,)
            residuals of the synthetic spectrum
        """
        update = not isJacobian
        save = not isJacobian and self.filename is not None
        reuse_wavelength_grid = isJacobian
        radial_velocity_mode = "robust" if not isJacobian else "fast"

        # change parameters
        for name, value in zip(self.parameter_names, param):
            if self.derived_param is not None:
                if name in self.derived_param.keys():
                    raise ValueError(f"fitting parameter {name} cannot be also in derived parameter.")
            sme[name] = value
        # change derived parameters
        if self.derived_param is not None:
            for name in self.derived_param.keys():
                if "abund" in name:
                    abund_name = name.split()[1]
                    sme.abund[abund_name] = self.derived_param[name](sme) - sme.monh
                    print(f"Changing derived parameter {name} to {sme.abund[abund_name]:.2f}.")
                else:
                    sme[name] = self.derived_param[name](sme)
                    print(f"Changing derived parameter {name} to {sme[name]:.2f}.")
        # run spectral synthesis
        try:
            result = self.synthesizer.synthesize_spectrum(
                sme,
                updateStructure=update,
                reuse_wavelength_grid=reuse_wavelength_grid,
                segments=segments,
                passLineList=True,
                updateLineList=self.update_linelist,
                radial_velocity_mode=radial_velocity_mode,
                linelist_mode=linelist_mode,
                smelib_lineinfo_mode=smelib_lineinfo_mode,
            )
        except AtmosphereError as ae:
            # Something went wrong (left the grid? Don't go there)
            # If returned value is not finite, the fit algorithm will not go there
            logger.debug(ae)
            return np.full(spec.size, np.inf)

        segments = Synthesizer.check_segments(sme, segments)

        # Get the correct results for the comparison
        synth = sme.synth if update else result[1]
        synth = synth[segments]
        synth = synth[mask] if mask is not None else synth

        if sme.telluric is not None:
            tell = sme.telluric[segments]
            tell = tell[mask] if mask is not None else tell
            synth = synth * tell

        resid = (synth - spec) / (uncs)
        resid = resid.ravel()
        resid = np.nan_to_num(resid, copy=False)

        # Update progress bars
        if isJacobian:
            self.progressbar_jacobian.update(1)
        else:
            self.progressbar.total += 1
            self.progressbar.update(1)

        if not isJacobian:
            # Update f_scale
            # self.f_scale.value = np.percentile(np.abs(resid), 95)
            # print(self.f_scale.value)
            # Save result for jacobian
            self._latest_residual = resid
            self.iteration += 1
        logger.debug("%s", {n: f"{v:.3f}" for n, v in zip(self.parameter_names, param)})

        # Also save intermediary results, because we can
        if save:
            self.backup(sme)

        return resid

    def _jacobian(
        self,
        param,
        *args,
        bounds=None,
        segments="all",
        step_sizes=None,
        method="2-point",
        linelist_mode='all',
        smelib_lineinfo_mode=0,
        **_,
    ):
        """
        Approximate the jacobian numerically
        The calculation is the same as "2-point"
        but we can tell residuals that we are within a jacobian.

        Note that when we reuse the wavelength grid, the results are
        slightly different for reasons(?). Therefore the step size should
        be larger than those differences, which is why we specify the
        step size for each parameter.
        """
        self.progressbar_jacobian.reset()

        # Here we replace the scipy version of approx_derivative with our own
        # The only difference being that we use Multiprocessing for the jacobian
        g = approx_derivative(
            self._residuals,
            param,
            method=method,
            # This feels pretty bad, passing the latest synthetic spectrum
            # by reference as a parameter of the residuals function object
            f0=self._latest_residual,
            abs_step=step_sizes,
            bounds=bounds,
            args=args,
            kwargs={
                "isJacobian": True,
                "segments": segments,
                "linelist_mode": linelist_mode,
                "smelib_lineinfo_mode": smelib_lineinfo_mode,
            },
        )

        if not np.all(np.isfinite(g)):
            g[~np.isfinite(g)] = 0
            logger.warning(
                "Some derivatives are non-finite, setting them to zero. "
                "Final uncertainties will be inaccurate. "
                "You might be running into the boundary of the grid"
            )
        self._latest_jacobian = np.copy(g)
        return g

    def get_bounds(self, sme):
        """
        Create Bounds based on atmosphere grid and general rules

        Note that bounds define by definition a cube in the parameter space,
        but the grid might not be a cube. I.e. Not all combinations of teff, logg, monh are valid
        This method will choose the outerbounds of that space as the boundary, which means that
        we can still run into problems when interpolating the atmospheres

        Parameters
        ----------
        param_names : array(str)
            names of the parameters to vary
        sme : SME_Structure
            sme structure to get bounds for

        Raises
        ------
        IOError
            If the atmosphere file can't be read, allowed types
            are IDL savefiles (.sav), and .krz files

        Returns
        -------
        bounds : dict
            Bounds for the given parameters
        """

        bounds = {}

        # Create bounds based on atmosphere grid
        if (
            "teff" in self.parameter_names
            or "logg" in self.parameter_names
            or "monh" in self.parameter_names
        ):
            if sme.atmo.method == "grid":
                atmo_source = sme.atmo.source
                _, ext = splitext(atmo_source)
                atmo_file = self.lfs_atmo.get(atmo_source)

                if ext == ".sav":
                    atmo_grid = SavFile(
                        atmo_file, source=atmo_source, lfs=self.lfs_atmo
                    )

                    teff = np.unique(atmo_grid.teff)
                    teff = np.min(teff), np.max(teff)
                    bounds["teff"] = teff

                    logg = np.unique(atmo_grid.logg)
                    logg = np.min(logg), np.max(logg)
                    bounds["logg"] = logg

                    monh = np.unique(atmo_grid.monh)
                    monh = np.min(monh), np.max(monh)
                    bounds["monh"] = monh
                elif ext == ".krz":
                    # krz atmospheres are fixed to one parameter set
                    # allow just "small" area around that
                    atmo = KrzFile(atmo_file, source=atmo_source)
                    bounds["teff"] = atmo.teff - 500, atmo.teff + 500
                    bounds["logg"] = atmo.logg - 1, atmo.logg + 1
                    bounds["monh"] = atmo.monh - 1, atmo.monh + 1
                else:
                    raise IOError(f"File extension {ext} not recognized")
            if sme.atmo.method == "embedded":
                atmo = sme.atmo
                bounds["teff"] = atmo.teff - 500, atmo.teff + 500
                bounds["logg"] = atmo.logg - 1, atmo.logg + 1
                bounds["monh"] = atmo.monh - 1, atmo.monh + 1
        # Add generic bounds
        bounds.update(
            {
                "vmic": [0, clight],
                "vmac": [0, clight],
                "vsini": [0, clight],
                "ipres": [1, 500_000],
            }
        )
        # bounds.update({"abund %s" % el: [-10, 11] for el in abund_elem})

        result = np.array([[-np.inf, np.inf]] * self.nparam)
        solar = Abund.solar()
        for i, name in enumerate(self.parameter_names):
            if name[:5].lower() == "abund":
                element = name[5:].strip().capitalize()
                if element in sme.nlte.elements:
                    fname = sme.nlte.grids[element]
                    fname = self.lfs_nlte.get(fname)
                    grid = DirectAccessFile(fname)
                    available = grid["abund"]
                    xmin, xmax = available.min(), available.max()
                    xmin += solar[element]
                    xmax += solar[element]
                    if xmin == xmax:
                        xmin -= 1
                        xmax += 1
                    result[i] = [xmin, xmax]
                else:
                    result[i] = [-10, 11]
            elif name[:8].lower() == "linelist":
                if "excit" in name.lower():
                    value = sme.linelist[int(name.split()[1])]["excit"][0]
                    result[i] = [value - 0.005, value + 0.005]
                else:
                    result[i] = [-np.inf, np.inf]
            else:
                result[i] = bounds[name]

        result = result.T

        if len(result) > 0:
            return result
        else:
            return [-np.inf, np.inf]

    def get_scale(self):
        """
        Returns scales for each parameter so that values are on order ~1

        Parameters
        ----------
        param_names : list(str)
            names of the parameters

        Returns
        -------
        scales : list(float)
            scales of the parameters in the same order as input array
        """

        # The only parameter we want to scale right now is temperature,
        # as it is orders of magnitude larger than all others
        scales = {"teff": 1000}
        scales = [
            scales[name] if name in scales.keys() else 1
            for name in self.parameter_names
        ]
        return scales

    def get_step_sizes(self, parameter_names):
        steps = {
            "teff": 10,
            "logg": 0.01,
            "monh": 0.01,
            "vmic": 0.01,
            "vmac": 0.1,
            "vsini": 0.05,
            "vrad": 0.05,
            "gam6": 0.02,
            "ipres": 1000,
        }
        step_sizes = []
        for param in parameter_names:
            if param in steps.keys():
                step_sizes += [steps[param]]
            elif param.startswith("Abund"):
                step_sizes += [0.01]
            elif param.startswith("linelist"):
                if "excit" in param:
                    step_sizes += [0.001]
                else:
                    step_sizes += [0.2]
            else:
                step_sizes += [0.001]
        step_sizes = np.asarray(step_sizes)
        return step_sizes

    def get_default_values(self, sme):
        """Default parameter values for each name"""
        d = {"teff": 5778, "logg": 4.4, "monh": 0, "vmac": 1, "vmic": 1}
        d.update({f"{el} abund": v for el, v in Abund.solar()().items()})

        def default(name):
            logger.info("No value for %s set, using default value %s", name, d[name])
            return d[name]

        values = [
            sme[s] if sme[s] is not None else default(s) for s in self.parameter_names
        ]
        return np.array(values)

    @staticmethod
    def estimate_uncertainties(resid, deriv):
        """
        Estimate the uncertainties by fitting the cumulative distribution of
        derivative / uncertainties vs. residual / derivative
        with the generalized normal distribution and use the 68% percentile
        as the 1 sigma approximation for a normally distributed variable

        Parameters
        ----------
        unc : array of shape (n,)
            uncertainties
        resid : array of shape (n,)
            residuals of the least squares fit
        deriv : array of shape (n, p)
            derivatives (jacobian) of the least squares fit for each parameter

        Returns
        -------
        freep_unc : array of shape (p,)
            uncertainties for each free paramater, in the same order as self.parameter_names
        """

        nparameters = deriv.shape[1]
        freep_unc = np.zeros(nparameters)

        for i in range(nparameters):
            pder = deriv[:, i]
            idx = pder != 0
            idx &= np.isfinite(pder)

            if np.count_nonzero(idx) <= 5:
                logger.warning(
                    "Not enough data points with a suitable derivative "
                    "to determine the uncertainties"
                )
                continue
            # Sort pixels according to the change of the i
            # parameter needed to match the observations
            idx_sort = np.argsort(resid[idx] / pder[idx])
            ch_x = resid[idx][idx_sort] / pder[idx][idx_sort]
            # Weights of the individual pixels also sorted
            # uncertainties are already included in pder / unc[idx][idx_sort]
            ch_y = np.abs(pder[idx][idx_sort])
            # Cumulative weights
            ch_y = np.cumsum(ch_y)
            # Normalized cumulative weights
            ch_y /= ch_y[-1]

            # hmed = np.interp(0.5, ch_y, ch_x)
            interval = np.interp([0.16, 0.84], ch_y, ch_x)
            sigma_estimate = (interval[1] - interval[0]) / 2
            freep_unc[i] = sigma_estimate

        return freep_unc

    def update_fitresults(self, sme, result, segments):
        # Update SME structure
        sme.fitresults.clear()

        popt = result.x
        sme.fitresults.values = popt
        sme.fitresults.parameters = self.parameter_names

        # Determine the covariance
        # hessian == fisher information matrix
        fisher = result.jac.T.dot(result.jac)
        covar = np.linalg.pinv(fisher)
        sig = np.sqrt(covar.diagonal())

        # Update fitresults
        sme.fitresults.covariance = covar
        sme.fitresults.gradient = result.grad
        sme.fitresults.derivative = result.jac
        sme.fitresults.residuals = (
            (sme.spec[segments] - sme.synth[segments]) / sme.uncs[segments]
        )[sme.mask_line[segments]]
        sme.fitresults.chisq = np.sum(sme.fitresults.residuals ** 2) / (
            result.fun.size - len(self.parameter_names)
        )
        sme.fitresults.iterations = self.iteration

        sme.fitresults.fit_uncertainties = [np.nan for _ in self.parameter_names]
        for i in range(len(self.parameter_names)):
            # Errors based on covariance matrix
            sme.fitresults.fit_uncertainties[i] = sig[i] * np.sqrt(sme.fitresults.chisq)

        try:
            sme.fitresults.uncertainties = self.estimate_uncertainties(
                sme.fitresults.residuals,
                sme.fitresults.derivative,
            )
        except:
            logger.warning(
                "Could not determine the uncertainties from the probability "
                "distribution, using fit uncertainties instead"
            )
            sme.fitresults.uncertainties = sme.fitresults.fit_uncertainties
        # sme.fitresults.uncertainties = sme.fitresults.fit_uncertainties

        return sme

    def sanitize_parameter_names(self, sme, param_names):
        # Sanitize parameter names
        param_names = [p.casefold() for p in param_names]
        param_names = [p.capitalize() if p[:5] == "abund" else p for p in param_names]

        param_names = [p if p != "grav" else "logg" for p in param_names]
        param_names = [p if p != "feh" else "monh" for p in param_names]

        # Parameters are unique
        # But keep the order the same
        param_names, index = np.unique(param_names, return_index=True)
        param_names = param_names[np.argsort(index)]
        param_names = list(param_names)

        if "vrad" in param_names:
            param_names.remove("vrad")
            if sme.vrad_flag in ["fix", "none"]:
                sme.vrad_flag = "whole"
                logger.info(
                    f"Removed fit parameter 'vrad', instead set radial velocity flag to {sme.vrad_flag}"
                )

        if "cont" in param_names:
            param_names.remove("cont")
            if sme.cscale_flag in ["fix", "none"]:
                sme.cscale_flag = "linear"
                logger.info(
                    "Removed fit parameter 'cont', instead set continuum flag to %s",
                    sme.cscale_flag,
                )
        return param_names

    def solve(
        self,
        sme,
        param_names=None,
        segments="all",
        bounds=None,
        step_sizes=None,
        derived_param=None,
        dynamic_param=None,
        linelist_mode="all",
        smelib_lineinfo_mode=0,
    ):
        """
        Find the least squares fit parameters to an observed spectrum

        NOTE: intermediary results will be saved in filename ("sme.npy")

        Parameters
        ----------
        sme : SME_Struct
            sme struct containing all input (and output) parameters
        param_names : list, optional
            the names of the parameters to fit (default: ["teff", "logg", "monh"])
        filename : str, optional
            the sme structure will be saved to this file, use None to suppress this behaviour (default: "sme.npy")

        Returns
        -------
        sme : SME_Struct
            same sme structure with fit results in sme.fitresults, and best fit spectrum in sme.smod
        """

        assert "wave" in sme, "SME Structure has no wavelength"
        assert "spec" in sme, "SME Structure has no observation"

        if derived_param is not None and dynamic_param is not None and derived_param is not dynamic_param:
            raise ValueError("Specify only one of derived_param or dynamic_param, not both.")
        if dynamic_param is not None:
            warnings.warn(
                "'dynamic_param' is deprecated and will be removed in a future release; use 'derived_param' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
        self.derived_param = derived_param if derived_param is not None else dynamic_param

        if self.restore and self.filename is not None:
            fname = self.filename.rsplit(".", 1)[0]
            fname = f"{fname}_iter.json"
            try:
                with open(fname) as f:
                    data = json.load(f)
                for fp in param_names:
                    sme[fp] = data[fp]
                logger.warning(f"Restoring existing backup data from {fname}")
            except:
                pass

        if "uncs" not in sme:
            sme.uncs = np.ones(sme.spec.size)
            logger.warning("SME Structure has no uncertainties, using all ones instead")
        if "mask" not in sme:
            sme.mask = np.full(sme.wave.size, MASK_VALUES.LINE)

        segments = Synthesizer.check_segments(sme, segments)

        # Clean parameter values
        if param_names is None:
            param_names = sme.fitparameters
        if param_names is None or len(param_names) == 0:
            logger.warning(
                "No Fit Parameters have been set. Using ('teff', 'logg', 'monh') instead."
            )
            param_names = ("teff", "logg", "monh")
        self.parameter_names = self.sanitize_parameter_names(sme, param_names)

        self.update_linelist = False
        for name in self.parameter_names:
            if name[:8] == "linelist":
                if self.update_linelist is False:
                    self.update_linelist = []
                try:
                    idx = int(name.split()[1])
                except IndexError:
                    raise ValueError(
                        f"Could not parse fit parameter {name}, expected a "
                        "linelist parameter like 'linelist n gflog'"
                    )
                self.update_linelist += [idx]
        if self.update_linelist:
            self.update_linelist = np.unique(self.update_linelist)

        # Create appropiate bounds
        if bounds is None:
            bounds = self.get_bounds(sme)
        # scales = self.get_scale()
        if step_sizes is None:
            step_sizes = self.get_step_sizes(self.parameter_names)
        # Starting values
        p0 = self.get_default_values(sme)
        if np.any((p0 < bounds[0]) | (p0 > bounds[1])):
            logger.warning(
                "Initial values are incompatible with the bounds, clipping initial values"
            )
            p0 = np.clip(p0, bounds[0], bounds[1])
        # Restore backup
        if self.restore:
            sme = self.restore_func(sme)

        # Get constant data from sme structure
        for seg in segments:
            sme.mask[seg, sme.uncs[seg] == 0] = MASK_VALUES.BAD
        mask = sme.mask_line[segments]
        spec = sme.spec[segments][mask]
        uncs = sme.uncs[segments][mask]

        # Divide the uncertainties by the spectrum, to improve the fit in the continuum
        # Just as in IDL SME, this increases the relative error for points inside lines
        # uncs /= np.abs(spec)

        # This is the expected range of the uncertainty
        # if the residuals are larger, they are dampened by log(1 + z)
        self.f_scale = 0.2 * np.nanmean(spec.ravel()) / np.nanmean(uncs.ravel())

        logger.info("Fitting Spectrum with Parameters: %s", ",".join(param_names))
        logger.debug("Initial values: %s", p0)
        logger.debug("Bounds: %s", bounds)

        if (
            sme.wran.min() * (1 - 100 / clight) > sme.linelist.wlcent.min()
            or sme.wran.max() * (1 + 100 / clight) < sme.linelist.wlcent.max()
        ):
            logger.warning(
                "The linelist extends far beyond the requested wavelength range."
                " This will slow down the calculation, consider using only relevant lines\n"
                f"Wavelength range: {sme.wran.min()} - {sme.wran.max()} Å"
                f" ; Linelist range: {sme.linelist.wlcent.min()} - {sme.linelist.wlcent.max()} Å"
            )

        # Setup LineList only once (Mingjie: not needed?)
        # dll = self.synthesizer.get_dll()
        # dll.SetLibraryPath()
        # _ = dll.InputLineList(sme.linelist)

        # Do the heavy lifting
        if self.nparam > 0:
            self.progressbar = tqdm(
                desc="Iteration", total=0, disable=~show_progress_bars
            )
            self.progressbar_jacobian = tqdm(
                desc="Jacobian", total=len(p0), disable=~show_progress_bars
            )
            with print_to_log():
                res = least_squares(
                    self._residuals,
                    jac=self._jacobian,
                    x0=p0,
                    bounds=bounds,
                    loss=sme.leastsquares_loss,
                    f_scale=self.f_scale,
                    method=sme.leastsquares_method,
                    x_scale=sme.leastsquares_xscale,
                    # These control the tolerance, for early termination
                    # since each iteration is quite expensive
                    xtol=sme.accxt,
                    ftol=sme.accft,
                    gtol=sme.accgt,
                    verbose=2,
                    args=(sme, spec, uncs, mask),
                    kwargs={
                        "bounds": bounds,
                        "segments": segments,
                        "step_sizes": step_sizes,
                        "method": sme.leastsquares_jac,
                        "linelist_mode": linelist_mode,
                        "smelib_lineinfo_mode": smelib_lineinfo_mode,
                    },
                )
                # The jacobian is altered by the loss function
                # This lets us keep the original for our uncertainty estimate
                res.jac = self._latest_jacobian

            self.progressbar.close()
            self.progressbar_jacobian.close()
            for i, name in enumerate(self.parameter_names):
                sme[name] = res.x[i]
            sme = self.update_fitresults(sme, res, segments)
            logger.debug("Reduced chi square: %.3f", sme.fitresults.chisq)
            try:
                for name, value, unc in zip(
                    self.parameter_names, res.x, sme.fitresults.fit_uncertainties
                ):
                    logger.info("%s\t%.5f +- %.5g", name.ljust(10), value, unc)
                logger.info("%s\t%s +- %s", "vrad".ljust(10), sme.vrad, sme.vrad_unc)
            except:
                pass
        elif len(param_names) > 0:
            # This happens when vrad and/or cscale are given as parameters but nothing else
            # We could try to reuse the already calculated synthetic spectrum (if it already exists)
            # However it is much lower resolution then the newly synthethized one (usually)
            # Therefore the radial velocity wont be as good as when redoing the whole thing
            sme = self.synthesizer.synthesize_spectrum(
                sme,
                segments,
                linelist_mode=linelist_mode,
                smelib_lineinfo_mode=smelib_lineinfo_mode,
            )
        else:
            raise ValueError("No fit parameters given")

        if self.filename is not None:
            sme.save(self.filename)

        return sme


def solve(
    sme,
    param_names=None,
    segments="all",
    filename=None,
    restore=False,
    derived_param=None,
    dynamic_param=None,
    **kwargs,
):
    solver = SME_Solver(filename=filename, restore=restore)
    return solver.solve(
        sme,
        param_names,
        segments,
        derived_param=derived_param,
        dynamic_param=dynamic_param,
        **kwargs,
    )
