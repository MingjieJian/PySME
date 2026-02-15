# -*- coding: utf-8 -*-
"""
Spectral Synthesis Module of SME
"""
import logging
import uuid
import warnings

import numpy as np
from scipy.constants import speed_of_light
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.ndimage import convolve
from tqdm import tqdm
from scipy.spatial.distance import cdist

from . import broadening
from .atmosphere.interpolation import AtmosphereInterpolator
from .continuum_and_radial_velocity import (
    apply_radial_velocity_and_continuum,
    match_rv_continuum,
    null_result,
)
from .iliffe_vector import Iliffe_vector
from .large_file_storage import setup_lfs
from .sme import MASK_VALUES
from .sme_synth import SME_DLL
from .util import show_progress_bars, boundary_vertices, safe_interpolation, interpolate_3DNLTEH_spectrum_RBF
from . import util
from .sme import SME_Structure

from contextlib import redirect_stdout
from copy import deepcopy
from pqdm.processes import pqdm
# from pqdm import threads

# from memory_profiler import profile

import os
import re
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

# Temp solution
import pandas as pd
pd.options.mode.chained_assignment = None  # None means no warning will be shown

logger = logging.getLogger(__name__)

clight = speed_of_light * 1e-3  # km/s

__DLL_DICT__ = {}
__DLL_IDS__ = {}

def _extract_params_from_fname(fname):
    match = re.match(r'teff([0-9.]+)_logg([0-9.]+)_monh([-0-9.]+)_vmic([0-9.]+)\.npz', fname)
    if match:
        return tuple(map(float, match.groups()))
    return None

def _load_all_grid_points(cdr_folder):
    param_list = []
    fname_map = {}
    for fname in os.listdir(cdr_folder):
        if fname.endswith('.npz'):
            key = _extract_params_from_fname(fname)
            if key:
                param_list.append(key)
                fname_map[key] = fname
    return np.array(param_list), fname_map

def _load_mask_and_ranges_from_compressed(full_path, wlcent, n_lines_total):
    """
    Load the strong-line mask and wavelength ranges from a compressed CDR grid file.

    The compressed file is expected to contain:
        - mask_bits     : 1D uint8 array, bit-packed boolean mask of strong lines
        - unique_widths : 1D float32 array, all distinct line_width values in this grid
        - codes         : 1D uint8/uint16/uint32 array, indices into unique_widths
                          corresponding to the strong lines. The order of `codes`
                          matches the order of strong-line indices given by:
                          strong_idx = np.nonzero(mask_full)[0]

    For each strong line i:
        width[i] = unique_widths[codes[i]]
        line_range_s[i] = wlcent[i] - 0.5 * width[i]
        line_range_e[i] = wlcent[i] + 0.5 * width[i]

    All non-strong lines are assigned NaN in line_range_s/e.

    Parameters
    ----------
    full_path : str
        Path to the compressed npz file for a single grid point.
    wlcent : array-like, shape (n_lines_total,)
        Central wavelengths of all lines in the global line list.
    n_lines_total : int
        Total number of lines in the line list.

    Returns
    -------
    mask_full : np.ndarray, dtype=bool, shape (n_lines_total,)
        Boolean array marking which lines are strong in this grid.
    line_range_s : np.ndarray, dtype=float32, shape (n_lines_total,)
        Start wavelength of the line range for each line (NaN for non-strong lines).
    line_range_e : np.ndarray, dtype=float32, shape (n_lines_total,)
        End wavelength of the line range for each line (NaN for non-strong lines).
    """
    data = np.load(full_path)
    mask_bits     = data["mask_bits"]
    unique_widths = data["unique_widths"]
    codes         = data["codes"]
    n_lines_total = int(data["n_lines_total"])

    # 1) Unpack bit-packed mask into a full boolean array
    mask_full = np.unpackbits(mask_bits).astype(bool)[:n_lines_total]
    if mask_full.size < n_lines_total:
        raise ValueError(
            f"{full_path}: unpacked mask length {mask_full.size} < n_lines_total={n_lines_total}. "
            "Please ensure the compressed database matches the current VALD linelist."
        )
    mask_full = mask_full[:n_lines_total]

    # 2) Indices of strong lines and their corresponding widths
    strong_idx = np.nonzero(mask_full)[0]
    widths     = unique_widths[codes]  # one-to-one with strong_idx

    # 3) Reconstruct line_range_s/e from central wavelengths and widths
    wlcent = np.asarray(wlcent)
    if wlcent.size != n_lines_total:
        raise ValueError(
            f"{full_path}: wlcent length {wlcent.size} != n_lines_total={n_lines_total}"
        )

    wl_strong = wlcent[strong_idx]
    s_vals = (wl_strong - 0.5 * widths).astype(np.float32)
    e_vals = (wl_strong + 0.5 * widths).astype(np.float32)

    line_range_s = np.full(n_lines_total, np.nan, dtype=np.float32)
    line_range_e = np.full(n_lines_total, np.nan, dtype=np.float32)
    line_range_s[strong_idx] = s_vals
    line_range_e[strong_idx] = e_vals

    return mask_full, line_range_s, line_range_e

class Synthesizer:
    def __init__(self, config=None, lfs_atmo=None, lfs_nlte=None, dll=None):
        self.config, self.lfs_atmo, self.lfs_nlte = setup_lfs(
            config, lfs_atmo, lfs_nlte
        )
        # dict: internal storage of the adaptive wavelength grid
        self.wint = {}
        # dll: the smelib object used for the radiative transfer calculation
        self.dll = dll if dll is not None else SME_DLL()
        self.dll = self.get_dll_id()
        self.atmosphere_interpolator = None
        # This stores a reference to the currently used sme structure, so we only log it once
        self.known_sme = None
        self.update_cdr_switch = False
        # logger.info("Don't forget to cite your sources. Use sme.citation()")

    def get_atmosphere(self, sme):
        """
        Return an atmosphere based on specification in an SME structure

        sme.atmo.method defines mode of action:
            "grid"
                interpolate on atmosphere grid
            "embedded"
                No change
            "routine"
                calls sme.atmo.source(sme, atmo)

        Parameters
        ---------
            sme : SME_Struct
                sme structure with sme.atmo = atmosphere specification

        Returns
        -------
        sme : SME_Struct
            sme structure with updated sme.atmo
        """

        # Handle atmosphere grid or user routine.
        atmo = sme.atmo

        if atmo.method == "grid":
            if self.atmosphere_interpolator is None:
                self.atmosphere_interpolator = AtmosphereInterpolator(
                    depth=atmo.depth,
                    interp=atmo.interp,
                    geom=atmo.geom,
                    lfs_atmo=self.lfs_atmo,
                )
            else:
                self.atmosphere_interpolator.depth = atmo.depth
                self.atmosphere_interpolator.interp = atmo.interp
                self.atmosphere_interpolator.geom = atmo.geom

            atmo = self.atmosphere_interpolator.interp_atmo_grid(
                atmo.source, sme.teff, sme.logg, sme.monh
            )
        elif atmo.method == "routine":
            atmo = atmo.source(sme, atmo)
        elif atmo.method == "embedded":
            # atmo structure already extracted in sme_main
            pass
        else:
            raise AttributeError("Source must be 'grid', 'routine', or 'embedded'")

        sme.atmo = atmo
        return sme

    @staticmethod
    def get_wavelengthrange(wran, vrad, vsini):
        """
        Determine wavelengthrange that needs to be calculated
        to include all lines within velocity shift vrad + vsini
        """
        # 30 km/s == maximum barycentric velocity
        vrad_pad = 30.0 + 0.5 * np.clip(vsini, 0, None)  # km/s
        vbeg = vrad_pad + np.clip(vrad, 0, None)  # km/s
        vend = vrad_pad - np.clip(vrad, None, 0)  # km/s

        wbeg = wran[0] * (1 - vbeg / clight)
        wend = wran[1] * (1 + vend / clight)
        return wbeg, wend

    @staticmethod
    def new_wavelength_grid(wint):
        """Generate new wavelength grid within bounds of wint"""
        # Determine step size for a new model wavelength scale, which must be uniform
        # to facilitate convolution with broadening kernels. The uniform step size
        # is the larger of:
        #
        # [1] smallest wavelength step in WINT_SEG, which has variable step size
        # [2] 10% the mean dispersion of WINT_SEG
        # [3] 0.05 km/s, which is 1% the width of solar line profiles

        wbeg, wend = wint[0], wint[-1]
        wmid = 0.5 * (wend + wbeg)  # midpoint of segment
        wspan = wend - wbeg  # width of segment
        diff = wint[1:] - wint[:-1]
        jmin = np.argmin(diff)
        vstep1 = diff[jmin] / wint[jmin] * clight  # smallest step
        # vstep2 = 0.1 * wspan / (len(wint) - 1) / wmid * clight  # 10% mean dispersion
        # vstep2 = np.median(diff / wint[:-1] * clight) # Median step
        vstep3 = 0.05  # 0.05 km/s step
        vstep = np.max([vstep1, vstep3])  # select the largest
        # vstep = 0.02

        # Generate model wavelength scale X, with uniform wavelength step.
        nx = int(
            np.abs(np.log10(wend / wbeg)) / np.log10(1 + vstep / clight) + 1
        )  # number of wavelengths
        if nx % 2 == 0:
            nx += 1  # force nx to be odd

        # Resolution
        # IDL way
        # resol_out = 1 / ((wend / wbeg) ** (1 / (nx - 1)) - 1)
        # vstep = clight / resol_out
        # x_seg = wbeg * (1 + 1 / resol_out) ** np.arange(nx)

        # Python way (not identical, as IDL endpoint != wend)
        # difference approx 1e-7
        x_seg = np.geomspace(wbeg, wend, num=nx)
        resol_out = 1 / np.diff(np.log(x_seg[:2]))[0]
        vstep = clight / resol_out
        return x_seg, vstep

    @staticmethod
    def check_segments(sme, segments):
        if isinstance(segments, str) and segments == "all":
            segments = range(sme.nseg)
        else:
            segments = np.atleast_1d(segments)
            if np.any(segments < 0) or np.any(segments >= sme.nseg):
                raise IndexError("Segment(s) out of range")
            segments = np.unique(segments)

        if sme.mask is not None:
            segments = [
                seg for seg in segments if not np.all(sme.mask[seg] == MASK_VALUES.BAD)
            ]
        return segments

    @staticmethod
    def apply_radial_velocity_and_continuum_synth(
        wave, spec, wmod, smod, cmod, vrad, cscale, cscale_type, segments
    ):
        smod = apply_radial_velocity_and_continuum(
            wave, wmod, smod, vrad, cscale, cscale_type, segments
        )
        cmod = apply_radial_velocity_and_continuum(
            wave, wmod, cmod, vrad, None, None, segments
        )
        return smod, cmod

    @staticmethod
    def integrate_flux(mu, inten, deltav, vsini, vrt, osamp=None, wt=None, tdnlteH=False):
        """
        Produces a flux profile by integrating intensity profiles (sampled
        at various mu angles) over the visible stellar surface.

        Intensity profiles are weighted by the fraction of the projected
        stellar surface they represent, apportioning the area between
        adjacent MU points equally. Additional weights (such as those
        used in a Gauss-Legendre quadrature) can not meaningfully be
        used in this scheme.  About twice as many points are required
        with this scheme to achieve the precision of Gauss-Legendre
        quadrature.
        DELTAV, VSINI, and VRT must all be in the same units (e.g. km/s).
        If specified, OSAMP should be a positive integer.

        Parameters
        ----------
        mu : array(float) of size (nmu,)
            cosine of the angle between the outward normal and
            the line of sight for each intensity spectrum in INTEN.
        inten : array(float) of size(nmu, npts)
            intensity spectra at specified values of MU.
        deltav : float
            velocity spacing between adjacent spectrum points
            in INTEN (same units as VSINI and VRT).
        vsini : float
            maximum radial velocity, due to solid-body rotation.
        vrt : float
            radial-tangential macroturbulence parameter, i.e.
            np.sqrt(2) times the standard deviation of a Gaussian distribution
            of turbulent velocities. The same distribution function describes
            the radial motions of one component and the tangential motions of
            a second component. Each component covers half the stellar surface.
            See 'The Observation and Analysis of Stellar Photospheres', Gray.
        osamp : int, optional
            internal oversampling factor for convolutions.
            By default convolutions are done using the input points (OSAMP=1),
            but when OSAMP is set to higher integer values, the input spectra
            are first oversampled by cubic spline interpolation.

        Returns
        -------
        value : array(float) of size (npts,)
            Disk integrated flux profile.

        Note
        ------------
            If you use this algorithm in work that you publish, please cite
            Valenti & Anderson 1996, PASP, currently in preparation.
        """
        """
        History
        -----------
        Feb-88  GM
            Created ANA version.
        13-Oct-92 JAV
            Adapted from G. Marcy's ANA routi!= of the same name.
        03-Nov-93 JAV
            Switched to annular convolution technique.
        12-Nov-93 JAV
            Fixed bug. Intensity compo!=nts not added when vsini=0.
        14-Jun-94 JAV
            Reformatted for "public" release. Heavily commented.
            Pass deltav instead of 2.998d5/deltav. Added osamp
            keyword. Added rebinning logic at end of routine.
            Changed default osamp from 3 to 1.
        20-Feb-95 JAV
            Added mu as an argument to handle arbitrary mu sampling
            and remove ambiguity in intensity profile ordering.
            Interpret VTURB as np.sqrt(2)*sigma instead of just sigma.
            Replaced call_external with call to spl_{init|interp}.
        03-Apr-95 JAV
            Multiply flux by pi to give observed flux.
        24-Oct-95 JAV
            Force "nmk" padding to be at least 3 pixels.
        18-Dec-95 JAV
            Renamed from dskint() to rtint(). No longer make local
            copy of intensities. Use radial-tangential instead
            of isotropic Gaussian macroturbulence.
        26-Jan-99 JAV
            For NMU=1 and VSINI=0, assume resolved solar surface#
            apply R-T macro, but supress vsini broadening.
        01-Apr-99 GMH
            Use annuli weights, rather than assuming ==ual area.
        07-Mar-12 JAV
            Force vsini and vmac to be scalars.
        """

        # Make local copies of various input variables, which will be altered below.
        # Force vsini and especially vmac to be scalars. Otherwise mu dependence fails.

        if np.size(vsini) > 1:
            vsini = vsini[0]
        if np.size(vrt) > 1:
            vrt = vrt[0]
        nmu = np.size(mu)  # number of radii

        # Convert input MU to projected radii, R, of annuli for a star of unit radius
        #  (which is just sine, rather than cosine, of the angle between the outward
        #  normal and the line of sight).
        rmu = np.sqrt(1 - mu ** 2)  # use simple trig identity
        if nmu > 1:
            r = np.sqrt(
                0.5 * (rmu[:-1] ** 2 + rmu[1:] ** 2)
            )  # area midpoints between rmu
            r = np.concatenate(([0], r, [1]))
        else:
            r = np.array([0, 1])

        # Determine oversampling factor.
        if osamp is None:
            if vsini == 0:
                os = 2
            else:
                os = deltav / (vsini * r[r != 0])
                os = np.max(os[np.isfinite(os)])
                os = int(np.ceil(os)) + 1
        else:
            os = osamp
        # force integral value > 1
        os = round(np.clip(os, 2, 10))

        # Sort the projected radii and corresponding intensity spectra into ascending
        #  order (i.e. from disk center to the limb), which is equivalent to sorting
        #  MU in descending order.
        isort = np.argsort(rmu)
        rmu = rmu[isort]  # reorder projected radii
        if nmu == 1 and vsini != 0:
            logger.warning(
                "Vsini is non-zero, but only one projected radius (mu value) is set. No rotational broadening will be performed."
            )
            vsini = 0  # ignore vsini if only 1 mu

        # Calculate projected radii for boundaries of disk integration annuli.  The n+1
        # boundaries are selected such that r(i+1) exactly bisects the area between
        # rmu(i) and rmu(i+1). The in!=rmost boundary, r(0) is set to 0 (disk center)
        # and the outermost boundary, r(nmu) is set to 1 (limb).
        if wt is None:
            if nmu > 1:  # really want disk integration
                # Calculate integration weights for each disk integration annulus.  The weight
                # is just given by the relative area of each annulus, normalized such that
                # the sum of all weights is unity.  Weights for limb darkening are included
                # explicitly in the intensity profiles, so they aren't needed here.
                wt = r[1:] ** 2 - r[:-1] ** 2  # weights = relative areas
            else:
                wt = np.array([1.0])  # single mu value, full weight

        # Generate index vectors for input and oversampled points. Note that the
        # oversampled indicies are carefully chosen such that every "os" finely
        # sampled points fit exactly into one input bin. This makes it simple to
        # "integrate" the finely sampled points at the end of the routine.
        npts = inten.shape[1]  # number of points
        xpix = np.arange(npts, dtype=float)  # point indices
        nfine = os * npts  # number of oversampled points
        xfine = (0.5 / os) * (
            2 * np.arange(nfine, dtype=float) - os + 1
        )  # oversampled points indices

        # Loop through annuli, constructing and convolving with rotation kernels.

        yfine = np.empty(nfine)  # init oversampled intensities
        flux = np.zeros(nfine)  # init flux vector
        for imu in range(nmu):  # loop thru integration annuli

            #  Use external cubic spline routine (adapted from Numerical Recipes) to make
            #  an oversampled version of the intensity profile for the current annulus.
            ypix = inten[isort[imu]]  # extract intensity profile
            if os == 1:
                # just copy (use) original profile
                yfine = np.copy(ypix)
            else:
                # spline onto fine wavelength scale
                try:
                    cs = CubicSpline(
                        xpix, ypix, extrapolate=True
                    )
                    yfine = cs(xfine)
                except ValueError:
                    yfine = interp1d(
                        xpix, ypix, kind="linear", fill_value="extrapolate"
                    )(xfine)

            # Construct the convolution kernel which describes the distribution of
            # rotational velocities present in the current annulus. The distribution has
            # been derived analytically for annuli of arbitrary thickness in a rigidly
            # rotating star. The kernel is constructed in two pieces: o!= piece for
            # radial velocities less than the maximum velocity along the inner edge of
            # the annulus, and one piece for velocities greater than this limit.
            if vsini > 0:
                # nontrivial case
                r1 = r[imu]  # inner edge of annulus
                r2 = r[imu + 1]  # outer edge of annulus
                dv = deltav / os  # oversampled velocity spacing
                maxv = vsini * r2  # maximum velocity in annulus
                nrk = 2 * int(maxv / dv) + 3  ## oversampled kernel point
                # velocity scale for kernel
                v = dv * (np.arange(nrk, dtype=float) - ((nrk - 1) / 2))
                rkern = np.zeros(nrk)  # init rotational kernel
                j1 = np.abs(v) < vsini * r1  # low velocity points
                rkern[j1] = np.sqrt((vsini * r2) ** 2 - v[j1] ** 2) - np.sqrt(
                    (vsini * r1) ** 2 - v[j1] ** 2
                )  # generate distribution

                j2 = (np.abs(v) >= vsini * r1) & (np.abs(v) <= vsini * r2)
                rkern[j2] = np.sqrt(
                    (vsini * r2) ** 2 - v[j2] ** 2
                )  # generate distribution

                rkern = rkern / np.sum(rkern)  # normalize kernel

                # Convolve the intensity profile with the rotational velocity kernel for this
                # annulus. Pad each end of the profile with as many points as are in the
                # convolution kernel. This reduces Fourier ringing.
                yfine = convolve(yfine, rkern, mode="nearest")

            # Calculate projected sigma for radial and tangential velocity distributions.
            muval = mu[isort[imu]]  # current value of mu
            sigma = os * vrt / np.sqrt(2) / deltav  # standard deviation in points
            sigr = sigma * muval  # reduce by current mu value
            sigt = sigma * np.sqrt(1.0 - muval ** 2)  # reduce by np.sqrt(1-mu**2)

            # Figure out how many points to use in macroturbulence kernel.
            nmk = int(10 * sigma)
            nmk = np.clip(nmk, 3, (nfine - 3) // 2)

            # Construct radial macroturbulence kernel with a sigma of mu*VRT/np.sqrt(2).
            if sigr > 0:
                xarg = np.linspace(-nmk, nmk, 2 * nmk + 1) / sigr
                xarg = np.clip(-0.5 * xarg ** 2, -20, None)
                mrkern = np.exp(xarg)  # compute the gaussian
                mrkern = mrkern / np.sum(mrkern)  # normalize the profile
            else:
                mrkern = np.zeros(2 * nmk + 1)  # init with 0d0
                mrkern[nmk] = 1.0  # delta function

            # Construct tangential kernel with a sigma of np.sqrt(1-mu**2)*VRT/np.sqrt(2).
            if sigt > 0:
                xarg = np.linspace(-nmk, nmk, 2 * nmk + 1) / sigt
                xarg = np.clip(-0.5 * xarg ** 2, -20, None)
                mtkern = np.exp(xarg)  # compute the gaussian
                mtkern = mtkern / np.sum(mtkern)  # normalize the profile
            else:
                mtkern = np.zeros(2 * nmk + 1)  # init with 0d0
                mtkern[nmk] = 1.0  # delta function

            # Sum the radial and tangential components, weighted by surface area.
            area_r = 0.5  # assume equal areas
            area_t = 0.5  # ar+at must equal 1
            mkern = area_r * mrkern + area_t * mtkern  # add both components

            # Convolve the total flux profiles, again padding the spectrum on both ends to
            # protect against Fourier ringing.
            yfine = convolve(yfine, mkern, mode="nearest")

            # Add contribution from current annulus to the running total.
            if tdnlteH:
                flux = flux + wt[imu] * yfine * mu[imu] # add profile to running total
            else:
                flux = flux + wt[imu] * yfine  # add profile to running total

        flux = np.reshape(flux, (npts, os))  # convert to an array
        if tdnlteH:
            flux = np.pi * np.sum(flux, axis=1) / os  # sum, normalize
        else:
            flux = np.pi * np.sum(flux, axis=1) / os  # sum, normalize

        return flux

    def get_dll_id(self, dll=None):
        if dll is None:
            dll = self.dll
        if dll in __DLL_IDS__:
            dll_id = __DLL_IDS__[dll]
        elif dll in __DLL_DICT__:
            dll_id = dll
        else:
            dll_id = uuid.uuid4()
            __DLL_DICT__[dll_id] = dll
            __DLL_IDS__[dll] = dll_id
        return dll_id

    def get_dll(self, dll_id=None):
        if dll_id is None:
            dll_id = self.dll
        if dll_id in __DLL_DICT__:
            return __DLL_DICT__[dll_id]
        else:
            return dll_id
    
    # @profile
    def synthesize_spectrum(
        self,
        sme,
        segments="all",
        passLineList=True,
        passAtmosphere=True,
        passNLTE=True,
        updateStructure=True,
        updateLineList=False,
        reuse_wavelength_grid=False,
        radial_velocity_mode="robust",
        dll_id=None,
        linelist_mode='all',
        get_opacity=False,
        cdr_database=None,
        cdr_create=False,
        keep_line_opacity=False,
        vbroad_expend_ratio=2,
        contribution_function=False,
        smelib_lineinfo_mode=0,
    ):
        """
        Calculate the synthetic spectrum based on the parameters passed in the SME structure
        The wavelength range of each segment is set in sme.wran
        The specific wavelength grid is given by sme.wave, or is generated on the fly if sme.wave is None

        Will try to fit radial velocity RV and continuum to observed spectrum, depending on vrad_flag and cscale_flag

        Other important fields:
        sme.iptype: instrument broadening type

        Parameters
        ----------
        sme : SME_Struct
            sme structure, with all necessary parameters for the calculation
        setLineList : bool, optional
            wether to pass the linelist to the c library (default: True)
        passAtmosphere : bool, optional
            wether to pass the atmosphere to the c library (default: True)
        passNLTE : bool, optional
            wether to pass NLTE departure coefficients to the c library (default: True)
        reuse_wavelength_grid : bool, optional
            wether to use sme.wint as the output grid of the function or create a new grid (default: False)

        Returns
        -------
        sme : SME_Struct
            same sme structure with synthetic spectrum in sme.smod
        """

        # Prepare 3D NLTE H profile corrections
        if sme.tdnlte_H:
        #     sme.tdnlte_H_correction = self.get_H_3dnlte_correction(sme)
        # if sme.tdnlte_H_new:
            sme.tdnlte_H_correction = self.get_H_3dnlte_correction_rbf(sme)

        if sme is not self.known_sme:
            logger.debug("Synthesize spectrum")
            logger.debug("%s", sme)
            self.known_sme = sme

        # Define constants
        n_segments = sme.nseg
        cscale_degree = sme.cscale_degree

        # fix impossible input
        if "spec" not in sme or sme.spec is None:
            sme.vrad_flag = "none"
            sme.cscale_flag = "none"
        else:
            if "uncs" not in sme or sme.uncs is None:
                sme.uncs = np.ones(sme.spec.size)
            if "mask" not in sme or sme.mask is None:
                sme.mask = np.full(sme.spec.size, MASK_VALUES.LINE)
            for i in range(sme.nseg):
                mask = ~np.isfinite(sme.spec[i])
                mask |= sme.uncs[i] == 0
                sme.mask[i][mask] = MASK_VALUES.BAD

        if radial_velocity_mode != "robust" and ("cscale" not in sme or "vrad" not in sme):
            radial_velocity_mode = "robust"

        segments = self.check_segments(sme, segments)

        # Prepare arrays
        vrad, _, cscale, _ = null_result(sme.nseg, sme.cscale_degree, sme.cscale_type)

        wave = [np.zeros(0) for _ in range(n_segments)]
        smod = [[] for _ in range(n_segments)]
        cmod = [[] for _ in range(n_segments)]
        wmod = [[] for _ in range(n_segments)]
        central_depth = [[] for _ in range(n_segments)]
        line_range = [[] for _ in range(n_segments)]
        opacity = [[] for _ in range(n_segments)]
        if contribution_function:
            sme.contribution_function = [[] for _ in range(n_segments)]
        if updateStructure:
            sme.linelist._lines['nlte_flag'] = -1

        # If wavelengths are already defined use those as output
        if "wave" in sme:
            wave = [w for w in sme.wave]

        dll = self.get_dll(dll_id)
        dll.SetLineInfoMode(int(smelib_lineinfo_mode))

        # "dynamic" is the preferred name in publications; keep "auto" as compatibility alias.
        if linelist_mode == "auto":
            warnings.warn(
                "'linelist_mode=\"auto\"' is deprecated; use 'linelist_mode=\"dynamic\"' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            linelist_mode = "dynamic"
        elif linelist_mode not in ("all", "dynamic"):
            raise ValueError("linelist_mode must be one of: 'all', 'dynamic'")

        # Calculate the line central depth and line range if necessary
        if linelist_mode == 'dynamic':
            # logger.info(f'linelist mode: {linelist_mode}')
            if sme.linelist.cdr_paras is None or not {'central_depth', 'line_range_s', 'line_range_e'}.issubset(sme.linelist._lines.columns) or np.abs(sme.linelist.cdr_paras[0]-sme.teff) >= sme.linelist.cdr_paras_thres['teff'] or (np.abs(sme.linelist.cdr_paras[1]-sme.logg) >= sme.linelist.cdr_paras_thres['logg']) or (np.abs(sme.linelist.cdr_paras[2]-sme.monh) >= sme.linelist.cdr_paras_thres['monh']) or cdr_create:
                logger.info(f'Updating linelist central depth and line range.')
                sme = self.update_cdr(sme, cdr_database=cdr_database, cdr_create=cdr_create, show_progress_bars=show_progress_bars)

        # Input Model data to C library
        dll.SetLibraryPath()
        if passLineList:
            linelist_for_smelib = sme.linelist
            if linelist_mode == 'dynamic':
                line_indices = sme.linelist['wlcent'] < 0
                v_broad = np.sqrt(sme.vmac**2 + sme.vsini**2 + (clight/sme.ipres)**2)
                for i in range(sme.nseg):
                    line_indices |= (sme.linelist['line_range_e'] > sme.wran[i][0] * (1 - vbroad_expend_ratio*v_broad/clight)) & (sme.linelist['line_range_s'] < sme.wran[i][1] * (1 + vbroad_expend_ratio*v_broad/clight))
                line_indices &= self.flag_strong_lines_by_bins(
                    sme.linelist['wlcent'],
                    sme.linelist['central_depth'],
                    threshold=sme.cdr_depth_thres,
                )
                sme.linelist._lines['use_indices'] = line_indices
                linelist_for_smelib = sme.linelist[line_indices]
            line_ion_mask = dll.InputLineList(linelist_for_smelib)
            sme.line_ion_mask = line_ion_mask

            if smelib_lineinfo_mode in (1, 2):
                cols = set(linelist_for_smelib._lines.columns)
                required_cols = {"line_range_s", "line_range_e"}
                missing_cols = sorted(required_cols - cols)

                # Build strong mask from existing column or central depth.
                if "strong" in cols:
                    strong_all = np.asarray(linelist_for_smelib["strong"], dtype=np.uint8)
                elif "central_depth" in cols:
                    strong_all = np.asarray(
                        self.flag_strong_lines_by_bins(
                            linelist_for_smelib["wlcent"],
                            linelist_for_smelib["central_depth"],
                            threshold=sme.cdr_depth_thres,
                        ),
                        dtype=np.uint8,
                    )
                else:
                    strong_all = None

                if strong_all is None:
                    missing_cols.append("strong/central_depth")

                if missing_cols:
                    msg = (
                        f"SMElib lineinfo mode={smelib_lineinfo_mode} requested but missing "
                        f"precomputed columns: {', '.join(missing_cols)}"
                    )
                    if smelib_lineinfo_mode == 2:
                        raise ValueError(msg)
                    logger.warning("%s. Falling back to mode 0.", msg)
                    dll.SetLineInfoMode(0)
                else:
                    keep_mask = ~np.asarray(line_ion_mask, dtype=bool)
                    range_s = np.asarray(linelist_for_smelib["line_range_s"], dtype=np.float64)[keep_mask]
                    range_e = np.asarray(linelist_for_smelib["line_range_e"], dtype=np.float64)[keep_mask]
                    strong_mask = np.asarray(strong_all, dtype=np.uint8)[keep_mask]

                    if "central_depth" in cols:
                        depth = np.asarray(linelist_for_smelib["central_depth"], dtype=np.float64)[keep_mask]
                        dll.InputLinePrecomputedInfo(range_s, range_e, strong_mask, depth)
                    else:
                        dll.InputLinePrecomputedInfo(range_s, range_e, strong_mask)
        if hasattr(updateLineList, "__len__") and len(updateLineList) > 0:
            # TODO Currently Updates the whole linelist, could be improved to only change affected lines
            dll.UpdateLineList(sme.atomic, sme.species, updateLineList)
        if passAtmosphere:
            sme = self.get_atmosphere(sme)
            dll.InputModel(sme.teff, sme.logg, sme.vmic, sme.atmo)
            dll.InputAbund(sme.abund)
            dll.Ionization(0)
            dll.SetVWscale(sme.gam6)
            dll.SetH2broad(sme.h2broad)
        if passNLTE:
            sme.nlte.update_coefficients(sme, dll, self.lfs_nlte)

        # Loop over segments
        #   Input Wavelength range and Opacity
        #   Calculate spectral synthesis for each
        #   Interpolate onto geomspaced wavelength grid
        #   Apply instrumental and turbulence broadening
        sme.first_segment = True
        for il in tqdm(segments, desc="Segments", leave=True, disable=not show_progress_bars):
            wmod[il], smod[il], cmod[il], central_depth[il], line_range[il], opacity[il] = self.synthesize_segment(
                sme,
                il,
                reuse_wavelength_grid,
                dll_id=dll_id,
                get_opacity=get_opacity,
                keep_line_opacity=keep_line_opacity,
                contribution_function=contribution_function
            )
        for il in segments:
            if "wave" not in sme or len(sme.wave[il]) == 0:
                # trim padding
                wbeg, wend = sme.wran[il]
                itrim = (wmod[il] > wbeg) & (wmod[il] < wend)
                # Force endpoints == wavelength range
                wave[il] = np.concatenate(([wbeg], wmod[il][itrim], [wend]))

        if sme.specific_intensities_only:
            return wmod, smod, cmod, central_depth

        # For testing wavegrid
        sme.wmod = wmod.copy()
        sme.smod = smod.copy()
        sme.comd = cmod.copy()
        sme.opacity = opacity.copy()

        # Fit continuum and radial velocity
        # And interpolate the flux onto the wavelength grid
        if radial_velocity_mode == "robust":
            cscale, cscale_unc, vrad, vrad_unc = match_rv_continuum(
                sme, segments, wmod, smod
            )
            logger.debug("Radial velocity: %s", str(vrad))
            logger.debug("Continuum coefficients: %s", str(cscale))
        elif radial_velocity_mode == "fast":
            cscale, vrad = sme.cscale, sme.vrad
        else:
            raise ValueError("Radial Velocity mode not understood")

        smod, cmod = self.apply_radial_velocity_and_continuum_synth(
            wave,
            sme.spec,
            wmod,
            smod,
            cmod,
            vrad,
            cscale,
            sme.cscale_type,
            segments,
        )

        # Merge all segments
        # if sme already has a wavelength this should be the same
        if updateStructure:
            if "wave" not in sme:
                # TODO: what if not all segments are there?
                sme.wave = wave
            if "synth" not in sme:
                sme.synth = smod
            if "cont" not in sme:
                sme.cont = cmod

            for s in segments:
                sme.wave[s] = wave[s]
                sme.synth[s] = smod[s]
                sme.cont[s] = cmod[s]
                # sme.central_depth[s] = central_depth[s]
                # sme.line_range[s] = line_range[s]


            if passLineList and (self.update_cdr_switch or linelist_mode == 'all'):
                s = 0
                if len(central_depth[s]) > 0:
                    sme.linelist._lines.loc[~sme.line_ion_mask, 'central_depth'] = central_depth[s]
                    sme.linelist._lines.loc[sme.line_ion_mask, 'central_depth'] = np.nan
                    sme.linelist._lines.loc[~sme.line_ion_mask, 'line_range_s'] = line_range[s][:, 0]
                    sme.linelist._lines.loc[~sme.line_ion_mask, 'line_range_e'] = line_range[s][:, 1]
                    sme.linelist._lines.loc[sme.line_ion_mask, 'line_range_s'] = np.nan
                    sme.linelist._lines.loc[sme.line_ion_mask, 'line_range_e'] = np.nan
                    sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])
                    # Manually change the 2000 line_range to 0.03.
                    indices = np.isclose(sme.linelist['line_range_e'] - sme.linelist['line_range_s'], 2000, rtol=1e-4, atol=5, equal_nan=False)
                    sme.linelist._lines.loc[indices, 'line_range_s'] = sme.linelist._lines.loc[indices, 'wlcent']-0.3
                    sme.linelist._lines.loc[indices, 'line_range_e'] = sme.linelist._lines.loc[indices, 'wlcent']+0.3

            if sme.cscale_type in ["spline", "spline+mask"]:
                sme.cscale = np.asarray(cscale)
                sme.cscale_unc = np.asarray(cscale_unc)
            elif sme.cscale_flag not in ["fix", "none"]:
                for s in np.arange(sme.nseg):
                    if s not in segments:
                        cscale[s] = sme.cscale[s]
                sme.cscale = np.asarray(cscale)
                sme.cscale_unc = np.asarray(cscale_unc)

            sme.vrad = np.asarray(vrad)
            sme.vrad_unc = np.asarray(vrad_unc)
            nlte_flags = dll.GetNLTEflags()
            if linelist_mode == 'dynamic':
                sme.linelist._lines.loc[sme.linelist._lines['use_indices'], 'nlte_flag'] = nlte_flags.astype(int)
            else:
                sme.nlte.flags = nlte_flags
                sme.linelist._lines.loc[~sme.line_ion_mask, 'nlte_flag'] = nlte_flags.astype(int)

        # Store the adaptive wavelength grid for the future

            result = sme
        else:
            wave = Iliffe_vector(values=wave)
            smod = Iliffe_vector(values=smod)
            cmod = Iliffe_vector(values=cmod)
            result = wave, smod, cmod

        # Cleanup
        return result

    # @profile
    def synthesize_segment(
        self,
        sme,
        segment,
        reuse_wavelength_grid=False,
        keep_line_opacity=False,
        dll_id=None,
        get_opacity=False,
        contribution_function=False
    ):
        """Create the synthetic spectrum of a single segment

        Parameters
        ----------
        sme : SME_Struct
            The SME strcuture containing all relevant parameters
        segment : int
            the segment to synthesize
        setLineList : bool, optional
            wether to pass the linelist to the c library (default: True)
        reuse_wavelength_grid : bool
            Whether to keep the current wavelength grid for the synthesis
            or create a new one, depending on the linelist. Default: False
        keep_line_opacity : bool
            Whether to reuse existing line opacities or not. This should be
            True if the opacities have been calculated in another segment.

        Returns
        -------
        wgrid : array of shape (npoints,)
            Wavelength grid of the synthesized spectrum
        flux : array of shape (npoints,)
            The Flux of the synthesized spectrum
        cont_flux : array of shape (npoints,)
            The continuum Flux of the synthesized spectrum
        """        
        logger.debug("Segment %i out of %i", segment, sme.nseg)
        dll = self.get_dll(dll_id)

        # Input Wavelength range and Opacity
        vrad_seg = sme.vrad[segment] if sme.vrad[segment] is not None else 0
        wbeg, wend = self.get_wavelengthrange(sme.wran[segment], vrad_seg, sme.vsini)

        # if passLineList:
            # if linelist_mode == 'all':
            #     line_ion_mask = dll.InputLineList(sme.linelist)
            #     sme.line_ion_mask = line_ion_mask
            # elif linelist_mode == 'auto':
            #     # Check if the current stellar parameters are within the range of that in the linelist
            #     if np.abs(sme.linelist.cdr_paras[0]-sme.teff) >= 500 or (np.abs(sme.linelist.cdr_paras[1]-sme.logg) >= 1) or (np.abs(sme.linelist.cdr_paras[2]-sme.monh) >= 0.5): 
            #         logger.warning(f'The current stellar parameters are out of the range (+-500K, +- 1 or +-0.5) of that used to calculate the central depth and ranges of the linelist. \n Current stellar parameters: Teff: {sme.teff}, logg: {sme.logg}, monh: {sme.monh}. Linelist parameters: Teff: {sme.linelist.cdr_paras[0]}, logg: {sme.linelist.cdr_paras[1]}, monh: {sme.linelist.cdr_paras[2]}.')
                
            #     v_broad = np.sqrt(sme.vmic**2 + sme.vmac**2 + sme.vsini**2)
            #     del_wav = v_broad * sme.linelist['wlcent'] / clight
            #     ipres_segment = sme.ipres if np.size(sme.ipres) == 1 else sme.ipres[segment]
            #     if ipres_segment != 0:
            #         del_wav += sme.linelist['wlcent'] / ipres_segment
            #     indices = (~((sme.linelist['line_range_e'] < wbeg - del_wav - line_margin) | (sme.linelist['line_range_s'] > wend + del_wav + line_margin))) & (sme.linelist['central_depth'] > sme.cdr_depth_thres)
            #     # logger.info(f"There are {len(sme.linelist[indices][np.char.find(sme.linelist[indices]['species'], 'Fe') >= 0])} Fe lines in sub linelist.")
            #     _ = dll.InputLineList(sme.linelist[indices])
            #     sme.linelist._lines['use_indices'] = indices
        # if hasattr(updateLineList, "__len__") and len(updateLineList) > 0:
        #     # TODO Currently Updates the whole linelist, could be improved to only change affected lines
        #     dll.UpdateLineList(sme.atomic, sme.species, updateLineList)

        # if passAtmosphere:
        #     dll.InputModel(sme.teff, sme.logg, sme.vmic, sme.atmo)
        #     dll.InputAbund(sme.abund)
        #     dll.Ionization(0)
        #     dll.SetVWscale(sme.gam6)
        #     dll.SetH2broad(sme.h2broad)

        # if passNLTE:
        #     sme.nlte.update_coefficients(sme, dll, self.lfs_nlte, sme.first_segment)        
        
        dll.InputWaveRange(wbeg-2, wend+2)
        dll.Opacity()

        # Priority for wavelength grid passed to SMElib:
        # 1) user-provided sme.wint for this segment
        # 2) internal cache when reuse_wavelength_grid=True
        # 3) None (let SMElib compute it)
        user_wint_seg = None
        if getattr(sme, "wint", None) is not None:
            try:
                user_wint_seg = sme.wint[segment]
                if user_wint_seg is not None and len(user_wint_seg) == 0:
                    user_wint_seg = None
            except (IndexError, KeyError, TypeError):
                user_wint_seg = None

        if user_wint_seg is not None:
            wint_seg = np.asarray(user_wint_seg, dtype=np.float64)
        elif reuse_wavelength_grid and segment in self.wint:
            wint_seg = self.wint[segment]
        else:
            wint_seg = None

        # Only calculate line opacities in the first segment
        #   Calculate spectral synthesis for each
        _, wint, sint, cint = dll.Transf(
            sme.mu,
            accrt=sme.accrt,  # threshold line opacity / cont opacity
            accwi=sme.accwi,
            keep_lineop=keep_line_opacity and not sme.first_segment,
            wave=wint_seg,
        )

        # Insert the new 3DNLTE correction
        if sme.tdnlte_H:
            interpolator = interp1d(util.lambda_H_3DNLTE, sme.tdnlte_H_correction, kind="linear", fill_value=1, bounds_error=False, assume_sorted=True)
            correction_3dnlte_H_interp = interpolator(wint)

            sint *= correction_3dnlte_H_interp

        # # Assign the nlte flags
        # nlte_flags = dll.GetNLTEflags()
        # sme.nlte.flags = nlte_flags

        # Store the adaptive wavelength grid for the future
        # if it was newly created
        if user_wint_seg is None and wint_seg is None:
            self.wint[segment] = wint

        if not sme.specific_intensities_only:
            # Create new geomspaced wavelength grid, to be used for intermediary steps
            wgrid, vstep = self.new_wavelength_grid(wint)

            logger.debug("Integrate specific intensities")
            # Radiative Transfer Integration
            # Continuum
            cint = self.integrate_flux(sme.mu, cint, 1, 0, 0)
            cint = np.interp(wgrid, wint, cint)

            # Broaden Spectrum
            y_integrated = np.empty((sme.nmu, len(wgrid)))
            for imu in range(sme.nmu):
                y_integrated[imu] = np.interp(wgrid, wint, sint[imu])

            # Turbulence broadening
            # Apply macroturbulent and rotational broadening while integrating intensities
            # over the stellar disk to produce flux spectrum Y.
            sint = self.integrate_flux(sme.mu, y_integrated, vstep, sme.vsini, sme.vmac)
            wint = wgrid

            # instrument broadening
            if "iptype" in sme:
                logger.debug("Apply detector broadening")
                ipres = sme.ipres.item() if np.size(sme.ipres) == 1 else sme.ipres[segment]
                sint = broadening.apply_broadening(ipres, wint, sint, type=sme.iptype, sme=sme)

            # Apply the correction on Ha, Hb and Hgamma line here.
            if sme.tdnlte_H:
                correction_resample = safe_interpolation(sme.tdnlte_H_correction[0], sme.tdnlte_H_correction[1], wint, fill_value=1)
                sint *= correction_resample

        # Divide calculated spectrum by continuum
        if sme.normalize_by_continuum:
            sint /= cint

        # Mingjie: run CentralDepth
        central_depth = dll.CentralDepth(sme.mu, sme.accrt)
        line_range = dll.GetLineRange()
        if get_opacity:
            opacity = []
            for wave_single in sme.wave[segment]:
                opacity.append(dll.GetLineOpacity(wave_single))
        else:
            opacity = None

        if contribution_function:
            cf = dll.GetContributionfunction(sme.mu, sme.wave[segment])[0]
            # Temporary: set the cf in outermost to be the same as the second one in rhox,
            #            to avoid impractical numbers
            cf[..., 0] = cf[..., 1]
            if not sme.specific_intensities_only:
                cf = np.array([self.integrate_flux(sme.mu, cf_mu, 1, 0, 0) for cf_mu in cf])
            sme.contribution_function[segment] = cf

        sme.first_segment = False
        return wint, sint, cint, central_depth, line_range, opacity
    
    def update_cdr(self, sme, cdr_database=None, cdr_create=False, cdr_grid_overwrite=False, mode='linear', dims=['teff', 'logg', 'monh'], show_progress_bars=show_progress_bars):
        '''
        Update or get the central depth and wavelength range of a line list. This version separate the parallel and non-parallel mode completely.
        Author: Mingjie Jian
        '''

        N_line_chunk, parallel, n_jobs, pysme_out = sme.cdr_N_line_chunk, sme.cdr_parallel, sme.cdr_n_jobs, sme.cdr_pysme_out
        self.update_cdr_switch = True

        # Decide how many chunks to be divided
        N_chunk = int(np.ceil(len(sme.linelist) / N_line_chunk))

        # Divide the line list to sub line lists
        sub_linelist = [sme.linelist[N_line_chunk*i:N_line_chunk*(i+1)] for i in range(N_chunk)]

        if sum(len(item) for item in sub_linelist) != len(sme.linelist):
            raise ValueError
        
        if cdr_database is not None:
            # cdr_database is provided, use it to update the central depth and line range
            self._interpolate_or_compute_and_update_linelist(sme, cdr_database, cdr_create=cdr_create, cdr_grid_overwrite=cdr_grid_overwrite, mode=mode, dims=dims, show_progress_bars=show_progress_bars)
            return sme

        logger.info('[cdr] Using calculation to update central depth and line range.')
        sub_sme_init = SME_Structure()
        exclude_keys = ['_wave', '_synth', '_spec', '_uncs', '_mask', '_SME_Structure__wran', '_normalize_by_continuum', '_specific_intensities_only', '_telluric', '__cont', '_linelist', '_fitparameters', '_fitresults']
        for key, value in sme.__dict__.items():
            if key not in exclude_keys and 'cscale' not in key and 'vrad' not in key:
                setattr(sub_sme_init, key, deepcopy(value))
        sub_sme_init.wave = np.arange(5000, 5010, 1)

        if not parallel:
            for i in tqdm(range(N_chunk), disable=not show_progress_bars):
                sub_sme_init.linelist = sub_linelist[i]
                sub_sme_init = self.synthesize_spectrum(sub_sme_init)
                if i == 0:
                    stack_linelist = deepcopy(sub_sme_init.linelist)
                else:
                    stack_linelist._lines = pd.concat([stack_linelist._lines, sub_sme_init.linelist._lines])
        else:
            sub_sme = []
            sub_sme_init.linelist = sme.linelist[:1]
            sub_sme_init = self.synthesize_spectrum(sub_sme_init)
            for i in range(N_chunk):
                sub_sme.append(deepcopy(sub_sme_init))
                sub_sme[i].linelist = sub_linelist[i]

            if pysme_out:
                sub_sme = pqdm(sub_sme, self.synthesize_spectrum, n_jobs=n_jobs, disable=not show_progress_bars)
            else:
                with redirect_stdout(open(f"/dev/null", 'w')):
                    sub_sme = pqdm(sub_sme, self.synthesize_spectrum, n_jobs=n_jobs, disable=not show_progress_bars)
            
            for i in range(N_chunk):
                sub_linelist[i] = sub_sme[i].linelist
            stack_linelist = deepcopy(sub_linelist[0])
            stack_linelist._lines = pd.concat([ele._lines for ele in sub_linelist])  
            # logger.info(f'{sub_linelist}')

        # Remove
        if len(stack_linelist) != len(sme.linelist):
            raise ValueError
        for column in ['central_depth', 'line_range_s', 'line_range_e']:
            if column in sme.linelist.columns:
                sme.linelist._lines = sme.linelist._lines.drop(column, axis=1)
        for column in ['central_depth', 'line_range_s', 'line_range_e']:
            sme.linelist._lines[column] = stack_linelist._lines[column]
        # pickle.dump([sme.linelist._lines, stack_linelist._lines], open('linelist.pkl', 'wb'))

        # Manually change the depth of all H 1 lines to 1, to include them back.
        sme.linelist._lines.loc[sme.linelist['species'] == 'H 1', 'central_depth'] = 1 

        # Manually change the 2000 line_range to 0.03.
        indices = np.isclose(sme.linelist['line_range_e'] - sme.linelist['line_range_s'], 2000, rtol=1e-4, atol=5, equal_nan=False)
        sme.linelist._lines.loc[indices, 'line_range_s'] = sme.linelist._lines.loc[indices, 'wlcent']-0.3
        sme.linelist._lines.loc[indices, 'line_range_e'] = sme.linelist._lines.loc[indices, 'wlcent']+0.3

        # Write the stellar parameters used here to the line list
        sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])

        self.update_cdr_switch = False

        return sme

    def _interpolate_or_compute_and_update_linelist(
        self, sme, cdr_database, cdepth_decimals=4, cdepth_thres=0,
        range_decimals=2, cdr_create=False, cdr_grid_overwrite=False,
        mode='linear', dims=['teff', 'logg', 'monh'], show_progress_bars=False
    ):
        teff, logg, monh, vmic = sme.teff, sme.logg, sme.monh, sme.vmic
        param = np.array([teff, logg, monh, vmic])

        param_grid, fname_map = _load_all_grid_points(cdr_database)

        if len(dims) == 3 and len(param_grid) > 0:
            # Remove vmic in the grid
            # vmic_mask = np.isclose(param_grid[:, -1], fixed_vmic)
            param_grid = param_grid[:, :-1]

            fname_map = {
                k[:-1]: v for k, v in fname_map.items()
            }
            param = param[:-1]

        if len(param_grid) > 0:

            thres = sme.linelist.cdr_paras_thres
            
            lower = np.array([teff - thres['teff'], logg - thres['logg'], monh - thres['monh'], vmic - thres['vmic']])
            upper = np.array([teff + thres['teff'], logg + thres['logg'], monh + thres['monh'], vmic + thres['vmic']])

            if len(dims) == 3:
                lower = lower[:-1]
                upper = upper[:-1]

            in_box_mask = np.all((param_grid >= lower) & (param_grid <= upper), axis=1)
            filtered_grid = param_grid[in_box_mask]

            if mode == 'linear' and len(filtered_grid) >= len(dims)+1 and not cdr_create:
                delaunay = Delaunay(filtered_grid)
                simplex_index = delaunay.find_simplex(param)
                if simplex_index >= 0:
                    logger.info('[cdr] Using linear interpolation from grid database.')
                    vertex_indices = delaunay.simplices[simplex_index]
                    vertices = filtered_grid[vertex_indices]

                    n_lines_total = len(sme.linelist)
                    interpolated_arrays = {
                        'central_depth': np.zeros(n_lines_total, dtype=np.float32),
                        'line_range_s':  sme.linelist['wlcent'] - 0.3,
                        'line_range_e':  sme.linelist['wlcent'] + 0.3,
                    }
                    vertex_arrays = {k: [] for k in interpolated_arrays}
                    for pt in vertices:
                        logger.info(f'Using {fname_map[tuple(pt)]}')
                        data = np.load(os.path.join(cdr_database, fname_map[tuple(pt)]))['line_info']
                        iloc = data[:, 0].astype(int)

                        arr_cdepth = np.zeros(n_lines_total, dtype=np.float32)
                        # arr_lrs    = np.full(n_lines_total, np.nan, dtype=np.float32)
                        # arr_lre    = np.full(n_lines_total, np.nan, dtype=np.float32)
                        arr_lrs    = sme.linelist['wlcent'] - 0.3
                        arr_lre    = sme.linelist['wlcent'] + 0.3
                        arr_cdepth[iloc] = data[:, 1]
                        arr_lrs[iloc]    = data[:, 2]
                        arr_lre[iloc]    = data[:, 3]

                        vertex_arrays['central_depth'].append(np.log10(arr_cdepth))
                        vertex_arrays['line_range_s'].append(arr_lrs)
                        vertex_arrays['line_range_e'].append(arr_lre)

                    # Do a rough normalization
                    vertices[:, 0] /= 1000
                    param[0] /= 1000
                    for key in interpolated_arrays:
                        stacked = np.stack(vertex_arrays[key], axis=0)
                        interp  = LinearNDInterpolator(vertices, stacked)
                        interpolated_arrays[key] = interp(param)[0]

                    sme.linelist._lines['central_depth'] = 10 ** interpolated_arrays['central_depth']
                    sme.linelist._lines['line_range_s']  = interpolated_arrays['line_range_s']
                    sme.linelist._lines['line_range_e']  = interpolated_arrays['line_range_e']
                    if len(dims) == 4:
                        sme.linelist.cdr_paras = np.array([teff, logg, monh, vmic])
                    elif len(dims) == 3:
                        sme.linelist.cdr_paras = np.array([teff, logg, monh])
                    return

            elif mode == 'nearest' and len(filtered_grid) > 0 and not cdr_create:
                # box 内找到最近的 grid 点
                dists = cdist(filtered_grid, param[None, :])
                i_nearest = np.argmin(dists)
                nearest_pt = filtered_grid[i_nearest]
                logger.info(f'[CDF] Using nearest grid point {nearest_pt} in box for direct assignment.')

                data = np.load(os.path.join(cdr_database, fname_map[tuple(nearest_pt)]))['line_info']
                iloc = data[:, 0].astype(int)

                n_lines_total = len(sme.linelist)
                arr_cdepth = np.zeros(n_lines_total, dtype=np.float32)
                arr_lrs    = np.full(n_lines_total, np.nan, dtype=np.float32)
                arr_lre    = np.full(n_lines_total, np.nan, dtype=np.float32)

                arr_cdepth[iloc] = data[:, 1]
                arr_lrs[iloc]    = data[:, 2]
                arr_lre[iloc]    = data[:, 3]

                sme.linelist._lines['central_depth'] = arr_cdepth
                sme.linelist._lines['line_range_s']  = arr_lrs
                sme.linelist._lines['line_range_e']  = arr_lre
                sme.linelist.cdr_paras = nearest_pt
                    
                return

        # ---------- fallback: compute and save ----------
        fname = f"teff{teff:.0f}_logg{logg:.1f}_monh{monh:.1f}_vmic{vmic:.1f}.npz"
        full_path = os.path.join(cdr_database, fname)
        if os.path.exists(full_path) and not cdr_grid_overwrite:
            logger.info(f"{fname} exists and cdr_grid_overwrite is false, skipping generating cdr grid.")
        else:
            logger.info("[CDF] Fallback: recomputing line properties.")
            self.update_cdr(sme, cdr_database=None, show_progress_bars=show_progress_bars)

            mask = sme.linelist['central_depth'] > cdepth_thres
            filtered_df = sme.linelist[['central_depth', 'line_range_s', 'line_range_e']][mask]
            filtered_iloc = np.where(mask)[0]

            depth = np.round(filtered_df[:, 0], cdepth_decimals)
            range_s = np.round(filtered_df[:, 1], range_decimals)
            range_e = np.round(filtered_df[:, 2], range_decimals)

            line_info = np.column_stack([filtered_iloc, depth, range_s, range_e])
            logger.info(f'Saving {fname} to database.')
            np.savez_compressed(full_path, line_info=line_info)

    def keep_mask_by_cumulative_depth(self, depths: np.ndarray, threshold=0.01):
        """
        Vectorised version for a single wavelength bin.
        """
        order = depths.argsort()                    # weakest → strongest
        sorted_depths = depths[order]
        cumsum = sorted_depths.cumsum()
        cut = np.searchsorted(cumsum, threshold, side="right")
        drop_sorted = np.zeros_like(sorted_depths, dtype=bool)
        drop_sorted[:cut] = True                    # drop the weakest until threshold reached
        drop_mask = np.zeros_like(depths, dtype=bool)
        drop_mask[order] = drop_sorted              # map back to input order
        return ~drop_mask                           # True = keep

    @staticmethod
    def flag_strong_lines_by_bins(wl, depth, bin_width=0.2, threshold=0.001):
        wl = np.asarray(wl)
        depth = np.asarray(depth)

        # 1) 计算 bin_idx
        wl_min, wl_max = wl.min(), wl.max()
        edges = np.arange(wl_min, wl_max + bin_width, bin_width)
        bin_idx = np.searchsorted(edges, wl, side="right") - 1
        # 假设 bin_idx ∈ [0, n_bins-1]；若有越界需先 clip

        # 2) 一次性分组排序：先按 bin，再按 depth（升序：弱→强）
        order = np.lexsort((depth, bin_idx))     # 主键：bin_idx；次键：depth
        bin_sorted = bin_idx[order]
        depth_sorted = depth[order]

        # 3) 找到每个 bin 的切片边界
        starts = np.r_[0, np.flatnonzero(np.diff(bin_sorted)) + 1]
        ends   = np.r_[starts[1:], len(depth_sorted)]

        # 4) 在“按强度升序”的每个切片内，丢弃前 cut 个（累计和<=threshold）
        keep_sorted = np.empty_like(depth_sorted, dtype=bool)
        csum = np.cumsum(depth_sorted)  # 全局前缀和，便于O(1)求任一切片的局部cumsum

        for s, e in zip(starts, ends):
            if s == e:  # 空bin（通常不会被lexsort产生，但留作保护）
                continue
            base = csum[s-1] if s > 0 else 0.0
            # 该bin内局部累计和
            local = csum[s:e] - base
            # 需要丢弃的“最弱”个数
            cut = np.searchsorted(local, threshold, side="right")
            # 弱的前 cut 个 False，其余 True
            keep_sorted[s:e] = True
            if cut > 0:
                keep_sorted[s:s+cut] = False

        # 5) 映射回原顺序
        keep = np.empty_like(keep_sorted)
        keep[order] = keep_sorted
        return keep

    def flag_strong_lines_by_bins_old(self, df, bin_width=0.2, threshold=0.01, wl_col="wlcent", depth_col="central_depth", out_col="keep_mask", show_progress_bars=show_progress_bars):
        """
        Add a boolean 'keep_mask' column to the VALD line list indicating strong / weak lines.

        Parameters
        ----------
        df : pd.DataFrame
            VALD line list containing wavelength and depth columns.
        bin_width : float, default 0.2
            Width of each wavelength bin (same unit as wl_col, typically Å).
        threshold : float, default 0.01
            Maximum cumulative depth to discard within each bin.
        wl_col : str, default "wlcent"
            Column name for line‐center wavelength.
        depth_col : str, default "central_depth"
            Column name for line‐center depth.
        out_col : str, default "keep_mask"
            Name of the boolean column to write into `df`.

        Returns
        -------
        pandas.Series (dtype=bool)
            Boolean mask aligned to `df.index`; True → strong line.
        """
        # ----- 0.  Prepare numpy views -----
        wl = df[wl_col].to_numpy()
        depth = df[depth_col].to_numpy()
        n_lines = len(df)

        # ----- 1.  Build wavelength bin edges -----
        wl_min, wl_max = wl.min(), wl.max()
        edges = np.arange(wl_min, wl_max + bin_width, bin_width)  # rightmost edge inclusive
        n_bins = len(edges) - 1

        # Bin index for each line: 0 ... n_bins-1
        bin_idx = np.searchsorted(edges, wl, side="right") - 1

        # ----- 2.  Allocate output mask -----
        keep_mask = np.zeros(n_lines, dtype=bool)

        # ----- 3.  Loop over bins, apply cumulative-depth filter -----
        for b in tqdm(range(n_bins), disable=not show_progress_bars):
            idx = np.where(bin_idx == b)[0]         # indices of lines in this bin
            if idx.size == 0:                       # empty bin -> skip
                continue
            keep_mask[idx] = self.keep_mask_by_cumulative_depth(
                depth[idx],
                threshold=threshold
            )

        # ----- 4.  Save into DataFrame & return -----
        df[out_col] = keep_mask
        return df[out_col]
    def flag_strong_lines_by_database(
        self, sme, cdr_negligibe_database, dims=['teff', 'logg', 'monh'], mode='or'
    ):
        """
        Flag strong lines for the current star using the compressed CDR strong-line database,
        and write both the 'strong' flag and the corresponding line_range_s / line_range_e
        into the SME linelist.

        The compressed CDR database is assumed to contain one npz file per grid point in
        `cdr_negligibe_database`. Each file should include:
            - mask_bits     : 1D uint8 array, bit-packed boolean mask for strong lines
            - unique_widths : 1D float32 array, distinct line_width values in this grid
            - codes         : 1D uint8/uint16/uint32 array indexing unique_widths,
                              ordered consistently with the strong-line indices.

        Interpolation modes
        -------------------
        mode : {'or', 'nearest'}
            'or':
                Within the parameter box around the current star, construct a Delaunay
                triangulation in the selected parameter dimensions (dims). If the
                current parameter lies inside a simplex, take all vertices of that simplex:
                  - Combine their strong masks using a logical OR across vertices.
                  - For line_range_s/e, take the union of ranges: s = nanmin(s_i),
                    e = nanmax(e_i) across vertices (ignoring NaNs).
            'nearest':
                Within the parameter box, find the single nearest grid point and directly
                use its strong mask and line_range_s/e.

        Parameters
        ----------
        sme : SME_Structure
            The current SME object (stellar parameters and linelist are taken from here).
        cdr_negligibe_database : str
            Directory containing the compressed CDR grid files.
        dims : list of str, default ['teff', 'logg', 'monh']
            Parameter dimensions used for interpolation. If length is 3, vmic is ignored.
        mode : {'or', 'nearest'}, default 'or'
            Interpolation / combination strategy as described above.
        """
        teff, logg, monh, vmic = sme.teff, sme.logg, sme.monh, sme.vmic
        param = np.array([teff, logg, monh, vmic], dtype=float)

        # Load all grid points and map (teff, logg, monh, vmic) -> filename
        param_grid, fname_map = _load_all_grid_points(cdr_negligibe_database)

        # In 3D mode, drop vmic dimension from grid and from the parameter vector
        if len(dims) == 3 and len(param_grid) > 0:
            param_grid = param_grid[:, :-1]
            fname_map = {k[:-1]: v for k, v in fname_map.items()}
            param = param[:-1]

        if len(param_grid) == 0:
            logger.warning("[cdr] Compressed strong-line database is empty; "
                           "flag_strong_lines_by_database will be skipped.")
            return

        # ----- Select grid points within the parameter box around the current star -----
        thres = sme.linelist.cdr_paras_thres
        lower = np.array(
            [teff - thres['teff'], logg - thres['logg'], monh - thres['monh'], vmic - thres['vmic']],
            dtype=float
        )
        upper = np.array(
            [teff + thres['teff'], logg + thres['logg'], monh + thres['monh'], vmic + thres['vmic']],
            dtype=float
        )

        if len(dims) == 3:
            lower = lower[:-1]
            upper = upper[:-1]

        in_box_mask = np.all((param_grid >= lower) & (param_grid <= upper), axis=1)
        filtered_grid = param_grid[in_box_mask]

        if filtered_grid.size == 0:
            logger.warning("[cdr] No grid points found within the parameter box; "
                           "cannot flag strong lines from database.")
            return

        n_lines_total = len(sme.linelist)
        wlcent = sme.linelist['wlcent']

        # ----- mode = 'or': use simplex vertices and OR-combine masks -----
        if mode == 'or' and filtered_grid.shape[0] >= len(dims) + 1:
            from scipy.spatial import Delaunay

            delaunay = Delaunay(filtered_grid)
            simplex_index = delaunay.find_simplex(param)

            if simplex_index >= 0:
                logger.info("[cdr] Combining strong-line masks and line ranges from "
                            "compressed database (simplex OR).")
                vertex_indices = delaunay.simplices[simplex_index]
                vertices = filtered_grid[vertex_indices]

                masks = []
                lrs_list = []
                lre_list = []

                for pt in vertices:
                    fname = fname_map[tuple(pt)]
                    full_path = os.path.join(cdr_negligibe_database, fname)
                    logger.info(f"[cdr] Using compressed strong-line file {fname} for OR-combination.")
                    mask_v, lrs_v, lre_v = _load_mask_and_ranges_from_compressed(
                        full_path, wlcent, n_lines_total
                    )
                    masks.append(mask_v)
                    lrs_list.append(lrs_v)
                    lre_list.append(lre_v)

                masks_arr = np.vstack(masks)               # (n_vertex, N_lines)
                strong_final = np.any(masks_arr, axis=0)   # logical OR across vertices

                # For line_range_s/e, take the union across vertices:
                # s = min over vertices (ignoring NaN), e = max over vertices.
                lrs_arr = np.vstack(lrs_list)
                lre_arr = np.vstack(lre_list)

                line_range_s_final = np.nanmin(lrs_arr, axis=0)
                line_range_e_final = np.nanmax(lre_arr, axis=0)

                sme.linelist._lines['strong']       = strong_final
                sme.linelist._lines['line_range_s'] = line_range_s_final
                sme.linelist._lines['line_range_e'] = line_range_e_final
                return

            # If the parameter is not inside any simplex, fall back to nearest strategy
            logger.warning("[cdr] Target parameter is not inside any simplex; "
                           "falling back to 'nearest' strategy.")
            mode = 'nearest'

        # ----- mode = 'nearest': use the nearest grid point within the box -----
        if mode == 'nearest':
            dists = cdist(filtered_grid, param[None, :])
            i_nearest = np.argmin(dists)
            nearest_pt = filtered_grid[i_nearest]
            logger.info(f"[cdr] Using nearest grid point {nearest_pt} in box for "
                        "strong-line mask and line ranges.")

            fname = fname_map[tuple(nearest_pt)]
            full_path = os.path.join(cdr_negligibe_database, fname)
            strong_final, line_range_s_final, line_range_e_final = _load_mask_and_ranges_from_compressed(
                full_path, wlcent, n_lines_total
            )

            sme.linelist._lines['strong']       = strong_final
            sme.linelist._lines['line_range_s'] = line_range_s_final
            sme.linelist._lines['line_range_e'] = line_range_e_final
            return

        logger.warning("[cdr] flag_strong_lines_by_database: no valid interpolation/nearest "
                       "strategy found; 'strong' and line_range_s/e were not set.")
        return

    def get_H_3dnlte_correction_rbf(self, sme):
        """
        Compute the 3D NLTE correction factor for hydrogen lines, using RBF interpolator and in intensities.
    
        """

        logger.info(f"Getting H 3dnlte correction using RBF")

        sme_H_only = SME_Structure()
        sme_H_only.teff, sme_H_only.logg, sme_H_only.monh, sme_H_only.vmic, sme_H_only.vmac, sme_H_only.vsini = sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini
        sme_H_only.iptype = sme.iptype
        sme_H_only.ipres = sme.ipres
        sme_H_only.specific_intensities_only = True
        # sme_H_only.normalize_by_continuum = False
        for i in range(len(sme.nlte.elements)):
            sme_H_only.nlte.set_nlte(sme.nlte.elements[i], sme.nlte.grids[sme.nlte.elements[i]])
        sme_H_only.linelist = sme.linelist[(sme.linelist['species'] == 'H 1')]
        sme_H_only.wave = np.arange(4000, 6700, 0.02)
        sme_H_only.tdnlte_H = False
        sme_H_only_res = self.synthesize_spectrum(sme_H_only)

        int_3dnlte_H = []
        interpolator = interp1d(sme_H_only_res[0][0], sme_H_only_res[1][0], kind="linear", fill_value=1, bounds_error=False, assume_sorted=True)
        int_1d_H = interpolator(util.lambda_H_3DNLTE)
        for mu in sme_H_only.mu:
            int_3dnlte_H_mu, in_boundary = interpolate_3DNLTEH_spectrum_RBF(sme_H_only.teff, sme_H_only.logg, sme_H_only.monh, mu, boundary_vertices)
            int_3dnlte_H.append(int_3dnlte_H_mu)
        
        if not in_boundary:
            logger.info(f"Outside the H 3dnlte grid, not performing correction.")
            sme.tdnlte_H = False
            return None
        
        int_3dnlte_H = np.array(int_3dnlte_H)
        int_1d_H = np.array(int_1d_H)
        correction = int_3dnlte_H / int_1d_H

        return correction

    # def get_H_3dnlte_correction(self, sme):
    #     """
    #     Compute the 3D NLTE correction factor for hydrogen lines.

    #     This function:
    #     - Extracts only H I lines from the linelist in `sme`
    #     - Synthesizes an H-only spectrum using LTE
    #     - Interpolates the precomputed 3D NLTE correction profiles to match 
    #     the current SME model parameters and mu-angle grid
    #     - Computes a correction factor as the ratio of 3D NLTE to LTE intensities
    #     - Applies this correction to the full synthetic spectrum (`sme_res`)
    #     - Returns the corrected spectrum and the full correction matrix

    #     Parameters
    #     ----------
    #     sme : SME_Structure
    #         The SME input object containing model parameters (Teff, logg, FeH, mu, linelist, etc.)

    #     Returns
    #     -------

    #     correction_all : list
    #         A list containing:
    #         [0] - Wavelength array used for correction
    #         [1] - 2D correction factor array: shape (n_mu, n_wavelength)
    #         [2] - 2D I/Ic intensity profile (3D NLTE): shape (n_mu, n_wavelength)

    #     Notes
    #     -----
    #     - Uses global `H_lineprof` as the precomputed hydrogen profile database.
    #     - `safe_interpolation` must support extrapolation with `fill_value=1.0` to avoid numerical issues.
    #     - The correction is only applied within the defined `boundary_vertices`; outside that region, correction = 1.

    #     """
    #     # Generate the synthetic spectra using only the H lines
    #     logger.info(f"Getting H 3dnlte correction")
    #     sme_H_only = SME_Structure()
    #     sme_H_only.teff, sme_H_only.logg, sme_H_only.monh, sme_H_only.vmic, sme_H_only.vmac, sme_H_only.vsini = sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini
    #     sme_H_only.iptype = sme.iptype
    #     sme_H_only.ipres = sme.ipres
    #     # sme_H_only.specific_intensities_only = True
    #     # sme_H_only.normalize_by_continuum = False
    #     for i in range(len(sme.nlte.elements)):
    #         sme_H_only.nlte.set_nlte(sme.nlte.elements[i], sme.nlte.grids[sme.nlte.elements[i]])
    #     # sme_H_only = deepcopy(sme)
    #     sme_H_only.linelist = sme.linelist[(sme.linelist['species'] == 'H 1') | ((sme.linelist['wlcent'] > 6562.8-0.5) & (sme.linelist['wlcent'] < 6562.8+0.5))]
    #     sme_H_only.wave = np.arange(4000, 6700, 0.02)
    #     sme_H_only.tdnlte_H = False
    #     sme_H_only_res = self.synthesize_spectrum(sme_H_only)

    #     # Get the correction
    #     resample_H_all, in_boundary = interpolate_H_spectrum(H_lineprof, sme.teff, sme.logg, sme.monh, boundary_vertices)

    #     if not in_boundary:
    #         logger.info(f"Outside H 3dnlte grid, not performing correction.")
    #         sme.tdnlte_H = False
    #         return None

    #     # Integrate the intensity
    #     mu_array = resample_H_all.loc[resample_H_all['wl'] == resample_H_all.loc[0, 'wl'], 'mu'].values
    #     wmu_array = resample_H_all.loc[resample_H_all['wl'] == resample_H_all.loc[0, 'wl'], 'wmu'].values
    #     nmu = len(mu_array)

    #     sint = []
    #     cint = []
    #     for mu in mu_array:
    #         indices = resample_H_all['mu'] == mu
    #         wint = resample_H_all.loc[indices, 'wl'].values
    #         sint.append(resample_H_all.loc[indices, 'I_interp'].values)
    #         cint.append(resample_H_all.loc[indices, 'Ic_interp'].values)
    #     sint = np.array(sint)
    #     cint = np.array(cint)

    #     wgrid_all = []
    #     sint_all = []
    #     cint_all = []

    #     for indices in [wint < 4500, (wint > 4500) & (wint < 6000), wint > 6000]:
    #         wint_single = wint[indices]
    #         sint_single = sint[:, indices]
    #         cint_single = cint[:, indices]
            
    #         wgrid, vstep = self.new_wavelength_grid(wint_single)
    #         cint_single = self.integrate_flux(mu_array, cint_single, 1, 0, 0, wt=wmu_array, tdnlteH=True)
    #         cint_single = np.interp(wgrid, wint_single, cint_single)
            
    #         y_integrated = np.empty((nmu, len(wgrid)))
    #         for imu in range(nmu):
    #             y_integrated[imu] = np.interp(wgrid, wint_single, sint_single[imu])
    #         sint_single = self.integrate_flux(mu_array, y_integrated, vstep, sme.vsini, sme.vmac, wt=wmu_array, tdnlteH=True)

    #         if "iptype" in sme:
    #             logger.debug("Apply detector broadening")
    #             # ToDo: fit for different resolution in segments (but not necessary?)
    #             ipres = sme.ipres if np.size(sme.ipres) == 1 else sme.ipres[0]
    #             sint_single = broadening.apply_broadening(ipres, wint_single, sint_single, type=sme.iptype, sme=sme)

    #         wgrid_all.append(wgrid)
    #         sint_all.append(sint_single)
    #         cint_all.append(cint_single)

    #     wgrid_all = np.concatenate(wgrid_all)
    #     sint_all = np.concatenate(sint_all)
    #     cint_all = np.concatenate(cint_all)

    #     interpolator = interp1d(wgrid_all, sint_all/cint_all, kind="linear", fill_value=1, bounds_error=False, assume_sorted=True)
    #     correction_all = interpolator(sme_H_only_res.wave[0])
    #     # if np.all(sint_all) == 1:
    #     correction_all = correction_all / sme_H_only_res.synth[0]
    #     interpolator = interp1d(wgrid_all, sint_all, kind="linear", fill_value=1, bounds_error=False, assume_sorted=True)
    #     sint_all = interpolator(sme_H_only_res.wave[0])
    #     interpolator = interp1d(wgrid_all, cint_all, kind="linear", fill_value=1, bounds_error=False, assume_sorted=True)
    #     cint_all = interpolator(sme_H_only_res.wave[0])
    #     return [sme_H_only_res.wave[0], correction_all, sint_all, cint_all, sme_H_only_res.synth[0]]

def synthesize_spectrum(sme, segments="all",**args):
    synthesizer = Synthesizer()
    return synthesizer.synthesize_spectrum(sme, segments, **args)
