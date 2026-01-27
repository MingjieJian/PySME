# -*- coding: utf-8 -*-
import logging

import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import convolve

logger = logging.getLogger(__name__)


def apply_broadening(ipres, x_seg, y_seg, type="gauss", sme=None):
    """
    Broaden a spectrum by instrument resolution, with a given broadening type

    Parameters
    ----------
    ipres : float
        instrument resolution
    x_seg : array (npoints,)
        x values (wavelength) of the spectrum to broaden
    y_seg : array (npoints,)
        y values (intensities) of the spectrum to broaden
    type : {str, None}, optional
        broadening type to apply. Options are "gauss", "sinc", "table", None.
        If None, will try to use sme.iptype to determine type. "table" requires
        keyword sme to be passed as well. See functions the respective for details.
        (default: "gauss")
    sme : SME_Struct, optional
        sme structure with instrument profile data, required only for
        type="table" or type=None (default: None)

    Raises
    ------
    AttributeError
        if type requires SME_Struct, but its missing
        passed type not recognized

    Returns
    -------
    y_seg : array (npoints,)
        broadened intensity spectrum
    """

    if sme is None and type in ["table", None]:
        raise AttributeError(f"SME structure needs to be passed when using type={type}")

    # Using the log-linear wavelength grid requires using the first point
    # for specifying the width of the instrumental profile
    hwhm = 0.5 * x_seg[0] / ipres if ipres > 0 else 0

    if type is None:
        type = sme.iptype
    type = type.casefold()

    if type == "table":
        y_seg = tablebroad(x_seg, y_seg, sme.ip_x, sme.ip_y)
    elif type == "gauss":
        y_seg = gaussbroad(x_seg, y_seg, hwhm)
    elif type == "sinc":
        y_seg = sincbroad(x_seg, y_seg, hwhm)
    else:
        raise AttributeError(f"Unknown instrument profile type - {type}")

    return y_seg


def tablebroad(w, s, xip, yip):
    """
    Convolves a spectrum with an arbitrary instrumental profile.

    Parameters
    -------
    w : array of size (n,)
        wavelength scale of spectrum to be smoothed
    s : array of size (n,)
        spectrum to be smoothed
    xip : array of size (m,)
        x points of the instrument profile
    yip : array of size (m,)
        y points of the instrument profile
    Returns
    -------
    sout: array[n]
        the smoothed spectrum.
    """
    """
    History
    -------
        22-May-92 JAV
            Switched instrumental profile from multiple gaussians
            to gaussian with power-law wings.
        04-Aug-92 JAV
            Renamed from ipsmo.pro# changed f/ procedure to function.
            Switched f/ 10 to 15 Hamilton pixels in each wing.
        20-Oct-92 JAV
            Switched from gpfunc to ipfun (4 to 5 par).
        23-Aug-94 JAV
            Switched to explicitly passed IPs.
        Oct-18  AW
            Python version
    """

    # Define sizes
    dsdh = np.abs(np.min(np.diff(xip)))
    nip = 2 * int(15 / dsdh) + 1  ## profile points

    # Generate instrumental profile on model pixel scale.
    x = (
        np.arange(nip, dtype=float) - (nip - 1) / 2
    ) * dsdh  # offset in Hamilton pixels
    ip = interp1d(xip, yip, kind="cubic")(x)
    # ip = bezier_interp(xip, yip, x)  # spline onto new scale
    ip = ip[::-1]  # reverse for convolution
    ip = ip / np.sum(ip)  # ensure unit area

    # Pad spectrum ends to minimize impact of Fourier ringing.
    sout = convolve(s, ip, mode="nearest")

    return sout  # return convolved spectrum


def gaussbroad(w, s, hwhm):
    """
    Smooths a spectrum by convolution with a gaussian of specified hwhm.

    Parameters
    -------
    w : array[n]
        wavelength scale of spectrum to be smoothed
    s : array[n]
        spectrum to be smoothed
    hwhm : float
        half width at half maximum of smoothing gaussian.

    Returns
    -------
    sout: array[n]
        the gaussian-smoothed spectrum.
    """
    """
    History
    --------
        Dec-90 GB,GM
            Rewrote with fourier convolution algorithm.
        Jul-91 AL
            Translated from ANA to IDL.
        22-Sep-91 JAV
            Relaxed constant dispersion check# vectorized, 50% faster.
        05-Jul-92 JAV
            Converted to function, handle nonpositive hwhm.
        Oct-18 AW
            Python version
    """
    # Ensure hwhm is a scalar
    hwhm = np.atleast_1d(hwhm).flat[0]

    # Warn user if hwhm is negative.
    if hwhm < 0:
        logger.warning("Forcing negative smoothing width to zero.")

    # Return input argument if half-width is nonpositive.
    if hwhm <= 0:
        return s  # true: no broadening

    # Calculate (uniform) dispersion.
    nw = len(w)  ## points in spectrum
    wrange = w[-1] - w[0]
    dw = wrange / (nw - 1)  # wavelength change per pixel

    # Make smoothing gaussian; extend to 4 sigma.
    # 4.0 / sqrt(2.0 * alog(2.0)) = 3.3972872
    # sqrt(alog(2.0)) = 0.83255461
    # sqrt(alog(2.0) / pi) = 0.46971864
    # (*1.0000632 to correct for >4 sigma wings)
    if hwhm >= 5 * wrange:
        return np.full(nw, np.sum(s) / nw)
    ## points in half gaussian
    nhalf = int(3.3972872 * hwhm / dw)
    ## points in gaussian (odd!)
    ng = 2 * nhalf + 1
    # wavelength scale of gaussian
    wg = dw * (np.arange(ng, dtype=float) - (ng - 1) / 2)
    # convenient absisca
    xg = (0.83255461 / hwhm) * wg
    # unit area gaussian w / FWHM
    gpro = (0.46974832 * dw / hwhm) * np.exp(-xg * xg)
    gpro = gpro / np.sum(gpro)

    # Pad spectrum ends to minimize impact of Fourier ringing.
    sout = convolve(s, gpro, mode="nearest")

    return sout


def sincbroad(w, s, hwhm):
    """
    Smooths a spectrum by convolution with a sinc function of specified hwhm.

    Parameters
    ------
    w : array of size (n,)
        wavelength scale of spectrum to be smoothed
    s : array of size (n,)
        spectrum to be smoothed
    hwhm : float
        half width at half maximum of smoothing gaussian.

    Returns
    -------
    sout : array of size (n,)
        the sinc-smoothed spectrum.
    """
    """
    History
    -------
    Dec-90 GB,GM
        Rewrote with fourier convolution algorithm.
    Jul-91 AL
        Translated from ANA to IDL.
    22-Sep-91 JAV
        Relaxed constant dispersion check# vectorized, 50% faster.
    05-Jul-92 JAV
        Converted to function, handle nonpositive hwhm.
    14-Nov-93 JAV
        Adapted from macbro.pro
    23-Apr-93 JAV
        Verified that convolution kernel has specified hwhm. For IR FTS
        spectra: hwhm=0.0759 Angstroms, max change in profile is 0.4% of continuum.
    Oct-18 AW
        Python Version
    """

    # Warn user if hwhm is negative.
    if hwhm < 0:
        logger.warning("Forcing negative smoothing width to zero.")

    # Return input argument if half-width is nonpositive.
    if hwhm <= 0:
        return s  # true: no broadening

    # Calculate (uniform) dispersion.
    nw = len(w)  ## points in spectrum
    dw = (w[-1] - w[0]) / (nw - 1)  # wavelength change per pixel

    # Make sinc function out to 20th zero-crossing on either side. Error due to
    # ignoring additional lobes is less than 0.2% of continuum. Reducing extent
    # to 10th zero-crossing doubles maximum error.
    fwhm = 2.0 * hwhm  # full width at half maximum
    rperfw = 0.26525  # radians per fwhm of sinc
    xrange = 20 * np.pi  # 20th zero of sinc (radians)
    wrange = xrange * fwhm * rperfw  # 20th zero of sinc (wavelength)
    nhalf = int(wrange / dw + 0.999)  ## points in half sinc
    nsinc = 2 * nhalf + 1  ## points in sinc (odd!)
    wsinc = (np.arange(nsinc, dtype=float) - nhalf) * dw  # absissca (wavelength)
    xsinc = wsinc / (fwhm * rperfw)  # absissca (radians)
    xsinc[nhalf] = 1.0  # avoid divide by zero
    sinc = np.sin(xsinc) / xsinc  # calculate sinc
    sinc[nhalf] = 1.0  # insert midpoint
    xsinc[nhalf] = 0.0  # fix xsinc
    sinc = sinc / np.sum(sinc)  # normalize sinc

    # Pad spectrum ends to minimize impact of Fourier ringing.
    sout = convolve(s, sinc, mode="nearest")

    return sout
