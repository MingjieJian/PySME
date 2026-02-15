# The SME Structure

The SME structure contains all the information necessary to run PySME.  
It is by design similar to the IDL structure of SME, but with a lot of
changes for ease of use. Note that some parameters depend on the number of
segments in the observation. When creating a new structure it is therefore
recommended to define the wavelengths first.

## Stellar Parameters

Stellar parameters describe the star in general and are usually what we want to fit.

- `teff`: The effective temperature of the star in **Kelvin**.
- `logg`: The surface gravity of the star in log (cgs).
- `monh`: The overall metallicity of the star in log₁₀ relative to the individual abundances (see {ref}`abund`).
- `vsini`: The (projected) rotation velocity in km s⁻¹. Describes the rotational broadening.  
  - If `vsini` is non-zero, you should set several values for `mu`.
- `vmic`: The micro-turbulence velocity in km s⁻¹. Describes turbulence on scales **smaller** than the photon mean free path, adding additional line broadening.
- `vmac`: The macro-turbulence velocity in km s⁻¹. Describes turbulence on scales **larger** than the photon mean free path, also contributing to broadening.
- `mu`: Limb-angle values ($\mu = \cos\theta$) at which the radiative-transfer calculation is performed. $\mu = 1$ corresponds to disk center, $0$ to the limb.


## Radial velocity and Continuum

The radial velocity and continuum shift the continuum
in wavelength and intensity direction respectively.
How to best handle them often depends on the situation
and the quality of the observation data. Therefore
PySME has many options to determine them.

- `cscale`: The polynomial parameters for each segment, that are applied to the synthetic spectrum to match the observation.
    - The polynomial is calculated using the wavelength grid of the observation, shifted so that the first point is 0. I.e. the polynomial is f(wave - wave[0]).
- `vrad`: The radial velocity in km/s that was applied to each segment of the synthethized spectrum to match the observation.
- `cscale_flag`: Determines how the continuum is fitted, or if it is fitted at all.
    - none: No continuum correction
    - fix: Use whatever continuum scale has been set, but don't change it.
    - constant: scale everything by a factor
    - linear: First order polynomial, i.e. a straight line
    - quadratic: Second order polynomial
- `vrad_flag`: Determines how the radial velocity is fitted, or if it is fitted at all.
    - none: No radial velocity fitting
    - fix: Use the set value for the radial velocity, but don't change it
    - each: Fit each wavelength segment individally
    - whole: Fit the whole spectrum at once

- `normalize_by_continuum`:
    A flag that determines, whether the synthetic flux should be normalized
    by the continous intensities or not. As long as you have a normalized
    observation this should be True, but if you have a flux calibrated
    spectrum this should be set to False. Note that even if this is False,
    you can still fit a continuum normally using cscale_flag.

Spectra
-------

Spectra are given as a list of arrays[^#], where each array represents
one wavelength segment of the spectrum. If there is only one segment,
the list will only have one element. For legacy reasons there is also
an interface to the 'old' system and names (e.g. smod instead of synth)
from IDL SME. It is recommend however to use the new variables.

:wave: The wavelength grid of the observation and/or the synthetic spectrum
:spec: The observed spectrum
:uncs:
    The uncertainties of the observed spectrum. If None will use
    uncertainty 1 for all points.
:mask:
    The pixel mask for the observation, with defination of 0 as bad pixel, 1 as line pixel, 2 as continuum pixel, and 4 as vrad pixel. The masks are additive, i.e., you can set mask value to 5 for line and vrad pixel. Only the good pixels will contribute to the fit, but the synthetic spectrum will still be calculated.
    Note that ``dtype`` of ``sme.mask`` must be int.
:synth: The synthetic spectrum

:wran:
    The first and last wavelength of each segment. You only need
    this if you dont have an observation and dont want to
    specify the exact wavelength grid of the synthetic
    observation. Note that this is not an Illiffe vector.
:wint:
    Optional wavelength grid passed to SMElib `Transf` as adaptive grid for each segment.
    Accepts the same segment-aware input styles as `wave`:
    `numpy.ndarray` (single segment), `list` of arrays (multi-segment),
    or `Iliffe_vector`.
    If provided, `sme.wint[segment]` is used directly for synthesis.
    If not provided, PySME can reuse an internal cached adaptive grid
    when `reuse_wavelength_grid=True`; otherwise SMElib computes a new
    adaptive grid.


Abundance
---------

The individal abundances are stored in a seperate Abundance
object, which shares the same metallicity as the overall structure.
For more detailed information see [](abundance).

:abund: The abundance object

Linelist
--------

The sme structure does contain the whole linelist in the linelist property.
For legacy reasons, it also provides direct access to
the 'species' and 'atomic' arrays. They refer directly to the linelist however.
For more detailed information see :ref:`linelist`.

:linelist: The linelist object
:species: Names of the species of each spectral line
:atomic:
    Atomic linelist data with columns "atom_number", "ionization",
    "wlcent", "excit", "gflog", "gamrad", "gamqst", "gamvw"

Atmosphere
----------

Unlike the linelist the atmosphere is stored in an external file,
that is only referenced by name in the structure.
For more detailed information see :ref:`atmosphere`.

:atmo: The atmosphere object

NLTE
----

Unlike the linelist, but similar to the atmosphere, the NLTE
parameters are stored in external tables, which are only referenced
by name. For more detailed information see :ref:`nlte`.

:nlte: The NLTE object

Instrument Parameters
---------------------

PySME can also model instrumental broadening as part of the
spectral synthesis. For this you need to specify the resolution
and the broadening method to use.

:ipres: The resolution of the instrument to simulate
:iptype:
    The broadening profile of the instrument.
    One of "gauss", "sinc", "table"
:ip_x:
    The x points of the instrument profile.
    Only relevant if iptype is 'table'.
:ip_y:
    The y points of the instrument profile.
    Only relevant if iptype is 'table'.

Fitresults
----------

:fitparameters:
    The fitparameters used for the fitting.
    See :ref:`fitparameters`.
:fitresults: The fitresults object. See :ref:`fitresults`.

System Information
------------------

The sme structure does contain information about the host system.
E.g. which operating system was used. This is mostly for
legacy reasons, and potential debugging information.
For more information see :ref:`system_info`.

:system_info: The system information object. It replaces the idlver object.

Other Parameters
----------------

:gam6: van der Waals scaling factor (usually 1)
:h2broad: flag determing whether to use H2 broadening or not (usually True)
:accrt:
    Minimum accuracy for synthethized spectrum at wavelength grid
    points in `sme.wint` (or SMElib adaptive grid if `sme.wint` is not set).
    Values below 1e-4 are not meaningful.
:accwi:
    Minimum accuracy for linear spectrum interpolation vs. wavelength.
    Values below 1e-4 are not meaningful.
:version: The version of sme used to create this structure and spectrum
:id:
    The date and time when this structure or the
    last synthetic spectrum was created

[^#] They are called Illiffe vectors in the code, and they were that in IDL. But they are technically not Illiffe vectors anymore, but just lists of individal numpy arrays.
