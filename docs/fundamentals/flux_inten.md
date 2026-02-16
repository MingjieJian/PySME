# Flux and intensity

This page summarizes two switches that control how PySME outputs synthetic spectra.

## `normalize_by_continuum`

`normalize_by_continuum` controls whether the synthetic spectrum is divided by
the synthetic continuum:

- `True` (default): output a continuum-normalized spectrum (usually for normalized observations)
- `False`: keep flux-level output (usually for flux-calibrated observations)

Even when this is `False`, continuum fitting can still be enabled through
`cscale_flag`.

## `specific_intensities_only`

`specific_intensities_only` controls whether PySME returns angle-dependent
specific intensities or disk-integrated flux:

- `False` (default): integrate specific intensities over the stellar disk and return flux
- `True`: return specific intensities directly, without flux-level post-processing
    - Note that in this case, `synthesize_spectrum` will only return the specific intensities, not the `sme` object.

This is useful when you want the radiative-transfer output itself (as a function
of angle), rather than only the final integrated spectrum.

## Typical usage

- Normalized stellar spectrum fitting:
  - `normalize_by_continuum = True`
  - `specific_intensities_only = False`
- Flux-calibrated analysis:
  - `normalize_by_continuum = False`
  - `specific_intensities_only = False`
- Intensity-level diagnostics / custom integration:
  - `specific_intensities_only = True`
