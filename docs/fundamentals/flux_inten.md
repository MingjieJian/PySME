# Flux and intensity

This page summarizes two switches that control how PySME outputs synthetic spectra.

## `normalize_by_continuum`

`normalize_by_continuum` controls whether the synthetic spectrum is divided by
the synthetic continuum:

- `True` (default): output a continuum-normalized spectrum (usually for normalized observations)
- `False`: keep flux-level output (usually for flux-calibrated observations)

Even when this is `False`, continuum fitting can still be enabled through
`cscale_flag`.

This flag does not redefine what `cscale` means. `cscale` is still the
multiplicative correction applied to the synthetic spectrum to match the
observation. See [](continuum.md) for the interaction between
`normalize_by_continuum`, `cscale`, and `cscale_type`.

## `specific_intensities_only`

`specific_intensities_only` controls whether PySME returns angle-dependent
specific intensities or disk-integrated flux:

- `False` (default): integrate specific intensities over the stellar disk and return flux
- `True`: return specific intensities directly, without flux-level post-processing
    - In this mode, `synthesize_spectrum` stores trimmed intensity-level outputs on `sme`:
      - `sme.wint`: transfer wavelength grid per segment
      - `sme.sint`: line+continuum specific intensities `(nmu, nwave)`
      - `sme.cint`: continuum specific intensities `(nmu, nwave)`

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
