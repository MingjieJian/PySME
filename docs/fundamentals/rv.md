# Radial velocity (vrad)

`vrad` in PySME is segment-aware and has two functions (with `vrad_flag`):

- fit RV from observation vs. synthetic spectrum
- apply RV shift to synthetic spectrum before comparison

## Core parameters

- `vrad`: radial velocity in km/s for each segment
- `vrad_flag`: controls how RV is determined
  - `none`: no RV fitting
  - `fix`: use the current `sme.vrad` to shift the synthetic spectra and do not change it
  - `each`: fit one RV per segment
  - `whole`: fit one shared RV from all selected segments