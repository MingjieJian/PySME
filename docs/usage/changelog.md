# Changelog

This page stores the change log for pysme since May 2024.

## In-development

- (model) modify the atmosphere geometry to '' for auto desicion instead of 'PP'.
- (sme structure) add `sme.wint` as a segment-aware optional transfer grid input for synthesis.
- (synthesize) when `sme.wint` is provided, synthesis now prefers user-provided `wint` over internal cached wavelength grids.
- (docs) update long-spectrum examples to prefer `linelist_mode='dynamic'` (keep `auto` as deprecated compatibility alias).

## In github-repo

## 0.4.189

- (synthesize) synthesize function now support calculating line central-depth and line range, and can use it to select only the relevant lines in synthesize (and thus in solve). 
    - Available in both LTE and NLTE calculation.
- (solve) dynamic parameter function now support abundances.

## v0.4.187

- NLTE grids (from Amarsi et al. 2020) set as default grids.

## v0.4.184

- (vald) add save and merge function to ValdFile.
    - Coupling information ('LS', 'JK' etc) are missed during reading. Now added back.
    - There is no `vmic` column in VALD short extract_all mode, but pysme have. This leads to I/O fail for this kind of line list. Fixed.

## v0.4.183

- (mask) revert mask function back. The modification in v0.4.180 caused a bug on mask and it is fixed now.
- (synthesize) fix the vstep issue when synthesizing sparce spectrum.
- (nlte) add ResetDepartureCoefficients function to ensure correct NLTE calculation
- (solve) use 'vrad' instead of 'v_rad' for the radial velocity fitting.
- (solve) add dynamic parameter function; sme parameter can changes according to other parameters while not being included in the fitting.
- New readthedocs theme (alabaster -> sphinx_book_theme).

## v0.4.180

- ~~fix mask bug. The continuum fitting function is now available.~~ (see v0.4.183)

## v0.4.179

- Add line depth result to sme_structure.
- Add line range result to sme_structure.
- Support python 3.12.
- Support new VALD3 line format.
