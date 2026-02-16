# Spectral segment

A segment is one wavelength chunk of the spectrum.  
PySME is segment-aware by design: most spectral arrays are stored as a list-like
object with one entry per segment.

## Why segments exist

Segments are useful when:

- your observation is naturally split into orders/chunks
- different chunks have different wavelength sampling/resolution/radial velocity.
- you want to process only part of the spectrum in synthesis/solve

## How segments are defined

PySME can get segment information from either:

- `sme.wave`: explicit wavelength arrays (recommended when available)
- `sme.wran`: only segment boundaries `[[w0_start, w0_end], [w1_start, w1_end], ...]`

If `sme.wave` exists, it effectively defines segment boundaries.  
If only `sme.wran` is given, PySME synthesizes each segment within those ranges.

## Segment-aware inputs

Typical segment-aware fields are:

- `wave`
- `spec`
- `uncs`
- `mask`
- `synth`
- `cont`
- `wint` (optional transfer grid)
- `vrad`
- `cscale`

You can pass these as:

- a single 1D array for one segment
- a list of arrays for multiple segments

## Parameters that are not segmented

Some model parameters are global (single value for all segments), e.g.:

- `vmic`
- `vmac`
- `vsini`

Some are scalar-or-per-segment, e.g.:

- `ipres`: can be one value or one value per segment
- `vrad`: handled per segment, depending on `vrad_flag`

## Choosing which segments to run

Both `synthesize_spectrum(...)` and `solve(...)` accept a `segments` argument.

- `segments="all"`: run all segments
- `segments=[0, 2, ...]`: run selected segments only

Invalid indices raise an error.  
Segments that are fully masked as bad pixels are skipped automatically.

## Practical notes

- Keep segment ordering consistent across all segment-aware inputs.
- If you provide per-segment arrays, lengths must match within each segment.
- For simple single-chunk workflows, using one 1D `wave` array is enough.
