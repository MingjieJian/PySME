# Line Filtering

For wide wavelength coverage (or many segments), using the full line list in
every segment is expensive. PySME provides dynamic line filtering to keep only
relevant lines per segment (see Jian et al. in prep).

## Core Option

Use `linelist_mode` in synthesis or solve:

- `"all"`: use all lines (default).
- `"dynamic"`: filter lines by precomputed line properties (recommended for long spectra).
- `"auto"`: legacy alias of `"dynamic"` (deprecated).

## How Dynamic Filtering Works

When `linelist_mode="dynamic"`:

1. PySME ensures line metadata exists (`central_depth`, `line_range_s`, `line_range_e`).
   If missing or stale, it updates them via `update_cdr(...)`.
2. For each segment, PySME keeps lines that overlap the segment range
   (with broadening margin) and pass a strength threshold.
3. Only this reduced line subset is sent to SMElib for that segment.

This can significantly reduce runtime for long or segmented spectra.

## Main Controls

- `sme.cdr_depth_thres`: minimum line-strength threshold used in filtering.
- `sme.cdr_N_line_chunk`: chunk size used in `update_cdr`.
- `sme.cdr_parallel`: enable/disable parallel `update_cdr`.
- `sme.cdr_n_jobs`: number of parallel jobs.
- `cdr_database` / `cdr_create` (function args): reuse or build a CDR grid on disk.

## Example 1: Dynamic Filtering in Synthesis

```py
from pysme.synthesize import Synthesizer, synthesize_spectrum

synth = Synthesizer()
sme = synth.update_cdr(sme)              # populate central_depth / line_range_* once
sme.cdr_depth_thres = 0.02               # keep only stronger lines

sme = synthesize_spectrum(sme, linelist_mode="dynamic")
```

## Example 2: Dynamic Filtering in Solve

```py
from pysme.solve import solve

fit = ["teff", "logg", "monh", "vmic"]
sme = solve(
    sme,
    fit,
    linelist_mode="dynamic",
    cdr_database="path/to/cdr_grid",     # optional on-disk CDR cache/grid
    cdr_create=False,                     # set True to force regeneration
)
```

## Practical Guidance

- Start with `sme.cdr_depth_thres = 0.0` and increase gradually if needed.
- Use `"all"` for short, narrow windows where filtering overhead may not help.
- Use `"dynamic"` for wide ranges or many segments.

