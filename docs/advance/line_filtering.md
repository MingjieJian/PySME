# Line Filtering

For wide wavelength coverage (or many segments), using the full line list in
every segment is expensive. PySME provides dynamic line filtering to keep only
relevant lines per segment (see Jian et al. in prep).

## Core Options

Use `linelist_mode` and `line_select_method` together in synthesis or solve:

- `"all"`: use all lines (default).
- `"dynamic"`: filter lines by precomputed line properties (recommended for long spectra).
- `"auto"`: legacy alias of `"dynamic"` (deprecated).
- `line_select_method="internal"`: no external preselection metadata.
- `line_select_method="cdr"`: use `central_depth` + `line_range_*`.
- `line_select_method="almax"`: use `almax_ratio` + `line_range_*`.

## How Dynamic Filtering Works

When `linelist_mode="dynamic"`:

1. PySME ensures line metadata exists for the selected method.
   - `cdr`: `central_depth`, `line_range_s`, `line_range_e`
   - `almax`: `almax_ratio`, `line_range_s`, `line_range_e`
2. If metadata is missing or stale, PySME updates it automatically.
2. For each segment, PySME keeps lines that overlap the segment range
   (with broadening margin) and pass a strength threshold.
3. Only this reduced line subset is sent to SMElib for that segment.

This can significantly reduce runtime for long or segmented spectra.

For a complete parameter-by-parameter reference, including deprecated aliases
and recommended replacements, see [](line_selection_reference.md).

## Main Controls

- Shared controls:
  - `sme.line_select_method`: `internal | cdr | almax`
  - `sme.line_select_policy`: `auto | strict`
  - `sme.line_select_parallel`: enable/disable parallel metadata update
  - `sme.line_select_n_jobs`: number of worker processes
  - `sme.line_select_chunk_size`: chunk size for metadata update
  - `sme.line_select_recompute`: `if_stale | always | never`
  - `sme.line_select_stale_thres`: stale thresholds for model changes
- CDR method controls:
  - `sme.line_select_cdr_strength_thres`
  - `sme.line_select_cdr_bin_width`
  - `line_precompute_database` / `cdr_create` (function args)
- ALMAX method controls:
  - `sme.line_select_almax_threshold`
  - `sme.line_select_almax_use_bins`
  - `sme.line_select_almax_bin_width`

`line_precompute_database` is a shared on-disk cache for both `cdr` and
`almax`. Cache entries are isolated by `(method, linelist_hash, stellar params)`,
so one folder can safely store multiple linelists and both methods together.
Legacy `cdr_database` is still accepted as a deprecated alias.

### Recompute vs. reuse

- `line_select_recompute` controls whether line metadata is recomputed when it
  is missing or stale:
  - `if_stale`: recompute only when needed
  - `always`: always recompute
  - `never`: require existing metadata or cache entries
- `line_select_reuse` is deprecated. Non-default values still enable a limited
  internal reuse path by keeping line opacity around, but this is not a fully
  developed cache policy and should not be treated as a stable public API.

### ALMAX strong-line rule

- `line_select_almax_use_bins=False`:
  - `strong = (almax_ratio >= line_select_almax_threshold)`
- `line_select_almax_use_bins=True`:
  - `strong = flag_strong_lines_by_bins(wl, almax_ratio, threshold=line_select_almax_threshold, bin_width=line_select_almax_bin_width)`

`line_select_almax_threshold` is the single ALMAX threshold parameter for both
rules. If it is `None`, it falls back to `sme.accrt` (legacy-compatible
behavior).

CDR does not have a separate `use_bins` switch because its current strong-line
selection already uses the bin-based helper internally.

## Example 1: Dynamic Filtering in Synthesis (CDR)

```py
from pysme.synthesize import Synthesizer, synthesize_spectrum

synth = Synthesizer()
sme = synth.update_cdr(sme)              # populate central_depth / line_range_* once
sme.line_select_method = "cdr"
sme.line_select_cdr_strength_thres = 0.02

sme = synthesize_spectrum(sme, linelist_mode="dynamic")
```

## Example 2: Dynamic Filtering in Synthesis (ALMAX)

```py
from pysme.synthesize import synthesize_spectrum

sme.line_select_method = "almax"
sme.line_select_almax_threshold = sme.accrt
sme.line_select_almax_use_bins = True
sme.line_select_almax_bin_width = 0.2

sme = synthesize_spectrum(sme, linelist_mode="dynamic")
```

## Example 3: Dynamic Filtering in Solve

```py
from pysme.solve import solve

fit = ["teff", "logg", "monh", "vmic"]
sme.line_select_method = "almax"
sme.line_select_almax_use_bins = True

sme = solve(
    sme,
    fit,
    linelist_mode="dynamic",
    line_precompute_database="path/to/line_precompute_db",
    cdr_create=False,                     # set True to force regeneration
)
```

## Practical Guidance

- Start with `line_select_method="cdr"` for continuity with existing CDR workflows.
- For ALMAX, start with `line_select_almax_threshold = sme.accrt`.
- Enable `line_select_almax_use_bins=True` when you want bin-wise cumulative pruning.
- Use `"all"` for short, narrow windows where filtering overhead may not help.
- Use `"dynamic"` for wide ranges or many segments.
