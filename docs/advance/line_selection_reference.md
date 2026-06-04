# Line-Selection Reference

This page is the parameter reference for PySME line selection and line-info
precomputation. It complements [](line_filtering.md), which focuses more on
workflow and examples.

## Mental model

Line selection in PySME is controlled in three layers:

1. `linelist_mode`
   - `all`: synthesize with the full line list
   - `dynamic`: synthesize with a per-segment filtered subset
   - `auto`: deprecated alias of `dynamic`
2. `line_select_method`
   - `internal`: no external CDR/ALMAX metadata path
   - `cdr`: use `central_depth` and `line_range_*`
   - `almax`: use `almax_ratio` and `line_range_*`
3. `line_select_policy`
   - `auto`: use method-dependent automatic line-info handling
   - `strict`: require explicit method-specific line-info handling

In practice:

- `linelist_mode` decides whether dynamic filtering is used
- `line_select_method` decides how line metadata is generated/interpreted
- `line_select_policy` decides how strictly that metadata is enforced

## Shared parameters

### `linelist_mode`

Function argument in `solve(...)` and `synthesize_spectrum(...)`.

- `all`: use the full line list
- `dynamic`: use dynamic line filtering
- `auto`: deprecated alias of `dynamic`

### `sme.line_select_method`

- `internal`
- `cdr`
- `almax`

Controls which metadata path is used for line preselection.

### `sme.line_select_policy`

- `auto`
- `strict`

Controls how method-specific line metadata is consumed by the synthesis path.

### `sme.line_select_parallel`

Boolean. Enable or disable parallel metadata updates.

### `sme.line_select_n_jobs`

Worker count for parallel metadata updates.

- `None`: infer automatically
- positive integer: explicit worker count

### `sme.line_select_chunk_size`

Chunk size used when splitting the line list for metadata updates.

### `sme.line_select_recompute`

- `if_stale`: recompute line metadata only when missing or stale
- `always`: recompute every time
- `never`: do not recompute; require existing metadata or cache entries

This is the main control for metadata regeneration.

### `sme.line_select_stale_thres`

Dictionary of stale thresholds, typically including:

- `teff`
- `logg`
- `monh`
- `vmic`
- `accrt`

Used to decide whether previously computed metadata is still valid.

### `sme.line_precompute_database`

Preferred cache-directory parameter.

This is the shared on-disk cache for line precompute products from both `cdr`
and `almax`.

## CDR-specific parameters

### `sme.line_select_cdr_strength_thres`

Threshold for strong-line selection in CDR mode.

### `sme.line_select_cdr_bin_width`

Bin width used in CDR strong-line selection.

CDR currently uses the bin-based strong-line helper directly, so there is no
separate `cdr_use_bins` switch.

## ALMAX-specific parameters

### `sme.line_select_almax_threshold`

Threshold used by ALMAX-based selection.

If `None`, it falls back to `sme.accrt`.

### `sme.line_select_almax_use_bins`

Boolean switch controlling which ALMAX strong-line rule is used:

- `False`: simple threshold rule using `almax_ratio >= line_select_almax_threshold`
- `True`: bin-wise cumulative rule using `flag_strong_lines_by_bins(...)`

### `sme.line_select_almax_bin_width`

Bin width used when `line_select_almax_use_bins=True`.

## Deprecated or legacy parameters

These are still accepted for backward compatibility, but should not be used in
new code.

### `cdr_database`

Deprecated alias of `line_precompute_database`.

### `sme.cdr_N_line_chunk`

Legacy alias of `sme.line_select_chunk_size`.

### `sme.cdr_parallel`

Legacy alias of `sme.line_select_parallel`.

### `sme.cdr_n_jobs`

Legacy alias of `sme.line_select_n_jobs`.

### `sme.strong_depth_thres`

Legacy alias of `sme.line_select_cdr_strength_thres`.

### `sme.strong_bin_width`

Legacy compatibility field for bin-width based strong-line selection.
Prefer:

- `sme.line_select_cdr_bin_width`
- `sme.line_select_almax_bin_width`

### `sme.line_select_reuse`

Deprecated.

Non-default values currently only trigger a limited internal reuse path by
keeping line opacity around. This is not a fully developed or stable public
cache policy, and new code should leave it at the default `none`.

### `cdr_create`

Legacy-style function argument still used to force regeneration of cached line
metadata products. It remains supported, but is not yet replaced by a clearer
unified name.

## Recommended usage

### CDR workflow

```python
sme.line_select_method = "cdr"
sme.line_select_policy = "strict"
sme.line_select_recompute = "if_stale"
sme.line_select_parallel = False
sme.line_select_chunk_size = 2000
sme.line_select_cdr_strength_thres = 0.001
sme.line_select_cdr_bin_width = 0.2
sme.line_precompute_database = "/path/to/cache"
```

### ALMAX workflow

```python
sme.line_select_method = "almax"
sme.line_select_policy = "strict"
sme.line_select_recompute = "if_stale"
sme.line_select_parallel = False
sme.line_select_chunk_size = 2000
sme.line_select_almax_threshold = sme.accrt
sme.line_select_almax_use_bins = False
sme.line_select_almax_bin_width = 0.2
sme.line_precompute_database = "/path/to/cache"
```

## See also

- [](line_filtering.md)
- [](how-to.md)
