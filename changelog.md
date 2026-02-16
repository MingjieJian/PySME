# Changelog

## 2026-02-10

- Added `derived_param` as the preferred name for derived-parameter callbacks in `solve()`.
- Kept backward compatibility with `dynamic_param`:
  - If `dynamic_param` is used, a `DeprecationWarning` is emitted.
  - If both `derived_param` and `dynamic_param` are provided (and differ), `ValueError` is raised.
- Updated internal solver logic and user-facing messages to use the "derived parameter" terminology.
- Added `smelib_lineinfo_mode` passthrough in `solve()` call paths (`_residuals` and `_jacobian`) so fitting runs can use SMElib precomputed line-info modes.
- Updated `linelist_mode` naming:
  - Preferred values are now `"all"` and `"dynamic"`.
  - `"auto"` is kept as a deprecated compatibility alias and maps to `"dynamic"` with a `DeprecationWarning`.
- Added segment-aware optional input `sme.wint` for synthesis transfer grids.
  - Priority is now `sme.wint[segment]` first, then internal cached grids (when enabled), then SMElib adaptive grid generation.
- Updated user docs accordingly (`sme_struct`, `quickstart`, and `how-to`).
