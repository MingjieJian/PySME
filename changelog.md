# Changelog

## 2026-02-10

- Added `derived_param` as the preferred name for derived-parameter callbacks in `solve()`.
- Kept backward compatibility with `dynamic_param`:
  - If `dynamic_param` is used, a `DeprecationWarning` is emitted.
  - If both `derived_param` and `dynamic_param` are provided (and differ), `ValueError` is raised.
- Updated internal solver logic and user-facing messages to use the "derived parameter" terminology.
- Added `smelib_lineinfo_mode` passthrough in `solve()` call paths (`_residuals` and `_jacobian`) so fitting runs can use SMElib precomputed line-info modes.
