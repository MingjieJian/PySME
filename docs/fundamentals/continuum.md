# Continuum

In PySME, the continuum level of the synthetic spectrum can be adjusted to match the observed spectrum using `cscale` and `cscale_flag`.

- `cscale`: Polynomial coefficients (per segment) applied to the synthetic spectrum.
  The polynomial is evaluated on the observed wavelength grid after shifting the first point to zero, i.e., `f(wave - wave[0])`.
- `cscale_flag`: Controls whether and how continuum correction is applied.
  - `none`: No continuum correction.
  - `fix`: Use the current `cscale` values without updating them.
  - `constant`: Scale by a constant factor.
  - `linear`: First-order polynomial (straight line).
  - `quadratic`: Second-order polynomial.
