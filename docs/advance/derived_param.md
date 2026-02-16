# Derived Parameters

`derived_param` lets you define SME parameters as functions of other parameters
during fitting, without adding them to the free-parameter list.

## What It Does

In `solve(...)`, PySME:

1. sets the current free parameters,
2. evaluates each function in `derived_param`,
3. writes the returned value back to `sme`,
4. runs synthesis and computes residuals.

So derived parameters are updated at every iteration.

## API

```py
from pysme.solve import solve

sme = solve(
    sme,
    param_names=[...],             # free parameters
    derived_param={...},           # derived parameters
)
```

- `derived_param` is a `dict[str, callable]`.
- Key: SME parameter name (for example `"vmic"` or `"abund Mg"`).
- Value: function `f(sme) -> float`.

## Rules and Notes

- A parameter cannot be both free and derived in the same run.
- `dynamic_param` is still accepted as a legacy alias, but it is deprecated.
- For abundance keys (for example `"abund Mg"`), return the **final abundance**
  you want in the usual abundance scale. PySME applies the internal
  `monh` conversion for you when writing into `sme.abund`.

## Example 1: Tie `vmic` to `teff` and `logg`

```py
import numpy as np
from pysme.solve import solve

derived = {
    "vmic": lambda s: np.clip(1.1 + 1e-4 * (s.teff - 5500.0) - 0.3 * (s.logg - 4.0), 0.2, 5.0)
}

fit = ["teff", "logg", "monh", "vsini"]
sme = solve(sme, fit, derived_param=derived)
```

Here `vmic` is never fitted directly; it is recomputed from the current model
state at each iteration.

## Example 2: Enforce fixed `[Mg/Fe]`

```py
from pysme.abund import Abund
from pysme.solve import solve

solar_mg = Abund.solar()["Mg"]

# Target relation: [Mg/Fe] = +0.20  ->  A(Mg) = A_sun(Mg) + [M/H] + 0.20
derived = {
    "abund Mg": lambda s: solar_mg + s.monh + 0.20
}

sme = solve(sme, ["teff", "logg", "monh"], derived_param=derived)
```

