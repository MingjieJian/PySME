# Abundance

`sme.abund` stores elemental abundances, while `sme.monh` stores the global metallicity offset.
In PySME, these two are combined at synthesis time.

## Core rule in PySME

For elements heavier than He, the abundance used by SMElib is:

```
effective abundance = pattern abundance + sme.monh
```

Hydrogen and helium are not shifted by `monh`.

## What this means when you modify abundances

When you update one element in `sme.abund`, you are writing the **pattern**
value (before the global metallicity shift).  
So if you want a target final abundance `A_target` (for example in `H=12`), use:

```
pattern_value = A_target - sme.monh
```

then assign that pattern value.

```{warning}
Note that when printing out `abund`, `monh` is applied to the abundance values.
```

## Practical examples

```py
from pysme.abund import Abund

sme.abund = Abund.solar()
sme.monh = -0.20
```

### 1. Set a target absolute abundance (H=12)

If you want final `A(Mg) = 7.40`:

```py
target = 7.40
sme.abund["Mg"] = target - sme.monh
```

Because PySME adds `monh` later, this produces the intended final Mg abundance.

### 2. Keep scaled-solar composition, only change metallicity

```py
sme.monh = -0.50
```

This shifts all metals (Z > 2) together by `-0.50 dex`, while keeping abundance ratios fixed.

### 3. Set an element enhancement at fixed metallicity

If you want `[Mg/Fe] = +0.30` at current `monh`, add `+0.30` to Mg in the
pattern relative to your base pattern.

<!-- ## Tip for fitting / derived parameters

If your derived function returns a target **final** abundance, subtract `sme.monh`
before writing to `sme.abund[...]`.  
This matches how PySME internally combines pattern + metallicity. -->
