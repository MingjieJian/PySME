# Continuum

In PySME, continuum handling is split into two related but different pieces:

- `cscale`: the multiplicative correction that PySME applies to the synthetic
  spectrum to match the observation
- `normalize_by_continuum`: whether the synthetic spectrum is first divided by
  its own synthetic continuum before the final matching step

The most important point is:

`cscale` is applied to the synthetic spectrum, not to the observed spectrum.

So `cscale` should usually be interpreted as:

- "the correction needed to bring the synthetic spectrum onto the observation"

and not as:

- "the true continuum estimate of the observed spectrum"

## `cscale`

`cscale` stores one correction per segment.

- For polynomial modes, each segment stores polynomial coefficients
  evaluated on the observed wavelength grid after shifting the first point to
  zero, i.e. `f(wave - wave[0])`
- For spline modes, each segment stores a spline-like multiplicative array on
  the wavelength grid

After radial-velocity matching, PySME multiplies the synthetic spectrum by
`cscale`.

## `cscale_flag`

`cscale_flag` controls the polynomial degree or whether continuum fitting is
performed at all.

- `none`: do not fit a continuum correction
- `fix`: use the current `cscale` values without updating them
- `constant`: zeroth-order multiplicative correction
- `linear`: first-order multiplicative correction
- `quadratic`: second-order multiplicative correction
- higher-order polynomial options continue the same pattern

## `cscale_type`

`cscale_type` controls how PySME determines the correction.

### `mask`

PySME fits the correction only on pixels marked as continuum in the mask.

Use this when:

- you trust your continuum mask
- you want a more traditional continuum normalization workflow

### `match`

PySME fits the correction by matching the synthetic spectrum to the observation
over all good pixels.

Use this when:

- you do not have a reliable continuum mask
- the spectrum is only approximately normalized

### `match+mask`

Like `match`, but restricted to continuum-mask pixels.

This is often useful when:

- you want a synthetic-to-observation match
- but only trust selected continuum regions

### `matchlines`

Variant of `match` that emphasizes line pixels more strongly than continuum
pixels.

### `matchlines+mask`

Like `matchlines`, with mask-based restrictions.

### `spline`

Uses a spline-based multiplicative correction instead of a low-order
polynomial.

### `spline+mask`

Like `spline`, but only over continuum-mask pixels.

### `mcmc`

Joint continuum and radial-velocity fit using the MCMC path.

## What changes with different `cscale_type`

The key difference is not "whether continuum correction exists", but how the
synthetic-to-observation correction is estimated.

| `cscale_type` | Main idea | Uses mask? | Best mental model |
|---|---|---:|---|
| `mask` | Fit continuum from continuum pixels | Yes | continuum normalization from selected regions |
| `match` | Match synthetic to observation on good pixels | No | global multiplicative synthetic correction |
| `match+mask` | Match synthetic to observation on continuum pixels | Yes | masked synthetic correction |
| `matchlines` | Match mainly on line structure | No | line-driven multiplicative correction |
| `matchlines+mask` | line-driven matching with mask restriction | Yes | masked line-driven correction |
| `spline` | flexible multiplicative correction | No | high-flexibility matching |
| `spline+mask` | flexible correction with mask restriction | Yes | masked high-flexibility matching |

## `normalize_by_continuum`

This flag is separate from `cscale`.

- `True`:
  - PySME first divides the synthetic spectrum by the synthetic continuum
  - then applies radial velocity and `cscale`
  - this is the usual mode for normalized or nearly normalized observations
- `False`:
  - PySME keeps the synthetic spectrum on its flux scale
  - then applies radial velocity and `cscale`
  - this is appropriate only when you really want flux-level output

So even with `normalize_by_continuum=False`, PySME can still fit a continuum
correction using `cscale`.

## Why plotting can look confusing

Because `cscale` is applied to the synthetic spectrum, a common mistake is:

```python
cont = np.polyval(coef, x)
plt.plot(spec / cont)
plt.plot(synth / cont)
```

This can be misleading, especially for `match` and `match+mask`, because:

- `spec / cont` treats `cscale` as if it were an observed-spectrum continuum estimate
- `synth / cont` may partially undo the correction PySME already applied to the synthetic spectrum

That means:

- `mask` mode often behaves more like a classical normalization coefficient
- `match` / `match+mask` behave more like a synthetic-to-observation matching correction

So the same plotting recipe does not always have the same meaning across
`cscale_type`.

## Practical plotting guidance

### If you want to inspect the final fitted comparison

Plot the stored outputs directly:

```python
plt.plot(sme.wave[seg], sme.spec[seg])
plt.plot(sme.wave[seg], sme.synth[seg])
```

This shows the actual fitted result produced by PySME.

### If you want to inspect continuum-normalized observations using `mask`

It is reasonable to inspect:

```python
x = sme.wave[seg] - sme.wave[seg][0]
cont = np.polyval(sme.cscale[seg], x)
plt.plot(sme.wave[seg], sme.spec[seg] / cont)
```

because in `mask` mode the fitted correction is close to a traditional
continuum estimate from continuum regions.

### If you use `match` or `match+mask`

Be careful interpreting `spec / cont` as "normalized observation", because the
correction is solving a matching problem, not necessarily estimating the true
continuum of the observed spectrum.

In these modes, the most trustworthy plot is usually still the direct fitted
comparison:

```python
plt.plot(sme.wave[seg], sme.spec[seg])
plt.plot(sme.wave[seg], sme.synth[seg])
```

## Recommended defaults

- Nearly normalized stellar spectra:
  - `normalize_by_continuum = True`
  - start with `cscale_type = "mask"` if you trust your continuum mask
  - otherwise try `cscale_type = "match+mask"` or `match`
- Flux-calibrated spectra:
  - `normalize_by_continuum = False`
  - use with care, because synthetic and observed outputs are no longer on a
    normalized scale by default
