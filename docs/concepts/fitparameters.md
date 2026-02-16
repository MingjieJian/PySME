# Fitparameters

PySME can use a wide variety of fitparameters, however there are
a few things to consider.

## How to specify fitparameters

Fitparameters can either be passed to `pysme.solve.solve` or as part
of the sme structure. When in doubt the explicitly passed values
overwrite the structure. Fitparameters are passed using their names
as strings (e.g. `'teff'` or `'logg'`), which are usually straightforward.

There are however a few more tricky ones for individual abundances and
linelist parameters:

```{admonition} Abundances
Individual abundances are passed using strings of the form `'abund {El}'`
where `'El'` is the short form element we want to fit (e.g. `'Mg'`).
PySME supports the fitting of all elements up to `'Es'`.
Alternatively you can also use the element number, i.e. the charge number,
instead; e.g. `'3'` for Lithium.
```

```{admonition} Linelist
For linelist parameters use `'linelist {nr} {p}'`, where `'nr'` is the index
in the linelist and `'p'` is the name of the parameter (e.g. `'gflog'`).
Parameters of the long format (e.g. lande_upper) can only be varied
for a long format linelist. For a list of the linelist parameters,
see [](linelist).
```

## Bounds

Most parameters have bounds. PySME gets the bounds of each parameter from the
specified atmosphere grid if possible.
It also only allows for positive turbulence velocities for example.

Continuum and Radial velocities
-------------------------------

It might be logical to think that the continuum and especially the
radial velocity is a regular fitparameters. However there are several
options about how to handle them (see [here](./sme_struct.md#radial-velocity-and-continuum)).
It is therefore recommended to use `'cscale_flag'` and `'vrad_flag'` to
specify the desired fitting method.
If `'vrad'` is passed as a fitparameter it is equivalent to `'vrad_flag' = 'each'`,
and if `'cont'` is passed it is the same as `'cscale_flag' = 'linear'`.
See the [radial velocity](../fundamentals/rv.md) and [continuum](../fundamentals/continuum.md) page for more information.