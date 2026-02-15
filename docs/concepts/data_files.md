Data Files You Need
===================

PySME uses local cache data under `~/.sme/`.

Typical files:

- Atmosphere grids
- NLTE grids
- Optional custom resources

```{admonition} Accessing data files
The atmosphere and nlte data files should be downloaded from the server automatically when used, so network connection is required when using PySME (not only during installation).
These files can also be downloaded manually by:
- Download data files as part of IDL SME from [here](http://www.stsci.edu/~valenti/sme.html).
-  copy them into their respective storage locations in ~/.sme/atmospheres and ~/.sme/nlte_grids
    - atmospheres: everything from SME/atmospheres
- Download the nlte_grids in [zenodo](https://doi.org/10.5281/zenodo.3888393).
- The files (mainly atmosphere models and NLTE departure coefficnent grids) required by PySME will be saved inside `~/.sme/`. These files can be large thus if your home directory is small, we recommend to create a softlink for `~/.sme`.
```
