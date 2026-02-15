# Installation

PySME can be installed through PyPI (recommended) or from github directly.

Currently PySME is tested with Python verion 3.9-3.13.

## (optional) set up virtual environment
- `conda env create --name pysme`
    - You can specify python verion using `python=3.12`
- `source activate pysme`

## Install 
- For the "stable" version (recommended):
    - `pip install pysme-astro`

```{warning}
- PySME requires the pre-compled C++/Fortran SME library to run. Currently we deliver SME library with Linux and Mac version; for Windows users, we recommend to use WSL and hence the Linux version. 
- The files (mainly atmosphere models and NLTE departure coefficnent grids) required by PySME will be saved inside `~/.sme/`. These files can be large thus if your home directory is small, we recommend to create a softlink for `~/.sme`.
```

## Running SME
- An simple minimum example is provided in the [examples directory](https://github.com/MingjieJian/SME/tree/master/examples). Make sure to also download the provided input structure.
- You can then run it with: `python minimum.py`

```{admonition} Accessing data files
The atmosphere and nlte data files should be downloaded from the server automatically when used, so network connection is required when using PySME (not only during installation).
These files can also be downloaded manually by:
- Download data files as part of IDL SME from [here](http://www.stsci.edu/~valenti/sme.html).
-  copy them into their respective storage locations in ~/.sme/atmospheres and ~/.sme/nlte_grids
    - atmospheres: everything from SME/atmospheres
- Download the nlte_grids in [zenodo](https://doi.org/10.5281/zenodo.3888393).
```

## Uninstall 

You can uninstall PySME by:
```sh
pip uninstall pysme-astro
```

Note that several files (data file, SMElib file etc) will remain after the uninstall. 
They are all list in the output of the pip command, and it is recommended to remove them manually.