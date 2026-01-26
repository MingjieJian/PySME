![Python application](https://github.com/ivh/PySME/workflows/Python%20application/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pysme-astro/badge/?version=latest)](https://pysme-astro.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5547527.svg)](https://doi.org/10.5281/zenodo.5547527)

# PySME

Spectroscopy Made Easy (SME) is a software tool that fits an observed
spectrum of a star with a model spectrum. Since its initial release in
[1996](http://adsabs.harvard.edu/abs/1996A%26AS..118..595V), SME has been a
suite of IDL routines that call a dynamically linked library, which is
compiled from C++ and Fortran. This classic IDL version of SME is available
for [download](http://www.stsci.edu/~valenti/sme.html).

In 2018, we began reimplementing the IDL part of SME in Python 3,
adopting an object oriented paradigm and continuous integration practices
(code repository, build automation, self-testing, frequent builds).

# Installation

Install from PyPI:
```bash
pip install pysme-astro
```

Pre-built wheels are available for Linux (x86_64), macOS (arm64), and Windows (AMD64).

See the [documentation](https://pysme-astro.readthedocs.io/en/latest/usage/installation.html) for more details.

# Poster

A poster about PySME can be found here: [Poster](http://sme.astro.uu.se/poster.html)

# GUI

A GUI for PySME is available in its own repository [PySME-GUI](https://github.com/MingjieJian/PySME-GUI).

# Development

Clone with submodules and build:
```bash
git clone --recurse-submodules https://github.com/ivh/PySME.git
cd PySME
uv sync
uv build
```

Run tests:
```bash
uv run pytest
```

After editing C/Fortran code in `smelib/`, rebuild with `uv build` or run `cmake --build .` in the `build/` directory.

See [CLAUDE.md](CLAUDE.md) for more development details.
