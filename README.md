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

## From PyPI (upstream)

The `pysme-astro` package on PyPI is published from the upstream repository
[MingjieJian/SME](https://github.com/MingjieJian/SME), not this fork:

```bash
pip install pysme-astro
```

## From this repository

Install a wheel from [GitHub Releases](https://github.com/ivh/PySME/releases):

```bash
pip install https://github.com/ivh/PySME/releases/download/vX.Y.Z/pysme_astro-X.Y.Z-cpXXX-cpXXX-PLATFORM.whl
```

Replace version and platform as needed. Pre-built wheels are available for Linux (x86_64), macOS (arm64), and Windows (AMD64).

## Build from source

Clone with submodules and install:

```bash
git clone --recurse-submodules https://github.com/ivh/PySME.git
cd PySME
uv sync
```

Then use `uv run` to run scripts or tests:

```bash
uv run pytest
uv run python your_script.py
```

After editing C/Fortran code in `smelib/`, re-run `uv sync` to rebuild.

See [CLAUDE.md](CLAUDE.md) for more development details.

# Documentation

See the [documentation](https://pysme-astro.readthedocs.io/en/latest/usage/installation.html) for usage details.

# Poster

A poster about PySME can be found here: [Poster](http://sme.astro.uu.se/poster.html)

# GUI

A GUI for PySME is available in its own repository [PySME-GUI](https://github.com/MingjieJian/PySME-GUI).
