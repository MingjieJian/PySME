## Archived Repository

PySME has moved to the new canonical repository:

<https://github.com/SpectroscopyMadeEasy/PySME>

This repository (`PySME-legacy`) is now legacy/read-only.
Please open new issues and pull requests in the new repository.

---

![Python application](https://github.com/SpectroscopyMadeEasy/PySME-legacy/workflows/Python%20application/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pysme-astro/badge/?version=latest)](https://pysme-astro.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5547527.svg)](https://doi.org/10.5281/zenodo.5547527)

# PySME

Spectroscopy Made Easy (SME) is a software tool that fits an observed
spectrum of a star with a model spectrum. Since its initial release in
[1996](http://adsabs.harvard.edu/abs/1996A%26AS..118..595V), SME has been a
suite of IDL routines that call a dynamically linked library, which is
compiled from C++ and fortran. This classic IDL version of SME is available
for [download](http://www.stsci.edu/~valenti/sme.html).

In 2018, we began began reimplmenting the IDL part of SME in python 3,
adopting an object oriented paradigm and continuous itegration practices
(code repository, build automation, self-testing, frequent builds).

# SMElib source code

PySME is the Python front-end for SME. The underlying SME C++/Fortran library
(`SMElib`) lives in a separate repository:

<https://github.com/SpectroscopyMadeEasy/SMElib>

Most users do not need to download `SMElib` separately because supported PySME
installations provide the compiled library. You may want the `SMElib` source
code if you need to inspect the low-level implementation, debug the native
library, or compile it manually for an unsupported platform.

# Installation

A stable version is available on pip `pip install pysme-astro`, and it is recommended to install this verion.

If you are interested in the latest version you can do so by cloning this git.
```bash
# Clone the git repository
git clone https://github.com/MingjieJian/SME.git
# Move to the new directory
cd SME
# Install this folder (as an editable module)
pip install -e .
```
See also the [documentation](https://pysme-astro.readthedocs.io/en/latest/usage/installation.html).

# Poster

A poster about PySME can be found here: [Poster](http://sme.astro.uu.se/poster.html)

# GUI

A GUI for PySME is available in its own repository [PySME-GUI](https://github.com/MingjieJian/PySME-GUI).

# Windows

Unfortunately PySME is not supported in Windows right now. While there is a SME C libary for Windows, it is not compatible with the Python C Extension inteface on Windows as it was compiled with a different compiler. Therefore if you want to use PySME you would need to compile the SME library with the same compiler.
