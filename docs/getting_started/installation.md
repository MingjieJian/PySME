# Installation

PySME can be installed through PyPI (recommended; stable release) or from github (the latest version in `develop` branch) directly.

```{admonition} Supported environments:
- Platforms: Linux, macOS
- Windows: supported via WSL2 (install/run PySME inside the Linux subsystem)
- Python versions: 
    - 3.9–3.13.
```

## Set up virtual environment

This step is optional but recommended.

### conda

```bash
conda create -n pysme python=3.12
conda activate pysme
```

### venv (alternative)

```bash
python -m venv .venv
source .venv/bin/activate
```

## Install PySME

### Stable release (recommended)

```bash
pip install pysme-astro
```

### From Github

#### Clone the repository

```bash
git clone https://github.com/MingjieJian/SME.git
cd SME
```

#### Install PySME from source

```bash
pip install -U pip
pip install .
```

## Verify installation

```bash
python -c "import pysme; print('PySME version:', pysme.__version__)"
```

You should see an output of `PySME version: [version]`.

## Uninstall 

You can uninstall PySME by:
```sh
pip uninstall pysme-astro
```

```{warning}
Note that several files (data file, SMElib file etc) will remain after the uninstall. 
They are all list in the output of the pip command, and it is recommended to remove them manually.
```

The content below to be removed.

## Running SME
- An simple minimum example is provided in the [examples directory](https://github.com/MingjieJian/SME/tree/master/examples). Make sure to also download the provided input structure.
- You can then run it with: `python minimum.py`

```{warning}
PySME requires a **pre-compiled C++/Fortran SME library**.
Wheels are currently provided for Linux and macOS.  
For Windows, we recommend using **WSL2** (install PySME inside the Linux environment).

- PySME requires the pre-compled C++/Fortran SME library to run. Currently we deliver SME library with Linux and Mac version; for Windows users, we recommend to use WSL and hence the Linux version. 

```