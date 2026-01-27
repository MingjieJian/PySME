# PySME Development Guide

## Quick Start

```bash
uv sync                    # Install in dev mode (compiles native code)
uv run pytest              # Run tests
```

## Project Overview

PySME (Spectroscopy Made Easy) is a Python package for stellar spectroscopy analysis. It includes a native library (`libsme`) written in Fortran/C++ that handles the computationally intensive spectral synthesis.

## Build System

Uses **scikit-build-core** with CMake to compile the native extension:
- `smelib/` - Git submodule containing Fortran/C++ source (from MingjieJian/SMElib)
- `CMakeLists.txt` - Builds `libsme` shared library and `_smelib` Python extension
- `pyproject.toml` - Build configuration and cibuildwheel settings

The native code is compiled at build time (not runtime), producing platform-specific wheels.

### smelib Submodule

The `smelib/` directory is a git submodule pointing to `ivh/SMElib` (forked from `MingjieJian/SMElib`, which is the active upstream). The original repo `AWehrhahn/SMElib` is abandoned.

PySME's CMakeLists.txt compiles these smelib sources:
- `smelib/src/sme/sme_synth_faster.cpp` - main C++ synthesis code
- `smelib/src/sme/*.f`, `smelib/src/eos/*.f` - Fortran routines
- `smelib/pymodule/_smelib.cpp` - Python extension wrapper
- `smelib/src/data/` - runtime data files

Into:
- `libsme.so/.dylib/.dll` - shared library with Fortran routines
- `_smelib` - Python extension module

#### smelib cleanup

The fork has been cleaned of files not needed by PySME:
- Removed: CI configs (`.github/`, `.travis.yml`, `travis/`), test suite, debug scripts, backup files
- Kept: autotools build system (`bootstrap`, `configure.ac`, `Makefile.am`, etc.) and `pymodule/smelib.py` for standalone/IDL use

#### Build system history

**Old approach** (before current build system):
1. `pip install pysme-astro` installed pure Python
2. At runtime, PySME detected platform and downloaded pre-built `libsme` binaries from SMElib GitHub releases
3. Those releases were built by SMElib's own CI workflows (now removed from our fork)

**Current approach**:
1. `pip install pysme-astro` installs a wheel with compiled `libsme` included
2. Compilation happens at wheel build time via PySME's CI (CMake/scikit-build-core)
3. No runtime downloads needed

The standalone SMElib releases (from MingjieJian/SMElib) are still produced for IDL users who use SME without Python.

## CI/CD

GitHub Actions workflow (`.github/workflows/python-app.yml`):
- Triggers on version tags (`v*`) or manual dispatch
- Builds wheels for:
  - Linux x86_64 (manylinux)
  - macOS arm64
  - Windows AMD64
- Publishes to PyPI via trusted publishing (OIDC)

### Platform-specific notes

- **macOS**: Uses Homebrew gcc/gfortran-14, MACOSX_DEPLOYMENT_TARGET=14.0
- **Windows**: Uses MSYS2 MinGW gfortran, Ninja generator, custom wheel repair script
- **Linux**: Uses manylinux container with gfortran

## Release Process

1. Update version in `pyproject.toml`:
   ```
   version = "X.Y.Z"      # stable release
   version = "X.Y.ZaN"    # alpha (not installed by default)
   version = "X.Y.ZbN"    # beta
   version = "X.Y.ZrcN"   # release candidate
   ```

2. Sync lockfile and commit:
   ```bash
   uv sync
   git add pyproject.toml uv.lock
   git commit -m "Release vX.Y.Z"
   git tag vX.Y.Z
   git push && git push origin vX.Y.Z
   ```

3. GitHub Actions builds wheels and publishes to PyPI automatically.

4. If PyPI publish fails (permissions), create GitHub release manually:
   ```bash
   gh run download <run-id> --dir release-artifacts
   gh release create vX.Y.Z --prerelease --title "vX.Y.Z" \
     release-artifacts/sdist/* release-artifacts/wheels-*/*
   ```

## Key Files

- `src/pysme/` - Main Python package
- `src/pysme/smelib/` - Native library interface
- `smelib/` - Fortran/C++ source (submodule)
- `CMakeLists.txt` - Native build configuration
- `pyproject.toml` - Package metadata, build config, cibuildwheel config

## Dependencies

Runtime: numpy, scipy, pandas, astropy, tqdm, emcee, plotly, matplotlib, etc.
Build: scikit-build-core, numpy, CMake, Fortran compiler (gfortran)

## Testing

```bash
uv run pytest                           # Run all tests
uv run python -c "from pysme.sme import SME_Structure"  # Quick import test
```

## Local Development Setup

After cloning, initialize submodules and install:

```bash
git submodule update --init --recursive
uv sync
```

This compiles the native code and installs in editable mode. After editing C/Fortran code in `smelib/`, re-run `uv sync` to rebuild.
