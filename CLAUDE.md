# PySME Development Guide

## Quick Start

```bash
uv sync                    # Install dependencies (includes dev group)
uv run pytest              # Run tests
uv build                   # Build wheel locally
```

## Project Overview

PySME (Spectroscopy Made Easy) is a Python package for stellar spectroscopy analysis. It includes a native library (`libsme`) written in Fortran/C++ that handles the computationally intensive spectral synthesis.

## Build System

Uses **scikit-build-core** with CMake to compile the native extension:
- `smelib/` - Git submodule containing Fortran/C++ source (from MingjieJian/SMElib)
- `CMakeLists.txt` - Builds `libsme` shared library and `_smelib` Python extension
- `pyproject.toml` - Build configuration and cibuildwheel settings

The native code is compiled at build time (not runtime), producing platform-specific wheels.

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

2. Commit and tag:
   ```bash
   git add pyproject.toml
   git commit -m "Release vX.Y.Z"
   git tag vX.Y.Z
   git push && git push --tags
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
