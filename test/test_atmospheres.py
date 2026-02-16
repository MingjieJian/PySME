# -*- coding: utf-8 -*-
# TODO implement atmosphere tests
import numpy as np
import pytest

from pysme.atmosphere.atmosphere import Atmosphere
from pysme.atmosphere.interpolation import AtmosphereInterpolator
from pysme.atmosphere.savfile import SavFile

from .test_largefilestorage import lfs_atmo, skipif_lfs


@pytest.fixture
def atmosphere_name():
    # TODO iterate over all possible options
    return "marcs2012p_t1.0.sav"


@pytest.fixture
def atmosphere(atmosphere_name):
    atmo = Atmosphere(source=atmosphere_name, method="grid", interp="TAU", depth="RHOX")
    return atmo


@pytest.fixture
def interpolator(atmosphere, lfs_atmo):
    interp = AtmosphereInterpolator(
        depth=atmosphere.depth,
        interp=atmosphere.depth,
        geom=atmosphere.geom,
        lfs_atmo=lfs_atmo,
    )
    return interp


@pytest.fixture
def atmosphere_grid(atmosphere_name, lfs_atmo):
    name = lfs_atmo.get(atmosphere_name)
    atmo = SavFile(name, source=name)
    return atmo


@skipif_lfs
@pytest.mark.usefixtures("lfs_atmo")
def test_grid_point(atmosphere_name, atmosphere_grid, lfs_atmo, interpolator):
    # TODO: get this values from the grid
    teff = 7000
    logg = 4
    monh = 0

    atmo_interp = interpolator.interp_atmo_grid(atmosphere_name, teff, logg, monh)
    atmo_grid = atmosphere_grid.get(teff, logg, monh)

    assert np.allclose(atmo_interp.temp, atmo_grid.temp[1:])
    assert np.allclose(atmo_interp.tau, atmo_grid.tau[1:])
    assert np.allclose(atmo_interp.rhox, atmo_grid.rhox[1:])
    assert np.allclose(atmo_interp.rho, atmo_grid.rho[1:])
    assert np.allclose(atmo_interp.xna, atmo_grid.xna[1:])
    assert np.allclose(atmo_interp.xne, atmo_grid.xne[1:])
