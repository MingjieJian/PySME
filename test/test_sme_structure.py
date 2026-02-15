# -*- coding: utf-8 -*-
from os import remove
from os.path import dirname

import numpy as np
import pytest

from pysme.iliffe_vector import Iliffe_vector
from pysme.sme import SME_Structure as SME_Struct


@pytest.fixture
def cwd():
    return dirname(__file__)


@pytest.fixture
def filename(cwd):
    fname = "{}/__test.sme".format(cwd)
    yield fname
    try:
        remove(fname)
    except:
        pass


def test_empty_structure():
    """Test that all properties behave well when nothing is set"""
    empty = SME_Struct()

    assert isinstance(empty.version, str)
    assert empty.teff is not None
    assert empty.logg is not None
    assert empty.vmic == 0
    assert empty.vmac == 0
    assert empty.vsini == 0

    assert empty.nseg == 0
    assert empty.wave is None
    assert empty.spec is None
    assert empty.uncs is None
    assert empty.synth is None
    assert empty.cont is None
    assert empty.mask is None
    assert empty.mask_good is None
    assert empty.mask_bad is None
    # assert empty.mask_line is None
    # assert empty.mask_continuum is None

    assert empty.cscale.shape == (0, 1)
    assert empty.vrad.shape == (0,)
    assert empty.cscale_flag == "none"
    assert empty.vrad_flag == "none"
    assert empty.cscale_degree == 0

    assert empty.mu is not None
    assert empty.nmu == 7

    # assert empty.md5 is not None

    assert empty.linelist is not None
    assert empty.species is not None
    assert len(empty.species) == 0
    assert empty.atomic is not None

    assert empty.monh == 0
    assert not np.isnan(empty["abund Fe"])
    assert empty.abund["H"] == 12
    assert not np.isnan(empty.abund()["Mg"])

    assert empty.system_info is not None
    assert empty.system_info.arch == ""

    assert len(empty.fitparameters) == 0
    assert empty.fitresults is not None
    assert empty.fitresults.covariance is None

    assert empty.atmo is not None
    assert empty.atmo.depth is None

    assert empty.nlte is not None
    assert empty.nlte.elements == []


def test_save_and_load_structure(filename):
    sme = SME_Struct()
    assert sme.teff is not None

    sme.teff = 5000
    sme.save(filename)
    del sme
    sme = SME_Struct.load(filename)
    assert sme.teff == 5000

    remove(filename)

    data = np.linspace(1000, 2000, 100)
    sme.wave = data
    sme.spec = data
    sme.save(filename)
    sme = SME_Struct.load(filename)
    assert np.all(sme.wave[0] == data)
    assert np.all(sme.spec[0] == data)
    assert sme.nseg == 1


def test_load_idl_savefile(cwd):
    filename = "{}/testcase1.inp".format(cwd)
    sme = SME_Struct.load(filename)

    assert sme.teff == 5770
    assert sme.wave is not None

    assert sme.nseg == 1
    assert sme.cscale_flag == "linear"
    assert sme.vrad_flag == "each"


def test_cscale_degree():
    sme = SME_Struct()
    sme.cscale = 1

    flags = ["none", "fix", "constant", "linear", "quadratic"]
    degrees = [0, 0, 0, 1, 2]

    for f, d in zip(flags, degrees):
        sme.cscale_flag = f
        assert sme.cscale_degree == d
        assert sme.cscale.shape[0] == 0
        assert sme.cscale.shape[1] == d + 1


def test_idlver():
    sme = SME_Struct()
    sme.system_info.update()
    # assert sme.idlver.arch == "x86_64"


def test_fitresults():
    sme = SME_Struct()
    sme.fitresults.chisq = 100
    sme.fitresults.clear()
    assert sme.fitresults.chisq is None


def test_wint_single_segment_array_input():
    sme = SME_Struct()
    wint = np.linspace(5000.0, 5001.0, 11)

    sme.wint = wint

    assert isinstance(sme.wint, Iliffe_vector)
    assert sme.wint.nseg == 1
    assert np.allclose(sme.wint[0], wint)


def test_wint_multi_segment_list_input():
    sme = SME_Struct()
    wint0 = np.linspace(5000.0, 5001.0, 5)
    wint1 = np.linspace(6000.0, 6002.0, 7)

    sme.wint = [wint0, wint1]

    assert isinstance(sme.wint, Iliffe_vector)
    assert sme.wint.nseg == 2
    assert np.allclose(sme.wint[0], wint0)
    assert np.allclose(sme.wint[1], wint1)
