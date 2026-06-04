# -*- coding: utf-8 -*-
from os.path import dirname, join

import numpy as np
import pytest
from scipy.constants import speed_of_light

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile

# TODO create various kinds of default sme structures
# then run test on all of the relevant ones


@pytest.fixture
def sme_empty():
    from pysme.sme import SME_Structure as SME_Struct

    sme = SME_Struct()
    return sme


@pytest.fixture
def testcase1():
    from pysme.sme import SME_Structure as SME_Struct
    from pysme.synthesize import synthesize_spectrum

    c_light = speed_of_light * 1e-3

    # TODO get better test case for this
    cwd = dirname(__file__)
    fname = join(cwd, "testcase1.inp")
    sme = SME_Struct.load(fname)
    # Build a 2-segment input first, then synthesize both segments.
    w0 = np.array(sme.wave[0], copy=True)
    s0 = np.array(sme.spec[0], copy=True)
    u0 = np.array(sme.uncs[0], copy=True)
    m0 = np.array(sme.mask[0], copy=True)

    sme.wave = [w0.copy(), w0.copy()]
    sme.spec = [s0.copy(), s0.copy()]
    sme.uncs = [u0.copy(), u0.copy()]
    sme.mask = [m0.copy(), m0.copy()]
    sme.wran = [[w0[0], w0[-1]], [w0[0], w0[-1]]]
    sme = synthesize_spectrum(sme)

    rv = 10
    x_syn = sme.wave[0] * (1 - rv / c_light)
    y_syn = sme.synth[0]

    x_syn = np.array([x_syn, x_syn])
    y_syn = np.array([y_syn, y_syn])

    return sme, x_syn, y_syn, rv


@pytest.fixture
def sme_2segments():
    from pysme.sme import SME_Structure as SME_Struct

    cwd = dirname(__file__)

    sme = SME_Struct()
    sme.teff = 5000
    sme.logg = 4.4
    sme.vmic = 1
    sme.vmac = 1
    sme.vsini = 1
    sme.abund = Abund(monh=0, pattern="asplund2009")
    sme.linelist = ValdFile("{}/testcase1.lin".format(cwd))
    sme.atmo.source = "marcs2012p_t2.0.sav"
    sme.atmo.method = "grid"

    sme.wran = [[6550, 6560], [6560, 6574]]

    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    return sme
