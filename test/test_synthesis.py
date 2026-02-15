# -*- coding: utf-8 -*-
# TODO implement synthesis tests
import numpy as np
import pytest

from pysme.iliffe_vector import Iliffe_vector
from pysme.sme import SME_Structure as SME_Struct
from pysme.synthesize import Synthesizer, synthesize_spectrum


def test_synthesis_simple(sme_2segments):
    sme = sme_2segments
    sme2 = synthesize_spectrum(sme)

    # Check if a result is there it has the expected data type
    assert sme2.synth is not None
    assert isinstance(sme2.synth, Iliffe_vector)
    assert np.all(sme2.synth.ravel() != 0)

    assert sme.wave is not None
    assert isinstance(sme2.wave, Iliffe_vector)
    assert np.issubdtype(sme2.wave.dtype, np.floating)
    assert np.all(sme2.wave.ravel() != 0)

    assert sme.spec is None


def test_synthesis_segment(sme_2segments):
    sme = sme_2segments
    # Out of range
    with pytest.raises(IndexError):
        synthesize_spectrum(sme, segments=[3])

    with pytest.raises(IndexError):
        synthesize_spectrum(sme, segments=[-1])

    sme2 = synthesize_spectrum(sme, segments=[0])
    assert len(sme2.synth[0]) != 0
    assert len(sme2.synth[1]) == 0

    assert len(sme2.wave[0]) != 0
    assert len(sme2.wave[1]) == 0

    assert sme2.wave.shape[0] == 2
    assert sme2.wave.shape[1][1] == 0

    orig = np.copy(sme2.synth[0])
    sme2 = synthesize_spectrum(sme2, segments=[1])

    assert sme2.wave.shape[1][0] != 0
    assert sme2.wave.shape[1][1] != 0

    assert np.all(sme2.synth[0] == orig)


class _DummyDLL:
    def __init__(self):
        self.last_wave = "unset"

    def InputWaveRange(self, *_):
        return None

    def Opacity(self):
        return None

    def Transf(self, mu, accrt, accwi, keep_lineop, wave=None):
        self.last_wave = wave
        if wave is None:
            wint = np.linspace(5000.0, 5001.0, 5)
        else:
            wint = np.asarray(wave, dtype=float)
        sint = np.ones((len(mu), len(wint)), dtype=float)
        cint = np.ones((len(mu), len(wint)), dtype=float)
        return len(wint), wint, sint, cint

    def CentralDepth(self, mu, accrt):
        return np.zeros(0, dtype=float)

    def GetLineRange(self):
        return np.zeros((0, 2), dtype=float)


def _minimal_sme():
    sme = SME_Struct()
    sme.wran = [[5000.0, 5001.0]]
    sme.vsini = 0.0
    sme.vmac = 0.0
    sme.vrad_flag = "none"
    sme.specific_intensities_only = True
    sme.normalize_by_continuum = True
    return sme


def test_synthesize_segment_prefers_user_wint_over_cache():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    user_wint = np.linspace(5000.0, 5001.0, 7)
    cached_wint = np.linspace(5000.0, 5001.0, 9)
    sme.wint = user_wint
    synth.wint[0] = cached_wint.copy()

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert np.allclose(dll.last_wave, user_wint)
    assert np.allclose(synth.wint[0], cached_wint)


def test_synthesize_segment_uses_cache_when_user_wint_missing():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    cached_wint = np.linspace(5000.0, 5001.0, 9)
    synth.wint[0] = cached_wint

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert np.allclose(dll.last_wave, cached_wint)


def test_synthesize_segment_populates_cache_when_no_wint_available():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert dll.last_wave is None
    assert 0 in synth.wint
    assert np.allclose(synth.wint[0], np.linspace(5000.0, 5001.0, 5))
