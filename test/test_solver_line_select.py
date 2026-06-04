# -*- coding: utf-8 -*-
from os.path import dirname, join

import numpy as np
import pytest

import pysme.solve as solve_mod
import pysme.synthesize as synth_mod
from pysme.sme import SME_Structure as SME_Struct
from pysme.solve import solve
from pysme.synthesize import Synthesizer, synthesize_spectrum


cwd = dirname(__file__)
filename = join(cwd, "testcase1.inp")


def _concat_segments(vec, nseg):
    return np.concatenate([np.asarray(vec[i], dtype=float) for i in range(nseg)])


def _prepare_sme():
    sme = SME_Struct.load(filename)
    wran = np.asarray(sme.wran, dtype=float).reshape(-1, 2)
    wl = np.asarray(sme.linelist["wlcent"], dtype=float)
    lo = float(np.min(wran[:, 0]) - 1.0)
    hi = float(np.max(wran[:, 1]) + 1.0)
    keep = (wl >= lo) & (wl <= hi)
    if np.any(keep):
        sme.linelist = sme.linelist[keep]
    return sme


def _run_solve(method):
    sme = _prepare_sme()
    sme.line_select_method = method
    sme.line_select_policy = "strict" if method in ("almax", "cdr") else "auto"
    sme.line_select_recompute = "always"
    sme.line_select_reuse = "none"
    sme.line_select_parallel = False

    # Use only vrad to keep runtime low while still exercising the solve path.
    out = solve(sme, ["vrad"], linelist_mode="all")
    wave = _concat_segments(out.wave, out.nseg)
    flux = _concat_segments(out.synth, out.nseg)
    return out, wave, flux


@pytest.mark.parametrize(
    "method,required_cols",
    [
        ("internal", []),
        ("almax", ["almax_ratio", "line_range_s", "line_range_e", "strong"]),
        ("cdr", ["central_depth", "line_range_s", "line_range_e", "strong"]),
    ],
)
def test_solve_line_select_methods_run(method, required_cols):
    out, _, flux = _run_solve(method)
    cols = set(out.linelist._lines.columns)

    assert out.synth is not None
    assert np.all(np.isfinite(flux))
    for col in required_cols:
        assert col in cols


def test_solve_line_select_methods_match_internal_flux():
    _, wave_ref, flux_ref = _run_solve("internal")

    for method in ("almax", "cdr"):
        _, wave, flux = _run_solve(method)
        flux_interp = np.interp(wave_ref, wave, flux)
        diff = flux_interp - flux_ref
        assert np.max(np.abs(diff)) < 1e-10


def test_line_precompute_database_separates_method_and_linelist_hash(tmp_path):
    db = str(tmp_path)

    sme_cdr = _prepare_sme()
    sme_cdr.line_select_method = "cdr"
    sme_cdr.line_select_policy = "strict"
    sme_cdr.line_select_recompute = "always"
    sme_cdr.line_select_parallel = False
    _ = synthesize_spectrum(
        sme_cdr,
        linelist_mode="all",
        line_precompute_database=db,
    )

    cdr_files = list(tmp_path.glob("cdr_*.npz"))
    assert len(cdr_files) > 0

    sme_almax = _prepare_sme()
    sme_almax.line_select_method = "almax"
    sme_almax.line_select_policy = "strict"
    sme_almax.line_select_recompute = "always"
    sme_almax.line_select_parallel = False
    sme_almax.line_select_almax_threshold = sme_almax.accrt
    _ = synthesize_spectrum(
        sme_almax,
        linelist_mode="all",
        line_precompute_database=db,
    )

    almax_files = list(tmp_path.glob("almax_*.npz"))
    assert len(almax_files) > 0

    # Mismatched linelist (different order/content) should be treated as cache miss.
    sme_mismatch = _prepare_sme()
    sme_mismatch.linelist = sme_mismatch.linelist[1:]
    sme_mismatch.line_select_method = "cdr"
    sme_mismatch.line_select_policy = "strict"
    sme_mismatch.line_select_recompute = "never"
    sme_mismatch.line_select_parallel = False
    with pytest.raises(ValueError, match="no matching entry"):
        _ = synthesize_spectrum(
            sme_mismatch,
            linelist_mode="all",
            line_precompute_database=db,
        )


class _StopLineSelect(RuntimeError):
    pass


class _FakeDll:
    def SetLineInfoMode(self, mode):
        self.mode = mode

    def SetLibraryPath(self):
        pass

    def InputLineList(self, linelist):
        raise _StopLineSelect


def _make_minimal_cdr_ready_linelist(sme):
    nlines = len(sme.linelist)
    sme.linelist._lines["central_depth"] = np.full(nlines, 0.02)
    wl = np.asarray(sme.linelist["wlcent"], dtype=float)
    sme.linelist._lines["line_range_s"] = wl - 0.05
    sme.linelist._lines["line_range_e"] = wl + 0.05
    sme.linelist._lines["strong"] = np.ones(nlines, dtype=bool)
    sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic], dtype=float)
    sme.linelist.cdr_paras_thres["strong_depth"] = 0.001
    sme.linelist.cdr_paras_thres["strong_bin_width"] = 0.2


def _make_minimal_synth():
    synth = Synthesizer.__new__(Synthesizer)
    synth.dll = _FakeDll()
    synth.wint = {}
    synth.known_sme = None
    synth.update_cdr_switch = False
    return synth


def test_cdr_update_uses_resolved_line_select_config(monkeypatch):
    sme = _prepare_sme()
    sme.line_select_method = "cdr"
    sme.line_select_policy = "strict"
    sme.line_select_recompute = "always"
    sme.line_select_parallel = False
    sme.line_select_n_jobs = 3
    sme.line_select_chunk_size = 17
    sme.line_select_cdr_strength_thres = 0.0123
    sme.line_select_cdr_bin_width = 0.45
    sme.cdr_parallel = True
    sme.cdr_n_jobs = 99
    sme.cdr_N_line_chunk = 2
    sme.strong_depth_thres = 0.99
    sme.strong_bin_width = 0.88

    captured = {}

    def fake_update_cdr(self, sme, **kwargs):
        captured.update(kwargs)
        raise _StopLineSelect

    monkeypatch.setattr(Synthesizer, "update_cdr", fake_update_cdr)
    synth = Synthesizer.__new__(Synthesizer)
    synth.dll = _FakeDll()
    synth.wint = {}
    synth.known_sme = None
    synth.update_cdr_switch = False

    with pytest.raises(_StopLineSelect):
        synth.synthesize_spectrum(sme, linelist_mode="all")

    assert captured["chunk_size"] == 17
    assert captured["parallel"] is False
    assert captured["n_jobs"] == 3
    assert np.isclose(sme.strong_depth_thres, 0.0123)
    assert np.isclose(sme.strong_bin_width, 0.45)
    assert sme.cdr_N_line_chunk == 17
    assert sme.cdr_n_jobs == 3
    assert sme.cdr_parallel is False


def test_almax_update_uses_resolved_line_select_config(monkeypatch):
    sme = _prepare_sme()
    sme.line_select_method = "almax"
    sme.line_select_policy = "strict"
    sme.line_select_recompute = "always"
    sme.line_select_parallel = False
    sme.line_select_n_jobs = 4
    sme.line_select_chunk_size = 19
    sme.line_select_almax_threshold = 0.0042
    sme.line_select_almax_use_bins = True
    sme.line_select_almax_bin_width = 0.31
    sme.cdr_parallel = True
    sme.cdr_n_jobs = 77
    sme.cdr_N_line_chunk = 5
    sme.strong_bin_width = 0.91

    captured = {}

    def fake_update_almax(self, sme, **kwargs):
        captured.update(kwargs)
        raise _StopLineSelect

    monkeypatch.setattr(Synthesizer, "update_almax", fake_update_almax)
    synth = Synthesizer.__new__(Synthesizer)
    synth.dll = _FakeDll()
    synth.wint = {}
    synth.known_sme = None
    synth.update_cdr_switch = False

    with pytest.raises(_StopLineSelect):
        synth.synthesize_spectrum(sme, linelist_mode="all")

    assert captured["chunk_size"] == 19
    assert captured["parallel"] is False
    assert captured["n_jobs"] == 4
    assert np.isclose(captured["threshold"], 0.0042)
    assert captured["use_bins"] is True
    assert np.isclose(captured["bin_width"], 0.31)


def test_line_select_reuse_is_deprecated(monkeypatch):
    sme = _prepare_sme()
    sme.line_select_method = "cdr"
    sme.line_select_recompute = "always"
    sme.line_select_reuse = "once"

    def fake_update_cdr(self, sme, **kwargs):
        raise _StopLineSelect

    monkeypatch.setattr(Synthesizer, "update_cdr", fake_update_cdr)
    synth = Synthesizer.__new__(Synthesizer)
    synth.dll = _FakeDll()
    synth.wint = {}
    synth.known_sme = None
    synth.update_cdr_switch = False

    with pytest.deprecated_call(match="line_select_reuse"):
        with pytest.raises(_StopLineSelect):
            synth.synthesize_spectrum(sme, linelist_mode="all")


def test_jacobian_scale_parameter_shift_within_stale_threshold_does_not_recompute_cdr(monkeypatch):
    sme = _prepare_sme()
    sme.line_select_method = "cdr"
    sme.line_select_policy = "strict"
    sme.line_select_recompute = "if_stale"
    sme.line_select_parallel = False
    sme.line_select_cdr_strength_thres = 0.001
    sme.line_select_cdr_bin_width = 0.2
    _make_minimal_cdr_ready_linelist(sme)

    # Small Jacobian-like perturbation: well within default stale thresholds.
    sme.teff += 10.0
    sme.logg += 0.01
    sme.monh += 0.01
    sme.vmic += 0.05

    called = {"count": 0}

    def fake_update_cdr(self, sme, **kwargs):
        called["count"] += 1
        raise AssertionError("update_cdr should not be called for in-threshold Jacobian perturbations")

    monkeypatch.setattr(Synthesizer, "update_cdr", fake_update_cdr)
    synth = _make_minimal_synth()

    with pytest.raises(_StopLineSelect):
        synth.synthesize_spectrum(sme, linelist_mode="all", updateStructure=False)

    assert called["count"] == 0


def test_jacobian_scale_parameter_shift_beyond_stale_threshold_recomputes_cdr(monkeypatch):
    sme = _prepare_sme()
    sme.line_select_method = "cdr"
    sme.line_select_policy = "strict"
    sme.line_select_recompute = "if_stale"
    sme.line_select_parallel = False
    sme.line_select_cdr_strength_thres = 0.001
    sme.line_select_cdr_bin_width = 0.2
    _make_minimal_cdr_ready_linelist(sme)

    # Large perturbation: exceeds default Teff stale threshold of 250 K.
    sme.teff += 300.0

    called = {"count": 0}

    def fake_update_cdr(self, sme, **kwargs):
        called["count"] += 1
        raise _StopLineSelect

    monkeypatch.setattr(Synthesizer, "update_cdr", fake_update_cdr)
    synth = _make_minimal_synth()

    with pytest.raises(_StopLineSelect):
        synth.synthesize_spectrum(sme, linelist_mode="all", updateStructure=False)

    assert called["count"] == 1


def test_synthesize_aliases_cdr_database_to_line_precompute_database(monkeypatch):
    sme = _prepare_sme()
    captured = {}

    class _FakeSynthesizer:
        def __init__(self):
            pass

        def synthesize_spectrum(self, sme, segments="all", **kwargs):
            captured.update(kwargs)
            return sme

    monkeypatch.setattr(synth_mod, "Synthesizer", _FakeSynthesizer)

    with pytest.deprecated_call(match="cdr_database"):
        synthesize_spectrum(sme, cdr_database="/tmp/cdr-cache")

    assert captured["line_precompute_database"] == "/tmp/cdr-cache"
    assert "cdr_database" not in captured or captured["cdr_database"] is None


def test_solve_aliases_cdr_database_to_line_precompute_database(monkeypatch):
    sme = _prepare_sme()
    captured = {}

    class _FakeSolver:
        def __init__(self, filename=None, restore=False):
            self.filename = filename
            self.restore = restore

        def solve(self, sme, param_names=None, segments="all", **kwargs):
            captured.update(kwargs)
            return sme

    monkeypatch.setattr(solve_mod, "SME_Solver", _FakeSolver)

    with pytest.deprecated_call(match="cdr_database"):
        solve(sme, ["vrad"], cdr_database="/tmp/cdr-cache")

    assert captured["line_precompute_database"] == "/tmp/cdr-cache"
    assert "cdr_database" not in captured or captured["cdr_database"] is None
