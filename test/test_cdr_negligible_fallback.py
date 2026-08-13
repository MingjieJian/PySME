# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

from pysme.synthesize import Synthesizer


class _MinimalLineList:
    def __init__(self):
        self._lines = pd.DataFrame(
            {
                "species": ["Fe 1", "Fe 1", "Fe 2"],
                "wlcent": [5000.0, 5001.0, 5002.0],
            }
        )
        self.cdr_paras_thres = {"teff": 250.0, "logg": 0.5, "monh": 0.5, "vmic": 1.0}

    def __len__(self):
        return len(self._lines)

    def __getitem__(self, key):
        return self._lines[key]


def _write_cdr_file(folder, teff, logg, monh, vmic, strong_index=1, width=0.4):
    mask = np.zeros(3, dtype=bool)
    mask[strong_index] = True
    np.savez_compressed(
        folder / f"teff{teff:g}_logg{logg:g}_monh{monh:g}_vmic{vmic:g}.npz",
        mask_bits=np.packbits(mask).astype(np.uint8),
        unique_widths=np.array([width], dtype=np.float32),
        codes=np.zeros(1, dtype=np.uint8),
        n_lines_total=np.int32(mask.size),
    )


def _make_sme():
    class _MinimalSME:
        pass

    sme = _MinimalSME()
    sme.teff = 6212.0
    sme.logg = 3.1
    sme.monh = -4.4
    sme.vmic = 1.5
    sme.linelist = _MinimalLineList()
    return sme


def test_cdr_database_or_mode_falls_back_to_nearest_for_sparse_grid(tmp_path):
    _write_cdr_file(tmp_path, teff=6250.0, logg=3.5, monh=-4.0, vmic=2.0)
    sme = _make_sme()

    Synthesizer().flag_strong_lines_by_database(sme, str(tmp_path))

    lines = sme.linelist._lines
    assert list(lines["strong"]) == [False, True, False]
    assert np.isnan(lines.loc[0, "line_range_s"])
    assert np.isclose(lines.loc[1, "line_range_s"], 5000.8)
    assert np.isclose(lines.loc[1, "line_range_e"], 5001.2)


def test_cdr_database_writes_safe_columns_when_no_grid_point_in_box(tmp_path):
    _write_cdr_file(tmp_path, teff=5000.0, logg=1.0, monh=0.0, vmic=2.0)
    sme = _make_sme()

    Synthesizer().flag_strong_lines_by_database(sme, str(tmp_path))

    lines = sme.linelist._lines
    assert list(lines["strong"]) == [False, False, False]
    assert np.all(np.isnan(lines["line_range_s"]))
    assert np.all(np.isnan(lines["line_range_e"]))
