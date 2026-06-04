# -*- coding: utf-8 -*-
from pathlib import Path
import hashlib
import json
from urllib.error import URLError

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.config import Config
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import synthesize_spectrum

pytestmark = pytest.mark.filterwarnings(
    "ignore:Covariance of the parameters could not be estimated:scipy.optimize.OptimizeWarning"
)


CASES = {
    "sun_lte": {
        "teff": 5772,
        "logg": 4.44,
        "monh": 0.0,
        "vmic": 1.0,
        "vmac": 2.5,
        "vsini": 1.8,
        "nlte": False,
    },
    "sun_nlte": {
        "teff": 5772,
        "logg": 4.44,
        "monh": 0.0,
        "vmic": 1.0,
        "vmac": 2.5,
        "vsini": 1.8,
        "nlte": True,
    },
    "arcturus_lte": {
        "teff": 4286,
        "logg": 1.64,
        "monh": -0.52,
        "vmic": 1.7,
        "vmac": 3.0,
        "vsini": 1.5,
        "nlte": False,
    },
    "arcturus_nlte": {
        "teff": 4286,
        "logg": 1.64,
        "monh": -0.52,
        "vmic": 1.7,
        "vmac": 3.0,
        "vsini": 1.5,
        "nlte": True,
    },
}

LINELIST = Path(__file__).with_name("halpha_window_cdr_union.lin")
BASELINE = Path(__file__).with_name("halpha_regression.npz")
WRAN = [6561.0, 6564.2]


def _has_local_h_nlte_grid():
    config = Config()
    nlte_root = Path(config["data.nlte_grids"]).expanduser()
    pointers = Path(config["data.pointers.nlte_grids"]).expanduser()
    if not pointers.is_absolute():
        pointers = Path(config.filename).resolve().parent / pointers

    if not pointers.exists():
        return False

    mapping = json.loads(pointers.read_text())
    relpath = mapping.get("nlte_H_pysme.grd")
    if relpath is None:
        return False

    url = f'{config["data.file_server"].rstrip("/")}/{relpath.lstrip("/")}'
    cache_dir = nlte_root / "download" / "url" / hashlib.md5(url.encode()).hexdigest()
    return (cache_dir / "contents").exists()


def _make_structure(case):
    params = CASES[case]
    sme = SME_Structure()
    sme.teff = params["teff"]
    sme.logg = params["logg"]
    sme.monh = params["monh"]
    sme.vmic = params["vmic"]
    sme.vmac = params["vmac"]
    sme.vsini = params["vsini"]
    sme.abund = Abund(monh=params["monh"], pattern="asplund2009")
    sme.linelist = ValdFile(str(LINELIST))
    sme.wran = [WRAN]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.normalize_by_continuum = True
    if params["nlte"]:
        sme.nlte.set_nlte("H", "nlte_H_pysme.grd")
    return sme


def _synthesize_case(case):
    if CASES[case]["nlte"] and not _has_local_h_nlte_grid():
        pytest.skip("H NLTE grid not available in local cache")

    try:
        out = synthesize_spectrum(_make_structure(case))
    except (FileNotFoundError, URLError) as exc:
        pytest.skip(f"Halpha regression data unavailable: {exc}")

    wave = np.asarray(out.wave[0], dtype=float)
    synth = np.asarray(out.synth[0], dtype=float)
    return wave, synth


@pytest.mark.parametrize(
    "case", ["sun_lte", "sun_nlte", "arcturus_lte", "arcturus_nlte"]
)
def test_halpha_regression(case):
    baseline = np.load(BASELINE)
    wave, synth = _synthesize_case(case)

    assert np.allclose(wave, baseline[f"{case}_wave"], rtol=0.0, atol=1e-8)
    assert np.allclose(synth, baseline[f"{case}_synth"], rtol=1e-6, atol=1e-7)
