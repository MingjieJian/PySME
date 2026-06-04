import importlib
import sys
import types


def test_util_import_does_not_require_sme_synth(monkeypatch):
    monkeypatch.delitem(sys.modules, "pysme.util", raising=False)

    blocker = types.ModuleType("pysme.sme_synth")

    def _block(name):
        raise AssertionError(f"pysme.util imported pysme.sme_synth.{name} too early")

    blocker.__getattr__ = _block
    monkeypatch.setitem(sys.modules, "pysme.sme_synth", blocker)

    util = importlib.import_module("pysme.util")

    assert hasattr(util, "air2vac")
