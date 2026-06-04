# -*- coding: utf-8 -*-
from os.path import dirname, join

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.linelist.linelist import LineList
from pysme.linelist.vald import ValdFile

species = "Fe 1"
wlcent = 5502.9931
excit = 0.9582
gflog = -3.047
gamrad = 7.19
gamqst = -6.22
gamvw = 239.249
linedata = [species, wlcent, excit, gflog, gamrad, gamqst, gamvw]


def test_line_init():
    """Test that property values equal line data passed to __init__()."""
    line = LineList(linedata)
    assert isinstance(line, LineList)
    assert line.species[0] == species
    assert line.wlcent[0] == wlcent
    assert line.excit[0] == excit
    assert line.gflog[0] == gflog
    assert line.gamrad[0] == gamrad
    assert line.gamqst[0] == gamqst
    assert line.gamvw[0] == gamvw


def test_linelist_add_and_len():
    """Test that len() returns the number of lines (including 0) in list."""
    linelist = LineList()
    assert isinstance(linelist, LineList)
    assert len(linelist) == 0
    for iline in range(3):
        print(len(linelist._lines))
        print(linelist._lines)
        assert len(linelist) == iline
        linelist.add(*linedata)


def test_linelist_properties():
    """Test that properties are lists with one item per line.
    Test that property value equal line data passed to add().
    """
    linelist = LineList()
    linelist.add(*linedata)
    proplist = [
        linelist.species,
        linelist.wlcent,
        linelist.excit,
        linelist.gflog,
        linelist.gamrad,
        linelist.gamqst,
        linelist.gamvw,
    ]
    for iprop, prop in enumerate(proplist):
        print(linelist._lines)
        assert isinstance(prop, np.ndarray)
        assert len(prop) == 1
        assert prop[0] == linedata[iprop]


def test_valdfile():
    """Test class to read a VALD line data file."""
    testdir = dirname(__file__)
    linelist = ValdFile(join(testdir, "testcase1.lin"))

    assert len(linelist) == 44
    assert linelist.lineformat == "short"
    assert linelist.medium == "air"
    assert linelist[0].species == "V 1"

    with pytest.raises(IOError):
        ValdFile(join(testdir, "testcase1.npy"))


def test_medium():
    """Test class to read a VALD line data file."""
    testdir = dirname(__file__)
    vf = ValdFile(join(testdir, "testcase1.lin"), medium="vac")

    assert vf.medium == "vac"


def test_short_format():
    linelist = ValdFile(join(dirname(__file__), "testcase1.lin"))

    assert linelist.lineformat == "short"
    assert len(linelist) == 44
    assert linelist.atomic is not None
    assert linelist.species is not None

    with pytest.raises(AttributeError):
        _ = linelist.lulande

    with pytest.raises(AttributeError):
        _ = linelist.extra

    assert isinstance(linelist.abund, Abund)
    assert isinstance(linelist.atmo, str)


def test_long_format():
    linelist = ValdFile(join(dirname(__file__), "testcase3.lin"))

    assert linelist.lineformat == "long"
    assert len(linelist) == 67
    assert linelist.atomic is not None
    assert linelist.species is not None
    assert linelist.lulande is not None
    assert linelist.extra is not None

    assert linelist.reference is not None
    assert linelist.error is not None
    assert linelist.term_lower is not None
    assert linelist.term_upper is not None

    assert isinstance(linelist.abund, Abund)
    assert isinstance(linelist.atmo, str)


def test_save_reload_long_subset_updates_header(tmp_path):
    linelist = ValdFile(join(dirname(__file__), "testcase3.lin"))
    subset = linelist[:5]

    path = tmp_path / "subset_long.lin"
    subset.save(path, overwrite=True)
    assert path.read_text().splitlines()[0].split(",")[2].strip() == "5"
    assert "species_len" not in subset._lines.columns
    reloaded = ValdFile(path)

    assert len(reloaded) == len(subset)
    assert reloaded.lineformat == "long"
    assert np.allclose(reloaded.wlcent, subset.wlcent)


def test_save_reload_short_subset_updates_header(tmp_path):
    linelist = ValdFile(join(dirname(__file__), "testcase1.lin"))
    subset = linelist[:5]

    path = tmp_path / "subset_short.lin"
    subset.save(path, overwrite=True)
    assert path.read_text().splitlines()[0].split(",")[2].strip() == "5"
    assert "species_len" not in subset._lines.columns
    reloaded = ValdFile(path)

    assert len(reloaded) == len(subset)
    assert reloaded.lineformat == "short"
    assert np.allclose(reloaded.wlcent, subset.wlcent)
