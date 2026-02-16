# -*- coding: utf-8 -*-
import pytest
import requests

from pysme.config import Config
from pysme.large_file_storage import setup_atmo, setup_nlte


def lfs_available():
    config = Config()
    r = requests.head(config["data.file_server"])
    return r.status_code == 200


skipif_lfs = pytest.mark.skipif(not lfs_available(), reason="LFS not available")


@pytest.fixture
def lfs_nlte():
    lfs_nlte = setup_nlte()
    yield lfs_nlte


@pytest.fixture
def lfs_atmo():
    lfs_atmo = setup_atmo()
    yield lfs_atmo
