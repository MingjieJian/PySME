# -*- coding: utf-8 -*-
import logging

import pytest

from pysme.init_config import ensure_user_config
from pysme.large_file_storage import LargeFileStorage
from pysme.smelib import libtools
from pysme.util import warn_with_log


def test_warn_with_log_emits_warning_and_logger(caplog):
    with pytest.warns(UserWarning, match="hello logger"):
        with caplog.at_level(logging.WARNING):
            warn_with_log("hello logger")

    assert any(record.getMessage() == "hello logger" for record in caplog.records)


def test_parse_needed_arch_from_error_logs_instead_of_print(caplog):
    msg = "dlopen failed: need 'arm64'"

    with caplog.at_level(logging.DEBUG, logger=libtools.__name__):
        arch = libtools._parse_needed_arch_from_error(msg)

    assert arch == "arm64"
    messages = [record.getMessage() for record in caplog.records]
    assert any("Incompatible SMElib architecture detected" in msg for msg in messages)
    assert any("Original architecture error" in msg for msg in messages)


def test_ensure_user_config_logs_without_print(monkeypatch, tmp_path, caplog):
    monkeypatch.setenv("HOME", str(tmp_path))

    def fail_print(*args, **kwargs):
        raise AssertionError("ensure_user_config should use logger, not print")

    monkeypatch.setattr("builtins.print", fail_print)

    with caplog.at_level(logging.INFO, logger="pysme.init_config"):
        ensure_user_config()

    messages = [record.getMessage() for record in caplog.records]
    assert any("Creating default PySME user config" in msg for msg in messages)
    assert any("Copying atmosphere datafile references" in msg for msg in messages)
    assert any("Copying NLTE datafile references" in msg for msg in messages)


def test_large_file_storage_logs_directory_creation_without_print(monkeypatch, tmp_path, caplog):
    storage = tmp_path / "cache"

    def fail_print(*args, **kwargs):
        raise AssertionError("LargeFileStorage should use logger, not print")

    monkeypatch.setattr("builtins.print", fail_print)

    with caplog.at_level(logging.INFO, logger="pysme.large_file_storage"):
        LargeFileStorage("https://example.invalid", {}, storage)

    assert storage.exists()
    messages = [record.getMessage() for record in caplog.records]
    assert any("LargeFileStorage cache directory does not exist" in msg for msg in messages)
