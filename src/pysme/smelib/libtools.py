# -*- coding: utf-8 -*-
"""
Utilities for locating the SME C/Fortran library and data files.
"""

import ctypes as ct
import os
import platform
from os.path import dirname, join


def get_lib_name():
    """Get the name of the SME C library"""
    system = platform.system().lower()

    if system == "windows":
        return "libsme-5.dll"
    elif system == "darwin":
        return "libsme.dylib"

    arch = platform.machine()
    bits = 64

    return "sme_synth.so.{system}.{arch}.{bits}".format(
        system=system, arch=arch, bits=bits
    )


def get_lib_directory():
    """
    Get the directory name of the library. I.e. 'lib' on all systems
    except windows, and 'bin' on windows
    """
    if platform.system() in ["Windows"]:
        dirpath = "bin"
    else:
        dirpath = "lib"
    return dirpath


def get_full_libfile():
    """Get the full path to the sme C library"""
    localdir = dirname(__file__)
    libfile = get_lib_name()
    dirpath = get_lib_directory()
    libfile = join(localdir, dirpath, libfile)
    return libfile


def load_library(libfile=None):
    """
    Load the SME library using ctypes.CDLL

    Parameters
    ----------
    libfile : str, optional
        filename of the library to load, by default use the SME library in
        this package

    Returns
    -------
    lib : CDLL
        library object of the SME library
    """
    if libfile is None:
        libfile = get_full_libfile()
    try:
        os.add_dll_directory(dirname(libfile))
    except AttributeError:
        newpath = dirname(libfile)
        if "PATH" in os.environ:
            newpath += os.pathsep + os.environ["PATH"]
        os.environ["PATH"] = newpath
    return ct.CDLL(str(libfile))


def get_full_datadir():
    """
    Get the filepath to the datafiles of the SME library
    """
    localdir = dirname(__file__)
    datadir = join(localdir, "data/")
    return datadir
