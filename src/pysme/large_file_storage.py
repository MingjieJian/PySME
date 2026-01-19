# -*- coding: utf-8 -*-
"""
System to store large data files on a server
Load them whem required by the user
Update the pointer file on github when new datafiles become available

Pro: Versioning is effectively done by Git
Con: Need to run server
"""

import gzip
import json
import logging
import os
from os.path import basename
from pathlib import Path
from tempfile import NamedTemporaryFile

from astropy.utils.data import (
    clear_download_cache,
    download_file,
    import_file_to_cache,
)
from tqdm.auto import tqdm
from tqdm.utils import CallbackIOWrapper

from .config import Config
from .util import show_progress_bars

logger = logging.getLogger(__name__)

# We are lazy and want a simple check if a file is in the Path
Path.__contains__ = lambda self, key: (self / key).exists()


class LargeFileStorage:
    """
    Download large data files from data server when needed
    New versions of the datafiles are indicated in a 'pointer' file
    that includes the hash of the newest version of the files

    Raises
    ------
    FileNotFoundError
        If the datafiles can't be located anywhere
    """

    def __init__(self, server, pointers, storage):
        #:Server: Large File Storage Server address
        self.server = server

        if isinstance(pointers, str):
            path = Path(__file__).parent / pointers
            pointers = LargeFileStorage.load_pointers_file(path)

        #:dict(fname:hash): points from a filename to the current newest object id, usually a hash
        self.pointers = pointers
        #:Directory: directory of the current data files
        self.current = Path(storage).expanduser().absolute()

        # set the folder to download the data file into
        # need to set environment variable because astropy will put things into home otherwise
        os.environ['XDG_CACHE_HOME'] = str(Path(storage).absolute())
        # if someone is using astropy along with pysme, it might mess with their astropy file storage
        # not threadsafe, but multiprocessing safe, because threads shares environment variables
        self.PKGNAME = ''
        
        if not os.path.exists(storage):
            print('folder to store data file does not exist, creating')
        cache_path = Path(storage)
        cache_path.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def load_pointers_file(filename):
        try:
            with open(str(filename), "r") as f:
                pointers = json.load(f)
        except FileNotFoundError:
            logger.error("Could not find LargeFileStorage reference file %s", filename)
            pointers = {}
        return pointers

    def get(self, key):
        """
        Request a datafile from the LargeFileStorage
        Assures that tracked files are at the specified version
        And downloads data from the server if necessary

        Parameters
        ----------
        key : str
            Name of the requested datafile

        Raises
        ------
        FileNotFoundError
            If the requested datafile can not be found anywhere

        Returns
        -------
        fullpath : str
            Absolute path to the datafile
        """
        url = self.get_url(key)
        # If its a direct file link, pass that directly to
        if url.startswith("file://"):
            return url[7:]

        fname = download_file(url, cache=True, pkgname=self.PKGNAME)

        compression = self._test_compression(fname)
        if compression is None:
            pass
        elif compression == "gzip":
            # If the file is compressed
            # Replace the cache file with the decompressed file
            self._unpack_gzip(fname, key, url)
        else:
            raise ValueError(
                "The file is compressed using %s, which is not supported" % compression
            )

        return fname

    def _test_compression(self, fname):
        """Check filetype using the magic string"""
        with open(fname, "rb") as f:
            magic = f.read(6)
            if magic[:2] == b"\x1f\x8b":
                return "gzip"
            if magic[:6] == b"\x37\x7A\xBC\xAF\x27\x1C":
                return "7z"
            if magic[:5] == b"\x50\x4B\x03\x04":
                return "zip"
            return None

    def _unpack_gzip(self, fname, key, url):
        logger.debug("Unpacking data file %s", key)

        # We have to use a try except block, as this will crash with
        # permissions denied on windows, when trying to copy an open file
        # here the temporary file
        # Therefore we close the file, after copying and then delete it manually
        try:
            with gzip.open(fname, "rb") as f_in:
                with NamedTemporaryFile("wb", delete=False) as f_out:
                    with tqdm(
                        # total=f_in.size,
                        desc="Unpack",
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        disable=~show_progress_bars,
                    ) as t:
                        fobj = CallbackIOWrapper(t.update, f_in, "read")
                        while True:
                            chunk = fobj.read(1024)
                            if not chunk:
                                break
                            f_out.write(chunk)
                        f_out.flush()
                        t.reset()
            import_file_to_cache(url, f_out.name, pkgname=self.PKGNAME)
        finally:
            try:
                os.remove(f_out.name)
            except:
                pass

    def get_url(self, key):
        key = str(key)

        # Check if the file is tracked and/or exists in the storage directory
        if key not in self.pointers:
            if key not in self.current:
                if not os.path.exists(key):
                    raise FileNotFoundError(
                        f"File {key} does not exist and is not tracked by the Large File system"
                    )
                else:
                    return Path(key).as_uri()
            else:
                return (self.current / key).as_uri()

        # Otherwise get it from the cache or online if necessary
        newest = self.pointers[key]
        url = self.server + "/" + newest
        return url

    def clean_cache(self):
        """Remove unused cache files (from old versions)"""
        clear_download_cache(pkgname=self.PKGNAME)

    def delete_file(self, fname):
        """Delete a file, including the cache file"""
        clear_download_cache(fname, pkgname=self.PKGNAME)

    def move_to_cache(self, fname, key=None):
        """Move currently used files into cache directory and use symlinks instead,
        just as if downloaded from a server"""
        if key is None:
            key = basename(fname)
        import_file_to_cache(key, fname, pkgname=self.PKGNAME)
        self.pointers[key] = key


def setup_atmo(config=None):
    if config is None:
        config = Config()
    server = config["data.file_server"]
    storage = config["data.atmospheres"]
    pointers = config["data.pointers.atmospheres"]
    lfs_atmo = LargeFileStorage(server, pointers, storage)
    return lfs_atmo


def setup_nlte(config=None):
    if config is None:
        config = Config()
    server = config["data.file_server"]
    storage = config["data.nlte_grids"]
    pointers = config["data.pointers.nlte_grids"]
    lfs_nlte = LargeFileStorage(server, pointers, storage)
    return lfs_nlte


def setup_lfs(config=None, lfs_atmo=None, lfs_nlte=None):
    if config is None:
        config = Config()
    if lfs_atmo is None:
        lfs_atmo = setup_atmo(config)
    if lfs_nlte is None:
        lfs_nlte = setup_nlte(config)
    return config, lfs_atmo, lfs_nlte


def get_available_atmospheres(config=None):
    if config is None:
        config = Config()
    pointers = config["data.pointers.atmospheres"]
    storage = config["data.atmospheres"]
    data = get_available_files(pointers, storage)
    return data


def get_available_nlte_grids(config=None):
    if config is None:
        config = Config()
    pointers = config["data.pointers.nlte_grids"]
    storage = config["data.nlte_grids"]
    data = get_available_files(pointers, storage)
    return data


def get_available_files(pointers, storage):
    pointers = Path(__file__).parent / pointers
    storage = Path(storage).expanduser()
    data = LargeFileStorage.load_pointers_file(pointers)
    files = list(data.keys())
    files_non_lfs = [
        f
        for f in os.listdir(storage)
        if f not in data and not os.path.isdir(storage / f)
    ]
    files += files_non_lfs
    return files
