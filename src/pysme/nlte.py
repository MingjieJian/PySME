# -*- coding: utf-8 -*-
"""
NLTE module of SME
reads and interpolates departure coefficients from library files
"""

import itertools
import logging
import warnings

import numpy as np
from scipy import interpolate
from tqdm import tqdm

from .abund import Abund
from .abund import elements as abund_elem
from .data_structure import Collection, CollectionFactory, array, astype, oneof, this
from .util import show_progress_bars

# from memory_profiler import profile

logger = logging.getLogger(__name__)

class DirectAccessFile:
    """
    This function reads a single record from binary file that has the following
    structure:
    Version string              - 64 byte long
    # of directory bocks        - short int
    directory block length      - short int
    # of used directory blocks  - short int
    1st directory block
    key      - string of up to 256 characters padded with ' '
    datatype - 32-bit int 23-element array returned by SIZE
    pointer  - 64-bit int pointer to the beginning of the record

    2nd directory block
    ...
    last directory block
    1st data record
    ...
    last data record
    """

    def __init__(self, filename):
        key, pointer, dtype, shape, version = DirectAccessFile.read_header(filename)
        self.file = filename
        self.version = version
        self.key = key
        self.shape = shape
        self.pointer = pointer
        self.dtype = dtype

    def __getitem__(self, key):
        # Access data via brackets
        value = self.get(key)
        if value is None:
            raise KeyError(f"Key {key} not found")
        return value

    def get(self, key: str, alt=None) -> np.memmap:
        """get field from file"""
        idx = np.where(self.key == key)[0]

        if idx.size == 0:
            return alt
        else:
            idx = idx[0]

        return np.memmap(
            self.file,
            mode="r",
            offset=self.pointer[idx],
            dtype=self.dtype[idx],
            shape=self.shape[idx][::-1],
        )

    @staticmethod
    def idl_typecode(i):
        """
        relevant IDL typecodes and corresponding Numpy Codes
        Most specify a byte size, but not all
        """
        typecode = {
            0: "V",
            1: "B",
            2: "i2",
            3: "i4",
            4: "f4",
            5: "f8",
            6: "c4",
            7: "S",
            8: "O",
            9: "c8",
            10: "i8",
            11: "O",
            12: "u2",
            13: "u4",
            14: "i8",
            15: "u8",
        }
        return typecode[i]

    @staticmethod
    def get_typecode(dt):
        """
        relevant IDL typecodes and corresponding Numpy Codes
        Most specify a byte size, but not all
        """
        typecode = {
            "V": 0,
            "B": 1,
            "U": 1,
            "i2": 2,
            "i4": 3,
            "f4": 4,
            "f8": 5,
            "c4": 6,
            "S": 7,
            "O": 8,
            "c8": 9,
            "i8": 10,
            "O": 11,
            "u2": 12,
            "u4": 13,
            "i8": 14,
            "u8": 15,
        }
        if isinstance(dt, np.dtype):
            dt = dt.str
        if len(dt) > 2:
            dt = dt[1:]

        return typecode[dt]

    @staticmethod
    def get_dtypes(version):
        major, minor = version
        if major == 1 and minor == 00:
            header_dtype = np.dtype(
                [("nblocks", "<u2"), ("dir_length", "<u2"), ("ndir", "<u2")]
            )
            dir_dtype = np.dtype(
                [("key", "S256"), ("size", "<i4", 23), ("pointer", "<i8")]
            )
        elif major == 1 and minor >= 10:
            header_dtype = np.dtype(
                [("nblocks", "<u8"), ("dir_length", "<i2"), ("ndir", "<u8")]
            )
            dir_dtype = np.dtype(
                [("key", "S256"), ("size", "<i4", 23), ("pointer", "<i8")]
            )
        else:
            raise ValueError("DirectAccess File Version '{version}' not understood.")
        return header_dtype, dir_dtype

    @staticmethod
    def read_version(version_string: str):
        if isinstance(version_string, bytes):
            version_string = version_string.decode()

        major, minor = int(version_string[26]), int(version_string[28:30])
        version = (major, minor)
        return version

    @classmethod
    def read_header(cls, fname: str):
        """parse Header data"""
        with open(fname, "rb") as file:
            version_dtype = "S64"
            version = np.fromfile(file, version_dtype, count=1)
            version = cls.read_version(version[0])
            header_dtype, dir_dtype = cls.get_dtypes(version)
            header = np.fromfile(file, header_dtype, count=1)
            ndir = int(header["ndir"][0])
            directory = np.fromfile(file, dir_dtype, count=ndir)

        # Decode bytes to strings
        key = directory["key"]
        key = np.char.strip(np.char.decode(key))
        # Get relevant info from size parameter
        # ndim, n1, n2, ..., typecode, size
        dtype = np.array(
            [cls.idl_typecode(d[1 + d[0]]) for d in directory["size"]],
            dtype="U5",
        )
        shape = np.empty(ndir, dtype=object)
        shape[:] = [tuple(d[1 : d[0] + 1]) for d in directory["size"]]
        # Pointer to data arrays
        pointer = directory["pointer"]

        # Bytes (which represent strings) get special treatment to get the dimensions right
        # And to properly convert them to strings
        # Also Null terminator is important to remember
        idx = dtype == "B"
        dtype[idx] = [f"S{s[0]}" for s in shape[idx]]
        shape2 = np.empty(np.count_nonzero(idx), dtype=object)
        shape2[:] = [tuple(s[1:]) for s in shape[idx]]
        shape[idx] = shape2

        return key, pointer, dtype, shape, version

    @classmethod
    def write(cls, fname, **kwargs):
        major, minor = 1, 10
        ndir = len(kwargs)
        version = f"DirectAccess file Version {major}.{minor} 2011-03-24"
        version = np.asarray([version], dtype="S64")
        version_length = version.itemsize

        header_dtype, dir_dtype = cls.get_dtypes((major, minor))
        header_length = header_dtype.itemsize
        dir_length = dir_dtype.itemsize

        header = np.zeros(1, header_dtype)
        header[0]["nblocks"] = len(kwargs)
        header[0]["dir_length"] = dir_length
        header[0]["ndir"] = ndir

        directory = np.zeros(ndir, dir_dtype)
        pointer = version_length + header_length + ndir * dir_length

        for i, (key, value) in enumerate(kwargs.items()):
            value = np.asarray(value)

            directory[i]["key"] = key
            shape = directory[i]["size"]
            if np.issubdtype(value.dtype, np.dtype("U")) or np.issubdtype(
                value.dtype, np.dtype("S")
            ):
                value = value.astype("S")
                shape[0] = value.ndim + 1
                shape[1 : 2 + value.ndim] = value.itemsize, *value.shape
                shape[2 + value.ndim] = 1
                shape[3 + value.ndim] = value.size * value.itemsize
            else:
                shape[0] = value.ndim
                shape[1 : 1 + value.ndim] = value.shape
                shape[1 + value.ndim] = cls.get_typecode(value.dtype)
                shape[2 + value.ndim] = value.size

            directory[i]["pointer"] = pointer
            pointer += value.nbytes
            kwargs[key] = value

        with open(fname, "wb") as file:
            version.tofile(file)
            header.tofile(file)
            directory.tofile(file)

            for i, value in enumerate(kwargs.values()):
                value.tofile(file)


class Grid:
    """NLTE Grid class that handles all NLTE data reading and interpolation"""

    # @profile
    def __init__(
        self,
        sme,
        elem,
        lfs_nlte,
        selection="energy",
        solar=None,
        abund_format="Fe=12",
        min_energy_diff=None,
    ):
        #:str: Element of the NLTE grid
        self.elem = elem

        if 'use_indices' in sme.linelist._lines.columns:
            #:LineList: Whole LineList that was passed to the C library
            self.linelist = sme.linelist[sme.linelist['use_indices']]
            #:array(str): Elemental Species Names for the linelist
            self.species = sme.linelist.species[sme.linelist['use_indices']]
            # sme.linelist._lines = sme.linelist._lines.drop(columns=['use_indices'])
        else:
            #:LineList: Whole LineList that was passed to the C library
            self.linelist = sme.linelist
            #:array(str): Elemental Species Names for the linelist
            self.species = sme.linelist.species

        #:str: Name of the grid
        self.grid_name = sme.nlte.grids[elem]
        #:str: complete filename of the NLTE grid data file
        self.fname = lfs_nlte.get(self.grid_name)
        #:{"levels", "energy"}: Selection algorithm to match lines in grid with linelist
        self.selection = selection
        #:float: Minimum difference between energy levels to match, only relevant for selection=='energy'
        self.min_energy_diff = min_energy_diff

        #:DirectAccessFile: The NLTE data file
        self.directory = DirectAccessFile(self.fname)
        self.version = self.directory.version

        # Check atmosphere compatibility
        if "atmosphere_grid" in self.directory.key:
            self.atmosphere_grid = self.directory["atmosphere_grid"][0]
            if self.atmosphere_grid != sme.atmo.source:
                logger.warning(
                    "The {elem} NLTE grid was created with {self.atmosphere_grid}, but we are using {sme.atmo.source} for the synthesis"
                )

        # The possible labels
        self._teff = self.directory["teff"]
        self._grav = self.directory["grav"]
        self._feh = self.directory["feh"]
        self._xfe = self.directory["abund"]
        # The position of the models in the datafile
        self._keys = self.directory["models"].astype("U")

        try:
            depth_name = str.lower(sme.atmo.interp)
        except TypeError:
            logger.warning(
                "No interpolation axis specified in the atmosphere, using 'rhox'."
            )
            depth_name = "rhox"

        try:
            depth = self.directory[depth_name]
        except KeyError:
            other_depth_name = "tau" if depth_name == "rhox" else "rhox"
            depth = self.directory[other_depth_name]
            logger.warning(
                f"No data for {depth_name} in NLTE grid for {self.elem} found, "
                f"using {other_depth_name} instead."
            )
            depth_name = other_depth_name

        #:str: Which parameter is used as the depth axis
        self.depth_name = depth_name
        #:array(float): depth points of the atmosphere that is passed to the C library (in log10 scale)
        self._depth = depth

        self._grid = None
        self._points = None

        #:list(int): number of points in the grid to cache for each parameter, order; abund, teff, logg, monh
        self.subgrid_size = sme.nlte.subgrid_size
        #:float: Solar Abundance of the element
        self.abund_format = abund_format
        if solar is None:
            solar = Abund.solar()
        elif isinstance(solar, Abund):
            pass
        else:
            solar = Abund(0, solar)
        self.solar = solar

        #:dict: upper and lower parameters covered by the grid
        self.limits = {}
        #:array: NLTE data array
        self.bgrid = None
        #:array: Depth points of the NLTE data
        self.depth = None

        #:array: Indices of the lines in the NLTE data
        self.linerefs = None
        #:array: Indices of the lines in the LineList
        self.lineindices = None
        #:array: Indices of the lines in the bgrid
        self.iused = None

        #:str: citations in bibtex format, if known
        self.citation_info = ""

        self.first_warning = True

        conf = self.directory["conf"].astype("U")
        term = self.directory["term"].astype("U")
        try:
            species = self.directory["spec"].astype("U")
        except KeyError:
            logger.warning(
                f"Could not find 'species' field in NLTE file {self.grid_name}. Assuming they are all '{self.elem} 1'."
            )
            species = np.full(conf.shape, "%s 1" % self.elem)
            pass
        rotnum = self.directory["J"]  # rotational number of the atomic state
        if self.version[0] == 1 and self.version[1] >= 10:
            energies = self.directory["energy"]  # energy in eV
            self.citation_info = self.directory["citation"][()].decode()
        else:
            self.citation_info = None
            if self.selection != "levels":
                # logger.warning(
                #     "NLTE grid file version %s only supports level selection, not %s",
                #     self.version,
                #     self.selection,
                # )
                self.selection = "levels"

        if self.selection == "levels":
            self.lineindices, self.linerefs, self.iused = self.select_levels(
                conf, term, species, rotnum
            )
        elif self.selection == "energy":
            self.lineindices, self.linerefs, self.iused = self.select_energies(
                conf, term, species, rotnum, energies
            )

    def renew_linelist(
        self,
        sme,
    ):

        if 'use_indices' in sme.linelist._lines.columns:
            #:LineList: Whole LineList that was passed to the C library
            self.linelist = sme.linelist[sme.linelist['use_indices']]
            #:array(str): Elemental Species Names for the linelist
            self.species = sme.linelist.species[sme.linelist['use_indices']]
            # sme.linelist._lines = sme.linelist._lines.drop(columns=['use_indices'])
        else:
            #:LineList: Whole LineList that was passed to the C library
            self.linelist = sme.linelist
            #:array(str): Elemental Species Names for the linelist
            self.species = sme.linelist.species

        #:dict: upper and lower parameters covered by the grid
        self.limits = {}
        #:array: NLTE data array
        self.bgrid = None
        #:array: Depth points of the NLTE data
        self.depth = None

        #:array: Indices of the lines in the NLTE data
        self.linerefs = None
        #:array: Indices of the lines in the LineList
        self.lineindices = None
        #:array: Indices of the lines in the bgrid
        self.iused = None

        #:str: citations in bibtex format, if known
        self.citation_info = ""

        self.first_warning = True

        conf = self.directory["conf"].astype("U")
        term = self.directory["term"].astype("U")
        try:
            species = self.directory["spec"].astype("U")
        except KeyError:
            logger.warning(
                "Could not find 'species' field in NLTE file %s. Assuming they are all '%s 1'.",
                self.grid_name,
                self.elem,
            )
            species = np.full(conf.shape, "%s 1" % self.elem)
            pass
        rotnum = self.directory["J"]  # rotational number of the atomic state
        if self.version[0] == 1 and self.version[1] >= 10:
            energies = self.directory["energy"]  # energy in eV
            self.citation_info = self.directory["citation"][()].decode()
        else:
            self.citation_info = None
            if self.selection != "levels":
                # logger.warning(
                #     "NLTE grid file version %s only supports level selection, not %s",
                #     self.version,
                #     self.selection,
                # )
                self.selection = "levels"

        if self.selection == "levels":
            self.lineindices, self.linerefs, self.iused = self.select_levels(
                conf, term, species, rotnum
            )
        elif self.selection == "energy":
            self.lineindices, self.linerefs, self.iused = self.select_energies(
                conf, term, species, rotnum, energies
            )

    def solar_rel_abund(self, abund, elem):
        """Get the abundance of elem relative to H, i.e. [X/H]"""
        if self.abund_format == "H=12":
            # Its a bit more efficient, to simply use the values as is, instead of converting everything
            # just to end up back here
            idx = abund.elem_dict[elem]
            a = abund._pattern[idx]
            s = self.solar._pattern[idx]
            h = sh = 12
        else:
            a = abund.get_element(elem, type=self.abund_format)
            h = abund.get_element("H", type=self.abund_format)
            s = self.solar.get_element(elem, type=self.abund_format)
            sh = self.solar.get_element("H", type=self.abund_format)

        rabund = (a - h) - (s - sh)
        return rabund

    def scaled_rel_abund(self, abund):
        """Get the abundance of self.elem relative to Fe, i.e. [X/Fe]"""
        # Hydrogen NLTE grids use H itself as the abundance reference.
        # Using the generic [X/Fe] path would turn H into -[Fe/H] when the
        # selected abundance pattern differs from the default solar pattern.
        if self.elem == "H":
            return 0.0
        sel = self.solar_rel_abund(abund, self.elem)
        sfe = self.solar_rel_abund(abund, "Fe")
        rabund = sel - sfe
        return rabund

    def get(self, abund, teff, logg, monh, atmo):
        rabund = self.scaled_rel_abund(abund)

        if (len(self.limits) == 0 or not (
            (self.limits["xfe"][0] <= rabund <= self.limits["xfe"][-1])
            and (self.limits["teff"][0] <= teff <= self.limits["teff"][-1])
            and (self.limits["grav"][0] <= logg <= self.limits["grav"][-1])
            and (self.limits["feh"][0] <= monh <= self.limits["feh"][-1])
        )):
            _ = self.read_grid(rabund, teff, logg, monh)

        return self.interpolate(rabund, teff, logg, monh, atmo)

    # @profile
    def read_grid(self, rabund, teff, logg, monh):
        """Read the NLTE coefficients from the nlte_grid files for the given element
        The class will cache subgrid_size points around the target values as well

        Parameters
        ----------
        rabund : float
            relative (to solar) abundance of the element
        teff : float
            temperature in Kelvin
        logg : float
            surface gravity in log(cgs)
        monh : float
            Metallicity in H=12

        Returns
        -------
        nlte_grid : dict
            collection of nlte coefficient data (memmapped)
        linerefs : array (nlines,)
            linelevel descriptions (Energy level terms)
        lineindices: array (nlines,)
            indices of the used lines in the linelist
        """
        # find the n nearest parameter values in the grid (n == subgrid_size)
        x = np.argsort(np.abs(rabund - self._xfe))[: self.subgrid_size[0]]
        t = np.argsort(np.abs(teff - self._teff))[: self.subgrid_size[1]]
        g = np.argsort(np.abs(logg - self._grav))[: self.subgrid_size[2]]
        f = np.argsort(np.abs(monh - self._feh))[: self.subgrid_size[3]]

        x = x[np.argsort(self._xfe[x])]
        t = t[np.argsort(self._teff[t])]
        g = g[np.argsort(self._grav[g])]
        f = f[np.argsort(self._feh[f])]

        # Read the models with those parameters, and store depth and level
        # Create storage array
        ndepths, nlevel = self.directory[self._keys[0, 0, 0, 0]].shape
        nabund = len(x)
        nteff = len(t)
        ngrav = len(g)
        nfeh = len(f)

        self.bgrid = np.zeros((ndepths, nlevel, nabund, nteff, ngrav, nfeh))

        for i, j, k, l in tqdm(
            np.ndindex(nabund, nteff, ngrav, nfeh),
            desc="Loading NLTE %s" % self.elem,
            total=nabund * nteff * ngrav * nfeh,
            disable=not show_progress_bars,
        ):
            model = self._keys[f[l], g[k], t[j], x[i]]
            try:
                self.bgrid[:, :, i, j, k, l] = self.directory[model]
            except KeyError:
                warnings.warn(
                    f"Missing Model for element {self.elem}: T={self._teff[t[j]]}, logg={self._grav[g[k]]}, feh={self._feh[f[l]]}, abund={self._xfe[x[i]]:.2f}"
                )
        mask = np.zeros(self._depth.shape[:-1], bool)
        for i, j, k in itertools.product(f, g, t):
            mask[i, j, k] = True
        self.depth = self._depth[mask, :]
        self.depth.shape = nfeh, ngrav, nteff, -1

        # Reduce the stored data to only relevant energy levels
        # Remap the previous indices into a collapsed sequence
        # level_labels = level_labels[iused]
        if self.iused is not None:
            self.bgrid = self.bgrid[self.iused, ...]
        else:
            self.bgrid = self.bgrid[False]

        self._points = (
            self._xfe[x],
            self._teff[t],
            self._grav[g],
            self._feh[f],
        )
        self.limits = {
            "teff": self._points[1][[0, -1]],
            "grav": self._points[2][[0, -1]],
            "feh": self._points[3][[0, -1]],
            "xfe": self._points[0][[0, -1]],
        }

        return self.bgrid

    def select_levels(self, conf, term, species, rotnum):
        """
        Match our NLTE terms to transitions in the vald3-format linelist.

        Level descriptions in the vald3 long format look like this:
        'LS                 2p6.3s               2S'
        'LS                 2p6.3p               2P*'
        These are stored in line_term_low and line_term_upp.
        The array line_extra has dimensions [3 x nlines]. It stores J_low, E_up, J_up
        The sme.atomic array stores:
        0) atomic number, 1) ionization state, 2) wavelength (in A),
        3) excitation energy of lower level (in eV), 4) log(gf), 5) radiative,
        6) Stark, 7) and van der Waals damping parameters

        Parameters
        ----------
        conf : array of shape (nl,)
            electronic configuration (for identification), e.g., 2p6.5s
        term : array of shape (nl,)
            term designations (for identification), e.g., a5F
        species : array of shape (nl,)
            Element and ion for each atomic level (for identification), e.g. Na 1.
        rotnum : array of shape (nl,)
            rotational number J of atomic state (for identification).

        Returns
        -------
        lineindices : array of shape (nl,)
            Indices of the used lines in the linelist
        linerefs : array of shape (2, nlines,)
            Cross references for the lower and upper level in each transition,
            to their indices in the list of atomic levels.
            Missing levels use indices of -1.
        iused : array of shape (nl,)
            Indices used in the subgrid
        """

        lineindices = np.asarray(self.species, "U")
        lineindices = np.char.startswith(lineindices, self.elem)
        if not np.any(lineindices):
            logger.warning(f"No NLTE transitions for {self.elem} found")
            return None, None, np.full(len(species), False)

        sme_species = self.species[lineindices]

        # Extract data from linelist
        low = self.linelist["term_lower"][lineindices]
        upp = self.linelist["term_upper"][lineindices]
        # Remove quotation marks (if any are there)
        low = low.astype(str)
        upp = upp.astype(str)
        low = np.char.replace(low, "'", "")
        upp = np.char.replace(upp, "'", "")
        # Get only the relevant part
        parts_low = np.char.partition(low, " ")[:, (0, 2)]
        parts_upp = np.char.partition(upp, " ")[:, (0, 2)]

        # Transform into term symbol J (2*S+1) ?
        extra = self.linelist.extra[lineindices]
        extra = extra[:, [0, 2]] * 2 + 1
        extra = np.rint(extra).astype("i8")

        # Transform into term symbol J (2*S+1) ?
        rotnum = np.rint(2 * rotnum + 1).astype(int)

        # Create record arrays for each set of labels
        dtype = [
            ("species", sme_species.dtype),
            ("configuration", parts_upp.dtype),
            ("term", parts_upp.dtype),
            ("J", extra.dtype),
        ]
        level_labels = np.rec.fromarrays((species, conf, term, rotnum), dtype=dtype)
        line_label_low = np.rec.fromarrays(
            (sme_species, parts_low[:, 0], parts_low[:, 1], extra[:, 0]),
            dtype=dtype,
        )
        line_label_upp = np.rec.fromarrays(
            (sme_species, parts_upp[:, 0], parts_upp[:, 1], extra[:, 1]),
            dtype=dtype,
        )

        # Prepare arrays
        nlines = parts_low.shape[0]
        linerefs = np.full((nlines, 2), -1)
        iused = np.zeros(len(species), dtype=bool)

        # Loop through the NLTE levels
        # and match line levels
        for i, level in enumerate(level_labels):
            idx_l = line_label_low == level
            linerefs[idx_l, 0] = i
            iused[i] |= np.any(idx_l)

            idx_u = line_label_upp == level
            linerefs[idx_u, 1] = i
            iused[i] |= np.any(idx_u)

        lineindices = np.where(lineindices)[0]

        # Remap the linelevel references
        for j, i in enumerate(np.where(iused)[0]):
            linerefs[linerefs == i] = j

        return lineindices, linerefs, iused

    def select_energies(self, conf, term, species, rotnum, energies):
        lineindices = np.asarray(self.species, "U")
        lineindices = np.char.startswith(lineindices, self.elem)
        if not np.any(lineindices):
            logger.warning("No NLTE transitions for %s found", self.elem)
            return None, None, np.full(len(species), False)
        sme_species = self.species[lineindices]

        # Extract data from linelist
        term_low = self.linelist["term_lower"][lineindices].astype(str)
        term_upp = self.linelist["term_upper"][lineindices].astype(str)
        # Remove quotation marks (if any are there)
        term_low = np.char.replace(term_low, "'", "")
        term_upp = np.char.replace(term_upp, "'", "")
        # Split the string into configuration and term
        term_low = np.char.partition(term_low, " ")
        term_upp = np.char.partition(term_upp, " ")
        # Remove whitespaces
        term_low = np.char.replace(term_low, " ", "")
        term_upp = np.char.replace(term_upp, " ", "")

        # Get only the relevant part
        conf_low, term_low = term_low[:, 0], term_low[:, 2]
        conf_upp, term_upp = term_upp[:, 0], term_upp[:, 2]

        # Energy levels in the linelist in eV
        energy_low = self.linelist["excit"][lineindices]
        energy_upp = self.linelist["e_upp"][lineindices]

        # Transform into term symbol J (2*S+1) ?
        extra = self.linelist.extra[lineindices]
        extra = extra[:, [0, 2]] * 2 + 1
        extra = np.rint(extra).astype(int)
        j_low, j_upp = extra[:, 0], extra[:, 1]

        # Transform into term symbol J (2*S+1) ?
        rotnum = np.rint(2 * rotnum + 1).astype(int)

        dtype = [
            ("species", sme_species.dtype),
            ("configuration", term_low.dtype),
            ("term", term_low.dtype),
            ("J", extra.dtype),
            ("energy", energy_low.dtype),
        ]

        level_labels = np.rec.fromarrays(
            (species, conf, term, rotnum, energies), dtype=dtype
        )
        line_label_low = np.rec.fromarrays(
            (sme_species, conf_low, term_low, j_low, energy_low), dtype=dtype
        )
        line_label_upp = np.rec.fromarrays(
            (sme_species, conf_upp, term_upp, j_upp, energy_upp), dtype=dtype
        )

        nlines = term_low.shape[0]  # np.count_nonzero(lineindices)
        linerefs = np.full((nlines, 2), -1)
        iused = np.zeros(len(species), dtype=bool)

        idx_map = np.arange(len(species))
        # Maximum energy in the grid
        max_energy = np.max(energies)
        # Maximum seperation between energy levels
        if self.min_energy_diff is None:
            diff = np.abs(np.diff(energies))
            energy_diff_limit = np.min(diff[diff != 0])
        else:
            energy_diff_limit = np.abs(self.min_energy_diff)

        def match(label, level):
            # Try to match the label and the level
            # Using various metrics in decreasing order of confidence
            idx_species = label.species == level.species
            idx_conf = label.configuration == level.configuration
            idx_term = label.term == level.term
            idx_j = label.J == level.J
            # 1. Try to match conf/term/spec/J as usual
            idx = idx_species & idx_conf & idx_term & idx_j
            if np.any(idx):
                return np.where(idx)[0][0]
            # 2. If it fails, try to match conf/term/spec, but ignore J
            # TODO: Unless there is another match!
            idx = idx_species & idx_conf & idx_term  # & ~idx_j
            if np.any(idx):
                return np.where(idx)[0][0]
            # 3. Try to match H
            # If spec is 'H 1', then try to match conf with conf,
            # *or* try to match term with term, *or* try to match conf with term,
            # *or* try to match term with conf
            if level.species == "H 1":
                # TODO
                idx_term_conf = label.term == level.configuration
                idx_conf_term = label.configuration == level.term
                idx = idx_conf | idx_term | idx_term_conf | idx_conf_term
                if np.any(idx):
                    return np.where(idx)[0][0]
            # 4. Try to match energies (including H, if step 3 failed)
            # Find the level in the nlte grid with the same spec,
            # and the closest energy; *provided* that the desired energy does
            # not exceed the smallest energy difference between levels in the grid
            idx = idx_species
            if np.any(idx):
                if level.energy > max_energy or level.energy <= 0:
                    # We exceed the maximum energy of the grid, ignore NLTE
                    return -1
                diff = np.abs(label[idx].energy - level.energy)
                idx2 = np.argmin(diff)
                mindiff = diff[idx2]
                # difference needs to be smaller than some limit?
                if mindiff < energy_diff_limit:
                    return np.where(np.atleast_1d(idx_map[idx][idx2]))[0][0]
                else:
                    return -1
            # 5. If everything fails return nothing
            return -1

        # Loop through the linelist line levels
        # and match to NLTE line levels
        # match() return -1 of no match is found, or the index otherwise
        for i, level in enumerate(line_label_low):
            idx_l = match(level_labels, level)
            linerefs[i, 0] = idx_l
            iused[idx_l] |= idx_l != -1

        for i, level in enumerate(line_label_upp):
            idx_u = match(level_labels, level)
            linerefs[i, 1] = idx_u
            iused[idx_u] |= idx_u != -1

        # Lineindices as integer pointers
        lineindices = np.where(lineindices)[0]
        # Remap the linelevel references
        for j, i in enumerate(np.where(iused)[0]):
            linerefs[linerefs == i] = j

        return lineindices, linerefs, iused

    def interpolate(self, rabund, teff, logg, monh, atmo):
        """
        interpolate nlte coefficients on the model grid

        Parameters
        ----------
        rabund : float
            relative (to solar) abundance of the element
        teff : float
            temperature in Kelvin
        logg : float
            surface gravity in log(cgs)
        monh : float
            Metallicity in H=12

        Returns
        -------
        subgrid : array (ndepth, nlines)
            interpolated grid values
        """

        assert self._points is not None
        if self.bgrid is None or self.bgrid.shape[0] == 0:
            return None

        # Interpolate the depth scale to the target depth, this is unstructured data
        # i.e. each combination of parameters has a different depth scale (given in depth)
        ndepths, _, *nparam = self.bgrid.shape
        target_depth = atmo[self.depth_name]
        target_depth = np.log10(target_depth)
        ntarget = len(target_depth)

        param = [rabund, teff, logg, monh]
        right = [
            (self._points[i][-1] == p) and len(self._points[i] > 1)
            for i, p in enumerate(param)
        ]
        iabund = np.digitize(rabund, self._points[0], right=right[0]) - 1
        iteff = np.digitize(teff, self._points[1], right=right[1]) - 1
        ilogg = np.digitize(logg, self._points[2], right=right[2]) - 1
        imonh = np.digitize(monh, self._points[3], right=right[3]) - 1

        # Interpolate on the grid
        # We only interpolate on the 2**4 points around the requested values
        # self._points and self._grid are interpolated when reading the data in read_grid
        target = (rabund, teff, logg, monh)
        points = (
            self._points[0][iabund : iabund + 2],
            self._points[1][iteff : iteff + 2],
            self._points[2][ilogg : ilogg + 2],
            self._points[3][imonh : imonh + 2],
        )
        npoints = (
            max(points[0].size, 1),
            max(points[1].size, 1),
            max(points[2].size, 1),
            max(points[3].size, 1),
        )

        grid = np.empty((*npoints, ntarget, ndepths), float)
        for l, x, t, g, f in np.ndindex(ndepths, *npoints):
            xp = self.depth[imonh + f, ilogg + g, iteff + t, :]
            yp = self.bgrid[l, :, iabund + x, iteff + t, ilogg + g, imonh + f]
            xp = np.log10(xp)
            grid[x, t, g, f, :, l] = interpolate.interp1d(
                xp,
                yp,
                bounds_error=False,
                fill_value="extrapolate",
                kind="cubic",
            )(target_depth)

        # Check if we need to extrapolate
        if any(
            [t < min(p) or t > max(p) for t, p in zip(target, points) if len(p) > 0]
        ):
            if self.first_warning:
                logger.warning(
                    f"Extrapolate on the {self.elem} NLTE grid. Requested values of {target} on grid {points}"
                )
                self.first_warning = False

        # Some grids have only one value in that direction
        # Usually in abundance. Then we need to remove that dimension
        # to avoid nan output
        mask = [len(p) > 1 for p in points]
        if not all(mask):
            points = [p for m, p in zip(mask, points) if m]
            target = [t for m, t in zip(mask, target) if m]
            idx = [slice(None, None) if m else 0 for m in mask]
            grid = grid[tuple(idx)]

        method = "order"
        if method == "grid":
            # TODO: Interpolate with splines
            # Possibly in order of importance, since scipy doesn't have spline interpolation on a grid
            subgrid = interpolate.interpn(
                points,
                grid,
                target,
                method="linear",
                bounds_error=False,
                fill_value=None,
            )
            subgrid = subgrid[0]
        elif method == "order":
            for p, t in zip(points, target):
                grid = interpolate.interp1d(
                    p,
                    grid,
                    axis=0,
                    bounds_error=False,
                    fill_value="extrapolate",
                )(t)
            subgrid = grid
        return subgrid


@CollectionFactory
class NLTE(Collection):
    # fmt: off
    _fields = Collection._fields + [
        ("elements", [], astype(list), this,
            "list: elements for which nlte calculations will be performed"),
        ("grids", {}, astype(dict), this,
            "dict: nlte grid datafiles for each element"),
        ("subgrid_size", [3, 3, 3, 3], array(4, int), this,
            "array of shape (4,): defines size of nlte grid cache."
            "Each entry is for one parameter abund, teff, logg, monh"),
        ("flags", None, array(None, np.bool_), this,
            "array: contains a flag for each line, whether it was calculated in NLTE (True) or not (False)"),
        ("solar", None, this, this, "str: defines which default to use as the solar metallicitiies"),
        ("abund_format", "H=12", astype(str), this, "str: which abundance format to use for comparison"),
        ("selection", "energy", oneof("energy", "levels"), this, "str: which selection algorithm to use to match linelist and departure coefficients"),
        ("min_energy_diff", None, this, this, "float: difference between energy levels that are still matched. If None will default to the smallest non zero difference between energy levels in the grid.")
    ]
    # fmt: on

    # For marcs2012 atmosphere
    _default_grids = {
        "H": "nlte_H_pysme.grd",
        "Li": "nlte_Li_pysme.grd",
        "C": "nlte_C_pysme.grd",
        "N": "nlte_N_pysme.grd",
        "O": "nlte_O_pysme.grd",
        "Na": "nlte_Na_pysme.grd",
        "Mg": "nlte_Mg_pysme.grd",
        "Al": "nlte_Al_pysme.grd",
        "Si": "nlte_Si_pysme.grd",
        "K": "nlte_K_pysme.grd",
        "Ca": "nlte_Ca_pysme.grd",
        "Ti": "nlte_Ti_pysme.grd",
        "Mn": "nlte_Mn_pysme.grd",
        "Fe": "nlte_Fe_pysme.grd",
        "Ba": "nlte_Ba_pysme.grd",        
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.first = False

        if "solar" in kwargs.keys():
            self.solar = kwargs["solar"]

        # Convert IDL keywords to Python
        if "nlte_elem_flags" in kwargs.keys():
            elements = kwargs["nlte_elem_flags"]
            self.elements = [abund_elem[i] for i, j in enumerate(elements) if j == 1]

        if "nlte_subgrid_size" in kwargs.keys():
            self.subgrid_size = kwargs["nlte_subgrid_size"]

        if "nlte_grids" in kwargs:
            grids = kwargs["nlte_grids"]
            if isinstance(grids, (list, np.ndarray)):
                grids = {
                    abund_elem[i]: name.decode()
                    for i, name in enumerate(grids)
                    if name != ""
                }
            self.grids = grids

        # TODO
        #:dict: the cached subgrid data for each element
        # This is NOT saved on sme.save
        # But maybe should be ?
        self.grid_data = {}

    def set_nlte(self, element, grid=None):
        """
        Add an element to the NLTE calculations

        Parameters
        ----------
        element : str
            The abbreviation of the element to add to the NLTE calculations
        grid : str, optional
            Filename of the NLTE data grid to use for this element
            the file must be in nlte_grids directory
            Defaults to a set of "known" files for some elements
        """
        if element in self.elements:
            # Element already in NLTE
            # Change grid if given
            if grid is not None:
                self.grids[element] = grid
            return

        if grid is None:
            # Use default grid
            if element not in NLTE._default_grids.keys():
                raise ValueError(f"No default grid known for element {element}")
            grid = NLTE._default_grids[element]
            logger.info("Using default grid %s for element %s", grid, element)

        if element not in self.elements:
            self.elements += [element]
        self.grids[element] = grid

    def remove_nlte(self, element):
        """
        Remove an element from the NLTE calculations

        Parameters
        ----------
        element : str
            Abbreviation of the element to remove from NLTE
        """
        if element not in self.elements:
            # Element not included in NLTE anyways
            return

        self.elements.remove(element)
        self.grids.pop(element)

    @property
    def _citation_info(self):
        citations = [
            grid.citation_info
            for el, grid in self.grid_data.items()
            if grid.citation_info is not None
        ]
        citations = "\n".join(citations)
        return citations

    @_citation_info.setter
    def _citation_info(self, value):
        pass

    # @profile
    def update_coefficients(self, sme, dll, lfs_nlte):
        """pass departure coefficients to C library;
            If bmat, linerefs, lineindices are given, use those instead of recalculating them 
            and the code will not check if the line level, energy is matched or not. 
        """

        # Reset the departure coefficient every time, just to be sure
        # It would be more efficient to just Update the values, but this doesn't take long
        dll.ResetDepartureCoefficients()

        # Only print "Running in NLTE" message on the first run each time
        if np.all(self.grids == "") or np.size(self.elements) == 0:
            # No NLTE to do
            if self.first:
                self.first = False
                logger.info("Running in LTE")
            return sme
        if sme.linelist.lineformat == "short":
            if self.first:
                self.first = False
                logger.warning(
                    "NLTE line formation was requested, but VALD3 long-format linedata\n"
                    "are required in order to relate line terms to NLTE level corrections!\n"
                    "Line formation will proceed under LTE."
                )
            return sme

        if self.first:
            self.first = False
            logger.info("Running in NLTE: %s", ", ".join(self.elements))

        # Call each element to update and return its set of departure coefficients
        marked_for_removal = []
        for elem in self.elements:
            # Call function to retrieve interpolated NLTE departure coefficients
            # the abundances for NLTE are handled in the H-12 format
            grid = self.get_grid(sme, elem, lfs_nlte)

            if not np.any(grid.iused):
                # No lines are found for this element
                # remove it from the elements after this loop
                logger.warning(
                    "No %s NLTE lines found, removing it from NLTE calculations",
                    elem,
                )
                marked_for_removal += [elem]
                continue
            bmat = grid.get(sme.abund, sme.teff, sme.logg, sme.monh, sme.atmo)
            if bmat is None or grid.linerefs.size == 0:
                logger.warning(
                    "No %s NLTE lines found, removing it from NLTE calculations",
                    elem,
                )
                marked_for_removal += [elem]
                continue
            else:
                # Put corrections into the nlte_b matrix, don't cache the data
                for lr, li in zip(grid.linerefs, grid.lineindices):
                    # loop through the list of relevant _lines_, substitute both their levels into the main b matrix
                    # Make sure both levels have corrections available
                    if lr[0] != -1 and lr[1] != -1:
                        dll.InputDepartureCoefficients(bmat[:, lr], li)
        for elem in marked_for_removal:
            self.remove_nlte(elem)

        return sme

    # @profile
    def get_grid(self, sme, elem, lfs_nlte):
        """Read and interpolate the NLTE grid for the current element and parameters"""
        if self.grids[elem] is None:
            raise ValueError(f"Element {elem} has not been prepared for NLTE")

        # The grids are cached in the NLTE object, but not saved
        if elem in self.grid_data.keys():
            grid = self.grid_data[elem]
            grid.renew_linelist(sme)
        else:
            grid = Grid(
                sme,
                elem,
                lfs_nlte,
                solar=self.solar,
                abund_format=self.abund_format,
                selection=self.selection,
                min_energy_diff=self.min_energy_diff,
            )
            self.grid_data[elem] = grid

        return grid


# These are kept here for compatibility, but it is recommended to use the functions of sme.nlte
def nlte(sme, dll, elem, lfs_nlte):
    """Read and interpolate the NLTE grid for the current element and parameters"""
    grid = sme.nlte.get_grid(sme, elem, lfs_nlte)
    subgrid = grid.get(sme.abund, sme.teff, sme.logg, sme.monh, sme.atmo)
    return subgrid, grid.linerefs, grid.lineindices


def update_nlte_coefficients(sme, dll, lfs_nlte):
    """pass departure coefficients to C library"""
    return sme.nlte.update_coefficients(sme, dll, lfs_nlte)
