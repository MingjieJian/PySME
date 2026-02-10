# -*- coding: utf-8 -*-
""" Wrapper for SME C library """

import logging
import os, sys
from ctypes import cdll
from os.path import abspath, dirname, join, normpath
from posixpath import expanduser
from .smelib import libtools
import numpy as np

from .smelib.libtools import get_full_datadir, get_full_libfile, _parse_needed_arch_from_error, download_smelib

logger = logging.getLogger(__name__)

# Load the library
# This ensures that we can find the library when we load _smelib
libfile = get_full_libfile()

try:
    cdll.LoadLibrary(libfile)
except OSError as e:
    msg = str(e)
    # 仅在 macOS 且架构不匹配时尝试自愈
    if sys.platform == "darwin" and ("incompatible architecture" in msg or "mach-o file" in msg):
        need = _parse_needed_arch_from_error(msg)
        print("Detected arch mismatch; error says need:", need or "(unknown)")
        # 如果错误里没写 need，就按当前解释器架构
        if need is not None:
            download_smelib(force_arch=need)
            libtools.compile_interface()
            cdll.LoadLibrary(libfile)

from .smelib import _smelib

CURRENT_LIB = libfile


def reload_lib(libfile):
    global _smelib
    global CURRENT_LIB
    if libfile != CURRENT_LIB:
        del _smelib
        cdll.LoadLibrary(libfile)
        from .smelib import _smelib
    CURRENT_LIB = libfile


class SME_DLL:
    """Object Oriented interface for the SME C library"""

    def __init__(self, libfile=None, datadir=None):
        self.libfile = libfile
        reload_lib(libfile)

        self.SetLibraryPath(datadir)
        self.check_data_files_exist()

    @property
    def datadir(self):
        """str: Expected directory of the data files"""
        return self.GetLibraryPath()

    @datadir.setter
    def datadir(self, value):
        self.SetLibraryPath(value)

    @property
    def file(self):
        """str: (Expected) Location of the library file"""
        return libfile

    def check_data_files_exist(self):
        """
        Checks if required data files for the SME library exist.
        If they dont exist, SME will just segfault, without any hint.

        Raises
        ------
        FileNotFoundError
            If any of the files don't exist
        """
        names = self.GetDataFiles()
        directory = self.GetLibraryPath()
        for name in names:
            n = os.path.join(directory, name)
            if not os.path.exists(n):
                raise FileNotFoundError(
                    "Could not find required data file {name} in library directory {directory}".format(
                        name=name, directory=directory
                    )
                )

    def SMELibraryVersion(self):
        """Return SME library version"""
        return _smelib.LibraryVersion()

    def GetLibraryPath(self):
        """Get the data file directory"""
        return _smelib.GetLibraryPath()

    def GetDataFiles(self):
        """Get the required data files"""
        files = _smelib.GetDataFiles()
        return files.split(";")

    def SetLibraryPath(self, libpath=None):
        """Set the path to the library"""
        if libpath is None:
            libpath = get_full_datadir()
        libpath = normpath(libpath) + os.sep
        _smelib.SetLibraryPath(libpath)

    def InputWaveRange(self, wfirst, wlast):
        """
        Read in Wavelength range

        Will raise an exception if wfirst is larger than wlast

        Parameters
        ----------
        wfirst : float
            first wavelength of the segment
        wlast : float
            last wavelength of the segment
        """
        _smelib.InputWaveRange(wfirst, wlast)

    def SetVWscale(self, gamma6):
        """
        Set van der Waals scaling factor

        Parameters
        ----------
        gamma6 : float
            van der Waals scaling factor
        """
        _smelib.SetVWscale(gamma6)

    def SetH2broad(self, h2_flag=True):
        """Set flag for H2 molecule"""
        if h2_flag:
            _smelib.SetH2broad()
        else:
            _smelib.ClearH2broad()

    def ClearH2broad(self):
        """Clear flag for H2 molecule"""
        self.SetH2broad(False)

    def InputLineList(self, linelist):
        """
        Read in line list

        Parameters
        ---------
        atomic : array of size (nlines, 8)
            atomic linelist data for each line
            fields are: atom_number, ionization, wlcent, excit, gflog, gamrad, gamqst, gamvw
        species : array(string) of size (nlines,)
            names of the elements (with Ionization level)
        """

        # Check for unsupported lines in the linelist
        # The limits are given by which data is available in XSAHA in eos.f
        # The one letter elements need the space, as this is what we get from the
        # linelist species
        max_ion = {
            "H ": 2,
            "He": 3,
            "Li": 4,
            "Be": 4,
            "B ": 4,
            "C ": 6,
            "N ": 6,
            "O ": 6,
            "F ": 6,
            "Ne": 6,
            "Na": 6,
            "Mg": 6,
            "Al": 6,
            "Si": 6,
            "P ": 6,
            "S ": 6,
            "Cl": 5,
            "Ar": 5,
            "K ": 5,
            "Ca": 5,
            "Sc": 5,
            "Ti": 5,
            "V ": 5,
            "Cr": 5,
            "Mn": 5,
            "Fe": 5,
            "Co": 5,
            "Ni": 5,
            "Cu": 3,
            "Zn": 3,
            "Ga": 3,
            "Ge": 3,
            "As": 3,
            "Se": 3,
            "Br": 3,
            "Kr": 3,
            "Rb": 3,
            "Sr": 3,
            "Y ": 3,
            "Zr": 3,
            "Nb": 3,
            "Mo": 3,
            "Tc": 3,
            "Ru": 3,
            "Rh": 3,
            "Pd": 3,
            "Ag": 3,
            "Cd": 3,
            "In": 3,
            "Sn": 3,
            "Sb": 3,
            "Te": 3,
            "I ": 3,
            "Xe": 3,
            "Cs": 3,
            "Ba": 3,
            "La": 3,
            "Ce": 4,
            "Pr": 4,
            "Nd": 4,
            "Pm": 3,
            "Sm": 3,
            "Eu": 4,
            "Gd": 3,
            "Tb": 3,
            "Dy": 3,
            "Ho": 3,
            "Er": 3,
            "Tm": 3,
            "Yb": 3,
            "Lu": 3,
            "Hf": 3,
            "Ta": 3,
            "W ": 3,
            "Re": 3,
            "Os": 3,
            "Ir": 3,
            "Pt": 3,
            "Au": 3,
            "Hg": 3,
            "Tl": 3,
            "Pb": 3,
            "Bi": 3,
            "Po": 3,
            "At": 3,
            "Rn": 3,
            "Fr": 3,
            "Ra": 3,
            "Ac": 3,
            "Th": 3,
            "Pa": 3,
            "U ": 3,
            "Np": 3,
            "Pu": 3,
            "Am": 3,
            "Cm": 3,
            "Bk": 3,
            "Cf": 3,
            "Es": 3,
        }

        elem = linelist["species"].astype("U2")
        ion = linelist["ionization"]

        # For each element check if we are above the maximum ionization
        mask = np.full(len(linelist), False)
        for key, value in max_ion.items():
            mask |= (elem == key) & (ion > value)
        nlines = np.count_nonzero(mask)
        # And remove those lines
        if nlines > 0:
            warning = ", ".join(np.unique(linelist[mask]["species"]))
            logger.warning(
                "Discarding %i high ionization lines, of species %s",
                nlines,
                warning,
            )
            linelist = linelist[~mask]

        # Prepare arrays for the C code
        atomic = linelist["atomic"].T
        species = linelist["species"]
        species = np.asarray(species, "S8")

        _smelib.InputLineList(species, atomic)

        return mask

    def OutputLineList(self):
        """
        Return line list

        Returns
        -------
        atomic : array of size (nlines, 6)
            relevant data of the linelist
            wlcent, excit, gflog, gamrad, gamqst, gamvw
        """
        return _smelib.OutputLineList()

    def UpdateLineList(self, atomic, species, index):
        """
        Change line list parameters

        Parameters
        ---------
        atomic : array of size (nlines, 8)
            atomic linelist data for each line
            fields are: atom_number, ionization, wlcent, excit, gflog, gamrad, gamqst, gamvw
        species : array(string) of size (nlines,)
            names of the elements (with Ionization level)
        index : array(int) of size (nlines,)
            indices of the lines to update relative to the overall linelist
        """
        index = np.asarray(index, dtype=np.int16)
        species = species[index].astype("S8")
        atomic = atomic[index].T

        _smelib.UpdateLineList(
            species,
            atomic,
            index,
        )

    def InputModel(self, teff, grav, vturb, atmo):
        """Read in model atmosphere

        Parameters
        ---------
        teff : float
            effective Temperature in Kelvin
        grav : float
            surface gravity in log10(cgs)
        vturb : float
            turbulence velocity in km/s
        atmo : Atmo
            atmosphere structure (see Atmo for details)
        """

        if teff <= 0:
            raise ValueError("Temperature must be positive (unit is Kelvin)")
        if vturb < 0:
            raise ValueError("Turbulence velocity must be positive or zero")

        motype = atmo["depth"]
        depth = atmo[motype]
        ndepth = len(depth)
        temp = atmo["temp"]
        xne = atmo["xne"]
        xna = atmo["xna"]
        rho = atmo["rho"]
        vt = np.full(ndepth, vturb) if np.size(vturb) == 1 else vturb
        wlstd = atmo["wlstd"]
        opflag = atmo["opflag"].astype(np.int16)

        kwargs = {}
        if atmo["geom"] == "SPH":
            kwargs["radius"] = atmo["radius"]
            kwargs["height"] = atmo["height"]
            motype = "SPH"

        _smelib.InputModel(
            teff, grav, wlstd, motype, opflag, depth, temp, xne, xna, rho, vt, **kwargs
        )

    def InputAbund(self, abund):
        """
        Pass abundances to radiative transfer code.

        Calculate elemental abundances from abundance pattern and metallicity.
        Metallicity adjustment is not applied to H or He.
        Renormalize abundances after applying metallicity.
        Introduced limiter in case the proposed step in abundance is too large.

        Parameters
        ---------
        abund : Abund
            abundance structure to be passed (see Abund for more details)
        """
        abund = abund("sme", raw=True)
        # TODO: does this make sense or should we fix the C code?
        abund[np.isnan(abund)] = -99
        abund[abund < -99] = -99

        _smelib.InputAbund(abund)

    def Opacity(self):
        """Calculate opacities"""
        _smelib.Opacity()

    def GetOpacity(self, switch, species=None, key=None):
        """
        Returns specific cont. opacity, different output depending on the input

        Parameters
        ----------
        switch : str
            one of [COPSTD, COPRED, COPBLU, AHYD, AH2P, AHMIN, SIGH, AHE1, AHE2,
            AHEMIN, SIGHE, ACOOL, ALUKE, AHOT, SIGEL, SIGH2]
        key : str, optional
            for ACOOL, one of [new, old, fraction]
        species : str, optional
            for ACOOL and ALUKE it specifies the element
            ACOOL: C1, Mg1, Al1, Si1, Fe1, CH, NH, OH
            ALUKE: N1, O1, Mg2, Si2, Ca2
        """
        kwargs = {}
        if key is not None:
            kwargs["key"] = key
        if species is not None:
            kwargs["species"] = species
        return _smelib.GetOpacity(switch, **kwargs)

    def Ionization(self, ion=0):
        """
        Calculate ionization balance for current atmosphere and abundances.
        Ionization state is stored in the external library.
        Ion is a bit flag with values (add them together to use multiple):

        :1: adopt particle number densities from EOS
        :2: adopt electron number densities from EOS
        :4: adopt gas densities (g/cm^3) from EOS

        instead of using values from model atmosphere. Different abundance patterns
        in the model atmosphere (usually scaled solar) and SME (may be non-solar)
        can affect line shape, e.g. shape of hydrogen lines.

        Parameters
        ----------
        ion : int
            flag that determines the behaviour of the C function
        """
        _smelib.Ionization(ion)

    def GetDensity(self):
        """
        Retrieve density in each layer

        Returns
        -------
        density : array of size (ndepth,)
            Density of the atmosphere in each layer
        """
        return _smelib.GetDensity()

    def GetNatom(self):
        """
        Get XNA

        Returns
        -------
        XNA : array of size (ndepth,)
            XNA in each layer
        """
        return _smelib.GetNatom()

    def GetNelec(self):
        """
        Get XNE (Electron number density) for each layer in the atmosphere

        Returns
        -------
        XNE : array of size (ndepth,)
            XNE in each layer
        """
        return _smelib.GetNelec()

    def GetFraction(self, species, mode=0):
        """
        Get species fraction/partition function/number density.

        Parameters
        ----------
        species : str
            Species name in SPLIST (e.g., 'Fe', 'Fe+', 'CO')
        mode : int, optional
            0: number density (FRACT * PF)
            1: partition function
            other: FRACT
        """
        return _smelib.GetFraction(species, mode=mode)

    def Transf(
        self,
        mu,
        wave=None,
        nwmax=400000,
        accrt=1e-3,
        accwi=3e-3,
        keep_lineop=False,
        long_continuum=True,
    ):
        """
        Radiative Transfer Calculation

        Perform the radiative transfer calculation thorugh the atmosphere
        Requires that all parameters have been set beforehand

        Parameters
        ---------
        mu : array of shape (nmu,)
            mu angles (1 - cos(phi)) of different limb points along the stellar surface
        accrt : float
            accuracy of the radiative transfer integration
        accwi : float
            accuracy of the interpolation on the wavelength grid
        keep_lineop : bool, optional
            if True do not recompute the line opacities (default: False)
        long_continuum : bool, optional
            if True the continuum is calculated at every wavelength (default: True)
        nwmax : int, optional
            maximum number of wavelength points if wavelength grid is not set with wave (default: 400000)
        wave : array, optional
            wavelength grid to use for the calculation,
            if not set will use an adaptive wavelength grid with no constant step size (default: None)

        Returns
        -------
        nw : int
            number of actual wavelength points, i.e. size of wint_seg
        wint_seg : array of shape (nwave,)
            wavelength grid, the number of wavelengthpoints is equal to the number of lines * 2 - 1
            One point in the center of each line + plus one between the next line
        sint_seg : array of shape (nmu, nwave)
            spectrum for each mu point
        cint_seg : array of shape (nmu, nwave)
            continuum for each mu point
        """
        keep_lineop = 1 if keep_lineop else 0
        long_continuum = 1 if long_continuum else 0

        # keywords = {"mu", "wave", "nwmax", "accrt", "accwi", "keep_lineop", "long_continuum"}
        nw, wave, sint, cint = _smelib.Transf(
            mu, wave, nwmax, accrt, accwi, keep_lineop, long_continuum
        )

        # Resize the arrays
        wave = wave[:nw]
        sint = sint[:nw, :].T
        cint = cint[:nw, :].T

        sint = np.nan_to_num(sint, copy=False)
        cint = np.nan_to_num(cint, copy=False)

        return nw, wave, sint, cint

    def CentralDepth(self, mu, accrt):
        """
        This subroutine explicitly solves the transfer equation
        for a set of nodes on the star disk in the centers of spectral
        lines. The results are specific intensities.

        Parameters
        ----------
        mu : array of size (nmu,)
            mu values along the stellar disk to calculate
        accrt : float
            precision of the radiative transfer calculation

        Returns
        -------
        table : array of size (nlines,)
            Centeral depth (i.e. specific intensity) of each line
        """

        return _smelib.CentralDepth(mu, accrt)

    def GetLineOpacity(self, wave):
        """
        Retrieve line opacity data from the C library

        Parameters
        ----------
        wave : float
            Wavelength of the line opacity to retrieve

        Returns
        ---------
        lop : array
            line opacity
        cop : array
            continuum opacity including scatter
        scr : array
            Scatter
        tsf : array
            Total source function
        csf : array
            Continuum source function
        """
        lop, cop, scr, tsf, csf = _smelib.GetLineOpacity(wave)
        return lop, cop, scr, tsf, csf

    def GetLineRange(self):
        """Get the effective wavelength range for each line
        i.e. the wavelengths for which the line has significant impact

        Returns
        -------
        linerange : array of size (nlines, 2)
            lower and upper wavelength for each spectral line
        """
        return _smelib.GetLineRange()

    def InputDepartureCoefficients(self, bmat, lineindex):
        """
        Input NLTE departure coefficients

        Parameters
        ----------
        bmat : array of size (2, ndepth)
            departure coefficient matrix
        lineindex : float
            index of the line in the linelist
        """
        bmat = np.atleast_2d(bmat)
        _smelib.InputDepartureCoefficients(bmat, lineindex)

    def GetDepartureCoefficients(self, line):
        """Get the NLTE departure coefficients as stored in the C library

        Parameters
        ----------
        line : int
            requested line index, i.e. between 0 and number of lines

        Returns
        -------
        bmat : array of size (2, nrhox)
            departure coefficients for the given line index
        """
        return _smelib.GetDepartureCoefficients(line)

    def ResetDepartureCoefficients(self):
        """Reset departure coefficients from any previous call, to ensure LTE as default"""
        _smelib.ResetDepartureCoefficients()

    def GetNLTEflags(self):
        """Get an array that tells us which lines have been used with NLTE correction

        Returns
        -------
        nlte_flags : array(bool) of size (nlines,)
            True if line was used with NLTE, False if line is only LTE
        """
        nlte_flags = _smelib.GetNLTEflags()
        return nlte_flags.astype(bool)

    def GetContributionfunction(self, mu, wave, nwmax=400000, long_continuum=True):
        """Get the contribution function for a wavelength grid.

        Returns
        -------

        """
        long_continuum_val = 1 if long_continuum else 0

        res = _smelib.ContributionFunctions(mu, wave, nwmax, long_continuum_val)
        table, ctable = res
        if long_continuum:
            return table, ctable
        else:
            return table
