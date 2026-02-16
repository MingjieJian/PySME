# -*- coding: utf-8 -*-
"""
Elemental abundance data handling module
"""
import io
import json
import logging

import numpy as np
from flex.extensions.bindata import BinaryDataExtension

from .persistence import IPersist

logger = logging.getLogger(__name__)

_citation_anders1989 = r"""
@ARTICLE{1989GeCoA..53..197A,
       author = {{Anders}, E. and {Grevesse}, N.},
        title = "{Abundances of the elements: Meteoritic and solar}",
      journal = {\gca},
     keywords = {Abundance, Chemical Analysis, Chondrites, Meteoritic Composition, Photosphere, Solar Corona, Energetic Particles, Isotopes, Solar Wind, METEORITES, ABUNDANCE, SAMPLES, METEORITE, CI CHONDRITES, SUN, CORONA, PHOTOSPHERE, REVIEW, ELEMENTS, FRACTIONATION, SOLAR SYSTEM, NUCLIDES, SOLAR WIND, CHARGED PARTICLES, RARE GASES, REFRACTORY ELEMENTS, COMETS, SIDEROPHILES, VOLATILES, NUCLEOSYNTHESIS, HALLEY, COMPARISONS},
         year = 1989,
        month = jan,
       volume = {53},
       number = {1},
        pages = {197-214},
          doi = {10.1016/0016-7037(89)90286-X},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1989GeCoA..53..197A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_grevesse1996 = r"""
@INPROCEEDINGS{1996ASPC...99..117G,
       author = {{Grevesse}, N. and {Noels}, A. and {Sauval}, A.~J.},
        title = "{Standard Abundances}",
    booktitle = {Cosmic Abundances},
         year = 1996,
       editor = {{Holt}, Stephen S. and {Sonneborn}, George},
       series = {Astronomical Society of the Pacific Conference Series},
       volume = {99},
        month = jan,
        pages = {117},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1996ASPC...99..117G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_grevesse1998 = r"""
@ARTICLE{1998SSRv...85..161G,
       author = {{Grevesse}, N. and {Sauval}, A.~J.},
        title = "{Standard Solar Composition}",
      journal = {\ssr},
     keywords = {Sun: abundances, Meteorites: abundances, Solar spectroscopy},
         year = 1998,
        month = may,
       volume = {85},
        pages = {161-174},
          doi = {10.1023/A:1005161325181},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_asplund2005 = r"""
@INPROCEEDINGS{2005ASPC..336...25A,
       author = {{Asplund}, M. and {Grevesse}, N. and {Sauval}, A.~J.},
        title = "{The Solar Chemical Composition}",
    booktitle = {Cosmic Abundances as Records of Stellar Evolution and Nucleosynthesis},
         year = 2005,
       editor = {{Barnes}, III, Thomas G. and {Bash}, Frank N.},
       series = {Astronomical Society of the Pacific Conference Series},
       volume = {336},
        month = sep,
        pages = {25},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2005ASPC..336...25A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_grevesse2007 = r"""
@ARTICLE{2007SSRv..130..105G,
       author = {{Grevesse}, N. and {Asplund}, M. and {Sauval}, A.~J.},
        title = "{The Solar Chemical Composition}",
      journal = {\ssr},
     keywords = {Sun: abundances, photosphere, corona},
         year = 2007,
        month = jun,
       volume = {130},
       number = {1-4},
        pages = {105-114},
          doi = {10.1007/s11214-007-9173-7},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2007SSRv..130..105G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_asplund2009 = r"""
@ARTICLE{2009ARA&A..47..481A,
    author = {{Asplund}, Martin and {Grevesse}, Nicolas and {Sauval}, A. Jacques and
    {Scott}, Pat},
    title = "{The Chemical Composition of the Sun}",
    journal = {\araa},
    keywords = {Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
    year = "2009",
    month = "Sep",
    volume = {47},
    number = {1},
    pages = {481-522},
    doi = {10.1146/annurev.astro.46.060407.145222},
    archivePrefix = {arXiv},
    eprint = {0909.0948},
    primaryClass = {astro-ph.SR},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2009ARA&A..47..481A},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_asplund2021 = r"""
@ARTICLE{2021A&A...653A.141A,
       author = {{Asplund}, M. and {Amarsi}, A.~M. and {Grevesse}, N.},
        title = "{The chemical make-up of the Sun: A 2020 vision}",
      journal = {\aap},
     keywords = {Sun: abundances, Sun: photosphere, Sun: atmosphere, line: formation, meteorites, meteors, meteoroids, Sun: helioseismology, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2021,
        month = sep,
       volume = {653},
          eid = {A141},
        pages = {A141},
          doi = {10.1051/0004-6361/202140445},
archivePrefix = {arXiv},
       eprint = {2105.01661},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021A&A...653A.141A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_lodders2003 = r"""
@ARTICLE{2003ApJ...591.1220L,
    author = {{Lodders}, Katharina},
    title = "{Solar System Abundances and Condensation Temperatures of the Elements}",
    journal = {\apj},
    keywords = {Astrochemistry, Meteors, Meteoroids, Solar System: Formation- Sun: Abundances, Sun: Photosphere},
    year = "2003",
    month = "Jul",
    volume = {591},
    number = {2},
    pages = {1220-1247},
    doi = {10.1086/375492},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2003ApJ...591.1220L},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_lodders2010 = r"""
@INPROCEEDINGS{2010ASSP...16..379L,
       author = {{Lodders}, Katharina},
        title = "{Solar System Abundances of the Elements}",
     keywords = {Physics, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
    booktitle = {Principles and Perspectives in Cosmochemistry},
         year = 2010,
       editor = {{Goswami}, Aruna and {Reddy}, B. Eswar},
       series = {Astrophysics and Space Science Proceedings},
       volume = {16},
        month = jan,
        pages = {379},
          doi = {10.1007/978-3-642-10352-0_8},
archivePrefix = {arXiv},
       eprint = {1010.2746},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2010ASSP...16..379L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

_citation_atomic_weights = r"""
@article{loss2003atomic,
  title={Atomic weights of the elements 2001 (IUPAC Technical Report)},
  author={Loss, RD},
  journal={Pure and Applied Chemistry},
  volume={75},
  number={8},
  pages={1107--1122},
  year={2003},
  publisher={De Gruyter}
}
"""

# fmt: off
elements = (
    "H", "He",
    "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar",
    "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce",
    "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" ,
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es",)

# index of each element in the pattern data array
elements_dict = {el:i for i, el in enumerate(elements)}

# Array of atomic weights
# From Loss, R. D. 2003, Pure Appl. Chem., 75, 1107, "Atomic Weights
# of the Elements 2001". The isotopic composition of the Sun and stars
# may differ from the terrestrial values listed here!
# For radionuclides that appear only in Table 3, we have adopted the
# atomic mass of the longest lived isotope, rounding to at most two
# digits after the decimal point. These elements are: 43_Tc_98,
# 61_Pm_145, 85_At_210, 86_Rn_222, 87_Fr_223, 88_Ra_226, 89_Ac_227,
# 93_Np_237, 94_Pu_244, 95_Am_243, 96_Cm_247, 97_Bk_247, 98_Cf_251,
# 99_Es_252

_atom_weight = (
    1.00794, 4.002602,
    6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    22.989770, 24.3050, 26.981538, 28.0855, 30.973761, 32.065, 35.453, 39.948,
    39.0983, 40.078, 44.955910, 47.867, 50.9415, 51.9961, 54.938049, 55.845,
    58.933200, 58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160, 78.96,
    79.904, 83.798, 85.4678, 87.62, 88.90585, 91.224, 92.90638,
    95.94, 97.91, 95.94, 101.07, 102.90550, 106.42, 107.8682, 112.411,
    114.818, 118.710, 121.760, 127.60, 126.90447, 131.293, 132.90545, 137.327,
    138.9055, 140.116, 140.90765, 144.24, 144.91, 150.36, 151.964, 157.25,
    158.92534, 162.500, 164.93032, 167.259, 168.93421, 173.04, 174.967, 178.49,
    180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.966, 200.59,
    204.3833, 207.2, 208.98038, 209.99, 222.02, 223.02, 226.03, 227.03,
    232.0381, 231.03588, 238.02891, 237.05, 244.06, 243.06, 247.07, 247.07,
    251.08, 252.08,)

# Anders, Grevesse, 1989, GeCo Acta  53, 197
_anders1989 = (
    12.00, 10.99,
    1.16, 1.15, 2.60, 8.56, 8.05, 8.93, 4.56, 8.09,
    6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.50, 6.56,
    5.12, 6.36, 3.10, 4.99, 4.00, 5.67, 5.39, 7.67,
    4.92, 6.25, 4.21, 4.60, 2.88, 3.41, None, None,
    None, None, 2.60, 2.90, 2.24, 2.60, 1.42, 1.92,
    None, 1.84, 1.12, 1.69, 0.94, 1.86, 1.66, 2.00,
    1.00, None, None, None, None, 2.13, 1.22, 1.55,
    0.71, 1.50, None, 1.00, 0.51, 1.12, -0.10, 1.10,
    0.26, 0.93, 0.00, 1.08, 0.76, 0.88, None, 1.11,
    None, 1.45, 1.35, 1.80, 1.01, None, 0.90, 1.85,
    None, None, None, None, None, None, None, 0.12,
    None, None, None, None, None, None, None, None,
    None,
)

# Grevesse, Noels, Sauval, 1996, ASPC 99, 177
_grevesse1996 = (
    12.00, 10.99,
    1.16, 1.15, 2.60, 8.55, 7.97, 8.87, 4.56, 8.08,
    6.33, 7.58, 6.47, 7.55, 5.45, 7.33, 5.50, 6.52,
    5.12, 6.36, 3.17, 5.02, 4.00, 5.67, 5.39, 7.50,
    4.92, 6.25, 4.21, 4.60, 2.88, 3.41, None, None,
    None, None, 2.60, 2.97, 2.24, 2.60, 1.42, 1.92,
    None, 1.84, 1.12, 1.69, 0.94, 1.77, 1.66, 2.00,
    1.00, None, None, None, None, 2.13, 1.17, 1.58,
    0.71, 1.50, None, 1.01, 0.51, 1.12, -0.10, 1.14,
    0.26, 0.93, 0.00, 1.08, 0.76, 0.88, None, 1.11,
    None, 1.45, 1.35, 1.80, 1.01, None, 0.90, 1.95,
    None, None, None, None, None, None, None, None,
    None, None, None, None, None, None, None, None,
    None,
)

# Grevesse, Sauval, 1998, SSR 85, 161
_grevesse1998 = (
    12.00, 10.93,
    1.10, 1.40, 2.55, 8.56, 7.92, 8.83, 4.56, 8.08,
    6.33, 7.58, 6.47, 7.55, 5.45, 7.33, 5.50, 6.40,
    5.12, 6.36, 3.17, 5.02, 4.00, 5.67, 5.39, 7.50,
    4.92, 6.25, 4.21, 4.60, 2.88, 3.41, None, None,
    None, None, 2.60, 2.97, 2.24, 2.60, 1.42, 1.92,
    None, 1.84, 1.12, 1.69, 0.94, 1.77, 1.66, 2.00,
    1.00, None, None, None, None, 2.13, 1.17, 1.58,
    0.71, 1.50, None, 1.01, 0.51, 1.12, -0.10, 1.14,
    0.26, 0.93, 0.00, 1.08, 0.06, 0.88, None, 1.11,
    None, 1.45, 1.35, 1.80, 1.01, None, 0.90, 1.95,
    None, None, None, None, None, None, None, None,
    None, None, None, None, None, None, None, None,
    None,
)

# Asplund, Grevesse, Sauval, 2005, ASPC 336, 252
_asplund2005 = (
    12.00, 10.93,
    1.05, 1.38, 2.70, 8.39, 7.78, 8.66, 4.56, 7.84,
    6.17, 7.53, 6.37, 7.51, 5.36, 7.14, 5.50, 6.18,
    5.08, 6.31, 3.05, 4.90, 4.00, 5.64, 5.39, 7.45,
    4.92, 6.23, 4.21, 4.60, 2.88, 3.58, None, None,
    None, 3.28, 2.60, 2.92, 2.21, 2.59, 1.42, 1.92,
    None, 1.84, 1.12, 1.69, 0.94, 1.77, 1.60, 2.00,
    1.00, None, None, 2.27, None, 2.17, 1.13, 1.58,
    0.71, 1.45, None, 1.01, 0.52, 1.12, 0.28, 1.14,
    0.51, 0.93, 0.00, 1.08, 0.06, 0.88, None, 1.11,
    None, 1.45, 1.38, None, 1.01, None, 0.90, 2.00,
    None, None, None, None, None, None, None, None,
    None, None, None, None, None, None, None, None,
    None,
)

# Grevesse, Asplund, Sauval (2007, Space Science Review, 130, 105)
_grevesse2007 = (
    12.00, 10.93,
    1.05, 1.38, 2.70, 8.39, 7.78, 8.66, 4.56, 7.84,
    6.17, 7.53, 6.37, 7.51, 5.36, 7.14, 5.50, 6.18,
    5.08, 6.31, 3.17, 4.90, 4.00, 5.64, 5.39, 7.45,
    4.92, 6.23, 4.21, 4.60, 2.88, 3.58, 2.29, 3.33,
    2.56, 3.25, 2.60, 2.92, 2.21, 2.58, 1.42, 1.92,
    None, 1.84, 1.12, 1.66, 0.94, 1.77, 1.60, 2.00,
    1.00, 2.19, 1.51, 2.24, 1.07, 2.17, 1.13, 1.70,
    0.58, 1.45, None, 1.00, 0.52, 1.11, 0.28, 1.14,
    0.51, 0.93, 0.00, 1.08, 0.06, 0.88, -0.17, 1.11,
    0.23, 1.25, 1.38, 1.64, 1.01, 1.13, 0.90, 2.00,
    0.65, None, None, None, None, None, None, 0.06,
    None, -0.52, None, None, None, None, None, None,
    None,
    )

# Asplund, Grevesse, Sauval, Scott (2009,  Annual Review of Astronomy
# and Astrophysics, 47, 481)
_asplund2009 = (
    12.00, 10.93,
    1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56,  7.93,
    6.24,  7.60,  6.45,  7.51,  5.41,  7.12,  5.50,  6.40,
    5.03,  6.34,  3.15,  4.95,  3.93,  5.64,  5.43,  7.50,
    4.99,  6.22,  4.19,  4.56,  3.04,  3.65,  2.30,  3.34,
    2.54,  3.25,  2.52,  2.87,  2.21,  2.58,  1.46,  1.88,
    None,  1.75,  0.91,  1.57,  0.94,  1.71,  0.80,  2.04,
    1.01,  2.18,  1.55,  2.24,  1.08,  2.18,  1.10,  1.58,
    0.72,  1.42,  None,  0.96,  0.52,  1.07,  0.30,  1.10,
    0.48,  0.92,  0.10,  0.84,  0.10,  0.85,  -0.12,  0.85,
    0.26,  1.40,  1.38,  1.62,  0.92,  1.17,  0.90,  1.75,
    0.65,  None,  None,  None,  None,  None,  None,  0.02,
    None,  -0.54,  None,  None,  None,  None,  None,  None,
    None,
    )

# Asplund, Grevesse, Sauval, Scott (2009,  Annual Review of Astronomy
# and Astrophysics, 47, 481)
_asplund2021 = (
    12.00, 10.914,
    0.96,  1.38,  2.70,  8.46,  7.83,  8.69,  4.40,  8.06,
    6.22,  7.55,  6.43,  7.51,  5.41,  7.12,  5.31,  6.38,
    5.07,  6.30,  3.14,  4.97,  3.90,  5.62,  5.42,  7.46,
    4.94,  6.20,  4.18,  4.56,  3.02,  3.62,  2.30,  3.34,
    2.54,  3.12,  2.32,  2.83,  2.21,  2.59,  1.47,  1.88,
    None,  1.75,  0.78,  1.57,  0.96,  1.71,  0.80,  2.02,
    1.01,  2.18,  1.55,  2.22,  1.08,  2.17,  1.11,  1.58,
    0.75,  1.42,  None,  0.95,  0.52,  1.08,  0.31,  1.10,
    0.48,  0.93,  0.10,  0.85,  0.10,  0.85,  -0.15,  0.79,
    0.26,  1.35,  1.32,  1.61,  0.91,  1.17,  0.92,  1.95,
    0.65,  None,  None,  None,  None,  None,  None,  0.03,
    None,  -0.54,  None,  None,  None,  None,  None,  None,
    None,
    )

# Lodders 2003 (ApJ, 591, 1220)
_lodders2003 = (
    12, 10.899,
    3.28, 1.41, 2.78, 8.39, 7.83, 8.69, 4.46, 7.87,
    6.30, 7.55, 6.46, 7.54, 5.46, 7.19, 5.26, 6.55,
    5.11, 6.34, 3.07, 4.92, 4.00, 5.65, 5.50, 7.47,
    4.91, 6.22, 4.26, 4.63, 3.10, 3.62, 2.32, 3.36,
    2.59, 3.28, 2.36, 2.91, 2.20, 2.60, 1.42, 1.96,
    1.82, 1.11, 1.70, 1.23, 1.74, 0.80, 2.11, 1.06,
    2.22, 1.54, 2.27, 1.10, 2.18, 1.18, 1.61, 0.75,
    1.46, 0.95, 0.85, 1.06, 0.31, 1.13, 0.49, 0.95,
    0.11, 0.94, 0.09, 0.77, -0.14, 0.65, 0.26, 1.37,
    1.35, 1.67, 0.72, 1.16, 0.81, 2.05, 0.68, 0.09,
    None, -0.49, None, None, None, None, None, None,
    None, None, None, None, None, None, None, None,
    None,
    )

# Lodders 2010, ASSP 16, 379
_lodders2010 = (
    12.00, 10.93,
    3.28, 1.32, 2.81, 8.39, 7.86, 8.73, 4.44, 8.05,
    6.29, 7.54, 6.46, 7.53, 5.45, 7.16, 5.25, 6.50,
    5.11, 6.31, 3.07, 4.93, 3.99, 5.65, 5.50, 7.46,
    4.90, 6.22, 4.27, 4.65, 3.10, 3.59, 2.32, 3.36,
    2.56, 3.28, 2.38, 2.90, 2.20, 2.57, 1.42, 1.94,
    None, 1.78, 1.10, 1.67, 1.22, 1.73, 0.78, 2.09,
    1.03, 2.20, 1.57, 2.27, 1.10, 2.18, 1.19, 1.60,
    0.77, 1.47, None, 0.96, 0.53, 1.09, 0.34, 1.14,
    0.49, 0.95, 0.14, 0.94, 0.11, 0.73, -0.14, 0.67,
    0.28, 1.37, 1.36, 1.64, 0.82, 1.19, 0.79, 2.06,
    0.67, None, None, None, None, None, None, 0.08,
    None, -0.52, None, None, None, None, None, None,
    None,
)


class Abund(IPersist):
    """Elemental abundance data and methods.
    Valid abundance pattern types are:

    'sme' - For hydrogen, the abundance value is the fraction of all
        nuclei that are hydrogen, including all ionization states
        and treating molecules as constituent atoms. For the other
        elements, the abundance values are log10 of the fraction of
        nuclei of each element in any form relative to the total for
        all elements in any form. For the Sun, the abundance values
        of H, He, and Li are approximately 0.92, -1.11, and -11.0.
    'n/nTot' - Abundance values are the fraction of nuclei
        of each element in any form relative to the total for all
        elements in any form. For the Sun, the abundance values of
        H, He, and Li are approximately 0.92, 0.078, and 1.03e-11.
    'n/nH' - Abundance values are the fraction of nuclei
        of each element in any form relative to the number of
        hydrogen nuclei in any form. For the Sun, the abundance
        values of H, He, and Li are approximately 1, 0.085, and
        1.12e-11.
    'H=12' - Abundance values are log10 of the fraction of nuclei of
        each element in any form relative to the number of hydrogen
        in any form plus an offset of 12. For the Sun, the nuclei
        abundance values of H, He, and Li are approximately 12,
        10.9, and 1.05.
    """

    def __init__(
        self,
        monh=0,
        pattern="solar",
        type="sme",
        citation_info=_citation_atomic_weights,
    ):
        self.monh = monh
        # The internal type is fixed to this value
        self._type_internal = "H=12"
        self.type = type
        self.citation_info = citation_info
        if isinstance(pattern, str):
            self.set_pattern_by_name(pattern)
        else:
            self.set_pattern_by_value(pattern, type)

    def __call__(self, type="H=12", raw=False):
        """Return abundances for all elements.
        Apply current [M/H] value to the current abundance pattern.
        Transform the resulting abundances to the requested abundance type.
        """
        pattern = np.copy(self._pattern)
        if self.monh is not None:
            pattern[2:] += self.monh
        return self.totype(pattern, type, raw=raw, copy=False)

    def __getitem__(self, elem):
        return self.get_element(elem)

    def __setitem__(self, elem, abund):
        self.update_pattern({elem: abund})

    def __str__(self):
        if self._pattern is None:
            if self._monh is None:
                return "[M/H] is not set. Abundance pattern is not set."
            else:
                return f"[M/H]={self._monh:.3f}. Abundance pattern is not set"
        else:
            a = self.get_pattern("H=12", raw=True)
            if self._monh is None:
                out = "[M/H] is not set. Values below are the abundance pattern.\n"
            else:
                a = np.copy(a)
                a[2:] += self._monh
                out = (
                    f" [M/H]={self._monh:.3f} applied to abundance pattern. "
                    "Values below are abundances.\n"
                )
            for i in range(9):
                for j in range(11):
                    out += "  {:<5s}".format(elements[11 * i + j])
                out += "\n"
                for j in range(11):
                    out += "{:7.3f}".format(a[11 * i + j])
                if i < 8:
                    out += "\n"
            return out

    def __repr__(self):
        return self.__str__()

    # These are there for the atmosphere interpolation
    # The internal format is always H=12, so that should be fine
    def __add__(self, other):
        result = Abund(
            monh=self.monh, pattern=np.copy(self._pattern), type=self._type_internal
        )
        if isinstance(other, Abund):
            result._pattern += other.get_pattern(self._type_internal, raw=True)
            result.monh += other.monh
        else:
            raise NotImplementedError
        result.type = self.type
        return result

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        result = Abund(
            monh=self.monh, pattern=np.copy(self._pattern), type=self._type_internal
        )
        if isinstance(other, Abund):
            raise NotImplementedError
        else:
            result._pattern *= other
            result.monh *= other
        result.type = self.type
        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def _get_index(self, elem):
        try:
            return elements_dict[elem]
        except KeyError:
            raise KeyError(
                f"Got element abreviation {elem}, expected one of {elements}"
            )

    def __copy__(self):
        pattern = self._pattern.copy()
        other = Abund(
            monh=self.monh,
            type=self.type,
            pattern=pattern,
            citation_info=self.citation_info,
        )
        return other

    @staticmethod
    def fromtype(pattern, fromtype, raw=False):
        """Return a copy of the input abundance pattern, transformed from
        the input type to the 'H=12' type. Valid abundance pattern types
        are 'sme', 'n/nTot', 'n/nH', 'n/nFe', and 'H=12'.
        """
        elem = elements

        if isinstance(pattern, dict):
            abund = np.array([pattern[el] for el in elem], dtype=float)
        else:
            abund = np.array(pattern, dtype=float)

        if np.isnan(abund[0]):
            raise ValueError("Pattern must define abundance of H")

        # print(abund)

        type = fromtype.lower()
        if type == "h=12":
            pass
        elif type == "sme":
            # sme -> H=12
            abund[1:] += 12 - np.log10(abund[0])
            abund[0] = 12
        elif type == "kurucz":
            # kurucz -> H=12
            abund[2:] = abund[2:] + np.log10(abund[1] / abund[0] + 1) + 12
            abund[1] = np.log10(abund[1]/abund[0]) + 12
            abund[0] = 12
        elif type == "n/ntot":
            abund /= abund[0]
            abund = 12 + np.log10(abund)
        elif type == "n/nh":
            abund = 12 + np.log10(abund)
        elif type == "n/nfe":
            abund /= abund[0]
            abund = 12 + np.log10(abund)
        elif type == "fe=12":
            abund = 10 ** (abund - 12)
            abund /= abund[0]
            abund = 12 + np.log10(abund)
        else:
            raise ValueError(
                "got abundance type '{}',".format(type)
                + " should be 'H=12', 'n/nH', 'n/nTot', 'n/nFe', 'Fe=12', or 'sme'"
            )

        if raw:
            return abund
        else:
            return {el: abund[elements_dict[el]] for el in elements}
        
    @staticmethod
    def totype(pattern, totype, raw=False, copy=True, X=None):
        """Return a copy of the input abundance pattern, transformed from
        the 'H=12' type to the output type. Valid abundance pattern types
        are 'sme', 'kurucz', 'n/nTot', 'n/nH', and 'H=12'.
        """
        if isinstance(pattern, dict):
            abund = [pattern[el] if el in pattern.keys() else np.nan for el in elements]
            abund = np.array(abund)
        elif copy:
            abund = np.copy(pattern)
        else:
            abund = pattern

        type = totype.lower()

        if type == "h=12":
            pass
        elif type == "sme":
            abund[0] = 1 / (1 + np.nansum(10**(abund[1:]-12)))
            abund[1:] = np.log10(abund[0] * 10**(abund[1:] - 12))
        elif type == "kurucz":
            # H=12 -> kurucz
            abund[0] = 1 / (1 + np.nansum(10**(abund[1:] - 12)))
            abund[1] = 1 / (10**(12-abund[1]) + 1 + 10**(12-abund[1]) * np.nansum(10**(abund[2:]-12)))
            abund[2:] = abund[2:] - 12 - np.log10(1 + 10**(abund[1]-12))
        elif type == "n/ntot":
            abund = 10 ** (abund - 12)
            abund /= np.nansum(abund)
        elif type == "n/nh":
            abund = 10 ** (abund - 12)
        elif type == "n/nfe":
            idx_fe = elements_dict["Fe"]
            abund = 10 ** (abund - 12)
            abund /= abund[idx_fe]
        elif type == "fe=12":
            idx_fe = elements_dict["Fe"]
            abund = 10 ** (abund - 12)
            abund /= abund[idx_fe]
            abund = np.log10(abund) + 12
        else:
            raise ValueError(
                "got abundance type '{}',".format(type)
                + " should be 'H=12', 'kurucz', 'n/nH', 'n/nTot', 'n/nFe', 'Fe=12', or 'sme'"
            )

        if raw:
            return abund
        else:
            return {el: abund[elements_dict[el]] for el in elements}

    _formats = ["H=12", "sme", "n/nTot", "n/nH", "n/nFe", "Fe=12"]

    @property
    def elem(self):
        """Return the standard abbreviation for each element.
        Use property so user will not redefine elements.
        """
        return elements

    @property
    def elem_dict(self):
        """Return the position of each element in the raw array"""
        return elements_dict

    @property
    def monh(self):
        """float: Metallicity, the logarithmic offset added to
        the abundance pattern for all elements except hydrogen and helium."""
        return self._monh

    @monh.setter
    def monh(self, monh):
        self._monh = float(monh)

    @property
    def pattern(self):
        """array: Abundance pattern in the initial format"""
        return self.get_pattern(self.type, raw=False)

    def set_pattern_by_name(self, pattern_name):
        """Set the abundance pattern to one of the predefined options

        Parameters
        ----------
        pattern_name : str
            Name of the predefined option to use. One of:
            'grevesse2007', 'asplund2009', 'asplund2021', 'asplund2005',
            'lodders2003', 'anders1989', 'grevesse1996', 'grevesse1998', ,'lodders2010',
            'empty'

        Raises
        ------
        ValueError
            If an undefined pattern_name was given
        """
        self.type = "H=12"
        self.citation_info = _citation_atomic_weights + "\n"
        if pattern_name.lower() in ["anders1989"]:
            self._pattern = np.array(_anders1989, dtype=float)
            self.citation_info += _citation_anders1989
        elif pattern_name.lower() in ["grevesse1996"]:
            self._pattern = np.array(_grevesse1996, dtype=float)
            self.citation_info += _citation_grevesse1996
        elif pattern_name.lower() in ["grevesse1998"]:
            self._pattern = np.array(_grevesse1998, dtype=float)
            self.citation_info += _citation_grevesse1998
        elif pattern_name.lower() in ["asplund2005"]:
            self._pattern = np.array(_asplund2005, dtype=float)
            self.citation_info += _citation_asplund2005
        elif pattern_name.lower() in ["grevesse2007", "solar"]:
            self._pattern = np.array(_grevesse2007, dtype=float)
            self.citation_info += _citation_grevesse2007
        elif pattern_name.lower() in ["asplund2009"]:
            self._pattern = np.array(_asplund2009, dtype=float)
            self.citation_info += _citation_asplund2009
        elif pattern_name.lower() in ["asplund2021"]:
            self._pattern = np.array(_asplund2021, dtype=float)
            self.citation_info += _citation_asplund2021
        elif pattern_name.lower() == "lodders2003":
            self._pattern = np.array(_lodders2003, dtype=float)
            self.citation_info += _citation_lodders2003
        elif pattern_name.lower() == "lodders2010":
            self._pattern = np.array(_lodders2010, dtype=float)
            self.citation_info += _citation_lodders2010
        elif pattern_name.lower() == "empty":
            self._pattern = self.empty_pattern()
        else:
            raise ValueError(
                f"Got abundance pattern name {pattern_name}. It should be one of: "
                "'grevesse2007' (or 'solar'), 'asplund2009', 'asplund2021', "
                "'asplund2005', 'lodders2003', 'lodders2010', 'anders1989', "
                "'grevesse1996', 'grevesse1998', or 'empty'."
            )

    def set_pattern_by_value(self, pattern, type):
        self._pattern = self.fromtype(pattern, type, raw=True)
        self.type = type

    def update_pattern(self, updates):
        """Update the abundance pattern for several elements at once

        The abundance is first converted into the initially specified format,
        before being converted back to the internal format

        Parameters
        ----------
        updates : dict{str:float}
            the elements to update

        Raises
        ------
        KeyError
            If any of the element keys is not valid
        """
        pattern = self.get_pattern(type=self.type, raw=True)
        for key in updates:
            pos = self._get_index(key)
            pattern[pos] = updates[key]
        self.set_pattern_by_value(pattern, self.type)

    def get_pattern(self, type=None, raw=False):
        """
        Transform the specified abundance pattern from type used
        internally by SME to the requested type. Valid abundance pattern
        types are: 'sme', 'n/nTot', 'n/nH', 'H=12'

        Parameters
        ----------
        type : str
            The pattern format to retrieve the pattern as. Defaults to the
            original input format.
        raw : bool
            If True will return the pattern as a numpy array, with indices as defined in elem dict.
            Otherwise the return value is a dictionary, with the elements as keys.
        """
        if type is None:
            type = self.type
        return self.totype(self._pattern, type, raw=raw)

    def get_element(self, elem, type=None):
        if type is None:
            type = self.type
        pattern = self(type, raw=True)
        i = self._get_index(elem)
        return pattern[i]

    def empty_pattern(self):
        """Return an abundance pattern with value None for all elements."""
        pattern = np.full(len(elements), np.nan)
        pattern[0] = 0
        return pattern

    def to_dict(self):
        data = {
            "monh": self.monh,
            "type_internal": self._type_internal,
            "type": self.type,
            "citation_info": self.citation_info,
            "data": self._pattern,
        }
        return data

    @classmethod
    def from_dict(cls, data):
        obj = cls(
            monh=data["monh"],
            pattern=data["data"],
            type=data["type_internal"],
            citation_info=data["citation_info"],
        )
        obj.type = data["type"]
        return obj

    def _save(self):
        header = {
            "monh": self.monh,
            "type_internal": self._type_internal,
            "type": self.type,
            "citation_info": self.citation_info,
        }
        data = self._pattern
        ext = BinaryDataExtension(header, data)
        return ext

    @classmethod
    def _load(cls, ext: BinaryDataExtension):
        header = ext.header
        pattern = ext.data
        abund = cls(
            monh=header["monh"],
            pattern=pattern,
            type=header["type"],
            citation_info=header["citation_info"],
        )
        return abund

    def _save_v1(self, file, folder="abund"):
        """Save the data to a file handler"""
        if folder != "" or folder[-1] != "/":
            folder += "/"

        monh = float(self.monh) if self.monh is not None else None
        info = {
            "format": self.type,
            "monh": monh,
            "citation_info": self.citation_info,
        }
        file.writestr(f"{folder}info.json", json.dumps(info))

        b = io.BytesIO()
        np.save(b, self.get_pattern(raw=True))
        file.writestr(f"{folder}pattern.npy", b.getvalue())

    @staticmethod
    def _load_v1(file, names, folder=""):
        """Load the data from a file handler"""
        for name in names:
            if name.endswith("info.json"):
                info = file.read(name)
                info = json.loads(info)

                abund_format = info["format"]
                monh = info["monh"]
                citation_info = info.get("citation_info", _citation_atomic_weights)

            elif name.endswith("pattern.npy"):
                b = io.BytesIO(file.read(name))
                pattern = np.load(b)

        return Abund(
            monh=monh,
            pattern=pattern,
            type=abund_format,
            citation_info=citation_info,
        )

    @staticmethod
    def solar():
        """Return solar abundances of Grevesse 2007"""
        solar = Abund(pattern="solar", monh=0)
        return solar
