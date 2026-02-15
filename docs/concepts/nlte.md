# NLTE grids

Non Local Thermal Equilibrium (NLTE) calculations are important to accuarately fit certain lines.
PySME supports them using pre-computed grids of NLTE departure coefficients, which need to be created for every element. 
For common elements PySME provides grids (see below) via the LFS (see [](lfs.md)). 
If any of these grids are used, please kindly take care to cite the papers describing the NLTE models and departure coefficient calculations.

NLTE calculations need to be specified for each element they are supposed to be used for individually using `sme.nlte.set_nlte(el, grid)` (the `grid` can be omitted if there is a grid in lfs).
Similarly they can be disabled for each element using `sme.nlte.remove_nlte(el)`, where sme is your SME structure.
If no element is set to NLTE in the structure PySME will perform
LTE calculations only.

## Fields of NLTE object

- `elements`: The elements for which NLTE has been activated
- `grids`: The grid file that is used for each active element
- `subgrid_size`: A small segment of the NLTE grid will be cached in memory
    to speed up calculations. This sets the size of that cache
    by defining the number of points in each
    axis (rabund, teff, logg, monh).
- `flags`: After the synthesis all lines are flaged if they used NLTE

## Grid interpolation

The grid has 6 dimensions.

- teff: Effective Temperature
- logg: Surface Gravity
- monh: Overall Metallicity
- rabund: relative abundance of that element
- depth: optical depth in the atmosphere
- departure coefficients:
    The NLTE departure coefficients describing how much
    it varies from the LTE calculation

We then perform linear interpolation to the stellar parameters
we want to model. And we then perform a cubic spline fit to the depth scale
of the model atmosphere we specified (See [](atmosphere)).

We then use the linelist to find only the relevant transitions in the grid,
and pass the departure coefficients for each line to the C library.

## NLTE flags in line list

PySME provides information on whehter a line is synthesized in NLTE through the `nlte_flag` column in the line list.

- `1`: this line was synthesized with in NLTE.
- `0`: this line was synthesized in LTE.
- `-1`: this line was not included in the current synthesis pass (see [](../advance/line_filtering.md)).

## Recommended and default grids

| Element | Grid name | Zenodo version | Citation |
|:---:|:---:|:---:|:---:|
| H | `nlte_H_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Li | `nlte_Li_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| C | `nlte_C_pysme.grd` | 2 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| N | `nlte_N_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| O | `nlte_O_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Na | `nlte_Na_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Mg | `nlte_Mg_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Al | `nlte_Al_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Si | `nlte_Si_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| S | `nlte_Si_pysme.grd` | 7 | [Amarsi et al. 2025](https://ui.adsabs.harvard.edu/abs/2025A%26A...703A..35A/abstract) |
| K | `nlte_K_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Ca | `nlte_Ca_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Ti | `nlte_Ti_pysme.grd` | 5 | [Mallinson et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...687A...5M) |
| Mn | `nlte_Mn_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Fe | `nlte_Fe_pysme.grd` | 4 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |
| Cu | `nlte_Cu_pysme.grd` | 6 | [Caliskan et al. 2025](https://ui.adsabs.harvard.edu/abs/2025A%26A...696A.210C/abstract) |
| Ba | `nlte_Ba_pysme.grd` | 3 | [Amarsi et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..62A) |

All recommended/default grids are available in [Zenodo](https://zenodo.org/records/17064337).

## Deprecated grids

  - H 
    - marcs2012_H2018.grd
  - Li
    - marcs2012_Li.grd [(Lind et al. 2009)](https://ui.adsabs.harvard.edu/abs/2009A%26A...503..541L)
    - marcs2012_Li2009.grd [(Lind et al. 2009)](https://ui.adsabs.harvard.edu/abs/2009A%26A...503..541L)
    - nlte_Li_multi.grd
  - C
    - marcs2012_C.grd
  - N 
    - marcs2012_N.grd
  - O
    - marcs2012p_t1.0_O.grd [(Sitnova et al. 2013)](https://ui.adsabs.harvard.edu/abs/2013AstL...39..126S)
    - marcs2012_O2015.grd [(Amarsi et al. 2016b)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.3735A)
    - marcs2012_O.grd
  - Na
    - marcs2012p_t1.0_Na.grd [(Mashonkina et al. 2001)](https://ui.adsabs.harvard.edu/abs/2000ARep...44..790M)
    - marcs2012_Na.grd [(Lind et al. 2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...528A.103L)
    - marcs2012_Na2011.grd [(Lind et al. 2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...528A.103L)
    - nlte_Na_multi_full.grd
    - nlte_Na_multi_sun.grd
  - Al 
    - marcs2012_Al.grd
    - marcs2012_Al2017.grd
  - Mg
    - marcs2012_Mg2016.grd [(Osorio et al. 2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...586A.120O)
    - marcs2012_Mg.grd
  - Si
    - marcs2012_Si2016.grd [(Amarsi & Asplund 2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464..264A)
    - marcs2012_Si.grd
  - K
    - marcs2012_K.grd
  - Ca
    - marcs2012s_t2.0_Ca.grd [(Mashonkina et al. 2007)](https://ui.adsabs.harvard.edu/abs/2007A%26A...461..261M)
    - marcs2012p_t1.0_Ca.grd [(Mashonkina et al. 2017)](https://ui.adsabs.harvard.edu/abs/2007A%26A...461..261M)
    - marcs2012_Ca.grd
  - Ti
    - marcs2012s_t2.0_Ti.grd [(Sitnova et al. 2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.461.1000S)
  - Mn
    - marcs2012_Mn.grd
  - Fe
    - marcs2012_Fe2016.grd [(Amarsi et al. 2016a)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.1518A)
    - marcs2012s_t2.0_Fe.grd [(Mashonkina et al. 2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...528A..87M)
    - nlte_Fe_multi_full.grd
    - marcs2012_Fe.grd
  - Ba
    - marcs2012p_t1.0_Ba.grd [(Mashonkina et al. 1999)](https://ui.adsabs.harvard.edu/abs/1999A%26A...343..519M)
  - Eu
    - nlte_Eu.grd




