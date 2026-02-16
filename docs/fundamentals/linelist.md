# Linelist

`linelist` is one of the core inputs for spectral synthesis. It stores the
physical parameters of each transition (element/species, wavelength, `log gf`,
damping parameters, etc.).

## Getting a Linelist from VALD

PySME is designed to work well with all types of line lists in [VALD](https://vald.astro.uu.se/) database. 
A common workflow is:

1. Download a line list from VALD (typically from `extract stellar`).
2. Load it in PySME with `ValdFile`.

```py
from pysme.linelist.vald import ValdFile

vald = ValdFile("my_linelist.lin")
sme.linelist = vald
```

## Linelist is essentially a DataFrame

Internally, PySME stores `LineList` data in a `pandas.DataFrame`
(`sme.linelist._lines`).
So you can use DataFrame-style operations directly: filtering, sorting,
editing columns, adding flags, and so on.

```py

# Filter by wavelength range
sub = sme.linelist[(sme.linelist["wlcent"] > 6436.0) & (sme.linelist["wlcent"] < 6444.0)]

# Select one species (for example, Fe I)
fe1 = sme.linelist[sme.linelist["species"] == "Fe 1"]
```

## Short vs Long Format

VALD provides two line list formats:

- short: fewer fields, suitable for LTE.
- long: includes additional upper/lower level and quantum-number information.

If you plan to run NLTE, you should use **VALD long format**.
In PySME, NLTE relies on the extra level information available in long format;
short format does not provide enough information for that workflow.

## Line parameters

The short format fields are

:`species`:
    A string identifier including the element
    and ionization state or the molecule
:`wlcent`: The central wavelength of the line in Angstrom
:`gflog`:
    The log of the product of the statistical weight of
    the lower level and the oscillator strength for the transition.
:`excit`: The excitation energy in the lower level
:`ionization`: The ionization state of the species, where 1 is neutral
:`gamrad`: The radiation broadening parameter
:`gamqst`: Stark broadening parameter
:`gamvw`: van der Waals broadening parameter
:`lande`: The lande factor
:`depth`: An arbitrary depth estimation of the line
:`reference`: A citation where this data came from
:`atom_number`:
    Identifies the species by the atomic number
    (i.e. the number of protons)

In addition the long format has the following fields

:`lande_lower`: The lower Lande factor
:`lande_upper`: The upper Lande factor
:`j_lo`: The spin of the lower level
:`j_up`: The spin of the upper level
:`e_upp`: The energy of the upper level
:`term_lower`: The electron configuration of the lower level
:`term_upper`: The electron configuration of the upper level
:`error`: An uncertainty estimate for this linedata

### Runtime columns

:`nlte_flag`: Per-line NLTE usage flag written after synthesis. See [NLTE flags in line list](../concepts/nlte.md#nlte-flags-in-line-list).
:`central_depth`: Estimated central line depth used for line selection/filtering. See [Line Filtering](../advance/line_filtering.md).
:`line_range_s`: Start wavelength of the estimated line-contribution range. See [Line Filtering](../advance/line_filtering.md).
:`line_range_e`: End wavelength of the estimated line-contribution range. See [Line Filtering](../advance/line_filtering.md).


### Important Note

As far as the radiative transfer code is concerned the ionization is defined as part of the species term.
I.e. a line with species = "Fe 2" will be calculated as ionization = 2.
If no number is set within the species field, SME will use an ionization of 1.

Also atom_number is ignored in the radiative transfer calculations and therefore does not need to be set.

## Example Linelist

Below is a typical `sme.linelist` excerpt (middle rows omitted):

```text
   species     wlcent  gflog     excit  j_lo     e_upp  j_up  lande_lower  \
0     Fe 1  6436.4055 -2.460  4.186364   2.0  6.112128   1.0         0.68
1     Eu 2  6437.6106 -1.242  1.319612   5.0  3.245008   5.0         1.73
2     Eu 2  6437.6121 -1.280  1.319636   5.0  3.245027   5.0         1.73
3     Eu 2  6437.6135 -2.473  1.319612   5.0  3.245008   5.0         1.73
4     Eu 2  6437.6194 -2.511  1.319636   5.0  3.245027   5.0         1.73
..     ...        ...    ...       ...   ...       ...   ...          ...
62    Mn 1  6443.4689 -2.818  3.772307   1.5  5.695960   2.5         1.21
63    Mn 1  6443.4726 -3.538  3.772307   1.5  5.695960   2.5         1.21
64    Mn 1  6443.4821 -2.275  3.772307   1.5  5.695960   2.5         1.21
65    Mn 1  6443.4887 -2.964  3.772307   1.5  5.695960   2.5         1.21
66    Mn 1  6443.4938 -3.918  3.772307   1.5  5.695960   2.5         1.21

    lande_upper  lande  ...        term_lower  couple_upper  \
0          0.56   0.73  ...           3d8 c3F            LS
1          1.79   1.76  ...  4f7.(8S).5d a9D*            JJ
2          1.79   1.76  ...  4f7.(8S).5d a9D*            JJ
3          1.79   1.76  ...  4f7.(8S).5d a9D*            JJ
4          1.79   1.76  ...  4f7.(8S).5d a9D*            JJ
..          ...    ...  ...               ...           ...
62         1.37   1.49  ...       3d5.4s2 b4D            LS
63         1.37   1.49  ...       3d5.4s2 b4D            LS
64         1.37   1.49  ...       3d5.4s2 b4D            LS
65         1.37   1.49  ...       3d5.4s2 b4D            LS
66         1.37   1.49  ...       3d5.4s2 b4D            LS

                         term_upper  error nlte_flag atom_number ionization  \
0        3d6.(3F2).4s.4p.(3P*) v3D*   0.50         0         1.0        1.0
1   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)   1.00         0         1.0        2.0
2   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)   1.00         0         1.0        2.0
3   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)   1.00         0         1.0        2.0
4   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)   1.00         0         1.0        2.0
..                              ...    ...       ...         ...        ...
62                 3d6.(5D).4p z4D*   0.25         0         1.0        1.0
63                 3d6.(5D).4p z4D*   0.25         0         1.0        1.0
64                 3d6.(5D).4p z4D*   0.25         0         1.0        1.0
65                 3d6.(5D).4p z4D*   0.25         0         1.0        1.0
66                 3d6.(5D).4p z4D*   0.25         0         1.0        1.0
```
