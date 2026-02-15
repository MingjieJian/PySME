# Quickstart

## SME Structure

The first step in each SME project is to create an SME structure

```py    
from pysme.sme import SME_Structure
```

This can be done in done in a few different ways:
- assign values manually: `sme = SME_Structure()`
- load an existing SME save file (from Python or IDL): `sme = SME_Structure.load("sme.inp")`
- load an .ech file spectrum: `sme = SME_Structure.load("obs.ech")`


## Synthesize a spectrum

Once the SME structure is created, it can be used to synthesize a spectrum.
First we need to define necessary parameters:
* Stellar parameters (effective temperature `teff`, surface gravity `logg`, metallicity `monh`, chemical abundance `abund`)
```py
from pysme.abund import Abund
sme.teff, sme.logg, sme.monh = 5700, 4.4, -0.1
sme.abund = Abund.solar()
```
* LineList (linelist), e.g. from [VALD database](https://vald.astro.uu.se/)
```py
from pysme.linelist.vald import ValdFile
vald = ValdFile("linelist.lin")
sme.linelist = vald
```
* Wavelength grid or wavelength range
```py
import numpy as np
sme.wave = [np.arange(6436, 6440, 0.1), np.arange(6442, 6443, 0.1)]
# Or
sme.wran = [[6436, 6440], [6442, 6443]]
```
* Optional: provide the transfer grid used by SMElib
```py
# single segment
sme.wint = np.arange(6436, 6440, 0.02)
# multi-segment
sme.wint = [np.arange(6436, 6440, 0.02), np.arange(6442, 6443, 0.02)]
```

Then use the `synthesize_spectrum` function:
```py
from pysme.synthesize import synthesize_spectrum
sme = synthesize_spectrum(sme)
```

The synthesized spectra are stored in `sme.wave` and `sme.synth`:
```py
Iliffe_vector([array([6436. , 6436.1, 6436.2, ..., 6439.8, 6439.9]), 
               array([6442. , 6442.1, 6442.2, ... 6442.8, 6442.9])]) 
Iliffe_vector([array([0.99986606, 0.99984309, 0.99980888, ..., 0.99754268, 0.99807395]), 
               array([0.99983848, 0.99985176, 0.99986281, ..., 0.99979064, 0.99471954])])
```

## Fit an observed spectrum

Assuming that we have an observed spectrum, with its wavelength, normalized flux, and uncertainties, as an array of `wave`, `flux` and `uncertainties`.

They can be inserted into the SME structure with:
```py
sme.wave = wave
sme.spec = flux
sme.uncs = uncertainties
sme.mask = np.ones(len(wave), dtype=int)
```

- The wavelength is always given in Angstrom. 
- Note that the observation may be split into segments (orders etc). 
        - Then Wavelength is a list of arrays [segment1, segment2, ...], and the same applies to spec, uncs, and mask.
- Segment-aware inputs include: `wave`, `spec`, `uncs`, `mask`, and `wint`.
- Global (non-segmented) inputs include: `vmic`, `vmac`, and `vsini`.
- The mask values are: 0: bad pixel, 1: line pixel, 2: continuum pixel, 4: vrad pixel
        - The masks are additive, i.e., you can set mask value to 5 for line and vrad pixel.
        - Note that the dtype of sme.mask must be int.

Then the `solve` function can be used to find the best fit solution:
```py
from pysme.solve import solve
# for more details on the fitparameter option, see fitparameters
fitparameters = ["teff", "logg", "monh", "abund Mg"]
sme = solve(sme, fitparameters)
```

## NLTE correction

The non-Local Thermal Equilibrium (NLTE) correction can be applied to the spectrum by setting the `nlte` attribute of the SME structure:
```py
# SME also comes with a few NLTE grids, see NLTE section
# The current NLTE grid is only applicable to the MARCS model!
sme.nlte.set_nlte("Ca")
```
```{note}
Upon the first use of the NLTE correction, the NLTE grid will be downloaded from the server and this may takes a while, depending on the size of the NLTE grid.
```

The results in the output sme structure can, for example, be plotted using the gui module (under development).
```py
from gui import plot_plotly
fig = plot_plotly.FinalPlot(sme)
fig.save(filename="sme.html")
```

```{raw} html
<iframe src="../_static/sun.html"
        width="100%" height="480"
        style="border:none;"></iframe>
```

## Save and Load

The SME structure can be loaded with
```py
sme = SME_Structure.load("in.sme")
```
or saved with
```py
sme.save("out.sme")
```
