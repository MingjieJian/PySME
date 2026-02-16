# Synthesize your first spectrum

Three parts of information is needed to synthesize a spectrum:
- Stellar atmosphere model;
- Stellar parameters;
- Line list.

The first two can be handeled in PySME:

```py
from pysme.sme import SME_Structure
from pysme.abund import Abund

sme = SME_Structure()
sme.teff, sme.logg, sme.monh = 5700, 4.4, 0
sme.abund = Abund.solar()
```

Here we defined effective temperature `teff`, surface gravity `logg`, metallicity `monh`, and solar [chemical abundance](../concepts/abundance.md) `abund`.
For a full list of stellar parameters, see [](../concepts/sme_struct.md).

LineList can be download from [VALD database](https://vald.astro.uu.se/).
Here we provide an example line list: [sun.lin](https://raw.githubusercontent.com/MingjieJian/PySME/master/examples/sun.lin).
Load the line list using:
```py
from pysme.linelist.vald import ValdFile
vald = ValdFile("linelist.lin")
sme.linelist = vald
```
Define wavelength grid or wavelength range of the synthetic spectrum:
```py
import numpy as np
sme.wave = np.arange(6436, 6440, 0.1)
# Or
sme.wran = [6436, 6440] # pysme will choose sampling automatically
```

Then use the `synthesize_spectrum` function:
```py
from pysme.synthesize import synthesize_spectrum
sme = synthesize_spectrum(sme)
```

The synthesized spectra are stored in `sme.wave` and `sme.synth`:
```py
Iliffe_vector([array([6436. , 6436.1, 6436.2, ..., 6439.8, 6439.9])]) 
Iliffe_vector([array([0.99986606, 0.99984309, 0.99980888, ..., 0.99754268, 0.99807395])])
```

and can be plot by (for example):

```py
import matplotlib.pyplot as plt

plt.plot(sme.wave[0], sme.synth[0])
```

`[0]` is required because PySME stores spectra as segments; here the result is recognized as segment 0.