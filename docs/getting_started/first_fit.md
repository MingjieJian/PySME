# Fit an observed spectrum

Assuming that we have an observed spectrum, with its wavelength, normalized flux, and uncertainties, as an array of `wave`, `flux` and `uncertainties`.

They can be inserted into the SME structure with:
```py
sme = SME_Structure()
sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini = 5777, 4.4, 0, 1.09, 4.19, 1.60
sme.abund = Abund.solar()
sme.linelist = vald
sme.iptype = 'gauss'
sme.ipres = 42000
sme.wave = wave
sme.spec = flux
sme.uncs = uncertainties
```

The new inputs are the [instrument resolution](../concepts/sme_struct.md#instrument-parameters) and [observed spectra](../concepts/sme_struct.md#spectra).

Then the `solve` function can be used to find the best fit solution:
```py
from pysme.solve import solve
fitparameters = ["teff", "logg", "monh", "abund Mg"]
sme = solve(sme, fitparameters)
```

The [fitresults](../concepts/sme_struct.md#fitresults) are stored in `sme.fitresults`.