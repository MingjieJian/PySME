# Abundance

The abundance contains information about the overall
metallicity of the star as well as the individual
abundances of each element.
SME can simulate and fit the first 99 elements of
the periodic table (from Hydrogen H to Einsteinium Es).

## Abundance formats

PySME supports a variety of formats to input the abundance.
However they are all internally converted to the 'H=12'
type before use.

### H=12

Abundance values are log10 of the fraction of nuclei of
each element in any form relative to the number of hydrogen
in any form plus an offset of 12. For the Sun, the nuclei
abundance values of H, He, and Li are approximately 12,
10.9, and 1.05.

### n/nTot
Abundance values are the fraction of nuclei
of each element in any form relative to the total for all
elements in any form. For the Sun, the abundance values of
H, He, and Li are approximately 0.92, 0.078, and 1.03e-11.

### n/nH
Abundance values are the fraction of nuclei
of each element in any form relative to the number of
hydrogen nuclei in any form. For the Sun, the abundance
values of H, He, and Li are approximately 1, 0.085, and
1.12e-11.

### sme
For hydrogen, the abundance value is the fraction of all
nuclei that are hydrogen, including all ionization states
and treating molecules as constituent atoms. For the other
elements, the abundance values are log10 of the fraction of
nuclei of each element in any form relative to the total for
all elements in any form. For the Sun, the abundance values
of H, He, and Li are approximately 0.92, -1.11, and -11.0.

### Fe=12
Abundance values are log10 of the fraction of nuclei of
each element in any form relative to the number of iron
in any form plus an offset of 12. For the Sun, the nuclei
abundance values of H, He, and Li are approximately 16.5,
15.43, and 5.55.

### n/nFe
Abundance values are the fraction of nuclei
of each element in any form relative to the number of
iron nuclei in any form. For the Sun, the abundance
values of H, He, and Li are approximately 3.16e4, 2.69e3, and
3.55e-7.

## Solar metallicity

PySME contains multiple pre defined sets of solar abundances
for you to choose from. They are:

From solar photosphere:
- `anders1989`: [Anders & Grevesse (1989, GeCoA)](https://ui.adsabs.harvard.edu/abs/1989GeCoA..53..197A)
- `grevesse1996`: [Grevesse, Noels & Sauval (1996, ASPC)](https://ui.adsabs.harvard.edu/abs/1996ASPC...99..117G)
- `grevesse1998`: [Grevesse, & Sauval (1998, SSRv)](https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G) 
- `asplund2005`: [Asplund, Grevesse & Sauval (2005, ASPC)](https://ui.adsabs.harvard.edu/abs/2005ASPC..336...25A)
- `grevesse2007` (alias: `solar`): [Grevesse, Asplund & Sauval (2007, SSRv)](https://ui.adsabs.harvard.edu/abs/2007SSRv..130..105G)
- `asplund2009`: [Asplund, Grevesse & Sauval (2009, ARA&A)](https://ui.adsabs.harvard.edu/abs/2009ARA&A..47..481A)
- `asplund2021`: [Asplund, Amarsi & Grevesse (2021, A&A)](https://ui.adsabs.harvard.edu/abs/2021A&A...653A.141A)

From CI chondrites:
- `lodders2003`: [Lodders(2003 ApJ)](https://ui.adsabs.harvard.edu/abs/2003ApJ...591.1220L)
- `lodders2010`: [Lodders(2010 ASSP)](https://ui.adsabs.harvard.edu/abs/2010ASSP...16..379L)

They can be initialized by passing their name during the
creation of the abundance. E.g. Abund("grevesse2007").
The default solar abundance is 'grevesse2007' and is also
available using Abund.solar().

## Conversions

### sme -> H=12

We have:
- $S_\mathrm{H} = \frac{N_\mathrm{H}}{N_\mathrm{all}}$ for Hydrogen in `sme` scale;
- $S_\mathrm{X} = \log{\frac{N_\mathrm{X}}{N_\mathrm{all}}}$ for element $\mathrm{X}$ other than Hydrogen in `sme` scale;
- $H_\mathrm{H} = 12$ for Hydrogen in `H=12` scale;
- $H_\mathrm{X} = \log{\frac{N_\mathrm{X}}{N_\mathrm{H}}} + 12$ for element $\mathrm{X}$ other than Hydrogen in `H=12` scale.

thus

- $H_\mathrm{X} = S_\mathrm{X} - \log{S_\mathrm{H}}$

### H=12 -> sme

- $S_\mathrm{H} = \left[1 + \sum_\mathrm{X}{10^{H_\mathrm{X}-12}} \right]^{-1}$
- $S_\mathrm{X} = \log{(S_\mathrm{H} \times 10^{H_\mathrm{X}-12})}$

### kurucz -> H=12

We have:
- $K_\mathrm{H} = \frac{N_\mathrm{H}}{N_\mathrm{all}}$ for Hydrogen in `kurucz` scale;
- $K_\mathrm{He} = \frac{N_\mathrm{He}}{N_\mathrm{all}}$ for Helium in `kurucz` scale;
- $K_\mathrm{X} = \log{\frac{N_\mathrm{X}}{N_\mathrm{H}+N_\mathrm{He}}}$ for element $\mathrm{X}$ other than H and He in `kurucz` scale;
- $H_\mathrm{H} = 12$ for Hydrogen in `H=12` scale;
- $H_\mathrm{X} = \log{\frac{N_\mathrm{X}}{N_\mathrm{H}}} + 12$ for element $\mathrm{X}$ other than Hydrogen in `H=12` scale.

thus 

- $H_\mathrm{He} = \log{\frac{K_\mathrm{He}}{K_\mathrm{H}}} + 12$
- $H_\mathrm{X} = K_\mathrm{X} + \log{\left(\frac{K_\mathrm{He}}{K_\mathrm{H}} + 1\right)} + 12$

### H=12 -> kurucz

- $K_\mathrm{H} = \left(1 + \sum_\mathrm{X}10^{H_\mathrm{X}-12} \right)^{-1}$
- $K_\mathrm{He} = \left(10^{12-H_\mathrm{He}} + 1 + 10^{12-H_\mathrm{He}}\times\sum_\mathrm{X}10^{H_\mathrm{X}-12} \right)^{-1}$
- $K_\mathrm{X} = H_\mathrm{X}-12-\log{(1+10^{H_\mathrm{He}-12})}$
