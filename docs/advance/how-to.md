# PySME how to 

This page describes some new (and in testing) function of PySME.

## How to get the atmosphere grid

```py
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer

# Set up sme structure
sme = SME_Structure()
# Define stellar parameters
sme.teff, sme.logg, sme.monh, sme.vmic = 5777, 4.485, 0, 1.46
# Select atmosphere source; optional, default is MARCS models ('marcs2014.sav')
# Available options are listed in https://pysme-astro.readthedocs.io/en/latest/usage/lfs.html
sme.atmo.source = 'marcs2014.sav'
synthesizer = Synthesizer()
# Interplate the atmopshere
atmo = synthesizer.get_atmosphere(sme).atmo
```

The atmosphere grid is stored in `atmo`:

```py
Atmosphere(teff=5777.0,
logg=4.485,
abund= [M/H]=0.000 applied to abundance pattern. Values below are abundances.
  H      He     Li     Be     B      C      N      O      F      Ne     Na   
 12.000 12.114  1.046  1.376  2.696  8.386  7.776  8.656  4.556  7.836  6.166
  Mg     Al     Si     P      S      Cl     Ar     K      Ca     Sc     Ti   
  7.526  6.366  7.506  5.356  7.136  5.496  6.176  5.076  6.306  3.166  4.896
  V      Cr     Mn     Fe     Co     Ni     Cu     Zn     Ga     Ge     As   
  3.996  5.636  5.386  7.446  4.916  6.226  4.206  4.596  2.876  3.576  2.286
  Se     Br     Kr     Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru   
  3.326  2.556  3.246  2.596  2.916  2.206  2.576  1.416  1.916 -7.964  1.836
  Rh     Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe     Cs   
  1.116  1.656  0.936  1.766  1.596  1.996  0.996  2.186  1.506  2.236  1.066
  Ba     La     Ce     Pr     Nd     Pm     Sm     Eu     Gd     Tb     Dy   
  2.166  1.126  1.696  0.576  1.446 -7.964  0.996  0.516  1.106  0.276  1.136
  Ho     Er     Tm     Yb     Lu     Hf     Ta     W      Re     Os     Ir   
  0.506  0.926 -0.004  1.076  0.056  0.876 -0.174  1.106  0.226  1.246  1.376
  Pt     Au     Hg     Tl     Pb     Bi     Po     At     Rn     Fr     Ra   
  1.636  1.006  1.126  0.896  1.996  0.646 -7.964 -7.964 -7.964 -7.964 -7.964
  Ac     Th     Pa     U      Np     Pu     Am     Cm     Bk     Cf     Es   
 -7.964  0.056 -7.964 -0.524 -7.964 -7.964 -7.964 -7.964 -7.964 -7.964 -7.964,
vturb=1.0,
lonh=1.5,
source='marcs2012.sav',
method='grid',
geom='PP',
radius=1.0,
height=array([-3.28991027e+10,  1.08125362e+01,  7.27886834e+16, ...,
       -3.07627183e+21,  1.08125000e+01,  6.88842390e+27]),
opflag=array([1, 1, 1, ..., 0, 0, 0]),
wlstd=5000.0,
depth='RHOX',
interp='TAU',
rhox=array([0.01206928, 0.0155857 , 0.02009963, ..., 6.73889351, 7.1708134 ,
       7.68721481]),
tau=array([1.89675940e-05, 2.95631869e-05, 4.58225743e-05, ...,
       2.44857945e+01, 3.79713691e+01, 5.92142993e+01]),
temp=array([4101.87925384, 4145.25143759, 4189.80397932, ..., 9422.35229943,
       9677.82956675, 9935.68732234]),
rho=array([1.36430747e-09, 1.74484534e-09, 2.22803405e-09, ...,
       3.24239421e-07, 3.33732002e-07, 3.45981912e-07]),
xna=array([6.48271628e+14, 8.29097296e+14, 1.05869036e+15, ...,
       1.54084995e+17, 1.58596875e+17, 1.64417830e+17]),
xne=array([5.05875737e+10, 6.47950015e+10, 8.29144374e+10, ...,
       4.09172378e+15, 5.26581938e+15, 6.73330283e+15]),
citation_info='\n                @ARTICLE{2008A&A...486..951G,\n                    author = {{Gustafsson}, B. and {Edvardsson}, B. and {Eriksson}, K. and\n                    {J{\\o}rgensen}, U.~G. and {Nordlund}, {\\r{A}}. and {Plez}, B.},\n                    title = "{A grid of MARCS model atmospheres for late-type stars. I. Methods and general properties}",\n                    journal = {Astronomy and Astrophysics},\n                    keywords = {stars: atmospheres, Sun: abundances, stars: fundamental parameters, stars: general, stars: late-type, stars: supergiants, Astrophysics},\n                    year = "2008",\n                    month = "Aug",\n                    volume = {486},\n                    number = {3},\n                    pages = {951-970},\n                    doi = {10.1051/0004-6361:200809724},\n                    archivePrefix = {arXiv},\n                    eprint = {0805.0554},\n                    primaryClass = {astro-ph},\n                    adsurl = {https://ui.adsabs.harvard.edu/abs/2008A&A...486..951G},\n                    adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n                }\n            ')
```

## How to update central depth and line range

Two new parameters for each line, central depth and line range, are added to PySME since 0.4.189. 

Central depth is the line depth **before** any kind of boradening.
It should **not** be used as an indicator of how deep a line is in the final synthetic spectra, since actual line depth also depands on the broadning mechanisms.
However, it can be used for screnning out the weak lines in the line list, thus acclerating the synthesis speed.

Line range indicates the wavelength range where a line has contribution, and it can be used to select the revelant lines for the synthetic wavelength, also for acclerating the synthesis speed.

These parameters can be calculated as follow:

```py
from pysme.sme import SME_Structure
from pysme.linelist.vald import ValdFile
from pysme.synthesize import synthesize_spectrum, Synthesizer
from pysme.solve import solve

linelist = ValdFile('SME/examples/sun.lin')
```

After we read in the line list, it looks like:

```py
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

    lande_upper  lande  ...  gamvw  depth  \
0          0.56   0.73  ... -7.810  0.112   
1          1.79   1.76  ...  0.000  0.011   
2          1.79   1.76  ...  0.000  0.010   
3          1.79   1.76  ...  0.000  0.011   
4          1.79   1.76  ...  0.000  0.010   
..          ...    ...  ...    ...    ...   
62         1.37   1.49  ... -7.511  0.021   
63         1.37   1.49  ... -7.511  0.021   
64         1.37   1.49  ... -7.511  0.021   
65         1.37   1.49  ... -7.511  0.021   
66         1.37   1.49  ... -7.511  0.021   

                                            reference  couple_lower  \
0   N D        Kurucz Fe I 2014Fe               1 ...            LS   
1   E 0.05     Wisconsin REE ex(153)Eu+         3 ...            LS   
2   E 0.05     Wisconsin REE ex(151)Eu+         3 ...            LS   
3   E 0.05     Wisconsin REE ex(153)Eu+         3 ...            LS   
4   E 0.05     Wisconsin REE ex(151)Eu+         3 ...            LS   
..                                                ...           ...   
62  NC         Kurucz MnI 2007 Mn               8 ...            LS   
63  NC         Kurucz MnI 2007 Mn               8 ...            LS   
64  NC         Kurucz MnI 2007 Mn               8 ...            LS   
65  NC         Kurucz MnI 2007 Mn               8 ...            LS   
66  NC         Kurucz MnI 2007 Mn               8 ...            LS   

          term_lower couple_upper                       term_upper error  \
0            3d8 c3F           LS       3d6.(3F2).4s.4p.(3P*) v3D*  0.50   
1   4f7.(8S).5d a9D*           JJ  4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00   
2   4f7.(8S).5d a9D*           JJ  4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00   
3   4f7.(8S).5d a9D*           JJ  4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00   
4   4f7.(8S).5d a9D*           JJ  4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00   
..               ...          ...                              ...   ...   
62       3d5.4s2 b4D           LS                 3d6.(5D).4p z4D*  0.25   
63       3d5.4s2 b4D           LS                 3d6.(5D).4p z4D*  0.25   
64       3d5.4s2 b4D           LS                 3d6.(5D).4p z4D*  0.25   
65       3d5.4s2 b4D           LS                 3d6.(5D).4p z4D*  0.25   
66       3d5.4s2 b4D           LS                 3d6.(5D).4p z4D*  0.25   

   atom_number  ionization  
0          1.0         1.0  
1          1.0         2.0  
2          1.0         2.0  
3          1.0         2.0  
4          1.0         2.0  
..         ...         ...  
62         1.0         1.0  
63         1.0         1.0  
64         1.0         1.0  
65         1.0         1.0  
66         1.0         1.0  

[67 rows x 22 columns]
```

Now we calculate the central depth and line range, in Teff, logg and monh of 5000, 4.0, 0:

```py
sme = SME_Structure()
sme.teff, sme.logg, sme.monh = 5000, 4.0, 0
sme.linelist = linelist

synthesizer = Synthesizer()
sme = synthesizer.update_cdr(sme)
```

Then `sme.linelist` becomes:

```py
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

    lande_upper  lande  ...  couple_lower        term_lower  couple_upper  \
0          0.56   0.73  ...            LS           3d8 c3F            LS   
1          1.79   1.76  ...            LS  4f7.(8S).5d a9D*            JJ   
2          1.79   1.76  ...            LS  4f7.(8S).5d a9D*            JJ   
3          1.79   1.76  ...            LS  4f7.(8S).5d a9D*            JJ   
4          1.79   1.76  ...            LS  4f7.(8S).5d a9D*            JJ   
..          ...    ...  ...           ...               ...           ...   
62         1.37   1.49  ...            LS       3d5.4s2 b4D            LS   
63         1.37   1.49  ...            LS       3d5.4s2 b4D            LS   
64         1.37   1.49  ...            LS       3d5.4s2 b4D            LS   
65         1.37   1.49  ...            LS       3d5.4s2 b4D            LS   
66         1.37   1.49  ...            LS       3d5.4s2 b4D            LS   

                         term_upper error atom_number ionization  \
0        3d6.(3F2).4s.4p.(3P*) v3D*  0.50         1.0        1.0   
1   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00         1.0        2.0   
2   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00         1.0        2.0   
3   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00         1.0        2.0   
4   4f7.(8S<7/2>).6p<3/2> (7/2,3/2)  1.00         1.0        2.0   
..                              ...   ...         ...        ...   
62                 3d6.(5D).4p z4D*  0.25         1.0        1.0   
63                 3d6.(5D).4p z4D*  0.25         1.0        1.0   
64                 3d6.(5D).4p z4D*  0.25         1.0        1.0   
65                 3d6.(5D).4p z4D*  0.25         1.0        1.0   
66                 3d6.(5D).4p z4D*  0.25         1.0        1.0   

   central_depth line_range_s  line_range_e  
0       0.796796    6435.9555     6436.8555  
1       0.250884    6437.3106     6437.9106  
2       0.234380    6437.3121     6437.9121  
3       0.018652    6437.3135     6437.9135  
4       0.017111    6437.3194     6437.9194  
..           ...          ...           ...  
62      0.149999    6443.1689     6443.7689  
63      0.031953    6443.1726     6443.7726  
64      0.388658    6443.1821     6443.7821  
65      0.111368    6443.1887     6443.7887  
66      0.013534    6443.1938     6443.7938  

[67 rows x 25 columns]
```

The last three columns are the central depth, start and end wavelength of line range.
The function calculated these parameters for each line in a batch of 2000 lines (`sme.cdr_N_line_chunk`), in parallel (`sme.cdr_parallel`) of 10 (`sme.cdr_n_jobs`) jobs. 
These vairalbes can be modified to change the behaviour of the function.

The variable `sme.linelist.cdr_paras` will also be updated to an array of [Teff, logg, monh, vmic] used for running `update_cdr`.

## How to synthesize a long spectra

The synthesis of long spectra, i.e., a spectra covering a wide wavelength range (>100AA) or a spectra divided into many chunks but still cover a wide wavelength range, would be a bit tricky.

The deafult synthesis mode of PySME takes all the lines in the `linelist` into account to synthesize the spectra, thus the synthesis takes a long time.
To make it even worse, the synthesis is repeated for each segment, thus if we have 10 segment covering from 3000-9000AA, all the lines (including those in ~9000AA) will be used for the first segment (3000-3600AA), and the total time cost is 10 times longer than just using one segment to cover the whole range.

This situation can be improved by using the central depth and line range parameters to select the relevant lines for each segment, by simply adding `linelist_mode='dynamic'` in `synthesize_spectrum`:

```py
from pysme.synthesize import synthesize_spectrum

sme = SME_Structure()
sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini = 5777, 4.0, 0.07, 1.86, 0, 9.61
sme.wave = np.arange(3000, 6000, 0.02).reshape(-1, 100)
sme.linelist = line_list
sme.nlte.set_nlte('Na')
sme.cdr_depth_thres = 0.1
sme = synthesize_spectrum(sme, linelist_mode='dynamic')
```

`linelist_mode='auto'` is still accepted as a deprecated compatibility alias
and maps to `linelist_mode='dynamic'`.

The detailed explanaiton of what `linelist_mode='dynamic'` triggers is as follow:

1. If the `sme.linelist.cdr_paras` is None (means `update_cdr` hasn't run for it), or either the Teff, logg, monh and vmic (optional) in `sme.linelist.cdr_paras` differ from the input parameters by 250K, 0.5, 0.5 or 1 (you can change these parameter in `sme.linelist.cdr_paras_thres` as a dictionary), then run `update_cdr` under current input stellar parameters.
2. Find the beginning and ending wavelength of each segment, wbeg and wend. 
3. For each segment, only input the lines with:
    - `central_depth` > sme.cdr_depth_thres (default 0), and 
    - `line_range_s` <= wend + del_wav + line_margin, and 
    - `line_range_e` >= wbeg - del_wav - line_margin
4. Finish the synthesis.

Here del_wav is defined as sqrt(sme.vmic^2 + sme.vmac^2 + sme.vsini^2)/c*lambda + lambda/sme.ipres, and line_margin is set to 2AA (just for adding some margin).
NLTE calculation is also availalbe in this mode, by setting `sme.nlte.set_nlte`.

Note that the input line list will not change - only the lines being input to each segment synthesis will differ.  

Removing the `linelist_mode` variable or changing it to `all` falls back to the default way to perform synthesis, described in the beginning of this section.

## How to fit parameter for a long spectra

The vairalbe for such process is as same as the one used in the last section:

```py
from pysme.solve import solve
sme = SME_Structure()
sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini = 5777, 4.0, 0.07, 1.86, 0, 9.61
sme.wave = obs_wave
sme.spec = obs_flux
sme.uncs = obs_flux_err
sme.linelist = line_list
sme.cdr_depth_thres = 0.1
sme = solve(sme, linelist_mode='dynamic')
```

It only triggers the `linelist_mode` variable in `synthesize_spectrum`.
