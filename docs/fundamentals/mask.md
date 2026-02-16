# Mask

Masking the observed spectrum is supported in `solve`. 
The code will then allocate the corresponding pixel to different mask, and perform different manupulation to them.

The mask, `sme.mask` is an array with `dtype=int` with the same length with `sme.wave`.
The types of mask are:

|mask value|mask name|Description|
|:--:|:--:|:--:|
|0|bad pixel|Pixels excluded|
|1|line pixel|Pixels included in the parameter fitting|
|2|continuum pixel|Pixels included in continuum fitting|
|4|vrad pixel|Pixels included in radial velocity fitting|

The masks are additive, i.e., you can set mask value to 5 for line and vrad pixel. Only the good pixels will be included in the fit, but the synthetic spectrum will still be calculated.

In default, mask is all 1.