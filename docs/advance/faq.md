# FAQ


## How do I change the default log file?

Call `util.start_logging(filename)`.

```py
from pysme import util
util.start_logging("your_log_file.log")
```

## Get output of `A module that was compiled using NumPy 1.x cannot be run in NumPy 2.2.0 as it may crash.`

This is because the `_smelib` (the 3rd part of [PySME components](fordev.md)) is compiled by NumPy 1.x but the numpy is upgraded to 2.x. 
PySME will autoamtically deal with this error and compile `_smelib` again. 
If you see similar output like:
```
running build_ext
building '_smelib' extension
```
then you should be good to use, and no need to restart python etc.

## I get an error "Derivatives in the starting point are not finite"

Make sure your initial stellar parameters are within the
atmosphere grid defined by the atmosphere file set in sme.atmo.source

## I get an error "lnGAS: DGESVX failed to solved for corrections the partial pressures."

The most possible reason would be the abvundance of the element in error
is too low or nan, thus the EOS code cannot compute its EOS. 
