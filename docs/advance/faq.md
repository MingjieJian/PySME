# FAQ

## Where is the SMElib source code?

The SME native library (`SMElib`) is maintained in a separate repository:

<https://github.com/SpectroscopyMadeEasy/SMElib>

Most PySME users do not need it separately, but it is the right place to look
if you want to inspect or compile the underlying C++/Fortran library.

## How do I change the default log file?

Call `util.start_logging(filename)`.

```py
from pysme import util
util.start_logging("your_log_file.log")
```

## Get output of `A module that was compiled using NumPy 1.x cannot be run in NumPy 2.2.0 as it may crash.`

This is an ABI mismatch between your NumPy and the installed `_smelib` binary.

- For `<= v0.6.26`, PySME may attempt a runtime rebuild/download path.
- For `>= v0.6.27`, the recommended fix is to reinstall PySME in the target environment so binaries are rebuilt/reinstalled consistently.

If you see output like:
```
running build_ext
building '_smelib' extension
```
the extension is being rebuilt.

## I get an error "Derivatives in the starting point are not finite"

Make sure your initial stellar parameters are within the
atmosphere grid defined by the atmosphere file set in sme.atmo.source

## I get an error "lnGAS: DGESVX failed to solved for corrections the partial pressures."

The most possible reason would be the abvundance of the element in error
is too low or nan, thus the EOS code cannot compute its EOS. 
