# FAQ

## Where is the SMElib source code?

The SME native library (`SMElib`) is maintained in a separate repository:

<https://github.com/SpectroscopyMadeEasy/SMElib>

Most PySME users do not need it separately, but it is the right place to look
if you want to inspect or compile the underlying C++/Fortran library.

## What is the default logger level, and how do I change it?

PySME now uses `WARNING` as the default console logger level.

That means:

- `WARNING` and `ERROR` messages are shown by default
- `INFO` and `DEBUG` messages are hidden by default
- Python warnings emitted via `warnings.warn(...)` still appear as usual

If you want more verbose output in the terminal, change both the PySME logger
and the console handler:

```py
import logging
import pysme

pysme.logger.setLevel(logging.INFO)
pysme.console.setLevel(logging.INFO)
```

Use `logging.DEBUG` instead if you want the most detailed output.

## How do I write PySME logs to a file?

Call `util.start_logging(...)`.

```py
from pysme import util
util.start_logging("your_log_file.log")
```

You can also choose the file log level explicitly:

```py
from pysme import util
util.start_logging("your_log_file.log", level="INFO")
```

This configures the PySME logger and adds a file handler. If you also want
`INFO` or `DEBUG` messages to appear in the terminal, change `pysme.console`
as shown above.

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
