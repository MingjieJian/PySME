Changes from IDL SME
=============

Changes from IDL SME to PySME:

SME structure:
    * Python Class object
    * all properties can be accessed as properties
      (sme.name) or as items (["name"])
    * all properties are case insensitive
    * Can read IDL sav files, but not write, uses numpy npy format instead
    * fitresult values in fitresult object (e.g. covariance matrix)
    * pname, ab_free, and glob_free are now combined in fitparameters
    * idlver still called that, but contains information about Python version
    * cscale is the polynomial coefficients describing
      the continuum (of any degree)
    * preferred access to spectra, through wave, spec, uncs, mask, synth
        * all of which represent the data as Iliffe vectors,
          i.e. a list of arrays of arbitrary length (The underlying data
          structure is the same, and is still accessable, wave is
          renamed to wob) also available are mask_good, mask_bad,
          mask_line, and mask_continuum, which are boolean
          Iliffe vectors of the masks
    * obsolete properties removed (smod_orig, cmod_orig, cmod,
      jint, sint, psig_l, psig_r, rchisq, crms, lrms, vmac_pro,
      cintb, cintr, obs_type, clim)
    * redundant properties are now dynamic (nmu, nseg, md5)
    * property "obs_name" renamed to object
    * property "wave" renamed to wob, but "wave" still exists
      but as an Iliffe vector (as described above)
    * flags, cscale_flag and vrad_flag now use strings instead of integers
        * cscale_flag : "none", "fix", "constant", "linear", "quadratic"
        * vrad_flag : "none", "whole", "each"
    * cscale and vrad are always the expected size
      (i.e. equal to the number of segments, and polynomial degrees)

Abundance:
    * Python Class object
    * set of multiple solar abundances available
    * can output abundances in different formats
      ('sme', 'n/nTot', 'n/nH', and 'H=12')
    * internal format is 'H=12', i.e. log(X/H) + 12
    * metallicity is included as part of the abundance

LineList:
    * Python Class object
    * uses a Pandas Dataframe under the hood
    * combines species, atomic, lulande, extra, depth,
      lineref, term_upp, and term_low
    * also provides relative error measure for each line
      (based on error code in lineref)

NLTE:
    * works similar to IDL
    * NLTE can be activated/deactivated with
      sme.nlte.set_nlte(element_name) / remove_nlte(element_name)
    * NLTE lines are flagged by the C library in sme.nlte.flags

Synthetic Spectrum:
    * uses the same C library as IDL SME
    * sme_func, gives almost the same results as in IDL
    * there are slight differences in the interpolation of the atmospheres
    * continuum and radial velocity are now determined
      at the same time using a least_squares fit (with method "trf")

Solver:
    * Switched solving algorithm to use scipy.optimize.least_squares with method="trf"
        * "trf" = Trust-region Reflective
        * robust loss function "soft_l1"
          resid = 2 * ((1 + z)**0.5 - 1) with z = model - obs / uncs
    * also includes bounds, but only for the outer edges of the atmosphere grid
    * will find different minima than IDL SME, but within a reasonable
      range and often with slightly better chi square
    * for parameters outside the atmosphere grid, the an infinite residual
      is returned, which causes the least squares to avoid that area

Logging:
    * log with standard logging module
    * output both to console and to file

C interface:
    * as mentioned above uses the same library as IDL SME
    * general purpose interface is idl_call_external in cwrapper.py
        * numerics are converted to ctypes elements
        * numpy arrays are converted to appropiate ctypes
          (memory stays the same if possible)
        * strings (and string arrays) are copied into IDL_string
          structures (memory will change)
        * output is written back into original arrays if necessary
    * specific functions are implemented in sme_synth.py
        * data types are converted to the necessary values
        * the size of input arrays is stored for the
          output functions (e.g. number of lines)
        * errors in the C code will raise Exceptions (ValueErrors)

Plotting:
    * two plotting modules
    * one for matplotlib.pyplot (creates a window)
    * one for plot.ly (creates a html file with javascript)
    * both provide functionality to move through the segments
    * can also change the mask from within the plot
      (requires a jupyter notebook session for plot.ly)
