#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/arrayobject.h>

// Header of the SME library
#include "sme_synth_faster.h"

// Everything else is considered an error
const char OK_response = '\0';

static char module_docstring[] =
    "This module provides a Python interface to the SME libary";

static char smelib_LibraryVersion_docstring[] = "Return SME library version";
static PyObject *smelib_LibraryVersion(PyObject *self, PyObject *args)
{
    const char *version;
    PyObject *ret;

    version = SMELibraryVersion(0, NULL);
    ret = Py_BuildValue("s", version);

    return ret;
}

static char smelib_GetDataFiles_docstring[] = "Return the required data files";
static PyObject *smelib_GetDataFiles(PyObject *self, PyObject *args)
{
    const char *datafiles;
    PyObject *ret;

    datafiles = GetDataFiles(0, NULL);
    ret = Py_BuildValue("s", datafiles);

    return ret;
}

static char smelib_GetLibraryPath_docstring[] = "Return the current data file directory";
static PyObject *smelib_GetLibraryPath(PyObject *self, PyObject *args)
{
    const char *path;
    PyObject *ret;

    path = GetLibraryPath(0, NULL);
    ret = Py_BuildValue("s", path);

    return ret;
}

static char smelib_SetLibraryPath_docstring[] = "Set the data file directory";
static PyObject *smelib_SetLibraryPath(PyObject *self, PyObject *args)
{
    char *path;
    IDL_STRING idl_path;
    const char *result = NULL;
    void *args_c[1];

    if (!PyArg_ParseTuple(args, "s", &path))
        return NULL;

    // Create IDL String
    idl_path.slen = strlen(path);
    idl_path.s = path;
    idl_path.stype = 0;

    // Pass to SMELIB
    args_c[0] = &idl_path;
    result = SetLibraryPath(1, args_c);

    // Check for errors
    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    // Return the Ok
    Py_RETURN_NONE;
}

static char smelib_InputWaveRange_docstring[] = "Read in Wavelength range";
static PyObject *smelib_InputWaveRange(PyObject *self, PyObject *args)
{
    const int n = 2;
    const char *result = NULL;
    void *args_c[n];
    double wmin, wmax;

    if (!PyArg_ParseTuple(args, "dd", &wmin, &wmax))
        return NULL;

    args_c[0] = &wmin;
    args_c[1] = &wmax;
    result = InputWaveRange(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_SetVWscale_docstring[] = "Set van der Waals scaling factor";
static PyObject *smelib_SetVWscale(PyObject *self, PyObject *args)
{
    const int n = 1;
    const char *result = NULL;
    void *args_c[n];
    double vwscale;

    if (!PyArg_ParseTuple(args, "d", &vwscale))
        return NULL;

    args_c[0] = &vwscale;
    result = SetVWscale(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_SetH2broad_docstring[] = "Set flag for H2 molecule";
static PyObject *smelib_SetH2broad(PyObject *self, PyObject *args)
{
    const int n = 0;
    const char *result = NULL;
    void **args_c = NULL;

    result = SetH2broad(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_ClearH2broad_docstring[] = "Clear flag for H2 molecule";
static PyObject *smelib_ClearH2broad(PyObject *self, PyObject *args)
{
    const int n = 0;
    const char *result = NULL;
    void **args_c = NULL;

    result = ClearH2broad(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_InputLineList_docstring[] = "Read in line list";
static PyObject *smelib_InputLineList(PyObject *self, PyObject *args)
{
    const int n = 3;
    const char *result = NULL;
    void *args_c[n];
    int nlines = 0, nchar = 0;
    PyObject *species_obj = NULL, *linelist_obj = NULL;
    PyArrayObject *species_array = NULL, *linelist_array = NULL;
    IDL_STRING *species = NULL;
    char *species_data = NULL;
    double *linelist = NULL;

    if (!PyArg_ParseTuple(args, "OO", &species_obj, &linelist_obj))
        return NULL;

    // Convert to Numpy arrays
    species_array = (PyArrayObject *)PyArray_FROM_OTF(species_obj, NPY_STRING, NPY_ARRAY_IN_ARRAY);
    linelist_array = (PyArrayObject *)PyArray_FROM_OTF(linelist_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    // Check if that worked
    if (species_array == NULL || linelist_array == NULL)
    {
        goto fail;
    }

    // Check dimensions
    if (PyArray_NDIM(species_array) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected species array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(linelist_array) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "Expected linelist of ndim == 2");
        goto fail;
    }
    // Get sizes
    nlines = PyArray_DIM(species_array, 0);
    nchar = PyArray_ITEMSIZE(species_array);

    if (PyArray_DIM(linelist_array, 0) != 8)
    {
        PyErr_SetString(PyExc_ValueError, "Expected linelist to have 8 values");
        goto fail;
    }
    if (PyArray_DIM(linelist_array, 1) != nlines)
    {
        PyErr_SetString(PyExc_ValueError, "Expected both arrays to have the same length");
        goto fail;
    }

    // Get the pointers to pass to SMELIB
    linelist = (double *)PyArray_DATA(linelist_array);
    species_data = (char *)PyArray_DATA(species_array);
    species = (IDL_STRING *)malloc(nlines * sizeof(IDL_STRING));

    for (int i = 0; i < nlines; i++)
    {
        species[i].slen = nchar;
        species[i].s = &(species_data[nchar * i]);
        species[i].stype = 0;
    }

    args_c[0] = &nlines;
    args_c[1] = species;
    args_c[2] = linelist;
    result = InputLineList(n, args_c);

    free(species);
    Py_DECREF(species_array);
    Py_DECREF(linelist_array);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    return NULL;
}

static char smelib_OutputLineList_docstring[] = "Return line list";
static PyObject *smelib_OutputLineList(PyObject *self, PyObject *args)
{
    const int n = 2;
    const char *result = NULL;
    void *args_c[n];
    int nlines;
    double *linelist = NULL;
    PyArrayObject *linelist_array;

    nlines = GetNLINES();

    npy_intp dims[2] = {nlines, 6};
    linelist_array = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    linelist = (double *)PyArray_DATA(linelist_array);

    args_c[0] = &nlines;
    args_c[1] = linelist;
    result = OutputLineList(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DecRef((PyObject *)linelist_array);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject *)linelist_array;
}

static char smelib_UpdateLineList_docstring[] = "Change line list parameters";
static PyObject *smelib_UpdateLineList(PyObject *self, PyObject *args)
{
    const int n = 4;
    const char *result = NULL;
    void *args_c[n];
    short nlines;
    int nchar;
    double *linelist = NULL;
    IDL_STRING *species = NULL;
    short *index = NULL;
    PyObject *linelist_obj = NULL, *species_obj = NULL, *index_obj = NULL;
    PyArrayObject *linelist_array = NULL, *species_array = NULL, *index_array = NULL;
    char *species_data = NULL;

    if (!PyArg_ParseTuple(args, "OOO", &species_obj, &linelist_obj, &index_obj))
        return NULL;

    // Convert to Numpy arrays
    species_array = (PyArrayObject *)PyArray_FROM_OTF(species_obj, NPY_STRING, NPY_ARRAY_IN_ARRAY);
    linelist_array = (PyArrayObject *)PyArray_FROM_OTF(linelist_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    index_array = (PyArrayObject *)PyArray_FROM_OTF(index_obj, NPY_SHORT, NPY_ARRAY_IN_ARRAY);

    if (species_array == NULL || linelist_array == NULL || index_array == NULL)
    {
        goto fail;
    }

    // Check dimensions
    if (PyArray_NDIM(species_array) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected species array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(linelist_array) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "Expected linelist of ndim == 2");
        goto fail;
    }
    if (PyArray_NDIM(index_array) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected index array of ndim == 1");
        goto fail;
    }

    // Get sizes
    nlines = PyArray_DIM(species_array, 0);
    nchar = PyArray_ITEMSIZE(species_array);

    // Check sizes
    if (PyArray_DIM(linelist_array, 0) != 8)
    {
        PyErr_SetString(PyExc_ValueError, "Expected linelist to have 8 values");
        goto fail;
    }
    if (PyArray_DIM(linelist_array, 1) != nlines)
    {
        PyErr_SetString(PyExc_ValueError, "Expected all arrays to have the same length");
        goto fail;
    }
    if (PyArray_DIM(index_array, 0) != nlines)
    {
        PyErr_SetString(PyExc_ValueError, "Expected all arrays to have the same length");
        goto fail;
    }

    // Get the pointers to pass to SMELIB
    linelist = (double *)PyArray_DATA(linelist_array);
    species_data = (char *)PyArray_DATA(species_array);
    index = (short *)PyArray_DATA(index_array);
    species = (IDL_STRING *)malloc(nlines * sizeof(IDL_STRING));

    for (int i = 0; i < nlines; i++)
    {
        species[i].slen = nchar;
        species[i].s = &species_data[nchar * i];
        species[i].stype = 0;
    }

    args_c[0] = &nlines;
    args_c[1] = species;
    args_c[2] = linelist;
    args_c[3] = index;
    result = UpdateLineList(n, args_c);

    free(species);
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    Py_XDECREF(index_array);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    Py_XDECREF(index_array);
    return NULL;
}

static char smelib_InputModel_docstring[] = "Read in model atmosphere";
static PyObject *smelib_InputModel(PyObject *self, PyObject *args, PyObject *kwds)
{
    int n = 14;
    void *args_c[14];
    const char *result = NULL;
    double teff, grav, wlstd, radius = NAN;
    char *motype;
    IDL_STRING motype_idl;
    short nrhox;
    PyObject *opflag_obj = NULL, *depth_obj = NULL, *temp_obj = NULL;
    PyObject *xne_obj = NULL, *xna_obj = NULL, *rho_obj = NULL, *vt_obj = NULL;
    PyObject *height_obj = NULL;

    PyArrayObject *opflag_arr = NULL, *depth_arr = NULL, *temp_arr = NULL;
    PyArrayObject *xne_arr = NULL, *xna_arr = NULL, *rho_arr = NULL, *vt_arr = NULL;
    PyArrayObject *height_arr = NULL;

    // Need to make this constant because C++ needs it
    // but we cast to non constant and trust Python
    static const char *keywords[] = {"teff", "grav", "wlstd", "motype", "opflag",
                                     "depth", "temp", "xne", "xna", "rho", "vt", "radius", "height", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddsOOOOOOO|dO", const_cast<char **>(keywords),
                                     &teff, &grav, &wlstd, &motype, &opflag_obj, &depth_obj, &temp_obj,
                                     &xne_obj, &xna_obj, &rho_obj, &vt_obj, &radius, &height_obj))
        return NULL;

    // Convert to Numpy arrays
    opflag_arr = (PyArrayObject *)PyArray_FROM_OTF(opflag_obj, NPY_SHORT, NPY_ARRAY_IN_ARRAY);
    depth_arr = (PyArrayObject *)PyArray_FROM_OTF(depth_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    temp_arr = (PyArrayObject *)PyArray_FROM_OTF(temp_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    xna_arr = (PyArrayObject *)PyArray_FROM_OTF(xna_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    xne_arr = (PyArrayObject *)PyArray_FROM_OTF(xne_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    rho_arr = (PyArrayObject *)PyArray_FROM_OTF(rho_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    vt_arr = (PyArrayObject *)PyArray_FROM_OTF(vt_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (opflag_arr == NULL || depth_arr == NULL || temp_arr == NULL || xna_arr == NULL ||
        xne_arr == NULL || rho_arr == NULL || vt_arr == NULL)
    {
        goto fail;
    }

    // Only if given
    if (height_obj != NULL && height_obj != Py_None)
    {
        height_arr = (PyArrayObject *)PyArray_FROM_OTF(height_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (height_arr == NULL)
        {
            goto fail;
        }
    }

    if (PyArray_NDIM(opflag_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected opflag array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(depth_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected depth of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(temp_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected temp array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(xna_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected XNA array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(xne_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected XNE array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(rho_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected rho array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(vt_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected v_turb array of ndim == 1");
        goto fail;
    }
    if (height_arr != NULL)
    {
        if (PyArray_NDIM(height_arr) != 1)
        {
            PyErr_SetString(PyExc_ValueError, "Expected height array of ndim == 1");
            goto fail;
        }
    }

    // Get size
    nrhox = PyArray_DIM(depth_arr, 0);

    // Check size
    if (PyArray_DIM(opflag_arr, 0) != 20)
    {
        PyErr_SetString(PyExc_ValueError, "Expected opflag array of size 20");
        goto fail;
    }
    if (PyArray_DIM(depth_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected depth of the same size");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected temp array of the same size");
        goto fail;
    }
    if (PyArray_DIM(xna_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected XNA array of the same size");
        goto fail;
    }
    if (PyArray_DIM(xne_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected XNE array of the same size");
        goto fail;
    }
    if (PyArray_DIM(rho_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected rho array of the same size");
        goto fail;
    }
    if (PyArray_DIM(vt_arr, 0) != nrhox)
    {
        PyErr_SetString(PyExc_ValueError, "Expected v_turb array of the same size");
        goto fail;
    }
    if (height_arr != NULL)
    {
        if (PyArray_DIM(height_arr, 0) != nrhox)
        {
            PyErr_SetString(PyExc_ValueError, "Expected height array of the same size");
            goto fail;
        }
    }

    // Convert to IDL string
    motype_idl.s = motype;
    motype_idl.slen = strlen(motype);
    motype_idl.stype = 0;

    if (strcmp(motype, "SPH") == 0)
    {
        if (isnan(radius) || height_arr == NULL)
        {
            PyErr_SetString(PyExc_ValueError, "Model type is SPH but no height and/or radius are given");
            goto fail;
        }
        args_c[0] = &nrhox;
        args_c[1] = &teff;
        args_c[2] = &grav;
        args_c[3] = &wlstd;
        args_c[4] = &motype_idl;
        args_c[5] = &radius;
        args_c[6] = (short *)PyArray_DATA(opflag_arr);
        args_c[7] = (double *)PyArray_DATA(depth_arr);
        args_c[8] = (double *)PyArray_DATA(temp_arr);
        args_c[9] = (double *)PyArray_DATA(xne_arr);
        args_c[10] = (double *)PyArray_DATA(xna_arr);
        args_c[11] = (double *)PyArray_DATA(rho_arr);
        args_c[12] = (double *)PyArray_DATA(vt_arr);
        args_c[13] = (double *)PyArray_DATA(height_arr);
        n = 14;
    }
    else
    {
        args_c[0] = &nrhox;
        args_c[1] = &teff;
        args_c[2] = &grav;
        args_c[3] = &wlstd;
        args_c[4] = &motype_idl;
        args_c[5] = (short *)PyArray_DATA(opflag_arr);
        args_c[6] = (double *)PyArray_DATA(depth_arr);
        args_c[7] = (double *)PyArray_DATA(temp_arr);
        args_c[8] = (double *)PyArray_DATA(xne_arr);
        args_c[9] = (double *)PyArray_DATA(xna_arr);
        args_c[10] = (double *)PyArray_DATA(rho_arr);
        args_c[11] = (double *)PyArray_DATA(vt_arr);
        n = 12;
    }

    result = InputModel(n, args_c);

    Py_XDECREF(opflag_arr);
    Py_XDECREF(depth_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(xna_arr);
    Py_XDECREF(xne_arr);
    Py_XDECREF(rho_arr);
    Py_XDECREF(vt_arr);
    Py_XDECREF(height_arr);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(opflag_arr);
    Py_XDECREF(depth_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(xna_arr);
    Py_XDECREF(xne_arr);
    Py_XDECREF(rho_arr);
    Py_XDECREF(vt_arr);
    Py_XDECREF(height_arr);
    return NULL;
}

static char smelib_InputDepartureCoefficients_docstring[] = "Input Departure coefficients";
static PyObject *smelib_InputDepartureCoefficients(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;
    int linenumber;
    short nrhox;
    PyObject *bmatrix_obj = NULL;
    PyArrayObject *bmatrix_arr = NULL;

    if (!PyArg_ParseTuple(args, "Oi", &bmatrix_obj, &linenumber))
        return NULL;

    bmatrix_arr = (PyArrayObject *)PyArray_FROM_OTF(bmatrix_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (bmatrix_arr == NULL)
    {
        return NULL;
    }

    nrhox = GetNRHOX();

    if (PyArray_NDIM(bmatrix_arr) != 2)
    {
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with ndim == 2");
        return NULL;
    }
    if (PyArray_DIM(bmatrix_arr, 1) != 2)
    {
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with shape (nrhox, 2)");
        return NULL;
    }
    if (PyArray_DIM(bmatrix_arr, 0) != nrhox)
    {
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with shape (nrhox, 2)");
        return NULL;
    }

    args_c[0] = (double *)PyArray_DATA(bmatrix_arr);
    args_c[1] = &linenumber;
    result = InputDepartureCoefficients(n, args_c);

    Py_DECREF(bmatrix_arr);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_GetDepartureCoefficients_docstring[] = "Get NLTE b's for specific line";
static PyObject *smelib_GetDepartureCoefficients(PyObject *self, PyObject *args)
{
    const int n = 3;
    void *args_c[n];
    const char *result = NULL;
    int nrhox, linenumber;

    PyArrayObject *bmatrix_arr;
    double *bmatrix_data;

    if (!PyArg_ParseTuple(args, "i", &linenumber))
        return NULL;

    nrhox = GetNRHOX();
    npy_intp dims[2] = {2, nrhox};
    bmatrix_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    bmatrix_data = (double *)PyArray_DATA(bmatrix_arr);

    args_c[0] = bmatrix_data;
    args_c[1] = &nrhox;
    args_c[2] = &linenumber;
    result = GetDepartureCoefficients(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DecRef((PyObject *)bmatrix_arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject *)bmatrix_arr;
}

static char smelib_ResetDepartureCoefficients_docstring[] = "Reset to LTE";
static PyObject *smelib_ResetDepartureCoefficients(PyObject *self, PyObject *args)
{
    const int n = 0;
    void **args_c = NULL;
    const char *result = NULL;

    result = ResetDepartureCoefficients(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    Py_RETURN_NONE;
}

static char smelib_InputAbund_docstring[] = "Read in abundances";
static PyObject *smelib_InputAbund(PyObject *self, PyObject *args)
{
    const int n = 1;
    void *args_c[n];
    const char *result = NULL;

    PyObject *abund_obj = NULL;
    PyArrayObject *abund_arr = NULL;

    if (!PyArg_ParseTuple(args, "O", &abund_obj))
        return NULL;

    abund_arr = (PyArrayObject *)PyArray_FROM_OTF(abund_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (abund_arr == NULL)
    {
        return NULL;
    }

    if (PyArray_NDIM(abund_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected abundance array with ndim == 1");
        Py_XDECREF(abund_arr);
        return NULL;
    }

    if (PyArray_DIM(abund_arr, 0) != 99)
    {
        PyErr_SetString(PyExc_ValueError, "Expected abundance array with size 99");
        Py_XDECREF(abund_arr);
        return NULL;
    }

    args_c[0] = PyArray_DATA(abund_arr);
    result = InputAbund(n, args_c);

    Py_XDECREF(abund_arr);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    Py_RETURN_NONE;
}

static char smelib_Opacity_docstring[] = "Calculate opacities";
static PyObject *smelib_Opacity(PyObject *self, PyObject *args)
{
    const int n = 0;
    void **args_c = NULL;
    const char *result = NULL;

    result = Opacity(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    Py_RETURN_NONE;
}

static char smelib_GetOpacity_docstring[] = "Returns specific cont. opacity";
static PyObject *smelib_GetOpacity(PyObject *self, PyObject *args, PyObject *kwds)
{
    int n = 5;
    void *args_c[5];
    const char *result = NULL;
    char *choice = NULL, *species = NULL, *key = NULL;
    short number = -100;
    short length;
    IDL_STRING species_idl, key_idl;

    PyArrayObject *arr;

    static const char *keywords[] = {"flag", "species", "key", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss", const_cast<char **>(keywords),
                                     &choice, &species, &key))
    {
        return NULL;
    }

    if (strcmp(choice, "COPSTD") == 0)
        number = -3;
    if (strcmp(choice, "COPRED") == 0)
        number = -2;
    if (strcmp(choice, "COPBLU") == 0)
        number = -1;
    if (strcmp(choice, "AHYD") == 0)
        number = 0;
    if (strcmp(choice, "AH2P") == 0)
        number = 1;
    if (strcmp(choice, "AHMIN") == 0)
        number = 2;
    if (strcmp(choice, "SIGH") == 0)
        number = 3;
    if (strcmp(choice, "AHE1") == 0)
        number = 4;
    if (strcmp(choice, "AHE2") == 0)
        number = 5;
    if (strcmp(choice, "AHEMIN") == 0)
        number = 6;
    if (strcmp(choice, "SIGHE") == 0)
        number = 7;
    if (strcmp(choice, "ACOOL") == 0)
        number = 8;
    if (strcmp(choice, "ALUKE") == 0)
        number = 9;
    if (strcmp(choice, "AHOT") == 0)
        number = 10;
    if (strcmp(choice, "SIGEL") == 0)
        number = 11;
    if (strcmp(choice, "SIGH2") == 0)
        number = 12;
    if (number == -100)
    {
        PyErr_SetString(PyExc_ValueError, "Unrecognized Opacity option");
        return NULL;
    }

    if (species != NULL)
    {
        species_idl.slen = strlen(species);
        species_idl.s = species;
        species_idl.stype = 0;
    }
    if (key != NULL)
    {
        key_idl.slen = strlen(key);
        key_idl.s = key;
        key_idl.stype = 0;
    }

    if (number == 8)
    {
        if (species == NULL || key == NULL)
        {
            PyErr_SetString(PyExc_ValueError, "Both species and key keywords need to be set for flag 'ACOOL'");
            return NULL;
        }
        n = 5;
        args_c[3] = &species_idl;
        args_c[4] = &key_idl;
    }
    else if (number == 9)
    {
        if (species == NULL)
        {
            PyErr_SetString(PyExc_ValueError, "Species needs to be set for flag 'ALUKE'");
            return NULL;
        }
        n = 4;
        args_c[3] = &species_idl;
    }
    else
    {
        n = 3;
    }

    length = GetNRHOX();
    npy_intp dims[] = {length};
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    args_c[0] = &number;
    args_c[1] = &length;
    args_c[2] = PyArray_DATA(arr);

    result = GetOpacity(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    return (PyObject *)arr;
}

static char smelib_GetFraction_docstring[] = "Returns fraction/partition function for a species";
static PyObject *smelib_GetFraction(PyObject *self, PyObject *args, PyObject *kwds)
{
    int n = 4;
    void *args_c[4];
    const char *result = NULL;
    char *species = NULL;
    short mode = 2;
    short length;
    IDL_STRING species_idl;

    PyArrayObject *arr;

    static const char *keywords[] = {"species", "mode", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|h", const_cast<char **>(keywords),
                                     &species, &mode))
        return NULL;

    species_idl.slen = strlen(species);
    species_idl.s = species;
    species_idl.stype = 0;

    length = GetNRHOX();
    npy_intp dims[] = {length};
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    args_c[0] = &species_idl;
    args_c[1] = &mode;
    args_c[2] = &length;
    args_c[3] = PyArray_DATA(arr);

    result = GetFraction(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    return (PyObject *)arr;
}

static char smelib_Ionization_docstring[] = "Perform EOS calculations";
static PyObject *smelib_Ionization(PyObject *self, PyObject *args)
{
    const int n = 1;
    void *args_c[n];
    const char *result = NULL;
    short flag = 0;

    if (!PyArg_ParseTuple(args, "h", &flag))
        return NULL;

    args_c[0] = &flag;
    result = Ionization(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_WarnEx(PyExc_Warning, result, 2);
    }
    Py_RETURN_NONE;
}

static char smelib_GetDensity_docstring[] = "Returns density in g/cm^3";
static PyObject *smelib_GetDensity(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;

    short length = 0;
    PyArrayObject *arr = NULL;

    length = GetNRHOX();
    npy_intp dims[] = {length};
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    args_c[0] = &length;
    args_c[1] = PyArray_DATA(arr);
    result = GetDensity(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    return (PyObject *)arr;
}

static char smelib_GetNatom_docstring[] = "Returns atomic number density";
static PyObject *smelib_GetNatom(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;

    short length = 0;
    PyArrayObject *arr = NULL;

    length = GetNRHOX();
    npy_intp dims[] = {length};
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    args_c[0] = &length;
    args_c[1] = PyArray_DATA(arr);
    result = GetNatom(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    return (PyObject *)arr;
}

static char smelib_GetNelec_docstring[] = "Returns electron number density";
static PyObject *smelib_GetNelec(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;

    short length = 0;
    PyArrayObject *arr = NULL;

    length = GetNRHOX();
    npy_intp dims[] = {length};
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    args_c[0] = &length;
    args_c[1] = PyArray_DATA(arr);
    result = GetNelec(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }
    return (PyObject *)arr;
}

static char smelib_Transf_docstring[] = "Computes spectral synthesis";
static PyObject *smelib_Transf(PyObject *self, PyObject *args, PyObject *kwds)
{
    const int n = 12;
    void *args_c[n];
    const char *result = NULL;

    short nmu, keep_lineop = 1, long_continuum = 1;
    int nwmax = 40000, nw = 0;
    double accrt = 1e-4, accwi = 3e-3;
    npy_intp dims[1];
    npy_intp dims2[2];

    PyObject *mu_obj = NULL, *wave_obj = NULL;
    PyObject *return_tuple = NULL;

    PyArrayObject *mu_arr = NULL, *cint_arr = NULL, *cintr_arr = NULL;
    PyArrayObject *wave_arr = NULL, *sint_arr = NULL;

    static const char *keywords[] = {"mu", "wave", "nwmax", "accrt", "accwi", "keep_lineop", "long_continuum", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|Oiddhh", const_cast<char **>(keywords),
                                     &mu_obj, &wave_obj, &nwmax, &accrt, &accwi, &keep_lineop, &long_continuum))
        return NULL;

    mu_arr = (PyArrayObject *)PyArray_FROM_OTF(mu_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mu_arr == NULL)
        goto fail;

    if (PyArray_NDIM(mu_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected mu array of ndim == 1");
        goto fail;
    }

    if (wave_obj != NULL && wave_obj != Py_None)
    {
        // Reuse wavelength grid
        wave_arr = (PyArrayObject *)PyArray_FROM_OTF(wave_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (wave_arr == NULL)
            goto fail;

        if (PyArray_NDIM(wave_arr) != 1)
        {
            PyErr_SetString(PyExc_ValueError, "Expected wavelength array of ndim == 1");
            goto fail;
        }

        nw = PyArray_DIM(wave_arr, 0);
        nwmax = nw;
    }
    else
    {
        // Create a new wavelength grid
        nw = 0;
        dims[0] = nwmax;
        wave_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    }

    nmu = PyArray_DIM(mu_arr, 0);

    dims[0] = nmu;
    cintr_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    dims2[0] = nwmax;
    dims2[1] = nmu;
    cint_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims2, NPY_DOUBLE);
    sint_arr = (PyArrayObject *)PyArray_SimpleNew(2, dims2, NPY_DOUBLE);

    args_c[0] = &nmu;
    args_c[1] = PyArray_DATA(mu_arr);
    args_c[2] = PyArray_DATA(cint_arr);
    args_c[3] = PyArray_DATA(cintr_arr);
    args_c[4] = &nwmax;
    args_c[5] = &nw;
    args_c[6] = PyArray_DATA(wave_arr);
    args_c[7] = PyArray_DATA(sint_arr);
    args_c[8] = &accrt;
    args_c[9] = &accwi;
    args_c[10] = &keep_lineop;
    args_c[11] = &long_continuum;

    result = Transf(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        goto fail;
    }

    // we don't need this one
    Py_DECREF(cintr_arr);

    // return nw, wave, sint_arr, cint_arr
    return_tuple = Py_BuildValue("iOOO", nw, (PyObject *)wave_arr,
                                 (PyObject *)sint_arr, (PyObject *)cint_arr);
    return return_tuple;

fail:
    Py_XDECREF(mu_arr);
    Py_XDECREF(wave_arr);
    Py_XDECREF(cint_arr);
    Py_XDECREF(cintr_arr);
    Py_XDECREF(sint_arr);
    return NULL;
}

static char smelib_CentralDepth_docstring[] = "Computes line central depths";
static PyObject *smelib_CentralDepth(PyObject *self, PyObject *args, PyObject *kwds)
{
    const int n = 5;
    void *args_c[n];
    const char *result = NULL;
    npy_intp dims[1];

    int nmu, nwsize;
    double accrt = 1e-4;
    PyObject *mu_obj = NULL;
    PyArrayObject *mu_arr = NULL, *table_arr = NULL;

    static const char *keywords[] = {"mu", "accrt", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", const_cast<char **>(keywords),
                                     &mu_obj, &accrt))
        return NULL;

    mu_arr = (PyArrayObject *)PyArray_FROM_OTF(mu_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mu_arr == NULL)
        goto fail;

    if (PyArray_NDIM(mu_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected mu array with ndim == 1");
        goto fail;
    }

    nmu = PyArray_DIM(mu_arr, 0);
    nwsize = GetNLINES();

    dims[0] = nwsize;
    table_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT);

    args_c[0] = &nmu;
    args_c[1] = PyArray_DATA(mu_arr);
    args_c[2] = &nwsize;
    args_c[3] = PyArray_DATA(table_arr);
    args_c[4] = &accrt;
    result = CentralDepth(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        goto fail;
    }
    return (PyObject *)table_arr;

fail:
    Py_XDECREF(mu_arr);
    Py_XDECREF(table_arr);
    return NULL;
}

static char smelib_GetLineOpacity_docstring[] = "Returns specific line opacity";
static PyObject *smelib_GetLineOpacity(PyObject *self, PyObject *args)
{
    const int n = 7;
    void *args_c[n];
    const char *result = NULL;
    npy_intp dims[1];
    short depth;
    double wave;
    PyObject *return_tuple;
    PyArrayObject *lop = NULL, *cop = NULL, *scr = NULL, *tsf = NULL, *csf = NULL;

    if (!PyArg_ParseTuple(args, "d", &wave))
        return NULL;

    depth = GetNRHOX();
    dims[0] = depth;
    lop = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    cop = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    scr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    tsf = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    csf = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    // wave, nmu, lop, cop, scr, tsf, csf,
    args_c[0] = &wave;
    args_c[1] = &depth;
    args_c[2] = PyArray_DATA(lop);
    args_c[3] = PyArray_DATA(cop);
    args_c[4] = PyArray_DATA(scr);
    args_c[5] = PyArray_DATA(tsf);
    args_c[6] = PyArray_DATA(csf);
    result = GetLineOpacity(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(lop);
        Py_DECREF(cop);
        Py_DECREF(scr);
        Py_DECREF(tsf);
        Py_DECREF(csf);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return_tuple = Py_BuildValue("OOOOO", (PyObject *)lop, (PyObject *)cop,
                                 (PyObject *)scr, (PyObject *)tsf, (PyObject *)csf);
    return return_tuple;
}

static char smelib_GetLineRange_docstring[] = "Get validity range for every line";
static PyObject *smelib_GetLineRange(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;
    npy_intp dims[2];
    int n_lines;

    PyArrayObject *arr = NULL;

    n_lines = GetNLINES();
    dims[0] = n_lines;
    dims[1] = 2;
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    args_c[0] = PyArray_DATA(arr);
    args_c[1] = &n_lines;
    result = GetLineRange(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject *)arr;
}

static char smelib_GetNLTEflags_docstring[] = "Get line list NLTE flags";
static PyObject *smelib_GetNLTEflags(PyObject *self, PyObject *args)
{
    const int n = 2;
    void *args_c[n];
    const char *result = NULL;
    npy_intp dims[2];
    int n_lines;

    PyArrayObject *arr = NULL;

    n_lines = GetNLINES();
    dims[0] = n_lines;
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_SHORT);

    args_c[0] = PyArray_DATA(arr);
    args_c[1] = &n_lines;
    result = GetNLTEflags(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject *)arr;
}

static char smelib_ContributionFunctions_docstring[] = "Compute contribution functions.";
static PyObject *smelib_ContributionFunctions(PyObject *self, PyObject *args, PyObject *kwds)
{
    /* parse Python arguments here — same signature SME expects:
       int n_mu, double *mu, int nw_avail, double acc_rt, etc. */
    
    const int n = 8;
    void *args_c[n];
    const char *result = NULL;
    
    /* allocate NumPy arrays for the output */
    short nmu, keep_lineop = 1, long_continuum = 1;
    double accrt = 1e-4, accwi = 3e-3;
    int nw = 0, nrhox = GetNRHOX(), nwmax = 40000;

    // npy_intp dims[1], dims2[2], dims3[3];
    npy_intp dims3[3];

    PyObject *mu_obj = NULL, *wave_obj = NULL;
    PyArrayObject *mu_arr = NULL;
    PyArrayObject *wave_arr = NULL;
    PyArrayObject *table_arr = NULL, *ctable_arr = NULL;
    PyObject *ret_tuple = NULL;

    static const char *keywords[] = {"mu", "wave", "nwmax", "accrt", "accwi", "keep_lineop", "long_continuum", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|Oiddhh", const_cast<char **>(keywords),
                                     &mu_obj, &wave_obj, &nwmax, &accrt, &accwi, &keep_lineop, &long_continuum))
        return NULL;
    
    mu_arr = (PyArrayObject *)PyArray_FROM_OTF(mu_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mu_arr == NULL)
        goto fail;

    nmu = PyArray_DIM(mu_arr, 0);
    
    if (PyArray_NDIM(mu_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Expected mu array of ndim == 1");
        goto fail;
    }

    if (wave_obj != NULL && wave_obj != Py_None)
    {
        // Reuse wavelength grid
        wave_arr = (PyArrayObject *)PyArray_FROM_OTF(wave_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (wave_arr == NULL)
            goto fail;

        if (PyArray_NDIM(wave_arr) != 1)
        {
            PyErr_SetString(PyExc_ValueError, "Expected wavelength array of ndim == 1");
            goto fail;
        }

        nw = PyArray_DIM(wave_arr, 0);
        nwmax = nw;
    }
    else
    {
        goto fail;
        // // Create a new wavelength grid
        // nw = 0;
        // dims[0] = nwmax;
        // wave_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    }

    // dims2[0] = nw;
    // dims2[1] = nrhox;
    dims3[0] = nwmax;
    dims3[1] = nmu;
    dims3[2] = nrhox;
    table_arr = (PyArrayObject *)PyArray_SimpleNew(3, dims3, NPY_DOUBLE);
    ctable_arr = (PyArrayObject *)PyArray_SimpleNew(3, dims3, NPY_DOUBLE);

    args_c[0] = &nmu;
    args_c[1] = PyArray_DATA(mu_arr);
    args_c[2] = &nw;
    args_c[3] = &nw;
    args_c[4] = PyArray_DATA(wave_arr);
    args_c[5] = PyArray_DATA(table_arr);
    args_c[6] = PyArray_DATA(ctable_arr);
    args_c[7] = &long_continuum;
    
    result = Contribution_functions(n, args_c);
    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        goto fail;
    }
    ret_tuple = PyTuple_New(2);
    PyTuple_SET_ITEM(ret_tuple, 0, (PyObject *)table_arr);
    PyTuple_SET_ITEM(ret_tuple, 1, (PyObject *)ctable_arr);
    return ret_tuple;
    /* return the NumPy arrays in a tuple */
fail:
return NULL;

}

static PyMethodDef module_methods[] = {
    {"LibraryVersion", smelib_LibraryVersion, METH_NOARGS, smelib_LibraryVersion_docstring},
    {"GetDataFiles", smelib_GetDataFiles, METH_NOARGS, smelib_GetDataFiles_docstring},
    {"GetLibraryPath", smelib_GetLibraryPath, METH_NOARGS, smelib_GetLibraryPath_docstring},
    {"SetLibraryPath", smelib_SetLibraryPath, METH_VARARGS, smelib_SetLibraryPath_docstring},
    {"InputWaveRange", smelib_InputWaveRange, METH_VARARGS, smelib_InputWaveRange_docstring},
    {"SetVWscale", smelib_SetVWscale, METH_VARARGS, smelib_SetVWscale_docstring},
    {"SetH2broad", smelib_SetH2broad, METH_NOARGS, smelib_SetH2broad_docstring},
    {"ClearH2broad", smelib_ClearH2broad, METH_NOARGS, smelib_ClearH2broad_docstring},
    {"InputLineList", smelib_InputLineList, METH_VARARGS, smelib_InputLineList_docstring},
    {"OutputLineList", smelib_OutputLineList, METH_NOARGS, smelib_OutputLineList_docstring},
    {"UpdateLineList", smelib_UpdateLineList, METH_VARARGS, smelib_UpdateLineList_docstring},
    {"InputModel", (PyCFunction)(void (*)(void))smelib_InputModel, METH_VARARGS | METH_KEYWORDS, smelib_InputModel_docstring},
    {"InputDepartureCoefficients", smelib_InputDepartureCoefficients, METH_VARARGS, smelib_InputDepartureCoefficients_docstring},
    {"GetDepartureCoefficients", smelib_GetDepartureCoefficients, METH_VARARGS, smelib_GetDepartureCoefficients_docstring},
    {"ResetDepartureCoefficients", smelib_ResetDepartureCoefficients, METH_NOARGS, smelib_ResetDepartureCoefficients_docstring},
    {"InputAbund", smelib_InputAbund, METH_VARARGS, smelib_InputAbund_docstring},
    {"Opacity", smelib_Opacity, METH_NOARGS, smelib_Opacity_docstring},
    {"GetOpacity", (PyCFunction)(void (*)(void))smelib_GetOpacity, METH_VARARGS | METH_KEYWORDS, smelib_GetOpacity_docstring},
    {"GetFraction", (PyCFunction)(void (*)(void))smelib_GetFraction, METH_VARARGS | METH_KEYWORDS, smelib_GetFraction_docstring},
    {"Ionization", smelib_Ionization, METH_VARARGS, smelib_Ionization_docstring},
    {"GetDensity", smelib_GetDensity, METH_NOARGS, smelib_GetDensity_docstring},
    {"GetNatom", smelib_GetNatom, METH_NOARGS, smelib_GetNatom_docstring},
    {"GetNelec", smelib_GetNelec, METH_NOARGS, smelib_GetNelec_docstring},
    {"Transf", (PyCFunction)(void (*)(void))smelib_Transf, METH_VARARGS | METH_KEYWORDS, smelib_Transf_docstring},
    {"CentralDepth", (PyCFunction)(void (*)(void))smelib_CentralDepth, METH_VARARGS | METH_KEYWORDS, smelib_CentralDepth_docstring},
    {"GetLineOpacity", smelib_GetLineOpacity, METH_VARARGS, smelib_GetLineOpacity_docstring},
    {"GetLineRange", smelib_GetLineRange, METH_NOARGS, smelib_GetLineRange_docstring},
    {"GetNLTEflags", smelib_GetNLTEflags, METH_NOARGS, smelib_GetNLTEflags_docstring},
    {"ContributionFunctions", (PyCFunction)smelib_ContributionFunctions, METH_VARARGS | METH_KEYWORDS, smelib_ContributionFunctions_docstring},
    {NULL, NULL, 0, NULL}};

PyMODINIT_FUNC PyInit__smelib(void)
{
    static struct PyModuleDef smelibmodule = {
        PyModuleDef_HEAD_INIT,
        "smelib",         /* name of module */
        module_docstring, /* module documentation, may be NULL */
        -1,               /* size of per-interpreter state of the module,
                          or -1 if the module keeps state in global variables. */
        module_methods};

    PyObject *m = PyModule_Create(&smelibmodule);
    if (m == NULL)
        return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return m;
}
