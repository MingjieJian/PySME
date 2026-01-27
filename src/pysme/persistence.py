# -*- coding: utf-8 -*-
import io
import json
import logging
import os
import subprocess
import sys
import tempfile
from zipfile import ZIP_LZMA, ZIP_STORED, ZipFile

import numpy as np
from flex.flex import FlexExtension, FlexFile

from . import __version__

logger = logging.getLogger(__name__)


def to_flex(sme):
    header = {}
    extensions = {}

    for name in sme._names:
        value = sme[name]
        if isinstance(value, IPersist):
            extensions[name] = value._save()
        elif isinstance(value, FlexExtension):
            extensions[name] = value
        elif value is not None:
            header[name] = value

    ff = FlexFile(header, extensions)
    return ff


def from_flex(ff, sme):
    header = ff.header
    extensions = ff.extensions
    for name in sme._names:
        if name in updates.keys():
            name = updates[name]
        if name in header.keys():
            sme[name] = header[name]
        elif name in extensions.keys():
            if sme[name] is not None and isinstance(sme[name], IPersist):
                sme[name] = sme[name]._load(extensions[name])
            else:
                sme[name] = extensions[name]
    return sme


def save(filename, sme, format="flex", _async=False):
    """
    Create a folder structure inside a tarfile
    See flex-format for details

    Parameters
    ----------
    filename : str
        Filename of the final file
    sme : SME_Structure
        sme structure to save
    compressed : bool, optional
        whether to compress the output
    """
    ff = to_flex(sme)

    if format == "flex":
        file_ending = ".sme"
    else:
        file_ending = "." + format
    if not filename.endswith(file_ending):
        filename = filename + file_ending

    if format == "flex":
        if _async:
            ff.write_async(filename)
        else:
            ff.write(filename)
    elif format == "fits":
        ff.to_fits(filename, overwrite=True)
    elif format == "json":
        ff.to_json(filename)
    else:
        raise ValueError(
            "Format {!r} not understood, expected one of ['flex', 'fits', 'json'].".format(
                format
            )
        )


def load(fname, sme):
    """
    Load the SME Structure from disk

    Parameters
    ----------
    fname : str
        file to load
    sme : SME_Structure
        empty sme structure with default values set

    Returns
    -------
    sme :  SME_Structure
        loaded sme structure
    """
    try:
        ff = FlexFile.read(fname)
        try:
            sme = from_flex(ff, sme)
        finally:
            ff.close()
        return sme
    except Exception as ex:
        logger.error(ex)
        try:
            sme = load_v1(fname, sme)
        except:
            raise ex
        return sme


# Update this if the names in sme change
updates = {"idlver": "system_info"}


class IPersist:
    def _save(self):
        raise NotImplementedError

    @classmethod
    def _load(cls, ext):
        raise NotImplementedError

    def _save_v1(self, file, folder=""):
        saves_v1(file, self, folder)

    @classmethod
    def _load_v1(cls, file, names, folder=""):
        logger.setLevel(logging.INFO)
        data = cls()  # TODO Suppress warnings
        data = loads_v1(file, data, names, folder)
        logger.setLevel(logging.NOTSET)
        return data


# Version 1 IO (Deprecated)


def toBaseType(value):
    if value is None:
        return value
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.str):
        return str(value)

    return value


def save_v1(filename, data, folder="", compressed=False):
    """
    Create a folder structure inside a zipfile
    Add .json and .npy and .npz files with the correct names
    And subfolders for more complicated objects
    with the same layout
    Each class should have a save and a load method
    which can be used for this purpose

    Parameters
    ----------
    filename : str
        Filename of the final zipfile
    data : SME_struct
        data to save
    folder : str, optional
        subfolder to save data to
    compressed : bool, optional
        whether to compress the output
    """
    # We use LZMA for compression, since that yields the
    # smallest filesize of the existing compression algorithms
    if not compressed:
        compression = ZIP_STORED
    else:
        compression = ZIP_LZMA

    with ZipFile(filename, "w", compression) as file:
        saves_v1(file, data, folder=folder)


# TODO: this is specific for Collection type objects
# Move this to Collection, and not here
def saves_v1(file, data, folder=""):
    if folder != "" and folder[-1] != "/":
        folder = folder + "/"

    parameters = {}
    arrays = {}
    others = {}
    for key in data._names:
        value = getattr(data, key)
        if np.isscalar(value) or isinstance(value, dict):
            parameters[key] = value
        elif isinstance(value, (list, np.ndarray)):
            if np.size(value) > 20:
                arrays[key] = value
            else:
                parameters[key] = value
        else:
            others[key] = value

    info = json.dumps(parameters, default=toBaseType)
    file.writestr(f"{folder}info.json", info)

    for key, value in arrays.items():
        b = io.BytesIO()
        np.save(b, value)
        file.writestr(f"{folder}{key}.npy", b.getvalue())

    for key, value in others.items():
        if value is not None:
            value._save_v1(file, f"{folder}{key}")


def load_v1(filename, data):
    with ZipFile(filename, "r") as file:
        names = file.namelist()
        return loads_v1(file, data, names)


def loads_v1(file, data, names=None, folder=""):
    if folder != "" and folder[-1] != "/":
        folder = folder + "/"
    if names is None:
        names = file.namelist()

    subdirs = {}
    local = []
    for name in names:
        name_within = name[len(folder) :]
        if "/" not in name_within:
            local.append(name)
        else:
            direc, _ = name_within.split("/", 1)
            if direc not in subdirs.keys():
                subdirs[direc] = []
            subdirs[direc].append(name)

    for name in local:
        if name.endswith(".json"):
            info = file.read(name)
            info = json.loads(info)
            for key, value in info.items():
                key = updates.get(key, key)
                data[key] = value
        elif name.endswith(".npy"):
            b = io.BytesIO(file.read(name))
            key = name[len(folder) : -4]
            key = updates.get(key, key)
            data[key] = np.load(b)
        elif name.endswith(".npz"):
            b = io.BytesIO(file.read(name))
            key = name[len(folder) : -4]
            key = updates.get(key, key)
            value = np.load(b)
            data[key] = [value[f"arr_{i}"] for i in range(len(value))]

    for key, value in subdirs.items():
        data_key = updates.get(key, key)
        data[data_key] = data[data_key]._load_v1(file, value, folder=folder + key)

    return data


# IDL IO


def get_typecode(dtype):
    """Get the IDL typecode for a given dtype"""
    if dtype.name[:5] == "bytes":
        return "1"
    if dtype.name == "int16":
        return "2"
    if dtype.name == "int32":
        return "3"
    if dtype.name == "float32":
        return "4"
    if dtype.name == "float64":
        return "5"
    if dtype.name[:3] == "str":
        return dtype.name[3:]
    raise ValueError("Don't recognise the datatype")


temps_to_clean = []


def save_as_binary(arr):
    global temps_to_clean

    with tempfile.NamedTemporaryFile("w+", suffix=".dat", delete=False) as temp:
        if arr.dtype.name[:3] == "str" or arr.dtype.name == "object":
            arr = arr.astype(bytes)
            shape = (arr.dtype.itemsize, len(arr))
        elif np.issubdtype(arr.dtype, np.floating):
            # SME expects double precision, so we assure that here
            arr = arr.astype("float64")
            shape = arr.shape[::-1]
        else:
            shape = arr.shape[::-1]

        # Most arrays should be in the native endianness anyway
        # But if not we swap it to the native representation
        endian = arr.dtype.str[0]
        if endian == "<":
            endian = "little"
        elif endian == ">":
            endian = "big"
        elif endian == "|":
            endian = sys.byteorder

        if endian != sys.byteorder:
            arr = arr.newbyteorder().byteswap()
            endian = "native"

        arr.tofile(temp)
        value = [temp.name, str(list(shape)), get_typecode(arr.dtype), endian]
    temps_to_clean += [temp]
    return value


def clean_temps():
    global temps_to_clean
    for temp in temps_to_clean:
        try:
            os.remove(temp)
        except:
            pass

    temps_to_clean = []


def write_as_idl(sme):
    """
    Write SME structure into and idl format
    data arrays are stored in seperate temp files, and only the filename is passed to idl
    """

    vrad_flag = {"none": -2, "whole": -1, "each": 0, "fix": -2}[sme.vrad_flag]
    # cscale_flag = {"none": -3, "fix": -3, "constant": 0, "linear": 1, "quadratic": 1, }[
    #     sme.cscale_flag
    # ]
    # if not sme.normalize_by_continuum:
    #     cscale_flag = -2

    abund = sme.abund.get_pattern(type="sme", raw=True)
    abund[np.isnan(abund)] = -99

    fitvars = ["TEFF", "GRAV", "FEH", "VMIC", "VMAC", "VSINI", "GAM6", "VRAD"]
    fitvars = [s.upper() for s in sme.fitparameters if s.upper() in fitvars]
    if "logg" in sme.fitparameters:
        fitvars += ["GRAV"]
    if "monh" in sme.fitparameters:
        fitvars += ["FEH"]

    if sme.mask is None and sme.wave is not None:
        sme.mask = 1

    idl_fields = {
        "version": 5.1,
        "id": sme.id,
        "teff": sme.teff,
        "grav": sme.logg,
        "feh": sme.monh,
        "vmic": float(sme.vmic),
        "vmac": float(sme.vmac),
        "vsini": float(sme.vsini),
        "vrad": sme.vrad.tolist() if vrad_flag == 0 else sme.vrad[0],
        "vrad_flag": vrad_flag,
        "cscale": 1.0,
        "cscale_flag": 0,
        "gam6": sme.gam6,
        "h2broad": int(sme.h2broad),
        "accwi": sme.accwi,
        "accrt": sme.accrt,
        "clim": 0.01,
        "maxiter": 100,
        "chirat": 0.002,
        "nmu": sme.nmu,
        "nseg": sme.nseg,
        "abund": save_as_binary(abund),
        "species": save_as_binary(sme.species),
        "atomic": save_as_binary(sme.atomic),
        "lande": save_as_binary(sme.linelist.lande),
        "lineref": save_as_binary(sme.linelist.reference),
        "short_line_format": {"short": 1, "long": 2}[sme.linelist.lineformat],
        "wran": sme.wran.tolist(),
        "mu": sme.mu.tolist() if sme.nmu > 1 else sme.mu[0],
        "obs_name": "",
        "obs_type": 0,
        "glob_free": fitvars if len(fitvars) != 0 else "",
        "atmo": {
            "method": str(sme.atmo.method),
            "source": str(sme.atmo.source),
            "depth": str(sme.atmo.depth),
            "interp": str(sme.atmo.interp),
            "geom": str(sme.atmo.geom),
        },
    }

    if len(sme.nlte.elements) != 0:
        idl_fields["nlte"] = {}

        flags = np.zeros(99, dtype="int16")
        grids = ["" for _ in range(99)]
        for elem in sme.nlte.elements:
            flags[sme.abund.elem_dict[elem]] = 1
            grids[sme.abund.elem_dict[elem]] = sme.nlte.grids[elem]

        idl_fields["nlte"]["nlte_elem_flags"] = save_as_binary(flags)
        idl_fields["nlte"]["nlte_subgrid_size"] = save_as_binary(
            sme.nlte.subgrid_size.astype("int16")
        )
        idl_fields["nlte"]["nlte_grids"] = grids
        idl_fields["nlte"]["nlte_pro"] = "sme_nlte"

    if sme.iptype is not None:
        idl_fields["iptype"] = sme.iptype
        idl_fields["ipres"] = sme.ipres[0]
        # "ip_x": sme.ip_x,
        # "ip_y": sme.ip_y,
    else:
        idl_fields["iptype"] = "gauss"
        idl_fields["ipres"] = 0

    if sme.wave is not None:
        wind = np.cumsum(sme.wave.shape[1]) - 1
        idl_fields["wave"] = save_as_binary(sme.wave.ravel())
        idl_fields["wind"] = wind.tolist()
    if sme.spec is not None:
        idl_fields["sob"] = save_as_binary(sme.spec.ravel())
    if sme.uncs is not None:
        idl_fields["uob"] = save_as_binary(sme.uncs.ravel())
    if sme.mask is not None:
        idl_fields["mob"] = save_as_binary(sme.mask.ravel().astype("int16"))
    if sme.synth is not None:
        idl_fields["smod"] = save_as_binary(sme.synth.ravel())

    if "depth" in sme.linelist.columns:
        idl_fields["depth"] = save_as_binary(sme.linelist.depth)
    else:
        idl_fields["depth"] = save_as_binary(np.ones(len(sme.linelist)))

    if sme.linelist.lineformat == "long":
        idl_fields.update(
            {
                "line_extra": save_as_binary(sme.linelist.extra),
                "line_lulande": save_as_binary(sme.linelist.lulande),
                "line_term_low": save_as_binary(sme.linelist.term_lower),
                "line_term_upp": save_as_binary(sme.linelist.term_upper),
            }
        )

    sep = ""
    text = ""

    for key, value in idl_fields.items():
        if isinstance(value, dict):
            text += f"{sep}{key!s}:{{{key!s},$\n"
            sep = ""
            for key2, value2 in value.items():
                text += f"{sep}{key2!s}:{value2!r}$\n"
                sep = ","
            sep = ","
            text += "}$\n"
        else:
            text += f"{sep}{key!s}:{value!r}$\n"
            sep = ","
    return text


def save_as_idl(sme, fname):
    """
    Save the SME structure to disk as an idl save file

    This writes a IDL script to a temporary file, which is then run
    with idl as a seperate process. Therefore this reqires a working
    idl installation.

    There are two steps to this. First all the fields from the sme,
    structure need to be transformed into simple idl readable structures.
    All large arrays are stored in seperate binary files, for performance.
    The script then reads those files back into idl.
    """
    with tempfile.NamedTemporaryFile("w+", suffix=".pro") as temp:
        tempname = temp.name
        temp.write("print, 'Hello'\n")
        temp.write("sme = {sme,")
        # TODO: Save data as idl compatible data
        temp.write(write_as_idl(sme))
        temp.write("} \n")
        # This is the code that will be run in idl
        temp.write("print, 'there'\n")
        temp.write(
            """tags = tag_names(sme)
print, tags
new_sme = {}

for i = 0, n_elements(tags)-1 do begin
    arr = sme.(i)
    s = size(arr)
    if (s[0] eq 1) and (s[1] eq 4) then begin
        void = execute('shape = ' + arr[1])
        type = fix(arr[2])
        endian = string(arr[3])
        arr = read_binary(arr[0], data_dims=shape, data_type=type, endian=endian)
        if type eq 1 then begin
            ;string
            arr = string(arr)
        endif
    endif
    if (s[s[0]+1] eq 8) then begin
        ;struct
        tags2 = tag_names(sme.(i))
        new2 = {}
        tmp = sme.(i)

        for j = 0, n_elements(tags2)-1 do begin
            arr2 = tmp.(j)
            s = size(arr2)
            if (s[0] eq 1) and (s[1] eq 4) then begin
                void = execute('shape = ' + arr2[1])
                type = fix(arr2[2])
                endian = string(arr2[3])
                arr2 = read_binary(arr2[0], data_dims=shape, data_type=type, endian=endian)
                if type eq 1 then begin
                    ;string
                    arr2 = string(arr2)
                endif
            endif
            new2 = create_struct(temporary(new2), tags2[j], arr2)
        endfor
        arr = new2
    endif
    new_sme = create_struct(temporary(new_sme), tags[i], arr)
endfor

sme = new_sme\n"""
        )
        temp.write(f'save, sme, filename="{fname}"\n')
        temp.write("end\n")
        temp.flush()

        # with open(os.devnull, 'w') as devnull:
        print("IDL Script: ", tempname)
        subprocess.run(["idl", "-e", ".r %s" % tempname])
        # input("Wait for me...")
        clean_temps()
