"""GSLIB"""

# MIT License

# Copyright (c) 2022 Guillaume Rongier

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
import re
import importlib
import pathlib
import platform
import fileinput
import subprocess
import warnings
from dataclasses import dataclass
from typing import Union
import itertools
import numpy as np
import pandas as pd
import xarray as xr

from ..utils import (
    to_list,
    isnumbers,
    rename_duplicated_columns,
    get_ellipsoid_aabb,
    create_dataset,
)
from ..analysis.variography import Direction, VarioModel


################################################################################
# Base


def _convert_column(data, column):
    """
    Converts the names or indices of one or more columns of an array, a DataFrame,
    or a DataSet to GSLIB indices.
    """
    is_list = True
    if isinstance(column, (tuple, list, set, np.ndarray)) == False:
        is_list = False
        column = [column]

    indices = []
    for c in column:
        if pd.isnull(c):
            indices.append(0)
        elif isinstance(c, str):
            if data is not None:
                if hasattr(data, "keys") == True:
                    variable_names = list(data.keys())
                elif hasattr(data, "name") == True:
                    variable_names = [data.name]
                else:
                    variable_names = None
                indices.append(variable_names.index(c) + 1)
            else:
                indices.append(0)
        else:
            indices.append(c + 1)

    if len(indices) == 1 and is_list == False:
        return indices[0]
    return indices


def _convert_parameter(value):
    """
    Converts a parameter to a string for GSLIB.
    """
    value_str = value
    if isinstance(value, bool) == True:
        value_str = int(value)

    return str(value_str)


def _write_header(file, program_name, len_line=40):
    """
    Writes a header into GSLIB parameter file.
    """
    header = "Parameters for " + program_name
    file.write(int(len_line / 2 - len(header) / 2) * " " + header + "\n")
    file.write(
        int(len_line / 2 - len(header) / 2) * " " + "********************" + "\n"
    )
    file.write("\n")
    file.write("START OF PARAMETERS:\n")


def _extract_data_from_dataset(ds, x="X", y="Y", z="Z", r="Realization"):
    """
    Extracts the information needed to save a dataset in a GSLIB data file.

    TODO: Clean-up that function so that u, v, w and x, y, z are consistently
        handled.
    """
    shape = (
        ds[x].size if x in ds.coords else ds.shape[-1],
        ds[y].size if y in ds.coords else ds.shape[-2],
        ds[z].size if z in ds.coords else ds.shape[-3],
    )
    if "spacing" in ds.attrs:
        spacing = ds.attrs["spacing"]
    else:
        spacing = (
            (
                ds[x][1].item() - ds[x][0].item()
                if x in ds.coords and len(ds[x]) > 1
                else 1.0
            ),
            (
                ds[y][1].item() - ds[y][0].item()
                if y in ds.coords and len(ds[y]) > 1
                else 1.0
            ),
            (
                ds[z][1].item() - ds[z][0].item()
                if z in ds.coords and len(ds[z]) > 1
                else 1.0
            ),
        )
    if "origin" in ds.attrs:
        origin = ds.attrs["origin"]
    else:
        origin = (
            ds[x][0].item() if x in ds.coords else 0.5,
            ds[y][0].item() if y in ds.coords else 0.5,
            ds[z][0].item() if z in ds.coords else 0.5,
        )
    n_realizations = None
    dims = ("W", "V", "U")
    if r in ds.sizes:
        n_realizations = ds.sizes[r]
        dims = (r, "W", "V", "U")
    if isinstance(ds, xr.Dataset):
        extra_dims = tuple(set(ds.dims) - set(dims))
        ds = ds.drop_dims(extra_dims, errors="ignore")
    ds = ds.to_dataframe(dims).drop([x, y, z], axis=1).reset_index(drop=True)

    return ds, shape, spacing, origin, n_realizations


def _extract_data_from_array(x):
    """
    Extracts the information needed to save a GSLIB data file.
    """
    shape = None
    spacing = None
    origin = None
    n_realizations = None
    if x.ndim >= 3:
        shape = x.shape[-3:]
        spacing = (1.0, 1.0, 1.0)
        origin = (0.5, 0.5, 0.5)
        if x.ndim == 5:
            n_realizations = x.shape[1]
            x = np.swapaxes(x, 2, 4).reshape(x.shape[0], -1).T
        elif x.ndim == 4:
            x = np.swapaxes(x, 1, 3).reshape(x.shape[0], -1).T
        else:
            x = np.swapaxes(x, 0, 2).reshape(-1, 1)
    ds = pd.DataFrame(x)
    ds.columns = ["Variable " + str(i + 1) for i in range(x.shape[-1])]

    return ds, shape, spacing, origin, n_realizations


def _write_data(
    data, file_path, variable_names=None, shape=None, spacing=None, origin=None
):
    """
    Writes data to a GSLIB data file.
    """
    if isinstance(data, (xr.DataArray, xr.Dataset)) == True:
        data, shape, spacing, origin, n_realizations = _extract_data_from_dataset(data)
    elif isinstance(data, np.ndarray) == True:
        data, shape, spacing, origin, n_realizations = _extract_data_from_array(data)
    if isinstance(data, pd.DataFrame) and variable_names is None:
        variable_names = data.columns

    data = data.loc[:, ~data.map(type).eq(str).all()]

    if variable_names is not None and file_path.suffix != ".trn":
        header = file_path.stem + "\n" + str(len(variable_names))
        if shape is not None:
            for s in shape:
                header += " " + str(s)
            for o in origin:
                header += " " + str(o)
            for s in spacing:
                header += " " + str(s)
            # TODO: Find a more robust way to find the number of realizations.
            if n_realizations is not None:
                header += " " + str(n_realizations)
        header += "\n" + "\n".join(variable_names) + "\n"
    else:
        header = ""

    with open(file_path, "w") as file:
        file.write(header)
    data.to_csv(file_path, sep=" ", header=False, index=False, mode="a")


def _write_calibration(calibration, file_path):
    """
    Writes a calibration file.
    """
    with open(file_path, "w") as file:
        file.write("Thresholds for secondary variable\n")
        file.write(str(len(calibration["Thresholds"])) + "\n")
        for thres in calibration["Thresholds"]:
            file.write(str(thres) + "\n")
        file.write("The local prior distribution table:\n")
        for row in calibration["Calibration_table"]:
            file.write(" ".join(str(x) for x in row) + "\n")
        file.write("The calibration parameters B(i):\n")
        for param in calibration["Calibration_parameters"]:
            file.write(str(param) + "\n")


def _read_data(file, parameters, nan=-99999, var_name=None):
    """
    Reads the values from a GSLIB data file.
    """
    # Expand the parameters
    nvar = int(parameters[0])
    # Read the variable names
    _var_name = [file.readline().rstrip("\n").rstrip().lstrip() for i in range(nvar)]
    if isinstance(var_name, str):
        var_name = [var_name]
    if var_name is None or len(var_name) != len(_var_name):
        var_name = _var_name
    # Read the values
    df = np.loadtxt(file, ndmin=2)
    df[df == nan] = np.nan
    # Create the final dataframe
    df = pd.DataFrame(data=df, columns=var_name)
    # TODO: Modify this part when there's a Pandas' function handling duplicated
    # column names
    df = rename_duplicated_columns(df)

    return df.convert_dtypes()


def _read_regular_grid(file, parameters, nan=-99999, var_name=None):
    """
    Reads the values of a regular grid from a GSLIB data file.
    """
    # Expand the parameters
    nvar, nx, ny, nz = [int(s) for s in parameters[:4]]
    xmn, ymn, zmn, xsiz, ysiz, zsiz = [float(s) for s in parameters[4:10]]
    nrez = int(parameters[-1]) if len(parameters) == 11 else None
    # Reorganize the parameters
    shape = (nz, ny, nx, nvar) if nrez is None else (nrez, nz, ny, nx, nvar)
    spacing = (xsiz, ysiz, zsiz)
    origin = (xmn, ymn, zmn)
    # Read the variable names
    _var_name = [file.readline().rstrip("\n").rstrip().lstrip() for i in range(nvar)]
    if isinstance(var_name, str):
        var_name = [var_name]
    if var_name is None or len(var_name) != len(_var_name):
        var_name = _var_name
    # Read the values
    ds = np.loadtxt(file)
    ds[ds == nan] = np.nan
    ds = ds.reshape(shape)
    ds = ds.swapaxes(0, 2) if nrez is None else ds.swapaxes(1, 3)
    # Create the final dataset
    ds = create_dataset(ds, spacing, origin, n_realizations=nrez, var_name=var_name)

    return ds


def _read_variogram(file):
    """
    Reads the values from a GSLIB variogram file.
    """
    df = {
        "Type": [],
        "Tail": [],
        "Head": [],
        "Direction": [],
        "Lag": [],
        "Distance": [],
        "Value": [],
        "Pair_number": [],
        "Tail_mean": [],
        "Head_mean": [],
    }

    file.seek(0)
    line = file.readline()
    while line:
        line = line.rstrip("\n")
        if isnumbers(line.split()[0]) == False:
            description = re.split("tail:|head:|direction", line)
        else:
            line = line.split()
            df["Type"].append(description[0].rstrip())
            df["Tail"].append(description[1].rstrip().replace("  ", " "))
            df["Head"].append(description[2].rstrip().replace("  ", " "))
            df["Direction"].append(int(description[3]))
            df["Lag"].append(int(line[0]))
            df["Distance"].append(float(line[1]))
            df["Value"].append(float(line[2]))
            df["Pair_number"].append(int(line[3]))
            df["Tail_mean"].append(float(line[4]))
            df["Head_mean"].append(float(line[5]))
        line = file.readline()

    return pd.DataFrame(df)


def _read_variogram_model(file):
    """
    Reads the values from a GSLIB variogram file.
    """
    df = {
        "Direction": [],
        "Distance": [],
        "Semivariance": [],
        "Covariance": [],
        "Correlation": [],
    }

    file.seek(0)
    line = file.readline()
    while line:
        line = line.rstrip("\n")
        if isnumbers(line.split()[0]) == False:
            description = re.split(":", line)
        else:
            line = line.split()
            df["Direction"].append(int(description[-1]))
            df["Distance"].append(float(line[1]))
            df["Semivariance"].append(float(line[2]))
            df["Covariance"].append(float(line[4]))
            df["Correlation"].append(float(line[5]))
        line = file.readline()

    return pd.DataFrame(df)


def _read_calibration(file):
    """
    Reads the calibration parameters from the calibration file generated by BICALIB.
    """
    n_thres = int(file.readline())
    dic = {
        "Thresholds": np.empty(n_thres),
        "Calibration_table": [],
        "Calibration_parameters": np.empty(n_thres),
    }
    for i in range(n_thres):
        dic["Thresholds"][i] = float(file.readline())
    line = file.readline()
    line = file.readline()
    while "calibration parameters" not in line:
        line = [float(x) for x in line.rstrip("\n").split()]
        dic["Calibration_table"].append(line)
        line = file.readline()
    dic["Calibration_table"] = np.array(dic["Calibration_table"])
    for i in range(n_thres):
        dic["Calibration_parameters"][i] = float(file.readline())

    return dic


class DataManager:
    """
    Manager for the data files used as inputs to GSLIB programs or created as
    outputs by GSLIB programs.

    Parameters
    ----------
    output_file_name : str, default='gslib'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(self, output_file_name="gslib", output_dir_path="."):

        self.output_dir_path = pathlib.Path(output_dir_path)  # .resolve()
        self.output_dir_path.mkdir(parents=True, exist_ok=True)
        self.output_file_path = self.output_dir_path / output_file_name

    def write_data(
        self,
        data,
        file_name="data.dat",
        variable_names=None,
        shape=None,
        spacing=None,
        origin=None,
    ):
        """
        Writes some data to a GSLIB data file.

        Parameters
        ----------
        data : array-like of shape (-1, nvar)
            Values to write into the file.
        file_name : str, default='data.dat'
            Name of the data file. If no extension is provided, the extension
            `.dat` is used.
        variable_names : array-like of shape (nvar,)
            Name of the variables in `data`. If `data` is a dataframe or a
            dataset, this parameter is retrieved automatically.
        shape : array-like of shape (3,)
            Shape of the grid if the data is on a regular grid. If `data` is
            a dataset, this parameter is retrieved automatically.
        spacing : array-like of shape (3,)
            Spacing of the cells if the data is on a regular grid. If `data` is
            a dataset, this parameter is retrieved automatically.
        origin : array-like of shape (3,)
            Spacing of the cells if the data is on a regular grid. If `data` is
            a dataset, this parameter is retrieved automatically.
        """
        file_path = self.output_dir_path / file_name
        if file_path.suffix == "":
            file_path = file_path.with_suffix(".dat")

        if isinstance(data, dict) and "Calibration_parameters" in data:
            _write_calibration(data, file_path)
        else:
            _write_data(data, file_path, variable_names, shape, spacing, origin)

    def read_data(self, file_path=None, nan=-99999, var_name=None):
        """
        Reads some data from a GSLIB data file.

        Parameters
        ----------
        file_path : str, default=None
            Path to the data file. Uses `output_file_path` by default, which is
            only defined in the daughter classes.
        nan : float, default=-99999
            Value for NaN in the file.
        var_name : string or array-like of shape (n_vars,), default=None
            Name(s) of the variable(s) to use instead of those in the file, if
            any.

        Returns
        -------
        data : array-like
            The data from the file.
        """
        if file_path is None:
            file_path = self.output_file_path.with_suffix(".out")
        elif isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        data = None
        if file_path.exists():
            with open(file_path) as file:
                # Skip the title
                line = file.readline()
                if isnumbers(line) == True:
                    file.seek(0)
                    data = np.loadtxt(file_path)
                elif len(line.split()) > 1 and "tail" in line and "head" in line:
                    data = _read_variogram(file)
                elif "Model Variogram" in line:
                    data = _read_variogram_model(file)
                elif len(line.split()) > 1 and line.split()[0] == "Thresholds":
                    data = _read_calibration(file)
                else:
                    # Read the header
                    line = file.readline()
                    parameters = line.rstrip("\n").split()

                    if len(parameters) > 1:
                        data = _read_regular_grid(file, parameters, nan, var_name)
                    else:
                        data = _read_data(file, parameters, nan, var_name)
        else:
            warnings.warn("File " + str(file_path) + " does not exists.")

        return data

    def read_sgems_grid(self, file_path, spacing=None, origin=None, nan=-99999.0):
        """
        Reads some gridded data from a SGEMS data file.

        Parameters
        ----------
        file_path : str
            Path to the data file. Uses `output_file_path` by default, which is
            only defined in the daughter classes.
        spacing : array-like of shape (3,)
            Spacing of the cells if the data is on a regular grid.
        origin : array-like of shape (3,)
            Spacing of the cells if the data is on a regular grid.
        nan : float, default=-99999
            Value for NaN in the file.

        Returns
        -------
        data : array-like
            The data from the file.
        """
        if spacing is not None and origin is None:
            origin = tuple(s / 2.0 for s in spacing)
        elif spacing is None and origin is not None:
            spacing = tuple(o * 2.0 for o in origin)
        else:
            spacing = (1.0, 1.0, 1.0)
            origin = (0.5, 0.5, 0.5)

        ds = None
        if file_path.exists():
            with open(file_path) as file:
                parameters = file.readline().rstrip("\n").split()
                # Expands the parameters
                nx, ny, nz = [int(s) for s in parameters]
                nvar = int(file.readline().rstrip("\n").strip())
                # Read the variable names
                var_names = [
                    file.readline().rstrip("\n").rstrip().lstrip() for i in range(nvar)
                ]
                # Read the values
                ds = np.loadtxt(file)
                ds[ds == nan] = np.nan
                ds = ds.reshape((nz, ny, nx, nvar)).T
                ds = xr.Dataset(
                    data_vars={
                        var_name: (["U", "V", "W"], ds[i])
                        for i, var_name in enumerate(var_names)
                    },
                    coords={
                        "X": (
                            ["U"],
                            np.linspace(
                                origin[0], origin[0] + spacing[0] * (nx - 1), nx
                            ),
                        ),
                        "Y": (
                            ["V"],
                            np.linspace(
                                origin[1], origin[1] + spacing[1] * (ny - 1), ny
                            ),
                        ),
                        "Z": (
                            ["W"],
                            np.linspace(
                                origin[2], origin[2] + spacing[2] * (nz - 1), nz
                            ),
                        ),
                    },
                )
                ds.attrs["shape"] = (nx, ny, nz)
                ds.attrs["spacing"] = spacing
                ds.attrs["origin"] = origin
        else:
            warnings.warn("File " + str(file_path) + " does not exists.")

        return ds

    def convert_data_to_grid(self, data, shape, spacing=None, origin=None):
        """
        Converts a DataFrame to a DataSet containing a regular grid.

        Parameters
        ----------
        data : DataFrame
            Data to convert to a gridded format.
        shape : array-like (x, y, z)
            3D shape of the grid.
        spacing : array-like (x, y, z)
            Cell size of the grid.
        origin : array-like (x, y, z)
            Origin of the grid.

        Returns
        -------
        data : array-like
            The gridded data.
        """
        if spacing is not None and origin is None:
            origin = tuple(s / 2.0 for s in spacing)
        elif spacing is None and origin is not None:
            spacing = tuple(o * 2.0 for o in origin)
        else:
            spacing = (1.0, 1.0, 1.0)
            origin = (0.5, 0.5, 0.5)

        nvar = len(data.columns)
        ds = np.reshape(data.to_numpy(dtype=float), shape[::-1] + (nvar,)).T
        ds = xr.Dataset(
            data_vars={
                var_name: (["U", "V", "W"], ds[i])
                for i, var_name in enumerate(data.columns)
            },
            coords={
                "X": (
                    ["U"],
                    np.linspace(
                        origin[0], origin[0] + spacing[0] * (shape[0] - 1), shape[0]
                    ),
                ),
                "Y": (
                    ["V"],
                    np.linspace(
                        origin[1], origin[1] + spacing[1] * (shape[1] - 1), shape[1]
                    ),
                ),
                "Z": (
                    ["W"],
                    np.linspace(
                        origin[2], origin[2] + spacing[2] * (shape[2] - 1), shape[2]
                    ),
                ),
            },
        )
        ds.attrs["shape"] = tuple(shape)
        ds.attrs["spacing"] = tuple(spacing)
        ds.attrs["origin"] = tuple(origin)

        return ds


class GSLIB(DataManager):
    """
    Base wrapper class to run GSLIB programs in Python.

    Parameters
    ----------
    output_file_name : str, default='gslib'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self, output_file_name="gslib", output_dir_path=".", executable_path=None
    ):

        DataManager.__init__(self, output_file_name, output_dir_path)
        if executable_path is None:
            self.executable_path = pathlib.Path(
                importlib.util.find_spec("geomodpy").origin
            ).parents[0]
            self.executable_path /= "wrapper/bin"
        else:
            self.executable_path = executable_path
        self.parameters = []
        self.input_file_paths = []
        self.output_file_paths = [self.output_file_path.with_suffix(".out")]

        self.variable_name = None

    def _set_executable(self, executable_path, executable_name):
        """
        Sets or completes the executable path.
        """
        if self.executable_path is None:
            self.executable_path = pathlib.Path(executable_path).expanduser()

        if (self.executable_path / executable_name).exists() == False:
            if platform.system() == "Windows":
                self.executable_path /= "windows/" + executable_name
            elif platform.system() == "Linux":
                self.executable_path /= "linux/" + executable_name
            elif platform.system() == "Darwin":
                if platform.processor() == "i386":
                    self.executable_path /= "macos/intel/" + executable_name
                else:
                    self.executable_path /= "macos/silicon/" + executable_name
        else:
            self.executable_path /= executable_name

    def _write_parameter(self, file, value, description, file_name=None, len_line=40):
        """
        Writes a parameter value into GSLIB parameter file.
        """
        if isinstance(value, (tuple, list)) == True or (
            isinstance(value, np.ndarray) == True
            and value.ndim == 1
            and file_name is None
        ):
            value_str = "  ".join([_convert_parameter(i) for i in value])
        elif (
            isinstance(
                value, (np.ndarray, pd.DataFrame, xr.DataArray, xr.Dataset, dict)
            )
            == True
        ):
            self.write_data(value, file_name=file_name)
            value_str = self.output_dir_path / file_name
            self.input_file_paths.append(value_str)
            value_str = str(value_str)
        else:
            value_str = _convert_parameter(value)
        file.write(
            value_str + max(len_line - len(value_str), 1) * " " + description + "\n"
        )

    def _write_input_parameters(self):
        """
        Writes the input parameters into a parameter file.
        """
        with open(self.output_file_path.with_suffix(".par"), "w") as file:
            _write_header(file, type(self).__name__)
            for parameter in self.parameters:
                self._write_parameter(file, *parameter)

    def process_output_files(self):
        """
        Processes the output files before reading.
        """
        pass

    def clean_files(self):
        """
        Deletes the files generated by GSLIB programs.
        """
        for path in self.input_file_paths:
            path.unlink(missing_ok=True)
        for ext in [
            ".ap",
            ".dbg",
            ".cal",
            ".geo",
            ".out",
            ".par",
            ".rep",
            ".sum",
            ".trn",
            ".vp",
            ".wd",
        ]:
            self.output_file_path.with_suffix(ext).unlink(missing_ok=True)

    def run(self, parameter_file_path=None, print_log=True, clean_files=True):
        """
        Runs a GSLIB program.

        Parameters
        ----------
        parameter_file_path : str, default=None
            Path to the parameter file, which is defined in the daughter classes.
        print_log : bool, default=True
            If true, prints the log of the GSLIB program while it's running,
            otherwise doesn't print anything.
        clean_files : bool, default=True
            If true, delete the files created for and by the run, otherwise keep them.

        Returns
        -------
        output : array-like
            The output generated by the program.
        """
        if parameter_file_path is None:
            self._write_input_parameters()
            parameter_file_path = self.output_file_path.with_suffix(".par")

        if print_log == True:
            stdout = subprocess.PIPE
            stderr = subprocess.STDOUT
        else:
            stdout = subprocess.DEVNULL
            stderr = subprocess.DEVNULL
        with subprocess.Popen(
            [self.executable_path, parameter_file_path],
            stdout=stdout,
            stderr=stderr,
            cwd=self.output_dir_path,
            encoding="latin-1",
        ) as process:
            if print_log == True:
                while process.poll() is None:
                    out = process.stdout.readline()
                    if out:
                        print(out.strip())
                print(process.stdout.read())
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)

        self.process_output_files()
        if len(self.output_file_paths) > 1:
            output = [
                self.read_data(path, var_name=self.variable_name)
                for path in self.output_file_paths
            ]
        else:
            output = self.read_data(
                self.output_file_paths[0], var_name=self.variable_name
            )

        if clean_files == True:
            self.clean_files()

        return output


################################################################################
# Preprocessing


class DECLUS(GSLIB):
    """
    Cell declustering.

    Parameters
    ----------
    data : DataFrame
        Data to decluster.
    coord_columns : array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_column : int or str, default=2
        Index or name of the variable to decluster in `data`.
    y_anisotropy_factor : float, default=1.
        Anisotropy factor to consider rectangular cells. The cell size in the x
        direction is multiplied by this factor to get the cell size in the y
        direction; e.g., if a cell size of 10 is being considered and `y_anisotropy_factor`
        is 2, then the cell size in the y direction is 20.
    z_anisotropy_factor : float, default=1.
        Anisotropy factor to consider rectangular cells. The cell size in the x
        direction is multiplied by this factor to get the cell size in the z
        direction; e.g., if a cell size of 10 is being considered and `y_anisotropy_factor`
        is 3, then the cell size in the y direction is 30.
    use_min_mean : bool, default=False
        If True, a minimum mean value is looked for, otherwise a maximum mean
        value is looked for.
    n_cell_sizes : int, default=24
        Number of cell sizes to consider.
    cell_size_range : array-like (min, max), default=(1., 25.)
        Minimum and maximum cell size. These sizes apply directly to the x direction
        and the anisotropy factors adjust the sizes in the y and z directions.
    n_origin_offets : int, default=5
        Number of origin offsets. Each of the ncell cell sizes are considered with
        `n_origin_offets` different original starting points. This avoids erratic
        results caused by extreme values falling into specific cells. A good number
        is 4 in 2D and 8 in 3D.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='backtr'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        data_column=2,
        y_anisotropy_factor=1.0,
        z_anisotropy_factor=1.0,
        use_min_mean=False,
        n_cell_sizes=24,
        cell_size_range=(1.0, 25.0),
        n_origin_offets=5,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="declus",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "declus")

        self.parameters = [
            (data, "\\file with data", output_file_name + "_cluster.dat"),
            (
                _convert_column(
                    data, to_list(coord_columns, 3, np.nan) + [data_column]
                ),
                "\\ columns for X, Y, Z, and variable",
            ),
            (trimming_limits, "\\  triming limits"),
            (output_file_name + ".sum", "\\file for summary output"),
            (output_file_name + ".out", "\\file for output with data & weights"),
            (
                (y_anisotropy_factor, z_anisotropy_factor),
                "\\Y and Z cell anisotropy (Ysize=size*Yanis)",
            ),
            (use_min_mean, "\\0=look for minimum declustered mean (1=max)"),
            (
                [n_cell_sizes] + list(cell_size_range),
                "\\number of cell sizes, min size, max size",
            ),
            (n_origin_offets, "\\number of origin offsets"),
        ]

        self.output_file_paths += [self.output_file_path.with_suffix(".sum")]


class NSCORE(GSLIB):
    """
    Normal score transformation.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    data_columns : int, str, or array-like (var, weight), default=(2, np.nan)
        Indices or names of the variable and weight in `data`. If the weight is
        np.nan, pd.NA, or omitted, then equal weighting is applied.
    smoothed_distribution : array-like of shape (1,) or (2,), default=None
        A smoothed distribution for the transformation. This can be useful when
        the number of samples is limited and the original distribution is deemed
        too discrete so not representative.
    smoothed_distribution_columns : int, str, or array-like (var, weight), default=(1, np.nan)
        Indices or names of the variable and weight in `smoothed_distribution`.
        If the weight is np.nan, pd.NA, or omitted, then equal weighting is applied.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='nscore'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        data_columns=(2, np.nan),
        smoothed_distribution=None,
        smoothed_distribution_columns=(0, np.nan),
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="nscore",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "nscore")

        self.parameters = [
            (data, "\\file with data", output_file_name + "_distrib.dat"),
            (
                _convert_column(data, to_list(data_columns, 2, np.nan)),
                "\\  columns for variable and weight",
            ),
            (trimming_limits, "\\  triming limits"),
            (
                False if smoothed_distribution is None else True,
                "\\l=transform according to specified ref. dist.",
            ),
            (
                smoothed_distribution,
                "\\  file with reference dist.",
                output_file_name + "_smooth_distrib.dat",
            ),
            (
                _convert_column(smoothed_distribution, smoothed_distribution_columns),
                "\\  columns for variable and weight",
            ),
            (output_file_name + ".out", "\\file for output"),
            (output_file_name + ".trn", "\\file for output transformation table"),
        ]

        self.output_file_paths += [self.output_file_path.with_suffix(".trn")]


class BACKTR(GSLIB):
    """
    Normal score back transformation.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    transformation_table : array-like of shape (2,)
        Transformation table to back-transform the data, i.e., pairs of z and y
        values with z being the original data values or class bound values and y
        being the corresponding normal scores values.
    data_column : int, default=3
        Index or name of the variable in `data`.
    extrema : array-like (min, max), default=None
        Minimum and maximum values to be used for extrapolation in the tails.
    lower_tail : array-like (opt, param), default=None
        The first value specifies the back-transformation implementation in the
        lower tail of the distribution:
        1 = linear interpolation to the lower extrema;
        2 = power model interpolation to the lower extrema.
        The second value specifies the parameter of option 2.
    upper_tail : array-like (opt, param), default=None
        The first value specifies the back-transformation implementation in the
        upper tail of the distribution:
        1 = linear interpolation to the upper extrema;
        2 = power model interpolation to the upper extrema.
        4 = hyperbolic model extrapolation to the upper extrema.
        The second value specifies the parameter of option 2 or 4.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='backtr'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        transformation_table,
        data_column=3,
        extrema=None,
        lower_tail=None,
        upper_tail=None,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="backtr",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "backtr")

        vmin, vmax = transformation_table[0, 0], transformation_table[-1, 0]
        if extrema is None:
            extrema = (vmin, vmax)
        if lower_tail is None:
            lower_tail = (1, vmin)
        if upper_tail is None:
            upper_tail = (1, vmax)

        self.parameters = [
            (data, "\\file with data", output_file_name + "_nscore.dat"),
            (_convert_column(data, data_column), "\\  column with Gaussian variable"),
            (trimming_limits, "\\  triming limits"),
            (output_file_name + ".out", "\\file for output"),
            (
                transformation_table,
                "\\file with input transformation table",
                output_file_name + "_nscore.trn",
            ),
            (extrema, "\\minimum and maximum data value"),
            (lower_tail, "\\lower tail option and parameter"),
            (upper_tail, "\\upper tail option and parameter"),
        ]


################################################################################
# Variogram modeling


@dataclass(frozen=True)
class ExpVario:
    """
    Parameters for an experimental variogram.

    Parameters
    ----------
    tail_column : int or str, default=4
        Index or name of the tail variable from which to compute a variogram.
    head_column : int or str, default=4
        Index or name of the head variable from which to compute a variogram.
    vario_type : int or str, default=4
        Variogram type:
            1 = semivariogram;
            2 = cross semivariogram;
            3 = covariance;
            4 = correlogram;
            5 = general relative semivariogram;
            6 = pairwise relative semivariogram;
            7 = semivariogram of logarithms;
            8 = semimadogram;
            9 = indicator semivariogram (continuous variable);
            10 = indicator semivariogram (categorical variable).
    cutoff : float, default=None
        Only if type is 9 or 10, i.e., an indicator variogram.
    """

    tail_column: Union[int, str] = 4
    head_column: Union[int, str] = 4
    vario_type: int = 1
    cutoff: float = None


@dataclass(frozen=True)
class Offset:
    """
    Parameters for the direction along which to compute a variogram on a regular
    grid.

    Parameters
    ----------
    x_shift : int, default=1
        Number of grid nodes that must be shifted to move from a node on the grid
        to the next nearest node on the grid along the x axis.
    y_shift : int, default=1
        Number of grid nodes that must be shifted to move from a node on the grid
        to the next nearest node on the grid along the y axis.
    z_shift : int, default=0
        Number of grid nodes that must be shifted to move from a node on the grid
        to the next nearest node on the grid along the z axis.
    """

    x_shift: int = 1
    y_shift: int = 1
    z_shift: int = 0


class GAM(GSLIB):
    """
    Variogram calculation with regular data.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    exp_vario : ExpVario or array-like of shape (n_vario,), default=ExpVario()
        Experimental variogram(s) to compute.
    i_realization : int, default=0
        Grid or realization index
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    offset : Offset or array-like of shape (n_off,), default=Offset()
        Offset(s) along which to compute the variogram.
    n_lags : int, default=10
        Number of lags
    standardize_sill : bool, default=True
        If true, standardizes the sill.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='gam'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        exp_vario=ExpVario(),
        i_realization=0,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        offset=Offset(),
        n_lags=10,
        standardize_sill=True,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="gam",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "gam")

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        if isinstance(exp_vario, (tuple, list)) == False:
            exp_vario = [exp_vario]
        if isinstance(offset, (tuple, list)) == False:
            offset = [offset]
        data_columns = set()
        for vario in exp_vario:
            data_columns.add(vario.tail_column)
            data_columns.add(vario.head_column)
        data_columns = list(data_columns)

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                [len(data_columns)] + _convert_column(data, data_columns),
                "\\  number of variables, column numbers",
            ),
            (trimming_limits, "\\  triming limits"),
            (output_file_name + ".out", "\\file for variogram output"),
            (i_realization + 1, "\\grid or realization number"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            ((len(offset), n_lags), "\\number of directions, number of lags"),
        ]
        for _offset in offset:
            _offset = [_offset.x_shift, _offset.y_shift, _offset.z_shift]
            self.parameters.append((_offset, "\\ixd,iyd,izd"))
        self.parameters.append((standardize_sill, "\\standardize sill? (O=no, l=yes)"))
        self.parameters.append((len(exp_vario), "\\number of variograms"))
        for vario in exp_vario:
            _vario = [
                data_columns.index(vario.tail_column) + 1,
                data_columns.index(vario.head_column) + 1,
                vario.vario_type,
            ]
            if vario.cutoff is not None:
                _vario.append(vario.cutoff)
            self.parameters.append(
                (_vario, "\\tail variable, head variable, variogram type")
            )


class _GAMV(GSLIB):
    """
    Variogram calculation with irregular data. This implementation follows GSLIB
    inputs exactly: the number of lags and the lag distance are the same for each
    direction.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. In 1D or 2D,
        use np.nan or pd.NA to signal the coordinates to ignore. The last and
        second-to-last coordinates can also be omitted instead of set to NaN.
    exp_vario : ExpVario or array-like of shape (n_vario,), default=ExpVario()
        Experimental variogram(s) to compute.
    n_lags : int, default=10
        Number of lags.
    lag_separation_distance : float, default=5.
        Separation distance between each lag.
    lag_tolerance : float, default=3.
        Tolerance when determining which data ends in a lag.
    direction : Direction or array-like of shape (n_dir,), default=Direction()
        Direction(s) along which to compute the variogram.
    standardize_sill : bool, default=True
        If true, standardizes the sill.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='gamv'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        exp_vario=ExpVario(),
        n_lags=10,
        lag_separation_distance=5.0,
        lag_tolerance=3.0,
        direction=Direction(),
        standardize_sill=True,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="gamv",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "gamv")

        if isinstance(exp_vario, (tuple, list)) == False:
            exp_vario = [exp_vario]
        if isinstance(direction, (tuple, list)) == False:
            direction = [direction]
        data_columns = set()
        for vario in exp_vario:
            data_columns.add(vario.tail_column)
            data_columns.add(vario.head_column)
        data_columns = list(data_columns)

        n_directions = 0
        for _direction in direction:
            if isinstance(_direction.azimuth, float) == True:
                n_directions += 1
            else:
                n_directions += len(_direction.azimuth)

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(data, to_list(coord_columns, 3, np.nan)),
                "\\  columns for X, Y, Z coordinates",
            ),
            (
                [len(data_columns)] + _convert_column(data, data_columns),
                "\\  number of variables, column numbers",
            ),
            (trimming_limits, "\\  triming limits"),
            (output_file_name + ".out", "\\file for variogram output"),
            (n_lags, "\\number of lags"),
            (lag_separation_distance, "\\lag separation distance"),
            (lag_tolerance, "\\lag tolerance"),
            (n_directions, "\\number of directions"),
        ]
        for _direction in direction:
            azimuths = _direction.azimuth
            if isinstance(azimuths, float) == True:
                azimuths = [azimuths]
            for azimuth in azimuths:
                params = [
                    azimuth,
                    _direction.azimuth_tolerance,
                    _direction.horizontal_bandwidth,
                    _direction.dip,
                    _direction.dip_tolerance,
                    _direction.vertical_bandwidth,
                ]
                self.parameters.append((params, "\\azm, atol, bandh, dip, dtol, bandv"))
        self.parameters.append((standardize_sill, "\\standardize sill? (O=no, l=yes)"))
        self.parameters.append((len(exp_vario), "\\number of variograms"))
        for vario in exp_vario:
            _vario = [
                data_columns.index(vario.tail_column) + 1,
                data_columns.index(vario.head_column) + 1,
                vario.vario_type,
            ]
            if vario.cutoff is not None:
                _vario.append(vario.cutoff)
            self.parameters.append(
                (_vario, "\\tail variable, head variable, variogram type")
            )


class GAMV(DataManager):
    """
    Variogram calculation with irregular data. This implementation is more
    flexible than GSLIB: the number of lags and the lag distance can vary for
    each direction.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. In 1D or 2D,
        use np.nan or pd.NA to signal the coordinates to ignore. The last and
        second-to-last coordinates can also be omitted instead of set to NaN.
    exp_vario : ExpVario or array-like of shape (n_vario,), default=ExpVario()
        Experimental variogram(s) to compute.
    direction : Direction or array-like of shape (n_dir,), default=Direction()
        Direction(s) along which to compute the variogram.
    standardize_sill : bool, default=True
        If true, standardizes the sill.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='gamv'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        exp_vario=ExpVario(),
        direction=Direction(),
        standardize_sill=True,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="gamv",
        output_dir_path=".",
        executable_path=None,
    ):

        DataManager.__init__(self, output_file_name, output_dir_path)

        if isinstance(direction, (tuple, list)) == False:
            direction = [direction]

        groups = []
        for key, group in itertools.groupby(
            sorted(enumerate(direction), key=self.group_key), key=self.group_key
        ):
            groups.append(list(zip(*group)))

        self.gamv = []
        for i, group in enumerate(groups):
            gamv = _GAMV(
                data,
                coord_columns,
                exp_vario,
                group[1][0].n_lags,
                group[1][0].lag_separation_distance,
                group[1][0].lag_tolerance * group[1][0].lag_separation_distance,
                group[1],
                standardize_sill,
                trimming_limits,
                output_file_name + "_" + str(i + 1),
                output_dir_path,
                executable_path,
            )
            self.gamv.append((group[0], gamv))

    def group_key(self, direction):
        """
        Extracts a grouping key from a tuple with an index and a `Direction` object.
        """
        return (
            direction[1].n_lags,
            direction[1].lag_separation_distance,
            direction[1].lag_tolerance,
        )

    def run(self, parameter_file_path=None, print_log=True, clean_files=True):
        """
        Runs the GSLIB programs.

        Parameters
        ----------
        parameter_file_path : str, default=None
            Path to the parameter file, which is defined in the daughter classes.
        print_log : bool, default=True
            If true, prints the log of the GSLIB program while it's running,
            otherwise doesn't print anything.
        clean_files : bool, default=True
            If true, delete the files created for and by the run, otherwise keep them.

        Returns
        -------
        output : array-like
            The output generated by the program.
        """
        data = []
        for indices, gamv in self.gamv:
            _data = gamv.run(parameter_file_path, print_log, clean_files)
            for index, i_direction in zip(indices, _data["Direction"].unique()):
                _data.loc[_data["Direction"] == i_direction, "Direction"] = index + 1
            data.append(_data)

        return pd.concat(data)


class VMODEL(GSLIB):
    """
    Variogram values from a model.

    Parameters
    ----------
    vario_model : VarioModel, default=VarioModel()
        Variogram model specification.
    direction : Direction or array-like of shape (n_dir,), default=Direction()
        Direction(s) along which to compute the variogram.
    n_lags : int, default=10
        Number of lags.
    output_file_name : str, default='vmodel'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        vario_model=VarioModel(),
        direction=Direction(),
        n_lags=10,
        output_file_name="vmodel",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "vmodel")

        if isinstance(direction, (tuple, list)) == False:
            direction = [direction]

        self.parameters = [
            (output_file_name + ".out", "\\file for variogram output"),
            ((len(direction), n_lags), "\\number of directions and lags"),
        ]
        for _direction in direction:
            self.parameters += [
                (
                    (_direction.azimuth, _direction.dip, _direction.lag_distance),
                    "\\azm, dip, lag distance",
                )
            ]
        self.parameters += [
            (
                (len(vario_model.structures), vario_model.nugget_effect),
                "\\nst, nugget effect",
            )
        ]
        for structure in vario_model.structures:
            self.parameters += [
                (
                    [
                        structure.model,
                        structure.partial_sill,
                        structure.azimuth,
                        structure.dip,
                        structure.plunge,
                    ],
                    "\\it,cc,anql,ang2,ang3",
                )
            ]
            self.parameters += [
                (
                    [
                        "  ",
                        structure.range_hmax,
                        structure.range_hmin,
                        structure.range_vert,
                    ],
                    "\\a_hmax, a_hmin, a_vert",
                )
            ]


class VARMAP(GSLIB):
    """
    Variogram map calculation with regular or irregular data.

    Parameters
    ----------
    data : DataFrame
        Data to transform.
    shape : array-like (x, y, z), default=None
        If the data are on a regular grid, shape of the grid along the x, y, and
        z axes.
    spacing : array-like (x, y, z), default=None
        If the data are on a regular grid, cell size along the x, y, and z axes.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. In 1D or 2D,
        use np.nan or pd.NA to signal the coordinates to ignore. The last and
        second-to-last coordinates can also be omitted instead of set to NaN. To
        use instead of `shape` and `spacing` if the data are scattered.
    exp_vario : ExpVario or array-like of shape (n_vario,), default=ExpVario()
        Experimental variogram(s) to compute.
    n_lags : int or array-like (x, y, z), default=(10, 10, 0)
        Number of lags in the x, y, and z directions.
    lag_tolerence : float or array-like (x, y, z), default=(5., 5., 1.)
        Lag tolerances or “cell sizes” in the x, y, and z directions.
    min_n_pairs : int, default=5
        Minimum number of pairs needed to define a variogram value (set to
        missing if fewer than minpairs is found).
    standardize_sill : bool, default=True
        If true, standardizes the sill.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='varmap'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        shape=None,
        spacing=None,
        coord_columns=(0, 1, np.nan),
        exp_vario=ExpVario(),
        n_lags=(10, 10, 0),
        lag_tolerences=(5.0, 5.0, 1.0),
        min_n_pairs=5,
        standardize_sill=True,
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="varmap",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "varmap")

        is_regular_grid = False
        if shape is not None and spacing is not None:
            is_regular_grid = True
        else:
            shape = 1
            spacing = 1.0
        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        if isinstance(exp_vario, (tuple, list)) == False:
            exp_vario = [exp_vario]
        data_columns = set()
        for vario in exp_vario:
            data_columns.add(vario.tail_column)
            data_columns.add(vario.head_column)
        data_columns = list(data_columns)

        self.n_vars = len(data_columns)

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                [len(data_columns)] + _convert_column(data, data_columns),
                "\\  number of variables: column numbers",
            ),
            (trimming_limits, "\\  triming limits"),
            (is_regular_grid, "\\1=regular grid, 0=scattered values"),
            (shape, "\\if =1: nx, ny, nz"),
            (spacing, "\\        xsiz, ysiz, zsiz"),
            (
                _convert_column(data, to_list(coord_columns, 3, np.nan)),
                "\\if =0: columns for x,y, z coordinates",
            ),
            (output_file_name + ".out", "\\file for variogram output"),
            (to_list(n_lags), "\\nxlag, nylag, nzlag"),
            (to_list(lag_tolerences), "\\dxlag, dylag, dzlag"),
            (min_n_pairs, "\\minimum number of pairs"),
            (standardize_sill, "\\standardize sill? (O=no, l=yes)"),
            (len(exp_vario), "\\number of variograms"),
        ]
        for vario in exp_vario:
            _vario = [
                data_columns.index(vario.tail_column) + 1,
                data_columns.index(vario.head_column) + 1,
                vario.vario_type,
            ]
            if vario.cutoff is not None:
                _vario.append(vario.cutoff)
            self.parameters.append(
                (_vario, "\\tail variable, head variable, variogram type")
            )

    def process_output_files(self):
        """
        Processes the output files before reading.
        """
        for i, line in enumerate(
            fileinput.input(self.output_file_path.with_suffix(".out"), inplace=True)
        ):
            if i == 1:
                line = (
                    " ".join(line.strip("\n").split()[:-1])
                    + " "
                    + str(self.n_vars)
                    + "\n"
                )
            sys.stdout.write(line)


################################################################################
# Kriging


def _compute_search_ellipsoid(max_search_radii, search_ellipsoid_angles, vario_model):
    """
    Computes the search ellipsoid based on (a) variogram(s).
    """
    if isinstance(vario_model, (tuple, list)) == False:
        vario_model = [vario_model]

    vario_angles = np.array(
        [
            [model.structures[-1].azimuth for model in vario_model],
            [model.structures[-1].dip for model in vario_model],
            [model.structures[-1].plunge for model in vario_model],
        ]
    )
    if np.all(vario_angles == vario_angles[:, 0:1]):
        # If there's a single variogram model or all the models share the same
        # azimuth, then use the largest ranges for the search neighborhood.
        max_ranges = [
            max(
                [
                    struct.range_hmax
                    for model in vario_model
                    for struct in model.structures
                ]
            ),
            max(
                [
                    struct.range_hmin
                    for model in vario_model
                    for struct in model.structures
                ]
            ),
            max(
                [
                    struct.range_vert
                    for model in vario_model
                    for struct in model.structures
                ]
            ),
        ]
        _max_search_radii = [
            r if mr is None else mr
            for mr, r in zip(to_list(max_search_radii), max_ranges)
        ]
        _search_ellipsoid_angles = [
            va if a is None else a
            for va, a in zip(to_list(search_ellipsoid_angles), vario_angles[:, 0])
        ]
    else:
        # If the azimuths aren't aligned, compute the bounding box containing
        # all the variogram ellipsoids and use its dimensions for the search
        # neighborhood.
        _max_search_radii = get_ellipsoid_aabb(
            [
                max([struct.range_hmax for struct in model.structures])
                for model in vario_model
            ],
            [
                max([struct.range_hmin for struct in model.structures])
                for model in vario_model
            ],
            [
                max([struct.range_vert for struct in model.structures])
                for model in vario_model
            ],
            azimuth=vario_angles[0],
            dip=vario_angles[1],
            plunge=vario_angles[2],
        )
        _max_search_radii = np.max(_max_search_radii, axis=0)
        _search_ellipsoid_angles = [0.0, 0.0, 0.0]

    return _max_search_radii, _search_ellipsoid_angles


class KT3D(GSLIB):
    """
    2D or 3D kriging for points or blocks by simple kriging (SK), ordinary
    kriging (OK), or kriging with a polynomial trend model (KT) with up to nine
    monomial terms.

    Parameters
    ----------
    data : DataFrame
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : int, str, or array-like (var, sec_var, id), default=(2, np.nan, np.nan)
        Indices or names of the variable to be simulated, the secondary variable
        (e.g., for external drift if used), and the ids for cross-validation in
        `data`. The last and second-to-last variables can also be omitted instead
        of set to NaN.
    vario_model : VarioModel, default=VarioModel()
        Variogram model specification.
    kriging_type : int, default=0
        The kriging type:
            0 = simple kriging;
            1 = ordinary kriging;
            2 = simple kriging with a locally varying mean;
            3 = kriging with an external drift.
        SK is required by theory; only in cases where the number of original data
        found in the neighborhood is large enough can OK be used without the risk
        of spreading data values beyond their range of influence.
    simple_kriging_mean : float, default=None
        Mean value to use with simple kriging. If none, the mean value from the
        data is used.
    include_drift_terms : bool or array-like of shape (9,), default=False
        Indicators for those drift terms to be included in the trend model.
        include_drift_terms(i) is set to `True` if the drift term number i should
        be included, and is set to `False` if not. The nine drift terms correspond
        to the following:
            i = 1 linear drift in x;
            i = 2 linear drift in y;
            i = 3 linear drift in z;
            i = 4 quadratic drift in x;
            i = 5 quadratic drift in y;
            i = 6 quadratic drift in z;
            i = 7 cross quadratic drift in xy;
            i = 8 cross quadratic drift in xz;
            i = 9 cross quadratic drift in yz.
    estimate_trend : bool, default=False,
        If true, estimates the trend, otherwise estimates the variable. The trend
        may be kriged with ordinary kriging (all include_drift_terms(i) values set
        to 0) or with any combination of trend kriging (some include_drift_terms(i)
        terms set to 1).
    secondary_data : array-like of shape (n_var, n_z, n_y, n_x), default=None
        The locally varying mean or the external drift variable (the secondary
        variable must be gridded at the same resolution as the model being
        constructed by `kt3d`).
    secondary_data_column : int, default=4
        Index or name of the variable in `secondary_data`.
    use_case : int, default=0
        The use case:
            0 = kriging a grid of points or blocks;
            1 = cross validation with the data in `data`;
            2 = jackknifing with data in `jackknife_data`.
    jackknife_data : DataFrame, default=None
        Conditioning data.
    jackknife_coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `jackknife_data`. One
        or two of the coordinates can be set to np.nan or pd.NA, which indicates
        that the simulation is 2D or 1D. The last and second-to-last coordinates
        can also be omitted instead of set to NaN.
    jackknife_data_columns : int, str, or array-like (var, weight, sec_var), default=(2, np.nan)
        Indices or names of the variable to be simulated and the secondary
        variable (e.g., for external drift if used) in `jackknife_data`. The last
        and second-to-last variables can also be omitted instead of set to NaN.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    block_discretization : float or array-like (x, y, z), default=1
        Number of discretization points for a block. If the three values are all
        set to 1, then point kriging is performed.
    n_data_range : array-like (min, max), default=(4, 8)
        The minimum and maximum number of original data that should be used for
        kriging a block.
    max_data_per_octant : int, default=None
        The number of original data to use per octant. If None, then it is not used;
        otherwise, it overrides the max in `n_data_range` and the data are partitioned
        into octants and the closest `max_data_per_octant` data in each octant are
        retained for the simulation of a grid node.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If any of those
        directions are None, uses the maximum range of the variogram model.
    search_ellipsoid_angles : float or array-like (x, y, z), default=None
        The angle parameters that describe the orientation of the search ellipsoid.
        If any of those angles are None, uses the angles of the variogram model.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default=('Variable', 'Variance')
        Name of the output mean estimate and of the variance.
    output_file_name : str, default='kt3d'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        data_columns=(2, np.nan, np.nan),
        vario_model=VarioModel(),
        kriging_type=0,
        simple_kriging_mean=None,
        include_drift_terms=False,
        estimate_trend=False,
        secondary_data=None,
        secondary_data_column=4,
        use_case=0,
        jackknife_data=None,
        jackknife_coord_columns=(0, 1, np.nan),
        jackknife_data_columns=(2, np.nan),
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        block_discretization=1,
        n_data_range=(4, 8),
        max_data_per_octant=None,
        max_search_radii=None,
        search_ellipsoid_angles=None,
        trimming_limits=(-1.0e21, 1.0e21),
        debugging_level=0,
        variable_name=("Variable", "Variance"),
        output_file_name="kt3d",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "kt3d")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        data_columns = to_list(data_columns, 3, np.nan)
        if simple_kriging_mean is None:
            if isinstance(data_columns[0], str):
                simple_kriging_mean = np.mean(data[data_columns[0]])
            else:
                simple_kriging_mean = np.mean(data[:, data_columns[0]])
        max_data_per_octant = (
            max_data_per_octant if max_data_per_octant is not None else 0
        )
        max_search_radii, search_ellipsoid_angles = _compute_search_ellipsoid(
            max_search_radii, search_ellipsoid_angles, vario_model
        )

        self.spacing = spacing
        self.origin = origin
        self.use_case = use_case

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data,
                    data_columns[-1:]
                    + to_list(coord_columns, 3, np.nan)
                    + data_columns[:-1],
                ),
                "\\  columns for X, Y, Z, var, sec var",
            ),
            (trimming_limits, "\\  trimming limits"),
            (use_case, "\\option: 0=grid, 1=cross, 2=jackknife"),
            (jackknife_data, "\\file with data", output_file_name + "_xvk.dat"),
            (
                _convert_column(
                    jackknife_data,
                    to_list(jackknife_coord_columns, 3, np.nan)
                    + to_list(jackknife_data_columns, 2, np.nan),
                ),
                "\\  columns for X, Y, Z, var, sec var",
            ),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for kriged output"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (to_list(block_discretization), "\\x,y and z block discretization"),
            (n_data_range, "\\min, max data for kriging"),
            (max_data_per_octant, "\\max per octant (0-> not used)"),
            (max_search_radii, "\\maximum search radii"),
            (search_ellipsoid_angles, "\\angles for search ellipsoid"),
            ((kriging_type, simple_kriging_mean), "\\0=SK,1=OK,2=non-st SK,3=exdrift"),
            (to_list(include_drift_terms, 9), "\\drift: x,y,z,xx,yy,zz,xy,xz,zy"),
            (estimate_trend, "\\0, variable; 1, estimate trend"),
            (
                secondary_data,
                "\\gridded file with drift/mean",
                output_file_name + "_extdrift.dat",
            ),
            (
                _convert_column(secondary_data, secondary_data_column),
                "\\  column number in gridded file",
            ),
            (
                (len(vario_model.structures), vario_model.nugget_effect),
                "\\nst, nugget effect",
            ),
        ]
        for structure in vario_model.structures:
            self.parameters.append(
                (
                    [
                        structure.model,
                        structure.partial_sill,
                        structure.azimuth,
                        structure.dip,
                        structure.plunge,
                    ],
                    "\\it,cc,ang1,ang2,ang3",
                )
            )
            self.parameters.append(
                (
                    [
                        "  ",
                        structure.range_hmax,
                        structure.range_hmin,
                        structure.range_vert,
                    ],
                    "\\a_hmax, a_hmin, a_vert",
                )
            )

    def process_output_files(self):
        """
        Processes the output files before reading.
        """
        for i, line in enumerate(
            fileinput.input(self.output_file_path.with_suffix(".out"), inplace=True)
        ):
            if i == 1:
                if self.use_case == 0:
                    system = self.origin + self.spacing
                    line = line.strip("\n") + " " + " ".join(map(str, system)) + "\n"
            sys.stdout.write(line)


class COKB3D(GSLIB):
    """
    2D or 3D cokriging.

    Parameters
    ----------
    data : DataFrame
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : array-like of shape (n_vars,), default=(3, 4)
        Indices or names of the variables to be simulated in `data`.
    vario_models : array-like, default=VarioModel()
        Variogram models for each combination of variables in the right order.
        E.g., with 2 variables, the models represents variable 1 with 1,
        variable 1 with 2, variable 2 with 2.
    kriging_type : int, default=0
        The kriging type:
            0 = simple cokriging;
            1 = standardized ordinary cokriging with recentered variables and a
                single unbiasedness constraint;
            2 = traditional ordinary cokriging.
    simple_kriging_means : array-like of shape (n_vars,), default=None
        Mean values to use with simple kriging. If none, the mean values from the
        data are used.
    colocated_secondary_data : array-like of shape (n_var, n_z, n_y, n_x), default=None
        The gridded covariate (the secondary variable must be gridded at the same
        resolution as the model being constructed by `kt3d`).
    colocated_secondary_data_column : int, default=4
        Index or name of the variable in `colocated_secondary_data`.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    block_discretization : int or array-like (x, y, z), default=1
        Number of discretization points for a block. If the three values are all
        set to 1, then point kriging is performed.
    n_data_range : array-like (min, max), default=(1, 12)
        The minimum and maximum number of original data that should be used for
        kriging a block.
    max_secondary_data : int, default=8
        Maximum number of secondary data that should be used for kriging a block.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If any of those
        directions are None, uses the maximum range of the variogram model.
    max_secondary_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below) for the secondary
        data. If any of those directions are None, uses the maximum range of the
        variogram model.
    search_ellipsoid_angles : float or array-like (x, y, z), default=None
        The angle parameters that describe the orientation of the search ellipsoid.
        If any of those angles are None, uses the angles of the variogram model.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default=('Variable', 'Variance')
        Name of the output mean estimate and of the variance.
    output_file_name : str, default='cokb3d'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        data_columns=(3, 4),
        vario_models=(VarioModel(), VarioModel(), VarioModel()),
        kriging_type=0,
        simple_kriging_means=None,
        colocated_secondary_data=None,
        colocated_secondary_data_column=4,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        block_discretization=1,
        n_data_range=(1, 12),
        max_secondary_data=8,
        max_search_radii=None,
        max_secondary_search_radii=None,
        search_ellipsoid_angles=None,
        trimming_limits=(-1.0e21, 1.0e21),
        debugging_level=0,
        variable_name=("Variable", "Variance"),
        output_file_name="cokb3d",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "cokb3d")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        _simple_kriging_means = []
        for i in range(len(data_columns)):
            if simple_kriging_means is None or simple_kriging_means[i] is None:
                if isinstance(data_columns[i], str):
                    _simple_kriging_means.append(np.mean(data[data_columns[i]]))
                else:
                    _simple_kriging_means.append(np.mean(data[:, data_columns[i]]))
            else:
                _simple_kriging_means.append(simple_kriging_means[i])
        max_search_radii, search_ellipsoid_angles = _compute_search_ellipsoid(
            max_search_radii, search_ellipsoid_angles, vario_models
        )
        max_secondary_search_radii = to_list(max_secondary_search_radii)
        max_secondary_search_radii = [
            sr if sr is not None else r
            for sr, r in zip(max_secondary_search_radii, max_search_radii)
        ]

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (len(data_columns), "\\  number of variables primary+other"),
            (
                _convert_column(
                    data, to_list(coord_columns, 3, np.nan) + list(data_columns)
                ),
                "\\  columns for X,Y,Z and variables",
            ),
            (trimming_limits, "\\  trimming limits"),
            (
                True if colocated_secondary_data is not None else False,
                "\\co-located cokriging? (0=no, 1=yes)",
            ),
            (
                colocated_secondary_data,
                "\\  file with gridded covariate",
                output_file_name + "_extdrift.dat",
            ),
            (
                _convert_column(
                    colocated_secondary_data, colocated_secondary_data_column
                ),
                "\\  column for covariate",
            ),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for kriged output"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (to_list(block_discretization), "\\x, y, and z block discretization"),
            (
                list(n_data_range) + [max_secondary_data],
                "\\min primary,max primary,max all sec",
            ),
            (max_search_radii, "\\maximum search radii: primary"),
            (max_secondary_search_radii, "\\maximum search radii: all secondary"),
            (search_ellipsoid_angles, "\\angles for search ellipsoid"),
            (kriging_type, "\\kriging type (0=SK, 1=OK, 2=OK-trad)"),
            (_simple_kriging_means, "\\mean(i),i=1,nvar"),
        ]
        k = 0
        for i in range(len(data_columns)):
            for j in range(i, len(data_columns)):
                self.parameters += [
                    ((i + 1, j + 1), '''\\semivariogram for "i" and "j"'''),
                    (
                        (
                            len(vario_models[k].structures),
                            vario_models[k].nugget_effect,
                        ),
                        "\\  nst, nugget effect",
                    ),
                ]
                for structure in vario_models[k].structures:
                    self.parameters += [
                        (
                            [
                                structure.model,
                                structure.partial_sill,
                                structure.azimuth,
                                structure.dip,
                                structure.plunge,
                            ],
                            "\\  it,cc,ang1,ang2,ang3",
                        )
                    ]
                    self.parameters += [
                        (
                            [
                                "  ",
                                structure.range_hmax,
                                structure.range_hmin,
                                structure.range_vert,
                            ],
                            "\\  a_hmax, a_hmin, a_vert",
                        )
                    ]
                k += 1


class IK3D(GSLIB):
    """
    Ordinary or simple indicator kriging of either categorical or cdf-type
    indicators.

    Parameters
    ----------
    data : DataFrame
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : int, str, or array-like (var, id), default=(2, np.nan)
        Index or name of the variable to be simulated and of the ids for
        cross-validation in `data`. The ids can be set to NaN when not doing a
        cross-validation or omitted.
    vario_model : VarioModel, default=VarioModel()
        Variogram model specification.
    variable_type : str, default='categorical'
        Type of variable, either `continuous` or `categorical`.
    classes : array-like of shape (n_classes,), default=(0, 1)
        Threshold values or category codes.
    distribution : array-like of shape (n_classes,), default=(0.75, 0.25)
        Global cdf or pdf values.
    kriging_type : int, default=0
        The kriging type used throughout the loop over all node:
            0 = simple kriging;
            1 = ordinary kriging.
    use_full_ind_kriging : bool, default=True
        If True, a full indicator kriging is performed at each grid node location
        to establish the conditional distribution. Otherwise, the median approximation
        is used; i.e., a single variogram is used for all categories; therefore,
        only one kriging system needs to be solved and the computer time is
        significantly reduced.
    median_approx_class : int or float, default=0
        Class whose variogram is used when `use_full_ind_kriging` is False, i.e.,
        when the median approximation is used.
    indicator_data : DataFrame, default=None
        Already transformed soft indicator values. Missing values are identified as
        less than the minimum in `extrema`, which would correspond to a constraint
        interval. Otherwise, the cdf data should steadily increase from 0 to 1
        and soft categorical probabilities must be between 0 to 1 and sum to 1.
    indicator_coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `indicator_data`. One
        or two of the coordinates can be set to np.nan or pd.NA, which indicates
        that the simulation is 2D or 1D. The last and second-to-last coordinates
        can also be omitted instead of set to NaN.
    indicator_data_columns : array-like of shape (n_classes,), default=(2, 3)
        Indices or names of the indicator variables to be simulated in `indicator_data`.
    use_case : int, default=0
        The use case:
            0 = kriging a grid of points or blocks;
            1 = cross validation with the data in `data`;
            2 for jackknifing with data in `jackknife_data`.
    jackknife_data : DataFrame, default=None
        Conditioning data.
    jackknife_coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `jackknife_data`. One
        or two of the coordinates can be set to np.nan or pd.NA, which indicates
        that the simulation is 2D or 1D. The last and second-to-last coordinates
        can also be omitted instead of set to NaN.
    jackknife_data_column : int or str, default=2
        Index or name of the variable to be simulated in `jackknife_data`.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_data_range : array-like (min, max), default=(4, 8)
        The minimum and maximum number of original data that should be used for
        kriging a block.
    max_data_per_octant : int, default=None
        The number of original data to use per octant. If None, then it is not used;
        otherwise, it overrides the max in `n_data_range` and the data are partitioned
        into octants and the closest `max_data_per_octant` data in each octant are
        retained for the simulation of a grid node.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If any of those
        directions are None, uses the maximum range of the variogram model.
    search_ellipsoid_angles : float or array-like (x, y, z), default=None
        The angle parameters that describe the orientation of the search ellipsoid.
        If any of those angles are None, uses the angles of the variogram model.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default=None
        Name of the output variable for each category.
    output_file_name : str, default='ik3d'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        data,
        coord_columns=(0, 1, np.nan),
        data_columns=(2, np.nan),
        vario_model=(VarioModel(), VarioModel()),
        variable_type="categorical",
        classes=(0, 1),
        distribution=(0.75, 0.25),
        kriging_type=0,
        use_full_ind_kriging=True,
        median_approx_class=0,
        indicator_data=None,
        indicator_coord_columns=(0, 1, np.nan),
        indicator_data_columns=(2, 3),
        use_case=0,
        jackknife_data=None,
        jackknife_coord_columns=(0, 1, np.nan),
        jackknife_data_column=2,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        n_data_range=(4, 8),
        max_data_per_octant=None,
        max_search_radii=None,
        search_ellipsoid_angles=None,
        trimming_limits=(-1.0e21, 1.0e21),
        debugging_level=0,
        variable_name=None,
        output_file_name="ik3d",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "ik3d")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        data_columns = to_list(data_columns, 2, np.nan)
        variable_types = {"categorical": 0, "continuous": 1}
        variable_type = variable_types[variable_type]
        if isinstance(vario_model, (tuple, list)) == False:
            vario_model = len(classes) * [vario_model]
        max_data_per_octant = (
            max_data_per_octant if max_data_per_octant is not None else 0
        )
        max_search_radii, search_ellipsoid_angles = _compute_search_ellipsoid(
            max_search_radii, search_ellipsoid_angles, vario_model
        )

        self.parameters = [
            (variable_type, "\\1=continuous(cdf), 0=categorical(pdf)"),
            (use_case, "\\option: 0=grid, 1=cross, 2=jackknife"),
            (jackknife_data, "\\file with data", output_file_name + "_jack.dat"),
            (
                _convert_column(
                    jackknife_data,
                    to_list(jackknife_coord_columns, 3, np.nan)
                    + [jackknife_data_column],
                ),
                "\\  columns for X,Y,Z,vr",
            ),
            (len(classes), "\\number thresholds/categories"),
            (classes, "\\  thresholds / categories"),
            (distribution, "\\  global cdf / pdf"),
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data,
                    data_columns[1:]
                    + to_list(coord_columns, 3, np.nan)
                    + list(data_columns[:1]),
                ),
                "\\  columns for X,Y,Z and variable",
            ),
            (indicator_data, "\\file with soft indicator input", "direct.ik"),
            (
                _convert_column(
                    indicator_data,
                    to_list(indicator_coord_columns, 3, np.nan)
                    + list(indicator_data_columns),
                ),
                "\\  columns for X,Y,Z, and indicators",
            ),
            (trimming_limits, "\\  trimming limits"),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for kriged output"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (n_data_range, "\\min, max data for kriging"),
            (max_search_radii, "\\maximum search radii"),
            (search_ellipsoid_angles, "\\angles for search ellipsoid"),
            (max_data_per_octant, "\\max per octant (0-> not used)"),
            (
                (not use_full_ind_kriging, median_approx_class),
                "\\0=full IK, 1=Median IK(threshold num)",
            ),
            (kriging_type, "\\0=SK, 1=OK"),
        ]
        for i, model in enumerate(vario_model):
            self.parameters += [
                (
                    (len(model.structures), model.nugget_effect),
                    "\\" + str(i + 1) + "  nst, nugget effect",
                )
            ]
            for structure in model.structures:
                self.parameters += [
                    (
                        [
                            structure.model,
                            structure.partial_sill,
                            structure.azimuth,
                            structure.dip,
                            structure.plunge,
                        ],
                        "\\   it,cc,anql,ang2,ang3",
                    )
                ]
                self.parameters += [
                    (
                        [
                            "  ",
                            structure.range_hmax,
                            structure.range_hmin,
                            structure.range_vert,
                        ],
                        "\\   a_hmax, a_hmin, a_vert",
                    )
                ]


################################################################################
# Simulations


class LUSIM(GSLIB):
    """
    LU Simulation. This program requires standard normal data and writes standard
    normal simulated values. Normal score transforms and back transforms are to
    be performed outside this program.

    Parameters
    ----------
    vario_model : VarioModel, default=VarioModel()
        Variogram model specification.
    data : DataFrame, default=None
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : int, str, or array-like (var, weight, sec_var), default=(2, np.nan, np.nan)
        Indices or names of the variable to be simulated, the declustering weight,
        and the secondary variable (e.g., for external drift if used) in `data`.
        The last and second-to-last variables can be set to NaN or omitted.
    shape : int or array-like (x, y, z), default=(50, 50, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(2., 2., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(1., 1., 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    seed : int, default=42
        Seed for random number generation.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default='Variable'
        Name of the simulated variable.
    output_file_name : str, default='sgsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        vario_model=VarioModel(),
        data=None,
        coord_columns=(0, 1, np.nan),
        data_column=2,
        shape=(50, 50, 1),
        spacing=(2.0, 2.0, 1.0),
        origin=(1.0, 1.0, 0.0),
        n_realizations=1,
        trimming_limits=(-1.0e21, 1.0e21),
        seed=42,
        debugging_level=0,
        variable_name="Variable",
        output_file_name="lusim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "lusim")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data, to_list(coord_columns, 3, np.nan) + [data_column]
                ),
                "\\  columns for X,Y,Z, normal scores",
            ),
            (trimming_limits, "\\  trimming limits"),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for realization(s)"),
            (n_realizations, "\\number of realizations"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (seed, "\\random number seed"),
            (
                (len(vario_model.structures), vario_model.nugget_effect),
                "\\nst, nugget effect",
            ),
        ]
        for structure in vario_model.structures:
            self.parameters.append(
                (
                    [
                        structure.model,
                        structure.partial_sill,
                        structure.azimuth,
                        structure.dip,
                        structure.plunge,
                    ],
                    "\\it,cc,ang1,ang2,ang3",
                )
            )
            self.parameters.append(
                (
                    [
                        "  ",
                        structure.range_hmax,
                        structure.range_hmin,
                        structure.range_vert,
                    ],
                    "\\a_hmax, a_hmin, a_vert",
                )
            )


class SGSIM(GSLIB):
    """
    Sequential Gaussian Simulation.

    Parameters
    ----------
    vario_model : VarioModel, default=VarioModel()
        Variogram model specification.
    data : DataFrame, default=None
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : int, str, or array-like (var, weight, sec_var), default=(2, np.nan, np.nan)
        Indices or names of the variable to be simulated, the declustering weight,
        and the secondary variable (e.g., for external drift if used) in `data`.
        The last and second-to-last variables can be set to NaN or omitted.
    kriging_type : int, default=0
        The kriging type:
            0 = simple kriging;
            1 = ordinary kriging;
            2 = simple kriging with a locally varying mean;
            3 = kriging with an external drift;
            4 = collocated cokriging with one secondary variable.
        SK is required by theory; only in cases where the number of original data
        found in the neighborhood is large enough can OK be used without the risk
        of spreading data values beyond their range of influence.
    secondary_data : array-like of shape (n_var, n_rez, n_x, n_y, n_z), default=None
        The locally varying mean, the external drift variable, or the secondary
        variable for collocated cokriging (the secondary variable must be gridded
        at the same resolution as the model being constructed by `SGSIM`).
    secondary_data_column : int, default=4
        Index or name of the variable in `secondary_data`.
    correlation_coefficient : float, default=0.6
        Correlation coefficient to use for collocated cokriging (used only if
        `kriging_type` is 4).
    variance_reduction_factor : float, default=1.
        Variance reduction factor to use for collocated cokriging (used only if
        `kriging_type` is 4). It modifies the kriging variance after collocated
        cokriging. The default should be 1.; however, depending on the continuity
        of the secondary variable realization, the variance of the resulting model
        can be too high. A reduction factor of about 0.6 may be required.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    n_data_range : array-like (min, max), default=(0, 8)
        The minimum and maximum number of original data that should be used to
        simulate a grid node. If there are fewer than ndmin data points, the node
        is not simulated.
    max_prev_sim_nodes : int, default=12
        The maximum number of previously simulated nodes to use for the simulation
        of another node.
    search_separately : bool, default=False
        If true, the data and previously simulated grid nodes are searched
        separately: The data are searched with a super block search, and the
        previously simulated nodes are searched with a spiral search. If false,
        the data are relocated to grid nodes and a spiral search is used and the
        parameter `max_prev_sim_nodes` is not considered.
    n_multigrids : int, default=0
        The number of multiple grid refinements to consider (used only if set to
        more than 0).
    max_data_per_octant : int, default=None
        The number of original data to use per octant. If None, then it is not used;
        otherwise, it overrides the max in `n_data_range` and the data are partitioned
        into octants and the closest `max_data_per_octant` data in each octant are
        retained for the simulation of a grid node.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If any of those
        directions are None, uses the maximum range of the variogram model.
    search_ellipsoid_angles : float or array-like (x, y, z), default=None
        The angle parameters that describe the orientation of the search ellipsoid.
        If any of those angles are None, uses the angles of the variogram model.
    shape_cov_lookup_table : int or array-like (x, y, z), default=361
        Shape of the of covariance lookup table.
    transform_data : bool, default=True
        If true, apply a normal score transformation to the data before simulation
        and back transform the realizations. Otherwise, no transformation is applied.
    smoothed_distribution : array-like of shape (1,) or (2,), default=None
        A smoothed distribution for the transformation. This can be useful when
        the number of samples is limited and the original distribution is deemed
        too discrete so not representative.
    smoothed_distribution_columns : array-like of shape (2,), default=(1, np.nan)
        Indices or names of the variable and weight in `smoothed_distribution`.
        If the weight is np.nan, pd.NA, or omitted, then equal weighting is applied.
    extrema : array-like (min, max), default=None
        Minimum and maximum values to be used for extrapolation in the tails.
    lower_tail : array-like of shape (2,), default=None
        The first value specifies the back-transformation implementation in the
        lower tail of the distribution:
        1 = linear interpolation to the lower extrema;
        2 = power model interpolation to the lower extrema.
        The second value specifies the parameter of option 2.
    upper_tail : array-like of shape (2,), default=None
        The first value specifies the back-transformation implementation in the
        upper tail of the distribution:
        1 = linear interpolation to the upper extrema;
        2 = power model interpolation to the upper extrema;
        4 = hyperbolic model extrapolation to the upper extrema.
        The second value specifies the parameter of option 2 or 4.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    seed : int, default=42
        Seed for random number generation.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default='Variable'
        Name of the simulated variable.
    output_file_name : str, default='sgsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        vario_model=VarioModel(),
        data=None,
        coord_columns=(0, 1, np.nan),
        data_columns=(2, np.nan, np.nan),
        kriging_type=0,
        secondary_data=None,
        secondary_data_column=4,
        correlation_coefficient=0.6,
        variance_reduction_factor=1.0,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        n_realizations=1,
        n_data_range=(0, 8),
        max_prev_sim_nodes=12,
        search_separately=False,
        n_multigrids=0,
        max_data_per_octant=None,
        max_search_radii=None,
        search_ellipsoid_angles=None,
        shape_cov_lookup_table=361,
        transform_data=False,
        smoothed_distribution=None,
        smoothed_distribution_columns=(0, np.nan),
        extrema=None,
        lower_tail=None,
        upper_tail=None,
        trimming_limits=(-1.0e21, 1.0e21),
        seed=42,
        debugging_level=0,
        variable_name="Variable",
        output_file_name="sgsim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "sgsim")
        self.variable_name = variable_name

        # kriging_type = {'simple': 0, 'ordinary': 1, 'varying mean': 2, 'external drift': 3, 'collocated': 4}

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        max_data_per_octant = (
            max_data_per_octant if max_data_per_octant is not None else 0
        )
        max_search_radii, search_ellipsoid_angles = _compute_search_ellipsoid(
            max_search_radii, search_ellipsoid_angles, vario_model
        )
        data_columns = to_list(data_columns, 3, np.nan)
        if data is None:
            vmin, vmax = -4.0, 4.0
        else:
            vmin, vmax = data[data_columns[0]].min(), data[data_columns[0]].max()
        if extrema is None:
            extrema = (vmin, vmax)
        if lower_tail is None:
            lower_tail = (1, vmin)
        if upper_tail is None:
            upper_tail = (1, vmax)
        if isinstance(secondary_data, xr.Dataset):
            secondary_data = secondary_data.drop_vars(
                [
                    var
                    for var in secondary_data.data_vars
                    if var != secondary_data_column
                ]
            )
        elif secondary_data is not None and isinstance(secondary_data_column, int):
            secondary_data = secondary_data[
                secondary_data_column : secondary_data_column + 1
            ]

        self.parameters = [
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(data, to_list(coord_columns, 3, np.nan) + data_columns),
                "\\  columns for X,Y,Z,vr,wt,sec.var.",
            ),
            (trimming_limits, "\\  trimming limits"),
            (transform_data, "\\transform the data (0=no, 1=yes)"),
            (output_file_name + ".trn", "\\  file for output trans table"),
            (
                False if smoothed_distribution is None else True,
                "\\  consider ref. distribution (0=no, 1=yes)",
            ),
            (
                smoothed_distribution,
                "\\  file with ref. distribution",
                output_file_name + "_smooth_distrib.dat",
            ),
            (
                _convert_column(
                    smoothed_distribution,
                    to_list(smoothed_distribution_columns, 2, np.nan),
                ),
                "\\  columns for vr and wt",
            ),
            (extrema, "\\  zmin,zmax(tail extrapolation)"),
            (lower_tail, "\\  lower tail option, parameter"),
            (upper_tail, "\\  upper tail option, parameter"),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for simulation output"),
            (n_realizations, "\\number of realizations to generate"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (seed, "\\random number seed"),
            (n_data_range, "\\min and max original data for sim"),
            (max_prev_sim_nodes, "\\number of simulated nodes to use"),
            (not search_separately, "\\assign data to nodes (O=no,l=yes)"),
            (
                (True if n_multigrids > 0 else False, n_multigrids),
                "\\multiple grid search (O=no,1-yes),num",
            ),
            (max_data_per_octant, "\\maximum data per octant (O=not used)"),
            (max_search_radii, "\\maximum search radii (hmax,hmin,vert)"),
            (search_ellipsoid_angles, "\\angles for search ellipsoid"),
            (to_list(shape_cov_lookup_table), "\\size of covariance lookup table"),
            (
                (kriging_type, correlation_coefficient, variance_reduction_factor),
                "\\ktype: 0=SK,1=OK,2=FLVM,3=EXDR,4=COLC",
            ),
            (
                secondary_data,
                "\\  file with LVM, EXDR, or COLC variable",
                output_file_name + "_secondary_data.dat",
            ),
            (
                _convert_column(secondary_data, secondary_data_column),
                "\\  column for secondary variable",
            ),
            (
                (len(vario_model.structures), vario_model.nugget_effect),
                "\\nst, nugget effect",
            ),
        ]
        for structure in vario_model.structures:
            self.parameters.append(
                (
                    [
                        structure.model,
                        structure.partial_sill,
                        structure.azimuth,
                        structure.dip,
                        structure.plunge,
                    ],
                    "\\it,cc,anql,ang2,ang3",
                )
            )
            self.parameters.append(
                (
                    [
                        "  ",
                        structure.range_hmax,
                        structure.range_hmin,
                        structure.range_vert,
                    ],
                    "\\a_hmax, a_hmin, a_vert",
                )
            )


class BICALIB(GSLIB):
    """
    Bivariate calibration.

    Parameters
    ----------
    secondary_data : DataFrame
        Full secondary data.
    calibration_data : DataFrame
        Calibration data containing pairs of primary and secondary data values.
    secondary_data_column : int or str, default=2
        Index or name of the secondary variable in `secondary_data`.
    calibration_data_columns : array-like (prim, sec, weight), default=(2, 3, 4)
        Index or name of the primary and secondary variables and declustering
        weight/bivariate frequency in `secondary_data`.
    thresholds_primary : array-like of shape (n_thres), default=(0.5, 1., 2.5, 5., 10.)
        Threshold values applied to the primary variable.
    thresholds_secondary : array-like of shape (n_thres), default=(0.5, 1., 2.5, 5., 10.)
        Threshold values applied to the secondary variable.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    output_file_name : str, default='bicalib'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        secondary_data,
        calibration_data,
        secondary_data_column=4,
        calibration_data_columns=(2, 3, 4),
        thresholds_primary=(0.5, 1.0, 2.5, 5.0, 10.0),
        thresholds_secondary=(0.5, 1.0, 2.5, 5.0, 10.0),
        trimming_limits=(-1.0e21, 1.0e21),
        output_file_name="bicalib",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "bicalib")

        self.parameters = [
            (
                secondary_data,
                "\\file with secondary data",
                output_file_name + "_ydata.dat",
            ),
            (
                _convert_column(secondary_data, secondary_data_column),
                "\\  column for secondary variable",
            ),
            (
                calibration_data,
                "\\file with calibration scatterplot",
                output_file_name + "_data.dat",
            ),
            (
                _convert_column(calibration_data, calibration_data_columns),
                "\\  columns of pri, sec, and weight",
            ),
            (trimming_limits, "\\  triming limits"),
            (output_file_name + ".out", "\\file for output data / distributions"),
            (output_file_name + ".cal", "\\file for output calibration (SISIM)"),
            (output_file_name + ".rep", "\\file for calibration report"),
            (len(thresholds_primary), "\\number of thresholds on primary"),
            (thresholds_primary, "\\  thresholds on primary"),
            (len(thresholds_secondary), "\\number of thresholds on secondary"),
            (thresholds_secondary, "\\  thresholds on secondary"),
        ]

        self.output_file_paths += [self.output_file_path.with_suffix(".cal")]


class SISIM(GSLIB):
    """
    Sequential Indicator Simulation.

    Parameters
    ----------
    vario_model : VarioModel or array-like of shape (n_models,), default=(VarioModel(), VarioModel())
        Variogram model specification for each class.
    variable_type : str, default='categorical'
        Type of variable, either `continuous` or `categorical`.
    class_codes : array-like of shape (n_classes,), default=(0, 1)
        Threshold values or category codes.
    proportions : array-like of shape (n_classes,), default=(0.75, 0.25)
        Global cdf or pdf values.
    data : DataFrame, default=None
        Conditioning data.
    coord_columns : int, str, array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_column : int or str, default=2
        Index or name of the variable to be simulated in `data`.
    kriging_type : int, default=0
        The kriging type used throughout the loop over all node:
            0 = simple kriging;
            1 = ordinary kriging.
        SK is required by theory, only in cases where the number of original data
        found in the neighborhood is large enough can OK be used without the risk
        of spreading data values beyond their range of influence. The global pdf
        values (specified with each category) are used for simple kriging.
    use_full_ind_kriging : bool, default=True
        If True, a full indicator kriging is performed at each grid node location
        to establish the conditional distribution. Otherwise, the median approximation
        is used; i.e., a single variogram is used for all categories; therefore,
        only one kriging system needs to be solved and the computer time is
        significantly reduced.
    median_approx_class : int or float, default=0
        Class whose variogram is used when `use_full_ind_kriging` is False, i.e.,
        when the median approximation is used.
    indicator_data : DataFrame, default=None
        Already transformed soft indicator values. Missing values are identified as
        less than the minimum in `extrema`, which would correspond to a constraint
        interval. Otherwise, the cdf data should steadily increase from 0 to 1
        and soft categorical probabilities must be between 0 to 1 and sum to 1.
    indicator_coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `indicator_data`. One
        or two of the coordinates can be set to np.nan or pd.NA, which indicates
        that the simulation is 2D or 1D. The last and second-to-last coordinates
        can also be omitted instead of set to NaN.
    indicator_data_columns : array-like of shape (n_classes,), default=(2, 3)
        Indices or names of the indicator variables to be simulated in `indicator_data`.
    use_markov_bayes_sim : bool, default=False
        If True, uses the Markov-Bayes option for cokriging with soft indicator
        data. Otherwise, doesn't use it.
    calibration : array-like of shape (n_classes,) or dict, default=(0.61, 0.54)
        Calibration required by the Markov-Bayes option if `use_markov_bayes_sim`
        is True. If `indicator_data` is used, `calibration` is an array-like
        with the calibration parameters; if `secondary_data` is used, `calibration`
        is a dict with the thresholds, the calibration table, and the calibration
        parameters.
    secondary_data : array-like of shape (n_var, n_x, n_y, n_z), default=None
        Gridded secondary variable to be used instead of `indicator_data`.
    secondary_data_column : int, default=4
        Index or name of the variable in `secondary_data`.
    prior_mean_data : array-like of shape (n_classes, n_x, n_y, n_z), default=None
        Prior cdf or pdf mean values to be used instead of `indicator_data`. Prior
        cdf or pdf information can be a valuable vehicle to incorporate soft knowledge
        regarding the phenomenon under study. They can also be used to incorporate
        in the simulation prior guesses for the lithofacies proportions, for example,
        stemming from geological interpretation or seismic imaging. These prior
        guesses could be then updated via IK in a Bayesian framework, i.e. the prior
        local ccdfs are updated to account for the hard data available in their
        neighborhoods.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    max_original_data : int, default=12
        The maximum number of original data that should be used to simulate a
        grid node.
    max_prev_sim_nodes : int, default=12
        The maximum number of previously simulated nodes to use for the simulation
        of another node.
    max_soft_ind_nodes : int, default=1
        The maximum soft indicator nodes for kriging, only used with secondary data.
    search_separately : bool, default=False
        If false, the data and previously simulated grid nodes are searched
        separately: The data are searched with a super block search, and the
        previously simulated nodes are searched with a spiral search. If true,
        the data are relocated to grid nodes and a spiral search is used and the
        parameter `max_prev_sim_nodes` is not considered.
    n_multigrids : int, default=0
        The number of multiple grid refinements to consider (used only if set to
        more than 0).
    max_data_per_octant : int, default=None
        The number of original data to use per octant. If None, then it is not used;
        otherwise, it overrides the max in `n_data_range` and the data are partitioned
        into octants and the closest `max_data_per_octant` data in each octant are
        retained for the simulation of a grid node.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If any of those
        directions are None, uses the maximum range of the variogram model.
    search_ellipsoid_angles : float or array-like (x, y, z), default=None
        The angle parameters that describe the orientation of the search ellipsoid.
        If any of those angles are None, uses the angles of the variogram model.
    shape_cov_lookup_table : int or array-like (x, y, z), default=361
        Shape of the of covariance lookup table.
    extrema : array-like (min, max), default=None
        Minimum and maximum values to be used for extrapolation in the tails.
    lower_tail : array-like of shape (2,), default=None
        The first value specifies the extrapolation in the lower tail of the distribution:
        1 = linear interpolation to the lower extrema;
        2 = power model interpolation to the lower extrema;
        3 = linear interpolation between tabulated quantiles (only for continuous
            variables).
        The second value specifies the parameter of option 2.
    middle_tail : array-like of shape (2,), default=None
        The first value specifies the interpolation within the middle of the
        distribution:
        1 = linear interpolation;
        2 = power model interpolation;
        3 = linear interpolation between tabulated quantile values (only for
            continuous variables).
        The second value specifies the parameter of option 2.
    upper_tail : array-like of shape (2,), default=None
        The first value specifies the interpolation in the upper tail of the distribution:
        1 = linear interpolation to the upper extrema;
        2 = power model interpolation to the upper extrema;
        3 = linear interpolation between tabulated quantiles;
        4 = hyperbolic model extrapolation to the upper extrema.
        The second value specifies the parameter of option 2 or 4.
    tabulated_data : DataFrame, default=None
        Simplified data if linear interpolation between tabulated values is the
        option selected for any of the three region. One legitimate choice is exactly
        the conditioning data, which is the option by default.
    tabulated_data_columns : int or str, default=(2, np.nan)
        Index or name of the variable to be used in `tabulated_data`.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    seed : int, default=42
        Seed for random number generation.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default='Class'
        Name of the simulated variable.
    output_file_name : str, default='sisim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        vario_model=(VarioModel(), VarioModel()),
        variable_type="categorical",
        class_codes=(0, 1),
        proportions=(0.75, 0.25),
        data=None,
        coord_columns=(0, 1, np.nan),
        data_column=2,
        kriging_type=0,
        use_full_ind_kriging=True,
        median_approx_class=0,
        indicator_data=None,
        indicator_coord_columns=(0, 1, np.nan),
        indicator_data_columns=(2, 3),
        use_markov_bayes_sim=False,
        calibration=(0.61, 0.54),
        secondary_data=None,
        secondary_data_column=4,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        n_realizations=1,
        prior_mean_data=None,
        max_original_data=12,
        max_prev_sim_nodes=12,
        max_soft_ind_nodes=1,
        search_separately=False,
        n_multigrids=0,
        max_data_per_octant=None,
        max_search_radii=None,
        search_ellipsoid_angles=None,
        shape_cov_lookup_table=361,
        extrema=None,
        lower_tail=None,
        middle_tail=None,
        upper_tail=None,
        tabulated_data=None,
        tabulated_data_columns=(2, np.nan),
        trimming_limits=(-1.0e21, 1.0e21),
        seed=42,
        debugging_level=0,
        variable_name="Class",
        output_file_name="sisim",
        output_dir_path=".",
        executable_path=None,
    ):

        executable = "sisim"
        if secondary_data is not None:
            executable += "_gs"
            output_file_name += "_gs"
        elif prior_mean_data is not None:
            executable += "_lm"
            output_file_name += "_lm"
        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, executable)
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        variable_types = {"categorical": 0, "continuous": 1}
        variable_type = variable_types[variable_type]
        if isinstance(vario_model, (tuple, list)) == False:
            vario_model = len(class_codes) * [vario_model]
        if tabulated_data is None:
            tabulated_data = data
            tabulated_data_columns = (data_column, np.nan)
        max_data_per_octant = (
            max_data_per_octant if max_data_per_octant is not None else 0
        )
        max_search_radii, search_ellipsoid_angles = _compute_search_ellipsoid(
            max_search_radii, search_ellipsoid_angles, vario_model
        )
        if data is None:
            vmin, vmax = min(class_codes), max(class_codes)
        else:
            vmin, vmax = data[data_column].min(), data[data_column].max()
        if extrema is None:
            extrema = (vmin, vmax)
        if lower_tail is None:
            lower_tail = (1, vmin)
        if middle_tail is None:
            middle_tail = (1, 0.0)
        if upper_tail is None:
            upper_tail = (1, vmax)

        self.spacing = spacing
        self.origin = origin
        self.n_realizations = n_realizations
        self.output_file_name = output_file_name

        self.parameters = [
            (variable_type, "\\1=continuous(cdf), 0=categorical(pdf)"),
            (len(class_codes), "\\number thresholds/categories"),
            (class_codes, "\\  thresholds / categories"),
            (np.asarray(proportions), "\\  global cdf / pdf"),
            (data, "\\file with data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data, to_list(coord_columns, 3, np.nan) + [data_column]
                ),
                "\\  columns for X,Y,Z, and variable",
            ),
        ]
        if secondary_data is not None:
            self.parameters += [
                (
                    secondary_data,
                    "\\file with gridded secondary variable",
                    output_file_name + "_ydata.dat",
                ),
                (
                    _convert_column(secondary_data, secondary_data_column),
                    "\\  column for secondary variable",
                ),
                (calibration, "\\file with calibration table", "bicalib.cal"),
            ]
        elif prior_mean_data is not None:
            self.parameters += [
                (
                    prior_mean_data,
                    "\\file with gridded indicator prior mean",
                    output_file_name + "_ydataprop.dat",
                ),
            ]
        else:
            self.parameters += [
                (indicator_data, "\\file with soft indicator input", "direct.ik"),
                (
                    _convert_column(
                        indicator_data,
                        to_list(indicator_coord_columns, 3, np.nan)
                        + list(indicator_data_columns),
                    ),
                    "\\  columns for X,Y,Z, and indicators",
                ),
                (use_markov_bayes_sim, "\\ Markov-Bayes simulation (0=no,1=yes)"),
                (calibration, "\\    calibration B(z) values"),
            ]
        self.parameters += [
            (trimming_limits, "\\trimming limits"),
            (extrema, "\\minimum and maximum data value"),
            (lower_tail, "\\  lower tail option and parameter"),
            (middle_tail, "\\  middle tail option and parameter"),
            (upper_tail, "\\  upper tail option and parameter"),
            (
                tabulated_data,
                "\\  file with tabulated values",
                output_file_name + "_tabulated_data.dat",
            ),
            (
                _convert_column(
                    tabulated_data, to_list(tabulated_data_columns, 2, np.nan)
                ),
                "\\    columns for variable, weight",
            ),
            (debugging_level, "\\debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "\\file for debugging output"),
            (output_file_name + ".out", "\\file for simulation output"),
            (n_realizations, "\\number of realizations to generate"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (seed, "\\random number seed"),
            (max_original_data, "\\maximum original data for each kriging"),
            (max_prev_sim_nodes, "\\maximum previous nodes for each kriging"),
        ]
        if prior_mean_data is None:
            self.parameters += [
                (max_soft_ind_nodes, "\\maximum soft indicator nodes for kriging")
            ]
        self.parameters += [
            (not search_separately, "\\assign data to nodes? (0=no,1=yes)"),
            (
                (True if n_multigrids > 0 else False, n_multigrids),
                "\\multiple grid search? (0=no,1=yes),num",
            ),
            (max_data_per_octant, "\\maximum per octant (0=not used)"),
            (max_search_radii, "\\maximum search radii"),
            (search_ellipsoid_angles, "\\angles for search ellipsoid"),
            (to_list(shape_cov_lookup_table), "\\size of covariance lookup table"),
            (
                (not use_full_ind_kriging, median_approx_class),
                "\\0=full IK, 1=median approx. (cutoff)",
            ),
        ]
        if prior_mean_data is None:
            self.parameters += [(kriging_type, "\\0=SK, 1=OK")]
        for i, model in enumerate(vario_model):
            self.parameters += [
                (
                    (len(model.structures), model.nugget_effect),
                    "\\" + str(i + 1) + "  nst, nugget effect",
                )
            ]
            for structure in model.structures:
                self.parameters += [
                    (
                        [
                            structure.model,
                            structure.partial_sill,
                            structure.azimuth,
                            structure.dip,
                            structure.plunge,
                        ],
                        "\\   it,cc,anql,ang2,ang3",
                    )
                ]
                self.parameters += [
                    (
                        [
                            "  ",
                            structure.range_hmax,
                            structure.range_hmin,
                            structure.range_vert,
                        ],
                        "\\   a_hmax, a_hmin, a_vert",
                    )
                ]

    def process_output_files(self):
        """
        Processes the output files before reading.
        """
        if self.output_file_name.split("_")[-1] == "lm":
            for i, line in enumerate(
                fileinput.input(self.output_file_path.with_suffix(".out"), inplace=True)
            ):
                if i == 1:
                    system = self.origin + self.spacing
                    line = (
                        line.strip("\n")
                        + " "
                        + " ".join(map(str, system))
                        + " "
                        + str(self.n_realizations)
                        + "\n"
                    )
                sys.stdout.write(line)


class GTSIM(GSLIB):
    """
    Truncated Gaussian simulation.

    Parameters
    ----------
    gaussian_realizations : array-like of shape (n_rez, n_z, n_y, n_x), default=None
        Gaussian realizations to truncate.
    proportions : array-like of shape (n_classes,), default=(0.25, 0.25, 0.5)
        Global facies proportions.
    local_proportions : array-like of shape (n_classes, n_z, n_y, n_x), default=None
        Local facies proportions.
    local_prop_columns : array-like of shape (n_classes,), default=(0, 1, 2)
        Index or name of the facies proportions in `local_proportions`.
    variable_name : string, default='Class'
        Name of the simulated variable.
    output_file_name : str, default='gtsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        gaussian_realizations,
        proportions=(0.25, 0.25, 0.5),
        local_proportions=None,
        local_prop_columns=(0, 1, 2),
        variable_name="Class",
        output_file_name="gtsim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "gtsim")
        self.variable_name = variable_name

        self.parameters = [
            (
                gaussian_realizations,
                "\\file with input Gaussian realizations",
                output_file_name + "_sgsim.out",
            ),
            (output_file_name + ".out", "\\file for output categorical realizations"),
            (
                (
                    tuple(gaussian_realizations.sizes.values())[0]
                    if isinstance(gaussian_realizations, xr.Dataset)
                    else gaussian_realizations.shape[0]
                ),
                "\\number of realizations",
            ),
            (
                (
                    tuple(gaussian_realizations.sizes.values())[1:][::-1]
                    if isinstance(gaussian_realizations, xr.Dataset)
                    else gaussian_realizations.shape[1:][::-1]
                ),
                "\\nx,ny,nz",
            ),
            (len(proportions), "\\number of categories"),
        ]
        for i, prop in enumerate(proportions):
            self.parameters += [
                (
                    (i + 1, prop),
                    "\\  cat("
                    + str(i + 1)
                    + ")  global proportion("
                    + str(i + 1)
                    + ")",
                )
            ]
        self.parameters += [
            (0 if local_proportions is None else 1, "\\proportion curves (0=no, 1=yes)")
        ]
        for i, col in enumerate(local_prop_columns):
            self.parameters += [
                (
                    local_proportions[col] if local_proportions is not None else None,
                    "\\  file with local proportion (" + str(i + 1) + ")",
                    output_file_name + "_propc" + str(i + 1) + ".dat",
                ),
                (0, "\\  column number for proportion"),
            ]


@dataclass(frozen=True)
class EllipsoidType:
    """
    Variogram model specification.

    Parameters
    ----------
    radius_1 : float, default=10.
        Radius of the ellipsoid along the direction of the azimuth angle.
    radius_2 : float, default=10.
        Radius of the ellipsoid along the direction of the dip angle.
    radius_3 : float, default=4.
        Radius of the ellipsoid along the direction of the plunge angle.
    azimuth : float, default=0.
        The first rotation angle (in degrees, clockwise) rotates the original
        Y axis (principal direction) in the horizontal plane.
    dip : float, default=0.
        The second rotation angle (in negative degrees down from horizontal)
        rotates the principal direction from the horizontal.
    plunge : float, default=0.
        The third rotation angle leaves the principal direction, defined by
        `angle_1` and `angle_2`, unchanged. The two directions orthogonal to
        that principal direction are rotated clockwise relative to the principal
        direction when looking toward the origin.
    weight : float, default=1.
         Relative proportion of the ellipsoid.
    """

    radius_1: float = 10.0
    radius_2: float = 10.0
    radius_3: float = 4.0
    azimuth: float = 0.0
    dip: float = 0.0
    plunge: float = 0.0
    weight: float = 1.0


class ELLIPSIM(GSLIB):
    """
    Boolean simulation of ellipsoids.

    Parameters
    ----------
    ellipsoid_type : EllipsoidType or array-like of shape (n_types,), default=(EllipsoidType(), EllipsoidType(2., 2., 2.), EllipsoidType(1., 1., 1.))
        Ellipsoid specification for each type.
    proportion : float, default=0.25
        Target proportion of ellipse volume.
    shape : int or array-like (x, y, z), default=(100, 100, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(1., 1., 1.)
        Cell size along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    origin : float or array-like (x, y, z), default=(0.5, 0.5, 0.)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    seed : int, default=42
        Seed for random number generation.
    variable_name : string, default='Class'
        Name of the simulated variable.
    output_file_name : str, default='ellipsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://www.statios.com/Quick/gslib.html

    References
    ----------
    Journel, A.G., Deutsch, C.V. (1992)
        GSLIB: Geostatistical Software Library and User's Guide
        http://claytonvdeutsch.com/wp-content/uploads/2019/03/GSLIB-Book-Second-Edition.pdf
    """

    def __init__(
        self,
        ellipsoid_type=(
            EllipsoidType(),
            EllipsoidType(2.0, 2.0, 2.0),
            EllipsoidType(1.0, 1.0, 1.0),
        ),
        proportion=0.25,
        shape=(100, 100, 1),
        spacing=(1.0, 1.0, 1.0),
        origin=(0.5, 0.5, 0.0),
        n_realizations=1,
        seed=42,
        variable_name="Class",
        output_file_name="ellipsim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "ellipsim")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)
        if isinstance(ellipsoid_type, (tuple, list)) == False:
            ellipsoid_type = [ellipsoid_type]

        self.parameters = [
            (output_file_name + ".out", "\\file for output realizations"),
            (n_realizations, "\\number of realizations"),
            ((shape[0], origin[0], spacing[0]), "\\nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "\\ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "\\nz,zmn,zsiz"),
            (seed, "\\random number seed"),
            (proportion, "\\target proportion (in ellipsoids)"),
        ]
        for i, etype in enumerate(ellipsoid_type):
            self.parameters += [
                (
                    (
                        etype.radius_1,
                        etype.radius_2,
                        etype.radius_3,
                        etype.azimuth,
                        etype.dip,
                        etype.plunge,
                        etype.weight,
                    ),
                    "\\radius[1,2,3],angle[1,2,3],weight",
                )
            ]
