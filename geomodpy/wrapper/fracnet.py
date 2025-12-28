"""FRACNET"""

# MIT License

# Copyright (c) 2023 Guillaume Rongier

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
import fileinput
import numpy as np

from ..utils import to_list
from .gslib import GSLIB, _convert_column


################################################################################
# Fracnet


class FRACNET(GSLIB):
    """
    Fracnet, a program for the stochastic simulation of fractures in layered systems.

    Parameters
    ----------
    fracture_strike : float, array-like, or DataFrame
        Histogram of fracture strikes. If a float, then there is a single strike
        without weight.
    fracture_strike_columns : int, str, or array-like (strike, weight), default=0
        Index or name of the strike variable and associated weights (if available)
        in `fracture_strike`.
    fracture_dip : float, array-like, or DataFrame
        Histogram of fracture dips. If a float, then there is a single dip without
        weight.
    fracture_dip_columns : int, str, or array-like (dip, weight), default=0
        Index or name of the dip variable and associated weights (if available)
        in `fracture_dip`.
    fracture_spacing : float, array-like, or DataFrame, default=None
        Histogram of fracture spacing. If a float, then there is a single spacing
        without weight.
    fracture_spacing_columns : int, str, or array-like (spacing, weight), default=0
        Index or name of the spacing variable and associated weights in (if available)
        `fracture_spacing`.
    min_fracture_spacing : int, default=None
        Minimum fracture spacing in cells (0 or None = not used).
    fracture_length : float, array-like, or DataFrame, default=None
        Histogram of fracture lengths. If a float, then there is a single length
        without weight.
    fracture_length_columns : array-like (length, weight), default=0
        Index or name of the length variable and associated weights in `length_spacing`.
    max_fracture_length : float, default=None
        Maximum fracture length (0 or None = not used; 1 = maximum length).
    density_field : int or array-like of shape (n_d_vars, n_d_z, n_d_y, n_d_x)
        Exhaustive "scanline" density map of the domain to be simulated. It
        contains the number of fractures per grid block. If an integer, then
        the whole domain is a single grid block, and this is the uniform density
        for the domain.
    density_field_column : int or str, default=0
        Index or name of the density variable in `density_field`.
    orientation_field : array-like of shape (n_o_vars, n_d_z, n_d_y, n_d_x), default=None
        Orientation field information (same dimensions as the density grid).
    orientation_field_column : int or str, default=0
        Index or name of the orientation variable in `orientation_field`.
    data : DataFrame, default=None
        Conditioning data containing the data location coordinates, the type
        (fracture = 1, matrix = 0), the strike and the dip.
    primary_fractures : array-like of shape (n_z, n_y, n_x), default=None
        Previously simulated fracture set to use as primary set.
    use_hierarchical_model : boolean, default=True
        If true, consider a hierarchical model.
    fracture_stopping_proba : float, default=1.
        Probability of stopping on a previously simulated fracture.
    use_max_grow : boolean, default=False
        If true, grow fractures to their maximum extent (it overrules fracture density).
    shale_laminaes : array-like of shape (n_d_z, n_d_y, n_d_x), default=None
        Previously simulated shale distribution (same dimensions as the density grid).
        1 means that the top boundary of the cell is made of shale, 0 means that
        there is no shale.
    shale_stopping_proba : float, default=1.
        Probability of stopping on a previously simulated fracture.
    shape : int or array-like (x, y, z), default=40
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=1.
        Cell size along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    seed : int, default=42
        Seed for random number generation.
    variable_name : string, default='Class'
        Name of the simulated variable.
    output_file_name : str, default='fracnet'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://doi.org/10.1016/S0098-3004(98)00071-5

    References
    ----------
    Gringarten, E. (1998)
        FRACNET: Stochastic simulation of fractures in layered systems
        https://doi.org/10.1016/S0098-3004(98)00071-5
    """

    def __init__(
        self,
        fracture_strike=0.0,
        fracture_strike_columns=0,
        fracture_dip=0.0,
        fracture_dip_columns=0,
        fracture_spacing=None,
        fracture_spacing_columns=0,
        min_fracture_spacing=0,
        fracture_length=None,
        fracture_length_columns=0,
        max_fracture_length=None,
        density_field=10,
        density_field_column=0,
        orientation_field=None,
        orientation_field_column=0,
        data=None,
        primary_fractures=None,
        use_hierarchical_model=True,
        fracture_stopping_proba=1.0,
        use_max_grow=False,
        shale_laminaes=None,
        shale_stopping_proba=1.0,
        shape=40,
        spacing=1.0,
        n_realizations=1,
        seed=42,
        variable_name="Class",
        output_file_name="fracnet",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "fracnet")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)

        if isinstance(density_field, (int, float)):
            density_field = np.full((1, 1, 1, 1), density_field)
            density_field_column = 0
        if isinstance(fracture_strike, (int, float)):
            fracture_strike = np.array([[fracture_strike]])
            fracture_strike_columns = (0, np.nan)
        if isinstance(fracture_dip, (int, float)):
            fracture_dip = np.array([[fracture_dip]])
            fracture_dip_columns = (0, np.nan)
        if isinstance(fracture_spacing, (int, float)):
            fracture_spacing = np.array([[fracture_spacing]])
            fracture_spacing_columns = (0, np.nan)
        if isinstance(fracture_length, (int, float)):
            fracture_length = np.array([[fracture_length]])
            fracture_length_columns = (0, np.nan)

        self.shape = shape
        self.spacing = spacing
        self.origin = [s / 2.0 for s in spacing]
        self.n_realizations = n_realizations

        self.parameters = [
            (False if data is None else True, "\\condition to data 1=yes, 0=no"),
            (
                data,
                "\\data file (x,y,z,frac/mat,strike,dip)",
                output_file_name + "_data.dat",
            ),
            (density_field, "\\density data file", output_file_name + "_density.dat"),
            (density_field[density_field_column].shape[::-1], "\\nxd, nyd, nzd"),
            (
                _convert_column(density_field, density_field_column),
                "\\column for density",
            ),
            (fracture_strike, "\\strikes data file", output_file_name + "_strike.dat"),
            (
                _convert_column(
                    fracture_strike, to_list(fracture_strike_columns, 2, np.nan)
                ),
                "\\columns: strike, weight",
            ),
            (fracture_dip, "\\dips data file", output_file_name + "_dip.dat"),
            (
                _convert_column(fracture_dip, to_list(fracture_dip_columns, 2, np.nan)),
                "\\columns: dip, weight",
            ),
            (
                False if fracture_spacing is None else True,
                "\\condition to spacing cdf 1=yes, 0=no",
            ),
            (
                fracture_spacing,
                "\\spacing data file",
                output_file_name + "_spacing.dat",
            ),
            (
                _convert_column(
                    fracture_spacing, to_list(fracture_spacing_columns, 2, np.nan)
                ),
                "\\columns: spacing, weight",
            ),
            (
                False if min_fracture_spacing is None else min_fracture_spacing,
                "\\min. spacing between frac. (pixels)",
            ),
            (
                (
                    False if max_fracture_length is None else True,
                    0 if max_fracture_length is None else max_fracture_length,
                ),
                "\\max length spec.: 1=yes, 0=no lmax",
            ),
            (
                False if fracture_length is None else True,
                "\\condition to length cdf 1=yes, 0=no",
            ),
            (fracture_length, "\\length data file", output_file_name + "_length.dat"),
            (
                _convert_column(
                    fracture_length, to_list(fracture_length_columns, 2, np.nan)
                ),
                "\\columns: length, weight",
            ),
            (
                False if primary_fractures is None else True,
                "\\existing primary set 1=yes, 0=no",
            ),
            (use_hierarchical_model, "\\hierarchical model 1=yes, 0=no"),
            (
                primary_fractures,
                "\\primary set data file",
                output_file_name + "_primary_set.dat",
            ),
            (
                (
                    False if fracture_stopping_proba is None else True,
                    1.0 if fracture_stopping_proba is None else fracture_stopping_proba,
                ),
                "\\prob.of abuting: 1=yes, 0=no; prob",
            ),
            (use_max_grow, "\\grow 2nd set to max extend"),
            (
                False if shale_laminaes is None else True,
                "\\existing shales 1=yes, 0=no",
            ),
            (shale_laminaes, "\\shale data file", output_file_name + "_shales.dat"),
            (
                (
                    False if shale_stopping_proba is None else True,
                    1.0 if shale_stopping_proba is None else shale_stopping_proba,
                ),
                "\\prob.of stopping: 1=yes, 0=no; prob",
            ),
            (
                False if orientation_field is None else True,
                "\\orientation field 1=yes, 0=no",
            ),
            (
                orientation_field,
                "\\orientation data file",
                output_file_name + "_orientation.dat",
            ),
            (
                _convert_column(orientation_field, orientation_field_column),
                "\\column for orientation",
            ),
            (output_file_name + ".out", "\\output file for simulation"),
            (seed, "\\random number seed"),
            (n_realizations, "\\number of simulations"),
            ((shape[0], spacing[0]), "\\nx,xsiz"),
            ((shape[1], spacing[1]), "\\ny,ysiz"),
            ((shape[2], spacing[2]), "\\nz,zsiz"),
        ]

    def process_output_files(self):
        """
        Processes the output files before reading.
        """
        for i, line in enumerate(
            fileinput.input(self.output_file_path.with_suffix(".out"), inplace=True)
        ):
            if i == 1:
                system = self.shape + self.origin + self.spacing
                line = (
                    line.strip("\n")
                    + " "
                    + " ".join(map(str, system))
                    + " "
                    + str(self.n_realizations)
                    + "\n"
                )
            sys.stdout.write(line)
