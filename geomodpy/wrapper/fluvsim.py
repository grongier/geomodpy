"""FLUVSIM"""

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
import fileinput
from copy import copy
import numpy as np
import pandas as pd

from ..utils import to_list
from .gslib import GSLIB, _convert_column


################################################################################
# Data management


def _complete_data_coord(
    data, coord_columns, data_columns, shape, spacing, origin, is_fluvsim=False
):
    """
    Completes the (x,y,z) coordinates of a data set if any are missing.
    """
    _data = copy(data)
    if _data is not None:
        if len(data_columns) > 1 and is_fluvsim == True:
            if isinstance(_data, pd.DataFrame):
                _data.loc[_data[data_columns[1]] == 2, data_columns[1]] = 1
            else:
                _data[_data[data_columns[1]] == 2, data_columns[1]] = 1
        if (
            isinstance(coord_columns[0], str) == False
            and np.isnan(coord_columns[0]) == True
        ):
            if isinstance(_data, pd.DataFrame):
                _data["GSLIBX"] = origin[0]
                coord_columns[0] = "GSLIBX"
            else:
                _data = np.c_[_data, np.zeros(len(_data))]
                coord_columns[0] = _data.shape[1] - 1
        if (
            isinstance(coord_columns[1], str) == False
            and np.isnan(coord_columns[1]) == True
        ):
            if isinstance(_data, pd.DataFrame):
                _data["GSLIBY"] = origin[1]
                coord_columns[1] = "GSLIBY"
            else:
                _data = np.c_[_data, np.zeros(len(_data))]
                coord_columns[1] = _data.shape[1] - 1
        if (
            isinstance(coord_columns[2], str) == True
            or np.isnan(coord_columns[2]) == False
        ):
            _data[coord_columns[2]] -= origin[2] - spacing[2] / 2.0
            _data[coord_columns[2]] /= spacing[2] * shape[2]
        else:
            if isinstance(_data, pd.DataFrame):
                _data["GSLIBZ"] = origin[2]
                coord_columns[2] = "GSLIBZ"
            else:
                _data = np.c_[_data, np.zeros(len(_data))]
                coord_columns[2] = _data.shape[1] - 1

    return _data


################################################################################
# Fluvsim


class FLUVSIM(GSLIB):
    """
    Fluvsim, a program for object-based stochastic modeling of fluvial
    depositional systems.

    Parameters
    ----------
    channel_orientation : float or array-like (min, mode, max), default=(-30., 0., 30.)
        Channel orientation (angle in degrees measured clockwise from North/Y-axis).
        The angles can be negative or positive (e.g., -10., 0., 10. would orient
        the channels parallel to the Y-axis with a 10 degree deviation in the
        channel direction).
    channel_amplitude : float or array-like (min, mode, max), default=200.
        Average departure from channel center line (horizontal distance units)
        representing the half width of the meander.
    channel_wavelength : float or array-like (min, mode, max), default=800.
        Horizontal correlation length for sinusoidal departure from channel center
        line (horizontal distance units).
    channel_thickness : float or array-like (min, mode, max), default=(1., 3., 5.)
        Channel thickness (vertical distance units).
    channel_thickness_undulation : float or array-like (min, mode, max), default=1.
        Average magnitude of channel thickness undulation (fraction relative to
        channel thickness). The channels will have a constant thickness if this
        parameter is set to either 0.0 or 1.0. The channels will vary between 0.8
        and 1.2 times the average thickness when this parameter is set to 0.2.
    channel_thickness_undulation_wavelength : float or array-like (min, mode, max), default=(250., 400., 450.)
        Horizontal correlation length for channel thickness undulation (horizontal
        distance units). The channel thickness undulation follows a Gaussian
        histogram; the correlation length is the distance range of correlation
        along the axis of the channel.
    channel_width_thickness_ratio : float or array-like (min, mode, max), default=(150., 200., 250.)
        Channel width to thickness ratio (horizontal to vertical distance units).
    channel_width_undulation : float or array-like (min, mode, max), default=1.
        Channel width undulation (fraction relative to channel width). See explanation
        for thickness undulation (line 34) for more details.
    channel_width_undulation_wavelength : float or array-like (min, mode, max), default=250.
         Horizontal correlation length for channel width undulation (horizontal
         distance units).
    levee_width : float or array-like (min, mode, max), default=(160., 240., 320.)
        Levee width (horizontal distance units).
    levee_height : float or array-like (min, mode, max), default=0.1
        Levee height (relative to thickness of channel).
    levee_depth_below_top : float or array-like (min, mode, max), default=(0.2, 0.3, 0.4)
        Levee depth below channel top (relative to thickness of channel).
    crevasse_attach_length : float or array-like (min, mode, max), default=80.
        Crevasse attachment length (horizontal distance units).
    crevasse_relative_thickness : float or array-like (min, mode, max), default=(0.25, 0.5, 0.75)
        Crevasse thickness next to channel (relative to channel thickness).
    crevasse_diameter : float or array-like (min, mode, max), default=500.
        Areal size of crevasse (diameter in horizontal distance units).
    max_channels : int, default=150
        Maximum number of channel objects.
    channel_prop : float, default=0.3
        Global proportions of channel sands.
    levee_prop : float, default=None
        Global proportions of levee sands.
    crevasse_prop : float, default=None
        Global proportions of crevasse sands; the shale proportion is 1.0 minus
        the sum of sand proportions.
    data : DataFrame, default=None
        Conditioning facies data.
    coord_columns : array-like (x, y, z), default=(0, 1, 2)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_columns : int or str, default=(3, 4)
        Index or name of the well number and the variable to be simulated in `data`.
    vert_prop_curve : array-like of shape (1 or 3, n_z), default=None
        Vertical facies proportion curves.
    is_vert_prop_net_to_gross : float, default=True
        If True, the vertical facies proportion curves correspond to the net-to-gross
        ratio (lumped proportion of all sand facies); otherwise, they correspond
        to each facies proportion.
    vert_prop_columns : array-like (channel, levee, crevasse), default=(0, np.nan, np.nan)
        Index or name of the net-to-gross or facies proportions in `vert_prop_curve`.
    areal_prop_map : array-like of shape (1 or 3, n_y, n_x), default=None
        Areal facies proportion maps.
    is_areal_prop_net_to_gross : float, default=False
        If True, the areal facies proportion curves correspond to the net-to-gross
        ratio (lumped proportion of all sand facies); otherwise, they correspond
        to each facies proportion.
    areal_prop_columns : array-like (channel, levee, crevasse), default=(0, 1, 2)
        Index or name of the net-to-gross or facies proportions in `areal_prop_map`.
    global_weight : float, default=1.
        Multiplicative weight for the global facies proportion. The program
        automatically determines weights for each component to allow convergence
        of all constraints; however, these multiplicative weights modify the
        automatically determined weights to place more importance on selected
        components.
    vert_weight : float, default=1.
        Multiplicative weight for the vertical facies proportion curve. The program
        automatically determines weights for each component to allow convergence
        of all constraints; however, these multiplicative weights modify the
        automatically determined weights to place more importance on selected
        components.
    areal_weight : float, default=1.
        Multiplicative weight for the areal facies proportion map. The program
        automatically determines weights for each component to allow convergence
        of all constraints; however, these multiplicative weights modify the
        automatically determined weights to place more importance on selected
        components.
    data_weight : float, default=1.
        Multiplicative weight for the conditioning facies data. The program
        automatically determines weights for each component to allow convergence
        of all constraints; however, these multiplicative weights modify the
        automatically determined weights to place more importance on selected
        components.
    shape : int or array-like (x, y, z), default=(100, 100, 50)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(40., 40., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(20., 20., 0.5)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    max_iter : int, default=100
        Maximum number of iterations.
    n_iter_no_change : int, default=10
        Maximum number of iterations without a change to the objective function.
    tol : float, default=0.05
        Minimum objective function (stopping criteria).
    init_temp : float, default=0.
        Initial temperature of the annealing schedule.
    reduct_factor : float, default=0.1
        Reduction factor of the annealing schedule.
    max_perturb : int, default=3
        Maximum number of perturbations at any one given temperature of the
        annealing schedule.
    n_perturb : int, default=1
        Target number of acceptable perturbations at a given temperature of the
        annealing schedule.
    n_max_perturb : int, default=8
        Stopping number (maximum number of times that `max_perturb` is reached).
    proba_channels_on_off : float, default=1.
        Probability of turning a channel entity on and one off as a perturbation
        mechanism.
    proba_channel_on : float, default=0.1
        Probability of turning a channel entity on as a perturbation mechanism.
    proba_channel_off : float, default=0.1
        Probability of turning a channel entity off as a perturbation mechanism.
    proba_fix_well : float, default=1.
        Probability of picking a well interval at random attempt to fix conditioning.
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
    output_file_name : str, default='fluvsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://doi.org/10.1016/S0098-3004(01)00075-9

    References
    ----------
    Deutsch, C.V., Tran, T.T. (2002)
        FLUVSIM: a program for object-based stochastic modeling of fluvial depositional systems
        https://doi.org/10.1016/S0098-3004(01)00075-9
    """

    def __init__(
        self,
        channel_orientation=(-30.0, 0.0, 30.0),
        channel_amplitude=200.0,
        channel_wavelength=800.0,
        channel_thickness=(1.0, 3.0, 5.0),
        channel_thickness_undulation=1.0,
        channel_thickness_undulation_wavelength=(250.0, 400.0, 450.0),
        channel_width_thickness_ratio=(150.0, 200.0, 250.0),
        channel_width_undulation=1.0,
        channel_width_undulation_wavelength=250.0,
        levee_width=(160.0, 240.0, 320.0),
        levee_height=0.1,
        levee_depth_below_top=(0.2, 0.3, 0.4),
        crevasse_attach_length=80.0,
        crevasse_relative_thickness=(0.25, 0.5, 0.75),
        crevasse_diameter=500.0,
        max_channels=150,
        channel_prop=0.3,
        levee_prop=None,
        crevasse_prop=None,
        data=None,
        coord_columns=(0, 1, 2),
        data_columns=(3, 4),
        vert_prop_curve=None,
        is_vert_prop_net_to_gross=True,
        vert_prop_columns=(0, np.nan, np.nan),
        areal_prop_map=None,
        is_areal_prop_net_to_gross=False,
        areal_prop_columns=(0, 1, 2),
        global_weight=1.0,
        vert_weight=1.0,
        areal_weight=1.0,
        data_weight=1.0,
        shape=(100, 100, 50),
        spacing=(40.0, 40.0, 1.0),
        origin=(20.0, 20.0, 0.5),
        n_realizations=1,
        max_iter=100,
        n_iter_no_change=10,
        tol=0.05,
        init_temp=0.0,
        reduct_factor=0.1,
        max_perturb=3,
        n_perturb=1,
        n_max_perturb=8,
        proba_channels_on_off=1.0,
        proba_channel_on=0.1,
        proba_channel_off=0.1,
        proba_fix_well=1.0,
        trimming_limits=(-1.0e21, 1.0e21),
        seed=42,
        debugging_level=0,
        variable_name="Class",
        output_file_name="fluvsim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "fluvsim")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.5)
        are_weights = [
            False if global_weight is None else True,
            False if vert_weight is None or vert_prop_curve is None else True,
            False if areal_weight is None or areal_prop_map is None else True,
            False if data_weight is None or data is None else True,
        ]
        weights = [
            0.0 if global_weight is None else global_weight,
            0.0 if vert_weight is None else vert_weight,
            0.0 if areal_weight is None else areal_weight,
            0.0 if data_weight is None else data_weight,
        ]
        are_facies = [
            True,
            False if levee_prop is None else True,
            False if crevasse_prop is None else True,
        ]
        proportions = [
            0.0 if channel_prop is None else channel_prop,
            0.0 if levee_prop is None else levee_prop,
            0.0 if crevasse_prop is None else crevasse_prop,
        ]
        if isinstance(vert_prop_columns, (tuple, list)) == False:
            vert_prop_columns = 3 * (vert_prop_columns,)
        elif len(vert_prop_columns) == 1:
            vert_prop_columns = 3 * vert_prop_columns
        if isinstance(areal_prop_columns, (tuple, list)) == False:
            areal_prop_columns = 3 * (areal_prop_columns,)
        elif len(areal_prop_columns) == 1:
            areal_prop_columns = 3 * areal_prop_columns

        coord_columns = to_list(coord_columns, 3, np.nan)
        data_columns = list(data_columns)
        data = _complete_data_coord(
            data, coord_columns, data_columns, shape, spacing, origin, is_fluvsim=True
        )

        self.shape = shape
        self.spacing = spacing
        self.origin = origin
        self.n_realizations = n_realizations

        self.parameters = [
            (data, "-file with well conditioning data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data, data_columns[0:1] + coord_columns + data_columns[1:]
                ),
                "-  columns for X, Y, Z, well #, facies",
            ),
            (trimming_limits, "-  trimming limits"),
            (debugging_level, "-debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "-file for debugging output"),
            (output_file_name + ".geo", "-file for geometric specification"),
            (output_file_name + ".out", "-file for simulation output"),
            (output_file_name + ".vp", "-file for vertical prop curve output"),
            (output_file_name + ".ap", "-file for areal prop map output"),
            (output_file_name + ".wd", "-file for well data output"),
            (n_realizations, "-number of realizations to generate"),
            (
                (shape[0], origin[0], spacing[0]),
                "-nx,xmn,xsiz - geological coordinates",
            ),
            (
                (shape[1], origin[1], spacing[1]),
                "-ny,ymn,ysiz - geological coordinates",
            ),
            (
                (shape[2], spacing[2] * shape[2]),
                "-nz, average thickness in physical units",
            ),
            (seed, "-random number seed"),
            (are_weights, "-1=on,0=off: global, vert, areal, wells"),
            (weights, "-weighting : global, vert, areal, wells"),
            (
                (max_iter, n_iter_no_change, tol),
                "-maximum iter, max no change, min. obj.",
            ),
            (
                (init_temp, reduct_factor, max_perturb, n_perturb, n_max_perturb),
                "-annealing schedule: t0,redfac,ka,k,num",
            ),
            (
                (
                    proba_channels_on_off,
                    proba_channel_on,
                    proba_channel_off,
                    proba_fix_well,
                ),
                "-Pert prob: 1on+1off, 1on, 1off, fix well",
            ),
            (are_facies, "-Facies(on): channel, levee, crevasse"),
            (proportions, "-Proportion: channel, levee, crevasse"),
            (
                vert_prop_curve,
                "-  vertical proportion curves",
                output_file_name + "_pcurve.dat",
            ),
            (not is_vert_prop_net_to_gross, "-     0=net-to-gross, 1=all facies"),
            (
                _convert_column(vert_prop_curve, vert_prop_columns),
                "-     column numbers",
            ),
            (
                areal_prop_map,
                "-  areal proportion map",
                output_file_name + "_arealprop.dat",
            ),
            (not is_areal_prop_net_to_gross, "-     0=net-to-gross, 1=all facies"),
            (
                _convert_column(areal_prop_map, areal_prop_columns),
                "-     column numbers",
            ),
            (max_channels, "-maximum number of channels"),
            (to_list(channel_orientation), "-channel:  orientation (degrees)"),
            (to_list(channel_amplitude), "-channel:  sinuosity: average departure"),
            (to_list(channel_wavelength), "-channel:  sinuosity: length scale"),
            (to_list(channel_thickness), "-channel:  thickness"),
            (to_list(channel_thickness_undulation), "-channel:  thickness undulation"),
            (
                to_list(channel_thickness_undulation_wavelength),
                "-channel:  thickness undul. length scale",
            ),
            (
                to_list(channel_width_thickness_ratio),
                "-channel:  width/thickness ratio",
            ),
            (to_list(channel_width_undulation), "-channel:  width: undulation"),
            (
                to_list(channel_width_undulation_wavelength),
                "-channel:  width: undulation length scale",
            ),
            (to_list(levee_width), "-levee:    average width"),
            (to_list(levee_height), "-levee:    average height"),
            (to_list(levee_depth_below_top), "-levee:    depth below top"),
            (to_list(crevasse_attach_length), "-crevasse: attachment length"),
            (
                to_list(crevasse_relative_thickness),
                "-crevasse: relative thickness by channel",
            ),
            (to_list(crevasse_diameter), "-crevasse: areal size (diameter)"),
        ]

        # self.output_file_paths += [self.output_file_path.with_suffix('.geo')]
        if vert_prop_curve is not None and vert_weight is not None:
            self.output_file_paths += [self.output_file_path.with_suffix(".vp")]
        if areal_prop_map is not None and areal_weight is not None:
            self.output_file_paths += [self.output_file_path.with_suffix(".ap")]
        # if data is not None and data_weight is not None:
        #     self.output_file_paths += [self.output_file_path.with_suffix('.wd')]

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
