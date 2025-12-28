"""ALLUVSIM"""

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
import numpy as np

from ..utils import to_list
from .gslib import GSLIB, _convert_column
from .fluvsim import _complete_data_coord


################################################################################
# Alluvsim


class ALLUVSIM(GSLIB):
    """
    Alluvsim, a program for event-based stochastic modeling of fluvial
    depositional systems.

    Parameters
    ----------
    channel_orientation : float or array-like (mean, stdev), default=(90., 1.)
        Primary azimuth of channel centerlines. Note current model assumes model
        proximal edge is x=0 and distal edge is x=xmax. Parameterized by mean
        and standard deviation of a Gaussian distribution.
    channel_source : float or array-like (mean, stdev), default=(500., 250.)
        Source location in Y coordinates. Source is located along the proximal
        edge of the model (e.g. x=0 for primary azimuth of 90°). Parameterized
        by mean and standard deviation of a Gaussian distribution.
    channel_thickness : float or array-like (mean, stdev), default=(4., 0.5, 0.2)
        Channel depth. Parameterized by mean and standard deviation of a Gaussian
        distribution.
    channel_width_thickness_ratio : float or array-like (mean, stdev), default=(10., 1.)
        Width-to-depth ratio. Parameterized by mean and standard deviation of a
        Gaussian distribution.
    channel_sinuosity : float or array-like (mean, stdev), default=(2., 0.2)
        Sinuosity. Based on a calibration of the Ferguson (1976) model.
        Parameterized by mean and standard deviation of a Gaussian distribution.
    levee_depth_below_top : float or array-like (mean, stdev), default=(0., 0.)
        Levee depth below top of channel fill. Parameterized by mean and standard
        deviation of a Gaussian distribution.
    levee_width : float or array-like (mean, stdev), default=(40., 5.)
        Levee width from edge of channel fill. Parameterized by mean and standard
        deviation of a Gaussian distribution.
    levee_height : float or array-like (mean, stdev), default=(0.5, 0.2)
        Levee height above top of channel fill. Parameterized by mean and standard
        deviation of a Gaussian distribution.
    levee_asymmetry_factor : float or array-like (mean, stdev), default=(0.3, 0.1)
        Factor for levee asymmetry on point bar and cut bank. For a value of 0,
        levees are symmetric and for a value of 1, levees are twice as wide on
        the cut bank side at the location of maximum curvature. Parameterized by
        mean and standard deviation of a Gaussian distribution.
    levee_distal_thinning_factor : float or array-like (mean, stdev), default=(0.3, 0.1)
        Factor for proximal to distal thinning along the centerline. For a value
        of 0, there is no thinning and for a value of 1 levee widths are doubled
        at the proximal edge and halved at the distal edge of the model.
        Parameterized by mean and standard deviation of a Gaussian distribution.
    n_splays_per_channel : int or array-like (mean, stdev), default=(0, 0)
        Number of crevasse splays along a single channel centerline.
        Parameterized by mean and standard deviation of a Gaussian distribution.
    n_lobes_per_splay : int or array-like (mean, stdev), default=(3, 2)
        Number of lobes within a single crevasse splay. Parameterized by mean
        and standard deviation of a Gaussian distribution.
    splay_source : float or array-like (mean, stdev), default=(50., 20.)
        Source location for a crevasse splay in X coordinate.
    lobe_length : float or array-like (mean, stdev), default=(200., 50.)
        Lobe length along the centerline. Parameterized by mean and standard
        deviation of a Gaussian distribution.
    lobe_max_width : float or array-like (mean, stdev), default=(30., 10.)
        Lobe maximum width. Parameterized by mean and standard deviation of a
        Gaussian distribution.
    lobe_length_max_width : float or array-like (mean, stdev), default=(100., 20.)
        Lobe length along centerline to position of maximum length.
        Parameterized by mean and standard deviation of a Gaussian distribution.
    lobe_attach_width : float or array-like (mean, stdev), default=(20., 10.)
        Lobe width at proximal edge. Parameterized by mean and standard
        deviation of a Gaussian distribution.
    lobe_height_width_ratio : float or array-like (mean, stdev), default=(0.03, 0.05)
        Lobe height-to-width ratio. Parameterized by mean and standard deviation
        of a Gaussian distribution.
    lobe_depth_width_ratio : float or array-like (mean, stdev), default=(0.02, 0.05)
        Lobe depth-to-width ratio. Parameterized by mean and standard deviation
        of a Gaussian distribution.
    fraction_abandoned_fill : float or array-like (mean, stdev), default=(0.3, 0.1)
        Fraction of abandoned channel fill assigned as fine grained. For a value
        of 0, the abandoned channel is coded as CH. For a value of 1, the entire
        abandoned channel is coded as FF(CH). For a value between 0 and 1, the
        contact between FF(CH) and CH elements are apparent.
        Parameterized by mean and standard deviation of a Gaussian distribution.
    max_migration_dist : float or array-like (mean, stdev), default=(25., 5.)
        Maximum meander migration distance. Translations calculated by the bank
        retreat model are standardized by this value. Parameterized by mean and
        standard deviation of a Gaussian distribution.
    n_time_steps: int, default=40
       Maximum number of centerlines: The algorithm terminates when this number
       of centerlines is generated or when net-to-gross is met.
    level_elevations : array-like of shape (n_levels,), default=(7., 13.1, 17.)
        List of associated levels: This is applied to define the vertical
        spacing of channels relative to Z = 0, where Zmax = (nz)(zsize).
    manning_n : float, default=0.0036
        Manning's n coefficient.
    scour_factor : float, default=10.
        Scour factor.
    gradient : float, default=0.001
        Topographic gradient.
    discharge : float, default=0.5
        Discharge.
    proba_proximal_avulsion : float, default=0.1
        Probability of avulsion proximal of the model (new centerline
        initialization).
    proba_avulsion : float, default=0.1
        Probability of avulsion within the model. One minus the sum of these
        probabilities is the probability of meander migration.
    net_to_gross : float, default=0.05
        Target net-to-gross ratio. The algorithm terminates when this
        net-to-gross ratio is exceeded.
    data : DataFrame, default=None
        Conditioning facies data. The wells are assumed to be vertical.
    coord_columns : array-like (x, y, z_top, z_bottom), default=(0, 1, 2, 3)
        Indices or names of the x, y, z_top, and z_bottom coordinates in `data`.
    data_columns : array-like (well, facies), default=(4, 5)
        Index or name of the well index and the facies in `data`.
    vert_prop_curve : array-like of shape (1,), default=None
        Relative vertical trend in channel density.
    vert_prop_column : int, default=0
        Index or name of the channel proportion in `vert_prop_curve`.
    areal_prop_map : array-like of shape (1 or 3, n_y, n_x), default=None
        Relative horizontal trend in channel density.
    areal_prop_column : int, default=0
        Index or name of the channel proportion in `areal_prop_map`.
    shape : float or array-like (x, y, z), default=(100, 100, 20)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(10., 10., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(5., 5., 0.5)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_channels_lookup : int, default=100
        Number of candidate centerlines calculated prior to model construction.
        Set this number several times larger than the maximum number of
        centerlines.
    n_nodes_fine_search : int, default=10
        Number of discretizations for spline interpolation between control nodes.
    n_nodes_corr : int, default=10
        Correlation length of the channel width Random Function (RF) (in node).
    anisotropy : float or array-like (x, y, z), default=1.
        Anisotropy ratios for nearest neighbor search along each direction.
    horiz_tol : float, default=10
        Buffer in number of control nodes. This prevents artifacts due to adjacent
        control nodes being set to honor different wells.
    vert_tol : float, default=3.
        Tolerance of net element interval thickness.
    facies_code_step : float, default=0.05
        Element code or category increment for differentiation of individual
        architectural elements.
    seed : int, default=42
        Seed for random number generation.
    variable_name : string, default='Class'
        Name of the simulated variable.
    output_file_name : str, default='alluvsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://doi.org/10.1016/j.cageo.2008.09.012

    References
    ----------
    Pyrcz, M.J., Boisvert, J.B., Deutsch, C.V. (2009)
        ALLUVSIM: A program for event-based stochastic modeling of fluvial depositional systems
        https://doi.org/10.1016/j.cageo.2008.09.012
    Zabel, F., Pyrcz, M.J. (2005)
        User’s Guide to Alluvsim Program
        https://ccg-server.engineering.ualberta.ca/CCG%20Publications/CCG%20Guidebooks/Vol%205%20User%20Guide%20to%20Alluvsim%20Program.pdf
    """

    def __init__(
        self,
        channel_orientation=(90.0, 1.0),
        channel_source=(500.0, 250.0),
        channel_thickness=(4.0, 0.5, 0.2),
        channel_width_thickness_ratio=(10.0, 1.0),
        channel_sinuosity=(2.0, 0.2),
        levee_depth_below_top=(0.0, 0.0),
        levee_width=(40.0, 5.0),
        levee_height=(0.5, 0.2),
        levee_asymmetry_factor=(0.3, 0.1),
        levee_distal_thinning_factor=(0.3, 0.1),
        n_splays_per_channel=(0, 0),
        n_lobes_per_splay=(3, 2),
        splay_source=(50.0, 20.0),
        lobe_length=(200.0, 50.0),
        lobe_max_width=(30.0, 10.0),
        lobe_length_max_width=(100.0, 20.0),
        lobe_attach_width=(20.0, 10.0),
        lobe_height_width_ratio=(0.03, 0.05),
        lobe_depth_width_ratio=(0.02, 0.05),
        fraction_abandoned_fill=(0.3, 0.1),
        max_migration_dist=(25.0, 5.0),
        n_time_steps=40,
        level_elevations=(7.0, 13.1, 17.0),
        manning_n=0.0036,
        scour_factor=10.0,
        gradient=0.001,
        discharge=0.5,
        proba_proximal_avulsion=0.1,
        proba_avulsion=0.1,
        net_to_gross=0.05,
        data=None,
        coord_columns=(0, 1, 2, 3),
        data_columns=(4, 5),
        vert_prop_curve=None,
        vert_prop_column=0,
        areal_prop_map=None,
        areal_prop_column=0,
        shape=(100, 100, 20),
        spacing=(10.0, 10.0, 1.0),
        origin=(5.0, 5.0, 0.5),
        n_channels_lookup=100,
        n_nodes_fine_search=10,
        n_nodes_corr=10,
        anisotropy=1.0,
        horiz_tol=10,
        vert_tol=3.0,
        facies_code_step=0.05,
        seed=42,
        variable_name="Class",
        output_file_name="alluvsim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "alluvsim")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)

        self.shape = shape
        self.spacing = spacing
        self.origin = origin
        self.n_realizations = 1

        self.parameters = [
            (data, "-file with well data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    data, data_columns[0:1] + coord_columns + data_columns[1:]
                ),
                "- wcol,xcol,ycol,ztcol,zbcol,fcol - removed nwell",
            ),
            (to_list(anisotropy), "-xanis,yanis,zanis"),
            ((horiz_tol, vert_tol), "- buffer, ztol"),
            (
                areal_prop_map,
                "-file with the horizontal trend",
                output_file_name + "_horizontal_trend.dat",
            ),
            (_convert_column(areal_prop_map, areal_prop_column), "- htcol"),
            (
                vert_prop_curve,
                "-file with the vertical trend",
                output_file_name + "_verttrend.dat",
            ),
            (_convert_column(vert_prop_curve, vert_prop_column), "- vtcol"),
            (
                (n_time_steps, n_time_steps, n_time_steps),
                "- ntime,max_assoc,max_assocwithin",
            ),
            (
                (len(level_elevations),) + tuple(level_elevations),
                "- nlevel, level elevations",
            ),
            (
                (net_to_gross, manning_n, scour_factor, gradient, discharge),
                "- NTGtarget,Cf,Scour_Factor,gradient,Q",
            ),
            (
                (n_channels_lookup, n_nodes_fine_search, n_nodes_corr),
                "- CHndraw,ndiscr,nCHcor",
            ),
            (
                (proba_proximal_avulsion, proba_avulsion),
                "- probAvulOutside,probAvulInside",
            ),
            (to_list(channel_orientation, 2, 0.0), "- CH element: mCHazi,stdevCHazi"),
            (to_list(channel_source, 2, 0.0), "-  mCHsource,stdevCHsource"),
            (
                to_list(channel_thickness, 3, 0.0),
                "-  mCHdepth,stdevCHdepth,stdevCHdepth2",
            ),
            (
                to_list(channel_width_thickness_ratio, 2, 0.0),
                "-  mCHwdratio,stdevCHwdratio",
            ),
            (to_list(channel_sinuosity, 2, 0.0), "-  mCHsinu,stdevCHsinu"),
            (
                to_list(levee_depth_below_top, 2, 0.0),
                "- LV Element: mLVdepth,stdevLVdepth",
            ),
            (to_list(levee_width, 2, 0.0), "-  mLVwidth,stdevLVwidth"),
            (to_list(levee_height, 2, 0.0), "-  mLVheight,stdevLVheight"),
            (to_list(levee_asymmetry_factor, 2, 0.0), "-  mLVasym,stdevLVasym"),
            (to_list(levee_distal_thinning_factor, 2, 0.0), "-  mLVthin,stdevLVthin"),
            (to_list(n_splays_per_channel, 2, 0), "- CS Element: mCSnum,stdevCSnum"),
            (to_list(n_lobes_per_splay, 2, 0), "-  mCSnumlobe,stdevCSnumlobe"),
            (to_list(splay_source, 2, 0.0), "-  mCSsource,stdevCSsource"),
            (to_list(lobe_length, 2, 0.0), "-  mCSLOLL,stdevCSLOLL"),
            (to_list(lobe_max_width, 2, 0.0), "-  mCSLOWW,stdevCSLOWW"),
            (to_list(lobe_length_max_width, 2, 0.0), "-  mCSLOl,stdevCSLOl"),
            (to_list(lobe_attach_width, 2, 0.0), "-  mCSLOw,stdevCSLOw"),
            (
                to_list(lobe_height_width_ratio, 2, 0.0),
                "-  mCSLO_hwratio,stdevCSLO_hwratio",
            ),
            (
                to_list(lobe_depth_width_ratio, 2, 0.0),
                "-  mCSLO_dwratio,stdevCSLO_dwratio",
            ),
            (
                to_list(fraction_abandoned_fill, 2, 0.0),
                "- FFCH Element: mFFCHprop,stdevFFCHprop",
            ),
            (to_list(max_migration_dist, 2, 0.0), "- mdistMigrate,stdevdistMigrate"),
            ((shape[0], origin[0], spacing[0]), "-nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "-ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "-nz,zmn,zsiz"),
            ((seed, facies_code_step), "-random number seed, color_incr"),
            (output_file_name + ".out", "-file for output facies file"),
            (output_file_name + ".geo", "-file for output updated streamlines"),
            (output_file_name + ".rep", "-well conditioning report"),
            (output_file_name + ".geoup", "-file for output updated streamlines"),
        ]

        self.output_file_paths += [self.output_file_path.with_suffix(".geo")]
        self.output_file_paths += [self.output_file_path.with_suffix(".rep")]
        # self.output_file_paths += [self.output_file_path.with_suffix('.geoup')]

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


################################################################################
# MAPSpp


class MAPSpp(GSLIB):
    """
    MAPSpp, a program to post-process object-based model to reproduce well data.

    Parameters
    ----------
    data : DataFrame
        Conditioning data.
    realization : array-like of shape (n_vars, n_z, n_y, n_x)
        Realization to post-process.
    coord_columns : int, str, or array-like (x, y, z), default=(1, 2, 3)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_column : int or str, default=4
        Index or name of the variable that was simulated in `data`.
    realization_column : int or str, default=0
        Index or name of the variable that was simulated in `realization`.
    shape : float or array-like (x, y, z), default=(100, 100, 20)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=(10., 10., 1.)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=(5., 5., 0.5)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    max_update_distance : array-like of shape (3,), default=(9, 9, 2)
        Maximum distance to update. This parameter is robust and will rarely need
        to be changed.
    weighting_exponent : float, default=2.
        Weighting exponent. This parameter is robust and will rarely need to be
        changed.
    maps_window : array-like of shape (3,), default=(9, 9, 2)
        MAPS window, in cells. This parameter is robust and will rarely need to
        be changed.
    trimming_limits : array-like (min, max), default=(-1.e21, 1.e21)
        Limits below and above which the values of the variable are ignored.
    seed : int, default=42
        Seed for random number generation.
    output_file_name : str, default='alluvsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://doi.org/10.1016/j.cageo.2008.09.012

    References
    ----------
    Pyrcz, M.J., Boisvert, J.B., Deutsch, C.V. (2009)
        ALLUVSIM: A program for event-based stochastic modeling of fluvial depositional systems
        https://doi.org/10.1016/j.cageo.2008.09.012
    """

    def __init__(
        self,
        data,
        realization,
        coord_columns=(1, 2, 3),
        data_column=4,
        realization_column=0,
        shape=(100, 100, 20),
        spacing=(10.0, 10.0, 1.0),
        origin=(5.0, 5.0, 0.5),
        max_update_distance=(9, 9, 2),
        weighting_exponent=2.0,
        maps_window=(2, 2, 2),
        trimming_limits=(-1.0e21, 1.0e21),
        seed=42,
        output_file_name="mapspp",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "MAPSpp")

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)

        coord_columns = to_list(coord_columns, 3, np.nan)
        data_column = [data_column]
        data = _complete_data_coord(data, coord_columns, data_column, origin)

        self.shape = shape
        self.spacing = spacing
        self.origin = origin

        self.parameters = [
            (data, "-file with well conditioning data", output_file_name + "_data.dat"),
            (
                _convert_column(data, to_list(coord_columns, 3, np.nan) + data_column),
                "-  columns for X, Y, Z, facies",
            ),
            (trimming_limits, "-  trimming limits"),
            (
                realization,
                "-file with inital image",
                output_file_name + "_realization.dat",
            ),
            (
                _convert_column(realization, realization_column),
                "-  column for categorical variable",
            ),
            (output_file_name + ".out", "-file for simulation output"),
            (1, "-number of realizations"),
            ((shape[0], origin[0], spacing[0]), "-nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "-ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "-nz,zmn,zsiz"),
            (seed, "-random number seed      |  These should be frozen and"),
            (
                max_update_distance,
                "-max distance to update  |  not made user adjustable",
            ),
            (weighting_exponent, "-weighting exponent      |"),
            (maps_window, "-MAPS window (cells)     |"),
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
                    + str(1)
                    + "\n"
                )
            sys.stdout.write(line)
