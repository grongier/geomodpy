"""SNESIM"""

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


import re
import fileinput
from copy import copy, deepcopy
import numpy as np

from ..utils import to_list
from .gslib import GSLIB, _convert_column


################################################################################
# Fluvsim


class SNESIM(GSLIB):
    """
    Snesim, a program for multiple-point simulation.

    Parameters
    ----------
    training_image : array-like of shape (n_vars, n_z, n_y, n_x)
        Training image.
    training_column : int or str, default=0
        Index or name of the variable to be simulated in `training_image`.
    class_codes : array-like of shape (n_classes,), default=None
        Category codes. If None, they are inferred from the training image.
    proportions : array-like of shape (n_classes,), default=None
        Target global pdf values. If None, they are inferred from the training image.
    servosystem_param : float, default=0.
        Parameter controlling the servosystem correction to “bend” the simulated
        marginal probability towards the target proportion.
    data : DataFrame, default=None
        Conditioning data.
    coord_columns : int, str, or array-like (x, y, z), default=(0, 1, np.nan)
        Indices or names of the x, y, and z coordinates in `data`. One or two of
        the coordinates can be set to np.nan or pd.NA, which indicates that the
        simulation is 2D or 1D. The last and second-to-last coordinates can also
        be omitted instead of set to NaN.
    data_column : int or str, default=3
        Index or name of the variable to be simulated in `data`.
    vert_prop_curve : array-like of shape (n_classes, n_z), default=None
        Vertical category proportion curves.
    secondary_data : array-like of shape (n_vars, n_z, n_y, n_z), default=None
        Secondary data.
    tau : array-like of shape (2,), default=None
        Two weighting factors to combine P(A|B) and P(A|C). If None, those factors
        are determined automatically.
    rotation_affinity : array-like of shape (n_vars, n_z, n_y, n_z), default=None
        Local rotation and affinity classes.
    affinity_factors : array-like n_factors*(x, y, z), default=((1., 1., 1.), (1., 0.5, 1.), (1., 2., 1.))
        Affinity factors for each affinity class.
    shape : int or array-like (x, y, z), default=(250, 250, 1)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z), default=1.
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z), default=0.5
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=1
        Number of realizations to simulate.
    n_multigrids : int, default=5
        The number of multiple grid refinements to consider.
    max_data : int, default=26
        Maximum number of primary data that should be used to simulate a grid node.
    min_replicates : int, default=10
        Minimum number of replicates for a data event to help get reliable statistics
        when using Bayesian-updating or servosystem to correct the marginal probability.
    max_search_radii : float or array-like (x, y, z), default=None
        The search radii in the maximum horizontal direction, minimum horizontal
        direction, and vertical direction (see angles below). If None, it is
        equal to a third of the grid dimensions.
    search_ellipsoid_angles : float or array-like (x, y, z), default=0.
        The angle parameters that describe the orientation of the search ellipsoid.
    seed : int, default=42
        Seed for random number generation.
    debugging_level : int, default=0
        Debugging level between 0 and 3. The larger the debugging level, the more
        information written out in an output file. By default, this file is
        automatically deleted (see `clean_files` in `run`).
    variable_name : string, default=('Class', 'Number_conditioning_data')
        Name of the simulated variable and of the number of conditioning data retained.
    output_file_name : str, default='fluvsim'
        Base name for the output files.
    output_dir_path : str, default='.'
        Path to the directory for the output files.
    executable_path : str, default=None
        Path to an executable for a GSLIB program. By default, looks for it in
        the `bin` directory of geomodpy.

    Source
    ------
    https://github.com/SCRFpublic/snesim-standalone

    References
    ----------
    Strebelle, S. (2002)
        Conditional Simulation of Complex Geological Structures Using Multiple-Point Statistics
        https://doi.org/10.1023/A:1014009426274
    """

    def __init__(
        self,
        training_image,
        training_column=0,
        class_codes=None,
        proportions=None,
        servosystem_param=0.0,
        data=None,
        coord_columns=(0, 1, np.nan),
        data_column=3,
        vert_prop_curve=None,
        secondary_data=None,
        tau=None,
        rotation_affinity=None,
        affinity_factors=((1.0, 1.0, 1.0), (1.0, 0.5, 1.0), (1.0, 2.0, 1.0)),
        shape=(250, 250, 1),
        spacing=1.0,
        origin=0.5,
        n_realizations=1,
        n_multigrids=5,
        max_data=26,
        min_replicates=10,
        max_search_radii=None,
        search_ellipsoid_angles=0.0,
        seed=42,
        debugging_level=0,
        variable_name=("Class", "Number_conditioning_data"),
        output_file_name="snesim",
        output_dir_path=".",
        executable_path=None,
    ):

        GSLIB.__init__(self, output_file_name, output_dir_path, executable_path)
        self._set_executable(executable_path, "snesim")
        self.variable_name = variable_name

        shape = to_list(shape, 3, 1)
        spacing = to_list(spacing, 3, 1.0)
        origin = to_list(origin, 3, 0.0)

        self.shape = shape
        self.spacing = spacing
        self.origin = origin
        self.n_realizations = n_realizations

        # TODO: Find a way to avoid the copy
        _data = copy(data)
        _training_image = deepcopy(training_image)
        self.is_zero = False
        if 0 in np.unique(_training_image[training_column]):
            self.is_zero = True
            _training_image[training_column] += 1
            if _data is not None:
                if isinstance(_data, np.ndarray) == True:
                    _data[:, data_column] += 1
                else:
                    _data[data_column] += 1
            if class_codes is not None:
                class_codes = tuple(c + 1 for c in class_codes)
        if class_codes is None:
            class_codes = np.unique(_training_image[training_column]).astype(int)
        if proportions is None:
            proportions = (
                np.unique(
                    np.asarray(_training_image[training_column]), return_counts=True
                )[1]
                / _training_image[training_column].size
            )
            sort_indices = np.argsort(class_codes)
            proportions = proportions[sort_indices]
        if max_search_radii is None:
            max_search_radii = (
                spacing[0] * shape[0] / 3.0,
                spacing[1] * shape[1] / 3.0,
                spacing[2] * shape[2] / 3.0,
            )

        self.parameters = [
            (_data, "- file with original data", output_file_name + "_data.dat"),
            (
                _convert_column(
                    _data, to_list(coord_columns, 3, np.nan) + [data_column]
                ),
                "- fcolumns for x, y, z, variable",
            ),
            (len(class_codes), "- number of categories"),
            (class_codes, "- category codes"),
            (proportions, "- (target) global pdf"),
            (
                0 if vert_prop_curve is None else 1,
                "-  use (target) vertical proportions (0=no, 1=yes)",
            ),
            (
                vert_prop_curve,
                "- file with target vertical proportions",
                output_file_name + "_vertprop.dat",
            ),
            (servosystem_param, "- servosystem parameter (0=no correction)"),
            (debugging_level, "- debugging level: 0,1,2,3"),
            (output_file_name + ".dbg", "- debugging file"),
            (output_file_name + ".out", "- file for simulation output"),
            (n_realizations, "- number of realizations to generate"),
            ((shape[0], origin[0], spacing[0]), "- nx,xmn,xsiz"),
            ((shape[1], origin[1], spacing[1]), "- ny,ymn,ysiz"),
            ((shape[2], origin[2], spacing[2]), "- nz,zmn,zsiz"),
            (seed, "- random number seed"),
            (max_data, "- max number of conditioning primary data"),
            (min_replicates, "- min. replicates number"),
            (
                (0 if secondary_data is None else 1, 1 if tau is None else 0),
                "- condition to LP (0=no, 1=yes), flag for iauto",
            ),
            (
                (1.0, 1.0) if tau is None else tau,
                "- two weighting factors to combine P(A|B) and P(A|C)",
            ),
            (
                secondary_data,
                "- file for local properties",
                output_file_name + "_localprop.dat",
            ),
            (
                0 if rotation_affinity is None else 1,
                "- condition to rotation and affinity (0=no, 1=yes)",
            ),
            (
                rotation_affinity,
                "- file for rotation and affinity",
                output_file_name + "_rot_aff.dat",
            ),
            (len(affinity_factors), "- number of affinity categories"),
        ]
        for i, factors in enumerate(affinity_factors):
            self.parameters += [(factors, "- affinity factors (X,Y,Z) icat=" + str(i))]
        self.parameters += [
            (n_multigrids, "- number of multiple grids"),
            (
                _training_image,
                "- file with training image",
                output_file_name + "_train.dat",
            ),
            (
                _training_image[training_column].shape,
                "- training image dimensions: nxtr, nytr, nztr",
            ),
            (
                _convert_column(_training_image, training_column),
                "- column for training variable",
            ),
            (to_list(max_search_radii), "- maximum search radii (hmax,hmin,hvert)"),
            (
                to_list(search_ellipsoid_angles),
                "- angles for search ellipsoid (amax,amin,avert)",
            ),
        ]

    def process_output_files(self):
        """
        Processes the output files before reading.

        TODO: It would be much better to modify the original Fortran code.
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
            elif self.is_zero == True and i > 3:
                match = re.match(r"^(\s*)(\S+)(\s+)(\S+)(\s+)(.*)$", line)
                if match:
                    column1 = match.group(2)
                    if column1 != "-99999":
                        column1 = str(int(column1) - 1)
                    line = f"{match.group(1)}{column1}{match.group(3)}{match.group(4)}{match.group(5)}{match.group(6)}"
            print(line, end="")
