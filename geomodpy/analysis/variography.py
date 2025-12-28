"""Variography"""

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


from dataclasses import dataclass
from typing import Union
import math
import numpy as np
import pandas as pd

from ..utils import set_rotation_matrix, set_rotation_matrix_2d
from .distances import squared_anisotropic_distance, distance_matrix


################################################################################
# Semivariogram model


@dataclass(frozen=True)
class Direction:
    """
    Parameters for the direction along which to compute a variogram.

    Parameters
    ----------
    azimuth : float or array-like, default=0.
        Azimuth or ensemble of azimuths.
    azimuth_tolerance : float, default=90.
        Tolerance around the azimuth.
    horizontal_bandwidth : float, default=50.
        Maximum horizontal deviation from the direction vector.
    dip : float, default=0.
        Dip.
    dip_tolerance : float, default=90.
        Tolerance around the dip.
    vertical_bandwidth : float, default=50.
        Maximum vertical deviation from the direction vector.
    n_lags : int, default=10
        Number of lags.
    lag_separation_distance : float, default=5.
        Separation distance between each lag.
    lag_tolerance : float, default=0.66
        Tolerance when determining which data ends in a lag, in fraction of
        `lag_separation_distance`.
    sampling_fraction : float, default=0.1
        Fraction of `lag_separation_distance` used when computing the values of
        a variogram model.
    """

    azimuth: Union[float, tuple, list, np.ndarray] = 0.0
    azimuth_tolerance: float = 90.0
    horizontal_bandwidth: float = 50.0
    dip: float = 0.0
    dip_tolerance: float = 90.0
    vertical_bandwidth: float = 50.0
    n_lags: int = 10
    lag_separation_distance: float = 5.0
    lag_tolerance: float = 0.66
    sampling_fraction: float = 0.1


@dataclass(frozen=True)
class VarioStruct:
    """
    Parameters for a variogram structure.

    Parameters
    ----------
    model : str, default='spherical'
        The type of structure, either `spherical`, `exponential`, `gaussian`,
        `power`, or `hole effect`.
    partial_sill : float, default=1.
        The partial sill, i.e., the variance contribution for the structure.
    azimuth : float, default=0.
        The first rotation angle (in degrees, clockwise) rotates the original
        Y axis (principal direction) in the horizontal plane.
    dip : float, default=0.
        The second rotation angle (in negative degrees down from horizontal)
        rotates the principal direction from the horizontal.
    plunge : float, default=0.
        The third rotation angle leaves the principal direction, defined by
        `azimuth` and `dip`, unchanged. The two directions orthogonal to that
        principal direction are rotated clockwise relative to the principal
        direction when looking toward the origin.
    range_hmax : float, default=10.
        The maximum horizontal range of the structure.
    range_hmin : float, default=10.
        The minimum horizontal range of the structure.
    range_vert : float, default=10.
        The vertical range of the structure.
    """

    model: str = "spherical"
    partial_sill: float = 1.0
    azimuth: float = 0.0
    dip: float = 0.0
    plunge: float = 0.0
    range_hmax: float = 10.0
    range_hmin: float = 10.0
    range_vert: float = 10.0

    def __post_init__(self):
        """
        Converts the variogram model to the right GSLIB index.
        """
        models = {
            "spherical": 1,
            "exponential": 2,
            "gaussian": 3,
            "power": 4,
            "hole effect": 5,
        }
        object.__setattr__(self, "model", models[self.model])


class VarioModel:
    """
    Variogram model specification.

    Parameters
    ----------
    nugget_effect : float, default=0.
        Isotropic nugget constant.
    structures : VarioStruct or array-like of shape (n_structs,), default=VarioStruct()
        Variogram structure or set of variogram structures.
    """

    def __init__(self, nugget_effect=0.0, structures=VarioStruct()):

        self.nugget_effect = nugget_effect
        self.structures = (
            structures if isinstance(structures, (tuple, list)) else [structures]
        )

    def max_covariance(self):
        """
        Computes the maximum covariance.
        """
        max_covariance = self.nugget_effect
        for structure in self.structures:
            if structure.model == 4:
                max_covariance += 999.0
            else:
                max_covariance += structure.partial_sill

        return max_covariance

    def rotation_matrix_2d(self, eps=1e-20):
        """
        Computes the 2D rotation matrix for this variogram model.

        Parameters
        ----------
        eps : float, default=1e-20
            Epsilon to avoid division by zero.
        """
        return set_rotation_matrix_2d(
            azimuth=[structure.azimuth for structure in self.structures],
            dip_anisotropy=[
                structure.range_hmin / max(structure.range_hmax, eps)
                for structure in self.structures
            ],
        )

    def rotation_matrix_3d(self, eps=1e-20):
        """
        Computes the 3D rotation matrix for this variogram model.

        Parameters
        ----------
        eps : float, default=1e-20
            Epsilon to avoid division by zero.
        """
        return set_rotation_matrix(
            azimuth=[structure.azimuth for structure in self.structures],
            dip=[structure.dip for structure in self.structures],
            plunge=[structure.plunge for structure in self.structures],
            dip_anisotropy=[
                structure.range_hmin / max(structure.range_hmax, eps)
                for structure in self.structures
            ],
            plunge_anisotropy=[
                structure.range_vert / max(structure.range_hmax, eps)
                for structure in self.structures
            ],
        )

    def _covariance(self, distance):
        """
        Calculates the covariance from a distance.
        """
        max_covariance = self.max_covariance()

        covariance = np.zeros(1 if distance.ndim == 1 else distance.shape[1:])
        for i, structure in enumerate(self.structures):
            # Spherical model
            if structure.model == 1:
                reduced_distance = distance[i] / structure.range_hmax
                covariance[reduced_distance < 1.0] += structure.partial_sill * (
                    1.0
                    - reduced_distance[reduced_distance < 1.0]
                    * (
                        1.5
                        - 0.5
                        * reduced_distance[reduced_distance < 1.0]
                        * reduced_distance[reduced_distance < 1.0]
                    )
                )
            # Exponential model
            elif structure.model == 2:
                covariance += structure.partial_sill * np.exp(
                    -3.0 * distance[i] / structure.range_hmax
                )
            # Gaussian model
            elif structure.model == 3:
                reduced_distance = distance[i] / structure.range_hmax
                covariance += structure.partial_sill * np.exp(
                    -3.0 * reduced_distance * reduced_distance
                )
            # Power model
            elif structure.model == 4:
                covariance += max_covariance - structure.partial_sill * (
                    distance[i] ** structure.range_hmax
                )
            # Hole effect model
            elif structure.model == 5:
                reduced_distance = distance[i] / structure.range_hmax
                d = 10.0 * structure.range_hmax
                covariance += (
                    structure.partial_sill
                    * np.exp(-3.0 * distance[i] / d)
                    * np.cos(reduced_distance * np.pi)
                )
                covariance += structure.partial_sill * np.cos(reduced_distance * np.pi)

        return covariance

    def covariance(self, point_1, point_2):
        """
        Calculates the covariance associated with a variogram model specified by
        a nugget effect and nested varigoram structures. The anisotropy definition
        can be different for each nested structure.

        Parameters
        ----------
        point_1 : array-like of shape (3,) or (n_points, 3)
            First point or set of points.
        point_2 : array-like of shape (3,) or (n_points, 3)
            Second point or set of points. If `point_1` and `point_2` contain
            multiple points, they must have the same shape and the covariance is
            computed for each pair of points.

        Returns
        -------
        covariance : ndarray of shape (1,) or (n_points,)
            Covariance.
        """
        if point_1.shape[-1] == 2:
            rotation_matrix = self.rotation_matrix_2d()
        else:
            rotation_matrix = self.rotation_matrix_3d()

        distance = np.sqrt(
            squared_anisotropic_distance(point_1, point_2, rotation_matrix)
        )

        return self._covariance(distance)

    def covariance_matrix(
        self,
        X,
        Y=None,
        vertices=None,
        discontinuities=None,
        visibility_graph=None,
        in_discontinuities=None,
    ):
        """
        Calculates the covariance associated with a variogram model specified by
        a nugget effect and nested variogram structures. The anisotropy definition
        can be different for each nested structure.

        Parameters
        ----------
        X : array-like of shape (n_points, 3)
            Coordinates of the data.
        Y : array-like of shape (n_points, 3), default=None
            Coordinates of the data.
        vertices : ndarray of shape (n_vertices, 3), default=None
            Coordinates of the vertices of the graph.
        discontinuities : list of ndarrays, default=None
            Set of discontinuities each represented by a list of tuple, each
            tuple symbolizing an edge of the polygon through the indices of its
            starting and ending vertices in `vertices`.
        visibility_graph: dict, default=None
            Visibility graph as a dictionary whose keys are the indices of the
            vertices forming an edge of the graph and values are the length of
            the edge.
        in_discontinuities : ndarray of shape (n_vertices,), default=None
            Indices of the polygon to which each vertex belongs.

        Returns
        -------
        covariance : ndarray of shape (n_points,)
            The covariance.
        """
        if X.shape[1] == 2:
            rotation_matrix = self.rotation_matrix_2d()
        else:
            rotation_matrix = self.rotation_matrix_3d()

        distance = distance_matrix(
            X,
            Y=Y,
            rotation_matrix=rotation_matrix,
            vertices=vertices,
            discontinuities=discontinuities,
            visibility_graph=visibility_graph,
            in_discontinuities=in_discontinuities,
        )

        return self._covariance(distance)

    def semivario(self, direction=Direction()):
        """
        Calculates a semivariogram model values at certain lags.

        Parameters
        ----------
        direction : Direction or array-like of shape (n_dir,), default=Direction()
            Direction(s) along which to compute the variogram.

        Returns
        -------
        values : DataFrame
            Values of the model at the lags.
        """
        if isinstance(direction, (tuple, list)) == False:
            direction = [direction]

        values = []
        for i, d in enumerate(direction):
            azimuths = d.azimuth
            if isinstance(azimuths, (int, float)) == True:
                azimuths = [azimuths]
            for azimuth in azimuths:

                n_lags = int(d.n_lags / d.sampling_fraction)
                lags = np.zeros((n_lags + 2, 3))
                lags[1] = 0.0001
                lags[2:] = np.tile(np.arange(1, n_lags + 1), (3, 1)).T

                offset = np.empty(3)
                offset[0] = (
                    math.sin(math.radians(azimuth))
                    * math.cos(math.radians(d.dip))
                    * d.lag_separation_distance
                    * d.sampling_fraction
                )
                offset[1] = (
                    math.cos(math.radians(azimuth))
                    * math.cos(math.radians(d.dip))
                    * d.lag_separation_distance
                    * d.sampling_fraction
                )
                offset[2] = (
                    math.sin(math.radians(d.dip))
                    * d.lag_separation_distance
                    * d.sampling_fraction
                )
                point = lags * offset

                cov = self.covariance(np.zeros(3), point)
                max_covariance = self.max_covariance()
                _values = {
                    "Azimuth": azimuth,
                    "Dip": d.dip,
                    "Distance": np.sqrt(np.sum(point * point, axis=1)),
                    "Semivariance": max_covariance - cov,
                    "Covariance": cov,
                    "Correlation": cov / max_covariance,
                }
                values.append(pd.DataFrame(_values))

        return pd.concat(values)

    def flatten(self):
        """
        Returns the variogram parameters as a single array.
        """
        parameters = [self.nugget_effect, len(self.structures)]
        for structure in self.structures:
            parameters += [
                structure.model,
                structure.partial_sill,
                structure.azimuth,
                structure.dip,
                structure.plunge,
                structure.range_hmax,
                structure.range_hmin,
                structure.range_vert,
            ]

        return np.array(parameters)
