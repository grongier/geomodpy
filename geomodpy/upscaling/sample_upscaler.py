"""Sample upscaler"""

# MIT License

# Copyright (c) 2025 Guillaume Rongier

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


import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_is_fitted

from ..gridding.structured import (
    calculate_cell_centers,
    compute_spiral_offsets,
    find_enclosing_cell_from_nearest_center,
    extract_coordinates,
)
from ..utils import groupby, iter_columns


################################################################################
# Well upscaling


def check_feature_names(X):
    """
    Gets the feature names of an array-like object, if any.
    """
    if isinstance(X, pd.DataFrame):
        return X.columns
    elif isinstance(X, pd.Series):
        return [X.name]
    else:
        return None


def apply_feature_names(X, feature_names):
    """
    Converts an array to a dataframe or a series if some feature names are
    available.
    """
    if feature_names is not None and X.ndim == 2:
        return pd.DataFrame(data=X, columns=feature_names)
    elif feature_names is not None and X.ndim == 1:
        return pd.Series(data=X, name=feature_names)
    else:
        return X


class SampleUpscaler(BaseEstimator):
    """
    Upscaling point samples into 3D structured grids.

    The upscaler uses a k-d tree to find the nearest cell center to each sample,
    then test for actual inclusion following a spiral search.

    Parameters
    ----------
    grid : xarray.Dataset, default=None
        Structured grid represented as a dataset. Either `grid` or `nodes` need
        to be provided.
    nodes : ndarray of shape (3, n_nodes_x, n_nodes_y, n_nodes_z), default=None
        3D coordinates of the nodes of the grid (i.e., the corners of the cells),
        which can be provided instead of `grid`.
    cell_centers : ndarray of shape (3, n_cells_x, n_cells_y, n_cells_z), default=None
        3D coordinates of the centers of the grid cells, calculated from the nodes
        if none are provided.
    x_nodes : str, default='X_nodes'
        Name of the coordinates of the nodes along the x axis in `grid`.
    y_nodes : str, default='Y_nodes'
        Name of the coordinates of the nodes along the y axis in `grid`.
    z_nodes : str, default='Z_nodes'
        Name of the coordinates of the nodes along the z axis in `grid`.
    x : str, default='X'
        Name of the coordinates of the cell centers along the x axis in `grid`.
        If None, the coordinates of the cell centers aren't returned.
    y : str, default='Y'
        Name of the coordinates of the cell centers along the y axis in `grid`.
        If None, the coordinates of the cell centers aren't returned.
    z : str, default='Z'
        Name of the coordinates of the cell centers along the z axis in `grid`.
        If None, the coordinates of the cell centers aren't returned.
    n_cell_radius : int or ndarray of shape (3,), default=5
        Maximum number of cells along each of the three axes of the grid to search
        for cell inclusion. 0 means that only the cell with the nearest center is
        tested, which is fine with regular grids. Irregular grids require a higher
        number.
    merge_outputs : bool, default=True
        If True, merges the `X` and `y` outputs of `resample` into a single
        object; otherwise keeps them as two separate objects.
    """

    def __init__(
        self,
        grid=None,
        nodes=None,
        cell_centers=None,
        x_nodes="X_nodes",
        y_nodes="Y_nodes",
        z_nodes="Z_nodes",
        x="X",
        y="Y",
        z="Z",
        n_cell_radius=5,
        merge_outputs=True,
    ):

        if grid is not None:
            self.nodes, self.cell_centers = extract_coordinates(
                grid, x_nodes=x_nodes, y_nodes=y_nodes, z_nodes=z_nodes, x=x, y=y, z=z
            )
        else:
            self.nodes = np.asarray(nodes)
            self.cell_centers = cell_centers
        if self.cell_centers is None:
            self.cell_centers = calculate_cell_centers(self.nodes)
        self.merge_outputs = merge_outputs

        self._shape = self.cell_centers.shape[1:]
        self._kdtree = KDTree(self.cell_centers.reshape(3, -1).T)
        self._spiral_offsets = compute_spiral_offsets(n_cell_radius)

    def fit(self, X=None, y=None):
        """
        Only here for compatibility.

        Parameters
        ----------
        X : None
            Ignored.
        y : None
            Ignored.

        Returns
        -------
        self : object
            Upscaler.
        """
        return self

    def resample(self, X, y, agg=None):
        """
        Upscales some point samples to the grid cells they lie into.

        Parameters
        ----------
        X : ndarray or DataFrame of shape (n_samples, 3)
            3D coordinates of the samples.
        y : ndarray or DataFrame of shape (n_samples, n_variables)
            Sample variables to upscale.
        agg : str or array-like of shape (n_variables,)
            Upscaling method to use for each variable of `y`:
                - 'mode': the mode.
                - 'gmean', 'geometric mean' or 'geometric': the geometric mean.
                - 'hmean', 'harmonic mean' or 'harmonic': the harmonic mean.
                - any other string (including the default): the arithmetic mean.

        Returns
        -------
        X_out : ndarray or DataFrame of shape (n_upscaled_samples, 3)
            3D coordinates of the upscaled samples.
        y_out : ndarray or DataFrame of shape (n_upscaled_samples, n_variables)
            Upscaled sample variables. Only returned if `merge_outputs` is True.
        """
        check_is_fitted(self, "_kdtree")
        X_names = check_feature_names(X)
        y_names = check_feature_names(y)
        if isinstance(agg, str):
            agg = [agg for i in iter_columns(y)]
        elif agg is None:
            agg = ["mean" for i in iter_columns(y)]

        _, U = self._kdtree.query(X)
        _U = find_enclosing_cell_from_nearest_center(
            np.asarray(X),
            np.stack(np.unravel_index(U, self._shape), axis=-1),
            self.nodes,
            self._spiral_offsets,
        )
        U[_U[:, 0] != -1] = np.ravel_multi_index(_U[_U[:, 0] != -1].T, self._shape)
        U[_U[:, 0] == -1] = -1
        y_out = []
        for c, a in zip(iter_columns(y), agg):
            _c, _u = groupby(c[U != -1], U[U != -1], agg=a)
            y_out.append(_c[:, np.newaxis])
        y_out = np.concatenate(y_out, axis=1)
        X_out = self.cell_centers.reshape(3, -1).T[_u]
        _U = np.unravel_index(_u, self._shape)

        X_out = apply_feature_names(X_out, X_names)
        if isinstance(X, pd.DataFrame):
            X_out["U"], X_out["V"], X_out["W"] = _U
        else:
            X_out = np.concatenate(
                (
                    X_out,
                    _U[0][:, np.newaxis],
                    _U[1][:, np.newaxis],
                    _U[2][:, np.newaxis],
                ),
                axis=-1,
            )
        y_out = apply_feature_names(y_out, y_names)

        if self.merge_outputs == True:
            if isinstance(X, pd.DataFrame):
                return pd.concat([X_out, y_out], axis=1)
            else:
                return np.concatenate([X_out, y_out], axis=1)
        else:
            if isinstance(y, pd.Series):
                y_out = y_out.iloc[:, 0]
            elif y.ndim == 1:
                y_out = y_out[:, 0]
            return X_out, y_out

    def fit_resample(self, X, y, agg=None):
        """
        Upscales some point samples to the grid cells they lie into.

        Parameters
        ----------
        X : ndarray or DataFrame of shape (n_samples, 3)
            3D coordinates of the samples.
        y : ndarray or DataFrame of shape (n_samples, n_variables)
            Sample variables to upscale.
        agg : str or array-like of shape (n_variables,)
            Upscaling method to use for each variable of `y`:
                - 'mode': the mode.
                - 'gmean', 'geometric mean' or 'geometric': the geometric mean.
                - 'hmean', 'harmonic mean' or 'harmonic': the harmonic mean.
                - any other string (including the default): the arithmetic mean.

        Returns
        -------
        X_out : ndarray or DataFrame of shape (n_upscaled_samples, 3)
            3D coordinates of the upscaled samples.
        y_out : ndarray or DataFrame of shape (n_upscaled_samples, n_variables)
            Upscaled sample variables. Only returned if `merge_outputs` is True.
        """
        return self.resample(X, y, agg=agg)
