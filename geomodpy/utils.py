"""Utils"""

# MIT License

# Copyright (c) 2022-2025 Guillaume Rongier

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
import numpy as np
import pandas as pd
import xarray as xr


################################################################################
# Basics


def is_sequence(seq):
    """
    Check if an object is a sequence. From numpy.distutils.misc_util.
    """
    if isinstance(seq, str):
        return False
    try:
        len(seq)
    except Exception:
        return False
    return True


def to_list(x, size=3, value=None):
    """
    Turns a single input into a list of size `size`. If the input is already a
    sequence but not a list, it is turned into a list.
    """
    x = [x] if is_sequence(x) == False else list(x)
    if value is None:
        x += [x[0]] * (size - len(x))
    else:
        x += [value] * (size - len(x))

    return x


def isnumbers(s):
    """
    Checks if a string contains numbers only.
    """
    s = re.sub(" {2,}", " ", s.lstrip()).split(" ")
    try:
        [float(i) for i in s]
    except ValueError:
        return False
    else:
        return True


def rename_duplicated_columns(df):
    """
    Checks if some columns have the same name and rename the last ones by adding
    an underscore and a number.

    TODO: Remove this function when there's a Pandas' function handling
    duplicated column names
    """
    col_counts = {}
    new_cols = []
    for col in df.columns:
        if col not in col_counts:
            col_counts[col] = 1
            new_cols.append(col)
        else:
            new_col = f"{col}_{col_counts[col]}"
            while new_col in df.columns:
                new_col = f"{col}_{col_counts[col]}"
                col_counts[col] += 1
            new_cols.append(new_col)
            col_counts[col] += 1
    df.columns = new_cols

    return df


################################################################################
# Rotation


def set_rotation_matrix_2d(azimuth=0.0, dip_anisotropy=1.0, eps=1e-20):
    """
    Sets up a 2D anisotropic rotation matrix.

    Parameters
    ----------
    azimuth : float or array-like of shape (n_matrices,), default=0.
        The first rotation angle (in degrees, clockwise) rotates the original
        Y axis (principal direction) in the horizontal plane.
    dip_anisotropy : float or array-like of shape (n_matrices,), default=1.
        Anisotropy for the dip angle.
    eps : float, default=1e-20
        Epsilon to avoid division by zero.

    Returns
    -------
    rotation_matrix : ndarray of shape (2, 2) or (n_matrices, 2, 2)
        Rotation matrix.
    """
    azimuth = np.array([azimuth] if isinstance(azimuth, (int, float)) else azimuth)
    alpha = np.empty(azimuth.shape)
    alpha[(azimuth >= 0.0) & (azimuth < 270.0)] = np.deg2rad(
        90.0 - azimuth[(azimuth >= 0.0) & (azimuth < 270.0)]
    )
    alpha[(azimuth < 0.0) | (azimuth >= 270.0)] = np.deg2rad(
        450.0 - azimuth[(azimuth < 0.0) | (azimuth >= 270.0)]
    )

    # Get the required sines and cosines
    sina = np.sin(alpha)
    cosa = np.cos(alpha)

    # Construct the rotation matrix
    afac1 = 1.0 / np.maximum(dip_anisotropy, eps)
    rotation_matrix = np.empty((alpha.shape[0], 2, 2))
    rotation_matrix[:, 0, 0] = cosa
    rotation_matrix[:, 0, 1] = sina
    rotation_matrix[:, 1, 0] = -afac1 * sina
    rotation_matrix[:, 1, 1] = afac1 * cosa

    return rotation_matrix


def set_rotation_matrix(
    azimuth=0.0,
    dip=0.0,
    plunge=0.0,
    dip_anisotropy=1.0,
    plunge_anisotropy=1.0,
    eps=1e-20,
):
    """
    Sets up an anisotropic rotation matrix.

    Parameters
    ----------
    azimuth : float or array-like of shape (n_matrices,), default=0.
        The first rotation angle (in degrees, clockwise) rotates the original
        Y axis (principal direction) in the horizontal plane.
    dip : float or array-like of shape (n_matrices,), default=0.
        The second rotation angle (in negative degrees down from horizontal)
        rotates the principal direction from the horizontal.
    plunge : float or array-like of shape (n_matrices,), default=0.
        The third rotation angle leaves the principal direction, defined by
        `azimuth` and `dip`, unchanged. The two directions orthogonal to that
        principal direction are rotated clockwise relative to the principal
        direction when looking toward the origin.
    dip_anisotropy : float or array-like of shape (n_matrices,), default=1.
        Anisotropy for the dip angle.
    plunge_anisotropy : float or array-like of shape (n_matrices,), default=1.
        Anisotropy for the plunge angle.
    eps : float, default=1e-20
        Epsilon to avoid division by zero.

    Returns
    -------
    rotation_matrix : ndarray of shape (3, 3) or (n_matrices, 3, 3)
        Rotation matrix.
    """
    # Convert the input angles to three angles that make more mathematical sense:
    #    alpha   angle between the major axis of anisotropy and the
    #            E-W axis. Note: Counter clockwise is positive.
    #    beta    angle between major axis and the horizontal plane.
    #            (The dip of the ellipsoid measured positive down)
    #    theta   Angle of rotation of minor axis about the major axis
    #            of the ellipsoid.
    azimuth = np.array([azimuth] if isinstance(azimuth, (int, float)) else azimuth)
    alpha = np.empty(azimuth.shape)
    alpha[(azimuth >= 0.0) & (azimuth < 270.0)] = np.deg2rad(
        90.0 - azimuth[(azimuth >= 0.0) & (azimuth < 270.0)]
    )
    alpha[(azimuth < 0.0) | (azimuth >= 270.0)] = np.deg2rad(
        450.0 - azimuth[(azimuth < 0.0) | (azimuth >= 270.0)]
    )
    beta = np.deg2rad(-1.0 * np.array(dip))
    theta = np.deg2rad(plunge)

    # Get the required sines and cosines
    sina = np.sin(alpha)
    sinb = np.sin(beta)
    sint = np.sin(theta)
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cost = np.cos(theta)

    # Construct the rotation matrix
    afac1 = 1.0 / np.maximum(dip_anisotropy, eps)
    afac2 = 1.0 / np.maximum(plunge_anisotropy, eps)
    rotation_matrix = np.empty((alpha.shape[0], 3, 3))
    rotation_matrix[:, 0, 0] = cosb * cosa
    rotation_matrix[:, 0, 1] = cosb * sina
    rotation_matrix[:, 0, 2] = -sinb
    rotation_matrix[:, 1, 0] = afac1 * (-cost * sina + sint * sinb * cosa)
    rotation_matrix[:, 1, 1] = afac1 * (cost * cosa + sint * sinb * sina)
    rotation_matrix[:, 1, 2] = afac1 * (sint * cosb)
    rotation_matrix[:, 2, 0] = afac2 * (sint * sina + cost * sinb * cosa)
    rotation_matrix[:, 2, 1] = afac2 * (-sint * cosa + cost * sinb * sina)
    rotation_matrix[:, 2, 2] = afac2 * (cost * cosb)

    return rotation_matrix


################################################################################
# Ellipsoid


def get_ellipsoid_aabb(radius_1, radius_2, radius_3, azimuth=0.0, dip=0.0, plunge=0.0):
    """
    Gets the axis-aligned bounding box of one or several ellipsoids.

    Parameters
    ----------
    radius_1 : float or array-like of shape (n_ellipsoids,)
        Radius of the ellipsoid along the azimuth.
    radius_2 : float or array-like of shape (n_ellipsoids,)
        Radius of the ellipsoid along the dip.
    radius_3 : float or array-like of shape (n_ellipsoids,)
        Radius of the ellipsoid along the plunge.
    azimuth : float or array-like of shape (n_ellipsoids,), default=0.
        The first rotation angle (in degrees, clockwise) rotates the original
        Y axis (principal direction) in the horizontal plane.
    dip : float or array-like of shape (n_ellipsoids,), default=0.
        The second rotation angle (in negative degrees down from horizontal)
        rotates the principal direction from the horizontal.
    plunge : float or array-like of shape (n_ellipsoids,), default=0.
        The third rotation angle leaves the principal direction, defined by
        `azimuth` and `dip`, unchanged. The two directions orthogonal to that
        principal direction are rotated clockwise relative to the principal
        direction when looking toward the origin.

    Returns
    -------
    d : ndarray of shape (n_ellipsoids, 3)
        Dimensions of the axis-aligned bounding boxes of each ellipsoid.

    Source
    ------
    https://tavianator.com/2014/ellipsoid_bounding_boxes.html
    """
    rotation_matrix = set_rotation_matrix(azimuth, dip, plunge, 1.0, 1.0)

    d = np.empty((rotation_matrix.shape[0], 3))
    m11 = rotation_matrix[:, 0, 0] * radius_1
    m12 = rotation_matrix[:, 0, 1] * radius_2
    m13 = rotation_matrix[:, 0, 2] * radius_3
    d[:, 0] = np.sqrt(m11 * m11 + m12 * m12 + m13 * m13)
    m21 = rotation_matrix[:, 1, 0] * radius_1
    m22 = rotation_matrix[:, 1, 1] * radius_2
    m23 = rotation_matrix[:, 1, 2] * radius_3
    d[:, 1] = np.sqrt(m21 * m21 + m22 * m22 + m23 * m23)
    m31 = rotation_matrix[:, 2, 0] * radius_1
    m32 = rotation_matrix[:, 2, 1] * radius_2
    m33 = rotation_matrix[:, 2, 2] * radius_3
    d[:, 2] = np.sqrt(m31 * m31 + m32 * m32 + m33 * m33)

    return d


################################################################################
# Data management


def _groupby_mode(values, unique_groups, inverse):
    """
    Computes the mode of a set of grouped values.
    """
    if len(values) > 0:
        vmin = values.min()
        counts = np.zeros((len(unique_groups), values.max() - vmin + 1), dtype=int)
        np.add.at(counts, (inverse, values - vmin), 1)

        return vmin + np.argmax(counts, axis=1)
    else:
        return np.array([])


def _groupby_mean(values, unique_groups, inverse, agg):
    """
    Computes the mean of a set of grouped values.
    """
    counts = np.zeros(len(unique_groups), dtype=int)
    sums = np.zeros(len(unique_groups), dtype=values.dtype)
    np.add.at(counts, inverse, 1)

    if agg == "geometric" or agg == "geometric mean" or agg == "gmean":
        np.add.at(sums, inverse, np.log(values))
        return np.exp(sums / counts)
    elif agg == "harmonic" or agg == "harmonic mean" or agg == "hmean":
        np.add.at(sums, inverse, 1.0 / values)
        return counts / sums
    else:
        np.add.at(sums, inverse, values)
        return sums / counts


def groupby(values, groups, agg="mean"):
    """
    Computes aggregated values based on some groups.

    Parameters
    ----------
    values : ndarray of shape (n_values,)
        Values to compute the aggregation on.
    groups : ndarray of shape (n_values,)
        Groups to which each value belongs to.
    agg : str, default='mean'
        The aggregation method to use, which can be:
            - 'mode': the mode.
            - 'gmean', 'geometric mean' or 'geometric': the geometric mean.
            - 'hmean', 'harmonic mean' or 'harmonic': the harmonic mean.
            - any other string (including the default): the arithmetic mean.

    Returns
    -------
    agg_values : ndarray of shape (n_unique_groups,)
        The aggregated values.
    unique_groups : ndarray of shape (n_unique_groups,)
        The unique groups corresponding to the aggregated values.
    """
    values = np.asarray(values)
    groups = np.asarray(groups)
    groups = groups[~np.isnan(values)]
    values = values[~np.isnan(values)]

    unique_groups, inverse = np.unique(groups, return_inverse=True)

    if agg == "mode":
        return _groupby_mode(values, unique_groups, inverse), unique_groups
    else:
        return _groupby_mean(values, unique_groups, inverse, agg), unique_groups


def iter_columns(x):
    """
    Iterator over the columns of an array-like object.

    Parameters
    ----------
    x : array-like of shape (n_values, n_variables)
        2D array-like object.

    Returns
    -------
    iter : iterator
        Iterator over the columns of the array-like object.
    """
    if isinstance(x, pd.DataFrame):
        for col in x.columns:
            yield x[col].values
    elif isinstance(x, pd.Series):
        yield x.values
    else:
        _x = np.asarray(x)
        if _x.ndim == 1:
            yield _x
        elif _x.ndim == 2:
            for col in _x.T:
                yield col
        else:
            raise ValueError("Input array must be 1D or 2D.")


def create_dataset(x, spacing, origin, n_realizations=None, var_name=None):
    """
    Creates a Xarray DataSet from an array.

    Parameters
    ----------
    x : array-like of shape (n_z, n_y, n_x), (n_z, n_y, n_x, n_vars), (n_realizations, n_z, n_y, n_x), or (n_realizations, n_z, n_y, n_x, n_vars)
        Array containing the data.
    spacing : float or array-like (x, y, z)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.
    n_realizations : int, default=None
        Number of realizations in `x`.
    var_name : string or array-like of shape (n_vars,), default=None
        Name(s) of the variable(s) in `x`.

    Returns
    -------
    ds : xarray.DataSet
        The DataSet containing `x`.
    """
    if x.ndim >= 4 and n_realizations is not None:
        shape = x.shape[1:4]
    else:
        shape = x.shape[:3]
    spacing = to_list(spacing)
    origin = to_list(origin)
    if x.ndim == 3 or (x.ndim == 4 and n_realizations is not None):
        x = x[..., np.newaxis]
    if var_name is None:
        var_name = ["Variable_" + str(i + 1) for i in range(x.shape[-1])]
    elif isinstance(var_name, str):
        var_name = [var_name]

    if n_realizations is None:
        dims = ["U", "V", "W"]
    else:
        dims = ["Realization", "U", "V", "W"]
    ds = xr.Dataset(
        data_vars={name: (dims, x[..., i]) for i, name in enumerate(var_name)},
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
    ds.attrs["shape"] = shape
    ds.attrs["spacing"] = tuple(spacing)
    ds.attrs["origin"] = tuple(origin)

    return ds
