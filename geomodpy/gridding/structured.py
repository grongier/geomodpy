"""Building structured grids"""

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
import numba as nb
import xarray as xr

from .geometry import centroid_hexahedron, is_point_in_hexahedron


################################################################################
# Grid elements


def get_cell_corners(nodes):
    """
    Get the corner nodes of the cells of a 3D structured grid.

    Parameters
    ----------
    nodes : ndarray of shape (3, n_nodes_x, n_nodes_y, n_nodes_z)
        3D coordinates of the nodes of the grid (i.e., the corners of the cells).

    Returns
    -------
    cell_corners : ndarray of shape (8, n_cells_x, n_cells_y, n_cells_z, 3)
        Cell corners.
    """
    return np.array(
        [
            nodes[:, :-1, :-1, :-1],  # (i,   j,   k)
            nodes[:, 1:, :-1, :-1],  # (i+1, j,   k)
            nodes[:, :-1, 1:, :-1],  # (i,   j+1, k)
            nodes[:, 1:, 1:, :-1],  # (i+1, j+1, k)
            nodes[:, :-1, :-1, 1:],  # (i,   j,   k+1)
            nodes[:, 1:, :-1, 1:],  # (i+1, j,   k+1)
            nodes[:, :-1, 1:, 1:],  # (i,   j+1, k+1)
            nodes[:, 1:, 1:, 1:],  # (i+1, j+1, k+1)
        ]
    )


def calculate_cell_centroids(nodes):
    """
    Calculates the cell centroids of a 3D structured grid from the grid nodes.

    Parameters
    ----------
    nodes : ndarray of shape (3, n_nodes_x, n_nodes_y, n_nodes_z)
        3D coordinates of the nodes of the grid (i.e., the corners of the cells).

    Returns
    -------
    cell_centroids : ndarray of shape (3, n_cells_x, n_cells_y, n_cells_z)
        Cell centroids.
    """
    return np.mean(get_cell_corners(nodes), axis=0)


@nb.njit
def calculate_cell_centers(nodes):
    """
    Calculates the cell centers of mass of a 3D structured grid from the grid nodes.

    Parameters
    ----------
    nodes : ndarray of shape (3, n_nodes_x, n_nodes_y, n_nodes_z)
        3D coordinates of the nodes of the grid (i.e., the corners of the cells).

    Returns
    -------
    cell_centers : ndarray of shape (3, n_cells_x, n_cells_y, n_cells_z)
        Cell centers.
    """
    cell_centers = np.empty(
        (3, nodes.shape[1] - 1, nodes.shape[2] - 1, nodes.shape[3] - 1)
    )

    for i in range(nodes.shape[1] - 1):
        for j in range(nodes.shape[2] - 1):
            for k in range(nodes.shape[3] - 1):
                cell_nodes = np.stack(
                    (
                        nodes[:, i, j, k],
                        nodes[:, i + 1, j, k],
                        nodes[:, i + 1, j + 1, k],
                        nodes[:, i, j + 1, k],
                        nodes[:, i, j, k + 1],
                        nodes[:, i + 1, j, k + 1],
                        nodes[:, i + 1, j + 1, k + 1],
                        nodes[:, i, j + 1, k + 1],
                    )
                )
                cell_centers[:, i, j, k] = centroid_hexahedron(cell_nodes)

    return cell_centers


################################################################################
# Xarray datasets as structured grids


def build_regular_structured_grid(shape, spacing, origin):
    """
    Builds a regular structured grid stored into a Xarray DataSet.

    Parameters
    ----------
    shape : int or array-like (x, y, z)
        Shape of the grid along the x, y, and z axes.
    spacing : float or array-like (x, y, z)
        Cell size along the x, y, and z axes.
    origin : float or array-like (x, y, z)
        Coordinates at the center of the cell at the bottom left of the grid
        along the x, y, and z axes.

    Returns
    -------
    grid : xarray.Dataset
        The grid.
    """
    nodes = np.stack(
        np.meshgrid(
            np.linspace(
                origin[0] - spacing[0] / 2,
                origin[0] + spacing[0] * (shape[0] - 0.5),
                shape[0] + 1,
            ),
            np.linspace(
                origin[1] - spacing[1] / 2,
                origin[1] + spacing[1] * (shape[1] - 0.5),
                shape[1] + 1,
            ),
            np.linspace(
                origin[2] - spacing[2] / 2,
                origin[2] + spacing[2] * (shape[2] - 0.5),
                shape[2] + 1,
            ),
            indexing="ij",
        ),
        axis=0,
    )

    cell_centers = calculate_cell_centers(nodes)

    coords = {
        "X": (["U", "V", "W"], cell_centers[0]),
        "Y": (["U", "V", "W"], cell_centers[1]),
        "Z": (["U", "V", "W"], cell_centers[2]),
        "X_nodes": (["I", "J", "K"], nodes[0]),
        "Y_nodes": (["I", "J", "K"], nodes[1]),
        "Z_nodes": (["I", "J", "K"], nodes[2]),
    }
    attrs = {
        "shape": tuple(shape),
        "spacing": tuple(spacing),
        "origin": tuple(origin),
    }

    return xr.Dataset(coords=coords, attrs=attrs)


def _extract_horizons(horizons, x, y, z):
    """
    Extracts horizons from a set of Xarray datasets and turn them into a NumPy
    array.
    """
    _horizons = np.empty((len(horizons), 3) + horizons[0].shape[:2])
    for i, horizon in enumerate(horizons):
        horizon_x = horizon[x]
        horizon_y = horizon[y]
        if horizon_x.ndim == 1 and horizon_y.ndim == 1:
            horizon_x, horizon_y = np.meshgrid(horizon_x, horizon_y, indexing="ij")
        _horizons[i, 0] = horizon_x
        _horizons[i, 1] = horizon_y
        _horizons[i, 2] = horizon[z][..., 0] if horizon[z].ndim == 3 else horizon[z]

    return _horizons


def _interpolate_nodes_from_horizons(horizons, nz):
    """
    Linearly interpolates horizons to get the nodes of a 3D irregular structured
    grid.
    """
    nodes = np.empty((3, horizons.shape[2], horizons.shape[3], np.sum(nz) + 1))
    start = 0
    for i, n in enumerate(nz):
        alpha = np.linspace(0.0, 1.0, num=n + 1)
        nodes[..., start : start + n + 1] = horizons[i, ..., np.newaxis] + alpha * (
            horizons[i + 1, ..., np.newaxis] - horizons[i, ..., np.newaxis]
        )
        start += n

    return nodes


def build_irregular_structured_grid(horizons, nz, x="X", y="Y", z="Z", eps=1e-8):
    """
    Builds an irregular structured grid stored into a Xarray DataSet. The grid
    follows a pillar-gridding approach, so all horizons must have the same shape.
    Sub-layers are linearly interpolated between the horizons, so the strata are
    proportional (no truncation or onlap).

    Parameters
    ----------
    horizons : array-like of xarray.Datasets
        List of horizons formated as Xarray datasets.
    nz : int or array-likeof shape (n_layers,)
        Number of cells for each layer along the z-axis. The number of layers
        must be equal to the number of horizons plus one.
    x : str, default='X'
        Name of the x-axis coordinates in the horizon datasets.
    y : str, default='Y'
        Name of the x-axis coordinates in the horizon datasets.
    z : str, default='Z'
        Name of the x-axis coordinates in the horizon datasets.
    eps : float, default=1e-8
        Small number used to define the formations from the horizons.

    Returns
    -------
    grid : xarray.Dataset
        The grid.
    """
    if isinstance(nz, int):
        nz = [nz]

    _horizons = _extract_horizons(horizons, x, y, z)
    nodes = _interpolate_nodes_from_horizons(_horizons, nz)

    cell_centers = calculate_cell_centers(nodes)

    # TODO: Make that operation more efficient
    formations = np.zeros(cell_centers.shape[1:], dtype=np.int32)
    for i in range(len(_horizons) - 1, 0, -1):
        top_horizon = _horizons[i].copy()
        top_horizon[2] += eps
        bottom_horizon = _horizons[i].copy()
        bottom_horizon[2] -= eps
        horizon_centers = calculate_cell_centers(
            np.stack((bottom_horizon, top_horizon), axis=-1)
        )
        formations[cell_centers[2] <= horizon_centers[2]] = i

    data_vars = {"Formation": (["U", "V", "W"], formations)}
    coords = {
        "X": (["U", "V", "W"], cell_centers[0]),
        "Y": (["U", "V", "W"], cell_centers[1]),
        "Z": (["U", "V", "W"], cell_centers[2]),
        "X_nodes": (["I", "J", "K"], nodes[0]),
        "Y_nodes": (["I", "J", "K"], nodes[1]),
        "Z_nodes": (["I", "J", "K"], nodes[2]),
    }
    attrs = {
        "shape": tuple(cell_centers.shape[1:]),
        # 'spacing': tuple(spacing),
        # 'origin': tuple(origin),
    }

    return xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)


################################################################################
# Grid operations


@nb.njit
def compute_spiral_offsets(n_cell_radius):
    """
    Computes 3D cell offsets in a structured grid as a spiral loop.

    Parameters
    ----------
    n_cell_radius : int or ndarray of shape (3,)
        Maximum number of cells along each of the three axes of the grid.

    Returns
    -------
    offsets : ndarray of shape (n_offsets, 3)
        3D offsets sorted by distance.
    """
    if isinstance(n_cell_radius, int):
        n_cell_radius = np.array((n_cell_radius, n_cell_radius, n_cell_radius))

    size = (
        (2 * n_cell_radius[0] + 1)
        * (2 * n_cell_radius[1] + 1)
        * (2 * n_cell_radius[2] + 1)
    )
    offsets = np.empty((size, 3), dtype=np.int64)
    distances = np.empty(size)
    i = 0
    for dx in range(-n_cell_radius[0], n_cell_radius[0] + 1):
        for dy in range(-n_cell_radius[1], n_cell_radius[1] + 1):
            for dz in range(-n_cell_radius[2], n_cell_radius[2] + 1):
                offsets[i, 0] = dx
                offsets[i, 1] = dy
                offsets[i, 2] = dz
                distances[i] = dx**2 + dy**2 + dz**2
                i += 1

    return offsets[np.argsort(distances)]


@nb.njit
def find_enclosing_cell_from_nearest_center(
    points, nearest_cell_centers, nodes, spiral_offsets
):
    """
    Find in which cells (if any) some points lie within a grid starting from the
    cells that have the center closest to the points.

    Parameters
    ----------
    points : ndarray of shape (n_points, 3)
        3D coodinates of the points.
    nearest_cell_centers : ndarray of shape (n_points, 3)
        3D indices of the nearest cells to the points.
    nodes : ndarray of shape (3, n_nodes_x, n_nodes_y, n_nodes_z)
        3D coordinates of the nodes of the grid (i.e., the corners of the cells).
    spiral_offsets : ndarray of shape (n_offsets, 3)
        3D offsets for spiral search.

    Returns
    -------
    _nearest_cell_centers : ndarray of shape (n_points,)
        Updated cell indices or -1 if the point lies outside the grid.
    """
    _nearest_cell_centers = -1 * np.ones_like(nearest_cell_centers)

    for i in range(len(points)):
        for j in range(len(spiral_offsets)):

            u = nearest_cell_centers[i, 0] + spiral_offsets[j, 0]
            v = nearest_cell_centers[i, 1] + spiral_offsets[j, 1]
            w = nearest_cell_centers[i, 2] + spiral_offsets[j, 2]

            if (
                u >= 0
                and v >= 0
                and w >= 0
                and u + 1 < nodes.shape[1]
                and v + 1 < nodes.shape[2]
                and w + 1 < nodes.shape[3]
            ):
                cell_nodes = np.stack(
                    (
                        nodes[:, u, v, w],
                        nodes[:, u + 1, v, w],
                        nodes[:, u + 1, v + 1, w],
                        nodes[:, u, v + 1, w],
                        nodes[:, u, v, w + 1],
                        nodes[:, u + 1, v, w + 1],
                        nodes[:, u + 1, v + 1, w + 1],
                        nodes[:, u, v + 1, w + 1],
                    )
                )
                if is_point_in_hexahedron(points[i], cell_nodes) == True:
                    _nearest_cell_centers[i, 0] = u
                    _nearest_cell_centers[i, 1] = v
                    _nearest_cell_centers[i, 2] = w
                    break

    return _nearest_cell_centers
