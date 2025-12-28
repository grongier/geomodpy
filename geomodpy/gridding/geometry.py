"""Geometry"""

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


################################################################################
# 3D volumes


@nb.njit
def volume_tetrahedron(tetrahedron):
    """
    Computes the volume of a tetrathedron.

    Parameters
    ----------
    tetrahedron : ndarray of shape (4, 3)
        3D coodinates of the vertices of the tetrahedron.

    Returns
    -------
    volume : float
        Volume of the tetrathedron.
    """
    tet_matrix = np.empty((3, 3))
    tet_matrix[:, 0] = tetrahedron[0] - tetrahedron[3]
    tet_matrix[:, 1] = tetrahedron[1] - tetrahedron[3]
    tet_matrix[:, 2] = tetrahedron[2] - tetrahedron[3]

    return np.linalg.det(tet_matrix) / 6.0


@nb.njit
def centroid_tetrahedron(tetrahedron):
    """
    Computes the centroid of a tetrathedron.

    Parameters
    ----------
    tetrahedron : ndarray of shape (4, 3)
        3D coodinates of the vertices of the tetrahedron.

    Returns
    -------
    centroid : ndarray of shape (3,)
        Centroid of the tetrathedron.
    """
    centroid = np.zeros(3)
    for i in range(len(tetrahedron)):
        centroid += tetrahedron[i]

    return centroid / 4.0


@nb.njit
def centroid_hexahedron(hexahedron):
    """
    Computes the centroid of a hexahedron.

    Parameters
    ----------
    hexahedron : ndarray of shape (8, 3)
        3D coodinates of the vertices of the hexahedron.

    Returns
    -------
    centroid : ndarray of shape (3,)
        Centroid of the hexahedron.
    """
    centroid = np.zeros(3)
    volume = 0.0
    tetrahedra = np.array(
        (
            (0, 4, 5, 6),
            (0, 1, 5, 6),
            (0, 1, 2, 6),
            (0, 4, 6, 7),
            (0, 3, 6, 7),
            (0, 2, 3, 6),
        )
    )
    for idx in tetrahedra:
        tet_volume = volume_tetrahedron(hexahedron[idx])
        volume += tet_volume
        tet_centroid = centroid_tetrahedron(hexahedron[idx])
        centroid += tet_volume * tet_centroid

    return centroid / volume


@nb.njit
def fractional_offsets_tetrahedron(point, tetrahedron):
    """
    Computes the fractional offsets of a point relative to a tetrahedron.

    Parameters
    ----------
    point : ndarray of shape (3,)
        3D coodinates of the point.
    tetrahedron : ndarray of shape (4, 3)
        3D coodinates of the vertices of the tetrahedron.

    Returns
    -------
    offsets : ndarray of shape (3,)
        Fractional offsets determining the coordinate of the point relative to
        the tetrahedron.

    Reference
    ---------
    Sadarjoen, I. A., de Boer, A. J., Post, F. H., & Mynett, A. E. (1998)
        Particle tracing in σ-transformed grids using tetrahedral 6-decomposition
        https://doi.org/10.1007/978-3-7091-7517-0_7
    """
    tet_matrix = np.empty((3, 3))
    tet_matrix[:, 0] = tetrahedron[0] - tetrahedron[3]
    tet_matrix[:, 1] = tetrahedron[1] - tetrahedron[3]
    tet_matrix[:, 2] = tetrahedron[2] - tetrahedron[3]

    return np.linalg.inv(tet_matrix) @ (point - tetrahedron[3])


@nb.njit
def is_point_in_tetrahedron(point, tetrahedron):
    """
    Checks if a point lies inside a tetrahedron.

    Parameters
    ----------
    point : ndarray of shape (3,)
        3D coodinates of the point.
    tetrahedron : ndarray of shape (4, 3)
        3D coodinates of the vertices of the tetrahedron.

    Returns
    -------
    is_inside : bool
        True if the point is inside, false otherwise.

    Reference
    ---------
    Sadarjoen, I. A., de Boer, A. J., Post, F. H., & Mynett, A. E. (1998)
        Particle tracing in σ-transformed grids using tetrahedral 6-decomposition
        https://doi.org/10.1007/978-3-7091-7517-0_7
    """
    offsets = fractional_offsets_tetrahedron(point, tetrahedron)

    return (
        (offsets[0] >= 0)
        and (offsets[1] >= 0)
        and (offsets[2] >= 0)
        and (offsets[0] + offsets[1] + offsets[2] <= 1)
    )


@nb.njit
def is_point_in_hexahedron(point, hexahedron):
    """
    Checks if a point lies inside a hexahedron using a tetrahedral 6-decomposition.

    Parameters
    ----------
    point : ndarray of shape (3,)
        3D coodinates of the point.
    hexahedron : ndarray of shape (8, 3)
        3D coodinates of the vertices of the hexahedron.

    Returns
    -------
    is_inside : bool
        True if the point is inside, false otherwise.

    Reference
    ---------
    Sadarjoen, I. A., de Boer, A. J., Post, F. H., & Mynett, A. E. (1998)
        Particle tracing in σ-transformed grids using tetrahedral 6-decomposition
        https://doi.org/10.1007/978-3-7091-7517-0_7
    """
    tetrahedra = np.array(
        (
            (0, 4, 5, 6),
            (0, 1, 5, 6),
            (0, 1, 2, 6),
            (0, 4, 6, 7),
            (0, 3, 6, 7),
            (0, 2, 3, 6),
        )
    )
    for idx in tetrahedra:
        if is_point_in_tetrahedron(point, hexahedron[idx]):
            return True
    return False
