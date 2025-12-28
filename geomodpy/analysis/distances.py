"""Distances"""

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


from itertools import chain
import math
import heapq
import numpy as np
import numba as nb

from ..utils import set_rotation_matrix_2d, set_rotation_matrix


################################################################################
# Basics


@nb.njit()
def isclose(a, b, rel_tol=1e-12, abs_tol=0.0):
    """
    Checks whether two values are close.
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


################################################################################
# Anisotropic distance


@nb.njit()
def anisotropic_distance(point_1, point_2, rotation_matrix):
    """
    Computes the anisotropic distance between two points given the coordinates
    of each point and a definition of the anisotropy.

    Parameters
    ----------
    point_1 : array-like of shape (2,) or (3,)
        Coordinates of the first point.
    point_2 : array-like of shape (2,) or (3,)
        Coordinates of the second point.
    rotation_matrix : array-like of shape (2, 2) or (3, 3)
        Rotation matrix.

    Returns
    -------
    distance : float
        The anisotropic distance.
    """
    point = rotation_matrix.dot(point_1 - point_2)
    distance = 0.0
    for u in range(len(point)):
        distance += point[u] * point[u]

    return math.sqrt(distance)


@nb.njit(parallel=True, nogil=True)
def _anisotropic_autodistance_matrix(X, rotation_matrix):
    """
    Calculates the anisotropic distance matrix the coordinates and a definition
    of the anisotropy.
    """
    distances = np.zeros((len(rotation_matrix), len(X), len(X)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(j + 1, len(X)):
                distances[k, j, i] = anisotropic_distance(
                    X[j], X[i], rotation_matrix[k]
                )
                distances[k, i, j] = distances[k, j, i]

    return distances


@nb.njit(parallel=True, nogil=True)
def _anisotropic_distance_matrix(X, Y, rotation_matrix):
    """
    Calculates the anisotropic distance matrix the coordinates and a definition
    of the anisotropy.
    """
    distances = np.empty((len(rotation_matrix), len(X), len(Y)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(len(Y)):
                distances[k, j, i] = anisotropic_distance(
                    X[j], Y[i], rotation_matrix[k]
                )

    return distances


def anisotropic_distance_matrix(X, Y=None, rotation_matrix=None):
    """
    Calculates the anisotropic distance matrix the coordinates and a definition
    of the anisotropy.

    Parameters
    ----------
    X : array-like of shape (n_points, 2) or (n_points, 3)
        Coordinates of the data.
    Y : array-like of shape (n_points, 2) or (n_points, 3), default=None
        Coordinates of the data.
    rotation_matrix : array-like of shape (n_rotations, 2, 2) or (n_rotations, 3, 3), default=None
        Rotation matrix.

    Returns
    -------
    distances : ndarray of shape (n_rotations, n_points)
        Squared anisotropic distance matrix.
    """
    if rotation_matrix is None:
        if X.shape[1] == 2:
            rotation_matrix = set_rotation_matrix_2d()
        else:
            rotation_matrix = set_rotation_matrix()

    if Y is None:
        return _anisotropic_autodistance_matrix(X, rotation_matrix)
    else:
        return _anisotropic_distance_matrix(X, Y, rotation_matrix)


def squared_anisotropic_distance(point_1, point_2, rotation_matrix):
    """
    Calculates the anisotropic distance between two points given the coordinates
    of each point and a definition of the anisotropy.

    Parameters
    ----------
    point_1 : array-like of shape (3,) or (n_points, 3)
        Coordinates of the first point.
    point_2 : array-like of shape (3,) or (n_points, 3)
        Coordinates of the second point.
    rotation_matrix : array-like of shape (3, 3) or (n_rotations, 3, 3)
        Rotation matrix.

    Returns
    -------
    distance : ndarray of shape (1,) or (n_points,) or (n_rotations,) or (n_rotations, n_points)
        Squared anisotropic distance.
    """
    return np.sum(np.dot(rotation_matrix, np.subtract(point_1, point_2).T) ** 2, axis=1)


################################################################################
# Discontinuity pre-processing


@nb.njit()
def _polygonize_discontinuity(discontinuity, eps):
    """
    Turns a discontinuity of more than one segment into polygons by
    duplicating and shifting its internal edges.
    """
    discontinuity_1 = np.empty((len(discontinuity), len(discontinuity[0])))
    discontinuity_1[0] = discontinuity[0]
    discontinuity_2 = np.empty((len(discontinuity) - 1, len(discontinuity[0])))
    discontinuity_2[0] = discontinuity[0]
    for i in range(1, len(discontinuity) - 1):
        u = discontinuity[i - 1] - discontinuity[i]
        u /= np.linalg.norm(u)
        v = discontinuity[i + 1] - discontinuity[i]
        v /= np.linalg.norm(v)
        w = u + v
        o = perp(u, v)
        w = eps * o * w / np.linalg.norm(w) / abs(o)
        discontinuity_1[i] = discontinuity[i] + w
        discontinuity_2[i] = discontinuity[i] - w
    discontinuity_1[-1] = discontinuity[-1]

    return np.concatenate((discontinuity_1, discontinuity_2[::-1]))


def polygonize_discontinuities(discontinuities, eps=1e-8):
    """
    Turns discontinuities of more than one segment into polygons by
    duplicating and shifting their internal edges.

    Parameters
    ----------
    discontinuities : array-like
        Set of discontinuities each represented by a list of 3D coordinates
        symbolizing the vertices of the discontinuities' edges.
    eps : float
        Shift to apply to the duplicated vertices to create a polygon.

    Returns
    -------
    discontinuities : list
        Set of discontinuities turned into polygons if needed.
    """
    _discontinuities = nb.typed.List()
    for discontinuity in discontinuities:
        if len(discontinuity) > 2 and np.any(discontinuity[0] != discontinuity[-1]):
            _discontinuities.append(
                _polygonize_discontinuity(np.asarray(discontinuity), eps)
            )
        else:
            _discontinuities.append(np.asarray(discontinuity))

    return _discontinuities


@nb.njit()
def _extract_polygons(vertices, discontinuities):
    """
    Extracts a list of polygons from discontinuities.
    """
    in_polygons = np.full(len(vertices), -1, dtype=nb.int64)
    polygons = nb.typed.List()
    for i, discontinuity in enumerate(discontinuities):
        polygon = nb.typed.List()
        for j, node in enumerate(discontinuity[:-1]):
            u = np.argmax(np.sum(vertices == node, axis=1))
            v = np.argmax(np.sum(vertices == discontinuity[j + 1], axis=1))
            polygon.append((u, v))
            in_polygons[u] = i
            in_polygons[v] = i
        polygons.append(polygon)

    return polygons, in_polygons


def extract_discontinuity_inputs(discontinuities, extra_vertices=None):
    """
    Extracts the inputs needed to compute distances around discontinuities.

    Parameters
    ----------
    discontinuities : array-like
        Set of discontinuities each represented by a list of 3D coordinates
        symbolizing the vertices of the discontinuities' edges.

    Returns
    -------
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices delimiting the discontinuities.
    polygons : list
        Discontinuities turned into polygons. Each polygon is a list of tuple,
        each tuple symbolizing an edge of the polygon through the indices
        of its starting and ending vertices in `vertices`.
    in_polygons : ndarray of shape (n_vertices,)
        Indices of the polygon to which each vertex belongs.
    """
    vertices = list(chain(*discontinuities))
    if extra_vertices is not None:
        vertices += extra_vertices
    vertices = np.unique(vertices, axis=0)

    polygons, in_polygons = _extract_polygons(vertices, discontinuities)

    return vertices, polygons, in_polygons


def format_discontinuities(discontinuities):
    """
    Formats some discontinuities to be usable by Numba.

    Parameters
    ----------
    discontinuities : list of lists
        Set of discontinuities each represented by a list of 3D coordinates
        symbolizing the vertices of the discontinuities' edges.

    Returns
    -------
    discontinuities : list of ndarrays
        Set of discontinuities each represented by a list of 3D coordinates
        symbolizing the vertices of the discontinuities' edges.
    """
    _discontinuities = nb.typed.List()
    for discontinuitie in discontinuities:
        _discontinuities.append(nb.typed.List(np.asarray(discontinuitie)))

    return _discontinuities


################################################################################
# Visibility graph


@nb.njit()
def perp(u, v):
    """
    Compares the orientation of two vectors.
    """
    return u[0] * v[1] - u[1] * v[0]


@nb.njit()
def intersect(vertex_1_1, vertex_1_2, vertex_2_1, vertex_2_2, include_extremities=True):
    """
    Checks whether two segments represented by the vertices at their extremities
    intersects.
    """
    u = vertex_1_2 - vertex_1_1
    v = vertex_2_2 - vertex_2_1
    w = vertex_1_1 - vertex_2_1
    D = perp(u, v)

    # Check if the segments are parallel
    if isclose(D, 0.0, abs_tol=1e-9):
        return None
    sI = perp(v, w) / D
    # Check if the potential intersection is outside of the segment
    if (sI < 0.0 or (isclose(sI, 0.0) and include_extremities == False)) or (
        sI > 1.0 or (isclose(sI, 1.0) and include_extremities == False)
    ):
        return None
    tI = perp(u, w) / D
    # Check if the potential intersection is outside of the segment
    if (tI < 0.0 or (isclose(tI, 0.0) and include_extremities == False)) or (
        tI > 1.0 or (isclose(tI, 1.0) and include_extremities == False)
    ):
        return None
    return vertex_1_1 + sI * u


@nb.njit()
def are_vertices_visible(vertices, i, j, polygons):
    """
    Checks whether two vertices are visible from one another or if a polygon
    stands in the way.
    """
    for polygon in polygons:
        for edge in polygon:
            intersection = intersect(
                vertices[i],
                vertices[j],
                vertices[edge[0]],
                vertices[edge[1]],
                include_extremities=False,
            )
            if intersection is not None:
                return False
    return True


@nb.njit()
def det(vertex_1, vertex_2, point):
    """
    Computes the determinant between a segment and a point.
    """
    return (vertex_1[0] - point[0]) * (vertex_2[1] - point[1]) - (
        vertex_2[0] - point[0]
    ) * (vertex_1[1] - point[1])


@nb.njit()
def is_point_in_polygon(point, vertices, polygon):
    """
    Checks if a point is inside a polygon. See:
    https://doi.org/10.1016/S0925-7721(01)00012-8
    """
    if (
        isclose(vertices[polygon[0][0]][1], point[1]) == True
        and isclose(vertices[polygon[0][0]][0], point[0]) == True
    ):
        return False

    winding_number = 0
    for edge in polygon:
        if isclose(vertices[edge[1]][1], point[1]) == True:
            if isclose(vertices[edge[1]][0], point[0]) == True:
                return False
            elif isclose(vertices[edge[0]][1], point[1]) == True and (
                vertices[edge[1]][0] > point[0]
            ) == (vertices[edge[0]][0] < point[0]):
                return False
        if (vertices[edge[0]][1] < point[1]) != (vertices[edge[1]][1] < point[1]):
            if vertices[edge[0]][0] >= point[0]:
                if vertices[edge[1]][0] > point[0]:
                    winding_number += (
                        2 * (vertices[edge[1]][1] > vertices[edge[0]][1]) - 1
                    )
                else:
                    d = det(vertices[edge[0]], vertices[edge[1]], point)
                    if isclose(d, 0.0, abs_tol=1e-9):
                        return False
                    elif (d > 0.0) == (vertices[edge[1]][1] > vertices[edge[0]][1]):
                        winding_number += (
                            2 * (vertices[edge[1]][1] > vertices[edge[0]][1]) - 1
                        )
            elif vertices[edge[1]][0] > point[0]:
                d = det(vertices[edge[0]], vertices[edge[1]], point)
                if isclose(d, 0.0, abs_tol=1e-9):
                    return False
                elif (d > 0.0) == (vertices[edge[1]][1] > vertices[edge[0]][1]):
                    winding_number += (
                        2 * (vertices[edge[1]][1] > vertices[edge[0]][1]) - 1
                    )

    if winding_number == 0:
        return False
    return True


@nb.njit()
def is_edge_inside(vertices, i, j, in_polygons, polygons):
    """
    Checks if the edge of a polygon is inside the polygon.
    """
    if (
        in_polygons[i] != -1
        and in_polygons[j] != -1
        and in_polygons[i] == in_polygons[j]
        and len(polygons[in_polygons[i]]) > 2
    ):
        middle = (vertices[i] + vertices[j]) / 2.0
        return is_point_in_polygon(middle, vertices, polygons[in_polygons[i]])
    return False


@nb.njit()
def is_vertex_reflex(vertices, i, polygon):
    """
    Checks if the vertex of a polygon is reflex.
    """
    for j in range(len(polygon)):
        if len(polygon) > 2 and i in polygon[j] and i in polygon[j - 1]:
            u = vertices[polygon[j][1]] - vertices[polygon[j][0]]
            v = vertices[polygon[j - 1][0]] - vertices[polygon[j - 1][1]]
            if perp(u, v) <= 0.0:
                return False
    return True


@nb.njit()
def is_edge_bitangent(vertices, i, j, in_polygons, polygons, eps):
    """
    Checks if the edge of a polygon is bitangent.
    """
    u = vertices[j] - vertices[i]
    u /= np.linalg.norm(u)

    is_bitangent = False
    if in_polygons[j] != -1 and len(polygons[in_polygons[j]]) > 2:
        extension = vertices[j] + eps * u
        is_bitangent = is_point_in_polygon(
            extension, vertices, polygons[in_polygons[j]]
        )
    if (
        in_polygons[i] != -1
        and len(polygons[in_polygons[i]]) > 2
        and is_bitangent == False
    ):
        extension = vertices[i] - eps * u
        is_bitangent = is_point_in_polygon(
            extension, vertices, polygons[in_polygons[i]]
        )

    return is_bitangent


@nb.njit()
def build_visibility_graph(
    vertices, polygons, in_polygons, rotation_matrix, reduced=True, eps=1e-5
):
    """
    Builds a visibility graph connecting some vertices around obstacles
    represented by polygons.

    Parameters
    ----------
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    polygons : list
        Discontinuities turned into polygons. Each polygon is a list of tuple,
        each tuple symbolizing an edge of the polygon through the indices
        of its starting and ending vertices in `vertices`.
    in_polygons : ndarray of shape (n_vertices,)
        Indices of the polygon to which each vertex belongs.
    rotation_matrix : array-like of shape (3, 3)
        Rotation matrix.
    reduced : bool, default=True
        If True, returns a reduced visibility graph, i.e., a graph without
        the edges that would never be used in a shortest path.
    eps: float, default=1e-5
        Extra length added to a segment to check whether it is bitangent.

    Returns
    -------
    visibility_graph: dict
        The visibility graph as a dictionary whose keys are the indices of
        the vertices forming an edge of the graph and values are the length
        of the edge.
    """
    is_reflex_vertices = None
    if reduced == True:
        is_reflex_vertices = np.empty(len(vertices), dtype=nb.int8)
        for i in range(len(vertices)):
            is_reflex_vertices[i] = is_vertex_reflex(
                vertices, i, polygons[in_polygons[i]]
            )

    visibility_graph = dict()
    for i in range(len(vertices)):
        if reduced == False or is_reflex_vertices[i] == 1:
            for j in range(i + 1, len(vertices)):
                if reduced == False or is_reflex_vertices[j] == 1:
                    if (
                        are_vertices_visible(vertices, i, j, polygons) == True
                        and is_edge_inside(vertices, i, j, in_polygons, polygons)
                        == False
                        and (
                            reduced == False
                            or is_edge_bitangent(
                                vertices, i, j, in_polygons, polygons, eps
                            )
                            == False
                        )
                    ):
                        visibility_graph[(i, j)] = anisotropic_distance(
                            vertices[i], vertices[j], rotation_matrix
                        )
                        visibility_graph[(j, i)] = visibility_graph[(i, j)]

    return visibility_graph


@nb.njit()
def update_visibility_graph(vertices, visibility_graph, rotation_matrix):
    """
    Updates the distances in a visibility graph based on a new rotation matrix.

    Parameters
    ----------
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    visibility_graph: dict
        Visibility graph as a dictionary whose keys are the indices of the
        vertices forming an edge of the graph and values are the length of
        the edge.
    rotation_matrix : array-like of shape (3, 3)
        Rotation matrix.

    Returns
    -------
    _visibility_graph: dict
        The visibility graph as a dictionary whose keys are the indices of
        the vertices forming an edge of the graph and values are the length
        of the edge.
    """
    _visibility_graph = dict()
    for key in visibility_graph:
        _visibility_graph[key] = anisotropic_distance(
            vertices[key[0]], vertices[key[1]], rotation_matrix
        )

    return _visibility_graph


@nb.njit()
def build_visibility_graphs(
    vertices, polygons, in_polygons, rotation_matrix, reduced=True, eps=1e-5
):
    """
    Builds visibility graphs connecting some vertices around obstacles
    represented by polygons based on several rotation matrices.

    Parameters
    ----------
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    polygons : list
        Discontinuities turned into polygons. Each polygon is a list of tuple,
        each tuple symbolizing an edge of the polygon through the indices
        of its starting and ending vertices in `vertices`.
    in_polygons : ndarray of shape (n_vertices,)
        Indices of the polygon to which each vertex belongs.
    rotation_matrix : array-like of shape (n_rotations, 3, 3)
        Rotation matrix.
    reduced : bool, default=True
        If True, returns a reduced visibility graph, i.e., a graph without
        the edges that would never be used in a shortest path.
    eps: float, default=1e-5
        Extra length added to a segment to check whether it is bitangent.

    Returns
    -------
    visibility_graphs: dict
        The visibility graphs as a list of dictionaries whose keys are the
        indices of the vertices forming an edge of the graph and values are the
        length of the edge.
    """
    visibility_graphs = nb.typed.List()
    visibility_graphs.append(
        build_visibility_graph(
            vertices,
            polygons,
            in_polygons,
            rotation_matrix[0],
            reduced=reduced,
            eps=eps,
        )
    )
    for i in range(1, len(rotation_matrix)):
        visibility_graphs.append(
            update_visibility_graph(vertices, visibility_graphs[0], rotation_matrix[i])
        )

    return visibility_graphs


@nb.njit()
def is_vertex_inside(vertices, i, polygons):
    """
    Checks whether a vertex is inside one of the polygons.
    """
    for polygon in polygons:
        if (
            len(polygon) > 2
            and is_point_in_polygon(vertices[i], vertices, polygon) == True
        ):
            return True
    return False


@nb.njit()
def add_to_visibility_graph(
    new_vertices,
    visibility_graph,
    vertices,
    polygons,
    in_polygons,
    rotation_matrix,
    reduced=True,
    eps=1e-5,
):
    """
    Adds new vertices to a visibility graph.

    Parameters
    ----------
    new_vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the new vertices to add.
    visibility_graph: dict
        Visibility graph as a dictionary whose keys are the indices of the
        vertices forming an edge of the graph and values are the length of
        the edge.
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    polygons : list
        Discontinuities turned into polygons. Each polygon is a list of tuple,
        each tuple symbolizing an edge of the polygon through the indices
        of its starting and ending vertices in `vertices`.
    in_polygons : ndarray of shape (n_vertices,)
        Indices of the polygon to which each vertex belongs.
    rotation_matrix : array-like of shape (3, 3)
        Rotation matrix.
    reduced : bool, default=True
        If True, returns a reduced visibility graph, i.e., a graph without
        the edges that would never be used in a shortest path.
    eps: float, default=1e-5
        Extra length added to a segment to check whether it is bitangent.

    Returns
    -------
    _visibility_graph: dict
        The visibility graph as a dictionary whose keys are the indices of
        the vertices forming an edge of the graph and values are the length
        of the edge.
    vertices : ndarray of shape (n_vertices, 3)
        The coordinates of the vertices of the graph, including the new ones.
    """
    vertices = np.concatenate((vertices, new_vertices))
    in_polygons = np.concatenate(
        (in_polygons, -1 * np.ones(len(new_vertices), dtype=nb.int64))
    )

    # TODO: Don't recompute the vertices from the original graph?
    is_reflex_vertices = None
    if reduced == True:
        is_reflex_vertices = np.empty(len(vertices), dtype=nb.int8)
        for i in range(len(vertices)):
            is_reflex_vertices[i] = is_vertex_reflex(
                vertices, i, polygons[in_polygons[i]]
            )

    _visibility_graph = dict()
    for key in visibility_graph:
        _visibility_graph[key] = visibility_graph[key]
    for i in range(len(vertices) - 1, len(vertices) - len(new_vertices) - 1, -1):
        if is_vertex_inside(vertices, i, polygons) == False and (
            reduced == False or is_reflex_vertices[i] == 1
        ):
            for j in range(i):
                if reduced == False or is_reflex_vertices[j] == 1:
                    if (
                        are_vertices_visible(vertices, i, j, polygons) == True
                        and is_edge_inside(vertices, i, j, in_polygons, polygons)
                        == False
                        and (
                            reduced == False
                            or is_edge_bitangent(
                                vertices, i, j, in_polygons, polygons, eps
                            )
                            == False
                        )
                    ):
                        _visibility_graph[(i, j)] = anisotropic_distance(
                            vertices[i], vertices[j], rotation_matrix
                        )
                        _visibility_graph[(j, i)] = _visibility_graph[(i, j)]

    return _visibility_graph, vertices


################################################################################
# Shortest path


@nb.njit()
def get_neighbors(vertices, visibility_graph):
    """
    Gets the neighbors of the vertices of a visibility graph.

    Parameters
    ----------
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    visibility_graph: dict
        Visibility graph as a dictionary whose keys are the indices of the
        vertices forming an edge of the graph and values are the length of
        the edge.

    Returns
    -------
    neighbors : list
        The indices in `vertices` of the neighbors for each vertex in
        `vertices`.
    """
    neighbors = nb.typed.List()
    for i in range(len(vertices)):
        neighbors.append(nb.typed.List.empty_list(nb.i8))

    for key in visibility_graph:
        neighbors[key[0]].append(key[1])

    return neighbors


@nb.njit()
def uniform_cost_search(source, target, graph, neighbors):
    """
    Finds the shortest distance in a graph using a uniform cost search.

    Parameters
    ----------
    source : int
        Index of the source vertex.
    target : int
        Index of the target vertex.
    graph: dict
        Graph as a dictionary whose keys are the indices of the vertices
        forming an edge of the graph and values are the length of the edge.
    neighbors : list
        Indices of the neighbors for each vertex.

    Returns
    -------
    distance : float
        The distance between the source and target vertices.
    """
    if len(neighbors[source]) == 0 or len(neighbors[target]) == 0:
        return np.inf

    node = source
    distance = 0.0
    frontier = [(distance, source)]
    visited = [source]
    while node != target and len(frontier) > 0:
        distance, node = heapq.heappop(frontier)
        visited.append(node)
        for neighbor in neighbors[node]:
            if neighbor not in visited:
                heapq.heappush(frontier, (distance + graph[(node, neighbor)], neighbor))

    return distance


################################################################################
# Discontinuous distance based on the shortest path


@nb.njit()
def discontinuous_distance(
    point_1, point_2, vertices, polygons, visibility_graph, in_polygons, rotation_matrix
):
    """
    Computes the discontinuous, anisotropic distance between two points given
    the coordinates of each point, the location of some discontinuities and a
    definition of the anisotropy. The discontinuous distance is based on the
    shortest distance around the discontinuities

    Parameters
    ----------
    point_1 : array-like of shape (3,)
        Coordinates of the first point.
    point_2 : array-like of shape (3,)
        Coordinates of the second point.
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    polygons : list
        Discontinuities turned into polygons. Each polygon is a list of tuple,
        each tuple symbolizing an edge of the polygon through the indices
        of its starting and ending vertices in `vertices`.
    visibility_graph: dict
        Visibility graph as a dictionary whose keys are the indices of the
        vertices forming an edge of the graph and values are the length of
        the edge.
    in_polygons : ndarray of shape (n_vertices,)
        Indices of the polygon to which each vertex belongs.
    rotation_matrix : array-like of shape (3, 3)
        Rotation matrix.

    Returns
    -------
    distance : float
        The discontinuous, anisotropic distance.
    """
    _visibility_graph, _vertices = add_to_visibility_graph(
        np.stack((point_1, point_2)),
        visibility_graph,
        vertices,
        polygons,
        in_polygons,
        rotation_matrix,
    )
    # TODO: Precompute the neighbors for the base graph?
    neighbors = get_neighbors(_vertices, _visibility_graph)

    return uniform_cost_search(
        len(_vertices) - 2, len(_vertices) - 1, _visibility_graph, neighbors
    )


@nb.njit(parallel=True, nogil=True)
def _discontinuous_autodistance_matrix(
    X, vertices, discontinuities, visibility_graph, in_discontinuities, rotation_matrix
):
    """
    Computes a discontinuous, anisotropic distance matrix given the coordinates
    of each points, the location of some discontinuities, and a definition of the
    anisotropy.
    """
    distances = np.zeros((len(rotation_matrix), len(X), len(X)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(j + 1, len(X)):
                distances[k, j, i] = discontinuous_distance(
                    X[j],
                    X[i],
                    vertices,
                    discontinuities,
                    visibility_graph[k],
                    in_discontinuities,
                    rotation_matrix[k],
                )
                distances[k, i, j] = distances[k, j, i]

    return distances


@nb.njit(parallel=True, nogil=True)
def _discontinuous_distance_matrix(
    X,
    Y,
    vertices,
    discontinuities,
    visibility_graph,
    in_discontinuities,
    rotation_matrix,
):
    """
    Computes a discontinuous, anisotropic distance matrix given the coordinates
    of each points, the location of some discontinuities, and a definition of the
    anisotropy.
    """
    distances = np.empty((len(rotation_matrix), len(X), len(Y)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(len(Y)):
                distances[k, j, i] = discontinuous_distance(
                    X[j],
                    Y[i],
                    vertices,
                    discontinuities,
                    visibility_graph[k],
                    in_discontinuities,
                    rotation_matrix[k],
                )

    return distances


################################################################################
# Discontinuous distance based on a penalty


@nb.njit()
def compute_penalty(
    intersection, vertices, discontinuity, j, edge_indices, rotation_matrix
):
    """
    Computes the penalty to go around a discontinuity.
    """
    penalty_1 = 2.0 * anisotropic_distance(
        intersection, vertices[discontinuity[j][0]], rotation_matrix
    )
    k = j - 1
    while k > -1 and k not in edge_indices:
        penalty_1 += 2.0 * anisotropic_distance(
            vertices[discontinuity[k][1]],
            vertices[discontinuity[k][0]],
            rotation_matrix,
        )
        k -= 1

    penalty_2 = 2.0 * anisotropic_distance(
        intersection, vertices[discontinuity[j][1]], rotation_matrix
    )
    l = j + 1
    while l < len(discontinuity) and l not in edge_indices:
        penalty_2 += 2.0 * anisotropic_distance(
            vertices[discontinuity[l][0]],
            vertices[discontinuity[l][1]],
            rotation_matrix,
        )
        l += 1

    return min(penalty_1, penalty_2)


@nb.njit()
def is_intersection_in(intersection, intersections):
    """
    Checks whether an intersection is was already encountered.
    """
    for intersection_2 in intersections:
        exist = True
        for i in range(len(intersection)):
            if isclose(intersection[i], intersection_2[i]) == False:
                exist = False
        if exist == True:
            return True
    return False


@nb.njit()
def get_intersections(point_1, point_2, vertices, discontinuity):
    """
    Computes all the intersections between a segment represented by two
    points and the edges representing a discontinuity.
    """
    intersections = []
    edge_indices = nb.typed.List()
    for i, edge in enumerate(discontinuity):
        intersection = intersect(point_1, point_2, vertices[edge[0]], vertices[edge[1]])
        if (
            intersection is not None
            and is_intersection_in(intersection, intersections) == False
        ):
            intersections.append(intersection)
            edge_indices.append(i)

    return intersections, edge_indices


@nb.njit()
def penalized_discontinuous_distance(
    point_1, point_2, vertices, discontinuities, rotation_matrix
):
    """
    Computes the discontinuous, anisotropic distance between two points given
    the coordinates of each point, the location of some discontinuities and a
    definition of the anisotropy. The discontinuous distance is based on a
    the distance to go around the discontinuities.

    Parameters
    ----------
    point_1 : array-like of shape (3,)
        Coordinates of the first point.
    point_2 : array-like of shape (3,)
        Coordinates of the second point.
    vertices : ndarray of shape (n_vertices, 3)
        Coordinates of the vertices of the graph.
    discontinuities : list of ndarrays
        Set of discontinuities each represented by a list of tuple, each
        tuple symbolizing an edge of the polygon through the indices of its
        starting and ending vertices in `vertices`.
    rotation_matrix : array-like of shape (3, 3)
        Rotation matrix.

    Returns
    -------
    distance : float
        The discontinuous, anisotropic distance.

    References
    ----------
    Pouzet , J. (1980)
        Estimation of a surface with known discontinuities for automatic contouring purposes
        https://doi.org/10.1007/BF01034744
    """
    distance = anisotropic_distance(point_1, point_2, rotation_matrix)
    for discontinuity in discontinuities:

        intersections, edge_indices = get_intersections(
            point_1, point_2, vertices, discontinuity
        )
        if len(intersections) > 0:
            for i, j in enumerate(edge_indices):
                distance += compute_penalty(
                    intersections[i],
                    vertices,
                    discontinuity,
                    j,
                    edge_indices,
                    rotation_matrix,
                )

    return distance


@nb.njit(parallel=True, nogil=True)
def _penalized_discontinuous_autodistance_matrix(
    X, vertices, discontinuities, rotation_matrix
):
    """
    Computes a discontinuous, anisotropic distance matrix given the coordinates
    of each points, the location of some discontinuities, and a definition of the
    anisotropy.
    """
    distances = np.zeros((len(rotation_matrix), len(X), len(X)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(j + 1, len(X)):
                distances[k, j, i] = penalized_discontinuous_distance(
                    X[j], X[i], vertices, discontinuities, rotation_matrix[k]
                )
                distances[k, i, j] = distances[k, j, i]

    return distances


@nb.njit(parallel=True, nogil=True)
def _penalized_discontinuous_distance_matrix(
    X, Y, vertices, discontinuities, rotation_matrix
):
    """
    Computes a discontinuous, anisotropic distance matrix given the coordinates
    of each points, the location of some discontinuities, and a definition of the
    anisotropy.
    """
    distances = np.empty((len(rotation_matrix), len(X), len(Y)))
    for k in range(len(rotation_matrix)):
        for j in nb.prange(len(X)):
            for i in range(len(Y)):
                distances[k, j, i] = penalized_discontinuous_distance(
                    X[j], Y[i], vertices, discontinuities, rotation_matrix[k]
                )

    return distances


################################################################################
# Distance


def distance_matrix(
    X,
    Y=None,
    rotation_matrix=None,
    vertices=None,
    discontinuities=None,
    visibility_graph=None,
    in_discontinuities=None,
):
    """
    Computes a (discontinuous) (anisotropic) distance matrix given the coordinates
    of each points, the location of some discontinuities, and a definition of the
    anisotropy. If none of the parameters after `Y` are provided, computes the
    Euclidean distance. If `rotation_matrix` is provided, computes the
    anisotropic distance. If `vertices` and `discontinuities` are provided,
    computes the discontinuous distance using a penalty. If `vertices`,
    `discontinuities`, `visibility_graph`, and `in_discontinuities` are provided,
    computes the discontinuous distance using the shortest path.

    Parameters
    ----------
    X : array-like of shape (n_points, 2) or (n_points, 3)
        Coordinates of the data.
    Y : array-like of shape (n_points, 2) or (n_points, 3), default=None
        Coordinates of the data. If None, computes the distances between the
        data of `X`.
    vertices : ndarray of shape (n_vertices, 2) or (n_vertices, 3), default=None
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
    rotation_matrix : array-like of shape (n_rotations, 2, 2) or (n_rotations, 3, 3), default=None
        Rotation matrix.

    Returns
    -------
    distances : ndarray of shape (n_rotations, n_points)
        Squared anisotropic distance matrix.
    """
    if rotation_matrix is None:
        if X.shape[1] == 2:
            rotation_matrix = set_rotation_matrix_2d()
        else:
            rotation_matrix = set_rotation_matrix()
    if vertices is None or discontinuities is None:
        if Y is None:
            return _anisotropic_autodistance_matrix(X, rotation_matrix)
        else:
            return _anisotropic_distance_matrix(X, Y, rotation_matrix)
    else:
        if visibility_graph is not None and in_discontinuities is not None:
            if Y is None:
                return _discontinuous_autodistance_matrix(
                    X,
                    vertices,
                    discontinuities,
                    visibility_graph,
                    in_discontinuities,
                    rotation_matrix,
                )
            else:
                return _discontinuous_distance_matrix(
                    X,
                    Y,
                    vertices,
                    discontinuities,
                    visibility_graph,
                    in_discontinuities,
                    rotation_matrix,
                )
        else:
            if Y is None:
                return _penalized_discontinuous_autodistance_matrix(
                    X, vertices, discontinuities, rotation_matrix
                )
            else:
                return _penalized_discontinuous_distance_matrix(
                    X, Y, vertices, discontinuities, rotation_matrix
                )
