################################################################################
# @file geometry.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import numpy as np
from .element import ElementType
from .math import len_vector, norm_vector, cross_product


def calc_area(all_vertices, element):
    """Calculate the area for the given element based on type.
    Returns 1 (1D), length between points (2D) and area within points (3D)."""
    vertices = all_vertices[element.vertice_ids]

    # 1D case with lines (volume) and points (area)
    if element.element_type == ElementType.POINT:
        area = 1.0
    # 2D case with quadrangles or triangles (volume) and lines (area)
    elif element.element_type == ElementType.LINE:
        a_v = vertices[1] - vertices[0]
        area = len_vector(a_v)
    # 3D case with heaxahedrons, prisms, pyramids or tetrahedrons (volume)
    # and quadrangles or triangles (area)
    elif element.element_type == ElementType.QUADRANGLE:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        c_v = vertices[3] - vertices[0]
        area = len_vector(cross_product(a_v, b_v))
        area += len_vector(cross_product(b_v, c_v))
        area *= 0.5
    elif element.element_type == ElementType.TRIANGLE:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        area = len_vector(cross_product(a_v, b_v))
        area *= 0.5
    else:
        raise TypeError(element.element_type)

    if area < 0.0:
        raise ValueError('Area', area)

    return area


def calc_volume(all_vertices, element):
    """Calculate the volume for the given element based on type.
    Returns length between points (2D) and volume within points (3D)."""
    vertices = all_vertices[element.vertice_ids]

    # 1D case with lines (volume) and points (area)
    if element.element_type == ElementType.LINE:
        a_v = vertices[1] - vertices[0]
        volume = len_vector(a_v)
    # 2D case with quadrangles or triangles (volume) and lines (area)
    elif element.element_type == ElementType.QUADRANGLE:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        c_v = vertices[3] - vertices[0]
        volume = len_vector(cross_product(a_v, b_v))
        volume += len_vector(cross_product(b_v, c_v))
        volume *= 0.5
    elif element.element_type == ElementType.TRIANGLE:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        volume = len_vector(cross_product(a_v, b_v))
        volume *= 0.5
    # 3D case with heaxahedrons, prisms, pyramids or tetrahedrons (volume)
    # and quadrangles or triangles (area)
    elif element.element_type == ElementType.HEXAEDER:
        volume = abs(np.dot(vertices[5] - vertices[4], cross_product(
            vertices[1] - vertices[4], vertices[6] - vertices[4])))
        volume += abs(np.dot(vertices[1] - vertices[4], cross_product(
            vertices[3] - vertices[4], vertices[6] - vertices[4])))
        volume += abs(np.dot(vertices[1] - vertices[4], cross_product(
            vertices[0] - vertices[4], vertices[3] - vertices[4])))
        volume += abs(np.dot(vertices[6] - vertices[4], cross_product(
            vertices[3] - vertices[4], vertices[7] - vertices[4])))
        volume += abs(np.dot(vertices[3] - vertices[1], cross_product(
            vertices[6] - vertices[1], vertices[2] - vertices[1])))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.PRISM:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        c_v = vertices[3] - vertices[0]
        d_v = vertices[4] - vertices[0]
        e_v = vertices[5] - vertices[0]
        volume = abs(np.dot(a_v, cross_product(b_v, e_v)))
        volume += abs(np.dot(a_v, cross_product(e_v, d_v)))
        volume += abs(np.dot(d_v, cross_product(e_v, c_v)))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.PYRAMID:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        c_v = vertices[3] - vertices[0]
        d_v = vertices[4] - vertices[0]
        volume = np.abs(np.dot(a_v, cross_product(b_v, d_v)))
        volume += np.abs(np.dot(b_v, cross_product(c_v, d_v)))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.TETRAEDER:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        c_v = vertices[3] - vertices[0]
        volume = np.abs(np.dot(a_v, cross_product(b_v, c_v)))
        volume *= 1.0/6.0
    else:
        raise TypeError(element.element_type)

    if volume < 0.0:
        raise ValueError('Volume', volume)

    return volume


def calc_distance(cells, faces, element):
    """Calculate the cell face distance for the given element."""
    d_x = cells[faces[element.face_id].cell_ids[0]].x_v - \
        faces[element.face_id].x_v
    tmp = np.dot(d_x, element.n_v) / \
        (len_vector(d_x) * len_vector(element.n_v))
    distance = len_vector(d_x) * max(min(tmp, 1.0), -1.0)

    if distance < 1e-16:
        raise ValueError('Distance', distance)

    return distance


def calc_weight(cells, element):
    """Calculate the interpolation weight for the given element.
    Returns 1 for boundary elements (cell[2] >= 0)."""
    if not element.is_boundary:
        d_1 = np.dot(element.n_v, element.x_v -
                     cells[element.cell_ids[0]].x_v) / element.area
        d_2 = np.dot(
            element.n_v, cells[element.cell_ids[1]].x_v - element.x_v) / element.area

        if d_1 + d_2 > 1e-16:
            weight = d_1 / (d_1 + d_2)
        else:
            raise ValueError('Weight', d_1 + d_2)
    else:
        weight = 1.0

    return weight


def calc_normal(all_vertices, cells, element):
    """Calculate the normal vector for the given element based on type."""
    vertices = all_vertices[element.vertice_ids]

    # 1D case with lines (volume) and points (area)
    if element.element_type == ElementType.POINT:
        n_v = norm_vector(vertices[0] - cells[element.cell_ids[0]].x_v)
    # 2D case with quadrangles or triangles (volume) and lines (area)
    elif element.element_type == ElementType.LINE:
        a_v = vertices[1] - vertices[0]
        n_v = norm_vector(np.array([-a_v[1], a_v[0], a_v[2]]))
    # 3D case with heaxahedrons, prisms, pyramids or tetrahedrons (volume)
    # and quadrangles or triangles (area)
    elif element.element_type == ElementType.QUADRANGLE:
        a_v = vertices[2] - vertices[0]
        b_v = vertices[3] - vertices[1]
        n_v = norm_vector(cross_product(a_v, b_v))
    elif element.element_type == ElementType.TRIANGLE:
        a_v = vertices[1] - vertices[0]
        b_v = vertices[2] - vertices[0]
        n_v = norm_vector(cross_product(a_v, b_v))
    else:
        raise TypeError(element.element_type)

    # check for normal vector direction (face/cell distance check)
    d_x = element.x_v - cells[element.cell_ids[0]].x_v
    cos_angle = np.dot(d_x, n_v) / (len_vector(d_x) * len_vector(n_v))
    angle = np.arccos(max(min(cos_angle, 1.0), -1.0)) * 180 / np.pi
    if 90.0 < angle < 270.0:
        n_v = -n_v

    return n_v


def calc_t1(element):
    """Calculate the first tangential vector for the given element
    based on type and normal vector."""
    n_0, n_1, n_2 = element.n_v
    t1_v = np.array([-n_2, 0.0, n_0]) if abs(n_0) > abs(n_1) \
        else np.array([0.0, -n_2, n_1])

    return norm_vector(t1_v)


def calc_t2(element):
    """Calculate the first tangential vector for the given element
    based on type and normal vector."""
    t2_v = cross_product(element.t1_v, element.n_v)
    return norm_vector(t2_v)
