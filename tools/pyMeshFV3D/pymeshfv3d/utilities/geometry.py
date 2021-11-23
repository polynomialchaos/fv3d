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


def calc_area(vertices, element):
    v = vertices[element.vertice_ids]

    if element.element_type == ElementType.POINT:
        area = 1.0
    elif element.element_type == ElementType.LINE:
        area = len_vector(v[1] - v[0])
    elif element.element_type == ElementType.TRIANGLE:
        area = len_vector(cross_product(v[1] - v[0], v[2] - v[0]))
        area *= 0.5
    elif element.element_type == ElementType.QUADRANGLE:
        area = len_vector(cross_product(v[1] - v[0], v[2] - v[0]))
        area += len_vector(cross_product(v[2] - v[0], v[3] - v[0]))
        area *= 0.5
    else:
        raise TypeError(element.element_type)

    if area < 0.0:
        raise(ValueError('Area', area))

    return area


def calc_volume(vertices, element):
    v = vertices[element.vertice_ids]

    if element.element_type == ElementType.LINE:
        volume = len_vector(v[1] - v[0])
    elif element.element_type == ElementType.TRIANGLE:
        volume = len_vector(cross_product(v[1] - v[0], v[2] - v[0]))
        volume *= 0.5
    elif element.element_type == ElementType.QUADRANGLE:
        volume = len_vector(cross_product(v[1] - v[0], v[2] - v[0]))
        volume += len_vector(cross_product(v[2] - v[0], v[3] - v[0]))
        volume *= 0.5
    elif element.element_type == ElementType.TETRAEDER:
        volume = abs(
            np.dot(v[1] - v[0], cross_product(v[2] - v[0], v[3] - v[0])))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.HEXAEDER:
        volume = abs(
            np.dot(v[5] - v[4], cross_product(v[1] - v[4], v[6] - v[4])))
        volume += abs(np.dot(v[1] - v[4],
                             cross_product(v[3] - v[4], v[6] - v[4])))
        volume += abs(np.dot(v[1] - v[4],
                             cross_product(v[0] - v[4], v[3] - v[4])))
        volume += abs(np.dot(v[6] - v[4],
                             cross_product(v[3] - v[4], v[7] - v[4])))
        volume += abs(np.dot(v[3] - v[1],
                             cross_product(v[6] - v[1], v[2] - v[1])))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.PRISM:
        volume = abs(
            np.dot(v[1] - v[0], cross_product(v[2] - v[0], v[5] - v[0])))
        volume += abs(np.dot(v[1] - v[0],
                             cross_product(v[5] - v[0], v[4] - v[0])))
        volume += abs(np.dot(v[4] - v[0],
                             cross_product(v[5] - v[0], v[3] - v[0])))
        volume *= 1.0/6.0
    elif element.element_type == ElementType.PYRAMID:
        volume = abs(
            np.dot(v[1] - v[0], cross_product(v[2] - v[0], v[4] - v[0])))
        volume += abs(np.dot(v[2] - v[0],
                             cross_product(v[3] - v[0], v[4] - v[0])))
        volume *= 1.0/6.0
    else:
        raise TypeError(element.element_type)

    if volume < 0.0:
        raise(ValueError('Volume', volume))

    return volume


def calc_distance(cells, faces, element):
    dx = cells[faces[element.face_id].cell_ids[0]].x - faces[element.face_id].x
    tmp = np.dot(dx, element.n) / (len_vector(dx) * len_vector(element.n))
    distance = len_vector(dx) * max(min(tmp, 1.0), -1.0)
    if distance < 1e-16:
        raise(ValueError('Distance', distance))
    return distance


def calc_weight(cells, element):
    if element.cell_ids[1] >= 0:
        d1 = np.dot(element.n, element.x -
                    cells[element.cell_ids[0]].x) / element.area
        d2 = np.dot(
            element.n, cells[element.cell_ids[1]].x - element.x) / element.area

        if d1 + d2 > 1e-16:
            weight = d1 / (d1 + d2)
        else:
            raise(ValueError('Weight', d1 + d2))
    else:
        weight = 1.0

    return weight


def calc_normal(vertices, cells, element):
    v = vertices[element.vertice_ids]

    if element.element_type == ElementType.POINT:
        n = norm_vector(v[0] - cells[element.cell_ids[0]].x)
    elif element.element_type == ElementType.LINE:
        a = v[1] - v[0]
        n = norm_vector(np.array([-a[1], a[0], a[2]]))
    elif element.element_type == ElementType.TRIANGLE:
        n = norm_vector(cross_product(v[1] - v[0], v[2] - v[0]))
    elif element.element_type == ElementType.QUADRANGLE:
        n = norm_vector(cross_product(v[2] - v[0], v[3] - v[1]))
    else:
        raise TypeError(element.element_type)

    dx = element.x - cells[element.cell_ids[0]].x
    tmp = np.dot(dx, n) / (len_vector(dx) * len_vector(n))
    angle = np.arccos(max(min(tmp, 1.0), -1.0)) * 180 / np.pi
    if angle > 90.0 and angle < 270.0:
        n = -n

    return n


def calc_t1(element):
    n = element.n
    t1 = np.array([-n[2], 0.0, n[0]] if abs(n[0]) >
                  abs(n[1]) else [0.0, -n[2], n[1]])
    return norm_vector(t1)


def calc_t2(element):
    t2 = cross_product(element.t1, element.n)
    return norm_vector(t2)
