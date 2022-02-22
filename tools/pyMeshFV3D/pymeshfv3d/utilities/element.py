################################################################################
# @file element.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
from enum import Enum, unique
import numpy as np


@unique
class ElementType(Enum):
    """MeshFV3D element types (sorted by dimension)"""
    POINT = 1
    LINE = 2
    TRIANGLE = 3
    QUADRANGLE = 4
    TETRAEDER = 5
    HEXAEDER = 6
    PRISM = 7
    PYRAMID = 8


class Element():
    """Element class."""

    def __init__(self, element_type, vertice_ids, region_id):
        self.element_type = element_type
        self.vertice_ids = list(np.atleast_1d(vertice_ids))
        self.region_id = region_id

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return '{:}'.format(self.__dict__)

    @property
    def n_vertice_ids(self):
        """Number of vertice IDs."""
        return len(self.vertice_ids)


class Cell(Element):
    """Cell class."""

    def __init__(self, element_type, vertice_ids, region_id):
        super().__init__(element_type, vertice_ids, region_id)
        self.face_ids = []
        self.x_v = None
        self.volume = None
        self.dx = None
        self.partition_id = None

    def get_faces_elements(self):
        """Generate list of face elements based on element type."""
        vertices = self.vertice_ids

        tmp_faces = []
        if self.element_type == ElementType.LINE:
            v_0, v_1 = vertices[0], vertices[1]
            tmp_faces.append(Face(ElementType.POINT, v_0))
            tmp_faces.append(Face(ElementType.POINT, v_1))
        elif self.element_type == ElementType.TRIANGLE:
            v_0, v_1 = vertices[0], vertices[1]
            v_2 = vertices[2]
            tmp_faces.append(Face(ElementType.LINE, (v_0, v_1)))
            tmp_faces.append(Face(ElementType.LINE, (v_1, v_2)))
            tmp_faces.append(Face(ElementType.LINE, (v_2, v_0)))
        elif self.element_type == ElementType.QUADRANGLE:
            v_0, v_1 = vertices[0], vertices[1]
            v_2, v_3 = vertices[2], vertices[3]
            tmp_faces.append(Face(ElementType.LINE, (v_0, v_1)))
            tmp_faces.append(Face(ElementType.LINE, (v_1, v_2)))
            tmp_faces.append(Face(ElementType.LINE, (v_2, v_3)))
            tmp_faces.append(Face(ElementType.LINE, (v_3, v_0)))
        elif self.element_type == ElementType.TETRAEDER:
            v_0, v_1 = vertices[0], vertices[1]
            v_2, v_3 = vertices[2], vertices[3]
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_1, v_2)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_3, v_1)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_3, v_2)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_2, v_3, v_1)))
        elif self.element_type == ElementType.HEXAEDER:
            v_0, v_1 = vertices[0], vertices[1]
            v_2, v_3 = vertices[2], vertices[3]
            v_4, v_5 = vertices[4], vertices[5]
            v_6, v_7 = vertices[6], vertices[7]
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_1, v_2, v_3)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_4, v_5, v_6, v_7)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_4, v_7, v_3)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_4, v_5, v_1)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_1, v_5, v_6, v_2)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_3, v_7, v_6, v_2)))
        elif self.element_type == ElementType.PRISM:
            v_0, v_1 = vertices[0], vertices[1]
            v_2, v_3 = vertices[2], vertices[3]
            v_4, v_5 = vertices[4], vertices[5]
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_1, v_4, v_3)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_1, v_2, v_5, v_4)))
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_2, v_5, v_3)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_1, v_2)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_3, v_4, v_5)))
        elif self.element_type == ElementType.PYRAMID:
            v_0, v_1 = vertices[0], vertices[1]
            v_2, v_3 = vertices[2], vertices[3]
            v_4 = vertices[4]
            tmp_faces.append(
                Face(ElementType.QUADRANGLE, (v_0, v_1, v_2, v_3)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_1, v_4)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_1, v_2, v_4)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_2, v_3, v_4)))
            tmp_faces.append(Face(ElementType.TRIANGLE, (v_0, v_3, v_1)))
        else:
            raise TypeError(self.element_type)

        return tmp_faces

    @property
    def n_face_ids(self):
        """Number of face IDs."""
        return len(self.face_ids)


class Boundary(Element):
    """Boundary class."""
    # pylint: disable=too-few-public-methods

    def __init__(self, element_type, vertice_ids, region_id):
        super().__init__(element_type, vertice_ids, region_id)
        self.face_id = None
        self.n_v = None
        self.t1_v = None
        self.t2_v = None
        self.distance = None
        self.partition_id = None


class Face(Element):
    """Face class."""
    # pylint: disable=too-many-instance-attributes

    def __init__(self, element_type, vertice_ids):
        super().__init__(element_type, vertice_ids, None)
        self.cell_ids = [-1, -1]
        self.boundary_id = None
        self.x_v = None
        self.n_v = None
        self.t1_v = None
        self.t2_v = None
        self.area = None
        self.weight = None

    @property
    def is_boundary(self):
        """Return flag if face is boundary (second cell index is below 0)."""
        return self.cell_ids[1] < 0


class Region():
    """Region class."""

    def __init__(self, name, is_boundary=False):
        self.name = name
        self.is_boundary = is_boundary

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return '{:}'.format(self.__dict__)
