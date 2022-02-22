################################################################################
# @file process.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import logging
import numpy as np
from pymeshfv3d.utilities import Mesh, ElementType, Cell, Boundary, Region
from pymeshfv3d.utilities import MetisLibrary
from pymeshfv3d.utilities import calc_area, calc_distance, calc_normal
from pymeshfv3d.utilities import calc_t1, calc_t2, calc_volume, calc_weight


class LoopContinue(Exception):
    """Outer loop continue exception."""
    pass


class ProcessMesh(Mesh):
    """Process mesh class."""

    def __init__(self, mesh, scale, n_partitions):
        super().__init__(mesh.name, mesh.elements, mesh.regions, mesh.vertices)
        self.scale = scale
        self.n_partitions = n_partitions

        self._set_regions()
        self._element_to_cell_boundaries()
        self._generate_cell_faces()
        self._link_boundaries()

        self._calc_geometry()

        if self.n_partitions >= 2:
            with MetisLibrary() as metis:
                metis.metis_part_mesh_nodal(self)

    def __str__(self):
        return '"{:}" (nb={:}, nc={:}, nf={:}, nv={:})'.format(
            self.name, self.n_boundaries, self.n_cells,
            self.n_faces, self.n_vertices)

    def _calc_geometry(self):
        self.vertices[:, 0] *= self.scale[0]
        self.vertices[:, 1] *= self.scale[1]
        self.vertices[:, 2] *= self.scale[2]

        for cell in self.cells:
            cell.x_v = np.mean(self.vertices[cell.vertice_ids], axis=0)
            cell.volume = calc_volume(self.vertices, cell)

        for face in self.faces:
            face.x_v = np.mean(self.vertices[face.vertice_ids], axis=0)
            face.area = calc_area(self.vertices, face)
            face.n_v = calc_normal(self.vertices, self.cells, face)
            face.t1_v = calc_t1(face)
            face.t2_v = calc_t2(face)
            face.weight = calc_weight(self.cells, face)

        for boundary in self.boundaries:
            boundary.n_v = -self.faces[boundary.face_id].n_v
            boundary.t1_v = calc_t1(boundary)
            boundary.t2_v = calc_t2(boundary)
            boundary.distance = calc_distance(self.cells, self.faces, boundary)

            region = self.regions[boundary.region_id]
            region.is_boundary = True

        for cell in self.cells:
            cell.dx = np.zeros(3)
            for idx in cell.face_ids:
                face = self.faces[idx]
                cell.dx += np.abs(face.n_v) * face.area

    def _element_to_cell_boundaries(self):
        face_type = self._get_face_type()

        self.cells = []
        self.boundaries = []
        for element in self.elements:
            if element.element_type.value > face_type.value:
                self.cells.append(
                    Cell(element.element_type, element.vertice_ids,
                         element.region_id)
                )
            else:
                self.boundaries.append(
                    Boundary(element.element_type, element.vertice_ids,
                             element.region_id)
                )

    def _generate_cell_faces(self):
        tmp_faces = []
        for i, cell in enumerate(self.cells):
            for face in cell.get_faces_elements():
                face.cell_ids[0] = i
                tmp_faces.append(face)

        checked = [False] * len(tmp_faces)
        vertice_faces = generate_vertice_faces(self.vertices, tmp_faces)

        self.faces = []
        for i, face in enumerate(tmp_faces):
            if checked[i]:
                continue
            i_sum = sum(face.vertice_ids)

            try:
                for vid in face.vertice_ids:
                    for j in vertice_faces[vid]:
                        if i == j:
                            continue

                        j_face = tmp_faces[j]
                        if not compare_faces(face, j_face, i_sum):
                            continue

                        checked[i] = True
                        checked[j] = True

                        face.cell_ids[1] = j_face.cell_ids[0]
                        j_face.cell_ids = None
                        self.faces.append(face)
                        raise LoopContinue
            except LoopContinue:
                continue

            # LoopContinue not raise -> it must be a boundary (only one face)
            if face.cell_ids[1] < 0:
                self.faces.append(face)

        # link cell face IDs
        for i, face in enumerate(self.faces):
            i_cell, j_cell = face.cell_ids
            self.cells[i_cell].face_ids.append(i)
            if not face.is_boundary:
                self.cells[j_cell].face_ids.append(i)
                face.boundary_id = -1

    def _link_boundaries(self):
        vertice_faces = generate_vertice_faces(self.vertices, self.faces)

        for i, boundary in enumerate(self.boundaries):
            i_sum = sum(boundary.vertice_ids)

            try:
                for vid in boundary.vertice_ids:
                    for j in vertice_faces[vid]:

                        j_face = self.faces[j]
                        if not j_face.is_boundary:
                            continue
                        if not compare_faces(boundary, j_face, i_sum):
                            continue

                        boundary.face_id = j
                        j_face.boundary_id = i
                        raise LoopContinue
            except LoopContinue:
                continue

    def _get_face_type(self):
        dim = self.dimension
        if dim == 3:
            return ElementType.QUADRANGLE
        if dim == 2:
            return ElementType.LINE

        return ElementType.POINT

    def _set_regions(self):
        self.regions = [Region(x) for x in self.regions]

    def to_numpy(self):
        """Change mesh lists to numpy arrays."""
        for cell in self.cells:
            cell.vertice_ids = np.array(cell.vertice_ids)
            cell.face_ids = np.array(cell.face_ids)

        for boundary in self.boundaries:
            boundary.vertice_ids = np.array(boundary.vertice_ids)

        for face in self.faces:
            face.vertice_ids = np.array(face.vertice_ids)
            face.cell_ids = np.array(face.cell_ids)

    @property
    def is_partitioned(self):
        """Flag if mesh is partitioned."""
        return self.n_partitions > 1

    @property
    def n_boundaries(self):
        """Return number of mesh boundaries."""
        return len(self.boundaries)

    @property
    def n_cells(self):
        """Return number of mesh cells."""
        return len(self.cells)

    @property
    def n_faces(self):
        """Return number of mesh faces."""
        return len(self.faces)

    @property
    def volume(self):
        """Return mesh volume."""
        return sum(x.volume for x in self.cells)


def compare_faces(f_1, f_2, f_1_sum):
    """Compare two faces."""
    if f_1.element_type != f_2.element_type:
        return False

    if f_1_sum != sum(f_2.vertice_ids):
        return False

    if set(f_1.vertice_ids) != set(f_2.vertice_ids):
        return False

    return True


def generate_vertice_faces(vertices, faces):
    """Return the list of faces belonging to the given vertices."""
    vertice_faces = [[] for _ in vertices]
    for i, face in enumerate(faces):
        for vertice in face.vertice_ids:
            vertice_faces[vertice].append(i)

    return vertice_faces


def process_mesh(mesh, scale, n_partitions):
    """Process the given mesh, apply scale and partition."""
    logging.info('Process mesh %s s=%s p=%s', mesh, scale, n_partitions)
    return ProcessMesh(mesh, scale, n_partitions)
