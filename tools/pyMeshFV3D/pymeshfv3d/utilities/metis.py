################################################################################
# @file metis.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import sys
import os
from ctypes import cdll, POINTER, byref, c_int, c_float, c_double
import numpy as np

c_int_p = POINTER(c_int)
c_float_p = POINTER(c_float)
c_double_p = POINTER(c_double)

lib_ending = {'darwin': 'dylib', 'win32': 'dll'}.get(sys.platform, 'so')

lib_path = os.path.join('/usr/lib/x86_64-linux-gnu')
if not os.path.exists(lib_path):
    lib_path = '/usr/local/Cellar/metis/5.1.0/lib'

lib_path = os.path.join(lib_path, 'libmetis.{:}'.format(lib_ending))


def wrap_vector(ctype, values):
    """Helper function for list wrapping"""
    return (ctype * len(values))(*values)


class MetisLibrary():
    """Class for handling METIS library interaction."""

    def __init__(self):
        self._lib = cdll.LoadLibrary(lib_path)

        self._lib.METIS_PartMeshNodal.argtypes = [c_int_p, c_int_p, c_int_p,
                                                  c_int_p, c_int_p, c_int_p,
                                                  c_int_p, c_float_p, c_int_p,
                                                  c_int_p, c_int_p, c_int_p]
        self._lib.METIS_PartMeshNodal.restype = c_int

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __del__(self):
        pass

    def metis_part_mesh_nodal(self, mesh):
        """Call METIS PartMeshNodal function."""
        n_total_vertices = sum(x.n_vertice_ids for x in mesh.cells)
        eptr, eind = [-1] * (mesh.n_cells+1), [-1] * (n_total_vertices)
        vwgt, vsize = [1] * (mesh.n_vertices), [1] * (mesh.n_vertices)
        tpwgt = [1.0 / mesh.n_partitions] * (mesh.n_partitions)
        epart, npart = (c_int * mesh.n_cells)(), (c_int * mesh.n_vertices)()

        k = 0
        for i in range(mesh.n_cells):
            eptr[i] = k
            for vid in mesh.cells[i].vertice_ids:
                eind[k] = vid
                k += 1

            eptr[i+1] = k

        i_error = self._lib.METIS_PartMeshNodal(
            byref(c_int(mesh.n_cells)), byref(c_int(mesh.n_vertices)),
            wrap_vector(c_int, eptr), wrap_vector(c_int, eind),
            wrap_vector(c_int, vwgt), wrap_vector(c_int, vsize),
            byref(c_int(mesh.n_partitions)), wrap_vector(c_float, tpwgt),
            wrap_vector(c_int, [-1] * 40), byref(c_int(0)),
            epart, npart
        )

        if i_error != 1:
            raise ValueError('Metis', i_error)

        for i, cell in enumerate(mesh.cells):
            cell.partition_id = epart[i]

        for i, boundary in enumerate(mesh.boundaries):
            i_cell = mesh.faces[boundary.face_id].cell_ids[0]
            boundary.partition_id = mesh.cells[i_cell].partition_id

        complete_parallel_data(mesh)

def complete_parallel_data(mesh):
    """Partition mesh based on METIS."""
    mesh.partition_cells = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_boundaries = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_faces = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_sends = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_sends_pid = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_receives = [[] for _ in range(mesh.n_partitions)]
    mesh.partition_receives_pid = [[] for _ in range(mesh.n_partitions)]

    for i, cell in enumerate(mesh.cells):
        mesh.partition_cells[cell.partition_id].append(i)

    for i, boundary in enumerate(mesh.boundaries):
        mesh.partition_boundaries[boundary.partition_id].append(i)

    is_sender = np.zeros((mesh.n_partitions, mesh.n_cells), dtype=np.int)
    is_receiver = np.zeros((mesh.n_partitions, mesh.n_cells), dtype=np.int)

    for i, face in enumerate(mesh.faces):
        i_cell = face.cell_ids[0]
        i_pid = mesh.cells[i_cell].partition_id
        mesh.partition_faces[i_pid].append(i)

        if not face.is_boundary:
            j_cell = face.cell_ids[1]
            j_pid = mesh.cells[j_cell].partition_id
            if i_pid == j_pid:
                continue

            mesh.partition_faces[j_pid].append(i)

            if is_sender[j_pid][i_cell] == 0:
                mesh.partition_sends[i_pid].append(i_cell)
                mesh.partition_sends_pid[i_pid].append(j_pid)
                is_sender[j_pid][i_cell] = 1

            if is_receiver[j_pid][i_cell] == 0:
                mesh.partition_receives[j_pid].append(i_cell)
                mesh.partition_receives_pid[j_pid].append(i_pid)
                is_receiver[j_pid][i_cell] = 1

            if is_sender[i_pid][j_cell] == 0:
                mesh.partition_sends[j_pid].append(j_cell)
                mesh.partition_sends_pid[j_pid].append(i_pid)
                is_sender[i_pid][j_cell] = 1

            if is_receiver[i_pid][j_cell] == 0:
                mesh.partition_receives[i_pid].append(j_cell)
                mesh.partition_receives_pid[i_pid].append(j_pid)
                is_receiver[i_pid][j_cell] = 1