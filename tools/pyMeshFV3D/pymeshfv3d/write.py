################################################################################
# @file write.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging
import datetime
import numpy as np
import h5py
from pymeshfv3d.utilities import ElementType

from .version import __version__

type_2_xdmf = [None] * (max([x.value for x in ElementType]) + 1)
type_2_xdmf[ElementType.POINT.value] = 1
type_2_xdmf[ElementType.LINE.value] = 2
type_2_xdmf[ElementType.TRIANGLE.value] = 4
type_2_xdmf[ElementType.QUADRANGLE.value] = 5
type_2_xdmf[ElementType.TETRAEDER.value] = 6
type_2_xdmf[ElementType.HEXAEDER.value] = 9
type_2_xdmf[ElementType.PRISM.value] = 8
type_2_xdmf[ElementType.PYRAMID.value] = 7


def write_mesh(file_name, mesh):
    """Write the mesh."""
    logging.info('Write mesh %s -> %s', mesh, file_name)
    mesh.to_numpy()

    with h5py.File(file_name, 'w') as h5_root:
        h5_root.attrs['date'] = np.string_(str(datetime.datetime.now()))
        h5_root.attrs['version'] = np.string_(str(__version__))
        h5_root.attrs['dimension'] = [mesh.dimension]
        h5_root.attrs['is_partitioned'] = [int(mesh.is_partitioned)]

        write_mesh_vertices(h5_root, mesh)
        write_mesh_cells(h5_root, mesh)
        write_mesh_faces(h5_root, mesh)
        write_mesh_boundaries(h5_root, mesh)
        write_mesh_regions(h5_root, mesh)
        if mesh.is_partitioned:
            write_mesh_partitions(h5_root, mesh)


def write_mesh_boundaries(h5_root, mesh):
    """Write the mesh boundaries."""
    max_boundary_vertices = max(x.n_vertice_ids for x in mesh.boundaries)
    max_boundary_vertices = max_boundary_vertices + 1
    if mesh.dimension < 3:
        max_boundary_vertices += 1

    fp_g = h5_root.create_group('BOUNDARIES')
    fp_g.attrs['n_boundaries'] = [len(mesh.boundaries)]
    fp_g.attrs['max_boundary_vertices'] = [max_boundary_vertices]

    fp_g.create_dataset('id', data=[x.region_id for x in mesh.boundaries])
    fp_g.create_dataset(
        'type', data=[x.element_type.value for x in mesh.boundaries])

    fp_g.create_dataset('n_vertices', data=[
                        x.n_vertice_ids for x in mesh.boundaries])
    fp_g.create_dataset('vertices', data=fill_array(
        [x.vertice_ids for x in mesh.boundaries], -1,
        np.int, size_y=max_boundary_vertices))
    fp_g.create_dataset('face', data=[x.face_id for x in mesh.boundaries])

    fp_g.create_dataset('n', data=np.array([x.n_v for x in mesh.boundaries]))
    fp_g.create_dataset('t1', data=np.array([x.t1_v for x in mesh.boundaries]))
    fp_g.create_dataset('t2', data=np.array([x.t2_v for x in mesh.boundaries]))
    fp_g.create_dataset('distance', data=np.array(
        [x.distance for x in mesh.boundaries]))

    if mesh.is_partitioned:
        fp_g.create_dataset(
            'pid', data=[x.partition_id for x in mesh.boundaries])

    fp_g.create_dataset('xdmf', data=gen_xdmf_2d(
        mesh.boundaries, max_boundary_vertices, mesh.dimension))


def write_mesh_cells(h5_root, mesh):
    """Write the mesh boundaries."""
    max_cell_vertices = max(x.n_vertice_ids for x in mesh.cells)
    max_cell_vertices = max_cell_vertices + 1
    if mesh.dimension < 2:
        max_cell_vertices += 1

    max_cell_faces = max(x.n_face_ids for x in mesh.cells)

    fp_g = h5_root.create_group('CELLS')
    fp_g.attrs['n_cells'] = [mesh.n_cells]
    fp_g.attrs['max_cell_vertices'] = [max_cell_vertices]
    fp_g.attrs['max_cell_faces'] = [max_cell_faces]

    fp_g.create_dataset('id', data=[x.region_id for x in mesh.cells])
    fp_g.create_dataset(
        'type', data=[x.element_type.value for x in mesh.cells])

    fp_g.create_dataset('n_vertices', data=[
                        x.n_vertice_ids for x in mesh.cells])
    fp_g.create_dataset('vertices', data=fill_array(
        [x.vertice_ids for x in mesh.cells], -1,
        np.int, size_y=max_cell_vertices))
    fp_g.create_dataset('n_faces', data=[x.n_face_ids for x in mesh.cells])
    fp_g.create_dataset('faces', data=fill_array(
        [x.face_ids for x in mesh.cells], -1, np.int, size_y=max_cell_faces))

    fp_g.create_dataset('volume', data=np.array(
        [x.volume for x in mesh.cells]))
    fp_g.create_dataset('x', data=np.array([x.x_v for x in mesh.cells]))
    fp_g.create_dataset('dx', data=np.array([x.dx for x in mesh.cells]))

    if mesh.is_partitioned:
        fp_g.create_dataset(
            'pid', data=[x.partition_id for x in mesh.cells])

    fp_g.create_dataset('xdmf', data=gen_xdmf_3d(
        mesh.cells, max_cell_vertices, mesh.dimension))


def write_mesh_faces(h5_root, mesh):
    """Write the mesh faces."""
    max_face_vertices = max(x.n_vertice_ids for x in mesh.faces)
    max_face_vertices = max_face_vertices + 1
    if mesh.dimension < 3:
        max_face_vertices += 1

    fp_g = h5_root.create_group('FACES')
    fp_g.attrs['n_faces'] = [len(mesh.faces)]
    fp_g.attrs['max_face_vertices'] = [max_face_vertices]

    fp_g.create_dataset(
        'type', data=[x.element_type.value for x in mesh.faces])

    fp_g.create_dataset('cells', data=[x.cell_ids for x in mesh.faces])
    fp_g.create_dataset(
        'boundary', data=[x.boundary_id for x in mesh.faces])

    fp_g.create_dataset('n_vertices', data=[
                        x.n_vertice_ids for x in mesh.faces])
    fp_g.create_dataset('vertices', data=fill_array(
        [x.vertice_ids for x in mesh.faces], -1,
        np.int, size_y=max_face_vertices))

    fp_g.create_dataset('x', data=np.array([x.x_v for x in mesh.faces]))
    fp_g.create_dataset('n', data=np.array([x.n_v for x in mesh.faces]))
    fp_g.create_dataset('t1', data=np.array([x.t1_v for x in mesh.faces]))
    fp_g.create_dataset('t2', data=np.array([x.t2_v for x in mesh.faces]))
    fp_g.create_dataset('area', data=np.array(
        [x.area for x in mesh.faces]))
    fp_g.create_dataset('lambda', data=np.array(
        [x.weight for x in mesh.faces]))

    fp_g.create_dataset('xdmf', data=gen_xdmf_2d(
        mesh.faces, max_face_vertices, mesh.dimension))


def write_mesh_partitions(h5_root, mesh):
    """Write the mesh boundaries."""
    n_partition_cells = [len(x) for x in mesh.partition_cells]
    n_partition_boundaries = [len(x) for x in mesh.partition_boundaries]
    n_partition_faces = [len(x) for x in mesh.partition_faces]
    n_partition_sends = [len(x) for x in mesh.partition_sends]
    n_partition_receives = [len(x) for x in mesh.partition_receives]

    fp_g = h5_root.create_group('PARTITIONS')
    fp_g.attrs['n_partitions'] = [mesh.n_partitions]
    fp_g.attrs['max_partition_cells'] = [max(n_partition_cells)]
    fp_g.attrs['max_partition_boundaries'] = [
        max(n_partition_boundaries)]
    fp_g.attrs['max_partition_faces'] = [max(n_partition_faces)]
    fp_g.attrs['max_partition_sends'] = [max(n_partition_sends)]
    fp_g.attrs['max_partition_receives'] = [max(n_partition_receives)]

    fp_g.create_dataset('n_partition_cells', data=n_partition_cells)
    fp_g.create_dataset('n_partition_boundaries',
                        data=n_partition_boundaries)
    fp_g.create_dataset('n_partition_faces', data=n_partition_faces)
    fp_g.create_dataset('n_partition_sends', data=n_partition_sends)
    fp_g.create_dataset('n_partition_receives',
                        data=n_partition_receives)

    fp_g.create_dataset('partition_cells', data=fill_array(
        mesh.partition_cells, -1, np.int))
    fp_g.create_dataset('partition_boundaries', data=fill_array(
        mesh.partition_boundaries, -1, np.int))
    fp_g.create_dataset('partition_faces', data=fill_array(
        mesh.partition_faces, -1, np.int))
    fp_g.create_dataset('partition_sends', data=fill_array(
        mesh.partition_sends, -1, np.int))
    fp_g.create_dataset('partition_sends_pid', data=fill_array(
        mesh.partition_sends_pid, -1, np.int))
    fp_g.create_dataset('partition_receives', data=fill_array(
        mesh.partition_receives, -1, np.int))
    fp_g.create_dataset('partition_receives_pid', data=fill_array(
        mesh.partition_receives_pid, -1, np.int))


def write_mesh_regions(h5_root, mesh):
    """Write the mesh regions."""
    fp_g = h5_root.create_group('REGIONS')
    fp_g.attrs['n_regions'] = [len(mesh.regions)]

    fp_g.create_dataset(
        'name', data=[np.string_(x.name) for x in mesh.regions])
    fp_g.create_dataset('is_boundary', data=[
                        int(x.is_boundary) for x in mesh.regions])


def write_mesh_vertices(h5_root, mesh):
    """Write the mesh regions."""
    fp_g = h5_root.create_group('VERTICES')
    fp_g.attrs['n_vertices'] = [mesh.n_vertices]

    fp_g.create_dataset('x', data=mesh.vertices)


def fill_array(datas, fill, dtype, size_y=None):
    """Generate a numpy array with fill values."""
    size_x = len(datas)
    size_y = max(len(x) for x in datas) if size_y is None else size_y

    tmp = np.full((size_x, size_y), fill, dtype=dtype)
    for i, data in enumerate(datas):
        tmp[i, :len(data)] = data

    return tmp


def gen_xdmf_3d(datas, size_y, dim):
    """Generate a XDMF output for cell elements."""
    size_x = len(datas)
    tmp = np.zeros((size_x, size_y), dtype=np.int)
    for i, data in enumerate(datas):
        if dim >= 2:
            tmp[i, 0] = type_2_xdmf[data.element_type.value]
            tmp[i, 1:1+len(data.vertice_ids)] = data.vertice_ids
        else:
            tmp[i, 0] = type_2_xdmf[data.element_type.value]
            tmp[i, 1] = 2
            tmp[i, 2:2+len(data.vertice_ids)] = data.vertice_ids

    return tmp


def gen_xdmf_2d(datas, size_y, dim):
    """Generate a XDMF output for face elements."""
    size_x = len(datas)
    tmp = np.zeros((size_x, size_y), dtype=np.int)
    for i, data in enumerate(datas):
        if dim > 2:
            tmp[i, 0] = type_2_xdmf[data.element_type.value]
            tmp[i, 1:1+len(data.vertice_ids)] = data.vertice_ids
        elif dim == 2:
            tmp[i, 0] = type_2_xdmf[data.element_type.value]
            tmp[i, 1] = 2
            tmp[i, 2:2+len(data.vertice_ids)] = data.vertice_ids
        else:
            tmp[i, 0] = type_2_xdmf[data.element_type.value]
            tmp[i, 1] = 1
            tmp[i, 2:2+len(data.vertice_ids)] = data.vertice_ids

    return tmp
