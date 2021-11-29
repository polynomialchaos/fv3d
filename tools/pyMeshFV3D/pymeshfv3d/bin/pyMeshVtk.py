################################################################################
# @file pyMeshVtk.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import sys
import argparse
import logging
from typing import Type
import h5py
from pymeshfv3d.utilities.element import Boundary
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from pymeshfv3d import LessThanFilter, ElementType, reader
from pychemistry.version import __version__

logging.getLogger().setLevel(logging.DEBUG)

err_handler = logging.StreamHandler(stream=sys.stderr)
err_handler.setLevel(logging.ERROR)
err_handler.setFormatter(logging.Formatter(
    '%(levelname)s:%(name)s:%(message)s'))
logging.getLogger().addHandler(err_handler)

out_handler = logging.StreamHandler(stream=sys.stdout)
out_handler.setLevel(logging.INFO)
out_handler.addFilter(LessThanFilter(logging.ERROR))
out_handler.setFormatter(logging.Formatter(
    '%(levelname)s:%(name)s:%(message)s'))
logging.getLogger().addHandler(out_handler)


vtk_constructor = [None] * 9
vtk_constructor[ElementType.POINT.value] = vtk.vtkVertex()
vtk_constructor[ElementType.LINE.value] = vtk.vtkLine()
vtk_constructor[ElementType.QUADRANGLE.value] = vtk.vtkQuad()
vtk_constructor[ElementType.TRIANGLE.value] = vtk.vtkTriangle()
vtk_constructor[ElementType.HEXAEDER.value] = vtk.vtkHexahedron()
vtk_constructor[ElementType.PRISM.value] = vtk.vtkWedge()
vtk_constructor[ElementType.PYRAMID.value] = vtk.vtkPyramid()
vtk_constructor[ElementType.TETRAEDER.value] = vtk.vtkTetra()


def main():
    """Main function entrance point."""

    # define the argument parser
    parser = argparse.ArgumentParser(
        description='pyMeshVtk - Mesh to Vtk converter')
    parser.add_argument('-d', '--debugging',
                        action='store_true', help='Debugging output')
    parser.add_argument('-m', '--mesh-file',
                        required=True, type=str, help='The mesh file')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s (version {:})'.format(__version__))
    # parser.add_argument('inifile', nargs='?',
    #                     help='The preprocessing input file')
    args = parser.parse_args()

    with h5py.File(args.mesh_file, 'r') as h5_root:
        date = h5_root.attrs['date']
        version = h5_root.attrs['version']
        dimension = h5_root.attrs['dimension']
        is_partitioned = h5_root.attrs['is_partitioned']

        # vertice data
        points = h5_root['VERTICES']['x'][:]

        # cell data
        cell_ids = h5_root['CELLS']['id'][:]
        cell_types = h5_root['CELLS']['type'][:]
        cell_n_vertices = h5_root['CELLS']['n_vertices'][:]
        cell_vertices = h5_root['CELLS']['vertices'][:]
        # cell_n_faces = h5_root['CELLS']['n_faces'][:]
        # cell_faces = h5_root['CELLS']['faces'][:]
        cell_volumes = h5_root['CELLS']['volume'][:]
        cell_xs = h5_root['CELLS']['x'][:]
        cell_dxs = h5_root['CELLS']['dx'][:]
        if is_partitioned:
            cell_pids = h5_root['CELLS']['pid'][:]

        # face data
        face_types = h5_root['FACES']['type'][:]
        face_cells = h5_root['FACES']['cells'][:]
        face_boundaries = h5_root['FACES']['boundary'][:]
        face_n_vertices = h5_root['FACES']['n_vertices'][:]
        face_vertices = h5_root['FACES']['vertices'][:]
        face_xs = h5_root['FACES']['x'][:]
        face_ns = h5_root['FACES']['n'][:]
        face_t1s = h5_root['FACES']['t1'][:]
        face_t2s = h5_root['FACES']['t2'][:]
        face_areas = h5_root['FACES']['area'][:]
        face_lambdas = h5_root['FACES']['lambda'][:]

        # boundary data
        boundary_ids = h5_root['BOUNDARIES']['id'][:]
        boundary_types = h5_root['BOUNDARIES']['type'][:]
        boundary_n_vertices = h5_root['BOUNDARIES']['n_vertices'][:]
        boundary_vertices = h5_root['BOUNDARIES']['vertices'][:]
        boundary_faces = h5_root['BOUNDARIES']['face'][:]
        boundary_ns = h5_root['BOUNDARIES']['n'][:]
        boundary_t1s = h5_root['BOUNDARIES']['t1'][:]
        boundary_t2s = h5_root['BOUNDARIES']['t2'][:]
        boundary_distances = h5_root['BOUNDARIES']['distance'][:]
        if is_partitioned:
            boundary_pids = h5_root['BOUNDARIES']['pid'][:]

    cells = generate_vtk_grid(
        points, cell_types, cell_n_vertices, cell_vertices)
    set_vtk_celldata(cells, 'id', cell_ids, vtk.VTK_INT)
    set_vtk_celldata(cells, 'volume', cell_volumes, vtk.VTK_DOUBLE)
    set_vtk_celldata(cells, 'x', cell_xs, vtk.VTK_DOUBLE)
    set_vtk_celldata(cells, 'dx', cell_dxs, vtk.VTK_DOUBLE)
    if is_partitioned:
        set_vtk_celldata(cells, 'pid', cell_pids, vtk.VTK_INT)
    write_vtk_grid(args.mesh_file + '.cells.vtu', cells)

    faces = generate_vtk_grid(
        points, face_types, face_n_vertices, face_vertices)
    set_vtk_celldata(faces, 'cells', face_cells, vtk.VTK_INT)
    set_vtk_celldata(faces, 'boundary', face_boundaries, vtk.VTK_INT)
    set_vtk_celldata(faces, 'x', face_xs, vtk.VTK_DOUBLE)
    set_vtk_celldata(faces, 'n', face_ns, vtk.VTK_DOUBLE)
    set_vtk_celldata(faces, 't1', face_t1s, vtk.VTK_DOUBLE)
    set_vtk_celldata(faces, 't2', face_t2s, vtk.VTK_DOUBLE)
    set_vtk_celldata(faces, 'area', face_areas, vtk.VTK_DOUBLE)
    set_vtk_celldata(faces, 'lambda', face_lambdas, vtk.VTK_DOUBLE)
    write_vtk_grid(args.mesh_file + '.faces.vtu', faces)

    boundaries = generate_vtk_grid(
        points, boundary_types, boundary_n_vertices, boundary_vertices)
    set_vtk_celldata(boundaries, 'id', boundary_ids, vtk.VTK_INT)
    set_vtk_celldata(boundaries, 'face', boundary_faces, vtk.VTK_INT)
    set_vtk_celldata(boundaries, 'n', boundary_ns, vtk.VTK_DOUBLE)
    set_vtk_celldata(boundaries, 't1', boundary_t1s, vtk.VTK_DOUBLE)
    set_vtk_celldata(boundaries, 't2', boundary_t2s, vtk.VTK_DOUBLE)
    set_vtk_celldata(boundaries, 'distance',
                     boundary_distances, vtk.VTK_DOUBLE)
    if is_partitioned:
        set_vtk_celldata(boundaries, 'pid', boundary_pids, vtk.VTK_INT)
    write_vtk_grid(args.mesh_file + '.boundaries.vtu', boundaries)


def generate_vtk_grid(points, element_types,
                      n_element_vertics, element_vertics):
    """Generate a unstructured grid for the given points and elements."""
    # define the points
    vtk_points = vtk.vtkPoints()
    for point in points:
        vtk_points.InsertNextPoint(point)

    # define the unstructured grid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(vtk_points)

    # append the cells
    for e_t, e_nv, e_v in zip(
            element_types, n_element_vertics, element_vertics):

        cell = vtk_constructor[e_t]
        cell_pids = cell.GetPointIds()
        for i in range(e_nv):
            cell_pids.SetId(i, e_v[i])

        grid.InsertNextCell(cell.GetCellType(), cell_pids)

    return grid


def set_vtk_celldata(grid, key, data, vtk_type):
    """Append cell data to the grid."""
    vtkArray = numpy_to_vtk(data, array_type=vtk_type)
    vtkArray.SetName(key)

    grid.GetCellData().AddArray(vtkArray)


def write_vtk_grid(path, grid):
    """Write the grid to the given path."""
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(grid)
    writer.Write()


################################################################################
# CALL BY SCRIPT
# ------------------------------------------------------------------------------
if __name__ == "__main__":

    main()
