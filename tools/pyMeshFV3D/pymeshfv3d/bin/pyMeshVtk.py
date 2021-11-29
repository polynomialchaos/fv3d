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
        points = h5_root['VERTICES']['x'][:]
        cell_types = h5_root['CELLS']['type'][:]
        cell_n_vertices = h5_root['CELLS']['n_vertices'][:]
        cell_vertices = h5_root['CELLS']['vertices'][:]

        face_types = h5_root['FACES']['type'][:]
        face_n_vertices = h5_root['FACES']['n_vertices'][:]
        face_vertices = h5_root['FACES']['vertices'][:]

        boundary_types = h5_root['BOUNDARIES']['type'][:]
        boundary_n_vertices = h5_root['BOUNDARIES']['n_vertices'][:]
        boundary_vertices = h5_root['BOUNDARIES']['vertices'][:]

    cells = generate_vtk_grid(
        points, cell_types, cell_n_vertices, cell_vertices)
    write_vtk_grid(args.mesh_file + '.cells.vtu', cells)

    faces = generate_vtk_grid(
        points, face_types, face_n_vertices, face_vertices)
    write_vtk_grid(args.mesh_file + '.faces.vtu', faces)

    boundaries = generate_vtk_grid(
        points, boundary_types, boundary_n_vertices, boundary_vertices)
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
