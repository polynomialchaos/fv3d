################################################################################
# @file vtk_helper.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import vtk
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
from .element import ElementType

vtk_constructor = [None] * 9
vtk_constructor[ElementType.POINT.value] = vtk.vtkVertex()
vtk_constructor[ElementType.LINE.value] = vtk.vtkLine()
vtk_constructor[ElementType.QUADRANGLE.value] = vtk.vtkQuad()
vtk_constructor[ElementType.TRIANGLE.value] = vtk.vtkTriangle()
vtk_constructor[ElementType.HEXAEDER.value] = vtk.vtkHexahedron()
vtk_constructor[ElementType.PRISM.value] = vtk.vtkWedge()
vtk_constructor[ElementType.PYRAMID.value] = vtk.vtkPyramid()
vtk_constructor[ElementType.TETRAEDER.value] = vtk.vtkTetra()


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


def set_vtk_celldata(grid, key, data):
    """Append cell data to the grid."""
    vtk_type = vtk.VTK_INT if data.dtype == np.int64 else vtk.VTK_DOUBLE
    vtkArray = numpy_to_vtk(data, array_type=vtk_type)
    vtkArray.SetName(key)

    grid.GetCellData().AddArray(vtkArray)


def shallow_copy(source):
    """Shallow copy of provided grid."""
    grid = vtk.vtkUnstructuredGrid()
    grid.ShallowCopy(source)
    return grid


def write_vtk_grid(path, grid):
    """Write the grid to the given path."""
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(grid)
    writer.Write()
