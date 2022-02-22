################################################################################
# @file __init__.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
from .utilities import LessThanFilter
from .element import ElementType, Element, Cell, Boundary, Face, Region
from .geometry import calc_area, calc_distance, calc_normal
from .geometry import calc_t1, calc_t2, calc_volume, calc_weight
from .mesh import Mesh
from .metis import MetisLibrary
from .vtk_helper import generate_vtk_grid, set_vtk_celldata
from .vtk_helper import write_vtk_grid, shallow_copy