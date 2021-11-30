################################################################################
# @file __init__.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .parameter import get_param, read_param_file, write_param_file
from .process import process_mesh
from .reader import read_handler_lin1d, read_handler_gmsh
from .utilities import LessThanFilter
from .utilities import generate_vtk_grid, set_vtk_celldata
from .utilities import write_vtk_grid, shallow_copy
from .write import write_mesh
