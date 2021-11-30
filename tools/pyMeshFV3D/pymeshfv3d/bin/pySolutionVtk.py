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
import h5py
from pymeshfv3d import LessThanFilter
from pymeshfv3d import generate_vtk_grid, set_vtk_celldata
from pymeshfv3d import write_vtk_grid, shallow_copy
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


def main():
    """Main function entrance point."""

    # define the argument parser
    parser = argparse.ArgumentParser(
        description='pySolutionVtk - Solution to Vtk converter')
    parser.add_argument('-d', '--debugging',
                        action='store_true', help='Debugging output')
    parser.add_argument('-m', '--mesh-file',
                        required=True, type=str, help='The mesh file')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s (version {:})'.format(__version__))
    parser.add_argument('solution_files', nargs='+',
                        help='The solution file')
    args = parser.parse_args()

    with h5py.File(args.mesh_file, 'r') as h5_root:
        # vertice data
        points = h5_root['VERTICES']['x'][:]

        # cell data
        cell_types = h5_root['CELLS']['type'][:]
        cell_n_vertices = h5_root['CELLS']['n_vertices'][:]
        cell_vertices = h5_root['CELLS']['vertices'][:]

    cells = generate_vtk_grid(
        points, cell_types, cell_n_vertices, cell_vertices)

    for solution_file in args.solution_files:
        with h5py.File(solution_file, 'r') as h5_root:
            variables = [str(x.decode()) for x in h5_root['tot_variables']]
            phi_total = h5_root['phi_total'][:]

            grid = shallow_copy(cells)
            for i, variable in enumerate(variables):
                set_vtk_celldata(
                    grid, variable, phi_total[:, i])

            file_name = solution_file.replace('.h5', '.vtu')
            write_vtk_grid(file_name, grid)


################################################################################
# CALL BY SCRIPT
# ------------------------------------------------------------------------------
if __name__ == "__main__":

    main()
