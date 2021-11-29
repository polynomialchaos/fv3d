################################################################################
# @file pyMeshFV3D.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import sys
import argparse
import logging
from pymeshfv3d import read_handler_lin1d, read_handler_gmsh
from pymeshfv3d import get_param, read_param_file, write_param_file
from pymeshfv3d import process_mesh, write_mesh, LessThanFilter
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
        description='pyMeshFV3D - Finite volume solver (FV3D) mesh preprocessing')
    parser.add_argument('-d', '--debugging',
                        action='store_true', help='Debugging output')
    parser.add_argument('-g', '--generate', dest='generate',
                        action='store_true', help='The preprocessing input file')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s (version {:})'.format(__version__))
    parser.add_argument('inifile', nargs='?',
                        help='The preprocessing input file')
    args = parser.parse_args()

    # setup additional output for debugging
    if args.debugging:
        debug_handler = logging.FileHandler(
            '{:}.log'.format(args.inifile), mode='w')
        debug_handler.setLevel(logging.DEBUG)
        debug_handler.setFormatter(logging.Formatter(
            '%(levelname)s:%(name)s:%(message)s'))
        logging.getLogger().addHandler(debug_handler)

    # modify arguments
    if args.generate:
        write_param_file(args.inifile)
        raise SystemExit

    # read the input file
    arguments = read_param_file(args.inifile)

    # define/read the mesh
    mesh_reader = get_param(arguments, 'Reader/mesh_reader')
    if mesh_reader == 'Gmsh':
        mesh = read_handler_gmsh(
            get_param(arguments, 'Reader/mesh_file'))
    elif mesh_reader == 'Lin1D':
        x_left = get_param(arguments, 'Reader/Lin1D/x_left')
        x_right = get_param(arguments, 'Reader/Lin1D/x_right')
        n_elements = get_param(arguments, 'Reader/Lin1D/n_elements')
        mesh = read_handler_lin1d(x_left, x_right, n_elements)
    else:
        raise KeyError('Reader/mesh_reader', mesh_reader)

    # process the mesh
    mesh_scale = get_param(arguments, 'Process/mesh_scale')
    n_partitions = get_param(arguments, 'Process/n_partitions')
    mesh = process_mesh(mesh, mesh_scale, n_partitions)

    # write the mesh
    file_name = '{:}.mesh.h5'.format(
        get_param(arguments, 'General/title'))
    write_mesh(file_name, mesh)


################################################################################
# CALL BY SCRIPT
# ------------------------------------------------------------------------------
if __name__ == "__main__":

    main()
