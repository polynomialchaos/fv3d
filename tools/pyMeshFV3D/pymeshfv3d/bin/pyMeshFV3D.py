################################################################################
# @file pyMeshFV3D.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import os
import argparse
import json
import shutil
from pymeshfv3d.reader import read_handler_lin1d, read_handler_gmsh
from pychemistry.version import __version__
from pymeshfv3d import process_mesh, write_mesh

main_path = os.path.abspath(__file__)


def get_arg_by_path(arguments, path, value='/value'):
    root = arguments
    for x in (path + value).split('/'):
        root = root[x]

    return root


def main():

    # define the argument parser
    parser = argparse.ArgumentParser(
        description='pyMeshFV3D - Finite volume solver (FV3D) mesh preprocessing')
    parser.add_argument('-g', '--generate', dest='generate',
                        action='store_true', help='The preprocessing input file')
    parser.add_argument('inifile', nargs='+',
                        help='The preprocessing input file')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s (version {:})'.format(__version__))
    args = parser.parse_args()

    # modify arguments
    args.inifile = args.inifile[0]
    if args.generate:
        tpl_file = os.path.join(os.path.dirname(main_path), 'template.json')
        shutil.copyfile(tpl_file, args.inifile)
        raise SystemExit

    # read the input file
    with open(args.inifile, 'r') as fp:
        arguments = json.load(fp)

    # define/read the mesh
    mesh_reader = get_arg_by_path(arguments, 'Reader/mesh_reader')
    if mesh_reader == 'Gmsh':
        mesh = read_handler_gmsh(
            get_arg_by_path(arguments, 'Reader/mesh_file'))
    elif mesh_reader == 'Lin1D':
        x_left = get_arg_by_path(arguments, 'Reader/Lin1D/x_left')
        x_right = get_arg_by_path(arguments, 'Reader/Lin1D/x_right')
        n_elements = get_arg_by_path(arguments, 'Reader/Lin1D/n_elements')
        mesh = read_handler_lin1d(x_left, x_right, n_elements)
    else:
        raise(KeyError('Reader/mesh_reader', mesh_reader))

    # process the mesh
    mesh_scale = get_arg_by_path(arguments, 'Process/mesh_scale')
    n_partitions = get_arg_by_path(arguments, 'Process/n_partitions')
    mesh = process_mesh(mesh, mesh_scale, n_partitions)

    # write the mesh
    file_name = '{:}.mesh.h5'.format(
        get_arg_by_path(arguments, 'General/title'))
    write_mesh(file_name, mesh)


####################################################################################################################################
# CALL BY SCRIPT
# -----------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    main()
