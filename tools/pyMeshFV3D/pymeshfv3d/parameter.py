################################################################################
# @file parameter.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import json

default_parameters = {
    "General": {
        "title": {
            "value": "untitled",
            "The project title": None
        }
    },
    "Reader": {
        "mesh_reader": {
            "value": "Gmsh",
            "The mesh reader": [
                "Gmsh",
                "Lin1D"
            ]
        },
        "mesh_file": {
            "value": "(NDEF)",
            "The mesh file": None
        },
        "Lin1D": {
            "x_left": {
                "value": 0.00000E+00,
                "The left mesh position": None
            },
            "x_right": {
                "value": 1.00000E+00,
                "The right mesh position": None
            },
            "n_elements": {
                "value": 10,
                "The number of elements": None
            }
        }
    },
    "Process": {
        "mesh_scale": {
            "value": [
                1.00000E+00,
                1.00000E+00,
                1.00000E+00
            ],
            "The mesh scale factor for each dimension": None
        },
        "n_partitions": {
            "value": 0,
            "The number of partitions": None
        }
    }
}


def get_param(parameters, path):
    """Get the parametery by path."""
    root = parameters
    for key in path.split('/'):
        root = root[key]

    if 'value' in root:
        return root['value']

    return root


def read_param_file(path):
    """Read the parameter file."""
    with open(path, 'r') as fp_ptr:
        data = json.load(fp_ptr)

    return data


def write_param_file(path):
    """Write a parameter file."""
    with open(path, 'w') as fp_ptr:
        json.dump(default_parameters, fp_ptr, indent=4)
