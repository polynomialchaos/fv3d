################################################################################
# @file math.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import math
import numpy as np


def len_vector(a):
    """Length of vector."""
    return math.sqrt(np.dot(a, a))


def norm_vector(a):
    """Normalize vector."""
    return a / len_vector(a)


def cross_product(a, b):
    """Cross-product of two vectors in 3D."""
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ])
