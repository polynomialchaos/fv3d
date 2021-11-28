################################################################################
# @file lin1d.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging
import numpy as np
from pymeshfv3d.utilities import ElementType, Element, Mesh


def read_handler_lin1d(x_left=0.0, x_right=1.0, n=10):
    """Generate a mesh (1D) by definition of left/right position and
    number of elements."""
    logging.info('Generate mesh (Lin1D) xL=%s, xR=%s, and n=%s',
                 x_left, x_right, n)

    regions = ['left', 'right', 'flow']
    vertices = [(x, 0.0, 0.0) for x in np.linspace(x_left, x_right, n+1)]

    elements = []
    elements.append(Element(ElementType.POINT, 0, regions.index('left')))

    for i in range(n):
        elements.append(
            Element(ElementType.LINE, (i, i+1), regions.index('flow')))

    elements.append(Element(ElementType.POINT, n, regions.index('right')))

    return Mesh('Lin1D', elements, regions, vertices)
