################################################################################
# @file mesh.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import numpy as np
from pymeshfv3d.utilities import ElementType


class Mesh():
    """Mesh object."""

    def __init__(self, name, elements, regions, vertices):
        self.name = name
        self.elements = elements
        self.regions = regions
        self.vertices = np.array(vertices)

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return '"{:}" (d={:}, ne={:}, nr={:}, nv={:})'.format(
            self.name, self.dimension, self.n_elements, self.n_regions, self.n_vertices)

    @property
    def dimension(self):
        """The mesh dimension by element type."""
        types = [x.element_type for x in self.elements]
        if any(x.value > ElementType.QUADRANGLE.value for x in types):
            return 3
        if any(x.value > ElementType.LINE.value for x in types):
            return 2

        return 1

    @property
    def n_elements(self):
        """The number of elements."""
        return len(self.elements)

    @property
    def n_regions(self):
        """The number of regions."""
        return len(self.regions)

    @property
    def n_vertices(self):
        """The number of vertices."""
        return len(self.vertices)
