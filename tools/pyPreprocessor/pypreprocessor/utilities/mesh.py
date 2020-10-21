####################################################################################################################################
# pyPreprocessor - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
#-----------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from pypreprocessor.utilities import ElementType, Cell, Boundary, Face

####################################################################################################################################
# Defintions
#-----------------------------------------------------------------------------------------------------------------------------------
class Mesh( object ):
    def __init__( self, vertices=None, regions=None, elements=None, title='untitled' ):
        self.vertices   = np.array( [] if vertices is None else vertices )
        self.regions    = [] if regions is None else regions
        self.elements   = [] if elements is None else elements
        self.title      = title

    @property
    def dimension( self ):
        types = [x.element_type for x in self.elements]
        if any( x.value > ElementType.QUADRANGLE.value for x in types ): return 3
        elif any( x.value > ElementType.LINE.value for x in types ): return 2
        else: return 1

    @property
    def n_vertices( self ):
        return len( self.vertices )

    @property
    def n_regions( self ):
        return len( self.regions )

    @property
    def n_elements( self ):
        return len( self.elements )

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return '"{:}" (nv={:}, nr={:}, ne={:}, d={:})'.format(
            self.title, self.n_vertices, self.n_regions, self.n_elements, self.dimension )

####################################################################################################################################
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------