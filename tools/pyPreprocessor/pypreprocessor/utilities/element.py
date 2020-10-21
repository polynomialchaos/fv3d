####################################################################################################################################
# pyPreprocessor - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
#-----------------------------------------------------------------------------------------------------------------------------------
from enum import Enum, unique
import numpy as np

####################################################################################################################################
# Defintions
#-----------------------------------------------------------------------------------------------------------------------------------
@unique
class ElementType( Enum ):
    POINT       = 1
    LINE        = 2
    TRIANGLE    = 3
    QUADRANGLE  = 4
    TETRAEDER   = 5
    HEXAEDER    = 6 
    PRISM       = 7
    PYRAMID     = 8

class Element( object ):
    def __init__( self, element_type, vertice_ids, region_id ):
        self.element_type   = element_type
        self.vertice_ids    = list( np.atleast_1d( vertice_ids ) )
        self.region_id      = region_id

    @property
    def n_vertice_ids( self ):
        return len( self.vertice_ids )

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return '{:}'.format( self.__dict__ )

class Cell( Element ):
    def __init__( self, element_type, vertice_ids, region_id, face_ids=None, x=None, volume=None, partition_id=None ):
        super().__init__( element_type, vertice_ids, region_id )
        self.face_ids       = [] if face_ids is None else list( np.atleast_1d( face_ids ) ) 
        self.x              = x
        self.volume         = volume
        self.partition_id   = partition_id

    def get_faces( self ):
        v = self.vertice_ids
        tmp_faces = []
    
        if self.element_type == ElementType.LINE:
            tmp_faces.append( Face( ElementType.POINT, v[0] ) )
            tmp_faces.append( Face( ElementType.POINT, v[1] ) )
        elif self.element_type == ElementType.TRIANGLE:
            tmp_faces.append( Face( ElementType.LINE, (v[0], v[1]) ) )
            tmp_faces.append( Face( ElementType.LINE, (v[1], v[2]) ) )
            tmp_faces.append( Face( ElementType.LINE, (v[2], v[0]) ) )
        elif self.element_type == ElementType.QUADRANGLE:
            tmp_faces.append( Face( ElementType.LINE, (v[0], v[1]) ) )
            tmp_faces.append( Face( ElementType.LINE, (v[1], v[2]) ) )
            tmp_faces.append( Face( ElementType.LINE, (v[2], v[3]) ) )
            tmp_faces.append( Face( ElementType.LINE, (v[3], v[0]) ) )
        elif self.element_type == ElementType.TETRAEDER:
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[1], v[2]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[3], v[1]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[3], v[2]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[2], v[3], v[1]) ) )
        elif self.element_type == ElementType.HEXAEDER:
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[1], v[2], v[3]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[4], v[5], v[6], v[7]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[4], v[7], v[3]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[4], v[5], v[1]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[1], v[5], v[6], v[2]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[3], v[7], v[6], v[2]) ) )
        elif self.element_type == ElementType.PRISM:
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[1], v[4], v[3]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[1], v[2], v[5], v[4]) ) )
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[2], v[5], v[3]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[1], v[2]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[3], v[4], v[5]) ) )
        elif self.element_type == ElementType.PYRAMID:
            tmp_faces.append( Face( ElementType.QUADRANGLE, (v[0], v[1], v[2], v[3]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[1], v[4]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[1], v[2], v[4]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[2], v[3], v[4]) ) )
            tmp_faces.append( Face( ElementType.TRIANGLE, (v[0], v[3], v[1]) ) )
        else:
            raise TypeError( self.element_type )
        
        return tmp_faces

    @property
    def n_face_ids( self ):
        return len( self.face_ids )

class Boundary( Element ):
    def __init__( self, element_type, vertice_ids, region_id, 
        face_id=None, n=None, t1=None, t2=None, distance=None, partition_id=None ):
        super().__init__( element_type, vertice_ids, region_id )
        self.face_id        = face_id
        self.n              = n
        self.t1             = t1
        self.t2             = t2
        self.distance       = distance
        self.partition_id   = partition_id

class Face( Element ):
    def __init__( self, element_type, vertice_ids, cell_ids=None, boundary_id=None,
        x=None, n=None, t1=None, t2=None, area=None, weight=None ):
        super().__init__( element_type, vertice_ids, None )
        self.cell_ids       = [-1, -1] if cell_ids is None else cell_ids
        self.boundary_id    = boundary_id
        self.x              = x
        self.n              = n
        self.t1             = t1
        self.t2             = t2
        self.area           = area
        self.weight         = weight

####################################################################################################################################
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------