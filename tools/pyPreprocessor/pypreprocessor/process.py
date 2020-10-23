####################################################################################################################################
# pyPreprocessor - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
#-----------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from pypreprocessor.utilities import Mesh, ElementType, Cell, Boundary, Cell, Region, metis_part_mesh_nodal
from pypreprocessor.utilities import calc_area, calc_distance, calc_normal, calc_t1, calc_t2, calc_volume, calc_weight

####################################################################################################################################
# Defintions
#-----------------------------------------------------------------------------------------------------------------------------------
class Continue( Exception ):
    pass

class ProcessMesh( Mesh ):
    def __init__( self, mesh, scale, n_partitions ):
        super().__init__( mesh.vertices, mesh.regions, mesh.elements )
        self.scale          = scale
        self.n_partitions   = n_partitions

        self.cells          = []
        self.boundaries     = []
        self.faces          = []
        self.regions        = [Region( x ) for x in self.regions]

        self._element_to_cell_boundaries()
        self._generate_cell_faces()
        self._calc_geometry()
        if self.n_partitions >= 2:
            metis_part_mesh_nodal( self )

    def _calc_geometry( self ):
        self.vertices[:,0] *= self.scale[0]
        self.vertices[:,1] *= self.scale[1]
        self.vertices[:,2] *= self.scale[2]

        for cell in self.cells:
            cell.x      = np.mean( self.vertices[cell.vertice_ids], axis=0 )
            cell.volume = calc_volume( self.vertices, cell )

        for face in self.faces:
            face.x      = np.mean( self.vertices[face.vertice_ids], axis=0 )
            face.area   = calc_area( self.vertices, face )
            face.n      = calc_normal( self.vertices, self.cells, face )
            face.t1     = calc_t1( face )
            face.t2     = calc_t2( face )
            face.weight = calc_weight( self.cells, face )

        for boundary in self.boundaries:
            boundary.n          = -self.faces[boundary.face_id].n
            boundary.t1         = calc_t1( boundary )
            boundary.t2         = calc_t2( boundary )
            boundary.distance   = calc_distance( self.cells, self.faces, boundary )

            region = self.regions[boundary.region_id]
            region.is_boundary = True

        for cell in self.cells:
            cell.dx = np.array( cell.dx )
            for idx in cell.face_ids:
                face = self.faces[idx]
                cell.dx += np.abs( face.n ) * face.area

    def _element_to_cell_boundaries( self ):
        face_type = self._get_face_type()

        for element in self.elements:
            if element.element_type.value > face_type.value:
                self.cells.append( Cell( element.element_type, element.vertice_ids, element.region_id ) )
            else:
                self.boundaries.append( Boundary( element.element_type, element.vertice_ids, element.region_id ) )

    def _generate_cell_faces( self ):
        tmp_faces = []
        for i, cell in enumerate( self.cells ):
            for face in cell.get_faces():
                face.cell_ids[0] = i
                tmp_faces.append( face )

        vertice_faces = generate_vertice_faces( self.vertices, tmp_faces )

        checked = [False] * len( tmp_faces )
        for i, face in enumerate( tmp_faces ):
            if checked[i]: continue
            i_sum = sum( face.vertice_ids )

            try:
                for vid in face.vertice_ids:
                    for j in vertice_faces[vid]:
                        if i == j: continue

                        j_face = tmp_faces[j]
                        if face.element_type != j_face.element_type: continue
                        if i_sum != sum( j_face.vertice_ids ): continue
                        if set( face.vertice_ids ) != set( j_face.vertice_ids ): continue

                        checked[i] = True
                        checked[j] = True

                        face.cell_ids[1] = j_face.cell_ids[0]
                        j_face.cell_ids  = None
                        self.faces.append( face )
                        raise Continue
            except Continue:
                continue

            if face.cell_ids[1] < 0:
                self.faces.append( face )

        vertice_faces = generate_vertice_faces( self.vertices, self.faces )

        for i, face in enumerate( self.faces ):
            self.cells[face.cell_ids[0]].face_ids.append( i )
            if face.cell_ids[1] >= 0:
                self.cells[face.cell_ids[1]].face_ids.append( i )
                face.boundary_id = -1

        for i, boundary in enumerate( self.boundaries ):
            i_sum = sum( boundary.vertice_ids )

            try:
                for vid in boundary.vertice_ids:
                    for j in vertice_faces[vid]:

                        j_face = self.faces[j]
                        if j_face.cell_ids[1] >= 0: continue
                        if boundary.element_type != j_face.element_type: continue
                        if i_sum != sum( j_face.vertice_ids ): continue
                        if set( boundary.vertice_ids ) != set( j_face.vertice_ids ): continue

                        boundary.face_id = j
                        j_face.boundary_id = i
                        raise Continue
            except Continue:
                continue

    def _get_face_type( self ):
        dim = self.dimension
        if dim == 3: return ElementType.QUADRANGLE
        elif dim == 2: return ElementType.LINE
        else: return ElementType.POINT

    @property
    def is_partitioned( self ):
        return self.n_partitions > 1

    @property
    def n_cells( self ):
        return len( self.cells )

    @property
    def n_boundaries( self ):
        return len( self.boundaries )

    @property
    def n_faces( self ):
        return len( self.faces )

    @property
    def volume( self ):
        return sum( x.volume for x in self.cells )

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return '"{:}" (nv={:}, nc={:}, nb={:}, nf={:}, v={:})'.format(
            self.title, self.n_vertices, self.n_cells, self.n_boundaries, self.n_faces, self.volume )

####################################################################################################################################
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------
def generate_vertice_faces( vertices, faces ):
    vertice_faces = [[] for _ in vertices]
    for i, face in enumerate( faces ):
        for vertice in face.vertice_ids:
            vertice_faces[vertice].append( i )

    return vertice_faces

def process_mesh( mesh, scale, n_partitions ):

    print( 'Processing mesh: {:}'.format( mesh ) )

    return ProcessMesh( mesh, scale, n_partitions )