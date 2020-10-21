####################################################################################################################################
# pyPostprocessor - Python script for FV3D postprocessing
# (c) Florian Eigentler | 2020
#-----------------------------------------------------------------------------------------------------------------------------------
import argparse, h5py
import numpy as np
from lxml.etree import ElementTree, Element, SubElement

####################################################################################################################################
# Defintions
#-----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------
def main():

    # define the argument parser
    parser = argparse.ArgumentParser( description='XDMF - Generate temporal collection' )
    parser.add_argument( '-m', '--mesh-file', dest='meshfile', required=True, type=str, help='The mesh file' )
    parser.add_argument( '-d', '--debug', dest='debug', action='store_true', help='Debug the mesh file' )
    parser.add_argument( 'inifiles', nargs='*', help='The solution files' )
    args = parser.parse_args()

    # modify arguments

    # read the mesh file
    with h5py.File( args.meshfile, 'r' ) as fp_m:
        n_vertices              = int( fp_m['VERTICES'].attrs['n_vertices'] )
        n_cells                 = int( fp_m['CELLS'].attrs['n_cells'] )
        max_cell_vertices       = int( fp_m['CELLS'].attrs['max_cell_vertices'] )
        n_boundaries            = int( fp_m['BOUNDARIES'].attrs['n_boundaries'] )
        max_boundary_vertices   = int( fp_m['BOUNDARIES'].attrs['max_boundary_vertices'] )
        n_faces                 = int( fp_m['FACES'].attrs['n_faces'] )
        max_face_vertices       = int( fp_m['FACES'].attrs['max_face_vertices'] )

        if args.debug:
            xdmf        = Element( 'Xdmf', {'Version': '2.0'} )
            domain      = SubElement( xdmf, 'Domain' )
            collection  = SubElement( domain, 'Grid', {'GridType': 'Collection', 'CollectionType': 'Spatial'} )

            for name, n, m, variables_sca, variables_vec in [
                ('CELLS', n_cells, max_cell_vertices, ['id', 'volume', 'pid'], ['x']),
                ('BOUNDARIES', n_boundaries, max_boundary_vertices, ['face', 'id', 'area', 'pid'], ['n', 't1', 't2']),
                ('FACES', n_faces, max_face_vertices, ['cell1', 'cell2', 'boundary', 'id', 'area', 'lambda'], ['x', 'n', 't1', 't2']),
            ]:
                pos             = SubElement( collection, 'Grid', {'Name': name, 'GridType': 'Uniform'} )
                p_top           = SubElement( pos, 'Topology', {'TopologyType': 'Mixed', 'NumberOfElements': str( n )} )
                p_top_di        = SubElement( p_top, 'DataItem', {'Dimensions': '{:} {:}'.format( n, m ), 'Format': 'HDF'} )
                p_top_di.text   = '{:}:/{:}/xdmf'.format( args.meshfile, name )

                p_geo           = SubElement( pos, 'Geometry', {'GeometryType': 'XYZ'} )
                p_geo_di        = SubElement( p_geo, 'DataItem', {'Dimensions': '{:} 3'.format( n_vertices ), 'Format': 'HDF'} )
                p_geo_di.text   = '{:}:/VERTICES/x'.format( args.meshfile )

                for variable in variables_sca:
                    if not variable in fp_m[name]: continue
                    p_att           = SubElement( pos, 'Attribute', {'Name': variable, 'AttributeType': 'Scalar', 'Center': 'Cell'} )
                    p_att_di        = SubElement( p_att, 'DataItem', {'Dimensions': str( n ), 'Format': 'HDF'} )
                    p_att_di.text   = '{:}:/{:}/{:}'.format( args.meshfile, name, variable )

                for variable in variables_vec:
                    if not variable in fp_m[name]: continue
                    p_att           = SubElement( pos, 'Attribute', {'Name': variable, 'AttributeType': 'Vector', 'Center': 'Cell'} )
                    p_att_di        = SubElement( p_att, 'DataItem', {'Dimensions': '{:} 3'.format( n ), 'Format': 'HDF'} )
                    p_att_di.text   = '{:}:/{:}/{:}'.format( args.meshfile, name, variable )

            tmp = ElementTree( xdmf )
            tmp.write( '{:}.xdmf'.format( args.meshfile ), encoding='UTF-8',
                doctype='<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>', xml_declaration=True, pretty_print=True )

    # read the solution file
    for inifile in args.inifiles:
        with h5py.File( inifile, 'r' ) as fp_h:
            all_variables   = [str( x.decode() ) for x in fp_h['variables']]
            variables       = [str( x.decode() ) for x in fp_h['variables']]
            vec_variables   = []
            if 'rho_u' in variables:
                vec_variables.append( ('rho_uvw', ('rho_u', 'rho_v', 'rho_w')) )
                vec_variables.append( ('uvw (prim)', ('u (prim)', 'v (prim)', 'w (prim)')) )
                for x in vec_variables:
                    for y in x[1]: variables.remove( y )

            solutions   = sorted( [key for key in fp_h['SOLUTIONS'] if key != 'last'] )

            xdmf        = Element( 'Xdmf', {'Version': '2.0'} )
            domain      = SubElement( xdmf, 'Domain' )
            collection  = SubElement( domain, 'Grid', {'GridType': 'Collection', 'CollectionType': 'Temporal'} )

            for i, iteration in enumerate( solutions ):
                pos = SubElement( collection, 'Grid', {'Name': iteration, 'GridType': 'Uniform'} )
                SubElement( pos, 'Time', {'Type': 'Single', 'Value': str( i )} )

                p_top           = SubElement( pos, 'Topology', {'TopologyType': 'Mixed', 'NumberOfElements': str( n_cells )} )
                p_top_di        = SubElement( p_top, 'DataItem', {'Dimensions': '{:} {:}'.format( n_cells, max_cell_vertices ), 'Format': 'HDF'} )
                p_top_di.text   = '{:}:/CELLS/xdmf'.format( args.meshfile )

                p_geo           = SubElement( pos, 'Geometry', {'GeometryType': 'XYZ'} )
                p_geo_di        = SubElement( p_geo, 'DataItem', {'Dimensions': '{:} 3'.format( n_vertices ), 'Format': 'HDF'} )
                p_geo_di.text   = '{:}:/VERTICES/x'.format( args.meshfile )

                for variable in variables:
                    j = all_variables.index( variable )

                    p_att           = SubElement( pos, 'Attribute', {'Name': variable, 'AttributeType': 'Scalar', 'Center': 'Cell'} )
                    p_att_di        = SubElement( p_att, 'DataItem', {'ItemType': 'HyperSlab', 'Dimensions': str( n_cells ), 'Type': 'HyperSlab'} )
                    p_att_di2       = SubElement( p_att_di, 'DataItem', {'Dimensions': '3 2', 'Format': 'XML'} )
                    p_att_di2.text  = '0 {:} 1 1 {:} 1'.format( j, n_cells )
                    p_att_di3       = SubElement( p_att_di, 'DataItem', {'Dimensions': '{:} {:}'.format( n_cells, len( variables ) ), 'Format': 'HDF'} )
                    p_att_di3.text  = '{:}:/SOLUTIONS/{:}/phi_total'.format( inifile, iteration )

                for variable, components in vec_variables:
                    j = all_variables.index( components[0] )

                    p_att           = SubElement( pos, 'Attribute', {'Name': variable, 'AttributeType': 'Vector', 'Center': 'Cell'} )
                    p_att_di        = SubElement( p_att, 'DataItem', {'ItemType': 'HyperSlab', 'Dimensions': '{:} {:}'.format( n_cells, len( components ) ), 'Type': 'HyperSlab'} )
                    p_att_di2       = SubElement( p_att_di, 'DataItem', {'Dimensions': '3 2', 'Format': 'XML'} )
                    p_att_di2.text  = '0 {:} 1 1 {:} {:}'.format( j, n_cells, len( components ) )
                    p_att_di3       = SubElement( p_att_di, 'DataItem', {'Dimensions': '{:} {:}'.format( n_cells, len( variables ) ), 'Format': 'HDF'} )
                    p_att_di3.text  = '{:}:/SOLUTIONS/{:}/phi_total'.format( inifile, iteration )

            tmp = ElementTree( xdmf )
            tmp.write( '{:}.xdmf'.format( inifile ), encoding='UTF-8',
                doctype='<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>', xml_declaration=True, pretty_print=True )

####################################################################################################################################
# CALL BY SCRIPT
#-----------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    main()