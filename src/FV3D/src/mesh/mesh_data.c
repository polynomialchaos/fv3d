//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "mesh_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void deallocate_partition( Partition_t **partition );
void print_partition( Partition_t *partition );
void deallocate_vertices( Vertices_t **vertices );
void print_vertices( Vertices_t *vertices );
void deallocate_cells( Cells_t **cells );
void print_cells( Cells_t *cells );
void deallocate_boundaries( Boundaries_t **boundaries );
void print_boundaries( Boundaries_t *boundaries );
void deallocate_faces( Faces_t **faces );
void print_faces( Faces_t *faces );
void deallocate_regions( Regions_t **regions );
void print_regions( Regions_t *regions );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
Mesh_t *allocate_mesh()
{
    Mesh_t *tmp = allocate( sizeof( Mesh_t ) );

    tmp->dimension      = 0;
    tmp->is_partitioned = 0;

    tmp->partition  = NULL;
    tmp->vertices   = NULL;
    tmp->cells      = NULL;
    tmp->boundaries = NULL;
    tmp->faces      = NULL;
    tmp->regions    = NULL;

    tmp->local_volume   = 0.0;
    tmp->global_volume  = 0.0;

    return tmp;
}

void print_mesh_info( Mesh_t *mesh )
{
    if (mesh == NULL) return;
    printf_r( "MESH\n" );
    printf_r( "dimension      = %d\n", mesh->dimension      );
    printf_r( "is_partitioned = %d\n", mesh->is_partitioned );

    print_partition( mesh->partition );
    print_vertices( mesh->vertices );
    print_cells( mesh->cells );
    print_boundaries( mesh->boundaries );
    print_faces( mesh->faces );
    print_regions( mesh->regions );

    printf_r( "local_volume  = %e\n", mesh->local_volume );
    printf_r( "global_volume = %e\n", mesh->global_volume );
}

void deallocate_mesh( Mesh_t **mesh )
{
    if ((*mesh) == NULL) return;

    deallocate_partition( &(*mesh)->partition );
    deallocate_vertices( &(*mesh)->vertices );
    deallocate_cells( &(*mesh)->cells );
    deallocate_boundaries( &(*mesh)->boundaries );
    deallocate_faces( &(*mesh)->faces );
    deallocate_regions( &(*mesh)->regions );

    deallocate( *mesh );
}

Partition_t *allocate_partition( Mesh_t *mesh, int n_partitions, int n_partition_cells, int n_partition_boundaries,
    int n_partition_faces, int n_partition_sends, int n_partition_receives )
{
    mesh->partition = allocate( sizeof( Partition_t ) );
    Partition_t *partition = mesh->partition;

    partition->n_partitions                 = n_partitions;
    partition->n_partition_cells            = n_partition_cells;
    partition->n_partition_boundaries       = n_partition_boundaries;
    partition->n_partition_faces            = n_partition_faces;
    partition->n_partition_sends            = n_partition_sends;
    partition->n_partition_receives         = n_partition_receives;

    partition->partition_cells              = allocate( sizeof( int ) * n_partition_cells );
    partition->partition_boundaries         = allocate( sizeof( int ) * n_partition_boundaries );
    partition->partition_faces              = allocate( sizeof( int ) * n_partition_faces );
    partition->partition_sends              = allocate( sizeof( int ) * n_partition_sends );
    partition->partition_sends_pid          = allocate( sizeof( int ) * n_partition_sends );
    partition->partition_receives           = allocate( sizeof( int ) * n_partition_receives );
    partition->partition_receives_pid       = allocate( sizeof( int ) * n_partition_receives );

    partition->n_partition_sends_to         = allocate( sizeof( int ) * n_partitions );
    partition->partition_sends_to           = allocate( sizeof( int ) * n_partitions * n_partition_sends );
    partition->n_partition_receives_from    = allocate( sizeof( int ) * n_partitions );
    partition->partition_receives_from      = allocate( sizeof( int ) * n_partitions * n_partition_receives );

    set_value_int_n( 0, partition->n_partition_sends_to, n_partitions );
    set_value_int_n( 0, partition->n_partition_receives_from, n_partitions );

    return partition;
}

Vertices_t *allocate_vertices( Mesh_t *mesh, int n_vertices )
{
    mesh->vertices = allocate( sizeof( Vertices_t ) );
    Vertices_t *vertices = mesh->vertices;

    vertices->n_vertices    = n_vertices;

    vertices->x             = allocate( sizeof( double ) * n_vertices * DIM );

    return vertices;
}

Cells_t *allocate_cells( Mesh_t *mesh, int n_cells, int max_cell_vertices, int max_cell_faces )
{
    mesh->cells = allocate( sizeof( Cells_t ) );
    Cells_t *cells = mesh->cells;

    cells->n_global_cells       = n_cells;
    cells->n_local_cells        = n_cells;
    cells->n_domain_cells       = n_cells;
    cells->max_cell_vertices    = max_cell_vertices;
    cells->max_cell_faces       = max_cell_faces;

    if (get_is_parallel())
    {
        cells->n_domain_cells   = mesh->partition->n_partition_cells;
        cells->n_local_cells    = mesh->partition->n_partition_cells + mesh->partition->n_partition_receives;
        n_cells                 = cells->n_local_cells;
    }

    cells->id           = allocate( sizeof( int ) * n_cells );
    cells->type         = allocate( sizeof( int ) * n_cells );
    cells->n_vertices   = allocate( sizeof( int ) * n_cells );
    cells->vertices     = allocate( sizeof( int ) * n_cells * max_cell_vertices );
    cells->n_faces      = allocate( sizeof( int ) * n_cells );
    cells->faces        = allocate( sizeof( int ) * n_cells * max_cell_faces );
    cells->x            = allocate( sizeof( double ) * n_cells * DIM );
    cells->volume       = allocate( sizeof( double ) * n_cells );
    cells->dx           = allocate( sizeof( double ) * n_cells * DIM );

    set_value_int_n( 0, cells->n_vertices, n_cells );
    set_value_int_n( 0, cells->n_faces, n_cells );

    return cells;
}

Boundaries_t *allocate_boundaries( Mesh_t *mesh, int n_boundaries, int max_boundary_vertices )
{
    mesh->boundaries = allocate( sizeof( Boundaries_t ) );
    Boundaries_t *boundaries = mesh->boundaries;

    boundaries->n_global_boundaries     = n_boundaries;
    boundaries->n_boundaries            = n_boundaries;
    boundaries->max_boundary_vertices   = max_boundary_vertices;

    if (get_is_parallel())
    {
        boundaries->n_boundaries    = mesh->partition->n_partition_boundaries;
        n_boundaries                = boundaries->n_boundaries;
    }

    boundaries->id          = allocate( sizeof( int ) * n_boundaries );
    boundaries->type        = allocate( sizeof( int ) * n_boundaries );
    boundaries->n_vertices  = allocate( sizeof( int ) * n_boundaries );
    boundaries->vertices    = allocate( sizeof( int ) * n_boundaries * max_boundary_vertices );
    boundaries->face        = allocate( sizeof( int ) * n_boundaries );
    boundaries->distance    = allocate( sizeof( double ) * n_boundaries );
    boundaries->n           = allocate( sizeof( double ) * n_boundaries * DIM );
    boundaries->t1          = allocate( sizeof( double ) * n_boundaries * DIM );
    boundaries->t2          = allocate( sizeof( double ) * n_boundaries * DIM );

    set_value_int_n( 0, boundaries->n_vertices, n_boundaries );

    return boundaries;
}

Faces_t *allocate_faces( Mesh_t *mesh, int n_faces, int max_face_vertices )
{
    mesh->faces = allocate( sizeof( Faces_t ) );
    Faces_t *faces = mesh->faces;

    faces->n_global_faces       = n_faces;
    faces->n_faces              = n_faces;
    faces->max_face_vertices    = max_face_vertices;

    if (get_is_parallel())
    {
        faces->n_faces  = mesh->partition->n_partition_faces;
        n_faces         = faces->n_faces;
    }

    faces->type         = allocate( sizeof( int ) * n_faces );
    faces->n_vertices   = allocate( sizeof( int ) * n_faces );
    faces->vertices     = allocate( sizeof( int ) * n_faces * max_face_vertices );
    faces->cells        = allocate( sizeof( int ) * n_faces * FACE_CELLS );
    faces->boundary     = allocate( sizeof( int ) * n_faces );
    faces->area         = allocate( sizeof( double ) * n_faces );
    faces->lambda       = allocate( sizeof( double ) * n_faces );
    faces->x            = allocate( sizeof( double ) * n_faces * DIM );
    faces->n            = allocate( sizeof( double ) * n_faces * DIM );
    faces->t1           = allocate( sizeof( double ) * n_faces * DIM );
    faces->t2           = allocate( sizeof( double ) * n_faces * DIM );

    set_value_int_n( 0, faces->n_vertices, n_faces );

    faces->n_internal_faces = 0;
    faces->n_boundary_faces = 0;
    faces->dist_cell_1      = NULL;
    faces->dist_cell_2      = NULL;
    faces->internal_faces   = NULL;
    faces->boundary_faces   = NULL;

    return faces;
}

Regions_t *allocate_regions( Mesh_t *mesh, int n_regions, int max_name_length )
{
    mesh->regions = allocate( sizeof( Regions_t ) );
    Regions_t *regions = mesh->regions;

    regions->n_regions          = n_regions;
    regions->max_name_length    = max_name_length;

    regions->name           = allocate_hdf5_string_buffer( n_regions, max_name_length );
    regions->is_boundary    = allocate( sizeof( int ) * n_regions );

    return regions;
}

void deallocate_partition( Partition_t **partition )
{
    if ((*partition) == NULL) return;

    deallocate( (*partition)->partition_cells );
    deallocate( (*partition)->partition_boundaries );
    deallocate( (*partition)->partition_faces );
    deallocate( (*partition)->partition_sends );
    deallocate( (*partition)->partition_sends_pid );
    deallocate( (*partition)->partition_receives );
    deallocate( (*partition)->partition_receives_pid );

    deallocate( (*partition)->n_partition_sends_to );
    deallocate( (*partition)->partition_sends_to );
    deallocate( (*partition)->n_partition_receives_from );
    deallocate( (*partition)->partition_receives_from );

    deallocate( *partition );
}

void print_partition( Partition_t *partition )
{
    if (partition == NULL) return;
    printf_r( "PARTITION\n" );

    printf_r( "n_partitions           = %d\n", partition->n_partitions           );
    printf_r( "n_partition_cells      = %d\n", partition->n_partition_cells      );
    printf_r( "n_partition_boundaries = %d\n", partition->n_partition_boundaries );
    printf_r( "n_partition_faces      = %d\n", partition->n_partition_faces      );
    printf_r( "n_partition_sends      = %d\n", partition->n_partition_sends      );
    printf_r( "n_partition_receives   = %d\n", partition->n_partition_receives   );
}

void deallocate_vertices( Vertices_t **vertices )
{
    if ((*vertices) == NULL) return;

    deallocate( (*vertices)->x );

    deallocate( *vertices );
}

void print_vertices( Vertices_t *vertices )
{
    if (vertices == NULL) return;
    printf_r( "VERTICES\n" );

    printf_r( "n_vertices = %d\n", vertices->n_vertices );
}

void deallocate_cells( Cells_t **cells )
{
    if ((*cells) == NULL) return;

    deallocate( (*cells)->id );
    deallocate( (*cells)->type );
    deallocate( (*cells)->n_vertices );
    deallocate( (*cells)->vertices );
    deallocate( (*cells)->n_faces );
    deallocate( (*cells)->faces );
    deallocate( (*cells)->x );
    deallocate( (*cells)->volume );
    deallocate( (*cells)->dx );

    deallocate( *cells );
}

void print_cells( Cells_t *cells )
{
    if (cells == NULL) return;
    printf_r( "CELLS\n" );

    printf_r( "n_global_cells    = %d\n", cells->n_global_cells    );
    printf_r( "n_local_cells     = %d\n", cells->n_local_cells     );
    printf_r( "n_domain_cells    = %d\n", cells->n_domain_cells    );
    printf_r( "max_cell_vertices = %d\n", cells->max_cell_vertices );
    printf_r( "max_cell_faces    = %d\n", cells->max_cell_faces    );
}

void deallocate_boundaries( Boundaries_t **boundaries )
{
    if ((*boundaries) == NULL) return;

    deallocate( (*boundaries)->id );
    deallocate( (*boundaries)->type );
    deallocate( (*boundaries)->n_vertices );
    deallocate( (*boundaries)->vertices );
    deallocate( (*boundaries)->face );
    deallocate( (*boundaries)->distance );
    deallocate( (*boundaries)->n );
    deallocate( (*boundaries)->t1 );
    deallocate( (*boundaries)->t2 );

    deallocate( *boundaries );
}

void print_boundaries( Boundaries_t *boundaries )
{
    if (boundaries == NULL) return;
    printf_r( "BOUNDARIES\n" );

    printf_r( "n_global_boundaries   = %d\n", boundaries->n_global_boundaries   );
    printf_r( "n_boundaries          = %d\n", boundaries->n_boundaries          );
    printf_r( "max_boundary_vertices = %d\n", boundaries->max_boundary_vertices );
}

void deallocate_faces( Faces_t **faces )
{
    if ((*faces) == NULL) return;

    deallocate( (*faces)->type );
    deallocate( (*faces)->n_vertices );
    deallocate( (*faces)->vertices );
    deallocate( (*faces)->cells );
    deallocate( (*faces)->boundary );
    deallocate( (*faces)->area );
    deallocate( (*faces)->lambda );
    deallocate( (*faces)->x );
    deallocate( (*faces)->n );
    deallocate( (*faces)->t1 );
    deallocate( (*faces)->t2 );

    deallocate( (*faces)->dist_cell_1 );
    deallocate( (*faces)->dist_cell_2 );
    deallocate( (*faces)->internal_faces );
    deallocate( (*faces)->boundary_faces );

    deallocate( *faces );
}

void print_faces( Faces_t *faces )
{
    if (faces == NULL) return;
    printf_r( "FACES\n" );

    printf_r( "n_global_faces    = %d\n", faces->n_global_faces    );
    printf_r( "n_faces           = %d\n", faces->n_faces           );
    printf_r( "max_face_vertices = %d\n", faces->max_face_vertices );
}

void deallocate_regions( Regions_t **regions )
{
    if ((*regions) == NULL) return;

    deallocate_hdf5_string_buffer( &(*regions)->name );
    deallocate( (*regions)->is_boundary );

    deallocate( *regions );
}

void print_regions( Regions_t *regions )
{
    if (regions == NULL) return;
    printf_r( "REGIONS\n" );

    printf_r( "n_regions   = %d\n", regions->n_regions   );
    printf_r( "flow_region = %d\n", regions->flow_region );

    for ( int i = 0; i < regions->n_regions; i++ )
        printf_r( "%d: %s\n", i, regions->name[i] );
}