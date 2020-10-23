//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#include <string.h>
#include "mesh_private.h"

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
void deallocate_vertices( Vertices_t **vertices );
void deallocate_cells( Cells_t **cells );
void deallocate_boundaries( Boundaries_t **boundaries );
void deallocate_faces( Faces_t **faces );
void deallocate_regions( Regions_t **regions );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
Mesh_t *allocate_mesh()
{
    Mesh_t *tmp = allocate( sizeof( Mesh_t ) );

    tmp->dimension      = 0;
    tmp->is_partitioned = 0;

    tmp->partition                          = allocate( sizeof( Partition_t ) );
    tmp->partition->n_partitions            = 0;
    tmp->partition->n_partition_cells       = 0;
    tmp->partition->n_partition_boundaries  = 0;
    tmp->partition->n_partition_faces       = 0;
    tmp->partition->n_partition_sends       = 0;
    tmp->partition->n_partition_receives    = 0;

    tmp->vertices               = allocate( sizeof( Vertices_t ) );
    tmp->vertices->n_vertices   = 0;

    tmp->cells                      = allocate( sizeof( Cells_t ) );
    tmp->cells->n_global_cells      = 0;
    tmp->cells->n_local_cells       = 0;
    tmp->cells->n_cells             = 0;
    tmp->cells->max_cell_vertices   = 0;
    tmp->cells->max_cell_faces      = 0;

    tmp->boundaries                         = allocate( sizeof( Boundaries_t ) );
    tmp->boundaries->n_global_boundaries    = 0;
    tmp->boundaries->n_local_boundaries     = 0;
    tmp->boundaries->n_boundaries           = 0;
    tmp->boundaries->max_boundary_vertices  = 0;

    tmp->faces                      = allocate( sizeof( Faces_t ) );
    tmp->faces->n_global_faces      = 0;
    tmp->faces->n_local_faces       = 0;
    tmp->faces->n_faces             = 0;
    tmp->faces->max_face_vertices   = 0;

    tmp->regions                = allocate( sizeof( Regions_t ) );
    tmp->regions->n_regions     = 0;
    tmp->regions->flow_region   = -1;

    tmp->total_volume = 0.0;

    return tmp;
}

void print_mesh_info( Mesh_t *mesh )
{
    printf_r( "dimension      = %d\n", mesh->dimension      );
    printf_r( "is_partitioned = %d\n", mesh->is_partitioned );

    printf_r( "n_partitions           = %d\n", mesh->partition->n_partitions           );
    printf_r( "n_partition_cells      = %d\n", mesh->partition->n_partition_cells      );
    printf_r( "n_partition_boundaries = %d\n", mesh->partition->n_partition_boundaries );
    printf_r( "n_partition_faces      = %d\n", mesh->partition->n_partition_faces      );
    printf_r( "n_partition_sends      = %d\n", mesh->partition->n_partition_sends      );
    printf_r( "n_partition_receives   = %d\n", mesh->partition->n_partition_receives   );

    printf_r( "n_vertices = %d\n", mesh->vertices->n_vertices );

    printf_r( "n_global_cells    = %d\n", mesh->cells->n_global_cells    );
    printf_r( "n_local_cells     = %d\n", mesh->cells->n_local_cells     );
    printf_r( "n_cells           = %d\n", mesh->cells->n_cells           );
    printf_r( "max_cell_vertices = %d\n", mesh->cells->max_cell_vertices );
    printf_r( "max_cell_faces    = %d\n", mesh->cells->max_cell_faces    );

    printf_r( "n_global_boundaries   = %d\n", mesh->boundaries->n_global_boundaries   );
    printf_r( "n_local_boundaries    = %d\n", mesh->boundaries->n_local_boundaries    );
    printf_r( "n_boundaries          = %d\n", mesh->boundaries->n_boundaries          );
    printf_r( "max_boundary_vertices = %d\n", mesh->boundaries->max_boundary_vertices );

    printf_r( "n_global_faces    = %d\n", mesh->faces->n_global_faces    );
    printf_r( "n_local_faces     = %d\n", mesh->faces->n_local_faces     );
    printf_r( "n_faces           = %d\n", mesh->faces->n_faces           );
    printf_r( "max_face_vertices = %d\n", mesh->faces->max_face_vertices );

    printf_r( "n_regions   = %d\n", mesh->regions->n_regions   );
    printf_r( "flow_region = %d\n", mesh->regions->flow_region );

    printf_r( "total_volume = %e\n", mesh->total_volume );
}

void deallocate_mesh( Mesh_t **mesh )
{
    deallocate_partition( &(*mesh)->partition );
    deallocate_vertices( &(*mesh)->vertices );
    deallocate_cells( &(*mesh)->cells );
    deallocate_boundaries( &(*mesh)->boundaries );
    deallocate_faces( &(*mesh)->faces );
    deallocate_regions( &(*mesh)->regions );

    deallocate( *mesh );
}

void allocate_partition( Partition_t *partition, int n_partitions,
    int n_cells, int n_boundaries, int n_faces, int n_sends, int n_receives )
{
    partition->partition_cells              = allocate( sizeof( int ) * n_cells );
    partition->partition_boundaries         = allocate( sizeof( int ) * n_boundaries );
    partition->partition_faces              = allocate( sizeof( int ) * n_faces );
    partition->partition_sends              = allocate( sizeof( int ) * n_sends );
    partition->partition_sends_pid          = allocate( sizeof( int ) * n_sends );
    partition->partition_receives           = allocate( sizeof( int ) * n_receives );
    partition->partition_receives_pid       = allocate( sizeof( int ) * n_receives );

    partition->n_partition_sends_to         = allocate( sizeof( int ) * n_partitions );
    partition->partition_sends_to           = allocate( sizeof( int ) * n_partitions * n_sends );
    partition->n_partition_receives_from    = allocate( sizeof( int ) * n_partitions );
    partition->partition_receives_from      = allocate( sizeof( int ) * n_partitions * n_receives );

    set_value_int_n( 0, partition->n_partition_sends_to, n_partitions );
    set_value_int_n( 0, partition->n_partition_receives_from, n_partitions );
}

void allocate_vertices( Vertices_t *vertices, int n_vertices )
{
    vertices->x = allocate( sizeof( double ) * n_vertices * DIM );
}

void allocate_cells( Cells_t *cells, int n_cells, int max_vertices, int max_faces )
{
    cells->id           = allocate( sizeof( int ) * n_cells );
    cells->type         = allocate( sizeof( int ) * n_cells );
    cells->n_vertices   = allocate( sizeof( int ) * n_cells );
    cells->vertices     = allocate( sizeof( int ) * n_cells * max_vertices );
    cells->n_faces      = allocate( sizeof( int ) * n_cells );
    cells->faces        = allocate( sizeof( int ) * n_cells * max_faces );
    cells->x            = allocate( sizeof( double ) * n_cells * DIM );
    cells->volume       = allocate( sizeof( double ) * n_cells );
    cells->dx           = allocate( sizeof( double ) * n_cells * DIM );

    set_value_int_n( 0, cells->n_vertices, n_cells );
    set_value_int_n( 0, cells->n_faces, n_cells );
}

void allocate_boundaries( Boundaries_t *boundaries, int n_boundaries, int max_vertices )
{
    boundaries->id          = allocate( sizeof( int ) * n_boundaries );
    boundaries->type        = allocate( sizeof( int ) * n_boundaries );
    boundaries->n_vertices  = allocate( sizeof( int ) * n_boundaries );
    boundaries->vertices    = allocate( sizeof( int ) * n_boundaries * max_vertices );
    boundaries->face        = allocate( sizeof( int ) * n_boundaries );
    boundaries->distance    = allocate( sizeof( double ) * n_boundaries );
    boundaries->n           = allocate( sizeof( double ) * n_boundaries * DIM );
    boundaries->t1          = allocate( sizeof( double ) * n_boundaries * DIM );
    boundaries->t2          = allocate( sizeof( double ) * n_boundaries * DIM );

    set_value_int_n( 0, boundaries->n_vertices, n_boundaries );
}

void allocate_faces( Faces_t *faces, int n_faces, int max_vertices )
{
    faces->type         = allocate( sizeof( int ) * n_faces );
    faces->n_vertices   = allocate( sizeof( int ) * n_faces );
    faces->vertices     = allocate( sizeof( int ) * n_faces * max_vertices );
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
}

void allocate_regions( Regions_t *regions, int n_regions, int max_name_length )
{
    regions->n_regions          = n_regions;
    regions->max_name_length    = max_name_length;

    regions->name           = allocate_hdf5_string_buffer( n_regions, max_name_length );
    regions->is_boundary    = allocate( sizeof( int ) * n_regions );
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

void deallocate_vertices( Vertices_t **vertices )
{
    if ((*vertices) == NULL) return;

    deallocate( (*vertices)->x );

    deallocate( *vertices );
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

void deallocate_regions( Regions_t **regions )
{
    if ((*regions) == NULL) return;

    deallocate_hdf5_string_buffer( &(*regions)->name );
    deallocate( (*regions)->is_boundary );

    deallocate( *regions );
}
