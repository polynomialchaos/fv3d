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
void deallocate_partition( Partition_t *partition );
void deallocate_vertices( Vertices_t *vertices );
void deallocate_cells( Cells_t *cells );
void deallocate_boundaries( Boundaries_t *boundaries );
void deallocate_faces( Faces_t *faces );
void deallocate_regions( Regions_t *regions );

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void free_mesh( Mesh_t *mesh )
{
    deallocate_partition( &mesh->partition );
    deallocate_vertices( &mesh->vertices );
    deallocate_cells( &mesh->cells );
    deallocate_boundaries( &mesh->boundaries );
    deallocate_faces( &mesh->faces );
    deallocate_regions( &mesh->regions );
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
}

void allocate_vertices( Vertices_t *vertices, int n_vertices )
{
    vertices->x = allocate( sizeof( double ) * n_vertices * 3 );
}

void allocate_cells( Cells_t *cells, int n_cells, int max_vertices, int max_faces )
{
    cells->id           = allocate( sizeof( int ) * n_cells );
    cells->type         = allocate( sizeof( int ) * n_cells );
    cells->n_vertices   = allocate( sizeof( int ) * n_cells );
    cells->vertices     = allocate( sizeof( int ) * n_cells * max_vertices );
    cells->n_faces      = allocate( sizeof( int ) * n_cells );
    cells->faces        = allocate( sizeof( int ) * n_cells * max_faces );
    cells->x            = allocate( sizeof( double ) * n_cells * 3 );
    cells->volume       = allocate( sizeof( double ) * n_cells );
    cells->ds           = allocate( sizeof( double ) * n_cells * 3 );
}

void allocate_boundaries( Boundaries_t *boundaries, int n_boundaries, int max_vertices )
{
    boundaries->id          = allocate( sizeof( int ) * n_boundaries );
    boundaries->type        = allocate( sizeof( int ) * n_boundaries );
    boundaries->n_vertices  = allocate( sizeof( int ) * n_boundaries );
    boundaries->vertices    = allocate( sizeof( int ) * n_boundaries * max_vertices );
    boundaries->face        = allocate( sizeof( int ) * n_boundaries );
    boundaries->distance    = allocate( sizeof( double ) * n_boundaries );
    boundaries->n           = allocate( sizeof( double ) * n_boundaries * 3 );
    boundaries->t1          = allocate( sizeof( double ) * n_boundaries * 3 );
    boundaries->t2          = allocate( sizeof( double ) * n_boundaries * 3 );
}

void allocate_faces( Faces_t *faces, int n_faces, int max_vertices )
{
    faces->type         = allocate( sizeof( int ) * n_faces );
    faces->n_vertices   = allocate( sizeof( int ) * n_faces );
    faces->vertices     = allocate( sizeof( int ) * n_faces * max_vertices );
    faces->cells        = allocate( sizeof( int ) * n_faces * 2 );
    faces->boundary     = allocate( sizeof( int ) * n_faces );
    faces->area         = allocate( sizeof( double ) * n_faces );
    faces->lambda       = allocate( sizeof( double ) * n_faces );
    faces->x            = allocate( sizeof( double ) * n_faces * 3 );
    faces->n            = allocate( sizeof( double ) * n_faces * 3 );
    faces->t1           = allocate( sizeof( double ) * n_faces * 3 );
    faces->t2           = allocate( sizeof( double ) * n_faces * 3 );
}

void allocate_regions( Regions_t *regions, int n_regions, int length )
{
    regions->name = allocate_hdf5_string_buffer( n_regions, length );
}

void deallocate_partition( Partition_t *partition )
{
    deallocate( partition->partition_cells );
    deallocate( partition->partition_boundaries );
    deallocate( partition->partition_faces );
    deallocate( partition->partition_sends );
    deallocate( partition->partition_sends_pid );
    deallocate( partition->partition_receives );
    deallocate( partition->partition_receives_pid );

    deallocate( partition->n_partition_sends_to );
    deallocate( partition->partition_sends_to );
    deallocate( partition->n_partition_receives_from );
    deallocate( partition->partition_receives_from );
}

void deallocate_vertices( Vertices_t *vertices )
{
    deallocate( vertices->x );
}

void deallocate_cells( Cells_t *cells )
{
    deallocate( cells->id );
    deallocate( cells->type );
    deallocate( cells->n_vertices );
    deallocate( cells->vertices );
    deallocate( cells->n_faces );
    deallocate( cells->faces );
    deallocate( cells->x );
    deallocate( cells->volume );
    deallocate( cells->ds );
}

void deallocate_boundaries( Boundaries_t *boundaries )
{
    deallocate( boundaries->id );
    deallocate( boundaries->type );
    deallocate( boundaries->n_vertices );
    deallocate( boundaries->vertices );
    deallocate( boundaries->face );
    deallocate( boundaries->distance );
    deallocate( boundaries->n );
    deallocate( boundaries->t1 );
    deallocate( boundaries->t2 );
}

void deallocate_faces( Faces_t *faces )
{
    deallocate( faces->type );
    deallocate( faces->n_vertices );
    deallocate( faces->vertices );
    deallocate( faces->cells );
    deallocate( faces->boundary );
    deallocate( faces->area );
    deallocate( faces->lambda );
    deallocate( faces->x );
    deallocate( faces->n );
    deallocate( faces->t1 );
    deallocate( faces->t2 );
}

void deallocate_regions( Regions_t *regions )
{
    deallocate_hdf5_string_buffer( &regions->name );
}
