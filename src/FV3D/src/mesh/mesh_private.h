//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef MESH_PRIVATE_H
#define MESH_PRIVATE_H

#include "fv3d_module.h"

//##################################################################################################################################
// DEFINES
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// MACROS
//----------------------------------------------------------------------------------------------------------------------------------

//##################################################################################################################################
// VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------
    // character(len=_STRLEN_)             :: mesh_file = 'untitled.mesh.h5'   !< The mesh file

    // integer                             :: dimension = 0                    !< mesh dimensionality
    // logical                             :: is_partitioned = .false.         !< partitioned mesh flag

    // integer                             :: n_partitions = 0                 !< number of partitions
    // integer                             :: n_partition_cells = 0            !< the number of cells per partition
    // integer                             :: n_partition_boundaries = 0       !< the number of boundaries per partition
    // integer                             :: n_partition_faces = 0            !< the number of faces per partition
    // integer                             :: n_partition_sends = 0            !< the number of sends
    // integer                             :: n_partition_receives = 0         !< the number of recieves

    // integer,                allocatable :: partition_cells(:)               !< the partition cells
    // integer,                allocatable :: partition_boundaries(:)          !< the partition boundaries
    // integer,                allocatable :: partition_faces(:)               !< the partition faces
    // integer,                allocatable :: partition_sends(:)               !< the partition sends
    // integer,                allocatable :: partition_sends_pid(:)           !< the partition sends partition id
    // integer,                allocatable :: partition_receives(:)            !< the partition receives
    // integer,                allocatable :: partition_receives_pid(:)        !< the partition receives partition id

    // integer,                allocatable :: n_partition_sends_to(:)          !< number of sends to
    // integer,                allocatable :: partition_sends_to(:,:)          !< sends to list
    // integer,                allocatable :: n_partition_receives_from(:)     !< number of receives from
    // integer,                allocatable :: partition_receives_from(:,:)     !< receives from list

    // integer                             :: n_vertices = 0                   !< number of vertices
    // real,   allocatable                 :: vertices(:,:)                    !< vertices array

    // type cell_t
    //     integer                         :: id = 0                           !< cell region id
    //     integer                         :: type = 0                         !< cell element type
    //     integer                         :: n_vertices = 0                   !< number of cell vertices
    //     integer,    allocatable         :: vertices(:)                      !< cell vertices
    //     integer                         :: n_faces = 0                      !< number of cell faces
    //     integer,    allocatable         :: faces(:)                         !< cell faces
    //     real                            :: x(3) = 0.0                       !< cell centre
    //     real                            :: volume = 0.0                     !< cell volume
    //     real                            :: ds(3) = 0.0                      !< cell metrics
    // end type cell_t

    // integer                             :: n_global_cells = 0               !< number of global cells
    // integer                             :: n_cells = 0                      !< number of cells
    // type(cell_t),   allocatable         :: cells(:)                         !< cells array

    // type boundary_t
    //     integer                         :: id = 0                           !< boundary region id
    //     integer                         :: type = 0                         !< boundary element type
    //     integer                         :: n_vertices = 0                   !< number of boundary vertices
    //     integer,    allocatable         :: vertices(:)                      !< boundary vertices
    //     integer                         :: face = 0                         !< boundary faces
    //     real                            :: distance = 0.0                   !< boundary distance
    //     real                            :: n(3) = 0.0                       !< boundary normal vector
    //     real                            :: t1(3) = 0.0                      !< boundary tangential vector 1
    //     real                            :: t2(3) = 0.0                      !< boundary tangential vector 2
    // end type boundary_t

    // integer                             :: n_global_boundaries = 0          !< number of global boundaries
    // integer                             :: n_boundaries = 0                 !< number of boundaries
    // type(boundary_t),   allocatable     :: boundaries(:)                    !< boundaries array

    // type face_t
    //     integer                         :: type = 0                         !< face element type
    //     integer                         :: cells(2) = 0                     !< face cells
    //     integer                         :: boundary = 0                     !< face boundary
    //     integer                         :: n_vertices = 0                   !< number of vertices
    //     integer,    allocatable         :: vertices(:)                      !< face vertices
    //     real                            :: area = 0.0                       !< face area
    //     real                            :: lambda = 0.0                     !< face weighting
    //     real                            :: x(3) = 0.0                       !< face centre
    //     real                            :: n(3) = 0.0                       !< face normal vector
    //     real                            :: t1(3) = 0.0                      !< face tangential vector 1
    //     real                            :: t2(3) = 0.0                      !< face tangential vector 2
    // end type face_t

    // integer                             :: n_global_faces = 0               !< number of global faces
    // integer                             :: n_faces = 0                      !< number of faces
    // type(face_t),   allocatable         :: faces(:)                         !< faces array

    // type region_t
    //     character(len=:),   allocatable :: name                             !< region name
    //     logical                         :: is_boundary = .false.            !< region is boundary flag
    //     real,   allocatable             :: phi(:)                           !< region vector of variables
    //     integer                         :: type = 0                         !< region boundary type
    //     integer                         :: function_id = 0                  !< region boundary function id
    // end type region_t

    // integer                             :: n_regions = 0                    !< number of regions
    // type(region_t), allocatable         :: regions(:)                       !< regions array

    // real,   allocatable,    target      :: dist_cell_1(:,:)                 !< faces to cell 1 distance
    // real,   allocatable,    target      :: dist_cell_2(:,:)                 !< faces to cell 2 distance

    // real                                :: total_volume = 0.0               !< total mesh volume

    // integer                             :: flow_region = 0                  !< flow region identifier

    // integer                             :: max_cell_vertices = 0            !< maximum number of cell vertices
    // integer                             :: max_cell_faces = 0               !< maximum number of cell faces
    // integer                             :: max_boundary_vertices = 0        !< maximum number of boundary vertices
    // integer                             :: max_face_vertices = 0            !< maximum number of face vertices

    // integer                             :: n_internal_faces = 0             !< number of internal faces
    // integer,    allocatable             :: internal_faces(:)                !< internal faces list
    // integer                             :: n_boundary_faces = 0             !< number of boundary faces
    // integer,    allocatable             :: boundary_faces(:)                !< boundary faces list

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define();

#endif /* MESH_PRIVATE_H */