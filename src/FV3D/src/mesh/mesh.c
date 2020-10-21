//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
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
string_t mesh_file = NULL;

int dimension = 0;
int is_partitioned = 0;

int n_partitions = 0;
int n_partition_cells = 0;
int n_partition_boundaries = 0;
int n_partition_faces = 0;
int n_partition_sends = 0;
int n_partition_receives = 0;

int *partition_cells = NULL;
int *partition_boundaries = NULL;
int *partition_faces = NULL;
int *partition_sends = NULL;
int *partition_sends_pid = NULL;
int *partition_receives = NULL;
int *partition_receives_pid = NULL;

int *n_partition_sends_to = NULL;
int *partition_sends_to = NULL;
int *n_partition_receives_from = NULL;
int *partition_receives_from = NULL;

int n_vertices = 0;
double *vertices = NULL;

int n_global_cells = 0;
int n_cells = 0;
Cell_t *cells = NULL;

int n_global_boundaries = 0;
int n_boundaries = 0;
Boundary_t *boundaries = NULL;

int n_global_faces = 0;
int n_faces = 0;
Face_t *faces = NULL;

int n_regions = 0;
Region_t *regions = NULL;

double *dist_cell_1 = NULL;
double *dist_cell_2 = NULL;

double total_volume = 0.0;

int flow_region = 0;

int n_internal_faces = 0;
int *internal_faces = NULL;
int n_boundary_faces = 0;
int *boundary_faces = NULL;

//##################################################################################################################################
// LOCAL FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_initialize();
void mesh_finalize();
void read_mesh_file();
void calc_mesh_metrics();

//##################################################################################################################################
// FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------
void mesh_define()
{
    register_initialize_routine( mesh_initialize );
    register_finalize_routine( mesh_finalize );

    string_t tmp = "untitled.mesh.h5";
    set_parameter( "Mesh/mesh_file", ParameterString, &tmp, "The mesh file", NULL, 0 );
}

void mesh_initialize()
{
    get_parameter( "Mesh/mesh_file", ParameterString, &mesh_file );

    read_mesh_file();
    calc_mesh_metrics();
}

void mesh_finalize()
{
    deallocate( mesh_file );

    deallocate( partition_cells );
    deallocate( partition_boundaries );
    deallocate( partition_faces );
    deallocate( partition_sends );
    deallocate( partition_sends_pid );
    deallocate( partition_receives );
    deallocate( partition_receives_pid );

    deallocate( n_partition_sends_to );
    deallocate( partition_sends_to );
    deallocate( n_partition_receives_from );
    deallocate( partition_receives_from );

    deallocate( vertices );

    deallocate( cells );
    deallocate( boundaries );
    deallocate( faces );
    deallocate( regions );

    deallocate( dist_cell_1 );
    deallocate( dist_cell_2 );

    deallocate( internal_faces );
    deallocate( boundary_faces );
}

void read_mesh_file()
{
    int is_parallel = get_is_parallel();
    int rank = get_rank();

    hsize_t dims[2];
    hsize_t offset[2];
    hsize_t count[2];

    hid_t file_id = open_hdf5_file( mesh_file );
        get_hdf5_attribute( file_id, "dimension", HDF5Int, &dimension );
        get_hdf5_attribute( file_id, "is_partitioned", HDF5Int, &is_partitioned );

        // the partition
        if (is_parallel)
        {
            check_error( (is_partitioned == 1) );

            hid_t group_id = open_hdf5_group( file_id, "PARTITIONS" );
                get_hdf5_attribute( group_id, "n_partitions", HDF5Int, &n_partitions );
                check_error( (n_partitions == get_n_procs()) );

                dims[0]     = 1;
                offset[0]   = rank;
                count[0]    = 1;

                get_hdf5_dataset_chunk_n( group_id, "n_partition_cells", HDF5Int,
                    &n_partition_cells, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_boundaries", HDF5Int,
                    &n_partition_boundaries, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_faces", HDF5Int,
                    &n_partition_faces, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_sends", HDF5Int,
                    &n_partition_sends, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n( group_id, "n_partition_receives", HDF5Int,
                    &n_partition_receives, 1, dims, NULL, offset, count, NULL, NULL );

                partition_cells         = allocate( sizeof( int ) * n_partition_cells );
                partition_boundaries    = allocate( sizeof( int ) * n_partition_boundaries );
                partition_faces         = allocate( sizeof( int ) * n_partition_faces );
                partition_sends         = allocate( sizeof( int ) * n_partition_sends );
                partition_sends_pid     = allocate( sizeof( int ) * n_partition_sends );
                partition_receives      = allocate( sizeof( int ) * n_partition_receives );
                partition_receives_pid  = allocate( sizeof( int ) * n_partition_receives );

                dims[0] = n_partition_cells; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_cells;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_cells", HDF5Int,
                    partition_cells, 1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_boundaries; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_boundaries;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_boundaries", HDF5Int,
                    partition_boundaries, 1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_faces; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_faces;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_faces", HDF5Int,
                    partition_faces, 1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_sends; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_sends;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends", HDF5Int,
                    partition_sends, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_sends_pid", HDF5Int,
                    partition_sends_pid, 1, dims, NULL, offset, count, NULL, NULL );

                dims[0] = n_partition_receives; offset[0] = rank;  count[0] = 1;
                dims[1] = 0; offset[1] = 0; count[1] = n_partition_receives;
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives", HDF5Int,
                    partition_receives, 1, dims, NULL, offset, count, NULL, NULL );
                get_hdf5_dataset_chunk_n_m( group_id, "partition_receives_pid", HDF5Int,
                    partition_receives_pid, 1, dims, NULL, offset, count, NULL, NULL );
            close_hdf5_group( group_id );
        }

        // the vertices
        hid_t group_id = open_hdf5_group( file_id, "VERTICES" );
            get_hdf5_attribute( group_id, "n_vertices", HDF5Int, &n_vertices );

            vertices = allocate( sizeof( double ) * n_vertices );
            dims[0] = 3 * n_vertices;
            dims[1] = 0;
            get_hdf5_dataset_chunk_n_m( group_id, "x", HDF5Int,
                vertices, 1, dims, NULL, NULL, NULL, NULL, NULL );
        close_hdf5_group( group_id );

    close_hdf5_file( file_id );
        // integer                             :: max_cell_vertices = 0            !< maximum number of cell vertices
    // integer                             :: max_cell_faces = 0               !< maximum number of cell faces
    // integer                             :: max_boundary_vertices = 0        !< maximum number of boundary vertices
    // integer                             :: max_face_vertices = 0            !< maximum number of face vertices


        //     ! the cells
        //     group_id = open_hdf5_group( file_id, "CELLS" )
        //         call get_hdf5_attribute( group_id, "n_cells", n_cells ); n_global_cells = n_cells
        //         call get_hdf5_attribute( group_id, "max_cell_vertices", max_cell_vertices )
        //         call get_hdf5_attribute( group_id, "max_cell_faces", max_cell_faces )
        //         if( is_parallel ) n_cells = n_partition_cells

        //         allocate( cells(n_cells+n_partition_receives) )

        //         if( is_parallel ) then
        //             allocate( tmp_stride(n_cells+n_partition_receives) )
        //             tmp_stride(:n_cells) = partition_cells
        //             tmp_stride(n_cells+1:) = partition_receives

        //             call get_hdf5_dataset( group_id, "id", cells(:)%id, stride=tmp_stride )
        //             call get_hdf5_dataset( group_id, "type", cells(:)%type, stride=tmp_stride )

        //             call get_hdf5_dataset( group_id, "n_vertices", cells(:)%n_vertices, stride=tmp_stride )
        //             allocate( tmp_i(max_cell_vertices,n_cells+n_partition_receives) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%vertices = tmp_i(1:cells(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "n_faces", cells(:)%n_faces, stride=tmp_stride )
        //             allocate( tmp_i(max_cell_faces,n_cells+n_partition_receives) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "faces", tmp_i, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%faces = tmp_i(1:cells(i)%n_faces,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             allocate( tmp_r(3,n_cells+n_partition_receives) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "x", tmp_r, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             call get_hdf5_dataset( group_id, "volume", cells(:)%volume, stride=tmp_stride )

        //             _DEALLOCATE(tmp_stride)
        //         else
        //             call get_hdf5_dataset( group_id, "id", cells(:)%id )
        //             call get_hdf5_dataset( group_id, "type", cells(:)%type )

        //             call get_hdf5_dataset( group_id, "n_vertices", cells(:)%n_vertices )
        //             allocate( tmp_i(max_cell_vertices,n_cells) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i )
        //             do i = 1, n_cells
        //                 cells(i)%vertices = tmp_i(1:cells(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "n_faces", cells(:)%n_faces )
        //             allocate( tmp_i(max_cell_faces,n_cells) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "faces", tmp_i )
        //             do i = 1, n_cells
        //                 cells(i)%faces = tmp_i(1:cells(i)%n_faces,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             allocate( tmp_r(3,n_cells) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "x", tmp_r )
        //             do i = 1, n_cells
        //                 cells(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             call get_hdf5_dataset( group_id, "volume", cells(:)%volume )
        //         end if
        //     call close_hdf5_group( group_id, "CELLS" )

        //     ! the boundaries
        //     group_id = open_hdf5_group( file_id, "BOUNDARIES" )
        //         call get_hdf5_attribute( group_id, "n_boundaries", n_boundaries ); n_global_boundaries = n_boundaries
        //         call get_hdf5_attribute( group_id, "max_boundary_vertices", max_boundary_vertices )
        //         if( is_parallel ) n_boundaries = n_partition_boundaries

        //         allocate( boundaries(n_boundaries) )

        //         if( is_parallel ) then
        //             call get_hdf5_dataset( group_id, "id", boundaries(:)%id, stride=partition_boundaries )
        //             call get_hdf5_dataset( group_id, "type", boundaries(:)%type, stride=partition_boundaries )

        //             call get_hdf5_dataset( group_id, "n_vertices", boundaries(:)%n_vertices, stride=partition_boundaries )
        //             allocate( tmp_i(max_boundary_vertices,n_boundaries) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%vertices = tmp_i(1:boundaries(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "face", boundaries(:)%face, stride=partition_boundaries )
        //             call get_hdf5_dataset( group_id, "distance", boundaries(:)%distance, stride=partition_boundaries )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "n", tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t1", tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t2", tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         else
        //             call get_hdf5_dataset( group_id, "id", boundaries(:)%id )
        //             call get_hdf5_dataset( group_id, "type", boundaries(:)%type )

        //             call get_hdf5_dataset( group_id, "n_vertices", boundaries(:)%n_vertices )
        //             allocate( tmp_i(max_boundary_vertices,n_boundaries) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%vertices = tmp_i(1:boundaries(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "face", boundaries(:)%face )
        //             call get_hdf5_dataset( group_id, "distance", boundaries(:)%distance )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "n", tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t1", tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t2", tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         end if
        //     call close_hdf5_group( group_id,"BOUNDARIES" )

        //     ! the faces
        //     group_id = open_hdf5_group( file_id, "FACES" )
        //         call get_hdf5_attribute( group_id, "n_faces", n_faces ); n_global_faces = n_faces
        //         call get_hdf5_attribute( group_id, "max_face_vertices", max_face_vertices )
        //         if( is_parallel ) n_faces = n_partition_faces

        //         allocate( faces(n_faces) )

        //         if( is_parallel ) then
        //             call get_hdf5_dataset( group_id, "type", faces(:)%type, stride=partition_faces )
        //             call get_hdf5_dataset( group_id, "cell1", faces(:)%cells(1), stride=partition_faces )
        //             call get_hdf5_dataset( group_id, "cell2", faces(:)%cells(2), stride=partition_faces )
        //             call get_hdf5_dataset( group_id, "boundary", faces(:)%boundary, stride=partition_faces )

        //             call get_hdf5_dataset( group_id, "n_vertices", faces(:)%n_vertices, stride=partition_faces )
        //             allocate( tmp_i(max_face_vertices,n_faces) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i, stride=partition_faces )
        //             do i = 1, n_cells
        //                 faces(i)%vertices = tmp_i(1:faces(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "area", faces(:)%area, stride=partition_faces )
        //             call get_hdf5_dataset( group_id, "lambda", faces(:)%lambda, stride=partition_faces )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "x", tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "n", tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t1", tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t2", tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         else
        //             call get_hdf5_dataset( group_id, "type", faces(:)%type )
        //             call get_hdf5_dataset( group_id, "cell1", faces(:)%cells(1) )
        //             call get_hdf5_dataset( group_id, "cell2", faces(:)%cells(2) )
        //             call get_hdf5_dataset( group_id, "boundary", faces(:)%boundary )

        //             call get_hdf5_dataset( group_id, "n_vertices", faces(:)%n_vertices )
        //             allocate( tmp_i(max_face_vertices,n_faces) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, "vertices", tmp_i )
        //             do i = 1, n_cells
        //                 faces(i)%vertices = tmp_i(1:faces(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, "area", faces(:)%area )
        //             call get_hdf5_dataset( group_id, "lambda", faces(:)%lambda )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "x", tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "n", tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t1", tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, "t2", tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         end if
        //     call close_hdf5_group( group_id, "FACES" )

        //     ! the regions
        //     group_id = open_hdf5_group( file_id, "REGIONS" )
        //         call get_hdf5_attribute( group_id, "n_regions", n_regions )

        //         allocate( regions(n_regions) )

        //         allocate( tmp_c(n_regions) ); tmp_c = ""
        //         call get_hdf5_dataset( group_id, "name", tmp_c )
        //         do i = 1, n_regions
        //             regions(i)%name = tmp_c(i)
        //         end do
        //         _DEALLOCATE( tmp_c )
        //     call close_hdf5_group( group_id, "REGIONS" )
        // call close_hdf5_file( file_id, mesh_file )

        // if( is_parallel ) call remap_local_mesh()
}

void remap_local_mesh()
{
        // global_cells        = 0
        // global_boundaries   = 0
        // global_faces        = 0

        // ! define the global index arrays
        // do i = 1, n_cells
        //     global_cells(partition_cells(i)) = i
        // end do

        // do i = 1, n_partition_receives
        //     global_cells(partition_receives(i)) = n_cells + i
        // end do

        // do i = 1, n_partition_sends
        //     partition_sends(i) = global_cells(partition_sends(i))
        // end do

        // do i = 1, n_partition_receives
        //     partition_receives(i) = global_cells(partition_receives(i))
        // end do

        // do i = 1, n_boundaries
        //     global_boundaries(partition_boundaries(i)) = i
        // end do

        // do i = 1, n_faces
        //     global_faces(partition_faces(i)) = i
        // end do

        // ! set the local cell indices
        // do i = 1, n_cells + n_partition_receives
        //     do j = 1, cells(i)%n_faces
        //         cells(i)%faces(j) = global_faces(cells(i)%faces(j))
        //     end do
        // end do

        // ! set the local boundary indices
        // do i = 1, n_boundaries
        //     boundaries(i)%face = global_faces(boundaries(i)%face)
        // end do

        // ! set the local face cell indices
        // do i = 1, n_faces
        //     faces(i)%cells(1) = global_cells(faces(i)%cells(1))

        //     if( faces(i)%cells(2) .gt. 0 ) then
        //         faces(i)%cells(2) = global_cells(faces(i)%cells(2))
        //     else
        //         faces(i)%boundary = global_boundaries(faces(i)%boundary)
        //     end if
        // end do

        // ! set the sends to and receives from list
        // allocate( n_partition_sends_to(n_partitions) ); n_partition_sends_to = 0
        // allocate( partition_sends_to(n_partition_sends,n_partitions) ); partition_sends_to = 0
        // allocate( n_partition_receives_from(n_partitions) ); n_partition_receives_from = 0
        // allocate( partition_receives_from(n_partition_receives,n_partitions) ); partition_receives_from = 0

        // do i = 1, n_partitions
        //     if( i-1 .eq. get_i_rank() ) cycle

        //     do j = 1, n_partition_sends
        //         if( partition_sends_pid(j) .ne. i ) cycle
        //         n_partition_sends_to(i) = n_partition_sends_to(i) + 1
        //         partition_sends_to(n_partition_sends_to(i),i) = partition_sends(j)
        //     end do

        //     do j = 1, n_partition_receives
        //         if( partition_receives_pid(j) .ne. i ) cycle
        //         n_partition_receives_from(i) = n_partition_receives_from(i) + 1
        //         partition_receives_from(n_partition_receives_from(i),i) = partition_receives(j)
        //     end do
        // end do
}

void calc_mesh_metrics()
{
        // total_volume_loc = sum( cells(:)%volume )
        // call allreduce_mpi( total_volume_loc, total_volume, MPI_SUM )

        // allocate( dist_cell_1(3,n_faces) ); dist_cell_1 = 0.0
        // allocate( dist_cell_2(3,n_faces) ); dist_cell_2 = 0.0

        // n_internal_faces = count( faces(:)%boundary .eq. 0 )
        // n_boundary_faces = n_faces - n_internal_faces
        // allocate(internal_faces(n_internal_faces)); internal_faces = 0
        // allocate(boundary_faces(n_boundary_faces)); boundary_faces = 0

        // ! define the cell metrics
        // do i = 1, n_cells
        //     do j = 1, cells(i)%n_faces
        //         cells(i)%ds = cells(i)%ds + &
        //             abs( faces(cells(i)%faces(j))%n ) * faces(cells(i)%faces(j))%area
        //     end do
        // end do

        // ! communication faces need to be updated with valid faces
        // do i = 1, n_faces
        //     fc = faces(i)%cells

        //     if( fc(1) .gt. n_cells ) then
        //         cells(fc(1))%n_faces = count( cells(fc(1))%faces .gt. 0 )
        //         cells(fc(1))%faces(1:cells(fc(1))%n_faces-1) = &
        //             pack( cells(fc(1))%faces, cells(fc(1))%faces .gt. 0 )
        //     endif

        //     if( fc(2) .gt. n_cells ) then
        //         cells(fc(2))%n_faces = count( cells(fc(2))%faces .gt. 0 )
        //         cells(fc(2))%faces(1:cells(fc(1))%n_faces-1) = &
        //             pack( cells(fc(2))%faces, cells(fc(2))%faces .gt. 0 )
        //     endif
        // end do

        // iif = 0; ibf = 0
        // do i = 1, n_faces
        //     fc = faces(i)%cells

        //     dist_cell_1(:,i) = faces(i)%x - cells(fc(1))%x

        //     if( fc(2) .gt. 0 ) then
        //         dist_cell_2(:,i)    = faces(i)%x - cells(fc(2))%x
        //         iif                 = iif + 1
        //         internal_faces(iif) = i
        //     else
        //         dist_cell_2(:,i)    = 0.0
        //         faces(i)%cells(2)   = n_cells + n_partition_receives + faces(i)%boundary
        //         ibf                 = ibf + 1
        //         boundary_faces(ibf) = i
        //     end if
        // end do

        // do i = 1, n_boundaries
        //     regions(boundaries(i)%id)%is_boundary = .true.
        // end do

        // do i = 1, n_regions
        //     if( regions(i)%is_boundary ) then
        //         flow_region = i
        //         exit
        //     end if
        // end do
}