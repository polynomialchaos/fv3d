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

    string_t tmp = "untitled_mesh.h5";
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

        //     _DEALLOCATE( partition_cells )
        // _DEALLOCATE( partition_boundaries )
        // _DEALLOCATE( partition_faces )
        // _DEALLOCATE( partition_sends )
        // _DEALLOCATE( partition_sends_pid )
        // _DEALLOCATE( partition_receives )
        // _DEALLOCATE( partition_receives_pid )

        // _DEALLOCATE( n_partition_sends_to )
        // _DEALLOCATE( partition_sends_to )
        // _DEALLOCATE( n_partition_receives_from )
        // _DEALLOCATE( partition_receives_from )

        // _DEALLOCATE( vertices )
        // _DEALLOCATE( cells )
        // _DEALLOCATE( faces )
        // _DEALLOCATE( boundaries )
        // _DEALLOCATE( regions )

        // _DEALLOCATE( dist_cell_1 )
        // _DEALLOCATE( dist_cell_2 )

        // _DEALLOCATE( internal_faces )
        // _DEALLOCATE( boundary_faces )
}

void read_mesh_file()
{
        // integer                             :: max_cell_vertices = 0            !< maximum number of cell vertices
    // integer                             :: max_cell_faces = 0               !< maximum number of cell faces
    // integer                             :: max_boundary_vertices = 0        !< maximum number of boundary vertices
    // integer                             :: max_face_vertices = 0            !< maximum number of face vertices

        // is_parallel = get_is_parallel()
        // i_rank = get_i_rank()

        // file_id = open_hdf5_file( mesh_file )
        //     ! the array sizes for visualization
        //     call get_hdf5_attribute( file_id, 'dim', dimension )
        //     call get_hdf5_attribute( file_id, 'is_partitioned', is_partitioned )

        //     if( (is_parallel) .and. (.not. is_partitioned) ) &
        //         call add_error(__LINE__, __FILE__, &
        //             'Parallel run requires a partitioned mesh!' )

        //     ! the partition
        //     if( is_parallel ) then
        //         group_id = open_hdf5_group( file_id, 'PARTITIONS' )
        //             call get_hdf5_attribute( group_id, 'n_partitions', n_partitions )

        //             if( get_n_procs() .ne. n_partitions ) &
        //                 call add_error(__LINE__, __FILE__, &
        //                     'Parallel run requires a specific amount of processors (req: ' &
        //                         // set_string( get_n_procs() ) // ', got: ' // set_string( n_partitions ) // ')!' )

        //             call get_hdf5_dataset( group_id, 'n_partition_cells', tmp_si, offset=[i_rank] )
        //             n_partition_cells = tmp_si(1)
        //             call get_hdf5_dataset( group_id, 'n_partition_boundaries', tmp_si, offset=[i_rank] )
        //             n_partition_boundaries = tmp_si(1)
        //             call get_hdf5_dataset( group_id, 'n_partition_faces', tmp_si, offset=[i_rank] )
        //             n_partition_faces = tmp_si(1)
        //             call get_hdf5_dataset( group_id, 'n_partition_sends', tmp_si, offset=[i_rank] )
        //             n_partition_sends = tmp_si(1)
        //             call get_hdf5_dataset( group_id, 'n_partition_receives', tmp_si, offset=[i_rank] )
        //             n_partition_receives = tmp_si(1)

        //             allocate( partition_cells(n_partition_cells) )
        //             allocate( partition_boundaries(n_partition_boundaries) )
        //             allocate( partition_faces(n_partition_faces) )
        //             allocate( partition_sends(n_partition_sends) )
        //             allocate( partition_sends_pid(n_partition_sends) )
        //             allocate( partition_receives(n_partition_receives) )
        //             allocate( partition_receives_pid(n_partition_receives) )

        //             call get_hdf5_dataset( group_id, 'partition_cells', partition_cells, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_boundaries', partition_boundaries, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_faces', partition_faces, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_sends', partition_sends, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_sends_pid', partition_sends_pid, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_receives', partition_receives, offset=[0,i_rank] )
        //             call get_hdf5_dataset( group_id, 'partition_receives_pid', partition_receives_pid, offset=[0,i_rank] )
        //         call close_hdf5_group( group_id, 'PARTITIONS' )
        //     end if

        //     ! the vertices
        //     group_id = open_hdf5_group( file_id, 'VERTICES' )
        //         call get_hdf5_attribute( group_id, 'n_vertices', n_vertices )

        //         allocate( vertices(3,n_vertices) )

        //         call get_hdf5_dataset( group_id, 'x', vertices )
        //     call close_hdf5_group( group_id, 'VERTICES' )

        //     ! the cells
        //     group_id = open_hdf5_group( file_id, 'CELLS' )
        //         call get_hdf5_attribute( group_id, 'n_cells', n_cells ); n_global_cells = n_cells
        //         call get_hdf5_attribute( group_id, 'max_cell_vertices', max_cell_vertices )
        //         call get_hdf5_attribute( group_id, 'max_cell_faces', max_cell_faces )
        //         if( is_parallel ) n_cells = n_partition_cells

        //         allocate( cells(n_cells+n_partition_receives) )

        //         if( is_parallel ) then
        //             allocate( tmp_stride(n_cells+n_partition_receives) )
        //             tmp_stride(:n_cells) = partition_cells
        //             tmp_stride(n_cells+1:) = partition_receives

        //             call get_hdf5_dataset( group_id, 'id', cells(:)%id, stride=tmp_stride )
        //             call get_hdf5_dataset( group_id, 'type', cells(:)%type, stride=tmp_stride )

        //             call get_hdf5_dataset( group_id, 'n_vertices', cells(:)%n_vertices, stride=tmp_stride )
        //             allocate( tmp_i(max_cell_vertices,n_cells+n_partition_receives) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%vertices = tmp_i(1:cells(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'n_faces', cells(:)%n_faces, stride=tmp_stride )
        //             allocate( tmp_i(max_cell_faces,n_cells+n_partition_receives) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'faces', tmp_i, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%faces = tmp_i(1:cells(i)%n_faces,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             allocate( tmp_r(3,n_cells+n_partition_receives) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'x', tmp_r, stride=tmp_stride )
        //             do i = 1, n_cells+n_partition_receives
        //                 cells(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             call get_hdf5_dataset( group_id, 'volume', cells(:)%volume, stride=tmp_stride )

        //             _DEALLOCATE(tmp_stride)
        //         else
        //             call get_hdf5_dataset( group_id, 'id', cells(:)%id )
        //             call get_hdf5_dataset( group_id, 'type', cells(:)%type )

        //             call get_hdf5_dataset( group_id, 'n_vertices', cells(:)%n_vertices )
        //             allocate( tmp_i(max_cell_vertices,n_cells) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i )
        //             do i = 1, n_cells
        //                 cells(i)%vertices = tmp_i(1:cells(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'n_faces', cells(:)%n_faces )
        //             allocate( tmp_i(max_cell_faces,n_cells) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'faces', tmp_i )
        //             do i = 1, n_cells
        //                 cells(i)%faces = tmp_i(1:cells(i)%n_faces,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             allocate( tmp_r(3,n_cells) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'x', tmp_r )
        //             do i = 1, n_cells
        //                 cells(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             call get_hdf5_dataset( group_id, 'volume', cells(:)%volume )
        //         end if
        //     call close_hdf5_group( group_id, 'CELLS' )

        //     ! the boundaries
        //     group_id = open_hdf5_group( file_id, 'BOUNDARIES' )
        //         call get_hdf5_attribute( group_id, 'n_boundaries', n_boundaries ); n_global_boundaries = n_boundaries
        //         call get_hdf5_attribute( group_id, 'max_boundary_vertices', max_boundary_vertices )
        //         if( is_parallel ) n_boundaries = n_partition_boundaries

        //         allocate( boundaries(n_boundaries) )

        //         if( is_parallel ) then
        //             call get_hdf5_dataset( group_id, 'id', boundaries(:)%id, stride=partition_boundaries )
        //             call get_hdf5_dataset( group_id, 'type', boundaries(:)%type, stride=partition_boundaries )

        //             call get_hdf5_dataset( group_id, 'n_vertices', boundaries(:)%n_vertices, stride=partition_boundaries )
        //             allocate( tmp_i(max_boundary_vertices,n_boundaries) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%vertices = tmp_i(1:boundaries(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'face', boundaries(:)%face, stride=partition_boundaries )
        //             call get_hdf5_dataset( group_id, 'distance', boundaries(:)%distance, stride=partition_boundaries )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'n', tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't1', tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't2', tmp_r, stride=partition_boundaries )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         else
        //             call get_hdf5_dataset( group_id, 'id', boundaries(:)%id )
        //             call get_hdf5_dataset( group_id, 'type', boundaries(:)%type )

        //             call get_hdf5_dataset( group_id, 'n_vertices', boundaries(:)%n_vertices )
        //             allocate( tmp_i(max_boundary_vertices,n_boundaries) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%vertices = tmp_i(1:boundaries(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'face', boundaries(:)%face )
        //             call get_hdf5_dataset( group_id, 'distance', boundaries(:)%distance )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'n', tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't1', tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_boundaries) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't2', tmp_r )
        //             do i = 1, n_boundaries
        //                 boundaries(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         end if
        //     call close_hdf5_group( group_id,'BOUNDARIES' )

        //     ! the faces
        //     group_id = open_hdf5_group( file_id, 'FACES' )
        //         call get_hdf5_attribute( group_id, 'n_faces', n_faces ); n_global_faces = n_faces
        //         call get_hdf5_attribute( group_id, 'max_face_vertices', max_face_vertices )
        //         if( is_parallel ) n_faces = n_partition_faces

        //         allocate( faces(n_faces) )

        //         if( is_parallel ) then
        //             call get_hdf5_dataset( group_id, 'type', faces(:)%type, stride=partition_faces )
        //             call get_hdf5_dataset( group_id, 'cell1', faces(:)%cells(1), stride=partition_faces )
        //             call get_hdf5_dataset( group_id, 'cell2', faces(:)%cells(2), stride=partition_faces )
        //             call get_hdf5_dataset( group_id, 'boundary', faces(:)%boundary, stride=partition_faces )

        //             call get_hdf5_dataset( group_id, 'n_vertices', faces(:)%n_vertices, stride=partition_faces )
        //             allocate( tmp_i(max_face_vertices,n_faces) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i, stride=partition_faces )
        //             do i = 1, n_cells
        //                 faces(i)%vertices = tmp_i(1:faces(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'area', faces(:)%area, stride=partition_faces )
        //             call get_hdf5_dataset( group_id, 'lambda', faces(:)%lambda, stride=partition_faces )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'x', tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'n', tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't1', tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't2', tmp_r, stride=partition_faces )
        //             do i = 1, n_faces
        //                 faces(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         else
        //             call get_hdf5_dataset( group_id, 'type', faces(:)%type )
        //             call get_hdf5_dataset( group_id, 'cell1', faces(:)%cells(1) )
        //             call get_hdf5_dataset( group_id, 'cell2', faces(:)%cells(2) )
        //             call get_hdf5_dataset( group_id, 'boundary', faces(:)%boundary )

        //             call get_hdf5_dataset( group_id, 'n_vertices', faces(:)%n_vertices )
        //             allocate( tmp_i(max_face_vertices,n_faces) ); tmp_i = 0
        //             call get_hdf5_dataset( group_id, 'vertices', tmp_i )
        //             do i = 1, n_cells
        //                 faces(i)%vertices = tmp_i(1:faces(i)%n_vertices,i)
        //             end do
        //             _DEALLOCATE( tmp_i )

        //             call get_hdf5_dataset( group_id, 'area', faces(:)%area )
        //             call get_hdf5_dataset( group_id, 'lambda', faces(:)%lambda )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'x', tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%x = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 'n', tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%n = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't1', tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%t1 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )

        //             allocate( tmp_r(3,n_faces) ); tmp_r = 0
        //             call get_hdf5_dataset( group_id, 't2', tmp_r )
        //             do i = 1, n_faces
        //                 faces(i)%t2 = tmp_r(:,i)
        //             end do
        //             _DEALLOCATE( tmp_r )
        //         end if
        //     call close_hdf5_group( group_id, 'FACES' )

        //     ! the regions
        //     group_id = open_hdf5_group( file_id, 'REGIONS' )
        //         call get_hdf5_attribute( group_id, 'n_regions', n_regions )

        //         allocate( regions(n_regions) )

        //         allocate( tmp_c(n_regions) ); tmp_c = ''
        //         call get_hdf5_dataset( group_id, 'name', tmp_c )
        //         do i = 1, n_regions
        //             regions(i)%name = tmp_c(i)
        //         end do
        //         _DEALLOCATE( tmp_c )
        //     call close_hdf5_group( group_id, 'REGIONS' )
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