/*******************************************************************************
 * @file mesh_private.h
 * @author Florian Eigentler
 * @brief
 * @version 1.0.0
 * @date 2021-11-15
 * @copyright Copyright (c) 2021
 ******************************************************************************/
#ifndef MESH_PRIVATE_H
#define MESH_PRIVATE_H

#include "solver/solver_module.h"

/*******************************************************************************
 * @brief Allocate the boundaries
 * @param mesh
 * @param n_boundaries
 * @param max_boundary_vertices
 * @return Boundaries_t*
 ******************************************************************************/
Boundaries_t *allocate_boundaries(Mesh_t *mesh, int n_boundaries,
                                  int max_boundary_vertices);

/*******************************************************************************
 * @brief Allocate the cells
 * @param mesh
 * @param n_cells
 * @param max_cell_vertices
 * @param max_cell_faces
 * @return Cells_t*
 ******************************************************************************/
Cells_t *allocate_cells(Mesh_t *mesh, int n_cells,
                        int max_cell_vertices, int max_cell_faces);

/*******************************************************************************
 * @brief Allocate the faces
 * @param mesh
 * @param n_faces
 * @param max_face_vertices
 * @return Faces_t*
 ******************************************************************************/
Faces_t *allocate_faces(Mesh_t *mesh, int n_faces, int max_face_vertices);

/*******************************************************************************
 * @brief Allocate the mesh
 * @return Mesh_t*
 ******************************************************************************/
Mesh_t *allocate_mesh();

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_partitions
 * @param n_partition_cells
 * @param n_partition_boundaries
 * @param n_partition_faces
 * @param n_partition_sends
 * @param n_partition_receives
 * @return Partition_t*
 ******************************************************************************/
Partition_t *allocate_partition(Mesh_t *mesh,
                                int n_partitions,
                                int n_partition_cells,
                                int n_partition_boundaries,
                                int n_partition_faces,
                                int n_partition_sends,
                                int n_partition_receives);

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_regions
 * @param max_name_length
 * @return Regions_t*
 ******************************************************************************/
Regions_t *allocate_regions(Mesh_t *mesh, int n_regions, int max_name_length);

/*******************************************************************************
 * @brief Allocate the mesh
 * @param mesh
 * @param n_vertices
 * @return Vertices_t*
 ******************************************************************************/
Vertices_t *allocate_vertices(Mesh_t *mesh, int n_vertices);

/*******************************************************************************
 * @brief Calculate the mesh metrics
 * @param mesh
 ******************************************************************************/
void calc_mesh_metrics(Mesh_t *mesh);

/*******************************************************************************
 * @brief Deallocate the boundaries
 * @param boundaries
 ******************************************************************************/
void deallocate_boundaries(Boundaries_t *boundaries);

/*******************************************************************************
 * @brief Deallocate the cells
 * @param cells
 ******************************************************************************/
void deallocate_cells(Cells_t *cells);

/*******************************************************************************
 * @brief Deallocate the faces
 * @param faces
 ******************************************************************************/
void deallocate_faces(Faces_t *faces);

/*******************************************************************************
 * @brief Deallocate the mesh
 * @param mesh
 ******************************************************************************/
void deallocate_mesh(Mesh_t *mesh);

/*******************************************************************************
 * @brief Deallocate the partition
 * @param partition
 ******************************************************************************/
void deallocate_partition(Partition_t *partition);

/*******************************************************************************
 * @brief Deallocate the regions
 * @param regions
 ******************************************************************************/
void deallocate_regions(Regions_t *regions);

/*******************************************************************************
 * @brief Deallocate the vertices
 * @param vertices
 ******************************************************************************/
void deallocate_vertices(Vertices_t *vertices);

/*******************************************************************************
 * @brief Define mesh
 ******************************************************************************/
void mesh_define();

/*******************************************************************************
 * @brief Finalize mesh
 ******************************************************************************/
void mesh_finalize();

/*******************************************************************************
 * @brief Initialize mesh
 ******************************************************************************/
void mesh_initialize();

/*******************************************************************************
 * @brief Print the boundaries
 * @param boundaries
 ******************************************************************************/
void print_boundaries(Boundaries_t *boundaries);

/*******************************************************************************
 * @brief Print the cellls
 * @param celss
 ******************************************************************************/
void print_cells(Cells_t *cells);

/*******************************************************************************
 * @brief Print the faces
 * @param faces
 ******************************************************************************/
void print_faces(Faces_t *faces);

/*******************************************************************************
 * @brief Print the mesh
 * @param mesh
 ******************************************************************************/
void print_mesh(Mesh_t *mesh);

/*******************************************************************************
 * @brief Print the partition
 * @param partition
 ******************************************************************************/
void print_partition(Partition_t *partition);

/*******************************************************************************
 * @brief Print the regions
 * @param regions
 ******************************************************************************/
void print_regions(Regions_t *regions);

/*******************************************************************************
 * @brief Print the vertices
 * @param vertices
 ******************************************************************************/
void print_vertices(Vertices_t *vertices);

/*******************************************************************************
 * @brief Read the mesh file
 * @param mesh_file
 * @return Mesh_t*
 ******************************************************************************/
Mesh_t *read_mesh_file(cstring_t mesh_file);

/*******************************************************************************
 * @brief Remap the mesh (global/local)
 * @param mesh
 ******************************************************************************/
void remap_local_mesh(Mesh_t *mesh);

#endif /* MESH_PRIVATE_H */