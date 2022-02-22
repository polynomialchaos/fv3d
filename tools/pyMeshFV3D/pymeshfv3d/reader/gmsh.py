################################################################################
# @file gmsh.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import os
import logging
from pymeshfv3d.utilities import ElementType, Element, Mesh

gmsh_to_ElementType = [None] * 15
gmsh_to_ElementType[14] = ElementType.POINT
gmsh_to_ElementType[1] = ElementType.LINE
gmsh_to_ElementType[2] = ElementType.TRIANGLE
gmsh_to_ElementType[3] = ElementType.QUADRANGLE
gmsh_to_ElementType[4] = ElementType.TETRAEDER
gmsh_to_ElementType[5] = ElementType.HEXAEDER
gmsh_to_ElementType[6] = ElementType.PRISM
gmsh_to_ElementType[7] = ElementType.PYRAMID

gmsh_elem_nodes = [None] * 15
gmsh_elem_nodes[14] = 1
gmsh_elem_nodes[1] = 2
gmsh_elem_nodes[2] = 3
gmsh_elem_nodes[3] = 4
gmsh_elem_nodes[4] = 4
gmsh_elem_nodes[5] = 8
gmsh_elem_nodes[6] = 6
gmsh_elem_nodes[7] = 5


def read_handler_gmsh(file_name):
    """Generate a mesh (1D) by reading Gmsh file."""
    logging.info('Generate mesh (Gmsh) p="%s"', file_name)

    with open(file_name, 'r') as fp:
        while True:
            line = fp.readline().strip()

            if line == '$MeshFormat':
                check_read_line(fp.readline(), '4.1 0 8')
                check_read_line(fp.readline(), '$EndMeshFormat')
            elif line == '$PhysicalNames':
                n_regions = int(fp.readline())

                regions = []
                for i in range(n_regions):
                    _, id_r, name = fp.readline().split(' ')
                    if i != int(id_r) - 1:
                        raise IndexError('Region', i, id_r)
                    regions.append(name.replace('"', '').strip())

                check_read_line(fp.readline(), '$EndPhysicalNames')
            elif line == '$Entities':
                n_entities = [int(x) for x in fp.readline().split(' ')]
                if len(n_entities) != 4:
                    raise ValueError('Number of Entities', n_entities)

                entities = []
                for i in range(4):
                    entities.append([])

                    for j in range(n_entities[i]):
                        tmp = fp.readline().split(' ')

                        if i == 0:
                            id_e, n_tags_1 = int(tmp[0]), int(tmp[4])
                            tags_1 = [int(x) for x in tmp[5:5+n_tags_1]]
                            # tags_2 = None
                        else:
                            id_e, n_tags_1 = int(tmp[0]), int(tmp[7])
                            tags_1 = [int(x) for x in tmp[8:8+n_tags_1]]
                            # n_tags_2 = int( tmp[8+n_tags_1] )
                            # tags_2 = [int( x ) for x in
                            # tmp[9+n_tags_1:9+n_tags_1+n_tags_2]]

                        r_id = tags_1[0]-1 if tags_1 else None
                        entities[i].append((id_e, r_id))

                check_read_line(fp.readline(), '$EndEntities')
            elif line == '$Nodes':
                n_blocks, n_vertices, _, _ = [
                    int(x) for x in fp.readline().split(' ')]
                vertices = [None] * (n_vertices)

                for i in range(n_blocks):
                    _, _, _, n_sub_blocks = \
                        [int(x) for x in fp.readline().split(' ')]

                    vertice_ids = []
                    for j in range(n_sub_blocks):
                        tmp = fp.readline()
                        vertice_ids.append(int(tmp))

                    for j in range(n_sub_blocks):
                        vertices[vertice_ids[j] - 1] = \
                            [float(x) for x in fp.readline().split(' ')]

                check_read_line(fp.readline(), '$EndNodes')
            elif line == '$Elements':
                n_blocks, n_elements, _, _ = [
                    int(x) for x in fp.readline().split(' ')]
                elements = [None] * (n_elements)

                for i in range(n_blocks):
                    id_d, id_e, e_type, n_sub_blocks = [
                        int(x) for x in fp.readline().split(' ')]

                    for j in range(n_sub_blocks):
                        tmp = [int(x)
                               for x in fp.readline().strip().split(' ')]
                        e_id, e_nodes = tmp[0], [
                            x-1 for x in tmp[1:1+gmsh_elem_nodes[e_type]]]
                        entitiy = get_entity_by_id(entities, id_d, id_e)
                        elements[e_id-1] = Element(
                            gmsh_to_ElementType[e_type], e_nodes, entitiy[1])

                check_read_line(fp.readline(), '$EndElements')
            else:
                if not line:
                    break
                raise KeyError('Identifier', line)

    return Mesh('Gmsh({:})'.format(os.path.basename(file_name)),
                elements, regions, vertices)


def get_entity_by_id(entities, id_d, id_e):
    for x in entities[id_d]:
        if x[0] == id_e:
            return x


def check_read_line(line, check):
    if line.strip() != check:
        raise Exception('Check', line, check)
