//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef OUTPUT_MODULE_H
#define OUTPUT_MODULE_H

#include "fv3d/fv3d_module.h"



extern int do_output_data;
extern int i_output_data;
extern string_t output_file;

void output_define();

void create_file_header();
void write_output(int iter, double t);

#endif /* OUTPUT_MODULE_H */