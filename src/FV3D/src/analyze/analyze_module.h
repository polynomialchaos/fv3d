//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef ANALYZE_MODULE_H
#define ANALYZE_MODULE_H

#include "fv3d/fv3d_module.h"



extern double *residual;

void analyze_define();

void calc_global_residual(double dt);

#endif /* ANALYZE_MODULE_H */