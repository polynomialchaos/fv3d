//##################################################################################################################################
// FV3D - Finite volume solver
// (c) 2020 | Florian Eigentler
//##################################################################################################################################
#ifndef IMPLICIT_MODULE_H
#define IMPLICIT_MODULE_H

#include "fv3d/fv3d_module.h"



extern int implicit_active;
extern int n_iter_inner;
extern int n_iter_lsoe;
extern int n_bdf_stages;
extern double **phi_old;

void implicit_define();

#endif /* IMPLICIT_MODULE_H */