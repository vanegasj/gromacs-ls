// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.h
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_H
#define VOROPP_COMMON_H

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "config.h"

namespace voro {

void check_duplicate(int n,double x,double y,double z,int id,double *qp);

void voro_fatal_error(const char *p,int status);
void voro_print_positions(std::vector<double> &v,FILE *fp=stdout);
FILE* safe_fopen(const char *filename,const char *mode);
void voro_print_vector(std::vector<int> &v,FILE *fp=stdout);
void voro_print_vector(std::vector<double> &v,FILE *fp=stdout);
void voro_print_face_vertices(std::vector<int> &v,FILE *fp=stdout);

}

#endif
