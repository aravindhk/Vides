// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#ifndef DEVICEMAPPING_H
#define DEVICEMAPPING_H
typedef struct {
  int nx;
  int ny;
  int nz;
  int n;
  int Nc;
  int *swap;
  double *boundary_conditions;
  double **dist;
  double **surf;
  double *surf2D;
  double *dVe;
} device_mapping;
#endif
