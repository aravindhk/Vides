// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#ifndef PHYSICALQUANTITIES_H
#define PHYSICALQUANTITIES_H
typedef struct {
  double *Phi;
  double *Na;
  double *Nd;
  double *Ef;
  double *rhof;
  double *vt;
  double *eps;
  double *free_charge;
  double *fixed_charge;
  double *Phiold;
} physical_quantities;
#endif
