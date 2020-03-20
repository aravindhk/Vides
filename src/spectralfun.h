// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "cmatrix.h"
#include "nrutil.h"
#include "constants.in"
#include "cfree_cmatrix.h"
#include "cdabs.h"
#include "complass.h"
#include "cmatmul.h"
#ifndef SPECTRALFUN_H
#define SPECTRALFUN_H
double *spectralfun(complex **G,complex **gamma1,int N);
#endif
