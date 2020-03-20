// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "complass.h"
#include "cmatrix.h"
#include "cmatnorm2diff.h"
#include "cmatinv.h"
#include "cmatmul.h"
#include "cmatmul3.h"
#include "cmatsum.h"
#include "cmatsub.h"
#include "constants.in"
#include "cfree_cmatrix.h"
#ifndef GZEROZERO_H
#define GZEROZERO_H
complex **Gzerozero(complex **wmH,complex **BETA,complex **BETADAGA,int N);
     //complex **Gzerozero(double E,complex **H,complex **BETA,complex **BETADAGA,int N)
#endif

