// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "cmatsub.h"
#include "cdabs.h"
#include "cmatsum.h"
#include <math.h>
#include "cfree_cmatrix.h"
#ifndef CMATNORM2DIFF_H
#define CMATNORM2DIFF_H
double cmatnorm2diff(complex **A,complex **B,int N);
#endif
