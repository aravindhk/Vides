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
#include "cmatmul.h"
#include "cmatRe.h"
#include "cmatmul3.h"
#include "cmatdaga.h"
#ifndef TRANSMISSION_H
#define TRANSMISSION_H
double transmission(complex **G,complex **gamma1,complex **gamma2,int N);
#endif
