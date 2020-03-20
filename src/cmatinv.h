// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "cmatrix.h"
#include "zgetri.h"
#include "ccvector.h"
#include "nrutil.h"
#include "zgetrf.h"
#include "cfree_ccvector.h"
#ifndef CMATINV_H
#define CMATINV_H
complex **cmatinv(complex **A,int N);
#endif
