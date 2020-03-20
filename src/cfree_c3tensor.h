// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#ifndef NR_END
#define NR_END 1
#endif
#ifndef FREE_ARG
#define FREE_ARG char*
#endif
#ifndef CFREE_C3TENSOR_H
#define CFREE_C3TENSOR_H
void cfree_c3tensor(complex ***t, int nrl, int nrh, int ncl, int nch,
		    int ndl, int ndh);
#endif
