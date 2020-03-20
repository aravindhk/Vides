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
#include "cfree_cmatrix.h"
#ifndef NR_END
#define NR_END 1
#endif
#ifndef FREE_ARG
#define FREE_ARG char*
#endif
#ifndef CFREE_CTENSOR4_H
#define CFREE_CTENSOR4_H
void cfree_ctensor4(complex ****m,int n1l,int n1h,int n2l,int n2h,int n3l,int n3h, int n4l, int n4h,int flagball);
#endif
