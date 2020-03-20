// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "constants.in"
#include "cmatrix.h"
#include "csum.h"
#include "complass.h"
#include "c3tensor.h"
#include "cfree_c3tensor.h"
#include "cmatsub.h"
#include "rgfblock_Lake.h"
#include "rgfblock.h"
#include "cIm.h"
#include "nrutil.h"
#include "spectralfun.h"
#include "cfree_ctensor4.h"
#include "cmatdaga.h"
#include "csub.h"
#include "transmission.h"
#ifndef LDOS_lAKE_H
#define LDOS_lAKW_H
void LDOS_Lake(double E,complex ***lowdiag,complex ***diag,complex ***updiag,
	       double ***A1,double ***A2,complex **sigmas,complex **sigmad,int n,
	       int Nc,int flagtrans,double *T,double thop,double eta);
#endif
