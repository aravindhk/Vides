// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include "complex.h"
#include "constants.in"
#include "cvectorm.h"
#include "cmatrix.h"
#include "csum.h"
#include "cdmul.h"
#include "complass.h"
#include "cfree_cvectorm.h"
#include "cmatsub.h"
#include "rgfblock.h"
#include "cIm.h"
#include "nrutil.h"
#include "spectralfunmode.h"
#include "cfree_ctensor4.h"
#include "cmatdaga.h"
#include "transmission.h"
#ifndef LDOSMODE_H
#define LDOSMODE_H
void LDOSmode(double E,complex ***lowdiag,complex ***diag,complex ***updiag,
	      double ***A1,double ***A2,complex **sigmas,complex **sigmad,int n,
	      int Nc,int flagtrans,double *T,int Nreal,int *order,double thop,
	      double eta);
#endif
