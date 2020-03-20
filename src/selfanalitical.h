// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "complass.h"
#include "create_beta1.h"
#include "cmatrix.h"
#include "constants.in"
#include "create_beta2mode.h"
#include "create_beta2transpmode.h"
#include "cfree_cmatrix.h"
#include "VAVdaga.h"
#include "ccmul.h"
#include "ccadd.h"
#include "ccsub.h"
#include "cdmul.h"
#include "ccsqrt.h"
#include "ccdiv.h"
#include "cIm.h"
#include "cRe.h"
#ifndef SELFANALITICAL_H
#define SELFANALITICAL_H
complex **selfanalitical(double E,double Em1,int N,
			 int Nm,int *order,double thop,double eta);
#endif

