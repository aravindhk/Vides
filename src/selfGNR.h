// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "complex.h"
#include "complass.h"
#include "create_beta1GNR.h"
#include "cmatrix.h"
#include "cmatinv.h"
#include "constants.in"
#include "create_beta2GNR.h"
#include "create_beta2transpGNR.h"
#include "cfree_cmatrix.h"
#include "create_updiagGNR.h"
#include "create_lowdiagGNR.h"
#include "Gzerozero.h"
#ifndef SELFgnr_H
#define SELFgnr_H
complex **selfGNR(double E,double *Em1,int N,double thop,double eta);
#endif

