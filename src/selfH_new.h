// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "complex.h"
#include "complass.h"
#include "cmatrix.h"
#include "cmatinv.h"
#include "cmatsub.h"
#include "cmatsum.h"
#include "cdabs.h"
#include "constants.in"
#include "cfree_cmatrix.h"
#include "cfree_ccvector.h"
#include "GzerozeroH_W.h"
#include "Gzerozero.h"
#include "eigenvalues_non_symmetric_matrix.h"
#include "cIm.h"
#include "cRe.h"
#include "cmul.h"
#include "csum.h"
#include "csub.h"
#include "cmatvect.h"
#include "domn.h"
#include "msolve.h"
#include "matvec.h"
#include "flip_cmatrix.h"
#ifndef SELFHNEW_H
#define SELFHNEW_H
complex **selfH_new(double E, complex ***diag, complex ***updiag, complex ***lowdiag, int N, int Nc, int lead, double eta);
#endif

