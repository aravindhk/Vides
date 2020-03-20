// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "cmatrix.h"
#include "cvectorm.h"
#include "ctensor4.h"
#include "cmatinv.h"
#include "cmatmul.h"
#include "cmatsub.h"
#include "cmatsum.h"
#include "cfree_cmatrix.h"
#include "cfree_cvectorm.h"
#include "cfree_ctensor4.h"
#include "cmatmul3.h"
#include "cdabs.h"
#ifndef RGFBLOCK_H
#define RGFBLOCK_H
complex ****rgfblock(complex ***lowdiag,complex ***diag,
		   complex ***updiag,int n,int Nc,int flagball);
#endif
