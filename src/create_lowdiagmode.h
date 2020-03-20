// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "complex.h"
#include "create_beta1.h"
#include "create_beta2mode.h"
#include "create_beta2transpmode.h"
#ifndef CREATE_LOWDIAGMODE_H
#define CREATE_LOWDIAGMODE_H
complex **create_lowdiagmode(int Nc,int n,int Nm,int *order,double thop);
#endif
