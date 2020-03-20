// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include "sign.h"
#include <math.h>
#include "constants.in"
//#include "soliduzzo.h"
#include "physical_quantities.h"
#include "device_mapping.h"
#ifndef FFUNPOISSON1d_H
#define FFUNPOISSON1d_H
void ffunpoisson1D(physical_quantities p,device_mapping d,
		   double *bb);
#endif
