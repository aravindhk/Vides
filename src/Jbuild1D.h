// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
//#include "soliduzzo.h"
#include <math.h>
#include "sign.h"
#include "constants.in"
#include "physical_quantities.h"
#include "device_mapping.h"
#ifndef jBUILD1d_H
#define jBUILD1d_H
void Jbuild1D(physical_quantities p,device_mapping d,
	      double ***J,int neq);
#endif
