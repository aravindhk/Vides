// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include <math.h>
#include "abbs.h"
#include "indice.h"
//#include "soliduzzo.h"
#include "Jbuild.h"
#include "physical_quantities.h"
#include "device_mapping.h"
//#include "dfunpoint.h"
#ifndef PREPAREJ_H
#define PREPAREJ_H
void preparej(physical_quantities p,device_mapping d,
	      double *a2,int *ia2,int *ja2,int *z12,int flagum,
	      double *****J,int neq);
#endif
