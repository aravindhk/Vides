// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include <stdio.h>
#include <math.h>
#include "abbs.h"
#include "indice2D.h"
//#include "soliduzzo.h"
#include "Jbuild2D.h"
#include "physical_quantities.h"
#include "device_mapping.h"
//#include "dfunpoint.h"
#ifndef PREPAREJ2d_H
#define PREPAREJ2d_H
void preparej2D(physical_quantities p,device_mapping d,
		double *a2,int *ia2,int *ja2,int *z12,int flagum,
		double ****J,int neq);
#endif
