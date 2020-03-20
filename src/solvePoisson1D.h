// ======================================================================
//  Copyright (c) 2004-2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Python.h"
#include "arrayobject.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "soliduzzo.h"
#include "ffunpoisson1D.h"
#include "norma2.h"
#include "msolve.h"
#include "matvec.h"
#include "domn.h"
#include "nrutil.h"
#include "preparej1D.h"
#include "constants.in"
//#include "stampaout.h"
#include "physical_quantities.h"
#include "device_mapping.h"
//#include "dfunpoint.h"
#ifndef SOLVEpOISSON1d_H
#define SOLVEpOISSON1d_H
static PyObject* py_solvePoisson1D(PyObject* self, PyObject* args);
//void solvesystem(physical_quantities p,device_mapping d,
//		    int *IERR,double tolldomn,
//		    int neq,soliduzzo *mat,
//		    double normapoisson,double sottoril,
//		    dfunpoint *nf,dfunpoint *pf,dfunpoint *accf,dfunpoint *donf,
//		    dfunpoint *dnf,dfunpoint *dpf,dfunpoint *daccf,dfunpoint *ddonf,int rank);
#endif
