// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
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
#include "preparej2D.h"
#include "ffunpoisson2D.h"
#include "norma2.h"
#include "msolve.h"
#include "matvec.h"
#include "domn.h"
#include "nrutil.h"
#include "create_J2D.h"
#include "preparej2D.h"
#include "destroy_J2D.h"
#include "constants.in"
//#include "stampaout.h"
#include "physical_quantities.h"
#include "device_mapping.h"

#if PY_MAJOR_VERSION >=3
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#endif

//#include "dfunpoint.h"
#ifndef SOLVEpOISSON2d_H
#define SOLVEpOISSON2d_H
static PyObject* py_solvePoisson2D(PyObject* self, PyObject* args);
//void solvesystem(physical_quantities p,device_mapping d,
//		    int *IERR,double tolldomn,
//		    int neq,soliduzzo *mat,
//		    double normapoisson,double sottoril,
//		    dfunpoint *nf,dfunpoint *pf,dfunpoint *accf,dfunpoint *donf,
//		    dfunpoint *dnf,dfunpoint *dpf,dfunpoint *daccf,dfunpoint *ddonf,int rank);
#endif
