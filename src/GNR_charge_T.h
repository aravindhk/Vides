// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Python.h"
#include "arrayobject.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "LDOS.h"
#include "nrutil.h"
#include "complex.h"
#include "complass.h"
#include "cmatrix.h"
#include "cvectorm.h"
#include "max.h"
#include "create_updiagGNR.h"
#include "create_lowdiagGNR.h"
#include "cfree_cmatrix.h"
#include "cfree_cvectorm.h"
#include "constants.in"
#include "min.h"
#include "selfschottky.h"
#include "selfGNR.h"
#include "Fermi_Dirac.h"
#include "gappo.h"

#if PY_MAJOR_VERSION >=3
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#endif
