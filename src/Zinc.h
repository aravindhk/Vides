// ======================================================================
//  Copyright (c) 2011, P. D'Amico, G. Fiori, G. Iannaccone  University of Pisa
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
#include "max.h"
#include "create_updiag.h"
#include "create_lowdiag.h"
#include "cfree_cmatrix.h"
#include "cfree_c3tensor.h"
#include "constants.in"
#include "min.h"
#include "selfH_W.h"
#include "Fermi_Dirac.h"
#include "Hamiltonian_py.h"

#if PY_MAJOR_VERSION >=3
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#endif
