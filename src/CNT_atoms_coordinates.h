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
#include "nanotube.h"
#include "nrutil.h"

#if PY_MAJOR_VERSION >=3
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#endif
