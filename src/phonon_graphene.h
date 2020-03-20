// ======================================================================
//  Copyright (c) 2009, A. Betti, G. Fiori, University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Python.h"
#include "arrayobject.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "constants.in"
#include "nrutil.h"
#include "complex.h"
#include "ceig.h"
#include "cmatrix.h"
#include "zhpevx.h"
#include "free_dvectorm.h"
#include "free_ivectorm.h"
