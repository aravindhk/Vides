// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cfree_ccvector.h"
void cfree_ccvector(complex *m, int nrl, int nrh)
{
  free((FREE_ARG) (m+nrl-NR_END));
}
