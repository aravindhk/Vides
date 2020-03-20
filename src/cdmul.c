// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cdmul.h"
complex cdmul(double a,complex b)
{
  complex c;
  c.r=a*b.r;
  c.i=a*b.i;
  return c;
}
