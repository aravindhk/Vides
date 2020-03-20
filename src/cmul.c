// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmul.h"
complex cmul(complex a,complex b)
{
  complex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.r*b.i+a.i*b.r;
  return c;
}
