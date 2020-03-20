// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cdabs.h"
double cdabs(complex a)
{
  double d;
  d=sqrt(pow(a.r,2)+pow(a.i,2));
  return d;
}
