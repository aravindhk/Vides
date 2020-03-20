// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "abbs.h"
double abbs(double x1)
{
  double ora;
  if (x1>0)
    ora=x1;
  else
    ora=-x1;
  return ora;
}
