// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "ccsub.h"
complex ccsub(complex a,complex b)
{
  complex out;
  double ar,br,ai,bi;
  ar=cRe(a);
  ai=cIm(a);
  br=cRe(b);
  bi=cIm(b);
  out=complass((ar-br),(ai-bi));
  return out;
}
