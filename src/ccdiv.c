// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "ccdiv.h"
complex ccdiv(complex a,complex b)
{
  complex out;
  double ar,br,ai,bi,r;
  ar=cRe(a);
  ai=cIm(a);
  br=cRe(b);
  bi=cIm(b);
  out=complass((ar*br+ai*bi),(ai*br-ar*bi));
  r=br*br+bi*bi;
  r=1/r;
  if (fabs(r)>1e-40)
    {
      out=cdmul(r,out);
      return out;
    }
  else
    {
      printf("Error in ccdiv called by selfanalitical \n");
      printf("division by zero \n");
      exit(1);
    }
}
