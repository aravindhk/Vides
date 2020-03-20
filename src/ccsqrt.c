// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "ccsqrt.h"
complex ccsqrt(complex z)
{
  complex out;
  double r,xtmp,x,y,ytmp;
  x=cRe(z);
  y=cIm(z);
  r=cdabs(z);
  if (fabs(r)<1e-40)
    {
      out=complass(0,0);
      return out;
    }
  xtmp=sqrt((r+fabs(x))*0.5);
  ytmp=y*0.5/xtmp;
  if (x>=0)
    out=complass(xtmp,ytmp);
  if (fabs(y)<1e-40)
    y=1;
  if (x<0) 
    {
      if (y<0)
	out=complass(fabs(ytmp),-1*fabs(xtmp));
      else
	out=complass(fabs(ytmp),fabs(xtmp));
    }
  return out;
}
	  
		     

