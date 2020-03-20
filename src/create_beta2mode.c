// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta2mode.h"
complex **create_beta2mode(int n,int Nm,int *order,double thop)
{
  complex zero,**out;
  int i,j;
  double t;
  t=-thop;
  zero=complass(0,0);
  out=cmatrix(0,Nm-1,0,Nm-1);
  for (i=0;i<Nm;i++)
    {
      for (j=0;j<Nm;j++)
	{
	  if (i==j)
	    out[i][j]=complass(2*t*cos(-pi*(order[i])/n)*cos(pi*(order[i])/n),
			       -2*t*sin(pi*(order[i])/n)*cos(pi*(order[i])/n));
	  else
	    out[i][j]=zero;
	}
    }
  return out;
}

