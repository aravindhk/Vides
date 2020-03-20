// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta2.h"
complex **create_beta2(int n,double thop)
{
  complex zero,t,**out;
  int i,j;
  
  t=complass(-thop,0);
  zero=complass(0,0);
  out=cmatrix(0,n-1,0,n-1);
  for (i=1;i<n;i++)
    {
      for (j=0;j<n;j++)
	{
	  if ((i==j)||(j==(i-1)))
	    out[i][j]=t;
	  else
	    out[i][j]=zero;
	}
    }
  for (j=0;j<n;j++)
    {
      if ((j==0)||(j==(n-1)))
	out[0][j]=t;
      else
	out[0][j]=zero;
    }
  return out;
}

