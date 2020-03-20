// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta1.h"
complex **create_beta1(int n,double thop)
{
  complex zero,t,**out;
  int i,j;
  
  t=complass(-thop,0);
  zero=complass(0,0);
  out=cmatrix(0,n-1,0,n-1);
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      {
	if (i==j)
	  out[i][j]=t;
	else
	  out[i][j]=zero;
      }
  return out;
}
