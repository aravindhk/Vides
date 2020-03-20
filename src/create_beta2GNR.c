// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta2GNR.h"
complex **create_beta2GNR(int n,double thop)
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
      //if ((j==0)||(j==(n-1)))
      if (j==0)
	out[0][j]=t;
      else
	out[0][j]=zero;
    }
  return out;
}

