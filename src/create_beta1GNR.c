// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta1GNR.h"
complex **create_beta1GNR(int n,int Nc,double thop)
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
  if (((Nc%4)>=0)&&((Nc%4)<=1))
    out[0][0]=complass(-thop*1.12,0);
  else
    out[n-1][n-1]=complass(-thop*1.12,0);
  return out;
}
