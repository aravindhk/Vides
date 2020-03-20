// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta2transpGNR.h"
complex **create_beta2transpGNR(int n,double thop)
{
  complex zero,t,**out,**out2;
  int i,j;
  
  out2=cmatrix(0,n-1,0,n-1);
  out=create_beta2GNR(n,thop);
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      out2[i][j]=out[j][i];
  cfree_cmatrix(out,0,n-1,0,n-1);
  return out2;
}

