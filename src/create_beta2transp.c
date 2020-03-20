// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_beta2transp.h"
complex **create_beta2transp(int n,double thop)
{
  complex zero,t,**out,**out2;
  int i,j;
  
  out2=cmatrix(0,n-1,0,n-1);
  out=create_beta2(n,thop);
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      out2[i][j]=out[j][i];
  cfree_cmatrix(out,0,n-1,0,n-1);
  return out2;
}

