// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatRe.h"
double **cmatRe(complex **A,int N)
{
  double **Re;
  int i,j;
  Re=dmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      Re[i][j]=A[i][j].r;
  return Re;
}
  
