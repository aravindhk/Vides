// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatnorm2diff.h"
double cmatnorm2diff(complex **A,complex **B,int N)
{
  double norma,num,den;
  complex **temp,**S,**D;
  int i,j;
  norma=0;
  num=0;
  den=0;
  temp=cmatsub(A,B,N);
  S=cmatsum(A,B,N);
  D=cmatsub(A,B,N);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	den+=cdabs(S[i][j]);
	num+=cdabs(D[i][j]);
	norma+=sqrt(temp[i][j].r*temp[i][j].r+temp[i][j].i*temp[i][j].i);
      }
  cfree_cmatrix(temp,0,N-1,0,N-1);
  cfree_cmatrix(D,0,N-1,0,N-1);
  cfree_cmatrix(S,0,N-1,0,N-1);
  //return norma;
  return num/den;
}
