// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatdaga.h"
complex **cmatdaga(complex **A,int N)
{
  int i,j,l;
  complex **C,zero,temp;
  C=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      C[i][j]=complass(A[j][i].r,-A[j][i].i);
  return C;
}
	
