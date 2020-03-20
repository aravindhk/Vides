// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatsum.h"
complex **cmatsum(complex **A,complex **B,int N)
{
  int i,j;
  complex **C;
  C=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	C[i][j].r=A[i][j].r+B[i][j].r;
	C[i][j].i=A[i][j].i+B[i][j].i;
      }
  return C;
}
	
