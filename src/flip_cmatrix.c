// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "flip_cmatrix.h"
void flip_cmatrix(complex **A,int N)
{
  int i,j;
  complex **change;
  change=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N-i;j++)
      {
	change[N-j-1][N-(j+i)-1]=complass(A[j][j+i].r,-A[j][j+i].i);
	change[N-(j+i)-1][N-j-1]=complass(A[j+i][j].r,-A[j+i][j].i);
	//	change[N-i][N-j]=complass(A[i][j].r,-A[i][j].i);
      }

  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      A[i][j]=change[i][j];
  
  cfree_cmatrix(change,0,N-1,0,N-1);
  return;
}
	
