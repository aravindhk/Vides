// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatinv.h"
complex **cmatinv(complex **A,int N)
{
  int i,j,k,*IPIV,INFO;
  complex **C,zero,*temp,*WORK;
  zero.r=0;
  zero.i=0;
  IPIV=ivector(0,N-1);
  WORK=ccvector(0,2*N-1);
  C=cmatrix(0,N-1,0,N-1);
  temp=ccvector(0,N*N-1);
  k=0;
  for (j=0;j<N;j++)
    for (i=0;i<N;i++)
      {
	temp[k]=A[i][j];
	k++;
      }
  zgetrf_(&N,&N,temp,&N,IPIV,&INFO);
  zgetri_(&N,temp,&N,IPIV,WORK,&N,&INFO);
  k=0;
  for (j=0;j<N;j++)
    for (i=0;i<N;i++)
      {
	C[i][j]=temp[k];
	k++;
      }
  //Free memory
  free_ivector(IPIV,0,N-1);
  cfree_ccvector(WORK,0,2*N-1);
  cfree_ccvector(temp,0,N*N-1);
  return C;
}
