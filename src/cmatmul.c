// ======================================================================
//  Copyright (c) 2004-2008, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatmul.h"
#include "zgemm.h"
#include <string.h>
#include "ccvector.h"
#include "cfree_ccvector.h"
complex **cmatmul(complex **A,complex **B,int N)
{
  int i,j,l,ix;
  complex **C,zero,temp,uno,*AA,*BB,*CC;
  char TRANSA[10],TRANSB[10];

  zero.r=0;
  zero.i=0;
  uno.r=1;
  uno.i=0;

  strcpy(TRANSA,"N");
  strcpy(TRANSB,"N");

  //  printf("N=%d \n",N);

  /*  for (i=0;i<N;i++)
    {
      for (j=0;j<N;j++)
	printf("%lg ",A[i][j].r);
      printf("\n");
    }

  for (i=0;i<N;i++)
    {
      for (j=0;j<N;j++)
	printf("%lg ",B[i][j].r);
      printf("\n");
    }

    printf("-------------------------------\n");*/

  //  AA=ccvector(0,N*N-1);
  //BB=ccvector(0,N*N-1);
  //  CC=ccvector(0,N*N-1);

  AA=&B[0][0];
  BB=&A[0][0];
  C=cmatrix(0,N-1,0,N-1);

  /*  for (i=0;i<N*N;i++)
    printf("%lg ",AA[i].r);

  for (i=0;i<N*N;i++)
  printf("%lg ",BB[i].r);*/
  //  exit(0);

  /*  for (i=0;i<N;i++){
      for (j=0;j<N;j++){
	ix=i+j*N;
	//	AA[ix].r=A[i][j].r;
	//	AA[ix].i=A[i][j].i;
	BB[ix].r=B[i][j].r;
        BB[ix].i=B[i][j].i;
      }
      }*/
  
  zgemm_(TRANSA,TRANSB,&N,&N,&N,&uno,AA,&N,BB,&N,&zero,&C[0][0],&N);
  
  /*  C=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      ix=i+j*N;
      C[i][j]=CC[ix];
      //      printf("%lg ",C[i][j].r);
    }
  }
  //  printf("-------------------------------\n");
  //  cfree_ccvector(AA,0,N*N-1);
  //cfree_ccvector(BB,0,N*N-1);
  cfree_ccvector(CC,0,N*N-1);
  */

  return C;
  

  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	C[i][j]=zero;
	for (l=0;l<N;l++)
	  if ((fabs(B[l][j].r)>0)||(fabs(B[l][j].i)>0))
	    {
	      C[i][j].r=C[i][j].r+A[i][l].r*B[l][j].r-A[i][l].i*B[l][j].i;
	      C[i][j].i=C[i][j].i+A[i][l].i*B[l][j].r+A[i][l].r*B[l][j].i;
	    }
      }
  return C;
}
	
