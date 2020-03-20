// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
// This subrountine perform a basis transformation, using the
// eigenvactors of the BETA2 matrix.
// A is the matrix to be transformed, N is the order of the
// matrix and Nm is the number of the considered eigenstates.
#include "VdagaAV.h"
complex **VdagaAV(complex **A,int N,int Nm,int *order)
{
  int i,j,l;
  complex **D,**M,zero,temp, **V,**Vdaga;
  zero.r=0;
  zero.i=0;
  V=cmatrix(0,N-1,0,Nm-1);
  Vdaga=cmatrix(0,Nm-1,0,N-1);
  // I create the matrix V
  for (i=0;i<N;i++)
    for (j=0;j<Nm;j++)
      {
	V[i][j]=complass(1/sqrt(N)*cos(2*pi*(order[j])/N*i),
			 1/sqrt(N)*sin(2*pi*(order[j])/N*i));
	Vdaga[j][i]=complass(V[i][j].r,-V[i][j].i);
      }

  //  // I create the matrix V
  //  for (i=0;i<N;i++)
  //    for (j=0;j<Nm;j++)
  //      {
  //	V[i][j]=complass(1/sqrt(N)*cos(2*pi*(j+1)/N*i),
  //			 1/sqrt(N)*sin(2*pi*(j+1)/N*i));
  //	Vdaga[j][i]=complass(V[i][j].r,-V[i][j].i);
  //      }



/*   for (i=0;i<N;i++) */
/*     { */
/*       for (j=0;j<Nm;j++) */
/* 	printf("%lg ",V[i][j].r); */
/*       printf("\n"); */
/*     } */

/*   printf("********************\n"); */
  
/*   for (i=0;i<N;i++) */
/*     { */
/*       for (j=0;j<Nm;j++) */
/* 	printf("%lg ",V[i][j].i); */
/*       printf("\n"); */
/*     } */
  
/*       exit(1); */
  

  //I compute M=A*V
  M=cmatrix(0,N-1,0,Nm-1);
  for (i=0;i<N;i++)
    for (j=0;j<Nm;j++)
      {
	M[i][j]=zero;
	for (l=0;l<N;l++)
	  {
	    M[i][j].r=M[i][j].r+A[i][l].r*V[l][j].r-A[i][l].i*V[l][j].i;
	    M[i][j].i=M[i][j].i+A[i][l].i*V[l][j].r+A[i][l].r*V[l][j].i;
	  }
      }

  //I compute D=Vdaga*M
  D=cmatrix(0,Nm-1,0,Nm-1);
  for (i=0;i<Nm;i++)
    for (j=0;j<Nm;j++)
      {
	D[i][j]=zero;
	for (l=0;l<N;l++)
	  {
	    D[i][j].r=D[i][j].r+Vdaga[i][l].r*M[l][j].r-Vdaga[i][l].i*M[l][j].i;
	    D[i][j].i=D[i][j].i+Vdaga[i][l].i*M[l][j].r+Vdaga[i][l].r*M[l][j].i;
	  }
      }
  cfree_cmatrix(M,0,N-1,0,Nm-1);
  cfree_cmatrix(V,0,N-1,0,Nm-1);
  cfree_cmatrix(Vdaga,0,Nm-1,0,N-1);
  return D;
}
	
