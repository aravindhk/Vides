// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
/* ************************************************************** */
/*      This subroutine compute the transmission coefficient      */
/*      given the Green matrix, the gamma functions               */
/*      and the order of the matrices (N)                         */
/* ************************************************************** */

#include "transmission.h"
double transmission(complex **G,complex **gamma1,complex **gamma2,int N)
{
  double **D1,T;
  complex zero;
  complex **C,**C1,**C2;
  int i,j,l;
  zero.r=0;
  zero.i=0;
  
  
/*    //I compute C=gamma1*G' where G' is the G conjugate/transpose */
  
/*    C=cmatrix(0,N-1,0,N-1); */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        { */
/*  	C[i][j]=zero; */
/*  	for (l=0;l<N;l++) */
/*  	  { */
/*  	    C[i][j].r=C[i][j].r+gamma1[i][l]*G[j][l].r; */
/*  	    C[i][j].i=C[i][j].i-gamma1[i][l]*G[j][l].i; */
/*  	  } */
/*        } */

/*    //I compute the matrix G*Gamma1*G' */
/*    C1=cmatmul(G,C,N); */
  
/*    //I compute the matrix Gamma2*G*Gamma1*G' */
/*    C2=cmatrix(0,N-1,0,N-1); */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        { */
/*  	C2[i][j]=zero; */
/*  	for (l=0;l<N;l++) */
/*  	  { */
/*  	    C2[i][j].r=C2[i][j].r+gamma2[i][l]*C1[l][j].r; */
/*  	    C2[i][j].i=C2[i][j].i+gamma2[i][l]*C1[l][j].i; */
/*  	  } */
/*        } */

  //I compute G'
  C=cmatdaga(G,N);
  //I compute the matrix G*Gamma1*G'
  C1=cmatmul3(G,gamma1,C,N);
  //I compute gamma2*G*gamma1*G'
  C2=cmatmul(gamma2,C1,N);

  //I take the real part
  D1=cmatRe(C2,N);

  //I compute Tr{Gamma2*G*Gamma1*G'}
  T=0;
  for (i=0;i<N;i++)
    T+=D1[i][i];
  cfree_cmatrix(C,0,N-1,0,N-1);
  cfree_cmatrix(C1,0,N-1,0,N-1);
  cfree_cmatrix(C2,0,N-1,0,N-1);
  free_dmatrix(D1,0,N-1,0,N-1);
  return T;
}
