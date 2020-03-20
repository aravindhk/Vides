// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
/* ************************************************************** */
/*      This subroutine compute the spectral function             */
/*      given the Green matrix, the gamma function                */
/*      and the order of the matrices (N),and stores the results  */
/*      in a double precision vector                              */
/* ************************************************************** */

#include "spectralfun.h"
double *spectralfun(complex **G,complex **gamma1,int N)
{
  double *SP,g1max,Gmax,**gamma1temp;
  complex zero,**Gtemp;
  complex **C;
  complex **Tr;
  int i,j,l;
  SP=dvector(0,N-1);
  zero.r=0;
  zero.i=0;


/*    gamma1temp=dmatrix(0,N-1,0,N-1); */
/*    Gtemp=cmatrix(0,N-1,0,N-1); */
/*    // I compute the maximum abs of gamma1 */
/*    g1max=-1e30; */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        if (g1max<=fabs(gamma1[i][j])) */
/*  	g1max=fabs(gamma1[i][j]); */
  
/*    // I compute the maximum abs of G */
/*    Gmax=-1e30; */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        if (cdabs(G[i][j])>=Gmax) */
/*  	Gmax=cdabs(G[i][j]); */
  
  //  printf("%lg %lg \n",g1max,Gmax);

/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        gamma1[i][j]=1; */
  

  
/*    g1max=1; */
/*    Gmax=1; */

/*    for (i=0;i<N;i++) */
/*      { */
/*        for (j=0;j<N;j++) */
/*  	printf("%lg ",gamma1[i][j]); */
/*        printf("\n"); */
/*      } */

/*    printf("------------------------------\n"); */

/*    // I normalize G and gamma1 */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        { */
/*  	Gtemp[i][j]=complass(G[i][j].r/Gmax,G[i][j].i/Gmax); */
/*  	gamma1temp[i][j]=gamma1[i][j]/g1max; */
/*        } */

/*    for (i=0;i<N;i++) */
/*      { */
/*        for (j=0;j<N;j++) */
/*  	printf("%lg ",gamma1[i][j]); */
/*        printf("\n"); */
/*      } */

/*    printf("--------------------------------\n"); */

/*    for (i=0;i<N;i++) */
/*      { */
/*        for (j=0;j<N;j++) */
/*  	printf("%lg ",G[i][j].r); */
/*        printf("\n"); */
/*      } */

/*    printf("--------------------------------\n"); */

/*    for (i=0;i<N;i++) */
/*      { */
/*        for (j=0;j<N;j++) */
/*  	printf("%lg ",G[i][j].i); */
/*        printf("\n"); */
/*      } */


/*    printf("******************************\n"); */
  
  
  
  Tr=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++) 
    for (j=0;j<N;j++)
      Tr[i][j]=complass(G[j][i].r,-G[j][i].i);

  C=cmatmul(gamma1,Tr,N);
  
/*    C=cmatrix(0,N-1,0,N-1); */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        { */
/*  	C[i][j]=zero; */
/*  	for (l=0;l<N;l++) */
/*  	  { */
/*  	    C[i][j].r=C[i][j].r+gamma1[i][l]*Tr[l][j].r; */
/*  	    C[i][j].i=C[i][j].i+gamma1[i][l]*Tr[l][j].i; */
/*  	  } */
/*        } */
  
  //I compute C=gamma1*G' where G' is the G conjugate/transpose
  
/*    C=cmatrix(0,N-1,0,N-1); */
/*    for (i=0;i<N;i++) */
/*      for (j=0;j<N;j++) */
/*        { */
/*  	C[i][j]=zero; */
/*  	for (l=0;l<N;l++) */
/*  	  { */
/*  	    C[i][j].r=C[i][j].r+gamma1temp[i][l]*Gtemp[j][l].r; */
/*  	    C[i][j].i=C[i][j].i-gamma1temp[i][l]*Gtemp[j][l].i; */
/*  	  } */
/*        } */

  //I compute Spectralfunction=Diagonal{G*C}/(2*pi)
  for (i=0;i<N;i++)
    {
      SP[i]=0;
      for (l=0;l<N;l++)
	SP[i]=SP[i]+G[i][l].r*C[l][i].r-G[i][l].i*C[l][i].i;
      SP[i]=SP[i]/(2*pi);
      //SP[i]=g1max*Gmax*Gmax*SP[i]/(2*pi);
    }
  cfree_cmatrix(C,0,N-1,0,N-1);
  cfree_cmatrix(Tr,0,N-1,0,N-1);
  //  free_dmatrix(gamma1temp,0,N-1,0,N-1);
  //  cfree_cmatrix(Gtemp,0,N-1,0,N-1);
/*    for (i=0;i<N;i++) */
/*      printf("%lg \n",SP[i]); */
  //  exit(1);
  return SP;
}
