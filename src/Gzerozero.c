// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Gzerozero.h"
// In wmH you put the EI-H, where H is the Hamiltonian, and 
// E is the energy. In BETA and BETADAGA you put the 
// hopping terms, i.e. the up and downdiagonal elements of the
// hamiltonian. OUTPUT : surface Green's function.
complex **Gzerozero(complex **wmH,complex **BETA,complex **BETADAGA,int N)
     //complex **Gzerozero(double E,complex **H,complex **BETA,complex **BETADAGA,int N);
{
  complex **ti,**titilde,**omega,**ID,**temp,**T,**ZERO,**temp1,**temp2,
    **temp3,**temp4,**temp5,**tinew,**titildenew,**temp6,**titildecycle,
    zero,**Gz;
  int i,j,k,l,ix,counter;
  double normati,eta,normatilde;
  FILE *fp;
  zero.r=0;
  zero.i=0;
  //  eta=1e-5;
  eta=1e-25;
  
  T=cmatrix(0,N-1,0,N-1);

  ID=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	if (i==j)
	  ID[i][j]=complass(1,0);
	else
	  ID[i][j]=zero;
      }
  
  ZERO=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      ZERO[i][j]=zero;
  
  // I compute t0
  temp=cmatinv(wmH,N);
  ti=cmatmul(temp,BETADAGA,N);
  titilde=cmatmul(temp,BETA,N);
  cfree_cmatrix(temp,0,N-1,0,N-1);
  
  // I compute T
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      T[i][j]=ti[i][j];
  
  titildecycle=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      titildecycle[i][j]=titilde[i][j];

/*    fp=fopen("prova","w"); */
/*    for (i=0;i<N;i++) */
/*      { */
/*        for (j=0;j<N;j++) */
/*  	printf("%lg ",titilde[i][j].r); */
/*        printf("\n"); */
/*      } */
/*    close(fp); */

/*    exit(1); */

  normati=1;
  normatilde=1;
  counter=1;
  while ((normati>=1e-10)||(normatilde>=1e-10))
    {
      temp=cmatmul(ti,titilde,N);
      temp1=cmatmul(titilde,ti,N);
      temp2=cmatsub(ID,temp,N);
      temp3=cmatsub(temp2,temp1,N);
      temp4=cmatinv(temp3,N);
      temp5=cmatmul(ti,ti,N);
      temp6=cmatmul(titilde,titilde,N);
      tinew=cmatmul(temp4,temp5,N);
      titildenew=cmatmul(temp4,temp6,N);
      normati=cmatnorm2diff(ti,ZERO,N);
      normatilde=cmatnorm2diff(titilde,ZERO,N);

      cfree_cmatrix(temp,0,N-1,0,N-1);
      cfree_cmatrix(temp1,0,N-1,0,N-1);
      cfree_cmatrix(temp2,0,N-1,0,N-1);
      cfree_cmatrix(temp3,0,N-1,0,N-1);
      cfree_cmatrix(temp4,0,N-1,0,N-1);
      cfree_cmatrix(temp5,0,N-1,0,N-1);
      cfree_cmatrix(temp6,0,N-1,0,N-1);
      
      temp=cmatmul(titildecycle,tinew,N);
      temp2=cmatsum(T,temp,N);
      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  T[i][j]=temp2[i][j];
      
      temp3=cmatmul(titildecycle,titildenew,N);
      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  titildecycle[i][j]=temp3[i][j];

      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  {
	    ti[i][j]=tinew[i][j];
	    titilde[i][j]=titildenew[i][j];
	  }
      cfree_cmatrix(tinew,0,N-1,0,N-1);
      cfree_cmatrix(titildenew,0,N-1,0,N-1);
      cfree_cmatrix(temp,0,N-1,0,N-1);
      cfree_cmatrix(temp2,0,N-1,0,N-1);
      cfree_cmatrix(temp3,0,N-1,0,N-1);
      counter++;
    }

  //printf("Iterazioni %d \n",counter);
  // I compute Gzerozero

  temp1=cmatmul(BETA,T,N);
  temp2=cmatsub(wmH,temp1,N);
  temp3=cmatinv(temp2,N);
  Gz=cmatmul3(BETA,temp3,BETADAGA,N);
    
  cfree_cmatrix(temp1,0,N-1,0,N-1);
  cfree_cmatrix(temp2,0,N-1,0,N-1);
  cfree_cmatrix(temp3,0,N-1,0,N-1);
  cfree_cmatrix(T,0,N-1,0,N-1);
  cfree_cmatrix(ID,0,N-1,0,N-1);
  cfree_cmatrix(ZERO,0,N-1,0,N-1);
  cfree_cmatrix(titildecycle,0,N-1,0,N-1);
  cfree_cmatrix(ti,0,N-1,0,N-1);
  cfree_cmatrix(titilde,0,N-1,0,N-1);
  return Gz; 
}
