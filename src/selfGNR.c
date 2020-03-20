// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "selfGNR.h"
// E is the energy (J), alpha is a parameter used by Jing Guo in
// his phenomenological approach, Em1 is the vector of the
// potential in the CNT at the contact, N is the number of 
// C atoms in each ring
complex **selfGNR(double E,double *Em1,int N,double thop,double eta)
{
  complex **Gs,**EImEm1,**beta1,zero,**temp,**temp1;
  complex **Gsnew,**beta,**betadaga;
  complex **BETA,**BETADAGA,**H1;
  int i,j,k,l,ix,counter;
  double norma;
  FILE *fp;
  zero.r=0;
  zero.i=0;

  //I create the E*I-Em1*I term
  EImEm1=cmatrix(0,4*N-1,0,4*N-1);

  // I create the diagonal the contact Hamiltonian
  for (i=0;i<4*N;i++)
    for (j=0;j<4*N;j++)
      {
	EImEm1[i][j]=zero;
      }

  // ****** I create the up_diag of the 4N x 4N matrix
  for (i=1;i<4;i++)
    {
      temp=create_updiagGNR(i+1,N,thop);
      for (j=0;j<N;j++)
	for (k=0;k<N;k++)
	  EImEm1[j+(i-1)*N][k+i*N]=temp[j][k];
      cfree_cmatrix(temp,0,N-1,0,N-1);
    }

  // ****** I create the low_diag of the 4N x 4N matrix
  for (i=0;i<3;i++)
    {
      temp=create_lowdiagGNR(i+2,N,thop);
      for (j=0;j<N;j++)
	for (k=0;k<N;k++)
	  EImEm1[j+(i+1)*N][k+i*N]=temp[j][k];
      cfree_cmatrix(temp,0,N-1,0,N-1);
    }

  // I create the diagonal element of the 4N x 4N matrix
  for (j=0;j<4;j++)
    for (i=0;i<N;i++)
      EImEm1[i+N*j][i+N*j]=complass(E+Em1[i+(3-j)*N],eta);
  

  // Initial guess for Gs
  //  Gs=cmatinv(EImEm1,4*N);


  // I create beta* the tight-binding matrices (coupling 
  // between rings
  beta=create_beta2GNR(N,thop);
  betadaga=create_beta2transpGNR(N,thop);
  BETA=cmatrix(0,4*N-1,0,4*N-1);
  BETADAGA=cmatrix(0,4*N-1,0,4*N-1);
  for (i=0;i<4*N;i++)
    for (j=0;j<4*N;j++)
      {
	BETA[i][j]=zero;
	BETADAGA[i][j]=zero;
      }
  
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	BETA[3*N+i][j]=beta[i][j];
	BETADAGA[i][j+3*N]=betadaga[i][j];
      }
  
  Gsnew=Gzerozero(EImEm1,BETA,BETADAGA,4*N);
  
  // we only need the small part of the Gsnew
  temp=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      temp[i][j]=Gsnew[(4-1)*N+i][(4-1)*N+j];

  
  //  cfree_cmatrix(Gs,0,4*N-1,0,4*N-1);
  cfree_cmatrix(Gsnew,0,4*N-1,0,4*N-1);
  cfree_cmatrix(EImEm1,0,4*N-1,0,4*N-1);
  cfree_cmatrix(BETA,0,4*N-1,0,4*N-1);
  cfree_cmatrix(BETADAGA,0,4*N-1,0,4*N-1);
  cfree_cmatrix(beta,0,N-1,0,N-1);
  cfree_cmatrix(betadaga,0,N-1,0,N-1);
  return temp; 
  //return Gsnew;
}
