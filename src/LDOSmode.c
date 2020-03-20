// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

/*      WE MAKE THE ASSUMPTION THAT THE HAMILTONIAN AS WELL AS    */
/*      THE SELF-ENERGY MATRICES ARE REPRESENTED IN THE MODE SPACE*/
/*      THE LDOS IS RETURNED IN THE REAL SPACE BASIS SET          */

/* ************************************************************** */
/*      This subroutine compute the LDOS                          */
/*      in case of ballistic transport.                           */
/*      diag, lowdiag, and updiag are the block diagonal matrices */
/*      of the threediagonal block matrix Hamiltonian.            */
/*      SigmaS and SigmaD are the self-consistent                 */
/*      energies of the source and of the drain,                  */
/*      respectively. A1,A2 are the spectral densities divided    */
/*      by 2pi, i.e. the Local Density Of States (LDOS) for       */
/*      the source and for the drain, respectively.               */
/*      In order to save computational time the -H is passed      */
/*      The spectral function A1 and A2 are returned in the       */
/* 	matrices A1 and A2, ordered by blocks. i.e. the           */
/* 	i-th row of A* is the spectral function associated        */
/*      with the i-th block and so on                             */
/* ************************************************************** */
#include "LDOSmode.h"
void LDOSMODE(double E,complex ***lowdiag,complex ***diag,complex ***updiag,
	      double ***A1,double ***A2,complex **sigmas,complex **sigmad,int n,
	      int Nc,int flagtrans,double *T,int Nreal,int *order,double thop,
	      double eta)
{
  int i,j,k,l,flaggo;
  complex ***d,****G,**temp1,**temp2;
  double *temp,xx,**transo,**helpo,**ID;
  complex **gamma1,**gamma2;
  FILE *fp;

  *A1=dmatrix(0,Nc-1,0,Nreal-1);
  *A2=dmatrix(0,Nc-1,0,Nreal-1);
  
  //Inizialization of useful vectors
  d=cvectorm(0,Nc-1);

  //  I compute the term E*I-H
  for (i=0;i<Nc;i++)
    {

      d[i]=cmatrix(0,n-1,0,n-1);
      for (j=0;j<n;j++)
	for (k=0;k<n;k++)
	  {
	    if (j==k)
	      d[i][j][k]=csum(complass(E,eta),diag[i][j][k]);
	    else
	      d[i][j][k]=diag[i][j][k];
	  }
    }


  
  // butto fuori Ginv
  
  d[0]=cmatsub(d[0],sigmas,n);
  d[Nc-1]=cmatsub(d[Nc-1],sigmad,n);
  
  //  I compute the term E*I-H-SigmaS-SigmaD
  
  G=rgfblock(lowdiag,d,updiag,n,Nc,1);
  
  // I compute gamma1 and gamma2
  // gamma=j*[sigma-sigma']
  
  // I compute sigma-sigma'
  temp1=cmatdaga(sigmas,n);
  gamma1=cmatsub(sigmas,temp1,n);
  cfree_cmatrix(temp1,0,n-1,0,n-1);
  
  temp1=cmatdaga(sigmad,n);
  gamma2=cmatsub(sigmad,temp1,n);
  cfree_cmatrix(temp1,0,n-1,0,n-1);
  
  // I compute j*[sigma-sigma']
  
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      {
	gamma1[i][j]=complass(-gamma1[i][j].i,gamma1[i][j].r);
	gamma2[i][j]=complass(-gamma2[i][j].i,gamma2[i][j].r);
      }

  for (i=0;i<Nc;i++)
    {
      temp=spectralfunmode(G[i][0],gamma1,n,Nreal,order);
      for (j=0;j<Nreal;j++)
	(*A1)[i][j]=temp[j];
      free_dvector(temp,0,Nreal-1);
      temp=spectralfunmode(G[i][1],gamma2,n,Nreal,order);
      for (j=0;j<Nreal;j++)
	(*A2)[i][j]=temp[j];
      free_dvector(temp,0,Nreal-1);
    }
  
  // if flagtrans eq 1 then I compute the transmission coefficient
  // in case of ballistic transport you need 
  if (flagtrans==1)
    {
      *T=transmission(G[Nc-1][0],gamma1,gamma2,n);
    }
  
  //Free the memory
  cfree_ctensor4(G,0,Nc-1,0,1,0,n-1,0,n-1,0);
  cfree_cvectorm(d,0,Nc-1,0,n-1,0,n-1);
  cfree_cmatrix(gamma1,0,n-1,0,n-1);
  cfree_cmatrix(gamma2,0,n-1,0,n-1);
  return;
}
