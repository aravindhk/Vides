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

#include "spectralfunmode.h"
double *spectralfunmode(complex **G,complex **gamma1,int N,int NReal,int *order)
{
  double *SP,g1max,Gmax,**gamma1temp;
  complex zero,**Gtemp;
  complex **C;
  complex **Tr,**TrReal,**GReal,**gamma1Real,**dummy;
  int i,j,l;
  SP=dvector(0,NReal-1);
  zero.r=0;
  zero.i=0;

  Tr=cmatrix(0,N-1,0,N-1);
  for (i=0;i<N;i++) 
    for (j=0;j<N;j++)
      Tr[i][j]=complass(G[j][i].r,-G[j][i].i);

/*   // I Anti-Transform G and Gdaga (Tr) and gamma1 */
/*   TrReal=VAVdaga(Tr,NReal,N,order); */
/*   GReal=VAVdaga(G,NReal,N,order); */
/*   gamma1Real=VAVdaga(gamma1,NReal,N,order); */

  C=cmatmul3(G,gamma1,Tr,N);
  dummy=VAVdaga(C,NReal,N,order);

  //I compute Spectralfunction=Diagonal{G*C}/(2*pi)
  for (i=0;i<NReal;i++)
    {
      SP[i]=dummy[i][i].r/(2*pi);
    }
  cfree_cmatrix(C,0,N-1,0,N-1);
  cfree_cmatrix(Tr,0,N-1,0,N-1);
  cfree_cmatrix(dummy,0,NReal-1,0,NReal-1);
  return SP;
}
