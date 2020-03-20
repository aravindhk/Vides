// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_J2D.h"
void create_J2D(double *****J,int nx,int ny,int neq)
{
  int i,j,k,n1,n2,n3,n4;
  double **Jel;
  n1=0;
  n2=0;
  n3=0;
  n4=0;
  // cerco il # massimo di autovalori per ogni zona
  //  printf("Memory allocation ...... \n");
  *J=dmatrixm(0,nx-1,0,ny-1);
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      {
	Jel=dmatrix(0,neq-1,0,7*neq-1);
	(*J)[i][j]=Jel;
      } 
  return;
}
