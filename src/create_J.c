// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_J.h"
void create_J(double ******J,int nx,int ny,int nz,int neq)
{
  int i,j,k,n1,n2,n3,n4;
  double **Jel,****JJ;
  n1=0;
  n2=0;
  n3=0;
  n4=0;
  // cerco il # massimo di autovalori per ogni zona
  //printf("Memory allocation ...... \n");
  *J=dvectormm(0,nx-1);
  for (i=0;i<nx;i++)
    {
      JJ=dmatrixm(0,ny-1,0,nz-1);
      for (k=0;k<nz;k++)
	for (j=0;j<ny;j++)
	  {
	    Jel=dmatrix(0,neq-1,0,7*neq-1);
	    JJ[j][k]=Jel;
	  }	 
      (*J)[i]=JJ;
    } 
  return;
}
