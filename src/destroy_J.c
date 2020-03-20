// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "destroy_J.h"
void destroy_J(double *****J,int nx,int ny,int nz,int neq)
{
  int i,j,k;
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
	free_dmatrix(J[i][j][k],0,neq-1,0,7*neq-1);
  free((FREE_ARG) (J[0][0]-NR_END));
  free((FREE_ARG) (J[0]-NR_END));
  free((FREE_ARG) (J-NR_END));
  return;
}
