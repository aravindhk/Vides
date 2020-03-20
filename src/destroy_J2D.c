// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "destroy_J2D.h"
void destroy_J2D(double ****J,int nx,int ny,int neq)
{
  int i,j,k;
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      free_dmatrix(J[i][j],0,neq-1,0,7*neq-1);
  free((FREE_ARG) (J[0]-NR_END));
  free((FREE_ARG) (J-NR_END));
  return;
}
