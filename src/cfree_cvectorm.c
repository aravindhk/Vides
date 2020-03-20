// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cfree_cvectorm.h"
void cfree_cvectorm(complex ***v,int nvl,int nvh,int nrl, int nrh,int ncl,int nch)
{
  int i,j,k;
  for(i=0;i<(nvh+1);i++)
    cfree_cmatrix(v[i],nrl,nrh,ncl,nch);
  free((FREE_ARG) (v-NR_END));
}
