// ======================================================================
//  Copyright (c) 2004-2008, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "free_dvectorm.h"
void free_dvectorm(double ***v,int nvl,int nvh,int nrl, int nrh,int ncl,int nch)
{
  int i;
  for(i=0;i<nvh;i++)
    free_dmatrix(v[i],nrl,nrh,ncl,nch);
  free((FREE_ARG) (v-NR_END));
}
