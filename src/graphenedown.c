// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "graphenedown.h"
void graphenedown(int i,int j,double *y,double *z,int n,double acc)
{
  int indz;
  if (((j%4)==0)||((j%4)==1))
    *y=i*sqrt(3)*acc;
  else
    *y=i*sqrt(3)*acc+sqrt(3)*acc*0.5;
  
  indz=floor(j/4);
  if ((j%4)==0)
    *z=indz*3*acc;
  else if ((j%4)==1)
    *z=indz*3*acc+acc;
  else if ((j%4)==2)
    *z=indz*3*acc+1.5*acc;
  else if ((j%4)==3)
    *z=indz*3*acc+2.5*acc;
  return;
}
