// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "indice2D.h"
int indice2D(int a,int nx,int ny)
{
  int ora,Np;
  Np=nx*ny;
  if ((a%5)==0)
    ora=0;
  else if ((a%5)==1)
    ora=-1;
  else if ((a%5)==2)
    ora=1;
  else if ((a%5)==3)
    ora=-nx;
  else if ((a%5)==4)
    ora=+nx;
  return ora+a/5*Np;
}
