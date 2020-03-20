// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "indice.h"
int indice(int a,int nx,int ny,int nz)
{
  int ora,Np;
  Np=nx*ny*nz;
  if ((a%7)==0)
    ora=0;
  else if ((a%7)==1)
    ora=-1;
  else if ((a%7)==2)
    ora=1;
  else if ((a%7)==3)
    ora=-nx;
  else if ((a%7)==4)
    ora=+nx;
  else if ((a%7)==5)
    ora=-nx*ny;
  else if ((a%7)==6)
    ora=+nx*ny;
  return ora+a/7*Np;
}
