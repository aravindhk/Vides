// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "indice1D.h"
int indice1D(int a,int nx)
{
  int ora,Np;
  Np=nx;
  if ((a%3)==0)
    ora=0;
  else if ((a%3)==1)
    ora=-1;
  else if ((a%3)==2)
    ora=1;
  return ora+a/3*Np;
}
