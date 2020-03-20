// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "norma2.h"
double norma2(double *vett,int nnn)
{
  int i;
  double norma;
  norma=0.0;
  for (i=0;i<nnn;i++)
    norma+=pow(vett[i],2);
  norma=sqrt(norma);
  return norma;
}
