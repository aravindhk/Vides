// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "stampaout.h"
void stampaout(char nome[20],double *ppp,int Np)
{
  FILE *fp3;
  int i;
  fp3=fopen(nome,"w");
  for (i=0;i<Np;i++)
    fprintf(fp3,"%.10lg \n",ppp[i]);
  fclose(fp3);
  return;
}

