// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "gappo.h"
#include <stdio.h>
double gappo(int n,double thop)
{
  double Eg3p,Eg3p_1,Eg3p_2,out,delta;
  int i,j,p;
  delta=0.12;
  p=2*n/3;
  Eg3p=thop*(4*cos(p*pi/(3*p+1))-2);
  Eg3p_1=thop*(2-4*cos((p+1)*pi/(3*p+2)));
  Eg3p_2=0;
  if ((2*n)==(3*p))
    out=Eg3p-8*delta*thop/(3*p+1)*sin(p*pi/(3*p+1))*sin(p*pi/(3*p+1));
  else if ((2*n)==(3*p+1))
    out=Eg3p_1+8*delta*thop/(3*p+2)*sin((p+1)*pi/(3*p+2))*sin((p+1)*pi/(3*p+2));
  else if ((2*n)==(3*p+2))
    out=Eg3p_2+2*delta*thop/(p+1);
  return out;
}
