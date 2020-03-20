// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "create_updiagGNR.h"
complex **create_updiagGNR(int Nc,int n,double thop)
{
  complex zero,t,**out;
  int i,j,tempo;
  
  tempo=Nc%4;
  if ((tempo==0)||(tempo==2))
    out=create_beta1GNR(n,Nc,thop);
  else if (tempo==1)
    out=create_beta2GNR(n,thop);
  else if (tempo==3)
    out=create_beta2transpGNR(n,thop);
  return out;
}
