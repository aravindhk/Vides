// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "create_lowdiag.h"
complex **create_lowdiag(int Nc,int n,double thop)
{
  complex zero,t,**out;
  int i,j,tempo;
  
  tempo=Nc%4;
  if ((tempo==0)||(tempo==2))
    out=create_beta1(n,thop);
  else if (tempo==1)
    out=create_beta2transp(n,thop);
  else if (tempo==3)
    out=create_beta2(n,thop);
  return out;
}
