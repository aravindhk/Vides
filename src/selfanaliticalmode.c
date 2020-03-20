// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "selfanaliticalmode.h"
// E is the energy (J), alpha is a parameter used by Jing Guo in
// his phenomenological approach, Em1 is the vector of the
// potential in the CNT at the contact, N is the number of 
// C atoms in each ring
complex **selfanaliticalmode(double E,double Em1,int N,int Nm,int *order,
			     double t,double eta)
{
  complex **EImEm1,zero,temp,**temp1,**g1,b,delta,a,**out;
  complex **beta,**betadaga,g1num,g1den;
  complex **beta1,**beta2,**beta2d,**H1;
  int i,j,k,l,ix,counter;
  double norma;
  FILE *fp;
  zero.r=0;
  zero.i=0;
  
  g1=cmatrix(0,Nm-1,0,Nm-1);
  EImEm1=cmatrix(0,Nm-1,0,Nm-1);
  
  for (i=0;i<Nm;i++)
    EImEm1[i][i]=complass(E-Em1,eta);
  
  beta1=create_beta1(N,t);
  beta2=create_beta2mode(N,Nm,order,t);
  beta2d=create_beta2transpmode(N,Nm,order,t);


  for (i=0;i<Nm;i++)
    for (j=0;j<Nm;j++)
      {
	if (i!=j)
	  g1[i][j]=zero;
	else
	  {
	    // g1=(b-sqrt(delta))/(2*a)
	    b=ccmul(EImEm1[i][i],EImEm1[i][i]);
	    temp=ccmul(beta2[i][i],beta2d[i][i]);
	    b=ccadd(b,temp);
	    temp=ccmul(beta1[i][i],beta1[i][i]);
	    b=ccsub(b,temp);
	    delta=ccmul(b,b);
	    temp=ccmul(EImEm1[i][i],EImEm1[i][i]);
	    temp=ccmul(temp,beta2[i][i]);
	    temp=ccmul(temp,beta2d[i][i]);
	    temp=cdmul(-4,temp);
	    delta=ccadd(delta,temp);
	    delta=ccsqrt(delta);
	    g1den=ccmul(EImEm1[i][i],beta2[i][i]);
	    g1den=ccmul(g1den,beta2d[i][i]);
	    g1den=cdmul(2,g1den);
	    temp=complass(0,1e-8);
	    g1den=ccadd(g1den,temp);
	    if ((cIm(delta)/cRe(g1den))<0)
	      g1num=ccadd(b,delta);
	    else
	      g1num=ccsub(b,delta);
	    g1[i][i]=ccdiv(g1num,g1den);
	  }
      }	

  // I return the selfenergy = beta2d*g1*beta2d
  for (i=0;i<Nm;i++)
    {
      temp=ccmul(beta2d[i][i],beta2[i][i]);
      g1[i][i]=ccmul(g1[i][i],temp);
    }
  
  cfree_cmatrix(beta1,0,Nm-1,0,Nm-1);
  cfree_cmatrix(beta2,0,Nm-1,0,Nm-1);
  cfree_cmatrix(beta2d,0,Nm-1,0,Nm-1);
  cfree_cmatrix(EImEm1,0,Nm-1,0,Nm-1);
  return g1; 
}
