// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 


#include "ffunpoisson2D.h"
void ffunpoisson2D(physical_quantities p,device_mapping d,
		   double *bb)
{
  FILE *fp;
  int i,Np,index;
  double ora;
  Np=d.nx*d.ny;
  // ricorda che bb e' cambiata di segno
  for (i=0;i<Np;i++)
    {
      index=(int)d.boundary_conditions[i];
      if (index>1004)
	{
	  bb[i]=(
		 +(d.dist[i][2]+d.dist[i][3])*0.5/d.dist[i][0]*eps0*(p.Phi[i-1]-p.Phi[i])*(p.eps[i-1]+p.eps[i])*0.5
		 +(d.dist[i][2]+d.dist[i][3])*0.5/d.dist[i][1]*eps0*(p.Phi[i+1]-p.Phi[i])*(p.eps[i+1]+p.eps[i])*0.5
		 +(d.dist[i][0]+d.dist[i][1])*0.5/d.dist[i][2]*eps0*(p.Phi[i-d.nx]-p.Phi[i])*(p.eps[i-d.nx]+p.eps[i])*0.5
		 +(d.dist[i][0]+d.dist[i][1])*0.5/d.dist[i][3]*eps0*(p.Phi[i+d.nx]-p.Phi[i])*(p.eps[i+d.nx]+p.eps[i])*0.5
		 );
	  bb[i]+=((q*(+p.free_charge[i]))*exp(-sign(p.free_charge[i])*(p.Phi[i]-p.Phiold[i])/(*p.vt))+q*p.fixed_charge[i]);
	  //bb[i]+=(q*ncntf(i,p.free_charge[i],p))+q*p.fixed_charge[i]*d.dVe[i];
	}
      else if (index==1001) 
	{
	  bb[i]=q*(p.Phi[i]-p.Phi[i+1]);
	}
      else if (index==1002) 
	{
	  bb[i]=q*(p.Phi[i]-p.Phi[i-1]);
	}
      else if (index==1003) 
	{
	  bb[i]=q*(p.Phi[i]-p.Phi[i+d.nx]);
	}
      else if (index==1004) 
	{
	  bb[i]=q*(p.Phi[i]-p.Phi[i-d.nx]);
	}
      else if (index<1000)
	bb[i]=q*(p.Phi[i]+d.boundary_conditions[i]);
      
      if (bb[i]!=0)
	bb[i]=-bb[i];
    }
  return;
}
