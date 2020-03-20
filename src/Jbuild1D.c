// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Jbuild1D.h"
void Jbuild1D(physical_quantities p,device_mapping d,
	      double ***J,int neq)
{
  double **Jel;
  int i,j,k,l,m,ix,jjj,index;
  for (i=0;i<d.nx;i++)
    {
      ix=i;
      for (l=0;l<neq;l++)
	for (m=0;m<3*neq;m++)
	  J[i][l][m]=0.;
	  
      index=d.boundary_conditions[ix];
      // ************************************************************************
      // ********************* PUNTI INTERNI ************************************
      // ************************************************************************
      if ((index==2000)||((index>=1003)))
	{
	  // POISSON 
	  J[i][0][0]=-(1.0/d.dist[ix][0]*eps0*(p.eps[ix-1]+p.eps[ix])*0.5+
		       1.0/d.dist[ix][1]*eps0*(p.eps[ix+1]+p.eps[ix])*0.5)
	    +q*(+p.free_charge[ix])*exp(-sign(p.free_charge[ix])*(p.Phi[ix]-p.Phiold[ix])/(*p.vt))/(*p.vt)*(-sign(p.free_charge[ix]));

	  J[i][0][1]=1.0/d.dist[ix][0]*eps0*(p.eps[ix-1]+p.eps[ix])*0.5;
	  J[i][0][2]=1.0/d.dist[ix][1]*eps0*(p.eps[ix+1]+p.eps[ix])*0.5;
	}
      
      // ************************************************************************
      // ****************** PUNTI DI NEUMANN ************************************
      // ************************************************************************
      
      
      if ((index>=1001)&&(index<=1002))
	{ 
	  J[i][0][0]=q; 
	  if (index==1001)
	    {
	      J[i][0][2]=-q; 
	    }
	  if (index==1002)
	    {
	      J[i][0][1]=-q;
	    }
	}
      
      // ************************************************************************
      // ****************** PUNTI DI DIRICHLET **********************************
      // ************************************************************************
      
      if (index<1000)
	{
	  J[i][0][0]=1;
	}
    }

  return;
}
