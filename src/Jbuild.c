// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Jbuild.h"
void Jbuild(physical_quantities p,device_mapping d,
	    double *****J,int neq)
{
  double **Jel,****JJ,dk,a1;
  int i,j,k,l,m,ix,jjj,index;
  for (i=0;i<d.nx;i++)
    {
      for (k=0;k<d.nz;k++)
	for (j=0;j<d.ny;j++)
	  {
	    ix=i+j*d.nx+k*d.nx*d.ny;
	    for (l=0;l<neq;l++)
	      for (m=0;m<7*neq;m++)
		J[i][j][k][l][m]=0.;
	     
	    index=d.boundary_conditions[ix];
	    // ************************************************************************
	    // ********************* PUNTI INTERNI ************************************
	    // ************************************************************************
	    if ((index==2000)||((index>=1007)))
	      {
		// POISSON 
		J[i][j][k][0][0]=-(d.surf[ix][0]/d.dist[ix][0]*eps0*(p.eps[ix-1]+p.eps[ix])*0.5+
				   d.surf[ix][1]/d.dist[ix][1]*eps0*(p.eps[ix+1]+p.eps[ix])*0.5+
				   d.surf[ix][3]/d.dist[ix][3]*eps0*(p.eps[ix+d.nx]+p.eps[ix])*0.5+
				   d.surf[ix][2]/d.dist[ix][2]*eps0*(p.eps[ix-d.nx]+p.eps[ix])*0.5+
				   d.surf[ix][4]/d.dist[ix][4]*eps0*(p.eps[ix-d.nx*d.ny]+p.eps[ix])*0.5+
				   d.surf[ix][5]/d.dist[ix][5]*eps0*(p.eps[ix+d.nx*d.ny]+p.eps[ix])*0.5)
		  +q*(+p.free_charge[ix])*exp(-sign(p.free_charge[ix])*(p.Phi[ix]-p.Phiold[ix])/(*p.vt))/(*p.vt)*(-sign(p.free_charge[ix]));
		//+q*(+p.free_charge[ix])*exp((p.Phi[ix]-p.Phiold[ix])/(*p.vt))/(*p.vt);
		  //+q*dncntf(ix,p.free_charge[ix],p);
		J[i][j][k][0][1]=d.surf[ix][0]/d.dist[ix][0]*eps0*(p.eps[ix-1]+p.eps[ix])*0.5;
		J[i][j][k][0][2]=d.surf[ix][1]/d.dist[ix][1]*eps0*(p.eps[ix+1]+p.eps[ix])*0.5;
		J[i][j][k][0][3]=d.surf[ix][2]/d.dist[ix][2]*eps0*(p.eps[ix-d.nx]+p.eps[ix])*0.5;
		J[i][j][k][0][4]=d.surf[ix][3]/d.dist[ix][3]*eps0*(p.eps[ix+d.nx]+p.eps[ix])*0.5;
		J[i][j][k][0][5]=d.surf[ix][4]/d.dist[ix][4]*eps0*(p.eps[ix-d.nx*d.ny]+p.eps[ix])*0.5;
		J[i][j][k][0][6]=d.surf[ix][5]/d.dist[ix][5]*eps0*(p.eps[ix+d.nx*d.ny]+p.eps[ix])*0.5;
	      }
	    
	    // ************************************************************************
	    // ****************** PUNTI DI NEUMANN ************************************
	    // ************************************************************************
	    
	    
	    if ((index>=1001)&&(index<=1006))
	      { 
		J[i][j][k][0][0]=q; 
		if (index==1001)
		  {
		    J[i][j][k][0][2]=-q; 
		  }
		if (index==1002)
		  {
		    J[i][j][k][0][1]=-q;
		  }
		if (index==1003)
		  {
		    J[i][j][k][0][4]=-q; 
		  }
		if (index==1004)
		  {
		    J[i][j][k][0][3]=-q;
		  }
		if (index==1005)
		  {
		    J[i][j][k][0][6]=-q;
		  }
		if (index==1006)
		  {
		    J[i][j][k][0][5]=-q;
		  }
	      }
	     
	    // ************************************************************************
	    // ****************** PUNTI DI DIRICHLET **********************************
	    // ************************************************************************
	     
	    if (index<1000)
	      {
		J[i][j][k][0][0]=q;
	      }
	  }
    }
  return;
}
