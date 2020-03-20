// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 


#include "preparej2D.h"
void preparej2D(physical_quantities p,device_mapping d,
		double *a2,int *ia2,int *ja2,int *z12,int flagum,
		double ****J,int neq)
{
  int i,j,k,ix,eq,l,m,jj,passuzzo,Np;
  double template;
  double **dummymatrix;
  double ****dummymatrix2;
  FILE *fpp;
  Np=d.nx*d.ny;
  // creo il vettore dummy J
  Jbuild2D(p,d,J,neq);
  
  m=0;
  *z12=0;
  jj=0;
  eq=flagum;
  if (flagum==0)
    {
      for (j=0;j<d.ny;j++)
	for (i=0;i<d.nx;i++)
	  {
	    ix=i+j*d.nx;
	    for (l=5*flagum;l<5*(flagum+1);l++)
	      {
		template=J[i][j][eq][l];
		if (abbs(template)>1e-40)
		  {
		    a2[jj]=J[i][j][eq][l];
		    ja2[jj]=ix+indice(l,d.nx,d.ny,d.nz)-l/5*Np+1;
		    *z12=*z12+1;
		    jj++;
		  }
	      }
	    ia2[m]=*z12+1;
	    m++;
	  }
    } 
  // adesso libero la memoria allocata in J dagli elementi Jel
  return;
}
