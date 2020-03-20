// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 


#include "preparej.h"
void preparej(physical_quantities p,device_mapping d,
	      double *a2,int *ia2,int *ja2,int *z12,int flagum,
	      double *****J,int neq)
{
  int i,j,k,ix,eq,l,m,jj,passuzzo,Np;
  double template;
  double **dummymatrix;
  double ****dummymatrix2;
  FILE *fpp;
  Np=d.nx*d.ny*d.nz;
  // creo il vettore dummy J
  Jbuild(p,d,J,neq);
  
  m=0;
  *z12=0;
  jj=0;
  eq=flagum;
  if (flagum==0)
    {
      for (k=0;k<d.nz;k++)
	for (j=0;j<d.ny;j++)
	  for (i=0;i<d.nx;i++)
	    {
	      ix=i+j*d.nx+k*d.nx*d.ny;
	      for (l=7*flagum;l<7*(flagum+1);l++)
		{
		  template=J[i][j][k][eq][l];
		  if (abbs(template)>1e-40)
		    {
		      a2[jj]=J[i][j][k][eq][l];
		      ja2[jj]=ix+indice(l,d.nx,d.ny,d.nz)-l/7*Np+1;
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
