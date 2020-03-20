// ======================================================================
//  Copyright (c) 2011, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 


#include "preparej1D.h"
void preparej1D(physical_quantities p,device_mapping d,
		double *a2,int *ia2,int *ja2,int *z12,int flagum,
		double ***J,int neq)
{
  int i,j,k,ix,eq,l,m,jj,passuzzo,Np;
  double template;
  double **dummymatrix;
  double ****dummymatrix2;
  FILE *fpp;
  Np=d.nx;
  // creo il vettore dummy J
  Jbuild1D(p,d,J,neq);
  
  m=0;
  *z12=0;
  jj=0;
  eq=flagum;
  if (flagum==0)
    {
      for (i=0;i<d.nx;i++)
	{
	  ix=i;
	  for (l=3*flagum;l<3*(flagum+1);l++)
	    {
	      template=J[i][eq][l];
	      if (abbs(template)>1e-40)
		{
		  a2[jj]=J[i][eq][l];
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
