// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "energysort.h"
void energysort(int **help,int n,double thop)
{
  int flag,*order,i;
  double *x,*y;
  flag=2;
  x=dvector(0,n-1);
  y=dvector(0,n-1);
  order=ivector(0,n-1);

  for (i=0;i<n;i++)
    {
      x[i]=fabs(2*fabs(thop)*cos(pi*(i)/n)-fabs(thop));
      y[i]=i;
    }

  dsort_(x,y,&n,&flag);
  for (i=0;i<n;i++)
    {
      order[i]=(int)(y[i]+0.5);
      (*help)[i]=order[i];
    }

  for (i=0;i<n;i++)
    printf("%d ",order[i]);
   printf("\n");

  for (i=0;i<(int)(n*0.5+0.5);i++)
    {
      (*help)[2*i]=order[i];
    }

  for (i=0;i<(int)(n*0.5);i++)
    {
      if (order[i]!=0)
	(*help)[2*i+1]=n-order[i];
      else
	(*help)[2*i+1]=(int)(n*0.5);
    }

  free_dvector(x,0,n-1);
  free_dvector(y,0,n-1);
  free_ivector(order,0,n-1);
  return;
}
  
