// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatvect.h"
complex *cmatvect(complex **A,complex *v,int N)
{
  int i,j,l;
  complex *C,zero,temp;
  zero.r=0;
  zero.i=0;
  C=ccvector(0,N-1);
  for (i=0;i<N;i++)
    {
      C[i]=zero;
      for (l=0;l<N;l++)
	if ((fabs(v[l].r)>0)||(fabs(v[l].i)>0))
	  {
	    C[i].r=C[i].r+A[i][l].r*v[l].r-A[i][l].i*v[l].i;
	    C[i].i=C[i].i+A[i][l].i*v[l].r+A[i][l].r*v[l].i;
	  }
    }
  return C;
}
	
