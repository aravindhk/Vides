// ======================================================================
//  Copyright (c) 2004-2008, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatmul3.h"
#include "cmatmul.h"
complex **cmatmul3(complex **A,complex **B,complex **C,int N)
{
  int i,j,l;
  complex **D,**M,zero,temp;
  zero.r=0;
  zero.i=0;

  M=cmatmul(B,C,N);
  D=cmatmul(A,M,N);
  cfree_cmatrix(M,0,N-1,0,N-1);
  return D;
}
	
