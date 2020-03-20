// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "ctensor4.h"
complex ****ctensor4(int n1l,int n1h,int n2l,int n2h,int n3l,int n3h, int n4l, int n4h)
{
  complex ****temp3,**temp2;
  int i,j;
  temp3=cmatrixm(n1l,n1h,n2l,n2h);
/*   for (i=0;i<=n1h;i++) */
/*     for (j=0;j<=n2h;j++) */
/*       { */
/* 	temp2=cmatrix(n3l,n3h,n4l,n4h); */
/* 	temp3[i][j]=temp2; */
/*       } */
  return temp3;
}
