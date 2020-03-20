// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
/* *********************************************** */
/*       This subroutine deallocate a complex      */
/*       4 tensor. If flagball==1 then deallocate  */
/*       the diagonal and the first and the last   */
/*       column, otherwise deallocate the whole    */
/*       matrix.                                   */
/* *********************************************** */
#include "cfree_ctensor4.h"
void cfree_ctensor4(complex ****m,int n1l,int n1h,int n2l,int n2h,int n3l,int n3h, int n4l, int n4h,int flagball)
/* free a complex matrix allocated by dmatrix() */
{
  int i,j;
  if (flagball!=1)
    {
      for (i=0;i<=n1h;i++)
	for (j=0;j<=n2h;j++)
	  cfree_cmatrix(m[i][j],n3l,n3h,n4l,n4h);
      free((FREE_ARG) (m[n1l]+n2l-NR_END));
      free((FREE_ARG) (m+n1l-NR_END));
    }
  else
    {
      for (i=0;i<=n1h;i++)
	{
	  cfree_cmatrix(m[i][i],n3l,n3h,n4l,n4h);
	  if (i!=n2l) cfree_cmatrix(m[i][n2l],n3l,n3h,n4l,n4h);
	  if (i!=n2h) cfree_cmatrix(m[i][n2h],n3l,n3h,n4l,n4h);
	}
      free((FREE_ARG) (m[n1l]+n2l-NR_END));
      free((FREE_ARG) (m+n1l-NR_END));
    }
  return;
}
