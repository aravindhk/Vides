// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cmatrixm.h"
complex ****cmatrixm(int nrl, int nrh, int ncl, int nch)
     /* allocate a complex matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  complex ****m;

  /* allocate pointers to rows */
  m=(complex ****) malloc((size_t)((nrow+NR_END)*sizeof(complex***)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(complex ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(complex**)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
