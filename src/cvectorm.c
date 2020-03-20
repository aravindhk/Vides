// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "cvectorm.h"
complex ***cvectorm(int nl, int nh)
/* allocate a complex vector with subscript range v[nl..nh] */
{
  complex ***v;
  v=(complex ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(complex**)));
  if (!v) printf("allocation failure in dvector()");
  return v-nl+NR_END;
}
