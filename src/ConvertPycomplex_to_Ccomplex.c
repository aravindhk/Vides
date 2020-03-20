// ======================================================================
//  Copyright (c) 2011, G. Fiori, P. D'Amico, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "ConvertPycomplex_to_Ccomplex.h"
complex ConvertPycomplex_to_Ccomplex(Py_complex c)
{
  complex out;
  out.r=c.real;
  out.i=c.imag;
  return out;
}

