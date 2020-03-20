// ======================================================================
//  Copyright (c) 2012, G. Fiori,  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "fphalf.h"
static PyObject* py_fphalf(PyObject* self, PyObject* args)
{
  PyObject *obj,*temp_obj;
  double eta,out;
  // I read the data from python class
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }
  eta=(double)PyFloat_AsDouble(obj);
  out=sub_fphalf_(&eta);

  return Py_BuildValue("d", out);
}
