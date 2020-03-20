// ======================================================================
//  Copyright (c) 2012, G. Fiori,  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "Fermi_Integrals.h"
static PyObject* py_Fermi_Integrals(PyObject* self, PyObject* args)
{
  PyObject *obj1,*obj2,*temp_obj;
  double eta,out,ord,relerr;
  int IERR;
  // I read the data from python class
  import_array();
  if (!PyArg_ParseTuple(args,"OO",&obj1,&obj2))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }
  relerr=1e-6;
  ord=(double)PyFloat_AsDouble(obj1);
  eta=(double)PyFloat_AsDouble(obj2);
  fermid_(&ord,&eta,&relerr,&out,&IERR);
  
  return Py_BuildValue("d", out);
}
