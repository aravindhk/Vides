// ======================================================================
//  Copyright (c) 2012, G. Fiori, P. Marconcini, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "writeout.h"
static PyObject* py_writeout(PyObject* self, PyObject* args)
{

  char *string2;
  PyObject *string;
  import_array();

  if (!PyArg_ParseTuple(args,"O",&string))
    {
      printf("NOT AN ARRAY! \n");
      exit(0);
    }
  #if PY_MAJOR_VERSION >=3
  string2=PyBytes_AsString(PyUnicode_AsEncodedString(string,"utf-8","Error"));
  #else
  string2=PyString_AsString(string);
  #endif
  printf(string2);
  printf("\n");
  char *s = "Hello from writeout!";
  return Py_BuildValue("s", s);
}
  
