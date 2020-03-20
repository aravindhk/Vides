// ======================================================================
//  Copyright (c) 2009, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "GNRgap.h"
static PyObject* py_GNRgap(PyObject* self, PyObject* args)
{
  PyObject *obj,*temp_obj;
  double Eg3p,Eg3p_1,Eg3p_2,out,delta;
  int i,j,p,n;
  double thop;
  // I read the data from python class
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }
  temp_obj=PyObject_GetAttrString(obj,"n");
  n=(int)PyInt_AsLong(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"thop");
  thop=(double)PyFloat_AsDouble(temp_obj);
  thop=fabs(thop);
  delta=0.12;
  p=2*n/3;
  Eg3p=thop*(4*cos(p*pi/(3*p+1))-2);
  Eg3p_1=thop*(2-4*cos((p+1)*pi/(3*p+2)));
  Eg3p_2=0;
  if ((2*n)==(3*p))
    out=Eg3p-8*delta*thop/(3*p+1)*sin(p*pi/(3*p+1))*sin(p*pi/(3*p+1));
  else if ((2*n)==(3*p+1))
    out=Eg3p_1+8*delta*thop/(3*p+2)*sin((p+1)*pi/(3*p+2))*sin((p+1)*pi/(3*p+2));
  else if ((2*n)==(3*p+2))
    out=Eg3p_2+2*delta*thop/(p+1);
  //  char *s = "Hello from CNT_charge_T";
  return Py_BuildValue("d", out);
}
