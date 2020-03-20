// ======================================================================
//  Copyright (c) 2011, G. Fiori, P. Marconcini, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "Simatrix.h"

static PyObject* py_Simatrix(PyObject* self, PyObject* args)
{

  int i, j, rank,Npar,Nh;
  double *par,dim[1],*h,cosx,cosy,cosz;

  PyObject *obj;
  PyArrayObject *temp, *out,*O_cosx,*O_cosy,*O_cosz;

  //  Py_complex cc;
  //  complex ccc;

  Npar=14;
  Nh=100;
  
  setvbuf(stdout,(char *) NULL,_IONBF,1);
  
  
/*
  Reads the data from python class.
*/
  
  import_array();
  if (!PyArg_ParseTuple(args,"OOOO",&O_cosx,&O_cosy,&O_cosz,&temp))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }
  
  cosx=PyFloat_AS_DOUBLE(O_cosx);
  cosy=PyFloat_AS_DOUBLE(O_cosy);
  cosz=PyFloat_AS_DOUBLE(O_cosz);

  par=dvector(0,Npar-1);
  for (i=0;i<Npar;i++)
    {
      par[i]=*(double *)(temp->data + i*temp->strides[0]);
      //if (i<3) printf("%lg ",par[i]);
    }
  
  h=dvector(0,100);
  elementsitosiavo_(&cosx,&cosy,&cosz,&par[0],&par[1],&par[2],&par[3],&par[4],&par[5],&par[6],&par[7],&par[8],&par[9],
		    &par[10],&par[11],&par[12],&par[13],h);

  obj=(PyObject *)PyArray_FromDimsAndData(1,&Nh,PyArray_DOUBLE,(char*)h);
  return Py_BuildValue("O", obj);
}  
