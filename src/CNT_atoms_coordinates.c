// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

//
//
// This subroutine computes the atoms coordinates in the 3D space
// of the carbon nanotubes. 
//
// Input variables 
//
// nanotube.acc [nm]: carbon-to-carbon distance
// nanotube.n: number of atoms along CNT ring
// nanotube.Nc: number of CNT rings
//
// Output variables
//
// N.B.: Atoms are ordered first within the ring and then along the 
// CNT axis. E.G.: if i runs from 0 to n-1 and j from 0 to Nc-1
// the ix-th atom with index (i,j) has index ix=i+j*n
//
// nanotube.x (vector n*Nc) [nm]: x coordinates of C atoms ordered as above
// nanotube.y (vector n*Nc) [nm]: y coordinates of C atoms ordered as above
// nanotube.z (vector n*Nc) [nm]: z coordinates of C atoms ordered as above
//
//

#include "CNT_atoms_coordinates.h"
static PyObject* py_CNT_atoms_coordinates(PyObject* self, PyObject* args)
{
  int i,j,count,Nc,n,ie;
  PyObject *obj,*temp_obj;
  PyArrayObject *Phi_array;
  double xt,yt,zt,*x,*y,*z,acc;

  // I read the data from python class
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }
  temp_obj=PyObject_GetAttrString(obj,"n");
  n=(int)PyInt_AsLong(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Nc");
  Nc=(int)PyInt_AsLong(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"acc");
  acc=(double)PyFloat_AsDouble(temp_obj);

  x=dvector(0,n*Nc-1);
  y=dvector(0,n*Nc-1);
  z=dvector(0,n*Nc-1);

  count=0;
  for (j=0;j<Nc;j++)
    {
      for (i=0;i<n;i++)
	{
	  nanotube(i,j,&xt,&yt,&zt,n,acc);
	  x[count]=xt;
	  y[count]=yt;
	  z[count]=zt;
	  count++;
	}
    }

  ie=n*Nc;
  // I set the attribute of the class to send data out
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)x);
  PyObject_SetAttrString(obj,"x",temp_obj);
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)y);
  PyObject_SetAttrString(obj,"y",temp_obj);
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)z);
  PyObject_SetAttrString(obj,"z",temp_obj);
  
  char *s = "Hello from CNT_atoms_coordinates";
  return Py_BuildValue("s", s);
}
