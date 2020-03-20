// ======================================================================
//  Copyright (c) 2011, G. Fiori, P. Marconcini, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "Sipassivation.h"
static PyObject* py_Sipassivation(PyObject* self, PyObject* args)
{

  int i, j, rank,Npar,Nh,temp_int1,temp_int2,temp_int3,temp_int4,numfn;
  double *par,dim[1],*h,x,y,z,*xfn,*yfn,*zfn,deltae;

  PyObject *obj;
  PyArrayObject *temp,*O_x,*O_y,*O_z,*O_numfn,*O_xfn,*O_yfn,*O_zfn,*O_deltae, *out;

  //  Py_complex cc;
  //  complex ccc;

  Nh=16;
  
  setvbuf(stdout,(char *) NULL,_IONBF,1);
  
  
/*
  Reads the data from python class.
*/
  
  import_array();
  if (!PyArg_ParseTuple(args,"OOOOOOOO",&O_x,&O_y,&O_z,&O_numfn,&O_xfn,&O_yfn,&O_zfn,&O_deltae))
    {
      //printf"NOT A CLASS! \n");
      exit(0);
    }
  
  x=PyFloat_AS_DOUBLE(O_x);
  y=PyFloat_AS_DOUBLE(O_y);
  z=PyFloat_AS_DOUBLE(O_z);
  deltae=PyFloat_AS_DOUBLE(O_deltae);
  #if PY_MAJOR_VERSION >=3
  numfn=PyLong_AS_LONG((PyObject*)O_numfn);
  #else
  numfn=PyInt_AS_LONG(O_numfn);
  #endif
  //  printf("%lg %lg %lg %lg %d \n",x,y,z,deltae,numfn);

  xfn=dvector(0,3);
  yfn=dvector(0,3);
  zfn=dvector(0,3);
  for (i=0;i<4;i++)
    {
      xfn[i]=*(double *)(O_xfn->data + i*O_xfn->strides[0]);
      yfn[i]=*(double *)(O_yfn->data + i*O_yfn->strides[0]);
      zfn[i]=*(double *)(O_zfn->data + i*O_zfn->strides[0]);
      //      printf("%lg %lg %lg \n",xfn[i],yfn[i],zfn[i]);
    }
  
  h=dvector(0,16);
  
  passivationsiavo_(&x,&y,&z,&numfn,xfn,yfn,zfn,&deltae,h);

  //  Py_DECREF(O_xfn);
  //  Py_DECREF(O_yfn);
  //  Py_DECREF(O_zfn);
  
  free_dvector(xfn,0,3);
  free_dvector(yfn,0,3);
  free_dvector(zfn,0,3);

  obj=(PyObject *)PyArray_FromDimsAndData(1,&Nh,PyArray_DOUBLE,(char*)h);
  return Py_BuildValue("O", obj);
}  
