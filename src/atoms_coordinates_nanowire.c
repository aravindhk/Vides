// ======================================================================
//  Copyright (c) 2011, G. Fiori, P. Marconcini, P. D'Amico, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "atoms_coordinates_nanowire.h"

static PyObject* py_atoms_coordinates_nanowire(PyObject* self, PyObject* args)
{

  int i, j, rank,Npar,Nh,Np,temp_int1,temp_int2;
  double par[5],dim[1],*h;

  PyObject *obj,*obj2,*obj3;
  PyArrayObject *temp, *out;

  int n_aux, Nc_aux, max, slices;
  int sqci, tilt;
  double a0,  edge, zmax,*ics, *ipsilon, *zeta;


  //  Py_complex cc;
  //  complex ccc;
  
  setvbuf(stdout,(char *) NULL,_IONBF,1);
  Npar=5;
  n_aux=100;
  Nc_aux=1000;
  
/*
  Reads the data from python class.
*/
  
  import_array();
  if (!PyArg_ParseTuple(args,"O",&temp))
    {
      //printf"NOT A CLASS! \n");
      exit(0);
    }
  
  for (i=0;i<Npar;i++)
    {
      par[i]=*(double *)(temp->data + i*temp->strides[0]);
      //      printf("%lg ",par[i]);
    }

  // par[0]=a0 in Angstrom
  // par[1]=sqci
  // par[2]=tilt
  // par[3]=edge
  // par[4]=zmax

  ics = dvector(0,n_aux*Nc_aux-1);
  ipsilon = dvector(0,n_aux*Nc_aux-1);
  zeta = dvector(0,n_aux*Nc_aux-1);
  for (i=0;i<n_aux*Nc_aux;i++)
    {
      ics[i]=-10;
      ipsilon[i]=-10;
      zeta[i]=-10;
    }

  
  temp_int1=par[1];
  temp_int2=par[2];
  nanowire_(&par[0], &temp_int1, &temp_int2, &par[3], &par[4], &slices, ics, ipsilon, zeta, &n_aux, &Nc_aux );

  Np=0;
  for (i=0;i<n_aux*Nc_aux;i++)
    {
      if (ics[i]>-10) Np++;
      zeta[i]-=par[0];
    }
  
  obj=(PyObject *)PyArray_FromDimsAndData(1,&Np,PyArray_DOUBLE,(char*)ics);
  obj2=(PyObject *)PyArray_FromDimsAndData(1,&Np,PyArray_DOUBLE,(char*)ipsilon);
  obj3=(PyObject *)PyArray_FromDimsAndData(1,&Np,PyArray_DOUBLE,(char*)zeta);
  return Py_BuildValue("OOO", obj,obj2,obj3);
  
  char *s = "Hello from CNT_charge_T";
  return Py_BuildValue("s", s);
  //  return 0;
}  
