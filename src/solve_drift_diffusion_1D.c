// ======================================================================
//  Copyright (c) 2012, G. Fiori  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "solve_drift_diffusion_1D.h"
static PyObject* py_solve_drift_diffusion_1D(PyObject* self, PyObject* args)
{

  time_t t1, t2;
  int i, j, ix, NE, nx, ny, Narray, rank,k;
  PyObject *obj,*temp_obj;
  PyArrayObject *mu_array,*ecs_array,*n1d_array,*y_array,*genric_array;
  double *y,*genric,*ecs;
  double Temperature,vt,nlims,nlimd,*mu,*tempo,tolljay,*n1d;
  int strides;
  
  setvbuf(stdout,(char *) NULL,_IONBF,1);
  
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      //printf"NOT A CLASS! \n");
      exit(0);
    }
  
  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);
  
  if (!rank) printf("*****************************************\n");
  if (!rank) printf("***** Computing Drift Diffusion 1D ******\n"); 
  if (!rank) printf("*****************************************\n\n");


  temp_obj=PyObject_GetAttrString(obj,"ny");
  ny=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("n # points in the longitudinal direction = %d \n",ny);
  temp_obj=PyObject_GetAttrString(obj,"tolljay");
  tolljay=(double)PyInt_AsLong(temp_obj);
  if (!rank) printf("Tollerance = %lg \n",tolljay);
  temp_obj=PyObject_GetAttrString(obj,"Temp");
  Temperature=(double)PyFloat_AsDouble(temp_obj);
  vt=kboltz*Temperature/q;
  if (!rank) printf("Vt = %lg \n",vt);
  if (!rank) printf("Temperature = %lg \n",Temperature);
  temp_obj=PyObject_GetAttrString(obj,"charge_left_contact");
  nlims=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Charge at the left contact = %lg \n",nlims);
  temp_obj=PyObject_GetAttrString(obj,"charge_right_contact");
  nlimd=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Charge at the right contact = %lg \n",nlimd);

/*
  USE THE FOLLOWING TWO LINES FOR CELL REPRESENTATION, OTHERWISE COMMENT IT!!
*/

/*
  ny = (ny)/4;
  n = n*4;
*/

  mu=dvector(0,ny-1);
  mu_array=(PyArrayObject*)PyObject_GetAttrString(obj,"mu");
  Narray=mu_array->dimensions[0];
  for (i=0;i<Narray;i++){
    tempo=(double *)(mu_array->data+i*mu_array->strides[0]);
    mu[i]=*tempo;
  }

  printf("mu \n");
  for (i=0;i<Narray;i++)
    printf("%lg ",mu[i]);
  printf("\n");

  genric=dvector(0,ny-1);
  genric_array=(PyArrayObject*)PyObject_GetAttrString(obj,"genric");
  Narray=genric_array->dimensions[0];
  for (i=0;i<Narray;i++){
    tempo=(double *)(genric_array->data+i*genric_array->strides[0]);
    genric[i]=*tempo;
  }

  printf("genric \n");
  for (i=0;i<Narray;i++)
    printf("%lg ",genric[i]);
  printf("\n");

  // I just get the first n points of the y vector
  y=dvector(0,ny-1);
  y_array=(PyArrayObject*)PyObject_GetAttrString(obj,"y");
  for (i=0;i<ny;i++){
    tempo=(double *)(y_array->data+i*y_array->strides[0]);
    y[i]=*tempo;
  }

  printf("y \n");
  for (i=0;i<Narray;i++)
    printf("%lg ",y[i]);
  printf("\n");


// I take the band vector
  ecs=dvector(0,ny-1);
  ecs_array=(PyArrayObject*)PyObject_GetAttrString(obj,"ecs");
  for (i=0;i<ny;i++){
    tempo=(double *)(ecs_array->data+i*ecs_array->strides[0]);
    ecs[i]=*tempo;
  }

  printf("ecs \n");
  for (i=0;i<Narray;i++)
    printf("%lg ",ecs[i]);
  printf("\n");


  // HERE THE COMPUTATION OF THE DRIFT DIFFUTION EQUATION STARTS
  
  n1d=dvector(0,ny-1);
  for (i=0;j<ny;j++)
    n1d[i]=0;

  sub_jay1d_(n1d,ecs,y,&ny,mu,genric,&tolljay,&vt,&nlims,&nlimd);

  n1d_array=(PyArrayObject*)PyObject_GetAttrString(obj,"n1d");
  for (i=0;i<ny;i++)
    *(double *)(n1d_array->data+i*n1d_array->strides[0])=n1d[i];
  free_dvector(n1d,0,ny-1);

  free_dvector(mu,0,ny-1);
  free_dvector(ecs,0,ny-1);
  free_dvector(genric,0,ny-1);
  free_dvector(y,0,ny-1);

  Py_DECREF(mu_array);
  Py_DECREF(genric_array);
  Py_DECREF(ecs_array);
  Py_DECREF(y_array);

  char *s = "Hello from CNT_charge_T";
  return Py_BuildValue("s", s);
 
}  
