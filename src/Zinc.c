// ======================================================================
//  Copyright (c) 2011, P. D'Amico, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "Zinc.h"

static PyObject* py_Zinc(PyObject* self, PyObject* args)
{

  int i, j, flag, rank, atoms, slices, n_aux, Nc_aux, max;
  int sqci, tilt;
  double a0,  edge, zmax, deltae, *ics, *ipsilon, *zeta, *H_aux, *skparameters;

  PyObject *obj, *temp_obj;
  PyArrayObject *temp_out, *zeta_out;


  //  Py_complex cc;
  //  complex ccc;
  
  setvbuf(stdout,(char *) NULL,_IONBF,1);
  
  
/*
  Reads the data from python class.
*/
  
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      //printf"NOT A CLASS! \n");
      exit(0);
    }
  
  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);
  printf("Rank = %d \n",rank);

  
  if (!rank) printf("*****************************************\n");
  if (!rank) printf("*** Computing NEGF in HAMILTONIAN *******\n"); 
  if (!rank) printf("*****************************************\n\n");

  temp_obj=PyObject_GetAttrString(obj,"sqci");
  sqci=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Section = %d \n",sqci);

  temp_obj=PyObject_GetAttrString(obj,"tilt");
  tilt=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Surface orientation = %d \n",tilt);

  temp_obj=PyObject_GetAttrString(obj,"edge");
  edge=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Diameter = %lg \n", edge);

  temp_obj=PyObject_GetAttrString(obj,"zmax");
  zmax=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Lenght of the wire = %lg \n", zmax);

  temp_obj=PyObject_GetAttrString(obj,"deltae");
  deltae=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Delta energy = %lg \n", deltae);

  temp_obj=PyObject_GetAttrString(obj,"n_aux");
  n_aux=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("n_aux = %d \n", n_aux);

  temp_obj=PyObject_GetAttrString(obj,"Nc_aux");
  Nc_aux = (int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Nc_aux = %d \n", Nc_aux);

  temp_obj=PyObject_GetAttrString(obj,"a0");
  a0 = (double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("a0 = %lg \n", a0);

  //if (flag == 0 ) a0 = 5.431;
  //  if (flag == 1 ) a0 =  5.6575;
 
  ics = dvector(0,n_aux*Nc_aux-1);
  ipsilon = dvector(0,n_aux*Nc_aux-1);
  zeta = dvector(0,n_aux*Nc_aux-1);
  H_aux = dvector(0,(Nc_aux*n_aux)*((Nc_aux*n_aux+1)/2)*(2+100) - 1);

  /*
  temp_out = (PyArrayObject*)PyObject_GetAttrString(obj,"skparameters");
  
  skparameters = dvector(0, 17);
  
  double *pass;
  int  strides;
  strides=sizeof(double);
  
  for(i=0; i<18; i++ ){
    pass =  (double *)(temp_out->data + i*strides);
    skparameters[i] = (double)(*pass);
  }
  
  for (i=0; i< 18; i++)
    printf("%lg ", skparameters[i]);
  
  printf("\n");
*/

  
  temp_obj=PyObject_GetAttrString(obj,"flag");
  flag = (int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Nc_aux = %d \n", flag);

  if (flag == 0 ){

    temp_out = (PyArrayObject*)PyObject_GetAttrString(obj,"skparameters");
    skparameters = dvector(0, 17);
    double *pass;
    int  strides;
    strides=sizeof(double);
    
    for(i=0; i<18; i++ ){
      pass =  (double *)(temp_out->data + i*strides);
      skparameters[i] = (double)(*pass);
    }
    
    for (i=0; i< 18; i++)
      printf("%lg ", skparameters[i]);
  
    printf("\n");

    nanowire_(&a0, &sqci, &tilt, &edge, &zmax, &slices, ics, ipsilon, zeta, &n_aux, &Nc_aux );
  
    hamiltonian_( &deltae, &atoms, H_aux, &n_aux, &Nc_aux, &max, skparameters );
  }

  if (flag == 1 ){

    temp_out = (PyArrayObject*)PyObject_GetAttrString(obj,"skparameters");
    skparameters = dvector(0, 28);
    double *pass;
    int  strides;
    strides=sizeof(double);
    
    for(i=0; i<29; i++ ){
      pass =  (double *)(temp_out->data + i*strides);
      skparameters[i] = (double)(*pass);
    }
    
    for (i=0; i< 29; i++)
      printf("%lg ", skparameters[i]);
    
    printf("\n");
    
    nanowire_( &a0, &sqci, &tilt, &edge, &zmax, &slices, ics, ipsilon, zeta, &n_aux, &Nc_aux );
  
    hamil_( &deltae, &atoms, H_aux, &n_aux, &Nc_aux, &max, skparameters );
  }

  //  for (j=0; j<102; j++)
  // printf("%lg, ", H_aux[j + 3*102]);
  //printf("\n");
  
  
  //exit(0);

  printf("n = %d \n", atoms);

  printf("Nc = %d \n", slices);

  printf("max = %d \n", max);

  //  exit(0);

  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"ics");
  for (i=0;i<n_aux*Nc_aux;i++)
    *(double *)(temp_out->data + i*temp_out->strides[0])=ics[i];
  free_dvector(ics, 0, n_aux*Nc_aux-1);

  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"ipsilon");
  for (i=0;i<n_aux*Nc_aux;i++)
    *(double *)(temp_out->data + i*temp_out->strides[0])=ipsilon[i];
  free_dvector(ipsilon, 0, n_aux*Nc_aux-1);

  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"zeta");
  for (i=0;i<n_aux*Nc_aux;i++)
    *(double *)(temp_out->data + i*temp_out->strides[0])=zeta[i];
  free_dvector(zeta, 0, n_aux*Nc_aux-1);

  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"H_aux");
  for (i=0;i<(Nc_aux*n_aux)*((Nc_aux*n_aux+1)/2)*(2+100);i++)
    *(double *)(temp_out->data + i*temp_out->strides[0])=H_aux[i];
  free_dvector(H_aux, 0, (Nc_aux*n_aux)*((Nc_aux*n_aux+1)/2)*(2+100) - 1);

  /*
  for (i=0; i<160; i++)
    printf("%lg, ", zeta[i]);

  printf("\n");
  */
  //  exit(0);


  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"atoms");
  *(double *)(temp_out->data) = atoms;


  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"slices");
  *(double *)(temp_out->data) = slices;

  temp_out=(PyArrayObject*)PyObject_GetAttrString(obj,"max");
  *(double *)(temp_out->data) = max;
  
  char *s = "Hello from CNT_charge_T";
  return Py_BuildValue("s", s);
  //  return 0;
}  
