// ======================================================================
//  Copyright (c) 2012, G. Fiori  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "solve_schroedinger_1D.h"
static PyObject* py_solve_schroedinger_1D(PyObject* self, PyObject* args)
{

  time_t t1, t2;
  int i, j, ix, NE, nx, ny, Narray, rank, Neig,k;
  PyObject *obj,*temp_obj;
  PyArrayObject *Ei_array,*Psi_array,*eig_array,*mass_array,*x_array,*Egap_array;
  double **Psi,**eig,**mass,*x,*m_temp,*eig_temp,*eigfun_temp,*Egap,*ecs;
  double Temperature,vt,*Ei,*tempo;
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
  if (!rank) printf("***** Computing Schroedinger 1D *********\n"); 
  if (!rank) printf("*****************************************\n\n");


  temp_obj=PyObject_GetAttrString(obj,"nx");
  nx=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("n # points in the vertical direction = %d \n",nx);
  temp_obj=PyObject_GetAttrString(obj,"ny");
  ny=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Number of slices along the channel = %d \n",ny);
  temp_obj=PyObject_GetAttrString(obj,"Temp");
  Temperature=(double)PyFloat_AsDouble(temp_obj);
  vt=kboltz*Temperature/q;
  if (!rank) printf("Vt = %lg \n",vt);
  if (!rank) printf("Temperature = %lg \n",Temperature);
  temp_obj=PyObject_GetAttrString(obj,"Neig");
  Neig=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("# considered eigenvalues = %d \n",Neig);

/*
  USE THE FOLLOWING TWO LINES FOR CELL REPRESENTATION, OTHERWISE COMMENT IT!!
*/

/*
  ny = (ny)/4;
  n = n*4;
*/

  Ei=dvector(0,nx*ny-1);
  Ei_array=(PyArrayObject*)PyObject_GetAttrString(obj,"Ei");
  Narray=Ei_array->dimensions[0];
  for (i=0;i<Narray;i++){
    tempo=(double *)(Ei_array->data+i*Ei_array->strides[0]);
    Ei[i]=*tempo;
  }

/*   printf("Ei \n"); */
/*   for (i=0;i<Narray;i++) */
/*     printf("%lg ",Ei[i]); */
/*   printf("\n"); */

  Egap=dvector(0,nx*ny-1);
  Egap_array=(PyArrayObject*)PyObject_GetAttrString(obj,"Egap");
  Narray=Egap_array->dimensions[0];
  for (i=0;i<Narray;i++){
    tempo=(double *)(Egap_array->data+i*Egap_array->strides[0]);
    Egap[i]=*tempo;
  }

/*   printf("Egap \n"); */
/*   for (i=0;i<Narray;i++) */
/*     printf("%lg ",Egap[i]); */
/*   printf("\n"); */

  // I just get the first n points of the x vector
  x=dvector(0,nx-1);
  x_array=(PyArrayObject*)PyObject_GetAttrString(obj,"x");
  for (i=0;i<nx;i++){
    tempo=(double *)(x_array->data+i*x_array->strides[0]);
    x[i]=*tempo;
  }
  
  Psi=dmatrix(0,nx*ny-1,0,Neig-1);
  eig=dmatrix(0,ny-1,0,Neig-1);  
  mass=dmatrix(0,nx-1,0,ny-1);
  
  mass_array=(PyArrayObject*)PyObject_GetAttrString(obj,"mass");
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      {
	tempo=(double *)(mass_array->data+i*mass_array->strides[0]+j*mass_array->strides[1]);
	mass[i][j]=*tempo;
      }

/*   printf("mass \n"); */
/*   for (i=0;i<nx;i++) */
/*     { */
/*       for (j=0;j<ny;j++) */
/* 	printf("%lg ",mass[i][j]); */
/*       printf("\n"); */
/*     } */
/*   printf("\n"); */
/*   exit(0); */


  // HERE THE COMPUTATION OF THE SCHROEDINGER EQUATION STARTS
  
  eig_temp=dvector(0,Neig-1);
  eigfun_temp=dvector(0,Neig*nx-1);
  m_temp=dvector(0,nx-1);
  ecs=dvector(0,nx-1);

  for (j=0;j<ny;j++)
    {
      for (i=0;i<nx;i++)
	{
	  m_temp[i]=mass[i][j];
	  ecs[i]=Ei[i+j*nx]+Egap[i+j*nx]*0.5;
	}

      //      printf("Eccomi qua \n");
	
      
      ris1d_(ecs,x,m_temp,&nx,eig_temp,eigfun_temp,&Neig);

      for (k=0;k<Neig;k++)
	{
	  for (i=0;i<nx;i++)
	    Psi[i+j*nx][k]=eigfun_temp[i+k*nx];
	  eig[j][k]=eig_temp[k];
	}

/*       FILE *fp; */
/*       fp=fopen("eigve2.out","w"); */
/*       for (i=0;i<nx;i++) */
/* 	fprintf(fp,"%lg %lg \n",x[i],Psi[i][0]); */
/*       fclose(fp); */
/*       exit(0); */
    
    }

    
/*
  give back the data
*/

  Psi_array=(PyArrayObject*)PyObject_GetAttrString(obj,"Psi");
  for (i=0;i<nx*ny;i++)
    for (j=0;j<Neig;j++)
      *(double *)(Psi_array->data+i*Psi_array->strides[0]+j*Psi_array->strides[1])=Psi[i][j];
  free_dmatrix(Psi,0,nx*ny-1,0,Neig-1);

  eig_array=(PyArrayObject*)PyObject_GetAttrString(obj,"eig");
  for (i=0;i<ny;i++)
    for (j=0;j<Neig;j++)
    *(double *)(eig_array->data+i*eig_array->strides[0]+j*eig_array->strides[1])=eig[i][j];
  free_dmatrix(eig,0,ny-1,0,Neig-1);

  /*  
      Free the memory.
  */

  free_dvector(eig_temp,0,Neig-1);
  free_dvector(eigfun_temp,0,Neig*nx-1);
  free_dvector(m_temp,0,nx-1);
  free_dvector(ecs,0,nx-1);
  free_dvector(Ei,0,nx*ny-1);
  free_dvector(Egap,0,nx*ny-1);
  free_dvector(x,0,nx-1);

  Py_DECREF(Ei_array);
  Py_DECREF(Egap_array);
  Py_DECREF(Psi_array);
  Py_DECREF(eig_array);
  Py_DECREF(mass_array);
  Py_DECREF(x_array);

  char *s = "Hello from CNT_charge_T";
  return Py_BuildValue("s", s);
 
}  
