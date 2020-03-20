#include "Python.h"
#include "arrayobject.h"
#include <stdio.h>
#include "nrutil.h"
#include "nonuniformgridmod.h"

#include "CNT_charge_T.c"
#include "CNTmode_charge_T.c"
#include "GNR_charge_T.c"
#include "CNT_atoms_coordinates.c"
#include "solvePoisson.c"
#include "solvePoisson2D.c"
#include "GNR_atoms_coordinates.c"
#include "GNRgap.c"
#include "H_charge_T.c"
#include "electron_GNR.c"
#include "phonon_GNR.c"
#include "phonon_graphene.c"
#include "rateACABS.c"
#include "rateACEM.c"
#include "rateOPTABS.c"
#include "rateOPTEM.c"
#include "Zinc.c"
#include "Simatrix.c"
#include "atoms_coordinates_nanowire.c"
#include "Sipassivation.c"
#include "writeout.c"
#include "solve_schroedinger_1D.c"
#include "solvePoisson1D.c"
#include "fphalf.c"
#include "Fermi_Integrals.c"

#if PY_MAJOR_VERSION >=3
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#endif


static PyObject* py_nonuniformgridmod(PyObject* self, PyObject* args)
{
  int n,nx,i,Np;
  PyArrayObject *meshdelta;
  PyArrayObject *b;
  double *outdata;
  double *mesh,*delta,*grid;
  int length,numregioni;
  // Import array e' importante altrimenti mi da
  // segmentation fault PyArray_FromDimsAndData
  
  import_array();

  if (!PyArg_ParseTuple(args,"O",&meshdelta))
    {
      printf("NOT AN ARRAY! \n");
      exit(0);
    }
  double e;
  e=3;
  nx=meshdelta->dimensions[0];
  Np=nx/2;
  mesh=dvector(0,Np-1);
  delta=dvector(0,Np-1);
  grid=dvector(0,10000);
  for (i=0;i<Np;i++)
    {
      outdata=(double *)(meshdelta->data+2*i*meshdelta->strides[0]);
      mesh[i]=*outdata;
      outdata=(double *)(meshdelta->data+(2*i+1)*meshdelta->strides[0]);
      delta[i]=*outdata;
      //      printf("%lg ",delta[i]);
    }
  
  //  printf("The dimension of the array is %d \n",nx);
  numregioni=Np-1;
  nonuniformgridmod_(grid,&Np,mesh,delta,&numregioni);

  b = (PyArrayObject *) PyArray_FromDimsAndData(1,&Np,PyArray_DOUBLE,(char*)grid);
  return  PyArray_Return(b);

  //  printf("Eccoci \n");
  //  char *s = "Hello from C!";
  //  return Py_BuildValue("s", s);
}

static PyObject* py_testing(PyObject* self, PyObject* args)
{
  int N,nx,i;
  double *Phi,*tempo;
  PyObject *obj,*temp_obj;
  PyArrayObject *Phi_array;
  
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  temp_obj=PyObject_GetAttrString(obj,"N");
  N=(int)PyInt_AsLong(temp_obj);
  printf("N is equal to %d \n",N);

  N=N+10;

  Phi_array=(PyArrayObject *)PyObject_GetAttrString(obj,"Phi");
  nx=Phi_array->dimensions[0];
  Phi=dvector(0,nx-1);
  printf("Length og the array = %d \n",nx);
  for (i=0;i<nx;i++)
    {
      tempo=(double *)(Phi_array->data+i*Phi_array->strides[0]);
      Phi[i]=(*tempo)*3;
    }

  // I set the attribute of the class to send data out
  temp_obj=  PyArray_FromDimsAndData(1,&nx,PyArray_DOUBLE,(char*)Phi);
  PyObject_SetAttrString(obj,"Phi",temp_obj);

  temp_obj=PyInt_AsLong(N);
  PyObject_SetAttrString(obj,"N",temp_obj);

  char *s = "Hello from testing!";
  return Py_BuildValue("s", s);
}



static PyMethodDef NanoTCAD_ViDESmod_methods[] = {
  {"nonuniformgridmod", py_nonuniformgridmod, METH_VARARGS},
  {"testing", py_testing, METH_VARARGS},
  {"CNT_charge_T", py_CNT_charge_T, METH_VARARGS},
  {"CNTmode_charge_T", py_CNTmode_charge_T, METH_VARARGS},
  {"GNR_charge_T", py_GNR_charge_T, METH_VARARGS},
  {"CNT_atoms_coordinates", py_CNT_atoms_coordinates, METH_VARARGS},
  {"solvePoisson", py_solvePoisson, METH_VARARGS},
  {"solvePoisson2D", py_solvePoisson2D, METH_VARARGS},
  {"GNR_atoms_coordinates", py_GNR_atoms_coordinates, METH_VARARGS},
  {"GNRgap", py_GNRgap, METH_VARARGS},
  {"H_charge_T", py_H_charge_T, METH_VARARGS},
  {"electron_GNR", py_electron_GNR, METH_VARARGS},
  {"phonon_GNR", py_phonon_GNR, METH_VARARGS},
  {"phonon_graphene", py_phonon_graphene, METH_VARARGS},
  {"rateACABS", py_rateACABS, METH_VARARGS},
  {"rateACEM", py_rateACEM, METH_VARARGS},
  {"rateOPTABS", py_rateOPTABS, METH_VARARGS},
  {"rateOPTEM", py_rateOPTEM, METH_VARARGS},
  {"Zinc", py_Zinc, METH_VARARGS},
  {"Simatrix", py_Simatrix, METH_VARARGS},
  {"atoms_coordinates_nanowire", py_atoms_coordinates_nanowire, METH_VARARGS},
  {"Sipassivation", py_Sipassivation, METH_VARARGS},
  {"writeout", py_writeout, METH_VARARGS},
  {"solve_schroedinger_1D", py_solve_schroedinger_1D, METH_VARARGS},
  {"solvePoisson1D", py_solvePoisson1D, METH_VARARGS},
  {"fphalf", py_fphalf, METH_VARARGS},
  {"Fermi_Integrals", py_Fermi_Integrals, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION <3
/*
 * Bind Python function names to our C functions
 */

/*
 * Python calls this to let us initialize our module
 */
void initNanoTCAD_ViDESmod()
{
	(void) Py_InitModule("NanoTCAD_ViDESmod", NanoTCAD_ViDESmod_methods);
}

#else

static struct PyModuleDef NanoTCAD_ViDESmod = {
  PyModuleDef_HEAD_INIT,
  "NanoTCAD_ViDESmod",
  "This is a module",
  -1,
  NanoTCAD_ViDESmod_methods,
  NULL,
  NULL
};

PyMODINIT_FUNC 
PyInit_NanoTCAD_ViDESmod(void)
{
  return PyModule_Create(&NanoTCAD_ViDESmod);
}

#endif
