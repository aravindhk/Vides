// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

/*  ****************************************************************************************  */
/*                 This subroutine computes the Poisson equation.                             */
/* 		p and d are the structures in which the physical                              */
/* 		quantities and the mapping of the structure are                               */
/* 		specified. IERR is a flag which                                               */
/* 		specify if the system converge (IERR=0) or                                    */
/* 		not (IERR=2). In this case the iteration stops.                               */
/* 		tolldomn is the tollerance for the linear system                              */
/* 		solver (DOMN). neq is the number of equation to be                            */
/* 		solved. In case you solve only the poisson equation                           */
/* 		eq=1. mat is the structure in which all the quantities                        */
/* 		relative to the material are specified. normapoisson is                       */
/* 		the tollerance for the poisson inner iteration loop.                          */
/* 		J is the 5 dimensional tensor in which the jacobian is                        */
/* 		specified. sottoril is the under/over relaxation coefficient.                 */
/* 		nf, pf, accf, and donf are the user specified array of                        */
/* 		point to function relative to the electron, hole                              */
/* 		ionized acceptors and ionized donors concentrations, respectively.            */
/* 		dnf, dpf, daccf, and ddonf are the user specified array of                    */
/* 		point to function relative to the derivative with respect                     */
/* 		to the electrostatic potential (eV) of the electron, hole                     */
/* 		ionized acceptors and ionized donors concentration, respectively.             */
/*  ****************************************************************************************  */		  


#include "solvePoisson.h"


static PyObject* py_solvePoisson(PyObject* self, PyObject* args)
{
  int i,j,k,ix,ciclom,ciclomm,ciclo,iterazione,Np,gumgum,neq,*temporaryint,rank,*IERR;
  double normab,normabold,normabb,normavp,*temporary,sottoril;
  void (*comp)(int *,double *,double *,int *,int *,int *,double *,int *);
  void (*comp2)(int *,double *,double *, int *,int *,int *,double *,int *,double *,int *);
  int NELT, *IA, *JA, NSAVE,*IWORK;
  double *B, *XLIN, *A, *R, *Z2,*****J; 
  double *DZ,*CSAV,*RWORK;
  //  double P2[50000][10],AP[50000][10],EMAP[50000][10];
  double *P2,*AP,*EMAP;
  int ISYM, ITOL, ITMAX;
  int  ITER, IUNIT,Z1;
  double TOL, ERR; 

  PyObject *grid,*interface,*temp_obj,*testo;
  int nx,ny,nz;
  PyArrayObject *eps_array,*boundary_conditions_array,
    *dVe_array,*surf_array,*dist_array,*Phi_array,*free_charge_array,*fixed_charge_array;
  device_mapping d;
  physical_quantities p;
  double tolldomn,normapoisson,vt;

  setvbuf(stdout,(char *) NULL,_IONBF,1);
  
  //************************************************************************
  // I read the data from python class
  //************************************************************************

  import_array();
  if (!PyArg_ParseTuple(args,"OO",&grid,&interface))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  temp_obj=PyObject_GetAttrString(interface,"rank");
  rank=(int)PyInt_AsLong(temp_obj);
  // I read the number of point along the three axis
  temp_obj=PyObject_GetAttrString(grid,"nx");
  d.nx=(int)PyInt_AsLong(temp_obj);
  temp_obj=PyObject_GetAttrString(grid,"ny");
  d.ny=(int)PyInt_AsLong(temp_obj);
  temp_obj=PyObject_GetAttrString(grid,"nz");
  d.nz=(int)PyInt_AsLong(temp_obj);

  // I set the number of total points in the 3D domain
  Np=d.nx*d.ny*d.nz;

  // I read the tollerance of the domn linear solver
  temp_obj=PyObject_GetAttrString(interface,"tolldomn");
  tolldomn=(double)PyFloat_AsDouble(temp_obj);

  // I read the norm of the variation of the potential for 
  // which the Poisson solver exits the loop
  temp_obj=PyObject_GetAttrString(interface,"normpoisson");
  normapoisson=(double)PyFloat_AsDouble(temp_obj);

  // I read the under relaxation coefficient
  temp_obj=PyObject_GetAttrString(interface,"underel");
  sottoril=(double)PyFloat_AsDouble(temp_obj);

  // I read the electrostatic potential along the 3D domain
  p.Phi=dvector(0,Np-1);
  Phi_array=(PyArrayObject *)PyObject_GetAttrString(interface,"Phi");

  for (i=0;i<Np;i++)
    {
      temporary=(double *)(Phi_array->data+i*Phi_array->strides[0]);
      p.Phi[i]=*temporary;
    }

  // I read the dielectric constant along the 3D domain
  p.eps=dvector(0,Np-1);
  eps_array=(PyArrayObject *)PyObject_GetAttrString(interface,"eps");
  for (i=0;i<Np;i++)
    {
      temporary=(double *)(eps_array->data+i*eps_array->strides[0]);
      p.eps[i]=*temporary;
    }

  // I read the volumes of the elements
  d.dVe=dvector(0,Np-1);
  dVe_array=(PyArrayObject *)PyObject_GetAttrString(grid,"dVe");
  for (i=0;i<Np;i++)
    {
      temporary=(double *)(dVe_array->data+i*dVe_array->strides[0]);
      d.dVe[i]=*temporary*1e-27;
    }

  // I read the surfaces of the elements
  d.surf=dmatrix(0,Np-1,0,5);
  surf_array=(PyArrayObject *)PyObject_GetAttrString(grid,"surf");
  for (i=0;i<Np;i++)
    for (j=0;j<6;j++)
      {
	//	printf("%d %d \n",i*surf_array->strides[0],j*8);
	temporary=(double *)(surf_array->data + i*surf_array->strides[0] + j*sizeof(double));
	d.surf[i][j]=*temporary*1e-18;;
      }

  //  printf("%d %d %d \n",surf_array->strides[0],surf_array->strides[0]/8,surf_array->nd);
  //  printf("%d %d \n",surf_array->dimensions[0],surf_array->dimensions[1]);

  // I read the distances between the first neighbours

  d.dist=dmatrix(0,Np-1,0,5);
  dist_array=(PyArrayObject *)PyObject_GetAttrString(grid,"dist");
  for (i=0;i<Np;i++)
    for (j=0;j<6;j++)
      {
	temporary=(double *)(dist_array->data + i*dist_array->strides[0] + j*sizeof(double));
	d.dist[i][j]=*temporary*1e-9;
      }

  // I read the boundary conditions along the 3D domain
  d.boundary_conditions=dvector(0,Np-1);
  boundary_conditions_array=(PyArrayObject *)PyObject_GetAttrString(interface,"boundary_conditions");
  for (i=0;i<Np;i++)
    {
      temporary=(double *)(boundary_conditions_array->data+i*boundary_conditions_array->strides[0]);
      d.boundary_conditions[i]=*temporary;
    }

  // I read the free charge term of the Poisson equation
  p.free_charge=dvector(0,Np-1);
  free_charge_array=(PyArrayObject *)PyObject_GetAttrString(interface,"free_charge");
  for (i=0;i<Np;i++)
    {
      temporary=(double *)(free_charge_array->data+i*free_charge_array->strides[0]);
      p.free_charge[i]=*temporary;
    }

  // I read the fixed charge term of the Poisson equation
  p.fixed_charge=dvector(0,Np-1);
  fixed_charge_array=(PyArrayObject *)PyObject_GetAttrString(interface,"fixed_charge");
  for (i=0;i<Np;i++)
    {
      temporary=(double *)(fixed_charge_array->data+i*fixed_charge_array->strides[0]);
      p.fixed_charge[i]=*temporary;
    }


  // ************************************************************************
  // Here it starts the solution of the 3D Poisson equation
  // ************************************************************************

  if (rank==0)
    {
      printf("*********************************\n");
      printf("******* Solving POISSON  ********\n"); 
      printf("*********************************\n\n");
    }


  p.vt=dvector(0,0);
  (*p.vt)=300*kboltz/q;
  ciclomm=2;
  normavp=1e60;
  iterazione=1;
  neq=1;

  // Initialization of the variables needed for the solution of the system
  
  // DOMN variables
  ISYM=0;
  IUNIT=10;
  ITMAX=20000;
  ITOL=1;
  
  comp=matvec_;
  comp2=msolve_;
  NELT=7*Np;
  NSAVE=10;
  IERR=ivector(0,0);

  // Since the following variables are workspace variables,
  // it is better to enlarge the dimension of these variables
  // even if is not explicitly required by DOMN

  P2=dvector(0,3*Np*NSAVE);
  AP=dvector(0,3*Np*NSAVE);
  EMAP=dvector(0,3*Np*NSAVE);

  // Initialization of variables
  
  create_J(&J,d.nx,d.ny,d.nz,neq);
  IA=ivector(0,NELT);
  JA=ivector(0,NELT);
  IWORK=ivector(0,NELT);
  B=dvector(0,Np);
  XLIN=dvector(0,Np);
  R=dvector(0,Np);
  A=dvector(0,NELT);
  DZ=dvector(0,Np);
  CSAV=dvector(0,NSAVE);
  RWORK=dvector(0,NELT);
  Z2=dvector(0,Np);

  p.Phiold=dvector(0,Np-1);

  if (!rank)
    printf("Initialization ...... \n");
  // vector initialization
  for (i=0;i<NELT;i++)
    {
      IA[i]=0;
      JA[i]=0;
      A[i]=0;
    }

  k=0;
  // initialization of the damned variables
  for (j=0;j<=NSAVE;j++)
    {
      for (i=0;i<Np;i++)
	{
	  P2[k]=0.;
	  AP[k]=0.;
	  EMAP[k]=0.;
	  k++;
	}
      CSAV[j]=0.;
    }


  for (i=0;i<Np;i++)
    {
      XLIN[i]=0;
      Z2[i]=0;
      RWORK[i]=0;
      DZ[i]=0;
      R[i]=0;
      IWORK[i]=0;
      p.Phiold[i]=p.Phi[i];
    }

  while (normavp>normapoisson)
    {
      ciclom=1;
      ciclo=1;
      // risolvo POISSON
      gumgum=0;
      ffunpoisson(p,d,B);

      normabold=norma2(B,Np);

      preparej(p,d,A,IA,JA,&Z1,gumgum,J,neq);

      msolve_(&Np,B,XLIN,&NELT,IA,JA,A,&ISYM,RWORK,IWORK);

      
      TOL=tolldomn;
      do{
	
	if (norma2(B,Np)>1e-40)
	  {
	    domn_(&Np,B,XLIN,&Z1,IA,JA,A,&ISYM,comp,comp2,
		  &NSAVE,&ITOL,&TOL,&ITMAX,&ITER,&ERR,IERR,&IUNIT,R,Z2,P2,
		  AP,EMAP,DZ,CSAV,RWORK,IWORK);
	  }
	else
	  {
	    *IERR=0;
	    for (k=0;k<Np;k++)
	      XLIN[k]=0.;
	  }     

	if (*IERR==2)
	  {
            if (!rank)
	      printf("Iterazion domn %d \n",ciclo);
	    TOL=TOL*5;
	    ciclo++;
	  }
      }
      while((ciclo<5)&&(*IERR==2));
      normavp=norma2(XLIN,Np);
      if (!rank)
        printf("NR internal iteration  %d ; norm = %lg \n",iterazione,normavp);
      iterazione++;
      for (i=0;i<Np;i++)
	{
	  p.Phi[i]=p.Phi[i]+(1-sottoril)*XLIN[i];
	}

      // control of the residual
      if (normavp>normapoisson)
	do
	  {
	    ffunpoisson(p,d,B);
	    normab=norma2(B,Np);
            if (!rank)
	      printf("Residual = %lg \n",normab);
	    if (normab>normabold)
	      {
		for (i=0;i<Np;i++)
		  p.Phi[i]-=XLIN[i]*(1-sottoril)/pow(2,ciclom);
                if (!rank)
		  printf("residuo = %lg \n",normab);
	      }
	    ciclom++;
	  }
	while((normab>normabold)&&(ciclom<ciclomm));
    }


  // I return back the Potential
  for (i=0;i<Np;i++)
    *(double *)(Phi_array->data+i*Phi_array->strides[0])=p.Phi[i];

  //  // I set the attribute of the class to send data out
  //  // In particular, I give back the computed potential
  //  temp_obj=  PyArray_FromDimsAndData(1,&Np,PyArray_DOUBLE,(char*)(p.Phi));
  //  PyObject_SetAttrString(interface,"Phi",temp_obj);

  Py_DECREF(surf_array);
  Py_DECREF(eps_array);
  Py_DECREF(boundary_conditions_array);
  Py_DECREF(dist_array);
  Py_DECREF(free_charge_array);
  Py_DECREF(fixed_charge_array);
  Py_DECREF(Phi_array);

  // Free some memory
  free_dvector(p.Phi,0,Np-1);
  free_dvector(p.Phiold,0,Np-1);
  free_dvector(p.eps,0,Np-1);
  free_dvector(d.dVe,0,Np-1);
  free_dmatrix(d.surf,0,Np-1,0,5);
  free_dmatrix(d.dist,0,Np-1,0,5);
  free_dvector(d.boundary_conditions,0,Np-1);
  destroy_J(J,d.nx,d.ny,d.nz,neq);
  free_dvector(P2,0,3*Np*NSAVE);
  free_dvector(AP,0,3*Np*NSAVE);
  free_dvector(EMAP,0,3*Np*NSAVE);
  free_ivector(IA,0,NELT);
  free_ivector(JA,0,NELT);
  free_ivector(IWORK,0,NELT);
  free_dvector(B,0,Np);
  free_dvector(XLIN,0,Np);
  free_dvector(R,0,Np);
  free_dvector(A,0,NELT);
  free_dvector(DZ,0,Np);
  free_dvector(CSAV,0,NSAVE);
  free_dvector(RWORK,0,NELT);
  free_dvector(Z2,0,Np);
  char *s = "Hello from solvePoisson";
  return Py_BuildValue("s", s);
}
