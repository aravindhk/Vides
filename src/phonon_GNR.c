// ======================================================================
//  Copyright (c) 2009, A. Betti, University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "phonon_GNR.h"
static PyObject* py_phonon_GNR(PyObject* self, PyObject* args)
{
 
  double K11[3][3],K12[3][3],K13[3][3];
  double K21[3][3],K22[3][3],K23[3][3],K24[3][3],K25[3][3],K26[3][3];
  double K31[3][3],K32[3][3],K33[3][3];
  double K41[3][3],K42[3][3],K43[3][3],K44[3][3],K45[3][3],K46[3][3];

  double Phi_r1,Phi_ti1,Phi_to1;
  double Phi_r2,Phi_ti2,Phi_to2;
  double Phi_r3,Phi_ti3,Phi_to3;
  double Phi_r4,Phi_ti4,Phi_to4;


  int i,j,r,n,m,s,**flag,**orderx,***ordery,cost1,cost2,l;
  complex **D,zero,**eigenvector;
  double kF,*ky1;
  double a,M,tempo,energyold;
  double *eigenvalue,*frequency,***energy;
  FILE *fp,*fp1,*fpAC,*fpOT,*fporder;
  
  PyObject *obj,*temp_obj;
  PyArrayObject *energyP_array,*kx_array,*objarray;

  int N,Nold,dimer,rank,Ny,Nyold;
  double *kx,*ky,*qy;
  double ***energyP,**energyP2D,**minAC;
  int n1,n2,n3,n4,dim[3],dimkx[1],dim2D[2],dimAC[2],dimqy[1];
  double *punttempo;
  int dim1,dim2,dim3;
  double aCC;


  import_array(); /* fundamental */

  /* I now read the python object, which is a class */

  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  printf("*************COMPUTATION OF PHONON DISPERSION CURVES*************\n");

  temp_obj=PyObject_GetAttrString(obj,"aCC");
  aCC=(double)PyFloat_AsDouble(temp_obj);

  /* now I extract the attribute N of the class through the command PyObject_GetAttrString */
 
  temp_obj=PyObject_GetAttrString(obj,"N");
  N=(int)PyInt_AsLong(temp_obj);

 /*  if (!rank) */
/*     printf("N is equal to %d \n",N); */

  /* now I extract the attribute dimer of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"dimer");
  dimer=(int)PyInt_AsLong(temp_obj);

 /*  if (!rank) */
/*     printf("dimer is equal to %d \n",dimer); */

  /* now I extract the attribute rank of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"dim1");
  dim1=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"dim2");
  dim2=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"dim3");
  dim3=(int)PyInt_AsLong(temp_obj);


  /* definition of force constant parameters */

  temp_obj=PyObject_GetAttrString(obj,"Phi_r1");
  Phi_r1=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_ti1");
  Phi_ti1=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_to1");
  Phi_to1=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_r2");
  Phi_r2=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_ti2");
  Phi_ti2=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_to2");
  Phi_to2=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_r3");
  Phi_r3=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_ti3");
  Phi_ti3=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_to3");
  Phi_to3=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_r4");
  Phi_r4=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_ti4");
  Phi_ti4=(double)PyFloat_AsDouble(temp_obj);
  temp_obj=PyObject_GetAttrString(obj,"Phi_to4");
  Phi_to4=(double)PyFloat_AsDouble(temp_obj);



  Ny=140; /* it would be sufficient for every GNR widths W<=10 nm: Ny=80 is the minimum for W= 1.12, 4.86 nm, Ny=140 is the minimum for W=2.62 nm and Ny=100 for W=10.10 nm */

  Nyold=Ny;
  
  
  /* if (!rank) */
/*     printf("/// number of grid points: N= %d ///\n",N);  */
    

  while ((((Ny-1)%(dimer+1))!=0) || ((Ny%2)==0)) 
    {
      Ny++;
    }

 /*  if (Ny!=Nyold) */
/*     { */
/*       if (!rank) */
/* 	printf("/// redefinition of the number of grid points: Ny= %d ///\n",Ny);  */
/*     } */
/*   else */
/*     { */
/*       if (!rank) */
/* 	printf("/// number of grid points: Ny= %d ///\n",Ny);  */
/*     } */

  M=weight/Nav;

  a=sqrt(3.0)*aCC;

  kx=dvector(0,N-1); /* it is defined equal to kxE. It is required for the computation of the scattering rate */
  
  ky=dvector(0,dimer+1); /* quantized vector for phonons */
  
  kx_array=(PyArrayObject *)PyObject_GetAttrString(obj,"qx");

  n1=kx_array->dimensions[0];

  dimkx[0]=n1;

  for (i=0;i<n1;i++)
    {

      punttempo=(double *)(kx_array->data+i*kx_array->strides[0]);

      kx[i]=*punttempo;
      
    }



  for (i=0;i<N;i++)
    {
      kx[i]= 2.0*pi/(sqrt(3.0)*a)*(i)/(N-1)-pi/(sqrt(3.0)*a);  /* for phonons: used for the computation of the scattering rate: kxP=kxE */
    }
  
  for (i=0;i<=(dimer+1);i++)
    {
      ky[i]=2.0*pi*(i)/((dimer+1.0)*a); /* dimer*6 is the number of phonon subbranches */ 

    }

  qy=dvector(0,dimer-1); /* transverse phonon wavevector */
  
  kx_array=(PyArrayObject *)PyObject_GetAttrString(obj,"qy");

  n1=kx_array->dimensions[0];

  dimqy[0]=n1;

  for (i=0;i<n1;i++)
    {

      punttempo=(double *)(kx_array->data+i*kx_array->strides[0]);

      qy[i]=*punttempo;
      
    }

  for (i=0;i<(dimer/2);i++)
    {
      qy[i]=ky[i]; /* dimer*6 is the number of phonon subbranches */ 

    }

  for (i=(dimer/2);i<(dimer);i++)
    {
      qy[i]=ky[i+2]; /* dimer*6 is the number of phonon subbranches */ 

    }

  ky1=dvector(0,Ny-1);

  flag=imatrix(0,5,0,5);
  orderx=imatrix(0,N-1,0,5); /*such matrix contains for each n=0,...,N-1 (kxn)  the order of the phonon branches for ky=0. In particular for each i the first three terms refer to number of the eigenvalues corresponding to acoustic branches. Instead the last three terms refer the optical branches */

  ordery=ivectorm(0,2); /* vector of int matrix: order of the phonon branches as a function of kx and ky: the first index is kx, the second one is ky, the third one is the number of the branch */

  energy=dvectorm(0,2); /* vector of double matrix: the first index is kx, the second is ky, and the third is the graphene branch (energy in cm^-1) */
 /*  energyP=dvectorm(0,N-1); */ /* it includes only the quantized vectors for the GNR */

  energyP=dvectorm(0,N-1); /* it includes only the quantized vectors for phonons (GNR) */
  
  for (n=0;n<N;n++)
    {
      
      energyP[n]=dmatrix(0,(dimer)-1,0,5);

    }


  energyP2D=dmatrix(0,N-1,0,dimer*6-1);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"energyP2D");

 /*  n1=energyP_array->dimensions[0]; */
/*   n2=energyP_array->dimensions[2]; */

  dim2D[0]=dim1;
  dim2D[1]=dim2*dim3;

  /* if (!rank) */
/*     printf("dim2D[0]= %d, dim2D[1]= %d\n",dim2D[0],dim2D[1]); */

  for (i=0;i<(dim2D[0]);i++)
    {
      for (j=0;j<(dim2D[1]);j++)
	{

	  punttempo=(double *)(energyP_array->data+i*energyP_array->strides[0]+j*energyP_array->strides[1]);
	  
	  energyP2D[i][j]=*punttempo;

	}
    }

  minAC=dmatrix(0,dimer-1,0,2);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"minAC");

  dimAC[0]=dimer;
  dimAC[1]=3;

  /* if (!rank) */
/*     printf("dimAC[0]= %d, dimAC[1]= %d\n",dimAC[0],dimAC[1]); */

  for (i=0;i<(dimAC[0]);i++)
    {
      for (j=0;j<(dimAC[1]);j++)
	{

	  punttempo=(double *)(energyP_array->data+i*energyP_array->strides[0]+j*energyP_array->strides[1]);
	  
	  minAC[i][j]=*punttempo;

	}
    }

  

  /* if (!rank) */
/*     printf("dim[0]= %d, dim[1]= %d, dim[2]= %d\n",dim[0],dim[1],dim[2]); */

  for (i=0;i<dim1;i++)
    {
      for (j=0;j<dim2;j++)
	{
	  for (n=0;n<dim3;n++)
	    {

	      energyP[i][j][n]=0.0;

	    }
	}
    }

  eigenvalue=dvector(0,5);
  frequency=dvector(0,5);


  eigenvector=cmatrix(0,5,0,5);
  D=cmatrix(0,5,0,5);

  for (n=0;n<=2;n++)
    {
      energy[n]=dmatrix(0,Ny-1,0,5);

      ordery[n]=imatrix(0,Ny-1,0,5); /* the initial order for each n (kxn) is the obvious order of the eigenvalues defined by the zhpevx function. It is the right order only if the branches do not intersect each other */
/*       energyP[n]=dmatrix(0,(dimer/2)-1,0,5); */

    }
  

  zero.r=0.0;
  zero.i=0.0;

  for (n=0;n<N;n++)
    {
      for (j=0;j<=5;j++)
	{
          orderx[n][j]=j; /* the initial order for each n (kxn) is the obvious order of the eigenvalues defined by the zhpevx function. It is the right order only if the branches do not intersect each other */
	}
    }

  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
          K11[i][j]=0.0;
          K12[i][j]=0.0;
          K13[i][j]=0.0;


          K21[i][j]=0.0;
          K22[i][j]=0.0;
          K23[i][j]=0.0;
          K24[i][j]=0.0;
          K25[i][j]=0.0;
          K26[i][j]=0.0;

          K31[i][j]=0.0;
          K32[i][j]=0.0;
          K33[i][j]=0.0;

          K41[i][j]=0.0;
          K42[i][j]=0.0;
          K43[i][j]=0.0;
          K44[i][j]=0.0;
          K45[i][j]=0.0;
          K46[i][j]=0.0;

	}
    }


 
  /* ALLOCATION K */

  /* first neighbours */

  K11[0][0]=Phi_r1;
  K11[1][1]=Phi_ti1;
  K11[2][2]=Phi_to1;

  
  K12[0][0]=1.0/4.0*Phi_r1+3.0/4.0*Phi_ti1;
  K12[0][1]=-sqrt(3.0)/4.0*Phi_r1+sqrt(3.0)/4.0*Phi_ti1;
  K12[1][0]=-sqrt(3.0)/4.0*Phi_r1+sqrt(3.0)/4.0*Phi_ti1;
  K12[1][1]=3.0/4.0*Phi_r1+1.0/4.0*Phi_ti1;
  K12[2][2]=Phi_to1;

  K13[0][0]=1.0/4.0*Phi_r1+3.0/4.0*Phi_ti1;
  K13[0][1]=sqrt(3.0)/4.0*Phi_r1-sqrt(3.0)/4.0*Phi_ti1;
  K13[1][0]=sqrt(3.0)/4.0*Phi_r1-sqrt(3.0)/4.0*Phi_ti1;
  K13[1][1]=3.0/4.0*Phi_r1+1.0/4.0*Phi_ti1;
  K13[2][2]=Phi_to1;

  /* second neighbours */

  K21[0][0]=3.0/4.0*Phi_r2+1.0/4.0*Phi_ti2;
  K21[0][1]=sqrt(3.0)/4.0*Phi_r2-sqrt(3.0)/4.0*Phi_ti2;
  K21[1][0]=sqrt(3.0)/4.0*Phi_r2-sqrt(3.0)/4.0*Phi_ti2;
  K21[1][1]=1.0/4.0*Phi_r2+3.0/4.0*Phi_ti2;
  K21[2][2]=Phi_to2;

  K22[0][0]=Phi_ti2;
  K22[1][1]=Phi_r2;
  K22[2][2]=Phi_to2;

  K23[0][0]=3.0/4.0*Phi_r2+1.0/4.0*Phi_ti2;
  K23[0][1]=-sqrt(3.0)/4.0*Phi_r2+sqrt(3.0)/4.0*Phi_ti2;
  K23[1][0]=-sqrt(3.0)/4.0*Phi_r2+sqrt(3.0)/4.0*Phi_ti2;
  K23[1][1]=1.0/4.0*Phi_r2+3.0/4.0*Phi_ti2;
  K23[2][2]=Phi_to2;

  K24[0][0]=3.0/4.0*Phi_r2+1.0/4.0*Phi_ti2;
  K24[0][1]=sqrt(3.0)/4.0*Phi_r2-sqrt(3.0)/4.0*Phi_ti2;
  K24[1][0]=sqrt(3.0)/4.0*Phi_r2-sqrt(3.0)/4.0*Phi_ti2;
  K24[1][1]=1.0/4.0*Phi_r2+3.0/4.0*Phi_ti2;
  K24[2][2]=Phi_to2;

  K25[0][0]=Phi_ti2;
  K25[1][1]=Phi_r2;
  K25[2][2]=Phi_to2;

  K26[0][0]=3.0/4.0*Phi_r2+1.0/4.0*Phi_ti2;
  K26[0][1]=-sqrt(3.0)/4.0*Phi_r2+sqrt(3.0)/4.0*Phi_ti2;
  K26[1][0]=-sqrt(3.0)/4.0*Phi_r2+sqrt(3.0)/4.0*Phi_ti2;
  K26[1][1]=1.0/4.0*Phi_r2+3.0/4.0*Phi_ti2;
  K26[2][2]=Phi_to2;

  /* third neighbours */

  K31[0][0]=1.0/4.0*Phi_r3+3.0/4.0*Phi_ti3;
  K31[0][1]=sqrt(3.0)/4.0*Phi_r3-sqrt(3.0)/4.0*Phi_ti3;
  K31[1][0]=sqrt(3.0)/4.0*Phi_r3-sqrt(3.0)/4.0*Phi_ti3;
  K31[1][1]=3.0/4.0*Phi_r3+1.0/4.0*Phi_ti3;
  K31[2][2]=Phi_to3;

  K32[0][0]=Phi_r3;
  K32[1][1]=Phi_ti3;
  K32[2][2]=Phi_to3;

  K33[0][0]=1.0/4.0*Phi_r3+3.0/4.0*Phi_ti3;
  K33[0][1]=-sqrt(3.0)/4.0*Phi_r3+sqrt(3.0)/4.0*Phi_ti3;
  K33[1][0]=-sqrt(3.0)/4.0*Phi_r3+sqrt(3.0)/4.0*Phi_ti3;
  K33[1][1]=3.0/4.0*Phi_r3+1.0/4.0*Phi_ti3;
  K33[2][2]=Phi_to3;

  /* fourth neighbours */

  K41[0][0]=25.0/28.0*Phi_r4+3.0/28.0*Phi_ti4;
  K41[0][1]=5.0*sqrt(3.0)/28.0*Phi_r4-5.0*sqrt(3.0)/28.0*Phi_ti4;
  K41[1][0]=5.0*sqrt(3.0)/28.0*Phi_r4-5.0*sqrt(3.0)/28.0*Phi_ti4;
  K41[1][1]=3.0/28.0*Phi_r4+25.0/28.0*Phi_ti4;
  K41[2][2]=Phi_to4;

  K42[0][0]=1.0/28.0*Phi_r4+27.0/28.0*Phi_ti4;
  K42[0][1]=-3.0*sqrt(3.0)/28.0*Phi_r4+3.0*sqrt(3.0)/28.0*Phi_ti4;
  K42[1][0]=-3.0*sqrt(3.0)/28.0*Phi_r4+3.0*sqrt(3.0)/28.0*Phi_ti4;
  K42[1][1]=27.0/28.0*Phi_r4+1.0/28.0*Phi_ti4;
  K42[2][2]=Phi_to4;

  K43[0][0]=4.0/7.0*Phi_r4+3.0/7.0*Phi_ti4;
  K43[0][1]=-2.0*sqrt(3.0)/7.0*Phi_r4+2.0*sqrt(3.0)/7.0*Phi_ti4;
  K43[1][0]=-2.0*sqrt(3.0)/7.0*Phi_r4+2.0*sqrt(3.0)/7.0*Phi_ti4;
  K43[1][1]=3.0/7.0*Phi_r4+4.0/7.0*Phi_ti4;
  K43[2][2]=Phi_to4;

  K44[0][0]=4.0/7.0*Phi_r4+3.0/7.0*Phi_ti4;
  K44[0][1]=2.0*sqrt(3.0)/7.0*Phi_r4-2.0*sqrt(3.0)/7.0*Phi_ti4;
  K44[1][0]=2.0*sqrt(3.0)/7.0*Phi_r4-2.0*sqrt(3.0)/7.0*Phi_ti4;
  K44[1][1]=3.0/7.0*Phi_r4+4.0/7.0*Phi_ti4;
  K44[2][2]=Phi_to4;

  K45[0][0]=1.0/28.0*Phi_r4+27.0/28.0*Phi_ti4;
  K45[0][1]=3.0*sqrt(3.0)/28.0*Phi_r4-3.0*sqrt(3.0)/28.0*Phi_ti4;
  K45[1][0]=3.0*sqrt(3.0)/28.0*Phi_r4-3.0*sqrt(3.0)/28.0*Phi_ti4;
  K45[1][1]=27.0/28.0*Phi_r4+1.0/28.0*Phi_ti4;
  K45[2][2]=Phi_to4;

  K46[0][0]=25.0/28.0*Phi_r4+3.0/28.0*Phi_ti4;
  K46[0][1]=-5.0*sqrt(3.0)/28.0*Phi_r4+5.0*sqrt(3.0)/28.0*Phi_ti4;
  K46[1][0]=-5.0*sqrt(3.0)/28.0*Phi_r4+5.0*sqrt(3.0)/28.0*Phi_ti4;
  K46[1][1]=3.0/28.0*Phi_r4+25.0/28.0*Phi_ti4;
  K46[2][2]=Phi_to4;



  /* allocation of the dynamical matrix D */

  /* n,m are the index for the cycle on the wavevector */



 /* graphene */


  fp=fopen("curveGM.out","w");
  fp1=fopen("orderky0.out","w");

  for (i=0;i<Ny;i++)
    {
/*       kx[i]= 2*pi/(sqrt(3)*a)*(i)/(N-1);  */

      ky1[i]=0.0;

    }


  kF=kx[N-1];



  for (n=0;n<N;n++)
    {

	  /* D_AA block */

      for (i=0;i<=2;i++)
	    {
	      for (j=0;j<=2;j++)
		{

		  D[i][j].r=(K11[i][j]+K12[i][j]+K13[i][j])+(K21[i][j]+K22[i][j]+                     K23[i][j]+K24[i][j]+K25[i][j]+K26[i][j])+(K31[i][j]+                            K32[i][j]+K33[i][j])+(K41[i][j]+K42[i][j]+K43[i][j]+                            K44[i][j]+K45[i][j]+K46[i][j]);  



		  D[i][j].i=0;

		  D[i][j].r=D[i][j].r-(K21[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])-sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0]))+K22[i][j]*cos(ky1[0]*a)+K23[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])+sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0]))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])-sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0]))+K25[i][j]*cos(ky1[0]*a)+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])+sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0])));

		  D[i][j].i=D[i][j].i-(K21[i][j]*(-sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])-cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0]))+K22[i][j]*(-sin(ky1[0]*a))+K23[i][j]*(-cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0])+sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0]))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0])+sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0]))+K25[i][j]*(sin(ky1[0]*a))+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[0])-sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[0])));

		}
	    }
    
	

	  /* D_AB block */

  
	  
	  for (i=0;i<=2;i++)
	    {
	      for (j=3;j<=5;j++)
		{

		  D[i][j].r=-(K11[i][j-3]*cos(a/sqrt(3.0)*kx[n])+K12[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])+sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K13[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])-sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K31[i][j-3]*(cos(a/sqrt(3.0)*kx[n])*cos(a*ky1[0])-sin(a/sqrt(3.0)*kx[n])*sin(a*ky1[0]))+K32[i][j-3]*cos(a*2.0/(sqrt(3.0))*kx[n])+K33[i][j-3]*(cos(a/sqrt(3.0)*kx[n])*cos(a*ky1[0])+sin(a/sqrt(3.0)*kx[n])*sin(a*ky1[0]))+K41[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K42[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[0])+sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[0]))+K43[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[0])+sin(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K44[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[0])-sin(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K45[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[0])-sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[0]))+K46[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])+sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0])));

		  D[i][j].i=-(K11[i][j-3]*(-sin(a/sqrt(3.0)*kx[n]))+K12[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])-cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K13[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])+cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K31[i][j-3]*(-sin(a/(sqrt(3.0))*kx[n])*cos(a*ky1[0])-cos(a/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K32[i][j-3]*sin(a*2.0/(sqrt(3.0))*kx[n])+K33[i][j-3]*(-sin(a/(sqrt(3.0))*kx[n])*cos(a*ky1[0])+cos(a/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K41[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])-cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0]))+K42[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[0])-cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[0]))+K43[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[0])-cos(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K44[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[0])+cos(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[0]))+K45[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[0])+cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[0]))+K46[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[0])+cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[0])));   

		 
		}
	    }

	  /* D_BA block */

	  for (i=3;i<=5;i++)
	    {
	      for (j=0;j<=2;j++)
		{
    
		 /*  D[i][j].r=D[i-3][j+3].r; */
/*                   D[i][j].i=-D[i-3][j+3].i; */

		  D[i][j].r=(D[j][i].r);
                  D[i][j].i=-(D[j][i].i);
		}
	    }


	  /* D_BB block */

	  for (i=3;i<=5;i++)
	    {
	      for (j=3;j<=5;j++)
		{
    
		  D[i][j].r=D[i-3][j-3].r;
                  D[i][j].i=D[i-3][j-3].i;
		  
		}
	    }

	  /* if (!rank) */
/* 	    if ((n%50)==0) */
/* 	      printf("kx=%lg, rank= %d\n",kx[n]/kF,rank); */


	  ceig(D,eigenvector,eigenvalue,6,6);


	  for (i=0;i<=5;i++)
	    {
	      if (eigenvalue[i]>=0)
		frequency[i]=sqrt(eigenvalue[i]/M);
              else
		frequency[i]=0;

	      if (n==0)
		{
		energy[0][0][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}
		else if (n==1)
		{
		  energy[1][0][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}
	      else if (n==2)
		{
		  energy[2][0][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}
	      else if (n>2)
		{
		  energy[0][0][i]=energy[1][0][i];
		  energy[1][0][i]=energy[2][0][i];
		  energy[2][0][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}


	    }             

          for (i=0;i<=5;i++)
	    {

	      for (j=i+1;j<=5;j++)
		{

		  
		  if (n>=2)
		    {
		      if((fabs(((energy[0][0][i]-energy[1][0][i])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][0][i])-energy[2][0][i])>fabs(((energy[0][0][i]-energy[1][0][i])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][0][i])-energy[2][0][j]))&&(fabs(((energy[0][0][j]-energy[1][0][j])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][0][j])-energy[2][0][j])>fabs(((energy[0][0][j]-energy[1][0][j])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][0][j])-energy[2][0][i])))
			{
			

                          flag[i][j]=1; /* flag[i][j]=1 means intersection between the eigenvalues i and j */
			  flag[j][i]=1;

                          orderx[n][i]=j; /* at the vector kxn the index of the eigenvalue corresponding to the branch i becomes j, and the index referred to the branch j becomes i */

			  orderx[n][j]=i;
			}
			
		    }
		}
	    }

		
		
	  if (!rank)
	    {
	      if (n>=2)
		{
		
		  fprintf(fp,"%lg %lg %lg %lg %lg %lg %lg \n",kx[n]/kF,energy[2][0][orderx[n][0]],energy[2][0][orderx[n][1]],energy[2][0][orderx[n][2]],energy[2][0][orderx[n][3]],energy[2][0][orderx[n][4]],energy[2][0][orderx[n][5]]);
		}
	      else if (n==0)
		{
		  fprintf(fp,"%lg %lg %lg %lg %lg %lg %lg \n",kx[n]/kF,energy[0][0][orderx[n][0]],energy[0][0][orderx[n][1]],energy[0][0][orderx[n][2]],energy[0][0][orderx[n][3]],energy[0][0][orderx[n][4]],energy[0][0][orderx[n][5]]);
		}
	       else if (n==1)
		{
		  fprintf(fp,"%lg %lg %lg %lg %lg %lg %lg \n",kx[n]/kF,energy[1][0][orderx[n][0]],energy[1][0][orderx[n][1]],energy[1][0][orderx[n][2]],energy[1][0][orderx[n][3]],energy[1][0][orderx[n][4]],energy[1][0][orderx[n][5]]);
		}
	      
	     
	      fprintf(fp1,"%lg %d %d %d %d %d %d  \n",kx[n]/kF,orderx[n][0],orderx[n][1],orderx[n][2],orderx[n][3],orderx[n][4],orderx[n][5]);
		
	    }


    }

  fclose(fp);

  fclose(fp1);
  /* at this point for each kxn in the n-th row of the matrix order the first three index are the numbers of the eigenvalues corresponding to AC branches, the last three index are referred instead to the OPT branches */ 


  kF=kx[N-1]; 

  for (n=0;n<Ny;n++)
    {


      ky1[n]=ky[1]*(dimer+1)*n/(Ny-1.0);

   
      for (i=0;i<=(dimer+1);i++)
	{
         
	  /* here we exploit the fact that (Ny-1)%(dimer-1)=0 for Ip. */

	  if (((dimer+1)*n/(Ny-1.0)-i)==0)
	  {
	      
	    ky1[n]=ky[i]; /* i valori i=dimer/2 e dimer/2+1 non sono buoni */

	    /* if (!rank) */
/* 	      printf("ky1=ky??? %lg, ky1= %lg, n=%d, i=%d\n", (dimer+1.0)*n/(Ny-1.0),ky1[n],n,i); /\* I impose the equality between ky1 and ky (ky is the quantized transverse vector for the GNR device *\/ */
	  }
	    
	}

    }


  for (n=0;n<N;n++)  /* cycle on kx */
    {
      
      /* for each n (kxn) we have to set to zero the matrix flag which defines the eventual intersection of the phonon branches (flag[i][j]=1 means intersection between the eigenvalues i and j) */
      

      for (m=0;m<Ny;m++)  /* cycle on ky */
	{
	  for (j=0;j<=5;j++)
	    {
	      if (n>2)
		{
		  ordery[0][m][j]=ordery[1][m][j];  /* aggiornamento dati vecchi */
		  ordery[1][m][j]=ordery[2][m][j];

		}

	      if (n==0)
		{
		  ordery[0][m][j]=orderx[n][j]; /* at fixed kxn the initial order of the eigenvalues along the ky direction is defined by the n-th row of the matrix orderx */
		}
	      else if (n==1)
		{
		  ordery[1][m][j]=orderx[n][j]; 
		}
	      else if (n>=2)
		{
		  ordery[2][m][j]=orderx[n][j]; 
		}

	    }
	}
    
     
	  

      for (m=0;m<Ny;m++)  /* cycle on ky */
	{

	  if (m>0)
	    {
	      for (r=m;r<Ny;r++)  /* cycle on ky */
		{
		  for (j=0;j<=5;j++)
		    {
		      if (n==0)
			{
			  ordery[0][r][j]=ordery[0][m-1][j]; /* at fixed kxn the initial order of the eigenvalues along the ky direction at kym is defined by the m-1 row of the matrix ordery */
			}
		      else if (n==1)
			{
			  ordery[1][r][j]=ordery[1][m-1][j]; 
			}
		      else if (n>=2)
			{
			  ordery[2][r][j]=ordery[2][m-1][j]; 
			}

		    }
		}
	    }

	  /* D_AA block */

	  for (i=0;i<=2;i++)
	    {
	      for (j=0;j<=2;j++)
		{

		  D[i][j].r=(K11[i][j]+K12[i][j]+K13[i][j])+(K21[i][j]+K22[i][j]+                     K23[i][j]+K24[i][j]+K25[i][j]+K26[i][j])+(K31[i][j]+                            K32[i][j]+K33[i][j])+(K41[i][j]+K42[i][j]+K43[i][j]+                            K44[i][j]+K45[i][j]+K46[i][j]);  



		  D[i][j].i=0.0;

		  D[i][j].r=D[i][j].r-(K21[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m])-sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m]))+K22[i][j]*cos(ky1[m]*a)+K23[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m])+sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m]))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m])-sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m]))+K25[i][j]*cos(ky1[m]*a)+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m])+sin(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m])));

		  D[i][j].i=(D[i][j].i-(K21[i][j]*(-sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m])-cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m]))+K22[i][j]*(-sin(ky1[m]*a))+K23[i][j]*(-cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m])+sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m]))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m])+sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m]))+K25[i][j]*(sin(ky1[m]*a))+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx[n])*sin(a/2.0*ky1[m])-sin(sqrt(3.0)/2.0*a*kx[n])*cos(a/2.0*ky1[m]))));

		}
	    }
    
	

	  /* D_AB block */

  
	  
	  for (i=0;i<=2;i++)
	    {
	      for (j=3;j<=5;j++)
		{


		  D[i][j].r=-(K11[i][j-3]*cos(a/sqrt(3.0)*kx[n])+K12[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])+sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K13[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])-sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K31[i][j-3]*(cos(a/sqrt(3.0)*kx[n])*cos(a*ky1[m])-sin(a/sqrt(3.0)*kx[n])*sin(a*ky1[m]))+K32[i][j-3]*cos(a*2.0/(sqrt(3.0))*kx[n])+K33[i][j-3]*(cos(a/sqrt(3.0)*kx[n])*cos(a*ky1[m])+sin(a/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K41[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K42[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[m])+sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[m]))+K43[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[m])+sin(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K44[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[m])-sin(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K45[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[m])-sin(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[m]))+K46[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])+sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m])));    

 

		  D[i][j].i=-(K11[i][j-3]*(-sin(a/sqrt(3.0)*kx[n]))+K12[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])-cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K13[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])+cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K31[i][j-3]*(-sin(a/(sqrt(3.0))*kx[n])*cos(a*ky1[m])-cos(a/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K32[i][j-3]*sin(a*2.0/(sqrt(3.0))*kx[n])+K33[i][j-3]*(-sin(a/(sqrt(3.0))*kx[n])*cos(a*ky1[m])+cos(a/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K41[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])-cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m]))+K42[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[m])-cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[m]))+K43[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[m])-cos(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K44[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx[n])*cos(a*ky1[m])+cos(a*2.0/(sqrt(3.0))*kx[n])*sin(a*ky1[m]))+K45[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx[n])*cos(a*3.0/2.0*ky1[m])+cos(a/(2.0*sqrt(3.0))*kx[n])*sin(a*3.0/2.0*ky1[m]))+K46[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx[n])*cos(a/2.0*ky1[m])+cos(a*5.0/(2.0*sqrt(3.0))*kx[n])*sin(a/2.0*ky1[m])));   

		}
	    }

	  /* D_BA block */

	  for (i=3;i<=5;i++)
	    {
	      for (j=0;j<=2;j++)
		{
    
		  D[i][j].r=(D[j][i].r);
                  D[i][j].i=-(D[j][i].i);
		}
	    }


	  /* D_BB block */

	  for (i=3;i<=5;i++)
	    {
	      for (j=3;j<=5;j++)
		{
    
		  D[i][j].r=D[i-3][j-3].r;
                  D[i][j].i=D[i-3][j-3].i;
		  
		}
	    }

	  if (!rank)
	    if ((m==0) & ((n%50)==0))
	      printf("kx=%lg  //GNR//\n",kx[n]/kF);


	  ceig(D,eigenvector,eigenvalue,6,6);

	  for (i=0;i<=5;i++)
	    {
	      if (eigenvalue[i]>=0)
		frequency[i]=sqrt(eigenvalue[i]/M);
              else
		frequency[i]=0;



	      if (n==0)
		{
		  energy[0][m][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		 
		}
	      else if (n==1)
		{
		  energy[1][m][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */

		}
	      else if (n==2)
		{
		  energy[2][m][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}
	      else if (n>2)
		{
		  energy[0][m][i]=energy[1][m][i];
		  energy[1][m][i]=energy[2][m][i];
		  energy[2][m][i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
		}

	    }             

          
	  for (i=0;i<=5;i++)
	    {
	      
	      for (j=i+1;j<=5;j++)
		{

		  
		  if (m>=2)
		    {
		      
		      if (n>=2)
			{
			  if((fabs(((energy[2][m-2][ordery[2][m-2][i]]-energy[2][m-1][ordery[2][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[2][m-1][ordery[2][m-1][i]])-energy[2][m][ordery[2][m][i]])>fabs(((energy[2][m-2][ordery[2][m-2][i]]-energy[2][m-1][ordery[2][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[2][m-1][ordery[2][m-1][i]])-energy[2][m][ordery[2][m][j]])) && (fabs(((energy[2][m-2][ordery[2][m-2][j]]-energy[2][m-1][ordery[2][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[2][m-1][ordery[2][m-1][j]])-energy[2][m][ordery[2][m][j]])>fabs(((energy[2][m-2][ordery[2][m-2][j]]-energy[2][m-1][ordery[2][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[2][m-1][ordery[2][m-1][j]])-energy[2][m][ordery[2][m][i]])))  
			    {
			  
			      /* eigenvalues exchange */
			      
			      tempo=ordery[2][m][j];
			      ordery[2][m][j]=ordery[2][m][i];
			      ordery[2][m][i]=tempo;

			  

			    }
			}
		      else if (n==0)
			{
			  if((fabs(((energy[0][m-2][ordery[0][m-2][i]]-energy[0][m-1][ordery[0][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[0][m-1][ordery[0][m-1][i]])-energy[0][m][ordery[0][m][i]])>fabs(((energy[0][m-2][ordery[0][m-2][i]]-energy[0][m-1][ordery[0][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[0][m-1][ordery[0][m-1][i]])-energy[0][m][ordery[0][m][j]])) && (fabs(((energy[0][m-2][ordery[0][m-2][j]]-energy[0][m-1][ordery[0][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[0][m-1][ordery[0][m-1][j]])-energy[0][m][ordery[0][m][j]])>fabs(((energy[0][m-2][ordery[0][m-2][j]]-energy[0][m-1][ordery[0][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[0][m-1][ordery[0][m-1][j]])-energy[0][m][ordery[0][m][i]])))  
			    {
			  
			      /* eigenvalues exchange */
			      
			      tempo=ordery[0][m][j];
			      ordery[0][m][j]=ordery[0][m][i];
			      ordery[0][m][i]=tempo;
			      
			  
			    }
			}
		      else if (n==1)
			{
			  if((fabs(((energy[1][m-2][ordery[1][m-2][i]]-energy[1][m-1][ordery[1][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[1][m-1][ordery[1][m-1][i]])-energy[1][m][ordery[1][m][i]])>fabs(((energy[1][m-2][ordery[1][m-2][i]]-energy[1][m-1][ordery[1][m-1][i]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[1][m-1][ordery[1][m-1][i]])-energy[1][m][ordery[1][m][j]])) && (fabs(((energy[1][m-2][ordery[1][m-2][j]]-energy[1][m-1][ordery[1][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[1][m-1][ordery[1][m-1][j]])-energy[1][m][ordery[1][m][j]])>fabs(((energy[1][m-2][ordery[1][m-2][j]]-energy[1][m-1][ordery[1][m-1][j]])/(ky1[m-2]-ky1[m-1])*(ky1[m]-ky1[m-1])+energy[1][m-1][ordery[1][m-1][j]])-energy[1][m][ordery[1][m][i]])))  
			    {
			  
			      /* eigenvalues exchange */
			      
			      tempo=ordery[1][m][j];
			      ordery[1][m][j]=ordery[1][m][i];
			      ordery[1][m][i]=tempo;
			      
			  
			    }
			}
		      



		    }

			 
		}
	    }

	  /* AAAA check branches along the kx direction AAAA*/


	 for (i=0;i<=5;i++)
	    {
	      
	      for (j=i+1;j<=5;j++)
		{

		  if (n>=2)
		    {
		      
		      if ((fabs(((energy[0][m][ordery[0][m][i]]-energy[1][m][ordery[1][m][i]])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][m][ordery[1][m][i]])-energy[2][m][ordery[2][m][i]])>fabs(((energy[0][m][ordery[0][m][i]]-energy[1][m][ordery[1][m][i]])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][m][ordery[1][m][i]])-energy[2][m][ordery[2][m][j]]))&&(fabs(((energy[0][m][ordery[0][m][j]]-energy[1][m][ordery[1][m][j]])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][m][ordery[1][m][j]])-energy[2][m][ordery[2][m][j]])>fabs(((energy[0][m][ordery[0][m][j]]-energy[1][m][ordery[1][m][j]])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energy[1][m][ordery[1][m][j]])-energy[2][m][ordery[2][m][i]])))
			{
			  tempo=ordery[2][m][j];
			  ordery[2][m][j]=ordery[2][m][i];
			  ordery[2][m][i]=tempo;

				  
			}
		    }
		}
	    }



	  for (i=0;i<(dimer/2);i++)
	    {
             

	      if ((ky1[m]==ky[i]) && ((n%50)==0))
		{
		  if (!rank)
		    printf("ky quantized (m=%d), quantum #: i= %d, ky= %lg, ky/kF=%lg\n",m,i,ky1[m],ky1[m]/ky1[Ny-1]);
		}


	      if (ky1[m]==ky[i])
		{ 
		  
		  if (n>=2)
		    {
		      energyP[n][i][0]=energy[2][m][ordery[2][m][0]];
		      energyP[n][i][1]=energy[2][m][ordery[2][m][1]];
		      energyP[n][i][2]=energy[2][m][ordery[2][m][2]];
		      energyP[n][i][3]=energy[2][m][ordery[2][m][3]];
		      energyP[n][i][4]=energy[2][m][ordery[2][m][4]];
		      energyP[n][i][5]=energy[2][m][ordery[2][m][5]];
		    }
		  else if (n==0)
		    {
		      energyP[n][i][0]=energy[0][m][ordery[0][m][0]];
		      energyP[n][i][1]=energy[0][m][ordery[0][m][1]];
		      energyP[n][i][2]=energy[0][m][ordery[0][m][2]];
		      energyP[n][i][3]=energy[0][m][ordery[0][m][3]];
		      energyP[n][i][4]=energy[0][m][ordery[0][m][4]];
		      energyP[n][i][5]=energy[0][m][ordery[0][m][5]];
		    }
		   else if (n==1)
		    {
		      energyP[n][i][0]=energy[1][m][ordery[1][m][0]];
		      energyP[n][i][1]=energy[1][m][ordery[1][m][1]];
		      energyP[n][i][2]=energy[1][m][ordery[1][m][2]];
		      energyP[n][i][3]=energy[1][m][ordery[1][m][3]];
		      energyP[n][i][4]=energy[1][m][ordery[1][m][4]];
		      energyP[n][i][5]=energy[1][m][ordery[1][m][5]];
		    }
		  

		  /* check AC branch with i=0 and corresponding to the zero branch of the graphene sheet */

		  if (N>=0)
		    {
		      if (n>=2)
			{
			  energyold=energyP[n][0][0];

			  if (kx[n]<0)
			    {
			      if(((energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1]))>((energyP[n-1][0][0]-energyP[n][0][0])/(kx[n-1]-kx[n])))
				{
				  energyP[n][0][0]=(energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energyP[n-1][0][0];
				  
				}
			      else if((energyP[n-1][0][0]-energyP[n][0][0])/(kx[n-1]-kx[n])>0)
				{
				  energyP[n][0][0]=(energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energyP[n-1][0][0];
				}
			      
			    }
			  else if (kx[n]>0)
			    {
			
			      if(energyP[n][0][0]!=energyP[N-1-n][0][0])
				{
				  energyP[n][0][0]=energyP[N-1-n][0][0];
				}
			    }
			}

		    }
		  
		
		}
	    }
	    

	  for (i=(dimer/2);i<(dimer);i++)
	    {
             

	      if ((ky1[m]==ky[i+2]) && ((n%50)==0))
		{
		  
		  if (!rank)
		    printf("ky quantized (m=%d), quantum #: i= %d, ky= %lg, ky/kF=%lg\n",m,i,ky1[m],ky1[m]/ky1[Ny-1]);

		}


	      if (ky1[m]==ky[i+2])
		{ 
		  
		  if (n>=2)
		    {
		      energyP[n][i][0]=energy[2][m][ordery[2][m][0]];
		      energyP[n][i][1]=energy[2][m][ordery[2][m][1]];
		      energyP[n][i][2]=energy[2][m][ordery[2][m][2]];
		      energyP[n][i][3]=energy[2][m][ordery[2][m][3]];
		      energyP[n][i][4]=energy[2][m][ordery[2][m][4]];
		      energyP[n][i][5]=energy[2][m][ordery[2][m][5]];
		    }
		  else if (n==0)
		    {
		      energyP[n][i][0]=energy[0][m][ordery[0][m][0]];
		      energyP[n][i][1]=energy[0][m][ordery[0][m][1]];
		      energyP[n][i][2]=energy[0][m][ordery[0][m][2]];
		      energyP[n][i][3]=energy[0][m][ordery[0][m][3]];
		      energyP[n][i][4]=energy[0][m][ordery[0][m][4]];
		      energyP[n][i][5]=energy[0][m][ordery[0][m][5]];
		    }
		  else if (n==1)
		    {
		      energyP[n][i][0]=energy[1][m][ordery[1][m][0]];
		      energyP[n][i][1]=energy[1][m][ordery[1][m][1]];
		      energyP[n][i][2]=energy[1][m][ordery[1][m][2]];
		      energyP[n][i][3]=energy[1][m][ordery[1][m][3]];
		      energyP[n][i][4]=energy[1][m][ordery[1][m][4]];
		      energyP[n][i][5]=energy[1][m][ordery[1][m][5]];
		    }
		  		    

		  /* check AC branch with i=0 and corresponding to the zero branch of the graphene sheet */

		  if (N>=0)
		    {
		      if (n>=2)
			{
			  energyold=energyP[n][0][0];

			  if (kx[n]<0)
			    {
			      if(((energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1]))>((energyP[n-1][0][0]-energyP[n][0][0])/(kx[n-1]-kx[n])))
				{
				  energyP[n][0][0]=(energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energyP[n-1][0][0];
				  
				}
			      else if((energyP[n-1][0][0]-energyP[n][0][0])/(kx[n-1]-kx[n])>0)
				{
				  energyP[n][0][0]=(energyP[n-2][0][0]-energyP[n-1][0][0])/(kx[n-2]-kx[n-1])*(kx[n]-kx[n-1])+energyP[n-1][0][0];
				}
			      
			    }
			  else if (kx[n]>0)
			    {
			
			      if(energyP[n][0][0]!=energyP[N-1-n][0][0])
				{
				  energyP[n][0][0]=energyP[N-1-n][0][0];
				}
			    }
			}

		    }
		  

		
		}
	    }
	 
	  
	}  /* cycle on ky1 */

    }  /* cycle on kx */


  /* fp=fopen("GNRsubbranches.out","w"); */

/*   for (i=0;i<N;i++) */
/*     { */

/*       fprintf(fp,"%lg  ",kx[i]/kF); */

/*       for (j=0;j<dimer;j++) */
/* 	{ */
/* 	  for (n=0;n<=5;n++) */
/* 	    { */
	      
/* 	      fprintf(fp,"%lg  ",energyP[i][j][n]); */
/* 	    } */
/* 	} */
      
/*       fprintf(fp,"\n"); */
/*     } */
  

/*   fclose(fp); */

  
  for (i=0;i<N;i++) 
    {

      for (j=0;j<dimer;j++)
	{
	  for (n=0;n<=5;n++)
	    {
	      m=j*6+n;

	      energyP2D[i][m]=energyP[i][j][n];
	    }
	}
   }


  /* determination of the minimum energy for each acoustic subbranches */
 
  for (m=0;m<(dimer);m++)
    {
      for (i=0;i<=2;i++)
	{
	  minAC[m][i]=1e6;
	}
    }


  for (m=0;m<(dimer);m++) /* cycle on kyP */
    {
      for (i=0;i<=2;i++)
	{
	  for (l=0;l<N;l++)  /* cycle on kxP */
	    {
	      if (minAC[m][i]>energyP[l][m][i])
		{
		  minAC[m][i]=energyP[l][m][i];
		}
	
	    }
	}
    }
  

  free_dvectorm(energy,0,2,0,Ny-1,0,5);

  free_ivectorm(ordery,0,2,0,Ny-1,0,5);

  free_imatrix(orderx,0,N-1,0,5);

  // I set the attribute of the class to send data out


  // I set the attribute of the class to send data out

  temp_obj=  PyArray_FromDimsAndData(1,dimkx,PyArray_DOUBLE,(char*)(&kx[0]));
  PyObject_SetAttrString(obj,"qx",temp_obj);
  
  temp_obj=  PyArray_FromDimsAndData(1,dimqy,PyArray_DOUBLE,(char*)(&qy[0]));
  PyObject_SetAttrString(obj,"qy",temp_obj);
  
  temp_obj=  PyArray_FromDimsAndData(2,dim2D,PyArray_DOUBLE,(char*)(&energyP2D[0][0]));
  PyObject_SetAttrString(obj,"energyP2D",temp_obj);
  
  temp_obj=  PyArray_FromDimsAndData(2,dimAC,PyArray_DOUBLE,(char*)(&minAC[0][0]));
  PyObject_SetAttrString(obj,"minAC",temp_obj);

 /*  energyP_array=  (PyArrayObject *) PyArray_FromDims(2,dim2D,PyArray_DOUBLE); */
  
/*   for (i=0;i<N;i++)  */
/*     { */

/*       for (j=0;j<dimer;j++) */
/* 	{ */
/* 	  for (n=0;n<=5;n++) */
/* 	    { */
/* 	      m=j*6+n; */

/* 	      (energyP_array->data+i*energyP_array->strides[0]+m*energyP_array->strides[1])= (char *) (&energyP2D[i][m]); */

/* 	    } */
/* 	} */
/*     } */

/*   PyObject_SetAttrString(obj,"energyP2D",temp_obj); */


  char *string = "Bye bye from phonon_GNR!";
  return Py_BuildValue("s", string);


}
