// ======================================================================
//  Copyright (c) 2009, A. Betti, University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "phonon_graphene.h"
static PyObject* py_phonon_graphene(PyObject* self, PyObject* args)
{
 
  double K11[3][3],K12[3][3],K13[3][3];
  double K21[3][3],K22[3][3],K23[3][3],K24[3][3],K25[3][3],K26[3][3];
  double K31[3][3],K32[3][3],K33[3][3];
  double K41[3][3],K42[3][3],K43[3][3],K44[3][3],K45[3][3],K46[3][3];

  double Phi_r1,Phi_ti1,Phi_to1;
  double Phi_r2,Phi_ti2,Phi_to2;
  double Phi_r3,Phi_ti3,Phi_to3;
  double Phi_r4,Phi_ti4,Phi_to4;


  int i,j,r,n,m,s,cost1,cost2,l;
  complex **D,zero,**eigenvector;
  double a,M,tempo;
  double *eigenvalue,*frequency;
  FILE *fp;
  
  PyObject *obj,*temp_obj;
  PyArrayObject *energy_array,*objarray;

  double kx,ky;
  double *energy;
  int n1,dimE[1];
  double *punttempo;
  int dim3,rank;
  double aCC;


  import_array(); /* fundamental */

  /* I now read the python object, which is a class */

  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  temp_obj=PyObject_GetAttrString(obj,"aCC");
  aCC=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"dim3");
  dim3=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"qx0");
  kx=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"qy0");
  ky=(double)PyFloat_AsDouble(temp_obj);

 if (!rank)
   printf("qx0= %lg, qy0= %lg\n",kx,ky);
  
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

  M=weight/Nav;

  a=sqrt(3.0)*aCC;

  energy=dvector(0,dim3-1); /* it include the graphene phonon energies for a specified (kx,ky) */
  
  energy_array=(PyArrayObject *)PyObject_GetAttrString(obj,"Egraphene");

  n1=energy_array->dimensions[0];

  dimE[0]=n1;

  for (i=0;i<n1;i++)
    {

      punttempo=(double *)(energy_array->data+i*energy_array->strides[0]);

      energy[i]=*punttempo;
      
    }


  eigenvalue=dvector(0,5);
  frequency=dvector(0,5);


  eigenvector=cmatrix(0,5,0,5);
  D=cmatrix(0,5,0,5);

 

  zero.r=0.0;
  zero.i=0.0;

  

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


    
  /* D_AA block */

  for (i=0;i<=2;i++)
    {
      for (j=0;j<=2;j++)
	{
	  
	  D[i][j].r=(K11[i][j]+K12[i][j]+K13[i][j])+(K21[i][j]+K22[i][j]+                     K23[i][j]+K24[i][j]+K25[i][j]+K26[i][j])+(K31[i][j]+                            K32[i][j]+K33[i][j])+(K41[i][j]+K42[i][j]+K43[i][j]+                            K44[i][j]+K45[i][j]+K46[i][j]);  
	  
	  
	  
	  D[i][j].i=0;
	  
	  D[i][j].r=D[i][j].r-(K21[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)-sin(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky))+K22[i][j]*cos(ky*a)+K23[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)+sin(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)-sin(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky))+K25[i][j]*cos(ky*a)+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)+sin(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky)));
	  
	  D[i][j].i=D[i][j].i-(K21[i][j]*(-sin(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)-cos(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky))+K22[i][j]*(-sin(ky*a))+K23[i][j]*(-cos(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky)+sin(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky))+K24[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky)+sin(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky))+K25[i][j]*(sin(ky*a))+K26[i][j]*(cos(sqrt(3.0)/2.0*a*kx)*sin(a/2.0*ky)-sin(sqrt(3.0)/2.0*a*kx)*cos(a/2.0*ky)));
	  
	}
    }
  
	

  /* D_AB block */

  
	  
  for (i=0;i<=2;i++)
    {
      for (j=3;j<=5;j++)
	{
	  
	  D[i][j].r=-(K11[i][j-3]*cos(a/sqrt(3.0)*kx)+K12[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)+sin(a/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K13[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)-sin(a/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K31[i][j-3]*(cos(a/sqrt(3.0)*kx)*cos(a*ky)-sin(a/sqrt(3.0)*kx)*sin(a*ky))+K32[i][j-3]*cos(a*2.0/(sqrt(3.0))*kx)+K33[i][j-3]*(cos(a/sqrt(3.0)*kx)*cos(a*ky)+sin(a/sqrt(3.0)*kx)*sin(a*ky))+K41[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)-sin(a*5.0/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K42[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx)*cos(a*3.0/2.0*ky)+sin(a/(2.0*sqrt(3.0))*kx)*sin(a*3.0/2.0*ky))+K43[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx)*cos(a*ky)+sin(a*2.0/(sqrt(3.0))*kx)*sin(a*ky))+K44[i][j-3]*(cos(a*2.0/(sqrt(3.0))*kx)*cos(a*ky)-sin(a*2.0/(sqrt(3.0))*kx)*sin(a*ky))+K45[i][j-3]*(cos(a/(2.0*sqrt(3.0))*kx)*cos(a*3.0/2.0*ky)-sin(a/(2.0*sqrt(3.0))*kx)*sin(a*3.0/2.0*ky))+K46[i][j-3]*(cos(a*5.0/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)+sin(a*5.0/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky)));
		  
	  D[i][j].i=-(K11[i][j-3]*(-sin(a/sqrt(3.0)*kx))+K12[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)-cos(a/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K13[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)+cos(a/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K31[i][j-3]*(-sin(a/(sqrt(3.0))*kx)*cos(a*ky)-cos(a/(sqrt(3.0))*kx)*sin(a*ky))+K32[i][j-3]*sin(a*2.0/(sqrt(3.0))*kx)+K33[i][j-3]*(-sin(a/(sqrt(3.0))*kx)*cos(a*ky)+cos(a/(sqrt(3.0))*kx)*sin(a*ky))+K41[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)-cos(a*5.0/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky))+K42[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx)*cos(a*3.0/2.0*ky)-cos(a/(2.0*sqrt(3.0))*kx)*sin(a*3.0/2.0*ky))+K43[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx)*cos(a*ky)-cos(a*2.0/(sqrt(3.0))*kx)*sin(a*ky))+K44[i][j-3]*(sin(a*2.0/(sqrt(3.0))*kx)*cos(a*ky)+cos(a*2.0/(sqrt(3.0))*kx)*sin(a*ky))+K45[i][j-3]*(sin(a/(2.0*sqrt(3.0))*kx)*cos(a*3.0/2.0*ky)+cos(a/(2.0*sqrt(3.0))*kx)*sin(a*3.0/2.0*ky))+K46[i][j-3]*(-sin(a*5.0/(2.0*sqrt(3.0))*kx)*cos(a/2.0*ky)+cos(a*5.0/(2.0*sqrt(3.0))*kx)*sin(a/2.0*ky)));   

		 
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
  
  
  ceig(D,eigenvector,eigenvalue,6,6);

  printf("*************GRAPHENE PHONON ENERGY*************\n");

  printf("energy (cm^-1):  ");

  for (i=0;i<=5;i++)
    {
      	 
      if (eigenvalue[i]>=0)
	frequency[i]=sqrt(eigenvalue[i]/M);
      else
	frequency[i]=0;

      energy[i]=frequency[i]/(2*pi*speed)*1e-2;  /* cm^-1 */
      
      printf("  %lg  ",energy[i]);

    }
  
  printf("\n");

  temp_obj=  PyArray_FromDimsAndData(1,dimE,PyArray_DOUBLE,(char*) (energy));
  PyObject_SetAttrString(obj,"Egraphene",temp_obj);

  char *string = "Bye bye from phonon_graphene!";
  return Py_BuildValue("s", string);


}
