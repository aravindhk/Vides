// ======================================================================
//  Copyright (c) 2009, A.Betti, University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "electron_GNR.h"
static PyObject* py_electron_GNR(PyObject* self, PyObject* args)
{
  int i,j,r,n,m;
  double kF,v;
  double a,M,tempo;
  FILE *fp;
  int N,Nold,dimer,rank;
  double *kx,*ky;
  PyObject *obj,*temp_obj;
  PyArrayObject *energyE_array,*kx_array,*ky_array,*objarray;
  double **energyE; /* electron energy as a function of kx and ky (first and second index, respectively) */
  int n1,dimkx[1],dimky[1],dim2D[2];
  double *punttempo,phi,cutoff,kxmin,kxmax,deltak,kxdown,kxup;
  int mmin,mmax,Nmin,Nmax,delta;
  double thop,aCC;

  import_array(); /* fundamental */

  /* I now read the python object, which is a class */

  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  printf("*************COMPUTATION OF ELECTRON DISPERSION CURVES*************\n");

  temp_obj=PyObject_GetAttrString(obj,"aCC");
  aCC=(double)PyFloat_AsDouble(temp_obj);
 
  temp_obj=PyObject_GetAttrString(obj,"thop");
  thop=(double)PyFloat_AsDouble(temp_obj);

  /* now I extract the attribute N of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"N");
  N=(int)PyInt_AsLong(temp_obj);
  
  if (!rank) 
    printf("total number of grid points N= %d \n",N);

  /* now I extract the attribute dimer of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"dimer");
  dimer=(int)PyInt_AsLong(temp_obj);
  
  if (!rank)
    printf("number of dimer lines= %d \n",dimer);

  /* now I extract the attribute rank of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"mmin");
  mmin=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"mmax");
  mmax=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"kxmin");
  kxmin=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"kxmax");
  kxmax=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"kxdown");
  kxdown=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"kxup");
  kxup=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"phi");
  phi=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"Ecutoff");
  cutoff=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"delta");
  delta=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"deltak");
  deltak=(double)PyFloat_AsDouble(temp_obj);

  if (!rank)
    printf("midgap phi= %lg, cutoff energy Ecutoff= %lg\n",phi,cutoff); 
    

  a=sqrt(3)*aCC;
  
  v=0.12*thop; /* energy correction (in eV) of the hopping parameter at the edges */

  kx=dvector(0,N-1); /* it is defined equal to kxE. It is required for the computation of the scattering rate */

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

  
  kF=kx[N-1]; 

  deltak=delta*(kx[1]-kx[0]); /* delta specifies the number of kx present between two kx points where rates and mobility are computed (sampled) */

  ky=dvector(0,dimer-1); /* quantized vector for electrons */

  ky_array=(PyArrayObject *)PyObject_GetAttrString(obj,"kyE");

  n1=ky_array->dimensions[0];

  dimky[0]=n1;

  for (i=0;i<n1;i++)
    {

      punttempo=(double *)(ky_array->data+i*ky_array->strides[0]);

      ky[i]=*punttempo;
      
    }


  for (i=0;i<(dimer);i++)
    {
      ky[i]=2.0*pi*(i+1)/((dimer+1.0)*a); /* dimer is the number of conduction subbands  */ 
    }

  energyE=dmatrix(0,N-1,0,2*dimer-1);

  energyE_array=(PyArrayObject *)PyObject_GetAttrString(obj,"energyE");

  dim2D[0]=N;
  dim2D[1]=2*dimer;

  /* if (!rank) */
/*     printf("dim2D[0]= %d, dim2D[1]= %d\n",dim2D[0],dim2D[1]); */

  for (i=0;i<(dim2D[0]);i++)
    {
      for (j=0;j<(dim2D[1]);j++)
	{

	  punttempo=(double *)(energyE_array->data+i*energyE_array->strides[0]+j*energyE_array->strides[1]);
	  
	  energyE[i][j]=*punttempo;

	}
    }

  
/*   fp=fopen("curveElectronGNR.out","w"); */

  for (i=0;i<N;i++)
    {

/*       fprintf(fp,"%lg  ",kx[i]/kF); */

      for (m=0;m<dimer;m++)
	{
	 
	  energyE[i][m]=thop*sqrt(1+4*cos(kx[i]*sqrt(3)/2.0*a)*cos(ky[m]*a/2.0)+4*cos(ky[m]*a/2.0)*cos(ky[m]*a/2.0));

/* 	  fprintf(fp,"%lg  ",energyE[i][m]); */
	  
	  /* EDGE RELAXATION */
	  
	  if (cos((pi*(m+1))/(dimer+1.0))>(-0.5))
	    {
	      energyE[i][m]=energyE[i][m]+(4*v/(dimer+1)*sin(((m+1)*pi)/(dimer+1.0))*sin(((m+1)*pi)/(dimer+1.0))*cos(kx[i]*aCC));
	    }
	  else
	    {
	      energyE[i][m]=energyE[i][m]-(4*v/(dimer+1)*sin(((m+1)*pi)/(dimer+1.0))*sin(((m+1)*pi)/(dimer+1.0))*cos(kx[i]*aCC));
	    }
	}

      for (m=dimer;m<(2*dimer);m++)
	{
	 
	  energyE[i][m]=-thop*sqrt(1+4*cos(kx[i]*sqrt(3)/2.0*a)*cos(ky[m-dimer]*a/2.0)+4*cos(ky[m-dimer]*a/2.0)*cos(ky[m-dimer]*a/2.0));

/* 	  fprintf(fp,"%lg  ",energyE[i][m]); */

	  /* EDGE RELAXATION: if the m-th conduction subband shifts upwards, then the m-th valence subband shift downwards */
	  
	  if (cos((pi*(m+1))/(dimer+1.0))>(-0.5))
	    {
	      energyE[i][m]=energyE[i][m]-(4*v/(dimer+1)*sin(((m+1)*pi)/(dimer+1.0))*sin(((m+1)*pi)/(dimer+1.0))*cos(kx[i]*aCC));
	    }
	  else
	    {
	      energyE[i][m]=energyE[i][m]+(4*v/(dimer+1)*sin(((m+1)*pi)/(dimer+1.0))*sin(((m+1)*pi)/(dimer+1.0))*cos(kx[i]*aCC));
	    }

	}


/*       fprintf(fp,"\n"); */
	  
    }
   
/*   fclose(fp); */

  /* case of channel potential different from 0 */

  for (i=0;i<N;i++)
    {
      for (m=0;m<(2*dimer);m++)
	{
	  
	  energyE[i][m]=energyE[i][m]-phi;  /* in eV */
	  
	}
    }
  
  

  Nmin=N-1;
  Nmax=0;

  mmin=dimer-1;
  mmax=0;

  for (m=(dimer/2);m<dimer;m++) /* kyE */
    {
      for (n=0;n<N;n++)  /* kxE */
	{
	  
	  /* only upwards subbands are taken into account */

	  if (n<((N-1)/2))
	    {
	      if (n>0)
		{
		  if ((energyE[n-1][m]>cutoff) & (energyE[n][m]<=cutoff))
		    {
		      
		      if (Nmin>n)
			{
			  Nmin=n;
			}
		    }		  
		}
	      
	    }
	  else if (n>=((N-1)/2))
	    {
	      if ((energyE[n-1][m]<=cutoff) & (energyE[n][m]>cutoff))
		{
		  if (Nmax<(n-1))
		    Nmax=n-1;

		}
		   
	    } 

	  if (energyE[n][m]<cutoff)
	    {
	      if (mmin>m)
		    mmin=m;
	      
	      if (mmax<m)
		mmax=m;
	      
	    }
	  

	}
    }

  if (!rank)
    {
      printf("quantum electron subband index considered in the computation of rates\n");

      printf("mmin= %d,mmax= %d, total # of subbands= %d\n",mmin,mmax,mmax-mmin+1);
    } 
  kxmin=kx[Nmin];

  kxmax=kx[Nmax];

  /* this part could be modified in the parallelized code */

  kxdown=kxmin;
  
  kxup=kxmax;

  if (!rank)
    {
      printf("minimum Nmin and maximum Nmax index of the kx range used for the computation of the scattering rates\n");
      printf("Nmin=%d, Nmax=%d, kxmin/kF=%lg, kxmax/kF=%lg\n",Nmin,Nmax,kx[Nmin]/kF,kx[Nmax]/kF);
      printf("total number of kx points to be computed= %lg\n",Nmax-Nmin);
    }


  temp_obj=PyInt_AsLong(mmin);
  PyObject_SetAttrString(obj,"mmin",temp_obj);

  temp_obj=PyInt_AsLong(mmax);
  PyObject_SetAttrString(obj,"mmax",temp_obj);

  temp_obj=PyFloat_FromDouble(kxmin);  /* double case */
  PyObject_SetAttrString(obj,"kxmin",temp_obj);

  temp_obj=PyFloat_FromDouble(kxmax);  /* double case */
  PyObject_SetAttrString(obj,"kxmax",temp_obj);

  temp_obj=PyFloat_FromDouble(kxdown);  /* double case */
  PyObject_SetAttrString(obj,"kxdown",temp_obj);

  temp_obj=PyFloat_FromDouble(kxup);  /* double case */
  PyObject_SetAttrString(obj,"kxup",temp_obj);

  temp_obj=PyFloat_FromDouble(deltak);  /* double case */
  PyObject_SetAttrString(obj,"deltak",temp_obj);

  temp_obj=  PyArray_FromDimsAndData(1,dimkx,PyArray_DOUBLE,(char*)(&kx[0]));
  PyObject_SetAttrString(obj,"qx",temp_obj);

  temp_obj=  PyArray_FromDimsAndData(1,dimky,PyArray_DOUBLE,(char*)(&ky[0]));
  PyObject_SetAttrString(obj,"kyE",temp_obj);
  
  temp_obj=  PyArray_FromDimsAndData(2,dim2D,PyArray_DOUBLE,(char*)(&energyE[0][0]));
  PyObject_SetAttrString(obj,"energyE",temp_obj);

  char *string = "Bye bye from electron_GNR!";
  return Py_BuildValue("s", string);

}
