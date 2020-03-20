// ======================================================================
//  Copyright (c) 2009, A. Betti, University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "rateOPTEM.h"
static PyObject* py_rateOPTEM(PyObject* self, PyObject* args)
{

  
  int n,m,j,states,i,l,iX,iY,jzf; /* jzf: j for the zone-folded subbranches */
  int abs,Heis,conv,cont,iter;

  double tot,max,maximum,sumOPT,gauss,norm;
  double tr,ti,tpr,tpi,factor; /* real and imaginary parts of t_k, factor due to the spinor */

  double sgn[3];
  double deltaE,deltaEold,deltag,deltag1; 
  double a,W,win;

  FILE *fp,*fp1,*fp2,*fcos;

  double G; /* form factor */
  int mtemp,iYtemp,mid,flag;
  double num1,sum,diff,den;

  int N,dimer,rank;
  double *kxE,*kxP,*kyE,*kyP;
  PyObject *obj,*temp_obj;
  PyArrayObject *energyE_array,*energyP_array,*kx_array,*ky_array,*objarray;
  double **energyE,**energyP2D,**rate,**minAC; /* electron energy as a function of kx and ky (first and second index, respectively) */
  int n1,dimkx[1],dimky[1],dimE2D[2],dimP2D[2],dimAC[2],dimrate[2],dimkyP[1];
  double *punttempo,Ecutoff,kxdown,kxup,deltak;
  int mmin,mmax,Nup,Ndown,deltan;
  double aCC,temp,thop,Dac;


  import_array(); /* fundamental */

  /* I now read the python object, which is a class */

  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  temp_obj=PyObject_GetAttrString(obj,"Dac");
  Dac=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"temp");
  temp=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"thop");
  thop=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"aCC");
  aCC=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"N");
  N=(int)PyInt_AsLong(temp_obj);
  
 /*  if (!rank)  */
/*     printf("N is equal to %d \n",N); */

  /* now I extract the attribute dimer of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"dimer");
  dimer=(int)PyInt_AsLong(temp_obj);
  
  /* if (!rank) */
/*     printf("dimer is equal to %d \n",dimer); */

  /* now I extract the attribute rank of the class through the command PyObject_GetAttrString */

  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"mmin");
  mmin=(int)PyInt_AsLong(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"mmax");
  mmax=(int)PyInt_AsLong(temp_obj);

 
  temp_obj=PyObject_GetAttrString(obj,"kxdown");
  kxdown=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"kxup");
  kxup=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"Ecutoff");
  Ecutoff=(double)PyFloat_AsDouble(temp_obj);

  temp_obj=PyObject_GetAttrString(obj,"deltak");
  deltak=(double)PyFloat_AsDouble(temp_obj);

  kxE=dvector(0,N-1); /* electron wavevector */

  kx_array=(PyArrayObject *)PyObject_GetAttrString(obj,"qx");

  n1=kx_array->dimensions[0];

  dimkx[0]=n1;

  for (i=0;i<n1;i++)
    {

      punttempo=(double *)(kx_array->data+i*kx_array->strides[0]);

      kxE[i]=*punttempo;
      
    }

  kxP=dvector(0,N-1); /* phonon wavevector */

  for (i=0;i<n1;i++)
    {
      
      kxP[i]=kxE[i];
      
    }

  /* conversion from kxdown (kxup) to Ndown (Nup) */


  for (i=0;i<n1;i++)
    {

      if (kxE[i]==kxdown)
	Ndown=i;

      if (kxE[i]==kxup)
	Nup=i;
      
    }

  kyE=dvector(0,dimer-1); /* electron wavevector */

  kx_array=(PyArrayObject *)PyObject_GetAttrString(obj,"kyE");

  dimky[0]=dimer;

  for (i=0;i<(dimky[0]);i++)
    {

      punttempo=(double *)(kx_array->data+i*kx_array->strides[0]);

      kyE[i]=*punttempo;
      
    }

  kyP=dvector(0,dimer-1); /* phonon wavevector */

  kx_array=(PyArrayObject *)PyObject_GetAttrString(obj,"qy");

  dimkyP[0]=dimer;

  for (i=0;i<(dimkyP[0]);i++)
    {

      punttempo=(double *)(kx_array->data+i*kx_array->strides[0]);

      kyP[i]=*punttempo;
      
    }
  
  energyE=dmatrix(0,N-1,0,2*dimer-1);

  energyE_array=(PyArrayObject *)PyObject_GetAttrString(obj,"energyE");

  dimE2D[0]=N;
  dimE2D[1]=2*dimer;

 /*  if (!rank) */
/*     printf("dimE2D[0]= %d, dimE2D[1]= %d\n",dimE2D[0],dimE2D[1]); */

  for (i=0;i<(dimE2D[0]);i++)
    {
      for (j=0;j<(dimE2D[1]);j++)
	{

	  punttempo=(double *)(energyE_array->data+i*energyE_array->strides[0]+j*energyE_array->strides[1]);
	  
	  energyE[i][j]=*punttempo;

	}
    }

  energyP2D=dmatrix(0,N-1,0,6*dimer-1);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"energyP2D");

  dimP2D[0]=N;
  dimP2D[1]=6*dimer;

  /* if (!rank) */
/*     printf("dimP2D[0]= %d, dimP2D[1]= %d\n",dimP2D[0],dimP2D[1]); */

  for (i=0;i<(dimP2D[0]);i++)
    {
      for (j=0;j<(dimP2D[1]);j++)
	{

	  punttempo=(double *)(energyP_array->data+i*energyP_array->strides[0]+j*energyP_array->strides[1]);
	  
	  energyP2D[i][j]=((2*pi*hbar*speed/q)*1e2*(*punttempo));   /* CONVERSION of the phonon energy from cm^-1 to eV */


	}
    }

  minAC=dmatrix(0,dimer-1,0,2);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"minAC");

  dimAC[0]=dimer;
  dimAC[1]=3;

 /*  if (!rank) */
/*     printf("dimAC[0]= %d, dimAC[1]= %d\n",dimAC[0],dimAC[1]); */

  for (i=0;i<(dimAC[0]);i++)
    {
      for (j=0;j<(dimAC[1]);j++)
	{

	  punttempo=(double *)(energyP_array->data+i*energyP_array->strides[0]+j*energyP_array->strides[1]);
	  
	  minAC[i][j]=*punttempo;

	}
    }

  rate=dmatrix(0,N-1,0,dimer-1);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"rateAE");

  dimrate[0]=N;
  dimrate[1]=dimer;

 /*  if (!rank) */
/*     printf("dimrate[0]= %d, dimrate[1]= %d\n",dimrate[0],dimrate[1]); */

  for (i=0;i<(dimrate[0]);i++)
    {
      for (j=0;j<(dimrate[1]);j++)
	{

	  punttempo=(double *)(energyP_array->data+i*energyP_array->strides[0]+j*energyP_array->strides[1]);
	  
	  rate[i][j]=*punttempo;

	}
    }

  deltan=deltak/(kxE[1]-kxE[0]);

  if (!rank)
    {
      printf("*************EMISSION of OPTICAL phonons*************\n");
  
      printf("rank 0: Ndown= %d, Nup= %d, deltan= %d, mmin= %d, mmax= %d \n",Ndown,Nup,deltan,mmin,mmax);
    }

  if (!rank)
    printf("discretization step (DS)= %lg, sampling deltak= %lg, deltan= deltak/DS= %d \n",kxE[1]-kxE[0],deltak,deltan);
 

  win=3.0; /* rectangular (1.0) or Gaussian (3.0) window function */

  a=sqrt(3)*aCC;

  W=(dimer/2.0-1.0/2.0)*a;

  mid=dimer/2;


  /* EMISSION OF OPTICAL PHONONS */


/*   fp1=fopen("momentumOPTEM.out","w"); */


      for (n=(Ndown);n<=(Nup);n=n+deltan)  /* cycle on kxE */
	{

	  if (!rank)
	    if ((n%1)==0)
	      if (!rank) 
		{
		  printf("\r");
		  printf("n= %d",n);
		}

	  cont=0;

	  for (l=0;l<(N-1);l++)  /* cycle on kxP */
	    {

	      /* we exploit N%2!=0 */

	      iX=n-l+(N-1)/2; /* OPT-EM */

	      if (iX>=(N-1))
		iX=iX-(N-1);
	      else if (iX<0)
		iX=iX+(N-1);
	      
	      /* check for backscattering event */

	      if (((kxE[n]>=0) & (kxE[iX]<=0)) || ((kxE[n]<=0) & (kxE[iX]>=0)))
		{

		  for (m=mmin;m<=mmax;m++) /* cycle on kyE */
		    {

		      flag=0;

		      if (energyE[n][m]>Ecutoff)
			flag=1;


		      if (flag==0)
			{
			 for (j=0;j<(dimer);j++)  /* cycle on kyP */
			  {

			  for (iY=0;iY<(dimer);iY++)  /* cycle on final kyE  */
			    {

			      /* the zone-folded subbranches have the same transverse quantum numbers of the original subbranches */

			      if (j>=(mid))
				{
				  jzf=j+2; /* zone-folded subbranches */

				}
			      else
				{
				  jzf=j;  /* original subbranches */
				}

			      
			       /* computation of the form factor */

			      /* 			   printf ("j= %d, jzf= %d\n",j,jzf); */
			  
			      mtemp=m+1;

			      iYtemp=iY+1;

			      num1=4.0*(pi*pi*pi)*jzf*mtemp*iYtemp;

			      sum=mtemp*mtemp+iYtemp*iYtemp;

			      diff=mtemp*mtemp-iYtemp*iYtemp;
			   
			      den=((pi*pi*pi*pi)*(jzf*jzf*jzf*jzf)-2.0*(pi*pi*pi*pi)*(jzf*jzf)*(sum)+(pi*pi*pi*pi)*(diff*diff));
			   
			   
			      if (((mtemp+iYtemp+jzf)%2)!=0) /* odd case */
				{
			       /* l'if a seguire dovrebbe essere sempre soddisfatto */
/* 			       printf ("sum= %d, resto= %d\n",mtemp+iYtemp+jzf,(mtemp+iYtemp+jzf)%2); */

				  if ((jzf!=(mtemp+iYtemp)) & (jzf!=(mtemp-iYtemp)) & (jzf!=(-mtemp+iYtemp)))
				    {
				      G=(2.0*num1*num1*2.0)/(den*den);
				    }
				}
			      else if (((mtemp+iYtemp+jzf)%2)==0)  /* even case */
				{

/* 			       printf ("sum= %d, resto= %d\n",mtemp+iYtemp+jzf,(mtemp+iYtemp+jzf)%2); */


				  if ((jzf==(mtemp-iYtemp)) & (jzf==0))
				    { 
				      G=1.0;
				    }
				  else if (jzf==(mtemp+iYtemp))
				    {
				      G=0.25;

				    }
				  else if ((jzf==(mtemp-iYtemp)) & (jzf!=0))
				    {
				      G=0.25;

				    }
				  else if ((jzf==(-mtemp+iYtemp)) & (jzf!=0))
				    {
				      G=0.25;
			       
				    }
				  else
				    {
				      G=0.0;
				    }
				}


/* 			    printf ("iY= %d\n",iY); */

/* 			   printf ("j= %d, jzf=  %d\n",j,jzf); */


			      if (G!=0.0)

				{

			  
				  /* we have to define the range of eigenvectors k where the delta function is different from zero */

				  for (i=4;i<=5;i++) /* 0-2: AC branches (i=0: ZA, i=1: TA, i=2: LA), 3-5: OPT branches (i=3: ZO, i=4: TO, i=5: LO) */
				    {
			       

				      if ((n>=5) & (n<=((N-1)/2-5))) 
					{
					  deltag1=fabs(energyE[n+5][m]-energyE[n-5][m]);
					}
				      else if ((n>=((N-1)/2+5)) & (n<=(N-6))) 
					{
					  deltag1=fabs(energyE[n+5][m]-energyE[n-5][m]);
					   
					}
				      else if (n>(N-6))
					{
					  deltag1=fabs(energyE[n][m]-energyE[n-10][m]);
					   
					}
				      else if (n<5)
					{
					  deltag1=fabs(energyE[n+10][m]-energyE[n][m]);
					   
					}
				      else if ((n>((N-1)/2-5)) & (n<((N-1)/2+5))) 
					{
					  deltag1=fabs(energyE[(N-1)/2][m]-energyE[(N-1)/2-10][m]);
					}
					     

				      /* check deltag1 for OPT */

				      if (i>=3)
					{
					  if (deltag1<(2e-3))
					    deltag1=2*1e-3;
					}
					     


				      /* finestra rettangolare */

				      deltag=1.0*deltag1;
				      deltaE=win*deltag;

			


			  /* computation of 1+cos(theta) due to the spinor */

			  /* in unit of the hopping parameter thop */


				    /* si deve stare attenti ad applicare lo zone-folding alle sottobranche piu' alte con concavita' verso il basso */

/* 				    if (m>=(mid)) */
				      {
					tr=1.0+2.0*cos((kxE[n])*aCC*3.0/2.0)*cos((kyE[m])*aCC*sqrt(3.0)/2.0);
					ti=-2.0*sin((kxE[n])*aCC*3.0/2.0)*cos((kyE[m])*aCC*sqrt(3.0)/2.0);
				
				      }
				   
/* 				    if (iY>=(mid)) */
				      {
					
					tpr=1.0+2.0*cos((kxE[iX])*aCC*3.0/2.0)*cos((kyE[iY])*aCC*sqrt(3.0)/2.0);
					tpi=2.0*sin((kxE[iX])*aCC*3.0/2.0)*cos((kyE[iY])*aCC*sqrt(3.0)/2.0); /*  conjugate complex */

				      }
				  



			    factor=1.0+(tr*tpr-ti*tpi)/(sqrt(tr*tr+ti*ti)*sqrt(tpr*tpr+tpi*tpi)); 

			  /* we have to define the range of eigenvectors k where the delta function is different from zero */

					 
/* 			    if (energyEG[g][iX][iY]<energyEG[g][n][m]) */
			      {
			

				if (fabs(energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i])<(deltaE))
				  {

				    norm=1.0/sqrt(pi);

				    gauss=norm*exp(-fabs((energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i])*(energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i]))/(fabs(((deltag)*(deltag)))));
					   

/* 				    gauss=1.0; */


				    sumOPT=(((hbar*Dac*Dac*deltak*pi*pi*2.0)/(a*a*rho*W)*1.0/((q*energyP2D[l][j*6+i])*(q*deltag))*(Bose_Einstein((q*energyP2D[l][j*6+i])/(kboltz*temp))+1.0))*(1.0-Fermi_Dirac((q*energyE[iX][iY])/(kboltz*temp)))/(1.0-Fermi_Dirac((q*energyE[n][m])/(kboltz*temp)))*gauss)*factor*G; /* OPT: we include the factor taking into account the case of degenerate semiconductor */


				    rate[n][m]+=sumOPT;  
					     

				   /*  if ((m==1) & (j==1) & (i==3))   */
/* 				      fprintf(fp1,"%lg %lg %lg\n",kxE[n]/kF,kxP[l]/kF,energyEG[g][n][m]); */

				  }  /* if |Efin-Ein+Eph|<deltaE */

			      }  /* if Efin<Ein */   
				   
			    } /* end cycle on i (graphene branches) */

			  } /* if G!=0.0 */

		       } /*  cycle on iY (allowed transverse vectors) */

		      }  /* cicle on kyP */

		    } /* if flag==0 */

		  }  /* cicle on kyE */

		} /* backscattering */

	    }  /* end cycle on kxP */

	} /* kxE */


/*        fclose(fp1); */

      fp=fopen("rateOPTEM.out","w");

      for (n=0;n<N;n++)  /* cycle on kxE */
	{
	  for (m=0;m<dimer;m++)  /* cycle on kyE */
	    {
	      
	      fprintf(fp," %lg %lg ",energyE[n][m],rate[n][m]);
	      
	    }
	  

	  fprintf(fp,"\n");

	}
		

      fclose(fp);

      temp_obj=  PyArray_FromDimsAndData(2,dimrate,PyArray_DOUBLE,(char*)(&rate[0][0]));
      PyObject_SetAttrString(obj,"rateOE",temp_obj);
     
      free_dmatrix(energyE,0,N-1,0,2*dimer-1);

      free_dmatrix(energyP2D,0,N-1,0,6*dimer-1);

      free_dmatrix(minAC,0,dimer-1,0,2);

      printf("\n");

      char *string = "Bye bye from rateOPTEM!";
      return Py_BuildValue("s", string);


}
