// ======================================================================
//  Copyright (c) 2009, A. Betti,  G. Fiori, University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
#include "rateACEM.h"
static PyObject* py_rateACEM(PyObject* self, PyObject* args)
{

  int n,m,j,states,i,l,iX,iY,r,rstar,jzf; /* jzf: j for the zone-folded subbranches */
  int abs,conv,cont,iter;

  double tot,max,maximum,sumAC,gauss;
  double tr,ti,tpr,tpi,factor; /* real and imaginary parts of t_k, factor due to the spinor */

  double sgn[4];
  double deltaE,deltaEold,deltag,deltag1,deltain; 
  double a,W,cutoff,cutoff1,cutoff2,win;
  double deltagmax,deltagmin,deltaEmax,deltaEmin;
  double totold,totmax,totmin,deltagold,number,x0,slope,broad,vs;
  double slopeold,slopestar,norm,c;

  FILE *fp,*fp1,*fp2;

  double Heis,tempG; /* coefficients for the Heisenberg's principle */
  int pippo,pluto;

  double G; /* form factor */
  int mtemp,iYtemp,mid,flag,band;
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
  int numberAC,imin,imax;


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
  
 /*  if (!rank) */
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

  temp_obj=PyObject_GetAttrString(obj,"numberAC");
  numberAC=(int)PyInt_AsLong(temp_obj);

  if (numberAC==2)
    {
      imin=1;
      imax=2;
    }
  else if (numberAC==1)
    {
      imin=2;
      imax=2;
    }
  else
    {
      if (!rank)
	printf("error in numberAC!\n");

      exit(1);
    }
  
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

  /* if (!rank) */
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

/*   if (!rank) */
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

  rate=dmatrix(0,N-1,0,dimer-1);

  energyP_array=(PyArrayObject *)PyObject_GetAttrString(obj,"rateAE");

  dimrate[0]=N;
  dimrate[1]=dimer;

/*   if (!rank) */
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
      printf("*************EMISSION of ACOUSTIC phonons*************\n");

      printf("rank 0: Ndown= %d, Nup= %d, deltan= %d, mmin= %d, mmax= %d \n",Ndown,Nup,deltan,mmin,mmax);
       
    }

  if (!rank)
    printf("discretization step (DS)= %lg, sampling deltak= %lg, deltan= deltak/DS= %d \n",kxE[1]-kxE[0],deltak,deltan);
 
 /*  win=1.0; /\* win is the factor related to the window function: 0.5 for the rectangular window, 3 for the Gaussian window *\/ */

  cutoff=10.0; /* in cm^-1: we apply the self-consistency only to phonon subbranches with the energy minimum smaller than the cutoff energy */

  cutoff1=10.0;

  cutoff2=50.0; /* in cm^-1: we apply the self-consistency only to phonon subbranches with the energy minimum smaller than the cutoff energy */

  a=sqrt(3)*aCC;

  W=(dimer/2.0-1.0/2.0)*a;

  mid=dimer/2;

 
/*       if (g%10==0) */

/*       fp=fopen("momentumACEM.out","w"); */


  for (n=(Ndown);n<=(Nup);n=n+deltan)  /* cycle on kxE */
    {


	  if (!rank)
	    if ((n%1)==0)
	      if (!rank) 
		{
		  printf("\r");
		  printf("n= %d",n);
		}

/*           band=0; */

	  for (m=mmin;m<=(mmax);m++) /* cycle on kyE */
	    {
	      /* assumo che lo step di qx sia uguale a quello di kx */

	      flag=0;

	      if (energyE[n][m]>Ecutoff)
		flag=1;

	   /*    else */
/* 		band++; */

	      /* assumo che lo step di qx sia uguale a quello di kx */
	      
/* 	      if ((flag==0) & (!iproc) & (band==1)) */
/* 		if ((n%1)==0) */
/* 		  printf("%d\n",n); */
	      
	      if (flag==0)
		{
		  for (j=(0);j<(dimer);j++)  /* cycle on kyP */
		    {

		      cont=0;

		      sumAC=0.0;

		      
		      for (iY=0;iY<(dimer);iY++)  /* cycle on final kyE  */
			  {

			    /* 			   iY=m; */


			      /* the zone-folded subbranches have transverse quantum numbers l/2+2,.....,l+1 */

			     
			    if (j>=(mid))
			      {

				jzf=j+2; /* zone-folded subbranches */

/* 			       printf ("j= %d, jzf= %d, dimer/2= %d\n",j,jzf,dimer/2); */

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

				  for (i=imin;i<=imax;i++) /* 0-2: AC branches (i=0: ZA, i=1: TA, i=2: LA), 3-5: OPT branches (i=3: ZO, i=4: TO, i=5: LO) */
				    {


				       /* 11 points kxE enclosed in the delta width */

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
					     

				      /* finestra rettangolare */

				       
				     /*  if ((i==0) & ((minAC[j][i])<=(cutoff1*0.124*1e-3))) */

				      if (((minAC[j][i])<=(cutoff*0.124*1e-3)))
					{

					  /* self-consistent ZA mode: Rect function */

					  win=0.5;

					  c=1.5180;  /* stop at the third zero of the Sinc function */

					  Heis=1.0/hbar*sqrt(2.0*c/pi); /* Gamma=Heis*DeltaE */

					  deltagmax=1.0*deltag1;
				       
					  if (deltagmax>5e-3)
					    deltagmax=5e-3; /* 2e-3 per gauss, 5e-3 per rect */

					  deltaEmax=win*deltagmax; /* gauss. */
					
					  deltagmin=1e-7;

					  deltaEmin=win*deltagmin;

					  deltag=deltagmax;

					  deltaE=win*deltag;
					}
				      else
					{
					  /* Gaussian function */

					  win=3.0;

					  Heis=2.0/hbar;

					  deltag=1.0*deltag1;
					  
					  deltaE=win*deltag;

					}


			
				abs=1;	/* ci sono energie con Efin>Ein */
				     
				iter=0;

				conv=0;  /* no convergenza di deltaE */
			
				tot=0;
				totold=0.0;

				totmax=0.0;
				totmin=0.0;

				deltagold=0.0;
				deltaEold=0.0;


				pippo=0; /* pippo=1 significa che ho un caso in cui il metodo della bisezione fallisce perche' Gamma_max-deltaEmax/(hbar*c)<0, pippo=0 significa tutto ok */

				pluto=0; /* serve per la stampa dei casi in cui ho la disuguaglianza precedente con il <0 (pippo=1), e fallisco oppure ho successo nella ricerca del deltaEmax giusto */


				while ((conv==0))
				  {

				    abs=0;
				    
				    tot=0;
				    
				    sumAC=0.0;

				    number=0;
				    

				    if ((i==0) & ((minAC[j][i])>(cutoff2*0.124*1e-3)))
				      {
					conv=1;
				      }

				    if ((i!=0) & ((minAC[j][i])>(cutoff*0.124*1e-3)))
				      {
					conv=1;

				      }

				    
				    for (l=0;l<(N-1);l++)  /* cycle on kxP */
				      {
					

					/* given kxE[n] and kxP[l], we have to compute the corresponding index iX of the final kxE[iX]=kxE[n]+kxP[l] */

					/* we exploit N%2!=0 */

					iX=n-l+(N-1)/2; /* AC-EM */

					if (iX>=(N-1))
					  iX=iX-(N-1);
					else if (iX<0)
					  iX=iX+(N-1);

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

					/* check for backscattering event, for a fixed kxE[n] and kxE[iX] */

					if (((kxE[n]>=0) & (kxE[iX]<=0)) || ((kxE[n]<=0) & (kxE[iX]>=0)))
					  {

					    /* we have to define the range of eigenvectors k where the delta function is different from zero */

/* 					    if (energyEG[g][iX][iY]<energyEG[g][n][m]) */
					      {

						if (fabs(energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i])<(deltaE))
						  {
					   
   
						    number++;

						    abs=1;


						   /*  if ((i==0) & ((minAC[j][i])<=(cutoff1*0.124*1e-3))) */

						    if (((minAC[j][i])<=(cutoff*0.124*1e-3)))
						      {
							gauss=1.0; /* RECT */
						      }
						    else
						      {
							norm=1.0/sqrt(pi); /* GAUSS */

							gauss=norm*exp(-fabs((energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i])*(energyE[iX][iY]-energyE[n][m]+energyP2D[l][j*6+i]))/(fabs(((deltag)*(deltag)))));
					   
						      }


						    /* momentum relaxation rate accounting for only backscattering */
				       
						    if ((j!=0) || ((j==0)  & (l!=((N-1)/2))))
						      {
							sumAC=((hbar*Dac*Dac*deltak)/(2.0*rho*W)*(kxP[l]*kxP[l]+kyP[j]*kyP[j])/((q*energyP2D[l][j*6+i])*(q*deltag))*(Bose_Einstein((q*energyP2D[l][j*6+i])/(kboltz*temp))+1.0)*(1.0-Fermi_Dirac((q*energyE[iX][iY])/(kboltz*temp)))/(1.0-Fermi_Dirac((q*energyE[n][m])/(kboltz*temp)))*gauss)*factor*G; 
						      }
						    else if ((j==0) & (l==((N-1)/2)))
						      {
							vs=q*(energyP2D[l+1][j*6+i]-energyP2D[l][j*6+i])/(kxP[l+1]-kxP[l])*1.0/hbar;

							sumAC=((hbar*Dac*Dac*deltak)/(2.0*rho*W)*1.0/((q*deltag))*(kboltz*temp/(hbar*hbar*vs*vs))*(1.0-Fermi_Dirac((q*energyE[iX][iY])/(kboltz*temp)))/(1.0-Fermi_Dirac((q*energyE[n][m])/(kboltz*temp)))*gauss)*factor*G; 	
						      } 


						    tot+=sumAC;



						  }
					      }

					  } /* backscattering */

				      } /* end cycle on kxP */
		

				    if ((i==0) & ((minAC[j][i])<=(cutoff1*0.124*1e-3)))
				      {

					deltagold=deltag;
					deltaEold=deltaE;

					if ((iter==0) & (abs==0))
					  {

					    conv=1;

					   /*  if (!iproc) */
/* 					      if ((pippo==1) & (abs==0)) */
/* 						printf("ERROR: pippo= %d, conv= %d, abs=%d, n= %d, m= %d, states= %d, j= %d, max=%lg, tot= %lg, totmax(prec)= %lg, tot-deltaEmax/(hbar*c)=%lg, sumAC= %lg, sgn[states]= %lg\n",pippo,conv,abs,n,m,states,j,deltagmax,tot,totmax,tot-deltagmax*q*Heis,sumAC,sgn[states]); */
					    
					    pippo=0;

					    pluto=0;

					  }

					if (G==0.0)
					  conv=1;


					if ((iter==0) & (abs==1))
					  {

					    if (G!=0.0)
					      {
						if ((tot-deltagmax*q*Heis)<0)
						  {

						    pippo=1;
						    pluto=1;

						    totmax=tot;

						    tempG=deltagmax; /* old */

						    deltagmax=deltagmax-0.02*tempG;
						    deltaEmax=win*deltagmax;
					    
						    deltag=deltagmax;

						    deltaE=win*deltag;

						   /*  if (!iproc) */
						   /*    printf("NO, abs=%d, n= %d, m= %d, states= %d, j= %d, max=%lg, totmax= %lg, difmax=%lg, sumAC= %lg, sgn[states]= %lg\n",abs,n,m,states,j,temp,totmax,tot-temp*q*Heis,sumAC,sgn[states]); */

						  }
						else
						  {
						   
						   
						    pippo=0;

						    totmin=0.0;
						    totmax=tot;

/* 						    if (!iproc) */
						    /*   if (pluto==1) */
/* 							{ */
/* 							  printf("SI, abs=%d, n= %d, m= %d, states= %d, j= %d, max=%lg, totmax= %lg, difmax=%lg, sumAC= %lg, sgn[states]= %lg\n",abs,n,m,states,j,deltagmax,totmax,totmax-deltagmax*q*Heis,sumAC,sgn[states]); */
/* 							} */

						    pluto=0;

						    x0=(deltagmin+deltagmax)/2.0;;

						    deltag=x0;
						    deltaE=win*deltag;
						  }
					      }
					  }
					else if (iter>0)
					  {

					    pippo=0;

					    pluto=0;

					    if ((tot-deltag*q*Heis)>0)
					      {
						deltagmax=deltag;
						deltaEmax=deltaE;
						totmax=tot;

						x0=(deltagmin+deltagmax)/2.0;

						deltag=x0;
						deltaE=win*deltag;
						
					      }
					    else
					      {
						deltagmin=deltag;
						deltaEmin=deltaE;
						totmin=tot;

						x0=(deltagmin+deltagmax)/2.0;

						deltag=x0;
						deltaE=win*deltag;

					      }
					  
    
					    if (((fabs(deltagmax-deltagmin))/deltagmax)<(1e-4))
					      {

						slope=((totmax-deltagmax*q*Heis)-(totmin-deltagmin*q*Heis))/(deltagmax*q-deltagmin*q);

						x0=(-(totmin-deltagmin*q*Heis)/slope+deltagmin*q)/q;
						
						tot=x0*q*Heis; /* stima */

						conv=1;
					      }


					  }
				
				      }
				    else if ((i==0) & ((minAC[j][i])<=(cutoff2*0.124*1e-3)) & ((minAC[j][i])>(cutoff1*0.124*1e-3)))
				      {
					if ((iter==0) & (abs==0))
					  conv=1;
					

					if (G==0.0)
					  conv=1;


					/* puo accadere che per iter>0 non ci siano piu punti sotto la finestra. Allora lo forzo ad uscire anche se non ha raggiunto la convergenza.  */

					if ((iter>0) & (abs==0))
					  {
					    conv=1;
					    
/* 					    tot=deltag*Heis*q; */
					  }

				
					if (abs==1)
					  {
						      
					    deltaEold=deltaE;

					    deltagold=deltag;
					    
					    deltag=(tot/(Heis*q))*0.9+0.1*deltagold;
					   
					    deltaE=win*deltag; 
					
					    if (((fabs(deltag-deltagold))/deltag)<1e-2)
					      {
						conv=1;
					      }
					  }
					
					
				      }
				    else if ((i!=0) & ((minAC[j][i])<=(cutoff*0.124*1e-3)))
				      {
					if ((iter==0) & (abs==0))
					  conv=1;
					

					if (G==0.0)
					  conv=1;


					/* puo accadere che per iter>0 non ci siano piu punti sotto la finestra. Allora lo forzo ad uscire anche se non ha raggiunto la convergenza.  */

					if ((iter>0) & (abs==0))
					  {
					    conv=1;
					    
					    tot=deltag*Heis*q;
					  }

				
					if (abs==1)
					  {
						      
					    deltaEold=deltaE;

					    deltagold=deltag;
					    
					    deltag=(tot/(Heis*q))*0.9+0.1*deltagold;
					   
					    deltaE=win*deltag; 
					
					    if (((fabs(deltag-deltagold))/deltag)<1e-2)
					      {
						conv=1;
					      }
					  }
					
					
				      }
				   


				    if ((i==0) & ((minAC[j][i])<=(cutoff1*0.124*1e-3)))
				      {
					if (pippo==1)
					  iter=0; /* problema: Gamma_max-deltaEmax/(hbar*c) < 0: sono ancora alla ricerca del deltaEmax giusto */
					else
					  iter++;
				      }
				    else
				      iter++;
				    
							    
			
				  } /* end self-consistency */
				 
   
				rate[n][m]+=tot;

				if (deltag>0.3)
				  printf ("i= %d, j= %d, deltag= %lg\n",i,j,deltag);


				/* if (iter1>0) */
/* 				    printf ("iter= %d, number= %d, deltag= %lg\n",iter1,number,deltag); */
		
			   } /* end cycle on i (graphene branches) */

			 }  /* if G!=0.0 */

		       } /*  cycle on iY (allowed transverse vectors) */
	      		   
   
		 } /* kyP */

		} /* if flag==0 */

	    } /* kyE */
	

	} /* kxE */
	

/*       fclose(fp); */

      fp=fopen("rateACEM.out","w");

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
      PyObject_SetAttrString(obj,"rateAE",temp_obj);
       
      free_dmatrix(energyE,0,N-1,0,2*dimer-1);

      free_dmatrix(energyP2D,0,N-1,0,6*dimer-1);
      
      free_dmatrix(minAC,0,dimer-1,0,2);

      printf("\n");

      char *string = "Bye bye from rateACEM!";
      return Py_BuildValue("s", string);


}
