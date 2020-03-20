// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 


//
// This subroutine computes the charge within the zig-zag (n,0) nanotube, 
// as well as the transmission coefficient, within the mode space approach.
//
// N.B.: Atoms are ordered first within the ring and then along the 
// CNT axis. E.G.: if i runs from 0 to n-1 (where n is the number of atoms
// along the ring) and j from 0 to Nc-1 (where Nc is the number of rings)
// the ix-th atom with index (i,j) has index ix=i+j*n
//
// INPUT:
//
// nanotube.n: number of atoms along CNT ring
// nanotube.Nc: number of CNT rings
// nanotube.Eupper [eV]: highest considered energy
// nanotube.Elower [eV]: lowest considered energy
// nanotube.dE [eV]: energy step
// nanotube.thop [eV]: hopping parameter
// nanotube.eta [eV]: broadening parameter
// nanotube.Temp [K]: CNT temperature
// nanotube.mu1 [eV]: Fermi level of the left reservoir
// nanotube.mu2 [eV]: Fermi level of the right reservoir
// nanotube.Nmodes: Number of considered modes
// nanotube.contact [char]: specifies if Schottky or doped reservoirs 
//                          have to be considered
// nanotube.Phi (vector n*Nc) [eV]: electrostatic potential in correspondence 
//                                  of each C atom ordered as above.
//
// OUTPUT:
//
// nanotube.charge (vector n*Nc): charge in correspondence of each atoms 
//                                (electrons-holes), i.e. if positive,
//                                then electrons, if negative then holes.
//                                Atoms are ordered as above.
// nanotube.E (vector) [eV]: computed energy levels stored in a vector
// nanotube.T (vector): nanotube.T[i] is the transmission coefficient 
//                      computed at nanotube.E[i].
//
//

#include "CNTmode_charge_T.h"
static PyObject* py_CNTmode_charge_T(PyObject* self, PyObject* args)
{
  int i,j,k,ix,ie,*order,Nm,icse,NE,n,Nc,nx,rank;
  complex ***DIAG,***UPDIAG,***LOWDIAG,***DIAGMODE;
  complex zero,**SIGMAS,**SIGMAD;
  double E,dE,Emin,Emax,*tempo,eta,
    **A1,**A2,*ncarcnt,T,
    sum,Eupper,Elower,*Phi,thop,*EE,*TE,mu1,mu2,vt,Egap,Temperature;
  PyObject *obj,*temp_obj;
  PyArrayObject *Phi_array;
  char *CNTboundary;

  setvbuf(stdout,(char *) NULL,_IONBF,1);

  // I read the data from python class
  import_array();
  if (!PyArg_ParseTuple(args,"O",&obj))
    {
      printf("NOT A CLASS! \n");
      exit(0);
    }

  temp_obj=PyObject_GetAttrString(obj,"rank");
  rank=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("*******************************************************\n");
  if (!rank) printf("*** Computing NEGF in CNT within the MODE SPACE *******\n"); 
  if (!rank) printf("*******************************************************\n\n");
  temp_obj=PyObject_GetAttrString(obj,"n");
  n=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("n of the nanotube = %d \n",n);
  temp_obj=PyObject_GetAttrString(obj,"Nc");
  Nc=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Number of rings of the nanotube = %d \n",Nc);
  temp_obj=PyObject_GetAttrString(obj,"Eupper");
  Eupper=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Upper Energy defined by the user = %lg \n",Eupper);
  temp_obj=PyObject_GetAttrString(obj,"Elower");
  Elower=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Lower Energy defined by the user = %lg \n",Elower);
  temp_obj=PyObject_GetAttrString(obj,"dE");
  dE=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Delta Energy = %lg \n",dE);
  temp_obj=PyObject_GetAttrString(obj,"thop");
  thop=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Hopping parameter = %lg \n",thop);
  temp_obj=PyObject_GetAttrString(obj,"eta");
  eta=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Eta = %lg \n",eta);
  temp_obj=PyObject_GetAttrString(obj,"Temp");
  Temperature=(double)PyFloat_AsDouble(temp_obj);
  vt=kboltz*Temperature/q;
  if (!rank) printf("Temperature = %lg \n",Temperature);
  temp_obj=PyObject_GetAttrString(obj,"mu1");
  mu1=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Electrochemical potential of the source = %lg \n",mu1);
  temp_obj=PyObject_GetAttrString(obj,"mu2");
  mu2=(double)PyFloat_AsDouble(temp_obj);
  if (!rank) printf("Electrochemical potential of the drain = %lg \n",mu2);
  temp_obj=PyObject_GetAttrString(obj,"contact");
  #if PY_MAJOR_VERSION >=3
  CNTboundary=PyBytes_AsString(PyUnicode_AsEncodedString(temp_obj,"utf-8","Error"));
  #else
  CNTboundary=PyString_AsString(temp_obj);
  #endif

  if (strcasecmp(CNTboundary,"doped")==0)
    if (!rank) printf("CNT contacts : doped reservoirs \n");
  else if (strcasecmp(CNTboundary,"schottky")==0)
    if (!rank) printf("CNT contacts : Schottky contacts \n");
  else
    {
      if (!rank) printf("Error in the definition of the CNT contact \n");
      if (!rank) printf("It can be either doped or schottky \n");
      exit(0);
    }
  temp_obj=PyObject_GetAttrString(obj,"Nmodes");
  Nm=(int)PyInt_AsLong(temp_obj);
  if (!rank) printf("Number of modes = %d \n",Nm);

  Phi=dvector(0,n*Nc-1);
  Phi_array=(PyArrayObject *)PyObject_GetAttrString(obj,"Phi");
  nx=Phi_array->dimensions[0];
  for (i=0;i<nx;i++)
    {
      tempo=(double *)(Phi_array->data+i*Phi_array->strides[0]);
      Phi[i]=*tempo;
    }
  
  Egap=fabs(2*thop*pi/(sqrt(3)*n));
  if (!rank) printf("Egap of the nanotube = %lg \n",Egap);
  if (!rank) printf("****************************** \n");

  // This is to order the eigenmodes to compute the 
  // analytical self-energy
  order=ivector(0,n-1);
  energysort(&order,n,thop);

  // Variable initialization and creation
  ie=0;
  zero.r=0;
  zero.i=0;
  

  Emin=300;
  Emax=-10;
  
  DIAGMODE=cvectorm(0,Nc-1);
  DIAG=cvectorm(0,Nc-1);
  LOWDIAG=cvectorm(0,Nc-1);
  UPDIAG=cvectorm(0,Nc-1);
  ncarcnt=dvector(0,n*Nc-1);
  for (i=0;i<n*Nc;i++)
    {
      ncarcnt[i]=0;
    }
  
  // Now I compute the energy range in which NEGF 
  // is solved if not specified by the user

  if ((Elower< -900)||(Eupper>900))

    for (i=0;i<n*Nc;i++)
      {
	if (Emax<=(-Phi[i]))
	  Emax=(-Phi[i]);
	if (Emin>=(-Phi[i]))
	  Emin=(-Phi[i]);
      }
  

  Emin=min(Emin,(mu2));
  Emin=min((mu1),Emin);
  Emax=max(Emax,(mu2));
  Emax=max((mu1),Emax);
  if (Eupper>900)
    Eupper=Emax+(Egap)*0.5+7*(vt);
  if (Elower<-900)
    Elower=Emin-(Egap)*0.5-7*(vt);

  if (!rank) printf("Eupper = %lg \n",Eupper);
  if (!rank) printf("Elower = %lg \n",Elower);
  if (!rank) printf("Current computed energy \n ");
  
  // NE is an extimation of the number of energy levels
  NE=(int)((Eupper-Elower)/dE);
  // ..... just to be sure to overestimate
  NE=NE+10;
  // Once known the number of energy levels I create the energy level and
  // the transmission coefficient matrix
  EE=dvector(0,NE-1);
  TE=dvector(0,NE-1);
  
  // I start to define the Hamiltonian for the carbon nanotube, changed by sign
  
  // I create the Diagonal block changed by sign
  for (i=0;i<Nc;i++)
    {
      DIAG[i]=cmatrix(0,n-1,0,n-1);
      for (j=0;j<n;j++)
	for (k=0;k<n;k++)
	  DIAG[i][j][k]=zero;
      
      for (j=0;j<n;j++)
	DIAG[i][j][j]=complass(-(-Phi[j+i*n]),0);
    }
  
  for (i=0;i<Nc;i++)
    DIAGMODE[i]=VdagaAV(DIAG[i],n,Nm,order);

  // I create the Updiagonal block changed by sign
  for (i=0;i<(Nc-1);i++)
    UPDIAG[i]=create_updiagmode(i,n,Nm,order,thop);

  
  UPDIAG[Nc-1]=cmatrix(0,Nm-1,0,Nm-1);
  for (i=0;i<Nm;i++)
    for (j=0;j<Nm;j++)
      UPDIAG[Nc-1][i][j]=zero;
  
  
  // I create the Lowdiag block changed by sign
  for (i=1;i<Nc;i++)
    LOWDIAG[i]=create_lowdiagmode(i-1,n,Nm,order,thop);
  
  LOWDIAG[0]=cmatrix(0,Nm-1,0,Nm-1);
  for (i=0;i<Nm;i++)
    for (j=0;j<Nm;j++)
      LOWDIAG[0][i][j]=zero;

  E=Elower;

  while (E<=(Eupper+dE*0.5))
    {
      
 
      
      // I create the Self-Energy in case of Schottky barrier contacts
      // or doped contacts

      if (strcasecmp(CNTboundary,"doped")==0)
	{
	  SIGMAS=selfanaliticalmode(E,-Phi[0],n,Nm,order,thop,eta);
	  SIGMAD=selfanaliticalmode(E,-Phi[(Nc-1)*n],n,Nm,order,thop,eta);
	}
      else
	{
	  SIGMAS=selfschottky(E,(mu1),Nm,thop,eta);
	  SIGMAD=selfschottky(E,(mu2),Nm,thop,eta);
	}


      // I print out the energy level for which 
      // the computation is performed every 50 energy levels
      if ((ie%50)==0)
	{      
	  if (!rank) printf("\r");
	  if (!rank) printf("%lg ",E);
	}
      
      LDOSMODE(E,LOWDIAG,DIAGMODE,UPDIAG,&A1,&A2,SIGMAS,SIGMAD,Nm,Nc,
	       1,&T,n,order,thop,eta);

      
      cfree_cmatrix(SIGMAS,0,Nm-1,0,Nm-1);
      cfree_cmatrix(SIGMAD,0,Nm-1,0,Nm-1);
      
      // I compute the charge
      for (i=0;i<Nc;i++)
	for (j=0;j<n;j++)
	  {
	    // electrons
	    ix=j+i*n;
	    if ((E>=-Phi[j+i*n]))
	      {
		ncarcnt[ix]=ncarcnt[ix]
		  -2*(
		      A1[i][j]*Fermi_Dirac((E-(mu1))/(vt))*dE);
	      }
	    
	    if ((E>=-Phi[j+i*n]))
	      {
		ncarcnt[ix]=ncarcnt[ix]
		  -2*(
		      A2[i][j]*Fermi_Dirac((E-(mu2))/(vt))*dE);
	      }
	    
	    //holes
	    if ((E<-Phi[j+i*n]))
	      {
		ncarcnt[ix]=ncarcnt[ix]
		  +2*(
		      A1[i][j]*(1-Fermi_Dirac((E-(mu1))/(vt)))*dE);
	      }
	    
	    if ((E<-Phi[j+i*n]))
	      {
		ncarcnt[ix]=ncarcnt[ix]
		  +2*(
		      A2[i][j]*(1-Fermi_Dirac((E-(mu2))/(vt)))*dE);
	      }
	    // last modifications
	    
	  }
      
      EE[ie]=E;
      TE[ie]=T;

      E+=dE;
      ie++;
      free_dmatrix(A1,0,Nc-1,0,n-1);
      free_dmatrix(A2,0,Nc-1,0,n-1);

    }

  if (!rank) printf("\n\n*******************************************************\n");
  if (!rank) printf("********* END NEGF in CNT within the MODE SPACE *******\n"); 
  if (!rank) printf("*******************************************************\n\n");

  // I set the attribute of the class to send data out
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)EE);
  PyObject_SetAttrString(obj,"E",temp_obj);
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)TE);
  PyObject_SetAttrString(obj,"T",temp_obj);
  ie=Nc*n;
  temp_obj=  PyArray_FromDimsAndData(1,&ie,PyArray_DOUBLE,(char*)ncarcnt);
  PyObject_SetAttrString(obj,"charge",temp_obj);
  
  // Free the memory
  free_ivector(order,0,n-1);
  cfree_cvectorm(DIAG,0,Nc-1,0,n-1,0,n-1);
  cfree_cvectorm(LOWDIAG,0,Nc-1,0,n-1,0,n-1);
  cfree_cvectorm(UPDIAG,0,Nc-1,0,n-1,0,n-1);
  //free_dvector(EE,0,NE-1);
  //free_dvector(TE,0,NE-1);
  char *s = "Hello from CNTmode_charge_T";
  return Py_BuildValue("s", s);
}
