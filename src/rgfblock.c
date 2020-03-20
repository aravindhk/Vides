// ======================================================================
//  Copyright (c) 2004-2010, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 
/* ****************************************************************** */
/*          This subroutine computes the inversion of                 */
/*          a tridiagonal block matrix, through the                   */
/*          Recursive Green's Function algorithm.                     */
/*          lowdiag,diag,and updiag are the block diagonal            */
/*          of the input matrix. lowdiag[0]=0 and                     */
/*          updiag[Nc]=0. n is the order of the block matrix          */
/*          and Nc is the number of blocks. In case flagball          */
/*          is equal to 1, only the first column and the last column  */
/*          of GREEN are computed.                                    */
/* ****************************************************************** */
#include "rgfblock.h"
complex ****rgfblock(complex ***lowdiag,complex ***diag,
		   complex ***updiag,int n,int Nc,int flagball)
{
  int i,j,k,l,ii,ij,Np,ix;
  complex zero,**zerom,uno,**ID;
  complex ***gl,***gr;
  complex **temp,****out,**temp1,**temp2;
  complex ****GREEN;
  double xx;
  FILE *fp;

  //Initialization of useful variables
  zero.r=0;
  zero.i=0;
  uno.r=1;
  uno.i=0;
  zerom=cmatrix(0,n-1,0,n-1);
  ID=cmatrix(0,n-1,0,n-1);
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      {
	if (i==j)
	  ID[i][j]=uno;
	else
	  ID[i][j]=zero;
      }
  
  Np=n*Nc;
  
  if (flagball==1)
    out=ctensor4(0,Nc-1,0,1,0,n-1,0,n-1);

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      zerom[i][j]=zero;
  gl=cvectorm(0,Nc-1);
  gr=cvectorm(0,Nc-1);
  
  GREEN=ctensor4(0,Nc-1,0,Nc-1,0,n-1,0,n-1);

  //I compute the left green's functions
  gl[0]=cmatinv(diag[0],n);
  for (i=1;i<Nc;i++)
    {
      temp=cmatmul3(lowdiag[i],gl[i-1],updiag[i-1],n);
      temp1=cmatsub(diag[i],temp,n);
      gl[i]=cmatinv(temp1,n);
      cfree_cmatrix(temp,0,n-1,0,n-1);
      cfree_cmatrix(temp1,0,n-1,0,n-1);
    }
  
  //I compute the right green's functions
  gr[Nc-1]=cmatinv(diag[Nc-1],n);
  
  for (i=Nc-2;i>=0;i--)
    {
      temp1=cmatmul3(updiag[i],gr[i+1],lowdiag[i+1],n);
      temp=cmatsub(diag[i],temp1,n);
      gr[i]=cmatinv(temp,n);
      cfree_cmatrix(temp,0,n-1,0,n-1);
      cfree_cmatrix(temp1,0,n-1,0,n-1);
    }
  
/*    fp=fopen("gr","w"); */
/*    for (i=0;i<Nc;i++) */
/*      { */
/*        for (j=0;j<n;j++) */
/*  	{ */
/*  	  for (k=0;k<n;k++) */
/*  	    fprintf(fp,"%lg ",gr[i][j][k].i); */
/*  	  fprintf(fp,"\n"); */
/*  	} */
/*      } */
/*    fclose(fp); */

/*    fp=fopen("gl","w"); */
/*    for (i=0;i<Nc;i++) */
/*      { */
/*        for (j=0;j<n;j++) */
/*  	{ */
/*  	  for (k=0;k<n;k++) */
/*  	    fprintf(fp,"%lg ",gl[i][j][k].i); */
/*  	  fprintf(fp,"\n"); */
/*  	} */
/*      } */
/*    fclose(fp); */

  // I compute the diagonal elements of the
  // GREEN matrix
  
  
  GREEN[Nc-1][Nc-1]=cmatrix(0,n-1,0,n-1);
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      GREEN[Nc-1][Nc-1][i][j]=gl[Nc-1][i][j];
  
  for (i=Nc-2;i>=0;i--)
    {
/*        temp=cmatmul3(GREEN[i+1][i+1],lowdiag[i+1],gl[i],n); */
/*        temp1=cmatmul(updiag[i],temp,n); */
/*        cfree_cmatrix(temp,0,n-1,0,n-1); */
/*        temp=cmatsum(ID,temp1,n); */
/*        GREEN[i][i]=cmatmul(gl[i],temp,n); */
/*        cfree_cmatrix(temp,0,n-1,0,n-1); */
/*        cfree_cmatrix(temp1,0,n-1,0,n-1); */
      temp=cmatmul3(updiag[i],gr[i+1],lowdiag[i+1],n);
      temp1=cmatmul(gl[i],temp,n);
      cfree_cmatrix(temp,0,n-1,0,n-1);
      temp=cmatsub(ID,temp1,n);
      cfree_cmatrix(temp1,0,n-1,0,n-1);
      temp1=cmatinv(temp,n);
      GREEN[i][i]=cmatmul(temp1,gl[i],n);
      cfree_cmatrix(temp,0,n-1,0,n-1);
      cfree_cmatrix(temp1,0,n-1,0,n-1);
    }

/*    fp=fopen("Gdiag","w"); */
/*    for (i=0;i<Nc;i++) */
/*      { */
/*        for (j=0;j<n;j++) */
/*  	{ */
/*  	  for (k=0;k<n;k++) */
/*  	    fprintf(fp,"%lg ",GREEN[i][i][j][k].i); */
/*  	  fprintf(fp,"\n"); */
/*  	} */
/*      } */
/*    fclose(fp); */

  //  exit(1);
  
  
  if (flagball==1)
    { 
      //I compute the first column of the GREEN matrix
      ix=0;
      for (i=1;i<Nc;i++)
	{
	  temp=cmatmul3(gr[i],lowdiag[i],GREEN[i-1][ix],n);
	  GREEN[i][ix]=cmatrix(0,n-1,0,n-1);
	  for (k=0;k<n;k++)
	    for (l=0;l<n;l++)
	      {
		if (fabs(temp[k][l].r)>1e-30)
		  GREEN[i][ix][k][l].r=-temp[k][l].r;
		else
		  GREEN[i][ix][k][l].r=temp[k][l].r;
		if (fabs(temp[k][l].i)>1e-30)
		  GREEN[i][ix][k][l].i=-temp[k][l].i;
		else
		  GREEN[i][ix][k][l].i=temp[k][l].i;
	      }
	  cfree_cmatrix(temp,0,n-1,0,n-1);
	}
      
      //I compute the last column of the GREEN matrix
      ix=Nc-1;
      for (i=Nc-2;i>=0;i--)
	{
	  temp=cmatmul3(gl[i],updiag[i],GREEN[i+1][ix],n);
	  GREEN[i][ix]=cmatrix(0,n-1,0,n-1);
	  for (k=0;k<n;k++)
	    for (l=0;l<n;l++)
	      {
		if (fabs(temp[k][l].r)>1e-30)
		  GREEN[i][ix][k][l].r=-temp[k][l].r;
		else
		  GREEN[i][ix][k][l].r=temp[k][l].r;
		if (fabs(temp[k][l].i)>1e-30)
		  GREEN[i][ix][k][l].i=-temp[k][l].i;
		else
		  GREEN[i][ix][k][l].i=temp[k][l].i;
	      }
	  cfree_cmatrix(temp,0,n-1,0,n-1);
	}
    }
  else
    {
      //I compute the whole matrix
      for (j=0;j<Nc;j++)
	{
	  for (i=j+1;i<Nc;i++)
	    {
	      temp1=GREEN[i-1][j];
	      temp=cmatmul(lowdiag[i],temp1,n);
	      temp=cmatmul(gr[i],temp,n);
	      GREEN[i][j]=cmatrix(0,n-1,0,n-1);
	      for (k=0;k<n;k++)
		for (l=0;l<n;l++)
		  {
		    GREEN[i][j][k][l].r=-temp[k][l].r;
		    GREEN[i][j][k][l].i=-temp[k][l].i;
		  }
	      
	    }
	  for (i=j-1;i>=0;i--)
	    {
	      temp1=GREEN[i+1][j];
	      temp=cmatmul(updiag[i],temp1,n);
	      temp=cmatmul(gl[i],temp,n);
	      GREEN[i][j]=cmatrix(0,n-1,0,n-1);
	      for (k=0;k<n;k++)
		for (l=0;l<n;l++)
		  {
		    GREEN[i][j][k][l].r=-temp[k][l].r;
		    GREEN[i][j][k][l].i=-temp[k][l].i;
		  }
	    }
	}
    }
    


  if (flagball==0)
    {
      return GREEN;
    }
  else
    {
      for (i=0;i<Nc;i++)
	for (j=0;j<2;j++)
	  {
	    out[i][j]=cmatrix(0,n-1,0,n-1);
	    for (k=0;k<n;k++)
	      for (l=0;l<n;l++)
		{
		  
		  out[i][j][k][l]=GREEN[i][j*(Nc-1)][k][l];
		}	
	  }
    }
  
  //Free the memory
  cfree_cmatrix(zerom,0,n-1,0,n-1);
  cfree_cmatrix(ID,0,n-1,0,n-1);
  cfree_cvectorm(gl,0,Nc-1,0,n-1,0,n-1);
  cfree_cvectorm(gr,0,Nc-1,0,n-1,0,n-1);
  if (flagball==1) cfree_ctensor4(GREEN,0,Nc-1,0,Nc-1,0,n-1,0,n-1,flagball);
  return out;
}
