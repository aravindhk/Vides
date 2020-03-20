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
#include "rgfblock_Lake.h"
#include "cmatmul_proc.h"
complex ****rgfblock_Lake(complex ***lowdiag,complex ***diag,
			  complex ***updiag,int n,int Nc,int flagball)
{
  int i,j,k,l,ii,ij,Np,ix,jj,kk;
  complex zero,**zerom,uno,**ID;
  complex ***gl,***gr;
  complex **temp,****out,**temp1,**temp2,**temp_old;
  complex ****GREEN;
  double xx;
  FILE *fp;

  temp=cmatrix(0,n-1,0,n-1);
  temp1=cmatrix(0,n-1,0,n-1);
  temp2=cmatrix(0,n-1,0,n-1);
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
  
  GREEN=ctensor4(0,Nc-1,0,Nc-1,0,n-1,0,n-1);

  //I compute the left green's functions
  gl[0]=cmatinv(diag[0],n);
  for (i=1;i<Nc;i++)
    {
      //      temp=cmatmul3(lowdiag[i],gl[i-1],updiag[i-1],n);
      cmatmul_proc(gl[i-1],updiag[i-1],temp1,n);
      cmatmul_proc(lowdiag[i],temp1,temp,n);
      // I do the following temp2=diag[i]-temp
      //temp_old=cmatsub(diag[i],temp,n);
      for (jj=0;jj<n;jj++)
	for (kk=0;kk<n;kk++)
	  temp2[jj][kk]=csub(diag[i][jj][kk],temp[jj][kk]);
      gl[i]=cmatinv(temp2,n);
      //      cfree_cmatrix(temp_old,0,n-1,0,n-1);
      //      cfree_cmatrix(temp1,0,n-1,0,n-1);
    }
  
/*   //I compute the right green's functions */
/*   gr[Nc-1]=cmatinv(diag[Nc-1],n); */
  
/*   for (i=Nc-2;i>=0;i--) */
/*     { */
/*       temp1=cmatmul3(updiag[i],gr[i+1],lowdiag[i+1],n); */
/*       temp=cmatsub(diag[i],temp1,n); */
/*       gr[i]=cmatinv(temp,n); */
/*       cfree_cmatrix(temp,0,n-1,0,n-1); */
/*       cfree_cmatrix(temp1,0,n-1,0,n-1); */
/*     } */
  
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
      //       temp=cmatmul3(GREEN[i+1][i+1],lowdiag[i+1],gl[i],n);
       cmatmul_proc(lowdiag[i+1],gl[i],temp1,n);
       cmatmul_proc(GREEN[i+1][i+1],temp1,temp,n);
       //temp1=cmatmul(updiag[i],temp,n);
       cmatmul_proc(updiag[i],temp,temp1,n);
       //       cfree_cmatrix(temp,0,n-1,0,n-1);
       // I perform ID+temp1
       //temp_old=cmatsum(ID,temp1,n);
       for (jj=0;jj<n;jj++)
	for (kk=0;kk<n;kk++)
	  temp2[jj][kk]=csum(ID[jj][kk],temp1[jj][kk]);
       GREEN[i][i]=cmatmul(gl[i],temp2,n);
       //cfree_cmatrix(temp_old,0,n-1,0,n-1);
       //       cfree_cmatrix(temp1,0,n-1,0,n-1);

/*       temp=cmatmul3(updiag[i],gr[i+1],lowdiag[i+1],n); */
/*       temp1=cmatmul(gl[i],temp,n); */
/*       cfree_cmatrix(temp,0,n-1,0,n-1); */
/*       temp=cmatsub(ID,temp1,n); */
/*       cfree_cmatrix(temp1,0,n-1,0,n-1); */
/*       temp1=cmatinv(temp,n); */
/*       GREEN[i][i]=cmatmul(temp1,gl[i],n); */
/*       cfree_cmatrix(temp,0,n-1,0,n-1); */
/*       cfree_cmatrix(temp1,0,n-1,0,n-1); */
    }

/*    fp=fopen("Gdiag","w"); */
/*    for (i=0;i<Nc;i++) */
/*      { */
/*        for (j=0;j<n;j++) */
/*  	{ */
/*  	  for (k=0;k<n;k++) */
/* 	    //	    fprintf(fp,"%lg ",GREEN[i][i][j][k].r); */
/* 	    fprintf(fp,"%lg ",gl[i][j][k].r); */
/* 	  fprintf(fp,"\n"); */
/*  	} */
/*        fprintf(fp,"------------------\n"); */
/*      } */
/*    fclose(fp); */

/*       //I compute the first column of the GREEN matrix */
/*       ix=0; */
/*       for (i=1;i<Nc;i++) */
/* 	{ */
/* 	  temp=cmatmul3(gr[i],lowdiag[i],GREEN[i-1][ix],n); */
/* 	  GREEN[i][ix]=cmatrix(0,n-1,0,n-1); */
/* 	  for (k=0;k<n;k++) */
/* 	    for (l=0;l<n;l++) */
/* 	      { */
/* 		if (fabs(temp[k][l].r)>1e-30) */
/* 		  GREEN[i][ix][k][l].r=-temp[k][l].r; */
/* 		else */
/* 		  GREEN[i][ix][k][l].r=temp[k][l].r; */
/* 		if (fabs(temp[k][l].i)>1e-30) */
/* 		  GREEN[i][ix][k][l].i=-temp[k][l].i; */
/* 		else */
/* 		  GREEN[i][ix][k][l].i=temp[k][l].i; */
/* 	      } */
/* 	  cfree_cmatrix(temp,0,n-1,0,n-1); */
/* 	} */

  // I also allocate the memory for the first column
  for (i=1;i<Nc;i++)
    GREEN[i][0]=cmatrix(0,n-1,0,n-1);
      
  //I compute the last column of the GREEN matrix
  //and return back the needed elements of G in out
  ix=Nc-1;
  for (i=Nc-2;i>=0;i--)
    {
      //      temp=cmatmul3(gl[i],updiag[i],GREEN[i+1][ix],n);
      cmatmul_proc(updiag[i],GREEN[i+1][ix],temp1,n);
      cmatmul_proc(gl[i],temp1,temp,n);
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
      //cfree_cmatrix(temp,0,n-1,0,n-1);
      
    }

  for (i=0;i<Nc;i++)
    {
      out[i][0]=cmatrix(0,n-1,0,n-1);
      out[i][1]=cmatrix(0,n-1,0,n-1);
      for (k=0;k<n;k++)
	for (l=0;l<n;l++)
	  {
	    out[i][1][k][l]=GREEN[i][Nc-1][k][l];
	    out[i][0][k][l]=GREEN[i][i][k][l];
	  }
    }
  
  //Free the memory
  cfree_cmatrix(temp,0,n-1,0,n-1);
  cfree_cmatrix(temp1,0,n-1,0,n-1);
  cfree_cmatrix(temp2,0,n-1,0,n-1);
  cfree_cmatrix(zerom,0,n-1,0,n-1);
  cfree_cmatrix(ID,0,n-1,0,n-1);
  cfree_cvectorm(gl,0,Nc-1,0,n-1,0,n-1);
  if (flagball==1) cfree_ctensor4(GREEN,0,Nc-1,0,Nc-1,0,n-1,0,n-1,flagball);
  return out;
}
