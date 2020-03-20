// ======================================================================
//  Copyright (c) 2011, G. Fiori  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

// Calculates the diagonal blocks in the form (E+i*eta)I-U, 
// where E is the energy and eta a small imaginary part.
// Result: lead self energy.

#include "selfH_new.h"
complex **selfH_new(double E, complex ***diag, complex ***updiag, complex ***lowdiag, int N, int Nc, int lead, double eta)
{
  complex **Gsnew, **Gz, **wmH, zero, **t, **selfenergy, **tdaga,**P,**temp,**temp2,**temp_min;
  complex *eigenvalues,**w01,**w23,*ld;
  complex **temp3,**temp4,**temp5,cc,cc2,**MM,*b,*v,**LAMBDA;
  int rank,rankR,rankC;
  double sum;

  double kvect[2000],num,den;
  int i, j, k,kk,*flag_right,ix,jj,ncomp,*index_zero_row,*index_zero_column,incr;

  zero.r = 0.0;
  zero.i = 0.0;

  wmH = cmatrix(0,4*N-1,0,4*N-1);
  P = cmatrix(0,4*N-1,0,4*N-1);
  t = cmatrix(0,N-1,0,N-1);
  //  tdaga = cmatrix(0,N-1,0,N-1);
  selfenergy = cmatrix(0,N-1,0,N-1);

  for(i=0; i<4*N; i++){
    for(j=0; j<4*N; j++){
      wmH[i][j]=zero;
      P[i][j]=zero;
    } 
  }


  //
  // This part is for the test with the python script script_GNR.py and self_W.py
  //

/*   // I take care of the diagonal block matrix */
/*   if (lead==0) */
/*     { */
/*       for(k=0; k<4; k++) */
/* 	for(i=0; i<N; i++) */
/* 	  for(j=0; j<N; j++) */
/* 	    if (i==j) */
/* 	      wmH[j+k*N][j+k*N]=complass(E+diag[k][j][j].r,+diag[k][j][j].i+eta); */
/* 	    else */
/* 	      wmH[i+k*N][j+k*N]=complass(diag[k][i][j].r,diag[k][i][j].i); */
/*     } */
/*   else */
/*     { */
/*       for(k=0; k<4; k++) */
/* 	for(i=0; i<N; i++) */
/* 	  for(j=0; j<N; j++) */
/* 	    if (i==j) */
/* 	      wmH[j+k*N][j+k*N]=complass(E+diag[Nc-4+k][j][j].r,+diag[Nc-4+k][j][j].i+eta); */
/* 	    else */
/* 	      wmH[i+k*N][j+k*N]=complass(diag[Nc-4+k][i][j].r,diag[Nc-4+k][i][j].i); */
/*     } */


/*   // I take care of the off-diagonal block matrix */
/*   for(k=0; k<3; k++){ */
/*     for(i=0; i<(1)*N; i++){ */
/*       for(j=0; j<N; j++){ */
/* 	wmH[i+k*N][j+(k+1)*N]=complass(+updiag[Nc-4+k][i][j].r,updiag[Nc-4+k][i][j].i); */
/*       } */
/*     } */
/*   } */
  
/*   for(k=1; k<4; k++){ */
/*     for(i=0; i<N; i++){ */
/*       for(j=0; j<N; j++){ */
/* 	wmH[i+k*N][j+(k-1)*N]=complass(+lowdiag[Nc-4+k][i][j].r,lowdiag[Nc-4+k][i][j].i); */
/*       } */
/*     } */
/*   } */


/*   FILE *fp; */
/*   fp=fopen("wmh.Re","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",wmH[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("wmh.Im","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",wmH[i][j].i); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("P.Re","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",P[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("P.Im","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",P[i][j].i); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   exit(0); */





  // I take care of the diagonal block matrix
/*   if (lead==0) */
/*     { */
/*       for(k=0; k<3; k++) */
/* 	for(i=0; i<N; i++) */
/* 	  for(j=0; j<N; j++) */
/* 	    if (i==j) */
/* 	      wmH[j+k*N][j+k*N]=complass(E+diag[k][j][j].r,+diag[k][j][j].i+eta); */
/* 	    else */
/* 	      wmH[i+k*N][j+k*N]=complass(diag[k][i][j].r,diag[k][i][j].i); */
/*     } */
/*   else */


  // I now create the wmH and P matrices for the left and right contact
  // left  => lead=1
  // right => lead=0

  if (lead==1)
    {
      {
	for(k=0; k<3; k++)
	  for(i=0; i<N; i++)
	    for(j=0; j<N; j++)
	      if (i==j)
		{
		  wmH[j+k*N][j+k*N]=complass(E+diag[Nc-4+k][j][j].r,+diag[Nc-4+k][j][j].i+eta);
		  /* if ((wmH[j+k*N][j+k*N].r<1e-20)&&(wmH[j+k*N][j+k*N].i<1e-20)) */
/* 		    wmH[j+k*N][j+k*N]=complass(2,0); */
		}
	
	      else
		wmH[i+k*N][j+k*N]=complass(diag[Nc-4+k][i][j].r,diag[Nc-4+k][i][j].i);
      }
      
      
      // I take care of the off-diagonal block matrix
      for(k=0; k<3; k++)
	{
	  for(i=0; i<(1)*N; i++)
	    {
	      for(j=0; j<N; j++)
		{
		  wmH[i+k*N][j+(k+1)*N]=complass(+updiag[Nc-4+k][i][j].r,updiag[Nc-4+k][i][j].i);
		}
	    }
	}
      
      for(k=1; k<3; k++){
	for(i=0; i<N; i++){
	  for(j=0; j<N; j++){
	    wmH[i+k*N][j+(k-1)*N]=complass(+lowdiag[Nc-4+k][i][j].r,lowdiag[Nc-4+k][i][j].i);
	  }
	}
      }
      
      // I take care of the T-10 element of wmH
      for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	  wmH[i+3*N][j]=complass(+updiag[Nc-5][i][j].r,updiag[Nc-5][i][j].i);
	}
      }
      
      // I take care of the P matrix
      for(i=0; i<N; i++)
	{
	  for(j=0; j<N; j++)
	    {
	      P[i][j+3*N]=complass(-lowdiag[Nc-4][i][j].r,-lowdiag[Nc-4][i][j].i);
	      P[i+3*N][j+2*N]=complass(-lowdiag[Nc-1][i][j].r,-lowdiag[Nc-1][i][j].i);
	      if (i==j)
		{
		  P[i+3*N][j+3*N]=complass(-(E+diag[Nc-1][i][j].r),-(+diag[Nc-1][i][j].i+eta));
		  /* if ((P[i+3*N][j+3*N].r<1e-20)&&(P[i+3*N][j+3*N].i<1e-20)) */
/* 		    P[i+3*N][j+3*N]=complass(2,0); */
		}
	      else
		P[i+3*N][j+3*N]=complass(-(diag[Nc-1][i][j].r),-(+diag[Nc-1][i][j].i));
	    }
	}
    }
  else
    {

      // I First compute the full wmH matrix

      // I take care of the diagonal block matrix
      for(k=0; k<4; k++)
	for(i=0; i<N; i++)
	  for(j=0; j<N; j++)
	    if (i==j)
	      wmH[j+k*N][j+k*N]=complass(E+diag[k][j][j].r,+diag[k][j][j].i+eta);
	    else
	      wmH[i+k*N][j+k*N]=complass(diag[k][i][j].r,diag[k][i][j].i);

      // I take care of the off-diagonal block matrix
      for(k=0; k<3; k++){
	for(i=0; i<(1)*N; i++){
	  for(j=0; j<N; j++){
	    wmH[i+k*N][j+(k+1)*N]=complass(+updiag[k][i][j].r,updiag[k][i][j].i);
	  } 
	}
      } 
      
      for(k=1; k<4; k++){
	for(i=0; i<N; i++){
	  for(j=0; j<N; j++){
	    wmH[i+k*N][j+(k-1)*N]=complass(+lowdiag[k][i][j].r,lowdiag[k][i][j].i);
	  } 
	}
      }    

      // I include the periodic part

      for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	  wmH[i+3*N][j]=complass(+updiag[3][i][j].r,updiag[3][i][j].i);
	  wmH[i][j+3*N]=complass(lowdiag[4][i][j].r,lowdiag[4][i][j].i);
	}
      }
      
      // I now flip the matrix
      flip_cmatrix(wmH,4*N);
      
      // I now create the H and P matrices 
      
      // I now work on the P matrix
      for (i=0;i<N;i++)
	for (j=0;j<N;j++)
	  {
	    P[i][j+3*N]=complass(-wmH[i][j+3*N].r,-wmH[i][j+3*N].i);
	    P[3*N+i][2*N+j]=complass(-wmH[3*N+i][2*N+j].r,-wmH[3*N+i][2*N+j].i);
	    P[3*N+i][3*N+j]=complass(-wmH[3*N+i][3*N+j].r,-wmH[3*N+i][3*N+j].i);
	  }


      // I create H from wmH, putting to zero the not-needed block matrices
      for (i=0;i<2*N;i++)
	for (j=0;j<N;j++)
	  {
	    wmH[i][3*N+j]=zero;
	    wmH[3*N+j][2*N+i]=zero;
	  }


    }


/*   // Now I check how many lines of T32 and D33  */
/*   rank=N; */
/*   index_zero_row=ivector(0,N-1); */
/*   for (i=0;i<N;i++) */
/*     { */
/*       sum=0; */
/*       for (j=0;j<4*N;j++) */
/* 	if (i==(j-3*N)) */
/* 	  { */
/* 	    sum=sum+cdabs(csub(P[i+3*N][j],complass(-E,-(2*lead-1)*eta))); */
/* 	  } */
/* 	else */
/* 	  sum=sum+cdabs(P[i+3*N][j]); */

/*       if (sum<=eta) */
/* 	{ */
/* 	  rank--; */
/* 	  index_zero_row[i]=1; */
/* 	} */
/*       else */
/* 	index_zero_row[i]=0; */
/*     } */
    
/*   printf("RANK %d \n",rank); */
    
/*   FILE *fp; */
/*   fp=fopen("wm.Re","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",wmH[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */
  
/*   fp=fopen("wm.Im","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",wmH[i][j].i); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */
  
/*   fp=fopen("P.Re","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",P[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */
  
/*   fp=fopen("P.Im","w"); */
/*   for (i=0;i<4*N;i++) */
/*     { */
/*       for (j=0;j<4*N;j++) */
/* 	fprintf(fp,"%lg ",P[i][j].i); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */
/*   exit(0); */
  
  
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      t[i][j]=wmH[3*N+i][j];


/*   for (i=0;i<N;i++) */
/*     for (j=0;j<N;j++) */
/*       { */
/* 	tdaga[i][j].r=t[j][i].r; */
/* 	tdaga[i][j].i=-t[j][i].i; */
/*       } */


//  FILE *fp;
/*   fp=fopen("t.Re","w"); */
/*   for (i=0; i<N; i++) */
/*     { */
/*       for (j=0; j<N; j++) */
/* 	{ */
/* 	  fprintf(fp,"%lg ",t[i][j].r); */
/* 	} */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("t.Im","w"); */
/*   for (i=0; i<N; i++) */
/*     { */
/*       for (j=0; j<N; j++) */
/* 	{ */
/* 	  fprintf(fp,"%lg ",t[i][j].i); */
/* 	} */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

  //   exit(0);

  temp=cmatsub(wmH,P,4*N);
  temp2=cmatinv(temp,4*N);
  cfree_cmatrix(temp,0,4*N-1,0,4*N-1);
  temp=cmatmul(temp2,P,4*N);
  cfree_cmatrix(temp2,0,4*N-1,0,4*N-1);
  temp_min=cmatrix(0,2*N-1,0,2*N-1);
  LAMBDA=cmatrix(0,2*N-1,0,2*N-1);
  for (i=0;i<2*N;i++)
    for (j=0;j<2*N;j++)
      {
	temp_min[i][j]=temp[i+2*N][j+2*N];
	LAMBDA[i][j]=zero;
      }


/*   FILE *fp; */
/*   fp=fopen("tmin.re","w"); */
/*   for (i=0; i<2*N; i++) */
/*     { */
/*       for (j=0; j<2*N; j++) */
/* 	{ */
/* 	  fprintf(fp,"%lg ",temp_min[i][j].r); */
/* 	} */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */


/*   fp=fopen("tmin.im","w"); */
/*   for (i=0; i<2*N; i++) */
/*     { */
/*       for (j=0; j<2*N; j++) */
/* 	{ */
/* 	  fprintf(fp,"%lg ",temp_min[i][j].i); */
/* 	} */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */
/*   //  exit(0); */


/*   printf("-------------------------\n"); */
/*   exit(0); */

  eigenvalues=ccvector(0,2*N-1);
  w23=cmatrix(0,2*N-1,0,2*N-1);

  eigenvalues_non_symmetric_matrix(temp_min,eigenvalues,w23,2*N,2*N,1);


  //  eigenvalues_non_symmetric_matrix_schur(temp_min,eigenvalues,w23,2*N,lead);

/*   FILE *fp; */
/*   fp=fopen("eigL.Re","w"); */
/*   for (i=0;i<2*N;i++) */
/*     fprintf(fp,"%lg \n",eigenvalues[i].r); */
/*   fclose(fp); */

/*   fp=fopen("eigL.Im","w"); */
/*   for (i=0;i<2*N;i++) */
/*     fprintf(fp,"%lg \n",eigenvalues[i].i); */
/*   fclose(fp); */
/*   exit(0); */


  // I now compute the eigenvalues lambda (ld) and the eigenvectors psi0 and psi1
  // and the left and right moving electrons (flag_right)
  //w01=cmatrix(0,2*N-1,0,2*N-1);
  ld=ccvector(0,2*N-1);
  flag_right=ivector(0,2*N-1);

  for (i=0;i<2*N;i++)
    for (j=0;j<2*N;j++)
      temp_min[i][j]=temp[i][j+2*N];

/*   for (i=0; i<2*N; i++) */
/*     { */
/*       for (j=0; j<2*N; j++) */
/* 	{ */
/* 	  printf("%lg ",temp_min[i][j].i); */
/* 	} */
/*       printf("\n"); */
/*     } */
/*   exit(0); */
  
/*   FILE *fp; */
/*   fp=fopen("w23.Re","w"); */
/*   for (i=0;i<2*N;i++) */
/*     { */
/*       for (j=0;j<2*N;j++) */
/* 	fprintf(fp,"%lg ",w23[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("w23.Im","w"); */
/*   for (i=0;i<2*N;i++) */
/*     { */
/*       for (j=0;j<2*N;j++) */
/* 	fprintf(fp,"%lg ",w23[i][j].i); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   fclose(fp); */

/*   fp=fopen("eigLabs","w"); */
/*   for (i=0;i<2*N;i++) */
/*     fprintf(fp,"%lg \n",cdabs(eigenvalues[i])); */
/*   fclose(fp); */

/*   exit(0); */


/*   printf("---------------\n"); */

  for (i=0;i<2*N;i++)
    {
      //      cc=eigenvalues[i];
      //      num=-cIm(cc)/cdabs(cc)/cdabs(cc);
      //      den=(cdabs(cc)*cdabs(cc)+cRe(cc))/(cdabs(cc)*cdabs(cc));
      //      kvect[i]=-atan(num/den);
      den=cdabs(complass(1+eigenvalues[i].r,eigenvalues[i].i))*cdabs(complass(1+eigenvalues[i].r,eigenvalues[i].i));
      ld[i]=complass((eigenvalues[i].r+cdabs(eigenvalues[i])*cdabs(eigenvalues[i]))/den,eigenvalues[i].i/den);
      
      //printf("%lg+j%lg \n",ld[i].r,ld[i].i);

      // note that |ld|<1 for the right contact means left travelling electrons

      if (cdabs(ld[i])<1)
	flag_right[i]=1;
      else
	flag_right[i]=0;

/*       for (kk=0;kk<2*N;kk++) */
/* 	{ */
/* 	  cc=zero; */
/* 	  for (j=0;j<2*N;j++) */
/* 	    cc=csum(cc,cmul(temp_min[kk][j],w23[j][i])); */
/* 	  cc2=complass(-eigenvalues[i].r/(cdabs(eigenvalues[i])*cdabs(eigenvalues[i])), */
/* 		       eigenvalues[i].i/(cdabs(eigenvalues[i])*cdabs(eigenvalues[i]))); */
/* 	  w01[kk][i]=cmul(cc2,cc); */
/* 	} */
      cc2=complass(-eigenvalues[i].r/(cdabs(eigenvalues[i])*cdabs(eigenvalues[i])),
		   eigenvalues[i].i/(cdabs(eigenvalues[i])*cdabs(eigenvalues[i])));
      LAMBDA[i][i]=cc2;
    }
  w01=cmatmul3(temp_min,w23,LAMBDA,2*N);




/*   printf("---------------\n"); */
/*   for (i=0;i<2*N;i++) */
/*     { */
/*       for (kk=0;kk<2*N;kk++) */
/* 	printf("%lg ",w01[kk][i].r); */
/*       printf("\n"); */
/*     } */
/*   printf("---------------\n"); */
/*   exit(0); */

/*   FILE *fp; */
/*   fp=fopen("eigL.Re","w"); */
/*   for (i=0;i<2*N;i++) */
/*     fprintf(fp,"%lg \n",ld[i].r); */
/*   fclose(fp); */

/*   fp=fopen("eigL.Im","w"); */
/*   for (i=0;i<2*N;i++) */
/*     fprintf(fp,"%lg \n",ld[i].i); */
/*   fclose(fp); */
    

/*   for (i=0;i<2*N;i++)  */
/*     printf("%d ",flag_right[i]); */
/*   printf("\n"); */



  // I find now the g00*T10
  MM=cmatrix(0,N-1,0,N-1);
/*   for (i=0;i<N;i++) */
/*     for (j=0;j<N;j++) */
/*       MM[i][j]=zero; */
  b=ccvector(0,N-1);
  cfree_cmatrix(temp,0,4*N-1,0,4*N-1);
  temp2=cmatrix(0,N-1,0,N-1);
  
  k=0;
  // I compute the MM matrix
  for (i=0;i<2*N;i++)
    {
      // left moving electrons
      if (flag_right[i]==1)
	{
	  for (jj=0;jj<N;jj++)
	    {
	      MM[k][jj]=w23[jj+N][i];
	    }
	  k=k+1;
	}
    }
  
  // Now I check the rank of MM, which can be smaller than N, due to the zero-padding.
  // If there is zero-padding, then I fill the N-rank
  // columns with the canonical vectors (since every vectors
  // is a solution)


  rankR=N;
  rankC=N;
  index_zero_row=ivector(0,N-1);
  index_zero_column=ivector(0,N-1);
  for (i=0;i<N;i++)
    {
      // check the rows
      sum=0;
      for (j=0;j<N;j++)
	sum=sum+cdabs(MM[i][j]);

      if (sum<=eta)
	{
	  rankR--;
	  index_zero_row[i]=1;
	}
      else
	index_zero_row[i]=0;

      // check the columns
      sum=0;
      for (j=0;j<N;j++)
	sum=sum+cdabs(MM[j][i]);

      if (sum<=eta)
	{
	  rankC--;
	  index_zero_column[i]=1;
	}
      else
	index_zero_column[i]=0;

    }
    
  if (rankR<=rankC)
    rank=rankR;
  else 
    rank=rankC;
  //  printf("RANK %d %d \n",rank,rankR);

/*   for (i=0;i<N;i++) */
/*     printf("%d %d \n",index_zero_row[i],index_zero_column[i]); */
/*   //  exit(0); */

  // I now fill the MM matrix with the canonical vectors
  incr=0;
  if (rank<N)
    if (rankR==N)
      {
	for (j=0;j<N;j++)
	  if (index_zero_column[j]==1)
	    {
	      for (i=0;i<N-rankC;i++)
		if (i==(j-rankC))
		  MM[i+rankC][j]=complass(1,0);
		else
		  MM[i+rankC][j]=zero;
	    }
      }
    else
      {
	for (j=0;j<N;j++)
	  if (index_zero_column[j]==1)
	    {
	      while (index_zero_row[incr]!=1)
		incr++;
	      for (i=0;i<N;i++)
		if (i==incr)
		  MM[i][j]=complass(1,0);
		else
		  MM[i][j]=zero;
	      incr++;
	    }
      }

/*     if (lead==1) */
/*       { */
/* 	for (i=0;i<N-rank;i++) */
/* 	  for (j=0;j<N-rank;j++) */
/* 	    if (i==j) */
/* 	      MM[rank+i][rank+j]=complass(1,0); */
/* 	    else */
/* 	      MM[rank+i][rank+j]=zero; */
/*       } */
/*     else */
/*       { */
/* 	for (i=0;i<N-rank;i++) */
/* 	  for (j=0;j<N-rank;j++) */
/* 	    if (i==j) */
/* 	      MM[rank+i][j]=complass(1,0); */
/* 	    else */
/* 	      MM[rank+i][j]=zero; */
/*       } */

/*   FILE *fp; */
/*   fp=fopen("MM","w"); */
   
/*   for (i=0;i<N;i++) */
/*     { */
/*       for (j=0;j<N;j++) */
/* 	fprintf(fp,"%lg ",MM[i][j].r); */
/*       fprintf(fp,"\n"); */
/*     } */
/*   exit(0); */

  temp=cmatinv(MM,N);

  for (ncomp=0;ncomp<N;ncomp++)
    {
      kk=0;
      for (i=0;i<2*N;i++)
	{
	  // left moving electrons
	  if (flag_right[i]==1)
	    {
	      b[kk]=cmul(ld[i],w01[ncomp][i]);
	      kk=kk+1;
	    }
	}
      
/*       printf("---------------\n"); */
/*       for (i=0;i<N;i++) */
/* 	{ */
/* 	  //for (j=0;j<N;j++) */
/* 	  //	    printf("%lg ",MM[i][j].r); */
/* 	  //printf("\n"); */
/* 	  printf("%lg+j%lg \n",b[i].r,b[i].i); */
/* 	} */
/*       /\*   exit(0); *\/ */

      v=cmatvect(temp,b,N);
      for (i=0;i<N;i++)
	{
	  temp2[ncomp][i]=v[i];
	  //	printf("%lg \n",v[ix]);
	}
      //  exit(0);
      cfree_ccvector(v,0,N-1);


    }
      
/*   for (i=0;i<N;i++) */
/*     { */
/*       for (j=0;j<N;j++) */
/* 	printf("%lg ",temp2[i][j].r); */
/*       printf("\n"); */
/*     } */
      
  //   exit(0);


  selfenergy=cmatmul(t,temp2,N);


/*   for (i=0;i<N;i++) */
/*     { */
/*       //      printf("%lg \n",b[i].i); */
/*       for (j=0;j<N;j++) */
/*       	printf("%lg+j%lg ",selfenergy[i][j].r,selfenergy[i][j].i); */
/*       printf("\n"); */
/*     } */


  if (lead==0)
    flip_cmatrix(selfenergy,N);

  // free vectors and matrices
  cfree_cmatrix(t,0,N-1,0,N-1);
  cfree_cmatrix(P,0,4*N-1,0,4*N-1);
  cfree_cmatrix(wmH,0,4*N-1,0,4*N-1);
  cfree_cmatrix(temp_min,0,2*N-1,0,2*N-1);
  cfree_cmatrix(LAMBDA,0,2*N-1,0,2*N-1);
  cfree_ccvector(eigenvalues,0,2*N-1);
  cfree_cmatrix(w23,0,2*N-1,0,2*N-1);
  cfree_cmatrix(w01,0,2*N-1,0,2*N-1);
  cfree_ccvector(ld,0,2*N-1);
  free_ivector(flag_right,0,2*N-1);
  free_ivector(index_zero_row,0,N-1);
  free_ivector(index_zero_column,0,N-1);
  cfree_cmatrix(MM,0,N-1,0,N-1);
  cfree_ccvector(b,0,N-1);
  cfree_cmatrix(temp2,0,N-1,0,N-1);
  cfree_cmatrix(temp,0,N-1,0,N-1);
  return selfenergy;
  //  exit(0);
}
