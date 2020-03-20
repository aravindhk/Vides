#include "ceig.h"
void ceig(complex **mat,complex **V,double *D,int n,int Nm)
{
  int i,j,k,IL,IU,M,*IWORK,INFO,LDZ,LWORK,LRWORK,LIWORK,*IFAIL;
  double VL,VU,ABSTOL,*W,*RWORK;
/*   char JOBZ[10],RANGE[10],UPLO[10]; */
  char *JOBZ,*RANGE,*UPLO;
  complex zero, *AP,*Z,*WORK;

 /*  printf("+++++++++++++++++\n"); */

/*   for (i=0;i<=5;i++) */
/*     { */
/*       for (j=0;j<=5;j++) */
/* 	{ */
/* 	  printf("%lg+%lgi ",mat[i][j].r,mat[i][j].i); */
/* 	} */
/*       printf("\n"); */
/*     } */
/*   printf("++++++++++++++++++\n"); */


  JOBZ=cvector(0,9);
  RANGE=cvector(0,9);
  UPLO=cvector(0,9);

  zero=complass(0.0,0.0);

  LWORK=2*n;
  LRWORK=7*n;
  LIWORK=5*n;

  LDZ=n;
  AP=ccvector(0,n*(n+1)-1);
  W=dvector(0,n-1);
  Z=ccvector(0,n*Nm-1);
  IFAIL=ivector(0,n-1);
  WORK=ccvector(0,LWORK-1);
  RWORK=dvector(0,LRWORK-1);
  IWORK=ivector(0,LIWORK-1);
  for (j=0;j<LIWORK;j++)
    IWORK[j]=0;
  for (j=0;j<LRWORK;j++)
    RWORK[j]=0;
  for (j=0;j<LWORK;j++)
    WORK[j]=zero;
  for (j=0;j<n;j++)
    W[j]=0;
  for (j=0;j<n*Nm;j++)
    Z[j]=zero;

/*   for (i=0;i<n;i++) */
/*     { */
/*       for (j=0;j<n;j++) */
/* 	printf("%lg ",mat[i][j].r); */
/*       printf("\n"); */
/*     } */

/*   for (i=0;i<n;i++) */
/*     { */
/*       for (j=0;j<n;j++) */
/* 	printf("%lg ",mat[i][j].i); */
/*       printf("\n"); */
/*     } */
  
  k=0;
  for (j=0;j<n;j++)
    for (i=0;i<=j;i++)
      {
	AP[k]=mat[i][j];
	k++;
      }

  VL=0;VU=0;
  strcpy(JOBZ,"N");   /* N: zhpevx computes only the eigenvalues. V: zhpevx computes also the eigenvectors */

  strcpy(RANGE,"I");
  strcpy(UPLO,"U");
  IL=1;
  IU=n;

  IL=n/2+1-Nm/2;
  IU=n/2+Nm/2;

  //  printf("autovaliamo %d %d \n",IL,IU);
  ABSTOL=1e-14; /* FUNDAMENTAL !!!!! (1/4/2010: originariamente era 1e-5) */
  M=Nm;
 
  zhpevx_(JOBZ,RANGE,UPLO,&n,AP,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z,&LDZ,	  
	  WORK,RWORK,IWORK,IFAIL,&INFO);
  
  if (INFO!=0)
    {
      printf("ERROR in CHEEVR %d \n",INFO);
      exit(1);
    }

 /*  if (INFO!=0) */
/*     printf("INFO = %d \n",INFO); */
  
  //printf("Matrix order %d \n",n);
  //printf("Number of eigenvalues %d \n",Nm);
  
  for (i=0;i<Nm;i++)
    {
      D[i]=W[i];
      //      printf("%lg \n",W[i]);
    }
  
/*   exit(0); */
  
  // free memory in free state!

  cfree_ccvector(AP,0,n*(n+1)-1);
  free_dvector(W,0,n-1);
  cfree_ccvector(Z,0,n*Nm-1);
/*   free_ccvector(WORK,0,8*n); */
  cfree_ccvector(WORK,0,LWORK-1);
  free_dvector(RWORK,0,LRWORK-1);
  free_ivector(IWORK,0,5*n);
  free_ivector(IFAIL,0,n);

  free_cvector(JOBZ,0,9);
  free_cvector(RANGE,0,9);
  free_cvector(UPLO,0,9);

  return;
}
