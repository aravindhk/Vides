#include "eigenvalues_non_symmetric_matrix.h"
//#include "selectfunction.h"
void eigenvalues_non_symmetric_matrix(complex **mat,complex *D,complex **Z,int n,int Nm,int lead)
{
  int i,j,k,IL,IU,M,*IWORK,INFO,LDZR,LDZL,LWORK,LRWORK,LIWORK,SDIM,LDA;
  double VL,VU,ABSTOL,*RWORK;
/*   char JOBVL[10],JOBVR[10],UPLO[10]; */
  char *JOBVL,*JOBVR,*UPLO;
  complex zero, *AP,*ZR,*WORK,*W,*ZL;

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

//  int (*select)(complex *a);
//  select=selectfunction_;
//  printf("%d \n",selectfunction_(&zero));
	 
  // If lead==1 then I compute the right eigenvectors
  // If lead==0 then I compute the left eigenvectors

  JOBVL=cvector(0,9);
  JOBVR=cvector(0,9);

  zero=complass(0.0,0.0);

  LWORK=2*n;
  LRWORK=2*n;

  LDA=n;
  LDZR=n;
  LDZL=n;
  AP=ccvector(0,n*n-1);
  W=ccvector(0,n-1);
  ZR=ccvector(0,n*Nm-1);
  ZL=ccvector(0,n*Nm-1);
  WORK=ccvector(0,LWORK-1);
  RWORK=dvector(0,LRWORK-1);
  for (j=0;j<LRWORK;j++)
    {
      RWORK[j]=0;
      WORK[j]=zero;
    }
  for (j=0;j<n;j++)
    W[j]=zero;
  for (j=0;j<n*Nm;j++)
    {
      ZR[j]=zero;
      ZL[j]=zero;
    }

/*   FILE *fp; */
/*   fp=fopen("matrix.out","w"); */
/*   for (i=0;i<n;i++) */
/*     { */
/*       for (j=0;j<n;j++) */
/* 	{ */
/* 	  fprintf(fp,"%lg ",mat[i][j].r); */
/* 	} */
/*       fprintf(fp,"\n"); */
/*     } */

  if (lead==1)
    {
      k=0;
      for (i=0;i<n;i++)
	for (j=0;j<n;j++)
	  {
	    AP[k]=mat[j][i];
	    k++;
	  }
    }
  else
    {
      k=0;
      for (i=0;i<n;i++)
	{
	  for (j=0;j<n;j++)
	    {
	      AP[k]=mat[i][j];
	      k++;
	    }
	}
    }

  strcpy(JOBVL,"N");   /* N: zhpevx computes only the eigenvalues. V: zhpevx computes also the eigenvectors */
  strcpy(JOBVR,"V");
  
  //  printf("sono qua prima \n");
  zgeev_(JOBVL,JOBVR,&n,AP,&LDA,W,ZL,&LDZL,ZR,&LDZR,WORK,&LWORK,RWORK,&INFO);
  //  printf("sono qua dopo %lg \n",WORK[0]);

  if (INFO!=0)
    {
      printf("ERROR in ZGEEV %d \n",INFO);
      exit(0);
    }

 /*  if (INFO!=0) */
/*     printf("INFO = %d \n",INFO); */
  
  //printf("Matrix order %d \n",n);
  //printf("Number of eigenvalues %d \n",Nm);
  

/*   if (lead==1) */
    {
      for (i=0; i<n; i++){
	for (j=0; j<n; j++){
	  Z[i][j] = ZR[i+j*n];
	}
      }
    }
/*   else */
/*     { */
/*       for (i=0; i<n; i++){ */
/* 	for (j=0; j<n; j++){ */
/* 	  Z[i][j] = ZR[i+j*n]; */
/* 	} */
/*       } */
/*     } */

  for (i=0;i<Nm;i++)
    {
      D[i]=W[i];
      //printf("%lg \n",D[i].r);
    }
  // free memory in free state!

  cfree_ccvector(AP,0,n*n-1);
  cfree_ccvector(W,0,n-1);
  cfree_ccvector(ZR,0,n*Nm-1);
  cfree_ccvector(ZL,0,n*Nm-1);
/*   free_ccvector(WORK,0,8*n); */
  cfree_ccvector(WORK,0,LWORK-1);
  free_dvector(RWORK,0,LRWORK-1);
  free_cvector(JOBVL,0,9);
  free_cvector(JOBVR,0,9);

  return;
}
