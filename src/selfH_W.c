// ======================================================================
//  Copyright (c) 2010, P.D'Amico, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

// Calculates the diagonal blocks in the form (E+i*eta)I-U, 
// where E is the energy and eta a small imaginary part.
// Result: lead self energy.

#include "selfH_W.h"
complex **selfH_W(double E, complex ***diag, complex ***updiag, complex ***lowdiag, int N, int Nc, int lead, double eta)
{
  complex **Gsnew, **Gz, **wmH, zero, **t, **selfenergy, **tdaga;
  int i, j, k;

  zero.r = 0.0;
  zero.i = 0.0;



  wmH = cmatrix(0,4*N-1,0,4*N-1);
  t = cmatrix(0,4*N-1,0,4*N-1);
  //  tdaga = cmatrix(0,4*N-1,0,4*N-1);
  selfenergy = cmatrix(0,N-1,0,N-1);

  for(i=0; i<4*N; i++){
    for(j=0; j<4*N; j++){
      wmH[i][j]=zero;
      t[i][j]=zero;
    } 
  }

  // I take care of the diagonal block matrix
  if (lead==0)
    {
      for(k=0; k<4; k++)
	for(i=0; i<N; i++)
	  for(j=0; j<N; j++)
	    if (i==j)
	      wmH[j+k*N][j+k*N]=complass(E+diag[k][j][j].r,+diag[k][j][j].i+eta);
	    else
	      wmH[i+k*N][j+k*N]=complass(diag[k][i][j].r,diag[k][i][j].i);
    }
  else
    {
      for(k=0; k<4; k++)
	for(i=0; i<N; i++)
	  for(j=0; j<N; j++)
	    if (i==j)
	      wmH[j+k*N][j+k*N]=complass(E+diag[Nc-4+k][j][j].r,+diag[Nc-4+k][j][j].i+eta);
	    else
	      wmH[i+k*N][j+k*N]=complass(diag[Nc-4+k][i][j].r,diag[Nc-4+k][i][j].i);
    }


  // I take care of the off-diagonal block matrix
  if (lead==1)
    {
      for(k=0; k<3; k++){
	for(i=0; i<(1)*N; i++){
	  for(j=0; j<N; j++){
	    wmH[i+k*N][j+(k+1)*N]=complass(+updiag[Nc-4+k][i][j].r,updiag[Nc-4+k][i][j].i);
	  } 
	}
      } 
      
      for(k=1; k<4; k++){
	for(i=0; i<N; i++){
	  for(j=0; j<N; j++){
	    wmH[i+k*N][j+(k-1)*N]=complass(+lowdiag[Nc-4+k][i][j].r,lowdiag[Nc-4+k][i][j].i);
	  } 
	}
      }    
    }
  else
    {
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
    }


  
  if (Nc<=4)
    {
      for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	  t[i+3*N][j]=complass(lowdiag[2][j][i].r,lowdiag[2][j][i].i);
    }
  else
    {
      if (lead==1)
	for(i=0; i<N; i++)
	  for(j=0; j<N; j++)
	    t[i+3*N][j]=complass(updiag[Nc-5][i][j].r,updiag[Nc-5][i][j].i);
      else
	for(i=0; i<N; i++)
	  for(j=0; j<N; j++)
	    t[i][j+3*N]=complass(lowdiag[4][i][j].r,lowdiag[4][i][j].i);
    }

/*   for (i=0;i<4*N;i++) */
/*     for (j=0;j<4*N;j++) */
/*       { */
/* 	tdaga[i][j].r=t[j][i].r; */
/* 	tdaga[i][j].i=-t[j][i].i; */
/*       } */


/*   if (lead==0){ */
/*   printf("LowDiag\n"); */
/*   for (i=0; i<4*N; i++) */
/*     { */
/*       for (j=0; j<4*N; j++) */
/* 	{ */
/* 	  printf("%lg+%lgi ",t[i][j].r,t[i][j].i); */
/* 	} */
/*       printf("\n"); */
/*     } */

/*   printf("wmH\n"); */
/*   for (i=0; i<4*N; i++) */
/*     { */
/*       for (j=0; j<4*N; j++) */
/* 	{ */
/* 	  printf("%lg+%lgi ",wmH[i][j].r,wmH[i][j].i); */
/* 	} */
/*       printf("\n"); */
/*     } */
/*   exit(0); */
/*   } */



  Gz = GzerozeroH_W(wmH, t, 4*N);

/*   printf("Self Total\n"); */
/*   for (i=0; i<4*N; i++){ */
/*     for (j=0; j<4*N; j++){ */
/*       printf("%lg+i%lg ", Gz[i][j].r,Gz[i][j].i); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   exit(0); */



//  if (lead==0)
//    Gz = Gzerozero(wmH, tdaga, t, 4*N);
//  else
/*   Gz = Gzerozero(wmH, t, tdaga, 4*N); */

/*   if (lead==1) */
/*     { */
/*       for (i=0;i<N;i++) */
/* 	for (j=0;j<N;j++) */
/* 	  selfenergy[i][j]=Gz[3*N+i][3*N+j]; */
/*     } */
/*   else */
/*     { */
/*       for (i=0;i<N;i++) */
/* 	for (j=0;j<N;j++) */
/* 	  //selfenergy[j][i]=Gz[4*N-1-i][4*N-1-j]; */
/* 	  selfenergy[i][j]=complass(Gz[3*N+N-1-i][3*N+N-1-j].r,Gz[3*N+N-1-i][3*N+N-1-j].i); */
/*     } */

/*   cfree_cmatrix(wmH,0,4*N-1,0,4*N-1); */
/*   cfree_cmatrix(Gz,0,4*N-1,0,4*N-1); */
/*   cfree_cmatrix(t,0,4*N-1,0,4*N-1); */
/*   cfree_cmatrix(tdaga,0,4*N-1,0,4*N-1); */

/* /\*   printf("Self Total\n"); *\/ */
/* /\*   for (i=0; i<N; i++){ *\/ */
/* /\*     for (j=0; j<N; j++){ *\/ */
/* /\*       printf("%lg+i%lg ", Gz[3*N+i][3*N+j].r,Gz[3*N+i][3*N+j].i); *\/ */
/* /\*     } *\/ */
/* /\*     printf("\n"); *\/ */
/* /\*   } *\/ */

/*   return selfenergy; */
  

  cfree_cmatrix(wmH,0,4*N-1,0,4*N-1);


  Gsnew=cmatmul(t, Gz, 4*N);


/*   if (lead==1) */
/*     { */
/*       printf("Self Total\n"); */
/*       for (i=0; i<4*N; i++){ */
/* 	for (j=0; j<4*N; j++){ */
/* 	  printf("(%lg+i%lg) ", Gsnew[i][j].r,Gsnew[i][j].i); */
/* 	  //printf("(%lg+i%lg) ", Gz[i][j].r,Gz[i][j].i); */
/* 	} */
/* 	printf("\n"); */
/*       } */
/*       exit(0); */
/*     } */


  if(lead == 1){
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
        selfenergy[i][j]=complass(Gsnew[i+3*N][j+3*N].r,Gsnew[i+3*N][j+3*N].i);
      } 
    }
  }
  else{
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
	//selfenergy[j][i]=complass(Gsnew[i][j].r,Gsnew[i][j].i);
        selfenergy[i][j]=complass(Gsnew[i][j].r,Gsnew[i][j].i);
      } 
    } 
  }


/*
  printf("Self\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      printf("(%f) ", selfenergy[i][j].r);
    }
    printf("\n");
  } 
*/




  cfree_cmatrix(Gz,0,4*N-1,0,4*N-1);
  cfree_cmatrix(Gsnew,0,4*N-1,0,4*N-1);
  cfree_cmatrix(t,0,4*N-1,0,4*N-1);

  return selfenergy;
}
