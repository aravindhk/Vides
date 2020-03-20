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

#include "bands.h"
void bands(complex ***diag, complex ***updiag, complex ***lowdiag, int N)
{
  int i, j, k;
  complex **matrix, **V;
  double x, *D;
  matrix = cmatrix(0,4*N-1,0,4*N-1);
  for (i=0;i<4*N;i++)
    for (j=0;j<4*N;j++)
      matrix[i][j]=complass(0,0);


  V = cmatrix(0,4*N-1,0,4*N-1);
  D = dvector(0,4*N-1);

  for (k=0;k<4;k++){
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        matrix[i+k*N][j+k*N].r = -diag[k][i][j].r;
      }
    }
  }

  for (k=0;k<3;k++){
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        matrix[i+k*N][j+(k+1)*N].r = -updiag[k][i][j].r;
      }
    }
  }

  for (k=1;k<4;k++){
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        matrix[i+k*N][j+(k-1)*N].r = -lowdiag[k][i][j].r;
      }
    }
  }

  for (i=0;i<4*N;i++)
    {
      for (j=0;j<4*N;j++)
	printf("%lg+i%lg ",matrix[i][j].r,matrix[i][j].i);
      printf("\n");
    }
  //  exit(0);


  x = 0.0;

  FILE *fp;
  fp=fopen("Bands.dat","w");
  
  for (x=0; x<3.14159; x+=0.1){
    
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
	matrix[i][j+3*N].r = -cos(x)*lowdiag[4][i][j].r;
	matrix[i][j+3*N].i = +sin(x)*lowdiag[4][i][j].r;
      }
    }
    
    
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
	matrix[i+3*N][j].r = -cos(x)*updiag[3][i][j].r;
	matrix[i+3*N][j].i = -sin(x)*updiag[3][i][j].r;
      }
    }
    
    ceig(matrix, V, D, 4*N, 4*N);
    
    fprintf(fp,"%f ", x);
    for (i=0;i<4*N;i++){
      fprintf(fp,"%f ", D[i]);
    }
    fprintf(fp,"\n");
    
  }
  
  fclose(fp);

//  exit(0);


/*
  complex diff;
  diff.i = 0.0;
  diff.r = 0.0;

  for(i=0; i<4*N; i++){
    for(j=i; j<4*N; j++){
      diff.r = diff.r + matrix[i][j].r - matrix[j][i].r;
      diff.i = diff.i + matrix[i][j].i + matrix[j][i].i;
    }  
  }

  printf("Test for Hermiticity\n"); 
  printf("(%1.3f,%1.3f)\n", diff.r, diff.i); 

  complex **V;
  double *D;
  complex  **test;
  V = cmatrix(0,4*N-1,0,4*N-1);
  D = dvector(0,4*N-1);
  test = cmatrix(0,3,0,3);

  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      test[i][j].r = 0.0;
      test[i][j].i = 0.0;
    }
  }

    for (i=0;i<4;i++){
      test[i][i].r = 1.0;
    }  
*/
/*
printf("Matrix\n");  
  printf("\n");
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      printf("(%1.0f) ", matrix[i+3*N][j].r); 
    }  
    printf("\n");
  }
*/
//  exit(0);


/*
  for(j=0; j<4*N; j++){
    printf("(%1.2f) ", D[j]); 
  } 


  exit(0);
*/

  return;
}
