// ======================================================================
//  Copyright (c) 2004-2010, P. D'Amico, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under the BSD license.
//  See the file "license.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

// In wmH you put the EI-H, where H is the Hamiltonian, and 
// E is the energy plus a small imaginary part. 
// t is the hopping matrix between two cells.
// 
// OUTPUT : g_{R}^{r} * t^{\dagger} = Z21*Z11^{-1},
// where g_{R}^{r} is the surface Green function.
// This result is then used to calculate the self-energy
// with the formula Sigma_R = t * g_{R}^{r} * t^{\dagger}
// in the calling routine self_H.c .

#include "GzerozeroH_W.h"
complex **GzerozeroH_W(complex **wmH, complex **t, int N)

{
  complex zero, **Gzeta, **A, **B, **temp1, **temp2, **temp3, **diagonalization;
  int i,j,k,l;

 
  FILE *fp;
  zero.r=0;
  zero.i=0;
 
  A = cmatrix(0,2*N-1,0,2*N-1);
  B = cmatrix(0,2*N-1,0,2*N-1);
  
  diagonalization = cmatrix(0,2*N-1,0,2*N-1);
  temp1 = cmatrix(0,N-1,0,N-1);
  temp2 = cmatrix(0,N-1,0,N-1);
  // temp3 = cmatrix(0,N-1,0,N-1);
/*
  THE FOLLOWING CICLES FILL THE A AND B MATRICES CORRECTLY
*/

/*
  Up-left matrix filling
*/

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      A[i][j]=zero; 
      B[i][j]=zero;
    }
  }

  for (i=0;i<N;i++){
      B[i][i].r=1.0;
  }

/*
  Up-right matrix filling
*/

  for (i=0;i<N;i++){
    for (j=N;j<2*N;j++){
      A[i][j]=zero; 
      B[i][j]=zero;
    }
  }

  for (i=0;i<N;i++){
      A[i][i+N].r=1.0;
  }

/*
  Down-left matrix filling
*/


  for (i=N;i<2*N;i++){
    for (j=0;j<N;j++){
      B[i][j]=zero;
    }
  }

  for (i=N;i<2*N;i++){
    for (j=0;j<N;j++){
      A[i][j].r=-t[j][i-N].r;
      A[i][j].i=t[j][i-N].i;  
    }  
  }


/*
  Down-right matrix filling. Those fillings can be done by reference also
*/

  for (i=N;i<2*N;i++){
    for (j=N;j<2*N;j++){
      A[i][j].r = wmH[i-N][j-N].r;
      A[i][j].i = wmH[i-N][j-N].i;
      B[i][j].r = t[i-N][j-N].r;
      B[i][j].i = t[i-N][j-N].i;
    }
  }


/*
  G(eneralized) E(igen)V(alue) P(roblem) with A and B as imputs: search solutions (lambda and v) for A.v = lambda B.v
*/

// zeig_W(A, B, diagonalization, 2*N); this line for the "standard" method


  zeig_S(A, B, diagonalization, 2*N);

  cfree_cmatrix(A,0,2*N-1,0,2*N-1);
  cfree_cmatrix(B,0,2*N-1,0,2*N-1);

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      temp1[i][j]=zero; 
      temp2[i][j]=zero;
      //temp3[i][j]=zero;
    }
  }

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      temp1[i][j]=diagonalization[i][j]; 
    }
  }

  for (i=N;i<2*N;i++){
    for (j=0;j<N;j++){
      temp2[i-N][j]=diagonalization[i][j]; 
    }
  }

  cfree_cmatrix(diagonalization,0,2*N-1,0,2*N-1);
  temp3 = cmatinv(temp1, N);


  Gzeta = cmatmul(temp2, temp3, N);

  cfree_cmatrix(temp1,0,N-1,0,N-1);
  cfree_cmatrix(temp2,0,N-1,0,N-1);
  cfree_cmatrix(temp3,0,N-1,0,N-1);
  return Gzeta;
}

