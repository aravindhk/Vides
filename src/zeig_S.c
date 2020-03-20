// ======================================================================
//  Copyright (c) 2010, P.D'Amico, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

// This rountine take as input the two matrices 
//
// A = |0           , 1  |  and   B = |1, 0|
//     |-t^{\dagger}, E-U|            |0, t| 
// 
// and solves the generalized problem A*v = lambda B*v
// using the Schur decomposition of the matrices A and B.
// Basically there exist two hermitian matrices, Q and Z,
// such that S = Q^{\dagger} * A * Z and T = Q^{\dagger} * B * Z, 
// where T and S are upper triangular matrices in the canonical 
// Schur form, i.e. the eigenvalues are on the diagonal.
// The eigenvalues lambda_i are given expressed in term of 
// the diagonal elements through S(i,i) / T(i,i).
// Then it reorders the eigenvectors in two groups, distinguishing
// between |lambda|<1 (call it lambda_<) and |lambda|>1 (call it lambda_>) and rearrange
// the Schur matrices in order to have the lambda_< on the left part of the
// triangular matrices.
// Notice that the matrices here are (2N x 2N) where N is the number of atoms 
// (eventually multiplied by the number of orbitals involved in the problem)
// that are in the unitary cell of the leads
// The method is built following the Michael Wimmer PhD Thesis (insert reference). 

#include "zeig_S.h"

void zeig_S(complex **A, complex**B, complex **diagonalization, int N)

{ 

  complex *AA, *BB, *ALPHA, *BETA, *VSL, *VSR, *WORK;
  int *BWORK, *DY, SELCTG;
  int i, j, ix, LDA, LDB, SDIM, INFO, LDVSL, LDVSR, LWORK;
  double *RWORK;
  char JOBVSL[10], JOBVSR[10], SORT[10];

  strcpy(JOBVSL,"V");   /* V (N): zgges does (not) computes the Left Schur eigenvectors*/
                       
  strcpy(JOBVSR,"V");   /* V (N): zgges does (not) computes the Right Schur eigenvectors*/

  strcpy(SORT,"N");  /* Specifies whether or not to order the eigenvalues on the diagonal of the generalized Schur form.*/
 		     /* S (N) Eigenvalues are (not) ordered. If S is selected, the eigenvectors are ordered on the diagonal according to SELCTG*/

  LDA = N;
  LDB = N;
  LDVSL = N;
  LDVSR = N;
  LWORK = 2*N;

  AA=ccvector(0,LDA*N-1);
  BB=ccvector(0,LDB*N-1);
  ALPHA = ccvector(0,N-1);
  BETA = ccvector(0,N-1);
  VSL = ccvector(0,LDVSL*N-1);
  VSR = ccvector(0,LDVSR*N-1);
  WORK = ccvector(0,LWORK-1);
  RWORK = dvector(0,8*N-1);

  BWORK = ivector(0,N-1);

  DY = ivector(0,N-1);

  for (i=0;i<N;i++){
      for (j=0;j<N;j++){
	ix=i+j*N;
	AA[ix].r=A[i][j].r;
	AA[ix].i=A[i][j].i;
	BB[ix].r=B[i][j].r;
        BB[ix].i=B[i][j].i;
      }
  }

  zgges_( JOBVSL, JOBVSR, SORT, &SELCTG, &N, AA, &LDA, BB, &LDB, &SDIM, ALPHA, BETA, VSL, &LDVSL, VSR, &LDVSR, WORK, &LWORK, RWORK, BWORK, &INFO );

  double numerator, denominator;
  numerator = 0.0;
  denominator = 0.0;


  for(i=0; i<N; i++){
    DY[i] = 0;
  }  


/*
  The following cicle assigns to DY[i] 1 if the corresponding eigenvalue has to be moved in the
  first columns (Z11 and Z21). 
*/

  for (i=0; i<N; i++){
    numerator = sqrt(ALPHA[i].r*ALPHA[i].r + ALPHA[i].i*ALPHA[i].i);
    denominator = sqrt(BETA[i].r*BETA[i].r + BETA[i].i*BETA[i].i);
    if (numerator < denominator){
      DY[i] = 1;
    }
  }
 

  int IJOB, WANTQ, WANTZ, M, LWORK1, *IWORK, LIWORK, INFO1;
  double  PL, PR, *DIF[2];
  complex *WORK1;

  IJOB = 0;
  WANTQ = 0;
  WANTZ = 1;
  LWORK1 = N;
  LIWORK = N;

  WORK1 = ccvector(0,LWORK1-1);
  IWORK = ivector(0,LIWORK-1);
  
  ztgsen_( &IJOB, &WANTQ, &WANTZ, DY, &N, AA, &LDA, BB, &LDB, ALPHA, BETA, VSL, &LDVSL, VSR, &LDVSR, &M, &PL, &PR, DIF, WORK1, &LWORK1, IWORK, &LIWORK, &INFO1 );

/*
  In "diagonalization" I store the matrix VSR (the rows are the generalized eigenvectors)
  sorted accordingly to the lambda_< and lambda_> criterion.
*/

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      diagonalization[i][j] = VSR[i+j*N];
    }
  }

  cfree_ccvector(AA, 0, LDA*N-1);
  cfree_ccvector(BB, 0, LDB*N-1);
  cfree_ccvector(ALPHA, 0, N-1);
  cfree_ccvector(BETA, 0, N-1);
  cfree_ccvector(VSL, 0, LDVSR*N-1);
  cfree_ccvector(VSR, 0, LDVSR*N-1);
  cfree_ccvector(WORK, 0, LWORK-1);
  free_dvector(RWORK, 0, 8*N-1);
  free_ivector(BWORK, 0, N-1);
  free_ivector(DY, 0, N-1);
  cfree_ccvector(WORK1, 0, LWORK1-1);
  free_ivector(IWORK, 0, LIWORK-1);
  
}

