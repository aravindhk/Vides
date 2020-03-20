#include "complex.h"
#ifndef ZGEMM_H
#define ZGEMM_H
void zgemm_(char *TRANSA,char *TRANSB,int *M,int *N,int *k,complex *ALPHA,
	    complex *A,int *LDA,complex *B,int *LDB,complex *BETA,complex *C,int *LDC);
#endif
