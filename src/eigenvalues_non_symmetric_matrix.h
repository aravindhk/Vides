#include "complex.h"
#include "cmatrix.h"
#include "complass.h"
#include "cvectorm.h"
#include "cfree_cmatrix.h"
#include "nrutil.h"
#include "zgeev.h"
#include "ccvector.h"
#include "cfree_ccvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef EIGENVALUES_NON_SYMMETRIC_MATRIX_H
#define EIGENVALUES_NON_SYMMETRIC_MATRIX_H
void eigenvalues_non_symmetric_matrix(complex **mat,complex *D,complex **Z,int n,int Nm,int lead);
#endif
