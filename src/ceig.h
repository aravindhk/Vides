#include "complex.h"
#include "cmatrix.h"
#include "complass.h"
#include "cvectorm.h"
#include "cfree_cmatrix.h"
#include "nrutil.h"
#include "zhpevx.h"
#include "ccvector.h"
#include "cfree_ccvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef CEIG_H
#define CEIG_H
void ceig(complex **mat,complex **V,double *D,int n,int Nm);
#endif
