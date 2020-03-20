#ifndef ZGEEV_H
#define ZGEEV_H
void zgeev_(char *jobvl, char *jobvr, int *n, complex *ap, int *LDA, complex *w, complex *zr,int *ldz,complex *zl,int *ldzl, complex *work, int *lwork, double *rwork, int *info);
#endif
