// ======================================================================
//  Copyright (c) 2011, P.D'Amico, G. Fiori, G. Iannaccone  University of Pisa
//  
//  This file is released under a RESTRICTED LICENSE.
//  See the file "licenseRESTRICTED.txt" for information on usage and
//  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
// ====================================================================== 

#include "Hamiltonian_py.h"

void py_hamiltonian(PyArrayObject* hamiltonian, int n, int M, int r, int c, double *Phi,complex ****diag, complex ****updiag, complex ****lowdiag)
{
  int i, j, k, l, s, ii, jj, N;
  int Num_dummy;

  N=n*c;

  import_array();

/*
  if (!PyArg_ParseTuple(args,"OII", &hamiltonian, &N, &M))
  {
    printf("NOT A CLASS! \n");
    exit(0);
  }
*/

//  printf("# of atoms * # of orbitals (in one slice) = %i\n", N); 
//  printf("# of slices = %i\n", M);

  /*Memory allocation*/  
  
  *diag=c3tensor(0,M-1,0,N-1,0,N-1);
  *lowdiag=c3tensor(0,M-1,0,N-1,0,N-1);
  *updiag=c3tensor(0,M-1,0,N-1,0,N-1);

  
  for(k=0; k<M; k++){
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
        (*diag)[k][i][j].r=0;
        (*diag)[k][i][j].i=0;
        (*lowdiag)[k][i][j].r=0;
        (*lowdiag)[k][i][j].i=0;
        (*updiag)[k][i][j].r=0;
        (*updiag)[k][i][j].i=0;
      }
    }
  }
  
  //  printf("QUI 2\n");
  
//  exit(0);
/*
  The #columns of the file that contains the Hamiltonian stored in the format i j t_{ij}^{1} ...... t_{ij}^{q},
  is given by the integer variable c.
*/
  int alpha, beta, gamma, delta,strides1;
  //  double **p = malloc((c*c) * sizeof (double*));
  complex *p;
  double *id, *jd;
  Py_complex *tempiamo;
//int *ptr = malloc(10 * sizeof (int));
  p=ccvector(0,c*c-1);

  //  printf("%d %d \n",hamiltonian->strides[1],sizeof(complex));
  //  exit(0);
  strides1=sizeof(complex);
  
  for(l=1; l<r; l++ ){
    for(s=0; s<c*c; s++){ 
      id =  (double *)(hamiltonian->data + l*hamiltonian->strides[0]); 
      jd =  (double *)(hamiltonian->data + l*hamiltonian->strides[0] + strides1);
      //      printf("%lg %lg \n",*id,*jd);
      //      p[s]  =  (double *)(hamiltonian->data + l*hamiltonian->strides[0] + (2 + s)*hamiltonian->strides[1]);
      tempiamo=(Py_complex *)(hamiltonian->data + l*hamiltonian->strides[0] + (2 + s)*strides1);
      p[s]  =  ConvertPycomplex_to_Ccomplex(*tempiamo);
      
      // I pass -H
      p[s].r=-p[s].r;
      p[s].i=-p[s].i;
	

      ii = (((int)(*id)) - 1)*(c) + 1 + s/c;
      jj = (((int)(*jd)) - 1)*(c) + 1 + s%c;

      alpha = (int)(((ii-1))/N);
      beta = (int)(((jj-1))/N);
      gamma = ((ii -1))%N;
      delta = ((jj -1))%N;
      //      printf("%d %d %d %d \n",alpha,beta,gamma,delta);

/*
      alpha = (int)(((ii-1)*c)/N);
      beta = (int)(((jj-1)*c)/N);
      gamma = ((ii -1)*(s+1))%N;
      delta = ((jj -1)*(s+1))%N;
*/

      if (alpha == beta){
        (*diag)[alpha][gamma][delta] = p[s];
	(*diag)[alpha][delta][gamma] = complass(p[s].r,-p[s].i);
        //(*diag)[alpha][delta][gamma].r = *p[s];
        //(*diag)[alpha][delta][gamma].i = 0;
      }
      else{
      //  if (alpha < beta){
          (*updiag)[beta-1][gamma][delta] = p[s];
          (*lowdiag)[alpha+1][delta][gamma] = complass(p[s].r,-p[s].i);

          //(*updiag)[beta-1][delta][gamma].r = *p[s];
          //(*updiag)[beta-1][delta][gamma].i = 0;
       // }
      /*  else{
          (*lowdiag)[alpha][gamma][delta].r = *p[s];
          (*lowdiag)[alpha][gamma][delta].i = 0;
          //(*lowdiag)[alpha][delta][gamma].r = *p[s];
          //(*lowdiag)[alpha][delta][gamma].i = 0;
        }*/  
      }
    }
  }


// Now I add the on-site energy
  Num_dummy=0;
  for(k=0; k<M; k++){
    for(i=0; i<n; i++)
      {     
	for(j=0; j<c; j++)
	  {
	    if (fabs(fabs((*diag)[k][i*c][i*c].r)-77777)>1e-3)
	      (*diag)[k][i*c+j][i*c+j].r+=Phi[k*n+i-Num_dummy];
	  }

	// I check the dummy atoms
	// right now it works only for 1 orbital
	if (fabs(fabs((*diag)[k][i*c][i*c].r)-77777)<1e-3)
	  Num_dummy++;
      }
  }
  //  printf("QUI 3\n");


/* FROM HERE ... 

    if ((ii)==(jj)){
      k=(int)((ii-1)/N);
      ii=(ii-1)%N;
      jj=(jj-1)%N;
      //printf("k= %d\n",k);
      (*diag)[k][ii][jj].r=*p;
      (*diag)[k][ii][jj].i=0;
    }



    else{

      if((ii)>(jj)){
        k=(int)((ii-1)/N);
        ii=(ii-1)%N;
        jj=(jj-1)%N;
        //printf("k= %d\n", k);
        (*lowdiag)[k][ii][jj].r=*p;
        (*lowdiag)[k][ii][jj].i=0;
      }


      else{
        k=(int)((jj-1)/N)-1;
        ii=(ii-1)%N;
        jj=(jj-1)%N;
        //printf("k= %d\n",k);
        (*updiag)[k][ii][jj].r=*p;
        (*updiag)[k][ii][jj].i=0;
      }   
    }
  }


.. TO HERE*/

/* Print on screen the matrix-system structure*/

//  printf("SYSTEM COMPOSED OF %i CELLS WITH %i ATOMS and %i ORBITALS PER ATOM \n", M, N, c); 

//exit(0);

/*  
  printf("DIAGONAL BLOCKS\n");  
  printf("\n");
  for(k=0; k<M; k++){
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
        printf("%f ", (*diag)[k][i][j].r); 
      }  
      printf("\n");
    }
    printf("\n");
    printf("\n");
  }
*/

/*   printf("UP-DIAGONAL BLOCKS\n");  */
/*   printf("\n"); */
/*   for(k=0; k<M; k++){ */
/*     for(i=0; i<N; i++){ */
/*       for(j=0; j<N; j++){ */
/*         printf("%f ", (*updiag)[k][i][j].r);  */
/*       }   */
/*       printf("\n"); */
/*     } */
/*     printf("\n"); */
/*     printf("\n"); */
/*   } */


/*   printf("LOW-DIAGONAL BLOCKS\n");  */
/*   printf("\n"); */
/*   for(k=0; k<M; k++){ */
/*     for(i=0; i<N; i++){ */
/*       for(j=0; j<N; j++){ */
/*         printf("%f ", (*lowdiag)[k][i][j].r);  */
/*       }   */
/*       printf("\n"); */
/*     } */
/*     printf("\n"); */
/*     printf("\n"); */
/*     } */

  //  exit(0);


//  printf("QUI 3\n");

//  char *stringa = "Hello from C!";
//  return Py_BuildValue("s", stringa);
}
