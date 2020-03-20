*DECK DOMN 
      SUBROUTINE DOMN(N,B,X,NELT,IA,JA,A,ISYM,  
     $     MATVEC,MSOLVE,  
     $     NSAVE,ITOL,TOL,ITMAX,ITER,ERR,IERR,IUNIT,R,Z,P,  
     $     AP,EMAP,DZ,CSAV,RWORK,IWORK ) 
C***BEGIN PROLOGUE  DOMN 
C***DATE WRITTEN   890404   (YYMMDD) 
C***REVISION DATE  890404   (YYMMDD) 
C***CATEGORY NO.  D2A4 
C***KEYWORDS  LIBRARY=SLATEC(SLAP), 
C             TYPE=DOUBLE PRECISION(DOMN-D), 
C             Non-Symmetric Linear system, Sparse,  
C             Iterative Precondition, Orthomin 
C***AUTHOR  Greenbaum, Anne, Courant Institute 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-300 
C             Livermore, CA 94550 (415) 423-3141 
C             seager@lll-crg.llnl.gov 
C***PURPOSE  Preconditioned Orthomin Sparse Iterative Ax=b Solver. 
C            Routine to solve a general linear system  Ax = b  using  
C            the Preconditioned Orthomin method. 
C***DESCRIPTION 
C *Usage: 
C     INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
C     INTEGER   ITER, IERR, IUNIT, IWORK(USER DEFINED) 
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N) 
C     DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE) 
C     DOUBLE PRECISION DZ(N), CSAV(NSAVE), RWORK(USER DEFIED) 
C     EXTERNAL MATVEC, MSOLVE 
C 
C     CALL DOMN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,  
C    $     NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R,  
C    $     Z, P, AP, EMAP, DZ, PSAV, APSV, QSAV, CSAV, RWORK, IWORK) 
C 
C *Arguments: 
C N      :IN       integer . 
C         Order of the Matrix. 
C B      :IN       Double Precision B(N). 
C         Right-hand side vector. 
C X      :INOUT    Double Precision X(N). 
C         On input X is your initial guess for solution vector. 
C         On output X is the final approximate solution. 
C NELT   :IN       integer . 
C         Number of Non-Zeros stored in A. 
C IA     :IN       integer  IA(NELT). 
C JA     :IN       integer  JA(NELT). 
C A      :IN       Double Precision A(NELT). 
C         These arrays contain the matrix data structure for A. 
C         It could take any form.  See "LONG DESCRIPTION", below 
C         for more late breaking details... 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all nonzero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C MATVEC :EXT      External. 
C         Name of a routine which performs the matrix vector multiply 
C         Y = A*X given A and X.  The name of the MATVEC routine must  
C         be declared external in the calling program.  The calling  
C         sequence to MATVEC is: 
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM ) 
C         Where N is the number of unknowns, Y is the product A*X 
C         upon return X is an input vector, NELT is the number of  
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.   
C         ISYM is a flag which, if non-zero, denotest that A is  
C         symmetric and only the lower or upper triangle is stored. 
C MSOLVE :EXT      External. 
C         Name of a routine which solves a linear system MZ = R for 
C         Z given R with the preconditioning matrix M (M is supplied via 
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  
C         be declared external in the calling program.  The calling  
C         sequence to MSOLVE is: 
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C         Where N is the number of unknowns, R is the right-hand side  
C         vector, and Z is the solution upon return.  RWORK is a  
C         double precision  
C         array that can be used to pass necessary preconditioning  
C         information and/or workspace to MSOLVE.  IWORK is an integer   
C         work array for the same purpose as RWORK. 
C NSAVE  :IN       integer . 
C         Number of  direction vectors to save and orthogonalize  
C         against.  NSAVE >= 0. 
C ITOL   :IN       integer . 
C         Flag to indicate type of convergence criterion. 
C         If ITOL=1, iteration stops when the 2-norm of the residual  
C         divided by the 2-norm of the right-hand side is less than TOL. 
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  
C         residual divided by the 2-norm of M-inv times the right hand  
C         side is less than TOL, where M-inv is the inverse of the  
C         diagonal of A. 
C         ITOL=11 is often useful for checking and comparing different  
C         routines.  For this case, the user must supply the "exact"  
C         solution or a very accurate approximation (one with an error  
C         much less than TOL) through a common block, 
C                     COMMON /SOLBLK/ SOLN(1) 
C         if ITOL=11, iteration stops when the 2-norm of the difference  
C         between the iterative approximation and the user-supplied 
C         solution divided by the 2-norm of the user-supplied solution  
C         is less than TOL.  Note that this requires the user to set up 
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine.  
C         The routine with this declaration should be loaded before the 
C         stop test so that the correct length is used by the loader.   
C         This procedure is not standard Fortran and may not work  
C         correctly on your system (although it has worked on every 
C         system the authors have tried).  If ITOL is not 11 then this 
C         common block is indeed standard Fortran. 
C TOL    :IN       Double Precision. 
C         Convergence criterion, as described above. 
C ITMAX  :IN       integer . 
C         Maximum number of iterations. 
C ITER   :OUT      integer . 
C         Number of iterations required to reach convergence, or  
C         ITMAX+1 if convergence criterion could not be achieved in  
C         ITMAX iterations. 
C ERR    :OUT      Double Precision. 
C         Error estimate of error in final approximate solution, as  
C         defined by ITOL. 
C IERR   :OUT      integer . 
C         Return error flag. 
C           IERR = 0 => All went well. 
C           IERR = 1 => Insufficient storage allocated  
C                       for WORK or IWORK. 
C           IERR = 2 => Method failed to converge in  
C                       ITMAX steps. 
C           IERR = 3 => Error in user input.  Check input 
C                       value of N, ITOL. 
C           IERR = 4 => User error tolerance set too tight. 
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded. 
C           IERR = 5 => Preconditioning matrix, M,  is not  
C                       Positive Definite.  $(r,z) < 0.0$. 
C           IERR = 6 => Breakdown of method detected. 
C                       $(p,Ap) < epsilon**2$. 
C IUNIT  :IN       integer . 
C         Unit number on which to write the error at each iteration,  
C         if this is desired for monitoring convergence.  If unit  
C         number is 0, no writing will occur. 
C R      :WORK     Double Precision R(N). 
C Z      :WORK     Double Precision Z(N). 
C P      :WORK     Double Precision P(N,0:NSAVE). 
C AP     :WORK     Double Precision AP(N,0:NSAVE). 
C EMAP   :WORK     Double Precision EMAP(N,0:NSAVE). 
C DZ     :WORK     Double Precision DZ(N). 
C CSAV   :WORK     Double Precision CSAV(NSAVE) 
C RWORK  :WORK     Double Precision RWORK(USER DEFINED). 
C         Double Precision array that can be used for workspace in  
C         MSOLVE. 
C IWORK  :WORK     integer  IWORK(USER DEFINED). 
C         integer  array that can be used for workspace in MSOLVE. 
C 
C *Precision:           Double Precision 
C *See Also: 
C         DSDOMN, DSLUOM, ISDOMN 
C 
C *Description 
C       This routine does  not care  what matrix data   structure is 
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE 
C       routines, with  the arguments as  described above.  The user 
C       could write any type of structure and the appropriate MATVEC 
C       and MSOLVE routines.  It is assumed  that A is stored in the 
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is 
C       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP 
C       routines DSDOMN and DSLUOM are examples of this procedure. 
C        
C       Two  examples  of  matrix  data structures  are the: 1) SLAP 
C       Triad  format and 2) SLAP Column format. 
C        
C       =================== S L A P Triad format =================== 
C       In  this   format only the  non-zeros are  stored.  They may 
C       appear  in *ANY* order.   The user  supplies three arrays of 
C       length NELT, where  NELT  is the number  of non-zeros in the 
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero 
C       the  user puts   the row  and  column index   of that matrix 
C       element in the IA and JA arrays.  The  value of the non-zero 
C       matrix  element is  placed in  the corresponding location of 
C       the A  array.  This is  an extremely easy data  structure to 
C       generate.  On  the other hand it  is  not too  efficient  on 
C       vector  computers   for the  iterative  solution  of  linear 
C       systems.  Hence, SLAP  changes this input  data structure to 
C       the SLAP   Column  format for the  iteration (but   does not 
C       change it back). 
C        
C       Here is an example of the  SLAP Triad   storage format for a 
C       5x5 Matrix.  Recall that the entries may appear in any order. 
C 
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. 
C                              1  2  3  4  5  6  7  8  9 10 11 
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       =================== S L A P Column format ================== 
C       This routine requires  that the  matrix  A be  stored in the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C        
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a  
C       column): 
C        
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C       1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C        
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  MATVEC, MSOLVE, ISDOMN,  
C                    DCOPY, DDOT, DAXPY, D1MACH 
C***END PROLOGUE  DOMN 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
      INTEGER   ITER, IERR, IUNIT, iwork(nelt) 
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)  
      DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE) 
      DOUBLE PRECISION DZ(N), CSAV(NSAVE), rwork(nelt)
      EXTERNAL MATVEC, MSOLVE 
C 
C         Check some of the input data. 
C***FIRST EXECUTABLE STATEMENT  DOMN 
      iunit=0
      open(unit=iunit,file='uscita.out',status='unknown')
      ITER = 0 
      IERR = 0 
      IF( N.LT.1 ) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
      EPS = D1MACH(3) 
      IF( TOL.LT.500.0*EPS ) THEN 
         TOL = 500.0*EPS 
         IERR = 4 
      ENDIF 
      FUZZ = EPS*EPS 
C          
C         Calculate initial residual and pseudo-residual, and check 
C         stopping criterion. 
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM) 
      DO 10 I = 1, N 
         R(I)  = B(I) - R(I)
 10   CONTINUE 
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C          
      IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     R, Z, P, AP, EMAP, DZ, CSAV, 
     $     RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200 
      IF( IERR.NE.0 ) RETURN 
C          
C          
C         ***** iteration loop ***** 
C          
CVD$R NOVECTOR 
CVD$R NOCONCUR 
      DO 100 K = 1, ITMAX 
         ITER = K 
         IP = MOD( ITER-1, NSAVE+1 ) 
C          
C         calculate direction vector p, a*p, and (m-inv)*a*p, 
C         and save if desired. 
         CALL DCOPY(N, Z, 1, P(1,IP), 1) 
         CALL MATVEC(N, P(1,IP), AP(1,IP), NELT, IA, JA, A, ISYM) 
         CALL MSOLVE(N, AP(1,IP), EMAP(1,IP), NELT, IA, JA, A, ISYM, 
     $        RWORK, IWORK) 
         IF( NSAVE.EQ.0 ) THEN 
            AKDEN = DDOT(N, EMAP, 1, EMAP, 1) 
         ELSE 
            IF( ITER.GT.1 ) THEN 
               LMAX = MIN( NSAVE, ITER-1 ) 
               DO 20 L = 1, LMAX 
                  IPO = MOD(IP+(NSAVE+1-L),NSAVE+1) 
                  BKL = DDOT(N, EMAP(1,IP), 1, EMAP(1,IPO), 1) 
                  BKL = BKL*CSAV(L) 
                  CALL DAXPY(N, -BKL,    P(1,IPO), 1,    P(1,IP), 1) 
                  CALL DAXPY(N, -BKL,   AP(1,IPO), 1,   AP(1,IP), 1) 
                  CALL DAXPY(N, -BKL, EMAP(1,IPO), 1, EMAP(1,IP), 1) 
 20            CONTINUE 
               IF( NSAVE.GT.1 ) THEN 
                  DO 30 L = NSAVE-1, 1, -1 
                     CSAV(L+1) = CSAV(L) 
 30               CONTINUE 
               ENDIF 
            ENDIF 
            AKDEN = DDOT(N, EMAP(1,IP), 1, EMAP(1,IP), 1) 
            IF( ABS(AKDEN).LT.EPS*EPS ) THEN 
               IERR = 6 
               RETURN 
            ENDIF 
            CSAV(1) = 1./AKDEN 
C          
C         calculate coefficient ak, new iterate x, new residual r, and 
C         new pseudo-residual z. 
         ENDIF 

         AKNUM = DDOT(N, Z, 1, EMAP(1,IP), 1) 
         AK = AKNUM/AKDEN 
         CALL DAXPY(N,  AK,    P(1,IP), 1, X, 1) 
         CALL DAXPY(N, -AK,   AP(1,IP), 1, R, 1) 
         CALL DAXPY(N, -AK, EMAP(1,IP), 1, Z, 1) 
C          
C         check stopping criterion. 
         IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE, 
     $        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $        R, Z, P, AP, EMAP, DZ, CSAV, 
     $        RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200 
C          
 100  CONTINUE 
C          
C         *****   end of loop  ***** 
C          
C         Stopping criterion not satisfied. 
      ITER = ITMAX + 1 
      IERR = 2 
C          
 200  continue
      close(iunit)
      RETURN 
C------------- LAST LINE OF DOMN FOLLOWS ---------------------------- 
      END 
*DECK DSDOMN 
      SUBROUTINE DSDOMN(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     RWORK, LENW, IWORK, LENIW ) 
C***BEGIN PROLOGUE  DSDOMN 
C***DATE WRITTEN   890404   (YYMMDD) 
C***REVISION DATE  890404   (YYMMDD) 
C***CATEGORY NO.  D2A4 
C***KEYWORDS  LIBRARY=SLATEC(SLAP), 
C             TYPE=DOUBLE PRECISION(SSDOMN-D), 
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition 
C***AUTHOR  Greenbaum, Anne, Courant Institute 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-300 
C             Livermore, CA 94550 (415) 423-3141 
C             seager@lll-crg.llnl.gov 
C***PURPOSE  Diagonally Scaled Orthomin Sparse Iterative Ax=b Solver. 
C            Routine to solve a general linear system  Ax = b using  
C            the Orthomin method with diagonal scaling. 
C***DESCRIPTION 
C *Usage: 
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
C     INTEGER  ITER, IERR, IUNIT, LENW, IWORK(10), LENIW 
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR 
C     DOUBLE PRECISION RWORK(7*N+3*N*NSAVE+NSAVE) 
C 
C     CALL DSDOMN(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL, 
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) 
C          
C *Arguments: 
C N      :IN       integer . 
C         Order of the Matrix. 
C B      :IN       Double Precision B(N). 
C         Right-hand side vector. 
C X      :INOUT    Double Precision X(N). 
C         On input X is your initial guess for solution vector. 
C         On output X is the final approximate solution. 
C NELT   :IN       integer . 
C         Number of Non-Zeros stored in A. 
C IA     :IN       integer  IA(NELT). 
C JA     :IN       integer  JA(NELT). 
C A      :IN       Double Precision A(NELT). 
C         These arrays should hold the matrix A in either the SLAP 
C         Triad format or the SLAP Column format.  See "LONG  
C         DESCRIPTION", below.  If the SLAP Triad format is chosen 
C         it is changed internally to the SLAP Column format. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all nonzero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C NSAVE  :IN       integer . 
C         Number of direction vectors to save and orthogonalize against. 
C ITOL   :IN       integer . 
C         Flag to indicate type of convergence criterion. 
C         If ITOL=1, iteration stops when the 2-norm of the residual  
C         divided by the 2-norm of the right-hand side is less than TOL. 
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  
C         residual divided by the 2-norm of M-inv times the right hand  
C         side is less than TOL, where M-inv is the inverse of the  
C         diagonal of A. 
C         ITOL=11 is often useful for checking and comparing different  
C         routines.  For this case, the user must supply the "exact"  
C         solution or a very accurate approximation (one with an error  
C         much less than TOL) through a common block, 
C         COMMON /SOLBLK/ SOLN( ) 
C         if ITOL=11, iteration stops when the 2-norm of the difference  
C         between the iterative approximation and the user-supplied 
C         solution divided by the 2-norm of the user-supplied solution  
C         is less than TOL. 
C TOL    :IN       Double Precision. 
C         Convergence criterion, as described above. 
C ITMAX  :IN       integer . 
C         Maximum number of iterations. 
C ITER   :OUT      integer . 
C         Number of iterations required to reach convergence, or  
C         ITMAX+1 if convergence criterion could not be achieved in  
C         ITMAX iterations. 
C ERR    :OUT      Double Precision. 
C         Error estimate of error in final approximate solution, as  
C         defined by ITOL. 
C IERR   :OUT      integer . 
C         Return error flag. 
C           IERR = 0 => All went well. 
C           IERR = 1 => Insufficient storage allocated  
C                       for WORK or IWORK. 
C           IERR = 2 => Method failed to converge in  
C                       ITMAX steps. 
C           IERR = 3 => Error in user input.  Check input 
C                       value of N, ITOL. 
C           IERR = 4 => User error tolerance set too tight. 
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded. 
C           IERR = 5 => Preconditioning matrix, M,  is not  
C                       Positive Definite.  $(r,z) < 0.0$. 
C           IERR = 6 => Breakdown of method detected. 
C                       $(p,Ap) < epsilon**2$. 
C IUNIT  :IN       integer . 
C         Unit number on which to write the error at each iteration,  
C         if this is desired for monitoring convergence.  If unit  
C         number is 0, no writing will occur. 
C RWORK  :WORK     Double Precision RWORK(LENW). 
C         Double Precision array used for workspace. 
C LENW   :IN       integer . 
C         Length of the double precision workspace, RWORK.   
C         LENW >= 7*N+NSAVE*(3*N+1). 
C IWORK  :WORK     integer  IWORK(LENIW). 
C         Used to hold pointers into the RWORK array. 
C LENIW  :IN       integer . 
C         Length of the double precision workspace, RWORK.  LENW >= 10. 
C          
C *Description: 
C       This routine  is simply a driver  for  the DOMN routine.  It 
C       calls the DSDS  routine  to set  up the  preconditioning and 
C       then   calls DOMN with the   appropriate   MATVEC and MSOLVE 
C       routines. 
C 
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix 
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP 
C       Column format.  The user can hand this routine either of the 
C       of these data structures and SLAP  will figure out  which on 
C       is being used and act accordingly. 
C        
C       =================== S L A P Triad format =================== 
C 
C       In  this   format only the  non-zeros are  stored.  They may 
C       appear  in *ANY* order.   The user  supplies three arrays of 
C       length NELT, where  NELT  is the number  of non-zeros in the 
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero 
C       the  user puts   the row  and  column index   of that matrix 
C       element in the IA and JA arrays.  The  value of the non-zero 
C       matrix  element is  placed in  the corresponding location of 
C       the A  array.  This is  an extremely easy data  structure to 
C       generate.  On  the other hand it  is  not too  efficient  on 
C       vector  computers   for the  iterative  solution  of  linear 
C       systems.  Hence, SLAP  changes this input  data structure to 
C       the SLAP   Column  format for the  iteration (but   does not 
C       change it back). 
C        
C       Here is an example of the  SLAP Triad   storage format for a 
C       5x5 Matrix.  Recall that the entries may appear in any order. 
C 
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. 
C                              1  2  3  4  5  6  7  8  9 10 11 
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       =================== S L A P Column format ================== 
C       This routine requires  that the  matrix  A be  stored in the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C        
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a  
C       column): 
C        
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C       1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C        
C *Precision:           Double Precision 
C *Side Effects: 
C       The SLAP Triad format (IA, JA, A)  is modified internally to 
C       be the   SLAP Column format.    See  the "LONG DESCRIPTION", 
C       below. 
C        
C *See Also: 
C         DOMN, DSLUOM 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  DS2Y, DCHKW, DSDS, DOMN, DSMV, DSDI 
C***END PROLOGUE  DSDOMN 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
      INTEGER  ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW 
      DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, RWORK(LENW) 
      EXTERNAL DSMV, DSDI 
      PARAMETER (LOCRB=1, LOCIB=11) 
C 
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format. 
C***FIRST EXECUTABLE STATEMENT  DSDOMN 
      iunit=0
      IERR = 0 
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
      CALL DS2Y( N, NELT, IA, JA, A, ISYM ) 
C          
C         Set up the workspace.  Compute the inverse of the  
C         diagonal of the matrix. 
      LOCIW = LOCIB 
C 
      LOCDIN = LOCRB 
      LOCR = LOCDIN + N 
      LOCZ = LOCR + N 
      LOCP = LOCZ + N 
      LOCAP = LOCP + N*(NSAVE+1) 
      LOCEMA = LOCAP + N*(NSAVE+1) 
      LOCDZ = LOCEMA + N*(NSAVE+1) 
      LOCCSA = LOCDZ + N 
      LOCW = LOCCSA + NSAVE 
C 
C         Check the workspace allocations. 
      CALL DCHKW('DSDOMN',LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR ) 
      IF( IERR.NE.0 ) RETURN 
C 
      IWORK(4) = LOCDIN 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
C 
      CALL DSDS(N, NELT, IA, JA, A, ISYM, RWORK(LOCDIN)) 
C          
C         Perform the Diagonally Scaled Orthomin iteration algorithm. 
      CALL DOMN(N, B, X, NELT, IA, JA, A, ISYM, DSMV,  
     $     DSDI, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP), RWORK(LOCAP), 
     $     RWORK(LOCEMA), RWORK(LOCDZ), RWORK(LOCCSA), 
     $     RWORK, IWORK ) 
      RETURN 
C------------- LAST LINE OF DSDOMN FOLLOWS ---------------------------- 
      END 
*DECK DSLUOM 
      SUBROUTINE DSLUOM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     RWORK, LENW, IWORK, LENIW ) 
C***BEGIN PROLOGUE  DSLUOM 
C***DATE WRITTEN   890404   (YYMMDD) 
C***REVISION DATE  890404   (YYMMDD) 
C***CATEGORY NO.  D2A4 
C***KEYWORDS  LIBRARY=SLATEC(SLAP), 
C             TYPE=DOUBLE PRECISION(SSLUOM-D), 
C             Non-Symmetric Linear system, Sparse,  
C             Iterative incomplete LU Precondition 
C***AUTHOR  Greenbaum, Anne, Courant Institute 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-300 
C             Livermore, CA 94550 (415) 423-3141 
C             seager@lll-crg.llnl.gov 
C***PURPOSE  Incomplete LU Orthomin Sparse Iterative Ax=b Solver. 
C            Routine to solve a general linear system  Ax = b  using  
C            the Orthomin method with Incomplete LU decomposition. 
C***DESCRIPTION 
C *Usage: 
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
C     INTEGER  ITER, IERR, IUNIT, LENW, IWORK(NEL+NU+4*N+2), LENIW 
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR 
C     DOUBLE PRECISION RWORK(NEL+NU+7*N+3*N*NSAVE+NSAVE) 
C 
C     CALL DSLUOM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL, 
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) 
C          
C *Arguments: 
C N      :IN       integer . 
C         Order of the matrix. 
C B      :IN       Double Precision B(N). 
C         Right-hand side vector. 
C X      :INOUT    Double Precision X(N). 
C         On input X is your initial guess for solution vector. 
C         On output X is the final approximate solution. 
C NELT   :IN       integer . 
C         Number of Non-Zeros stored in A. 
C IA     :INOUT    integer  IA(NELT). 
C JA     :INOUT    integer  JA(NELT). 
C A      :INOUT    Double Precision A(NELT). 
C         These arrays should hold the matrix A in either the SLAP 
C         Triad format or the SLAP Column format.  See "LONG  
C         DESCRIPTION", below.  If the SLAP Triad format is chosen 
C         it is changed internally to the SLAP Column format. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all nonzero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C NSAVE  :IN       integer . 
C         Number of direction vectors to save and orthogonalize against. 
C ITOL   :IN       integer . 
C         Flag to indicate type of convergence criterion. 
C         If ITOL=1, iteration stops when the 2-norm of the residual  
C         divided by the 2-norm of the right-hand side is less than TOL. 
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  
C         residual divided by the 2-norm of M-inv times the right hand  
C         side is less than TOL, where M-inv is the inverse of the  
C         diagonal of A. 
C         ITOL=11 is often useful for checking and comparing different  
C         routines.  For this case, the user must supply the "exact"  
C         solution or a very accurate approximation (one with an error  
C         much less than TOL) through a common block, 
C                     COMMON /SOLBLK/ SOLN(1) 
C         if ITOL=11, iteration stops when the 2-norm of the difference  
C         between the iterative approximation and the user-supplied 
C         solution divided by the 2-norm of the user-supplied solution  
C         is less than TOL.  Note that this requires the user to set up 
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine.  
C         The routine with this declaration should be loaded before the 
C         stop test so that the correct length is used by the loader.   
C         This procedure is not standard Fortran and may not work  
C         correctly on your system (although it has worked on every 
C         system the authors have tried).  If ITOL is not 11 then this 
C         common block is indeed standard Fortran. 
C TOL    :IN       Double Precision. 
C         Convergence criterion, as described above. 
C ITMAX  :IN       integer . 
C         Maximum number of iterations. 
C ITER   :OUT      integer . 
C         Number of iterations required to reach convergence, or  
C         ITMAX+1 if convergence criterion could not be achieved in  
C         ITMAX iterations. 
C ERR    :OUT      Double Precision. 
C         Error estimate of error in final approximate solution, as  
C         defined by ITOL. 
C IERR   :OUT      integer . 
C         Return error flag. 
C           IERR = 0 => All went well. 
C           IERR = 1 => Insufficient storage allocated  
C                       for WORK or IWORK. 
C           IERR = 2 => Method failed to converge in  
C                       ITMAX steps. 
C           IERR = 3 => Error in user input.  Check input 
C                       value of N, ITOL. 
C           IERR = 4 => User error tolerance set too tight. 
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded. 
C           IERR = 5 => Preconditioning matrix, M,  is not  
C                       Positive Definite.  $(r,z) < 0.0$. 
C           IERR = 6 => Breakdown of the method detected. 
C                       $(p,Ap) < epsilon**2$. 
C           IERR = 7 => Incomplete factorization broke down 
C                       and was fudged.  Resulting preconditioning may 
C                       be less than the best. 
C IUNIT  :IN       integer . 
C         Unit number on which to write the error at each iteration,  
C         if this is desired for monitoring convergence.  If unit  
C         number is 0, no writing will occur. 
C RWORK  :WORK     Double Precision RWORK(LENW). 
C         Double Precision array used for workspace.  NL is the  
C         number of non- 
C         zeros in the lower triangle of the matrix (including the 
C         diagonal).  NU is the number of nonzeros in the upper 
C         triangle of the matrix (including the diagonal). 
C LENW   :IN       integer . 
C         Length of the double precision workspace, RWORK.   
C         LENW >= NL+NU+4*N+NSAVE*(3*N+1) 
C IWORK  :WORK     integer  IWORK(LENIW) 
C         integer  array used for workspace.  NL is the number of non- 
C         zeros in the lower triangle of the matrix (including the 
C         diagonal).  NU is the number of nonzeros in the upper 
C         triangle of the matrix (including the diagonal). 
C         Upon return the following locations of IWORK hold information 
C         which may be of use to the user: 
C         IWORK(9)  Amount of integer  workspace actually used. 
C         IWORK(10) Amount of Double Precision workspace actually used. 
C LENIW  :IN       integer . 
C         Length of the double precision workspace, RWORK.   
C         LENW > NL+NU+4*N+12. 
C 
C *Description: 
C       This routine is  simply a driver  for  the DOMN routine.  It 
C       calls the DSILUS routine  to set  up the preconditioning and 
C       then  calls   DOMN  with the appropriate  MATVEC  and MSOLVE 
C       routines. 
C        
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix 
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP 
C       Column format.  The user can hand this routine either of the 
C       of these data structures and SLAP  will figure out  which on 
C       is being used and act accordingly. 
C        
C       =================== S L A P Triad format =================== 
C 
C       This routine requires that the  matrix A be   stored in  the 
C       SLAP  Triad format.  In  this format only the non-zeros  are 
C       stored.  They may appear in  *ANY* order.  The user supplies 
C       three arrays of  length NELT, where  NELT is  the number  of 
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For 
C       each non-zero the user puts the row and column index of that 
C       matrix element  in the IA and  JA arrays.  The  value of the 
C       non-zero   matrix  element is  placed  in  the corresponding 
C       location of the A array.   This is  an  extremely  easy data 
C       structure to generate.  On  the  other hand it   is  not too 
C       efficient on vector computers for  the iterative solution of 
C       linear systems.  Hence,   SLAP changes   this  input    data 
C       structure to the SLAP Column format  for  the iteration (but 
C       does not change it back). 
C        
C       Here is an example of the  SLAP Triad   storage format for a 
C       5x5 Matrix.  Recall that the entries may appear in any order. 
C 
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. 
C                              1  2  3  4  5  6  7  8  9 10 11 
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C        
C       =================== S L A P Column format ================== 
C       This routine requires  that the  matrix  A be  stored in the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C        
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a  
C       column): 
C        
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C       1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C *Precision:           Double Precision 
C *Side Effects: 
C       The SLAP Triad format (IA, JA,  A) is modified internally to 
C       be the  SLAP  Column format.  See  the   "LONG DESCRIPTION", 
C       below. 
C        
C *See Also: 
C         DOMN, DSDOMN 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  DS2Y, DCHKW, DSILUS, DOMN, DSMV, DSLUI 
C***END PROLOGUE  DSLUOM 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
      INTEGER  ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW 
      DOUBLE PRECISION B(N), X(N), A(N), RWORK(LENW) 
      EXTERNAL DSMV, DSLUI 
      PARAMETER (LOCRB=1, LOCIB=11) 
      iunit=0
C 
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format. 
C***FIRST EXECUTABLE STATEMENT  DSLUOM 
      IERR = 0 
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
      CALL DS2Y( N, NELT, IA, JA, A, ISYM ) 
C 
C         Count number of Non-Zero elements preconditioner ILU matrix. 
C         Then set up the work arrays. 
      NL = 0 
      NU = 0 
      DO 20 ICOL = 1, N 
C         Don't count diagonal. 
         JBGN = JA(ICOL)+1 
         JEND = JA(ICOL+1)-1 
         IF( JBGN.LE.JEND ) THEN 
CVD$ NOVECTOR 
            DO 10 J = JBGN, JEND 
               IF( IA(J).GT.ICOL ) THEN 
                  NL = NL + 1 
                  IF( ISYM.NE.0 ) NU = NU + 1 
               ELSE 
                  NU = NU + 1 
               ENDIF 
 10         CONTINUE 
         ENDIF 
 20   CONTINUE 
C          
      LOCIL = LOCIB 
      LOCJL = LOCIL + N+1 
      LOCIU = LOCJL + NL 
      LOCJU = LOCIU + NU 
      LOCNR = LOCJU + N+1 
      LOCNC = LOCNR + N 
      LOCIW = LOCNC + N 
C 
      LOCL   = LOCRB 
      LOCDIN = LOCL + NL 
      LOCU   = LOCDIN + N 
      LOCR   = LOCU + NU 
      LOCZ   = LOCR + N 
      LOCP   = LOCZ + N 
      LOCAP  = LOCP + N*(NSAVE+1) 
      LOCEMA = LOCAP + N*(NSAVE+1) 
      LOCDZ  = LOCEMA + N*(NSAVE+1) 
      LOCCSA = LOCDZ + N 
      LOCW   = LOCCSA + NSAVE 
C 
C         Check the workspace allocations. 
      CALL DCHKW('DSLUOM',LOCIW,LENIW, LOCW, LENW, IERR, ITER, ERR ) 
      IF( IERR.NE.0 ) RETURN 
C 
      IWORK(1) = LOCIL 
      IWORK(2) = LOCJL 
      IWORK(3) = LOCIU 
      IWORK(4) = LOCJU 
      IWORK(5) = LOCL 
      IWORK(6) = LOCDIN 
      IWORK(7) = LOCU 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
C 
C         Compute the Incomplete LU decomposition. 
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL), 
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU), 
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) ) 
C          
C         Perform the incomplete LU preconditioned OrthoMin algorithm. 
      CALL DOMN(N, B, X, NELT, IA, JA, A, ISYM, DSMV, 
     $     DSLUI, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP), RWORK(LOCAP), 
     $     RWORK(LOCEMA), RWORK(LOCDZ), RWORK(LOCCSA), 
     $     RWORK, IWORK ) 
      RETURN 
      END 
*DECK ISDOMN 
      FUNCTION ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     R, Z, P, AP, EMAP, DZ, CSAV, 
     $     RWORK, IWORK, AK, BNRM, SOLNRM) 
C***BEGIN PROLOGUE  ISDOMN 
C***REFER TO  DOMN, DSDOMN, DSLUOM 
C***DATE WRITTEN   890404   (YYMMDD) 
C***REVISION DATE  890404   (YYMMDD) 
C***CATEGORY NO.  D2A4 
C***KEYWORDS  LIBRARY=SLATEC(SLAP), 
C             TYPE=DOUBLE PRECISION(ISDOMN-D), 
C             Non-Symmetric Linear system, Sparse,  
C             Iterative Precondition, Stop Test, Orthomin 
C***AUTHOR  Greenbaum, Anne, Courant Institute 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-300 
C             Livermore, CA 94550 (415) 423-3141 
C             seager@lll-crg.llnl.gov 
C***PURPOSE  Preconditioned Orthomin Sparse Stop Test. 
C            This routine calculates the stop  test for the Orthomin 
C            iteration  scheme.  It returns a  nonzero if the  error 
C            estimate (the type of  which is  determined by ITOL) is 
C            less than the user specified tolerance TOL. 
C***DESCRIPTION 
C *Usage: 
C     INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
C     INTEGER   ITER, IERR, IUNIT, IWORK(nelt) 
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N) 
C     DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE) 
C     DOUBLE PRECISION DZ(N), CSAV(NSAVE), RWORK(nelt), AK 
C     DOUBLE PRECISION BNRM, SOLNRM 
C     EXTERNAL MSOLVE 
C 
C     IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE, 
C    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,  
C    $     EMAP, DZ, CSAV, RWORK, IWORK, AK, BNRM, SOLNRM) 
C    $     .NE.0 ) THEN ITERATION CONVERGED 
C 
C *Arguments: 
C N      :IN       integer . 
C         Order of the matrix. 
C B      :IN       Double Precision B(N). 
C         Right-hand side vector. 
C X      :IN       Double Precision X(N). 
C         On input X is your initial guess for solution vector. 
C         On output X is the final approximate solution. 
C NELT   :IN       integer . 
C         Number of Non-Zeros stored in A. 
C IA     :IN       integer  IA(NELT). 
C JA     :IN       integer  JA(NELT). 
C A      :IN       Double Precision A(NELT). 
C         These arrays should hold the matrix A in either the SLAP 
C         Triad format or the SLAP Column format.  See "LONG 
C         DESCRIPTION" in the DSDOMN or DSLUOM. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all nonzero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C MSOLVE :EXT      External. 
C         Name of a routine which solves a linear system MZ = R for 
C         Z given R with the preconditioning matrix M (M is supplied via 
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  
C         be declared external in the calling program.  The calling  
C         sequence to MSOLVE is: 
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C         Where N is the number of unknowns, R is the right-hand side  
C         vector, and Z is the solution upon return.  RWORK is a  
C         double precision  
C         array that can be used to pass necessary preconditioning  
C         information and/or workspace to MSOLVE.  IWORK is an integer   
C         work array for the same purpose as RWORK. 
C NSAVE  :IN       integer . 
C         Number of direction vectors to save and orthogonalize against. 
C ITOL   :IN       integer . 
C         Flag to indicate type of convergence criterion. 
C         If ITOL=1, iteration stops when the 2-norm of the residual  
C         divided by the 2-norm of the right-hand side is less than TOL. 
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  
C         residual divided by the 2-norm of M-inv times the right hand  
C         side is less than TOL, where M-inv is the inverse of the  
C         diagonal of A. 
C         ITOL=11 is often useful for checking and comparing different  
C         routines.  For this case, the user must supply the "exact"  
C         solution or a very accurate approximation (one with an error  
C         much less than TOL) through a common block, 
C                     COMMON /SOLBLK/ SOLN(1) 
C         if ITOL=11, iteration stops when the 2-norm of the difference  
C         between the iterative approximation and the user-supplied 
C         solution divided by the 2-norm of the user-supplied solution  
C         is less than TOL.  Note that this requires the user to set up 
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine.  
C         The routine with this declaration should be loaded before the 
C         stop test so that the correct length is used by the loader.   
C         This procedure is not standard Fortran and may not work  
C         correctly on your system (although it has worked on every 
C         system the authors have tried).  If ITOL is not 11 then this 
C         common block is indeed standard Fortran. 
C TOL    :IN       Double Precision. 
C         Convergence criterion, as described above. 
C ITMAX  :IN       integer . 
C         Maximum number of iterations. 
C ITER   :IN       integer . 
C         Number of iterations required to reach convergence, or  
C         ITMAX+1 if convergence criterion could not be achieved in  
C         ITMAX iterations. 
C ERR    :OUT      Double Precision. 
C         Error estimate of error in final approximate solution, as  
C         defined by ITOL. 
C IERR   :OUT      integer . 
C         Error flag.  IERR is set to 3 if ITOL is not on of the  
C         acceptable values, see above.  
C IUNIT  :IN       integer . 
C         Unit number on which to write the error at each iteration,  
C         if this is desired for monitoring convergence.  If unit  
C         number is 0, no writing will occur. 
C R      :IN       Double Precision R(N). 
C         The residual R = B-AX. 
C Z      :WORK     Double Precision Z(N). 
C P      :IN       Double Precision P(N,0:NSAVE). 
C         Workspace used to hold the conjugate direction vector(s). 
C AP     :IN       Double Precision AP(N,0:NSAVE). 
C         Workspace used to hold the matrix A times the P vector(s). 
C EMAP   :IN       Double Precision EMAP(N,0:NSAVE). 
C         Workspace used to hold M-inv times the AP vector(s). 
C DZ     :WORK     Double Precision DZ(N). 
C         Workspace. 
C CSAV   :DUMMY    Double Precision CSAV(NSAVE) 
C         Reserved for future use. 
C RWORK  :WORK     Double Precision RWORK(USER DEFINED). 
C         Double Precision array that can be used for workspace in  
C         MSOLVE. 
C IWORK  :WORK     integer  IWORK(USER DEFINED). 
C         integer  array that can be used for workspace in MSOLVE. 
C AK     :IN       Double Precision. 
C         Current iterate BiConjugate Gradient iteration parameter. 
C 
C *Function Return Values: 
C       0 : Error estimate (determined by ITOL) is *NOT* less than the  
C           specified tolerance, TOL.  The iteration must continue. 
C       1 : Error estimate (determined by ITOL) is less than the  
C           specified tolerance, TOL.  The iteration can be considered 
C           complete. 
C 
C *Precision:           Double Precision 
C *See Also: 
C         DOMN, DSDOMN, DSLUOM 
C 
C *Cautions: 
C     This routine will attempt to write to the fortran logical output  
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that 
C     this  logical  unit  must  be  attached  to  a  file or terminal 
C     before calling this routine with a non-zero  value  for   IUNIT. 
C     This routine does not check for the validity of a non-zero IUNIT 
C     unit number. 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  MSOLVE, DNRM2 
C***COMMON BLOCKS    SOLBLK 
C***END PROLOGUE  ISDOMN 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX 
      INTEGER   ITER, IUNIT, iwork(nelt) 
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N) 
      DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE) 
      DOUBLE PRECISION DZ(N), CSAV(NSAVE), rwork(nelt) 
      EXTERNAL MSOLVE 
      COMMON /SOLBLK/ SOLN(1) 
C          
C***FIRST EXECUTABLE STATEMENT  ISDOMN 
      ISDOMN = 0 
C          
      iunit=0
      IF( ITOL.EQ.1 ) THEN 
C         err = ||Residual||/||RightHandSide|| (2-Norms). 
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1) 
         ERR = DNRM2(N, R, 1)/BNRM 
      ELSE IF( ITOL.EQ.2 ) THEN 
C                  -1              -1 
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms). 
         IF(ITER .EQ. 0) THEN 
            CALL MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
            BNRM = DNRM2(N, DZ, 1) 
         ENDIF 
         ERR = DNRM2(N, Z, 1)/BNRM 
      ELSE IF( ITOL.EQ.11 ) THEN 
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms). 
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1) 
         DO 10 I = 1, N 
            DZ(I) = X(I) - SOLN(I) 
 10      CONTINUE 
         ERR = DNRM2(N, DZ, 1)/SOLNRM 
      ELSE 
C 
C         If we get here ITOL is not one of the acceptable values. 
         ERR = 1.0E10 
         IERR = 3 
      ENDIF 
C          
      IF(IUNIT .NE. 0) THEN 
         IF( ITER.EQ.0 ) THEN 
            WRITE(IUNIT,1000) NSAVE, N, ITOL 
         ENDIF 
         WRITE(IUNIT,1010) ITER, ERR, AK 
      ENDIF 
      IF(ERR .LE. TOL) ISDOMN = 1 
C          
      RETURN 
 1000 FORMAT(' Preconditioned Orthomin(',I3,') for ', 
     $     'N, ITOL = ',I5, I5, 
     $     /' ITER','   Error Estimate','            Alpha') 
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7) 
C------------- LAST LINE OF ISDOMN FOLLOWS ---------------------------- 
      END 
 
 
 

 



      DOUBLE PRECISION FUNCTION D1MACH(I) 
      INTEGER  I 
C 
C  DOUBLE-PRECISION MACHINE CONSTANTS 
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE. 
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE. 
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING. 
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING. 
C  D1MACH( 5) = LOG10(B) 
C 
      INTEGER  SMALL(2) 
      INTEGER  LARGE(2) 
      INTEGER  RIGHT(2) 
      INTEGER  DIVER(2) 
      INTEGER  LOG10(2) 
      INTEGER  SC, CRAY1(38), J 
      COMMON /D9MACH/ CRAY1 
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC 
      DOUBLE PRECISION DMACH(5) 
      EQUIVALENCE (DMACH(1),SMALL(1)) 
      EQUIVALENCE (DMACH(2),LARGE(1)) 
      EQUIVALENCE (DMACH(3),RIGHT(1)) 
      EQUIVALENCE (DMACH(4),DIVER(1)) 
      EQUIVALENCE (DMACH(5),LOG10(1)) 
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES. 
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF 
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR 
C  MANY MACHINES YET. 
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1 
C  ON THE NEXT LINE 
      DATA SC/0/ 
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW. 
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY 
C          mail netlib@research.bell-labs.com 
C          send old1mach from blas 
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com. 
C 
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES. 
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 / 
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 / 
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 / 
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 / 
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/ 
C 
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING 
C     32-BIT INTEGER S. 
C      DATA SMALL(1),SMALL(2) /    8388608,           0 / 
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 / 
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 / 
C      DATA DIVER(1),DIVER(2) /  620756992,           0 / 
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/ 
C 
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. 
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 / 
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 / 
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 / 
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 / 
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/ 
C 
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES. 
      IF (SC .NE. 987) THEN 
         DMACH(1) = 1.D13 
         IF (      SMALL(1) .EQ. 1117925532 
     $       .AND. SMALL(2) .EQ. -448790528) THEN 
*           *** IEEE BIG ENDIAN *** 
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2146435071 
            LARGE(2) = -1 
            RIGHT(1) = 1017118720 
            RIGHT(2) = 0 
            DIVER(1) = 1018167296 
            DIVER(2) = 0 
            LOG10(1) = 1070810131 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(2) .EQ. 1117925532 
     $       .AND. SMALL(1) .EQ. -448790528) THEN 
*           *** IEEE LITTLE ENDIAN *** 
            SMALL(2) = 1048576 
            SMALL(1) = 0 
            LARGE(2) = 2146435071 
            LARGE(1) = -1 
            RIGHT(2) = 1017118720 
            RIGHT(1) = 0 
            DIVER(2) = 1018167296 
            DIVER(1) = 0 
            LOG10(2) = 1070810131 
            LOG10(1) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. -2065213935 
     $       .AND. SMALL(2) .EQ. 10752) THEN 
*               *** VAX WITH D_FLOATING *** 
            SMALL(1) = 128 
            SMALL(2) = 0 
            LARGE(1) = -32769 
            LARGE(2) = -1 
            RIGHT(1) = 9344 
            RIGHT(2) = 0 
            DIVER(1) = 9472 
            DIVER(2) = 0 
            LOG10(1) = 546979738 
            LOG10(2) = -805796613 
         ELSE IF ( SMALL(1) .EQ. 1267827943 
     $       .AND. SMALL(2) .EQ. 704643072) THEN 
*               *** IBM MAINFRAME *** 
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2147483647 
            LARGE(2) = -1 
            RIGHT(1) = 856686592 
            RIGHT(2) = 0 
            DIVER(1) = 873463808 
            DIVER(2) = 0 
            LOG10(1) = 1091781651 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. 1120022684 
     $       .AND. SMALL(2) .EQ. -448790528) THEN 
*           *** CONVEX C-1 *** 
            SMALL(1) = 1048576 
            SMALL(2) = 0 
            LARGE(1) = 2147483647 
            LARGE(2) = -1 
            RIGHT(1) = 1019215872 
            RIGHT(2) = 0 
            DIVER(1) = 1020264448 
            DIVER(2) = 0 
            LOG10(1) = 1072907283 
            LOG10(2) = 1352628735 
         ELSE IF ( SMALL(1) .EQ. 815547074 
     $       .AND. SMALL(2) .EQ. 58688) THEN 
*           *** VAX G-FLOATING *** 
            SMALL(1) = 16 
            SMALL(2) = 0 
            LARGE(1) = -32769 
            LARGE(2) = -1 
            RIGHT(1) = 15552 
            RIGHT(2) = 0 
            DIVER(1) = 15568 
            DIVER(2) = 0 
            LOG10(1) = 1142112243 
            LOG10(2) = 2046775455 
         ELSE 
            DMACH(2) = 1.D27 + 1 
            DMACH(3) = 1.D27 
            LARGE(2) = LARGE(2) - RIGHT(2) 
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN 
               CRAY1(1) = 67291416 
               DO 10 J = 1, 20 
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J) 
 10               CONTINUE 
               CRAY1(22) = CRAY1(21) + 321322 
               DO 20 J = 22, 37 
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J) 
 20               CONTINUE 
               IF (CRAY1(38) .EQ. SMALL(1)) THEN 
*                  *** CRAY *** 
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0) 
                  SMALL(2) = 0 
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215) 
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214) 
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0) 
                  RIGHT(2) = 0 
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0) 
                  DIVER(2) = 0 
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215) 
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388) 
               ELSE 
                  WRITE(*,9000) 
                  STOP 779 
                  END IF 
            ELSE 
               WRITE(*,9000) 
               STOP 779 
               END IF 
            END IF 
         SC = 987 
         END IF 
*    SANITY CHECK 
      IF (DMACH(4) .GE. 1.0D0) STOP 778 
      IF (I .LT. 1 .OR. I .GT. 5) THEN 
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.' 
         STOP 
         END IF 
 
 
c	Modifica Pietro 
	d1mach=1e-90 
c      D1MACH = DMACH(I) 
      RETURN 
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/ 
     $' appropriate for your machine.') 
* /* Standard C source for D1MACH -- remove the * in column 1 */ 
*#include <stdio.h> 
*#include <float.h> 
*#include <math.h> 
*double d1mach_(long *i) 
*{ 
*	switch(*i){ 
*	  case 1: return DBL_MIN; 
*	  case 2: return DBL_MAX; 
*	  case 3: return DBL_EPSILON/FLT_RADIX; 
*	  case 4: return DBL_EPSILON; 
*	  case 5: return log10((double)FLT_RADIX); 
*	  } 
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i); 
*	exit(1); return 0; /* some compilers demand return values */ 
*} 
      END 
      SUBROUTINE I1MCRY(A, A1, B, C, D) 
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES **** 
      INTEGER  A, A1, B, C, D 
      A1 = 16777216*B + C 
      A = 16777216*A1 + D 
      END 
 
 
 
 
*DECK DS2Y 
      SUBROUTINE DS2Y (N, NELT, IA, JA, A, ISYM) 
C***BEGIN PROLOGUE  DS2Y 
C***PURPOSE  SLAP Triad to SLAP Column Format Converter. 
C            Routine to convert from the SLAP Triad to SLAP Column 
C            format. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D1B9 
C***TYPE      DOUBLE PRECISION (SS2Y-S, DS2Y-D) 
C***KEYWORDS  LINEAR SYSTEM, SLAP SPARSE 
C***AUTHOR  Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM 
C     DOUBLE PRECISION A(NELT) 
C 
C     CALL DS2Y( N, NELT, IA, JA, A, ISYM ) 
C 
C *Arguments: 
C N      :IN       integer  
C         Order of the Matrix. 
C NELT   :IN       integer . 
C         Number of non-zeros stored in A. 
C IA     :INOUT    integer  IA(NELT). 
C JA     :INOUT    integer  JA(NELT). 
C A      :INOUT    Double Precision A(NELT). 
C         These arrays should hold the matrix A in either the SLAP 
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is used, this format is 
C         translated to the SLAP Column format by this routine. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all non-zero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored. 
C 
C *Description: 
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix 
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP 
C       Column format.  The user can hand this routine either of the 
C       of these data structures.  If the SLAP Triad format is give 
C       as input then this routine transforms it into SLAP Column 
C       format.  The way this routine tells which format is given as 
C       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we 
C       have the SLAP Column format.  If that equality does not hold 
C       then it is assumed that the IA, JA, A arrays contain the 
C       SLAP Triad format. 
C 
C       =================== S L A P Triad format =================== 
C       This routine requires that the  matrix A be   stored in  the 
C       SLAP  Triad format.  In  this format only the non-zeros  are 
C       stored.  They may appear in  *ANY* order.  The user supplies 
C       three arrays of  length NELT, where  NELT is  the number  of 
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For 
C       each non-zero the user puts the row and column index of that 
C       matrix element  in the IA and  JA arrays.  The  value of the 
C       non-zero   matrix  element is  placed  in  the corresponding 
C       location of the A array.   This is  an  extremely  easy data 
C       structure to generate.  On  the  other hand it   is  not too 
C       efficient on vector computers for  the iterative solution of 
C       linear systems.  Hence,   SLAP changes   this  input    data 
C       structure to the SLAP Column format  for  the iteration (but 
C       does not change it back). 
C 
C       Here is an example of the  SLAP Triad   storage format for a 
C       5x5 Matrix.  Recall that the entries may appear in any order. 
C 
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left. 
C                              1  2  3  4  5  6  7  8  9 10 11 
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       =================== S L A P Column format ================== 
C 
C       This routine  requires that  the matrix A  be stored in  the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C 
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column): 
C 
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  QS2I1D 
C***REVISION HISTORY  (YYMMDD) 
C   871119  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF) 
C   920511  Added complete declaration section.  (WRB) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DS2Y 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT) 
      INTEGER  IA(NELT), JA(NELT) 
C     .. Local Scalars .. 
      DOUBLE PRECISION TEMP 
      INTEGER  I, IBGN, ICOL, IEND, ITEMP, J 
C     .. External Subroutines .. 
      EXTERNAL QS2I1D 
C***FIRST EXECUTABLE STATEMENT  DS2Y 
C 
C         Check to see if the (IA,JA,A) arrays are in SLAP Column 
C         format.  If it's not then transform from SLAP Triad. 
C 
      IF( JA(N+1).EQ.NELT+1 ) RETURN 
C 
C         Sort into ascending order by COLUMN (on the ja array). 
C         This will line up the columns. 
C 
      CALL QS2I1D( JA, IA, A, NELT, 1 ) 
C 
C         Loop over each column to see where the column indices change 
C         in the column index array ja.  This marks the beginning of the 
C         next column. 
C 
CVD$R NOVECTOR 
      JA(1) = 1 
      DO 20 ICOL = 1, N-1 
         DO 10 J = JA(ICOL)+1, NELT 
            IF( JA(J).NE.ICOL ) THEN 
               JA(ICOL+1) = J 
               GOTO 20 
            ENDIF 
 10      CONTINUE 
 20   CONTINUE 
      JA(N+1) = NELT+1 
C 
C         Mark the n+2 element so that future calls to a SLAP routine 
C         utilizing the YSMP-Column storage format will be able to tell. 
C 
      JA(N+2) = 0 
C 
C         Now loop through the IA array making sure that the diagonal 
C         matrix element appears first in the column.  Then sort the 
C         rest of the column in ascending order. 
C 
      DO 70 ICOL = 1, N 
         IBGN = JA(ICOL) 
         IEND = JA(ICOL+1)-1 
         DO 30 I = IBGN, IEND 
            IF( IA(I).EQ.ICOL ) THEN 
C 
C              Swap the diagonal element with the first element in the 
C              column. 
C 
               ITEMP = IA(I) 
               IA(I) = IA(IBGN) 
               IA(IBGN) = ITEMP 
               TEMP = A(I) 
               A(I) = A(IBGN) 
               A(IBGN) = TEMP 
               GOTO 40 
            ENDIF 
 30      CONTINUE 
 40      IBGN = IBGN + 1 
         IF( IBGN.LT.IEND ) THEN 
            DO 60 I = IBGN, IEND 
               DO 50 J = I+1, IEND 
                  IF( IA(I).GT.IA(J) ) THEN 
                     ITEMP = IA(I) 
                     IA(I) = IA(J) 
                     IA(J) = ITEMP 
                     TEMP = A(I) 
                     A(I) = A(J) 
                     A(J) = TEMP 
                  ENDIF 
 50            CONTINUE 
 60         CONTINUE 
         ENDIF 
 70   CONTINUE 
      RETURN 
C------------- LAST LINE OF DS2Y FOLLOWS ---------------------------- 
      END 
 
 
*DECK DCHKW 
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) 
C***BEGIN PROLOGUE  DCHKW 
C***SUBSIDIARY 
C***PURPOSE  SLAP WORK/IWORK Array Bounds Checker. 
C            This routine checks the work array lengths and interfaces 
C            to the SLATEC error handler if a problem is found. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  R2 
C***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D) 
C***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING 
C***AUTHOR  Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     CHARACTER*(*) NAME 
C     INTEGER  LOCIW, LENIW, LOCW, LENW, IERR, ITER 
C     DOUBLE PRECISION ERR 
C 
C     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR ) 
C 
C *Arguments: 
C NAME   :IN       Character*(*). 
C         Name of the calling routine.  This is used in the output 
C         message, if an error is detected. 
C LOCIW  :IN       integer . 
C         Location of the first free element in the integer  workspace 
C         array. 
C LENIW  :IN       integer . 
C         Length of the integer  workspace array. 
C LOCW   :IN       integer . 
C         Location of the first free element in the double precision 
C         workspace array. 
C LENRW  :IN       integer . 
C         Length of the double precision workspace array. 
C IERR   :OUT      integer . 
C         Return error flag. 
C               IERR = 0 => All went well. 
C               IERR = 1 => Insufficient storage allocated for 
C                           WORK or IWORK. 
C ITER   :OUT      integer . 
C         Set to zero on return. 
C ERR    :OUT      Double Precision. 
C         Set to the smallest positive magnitude if all went well. 
C         Set to a very large number if an error is detected. 
C 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  D1MACH, XERMSG 
C***REVISION HISTORY  (YYMMDD) 
C   880225  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI 
C           X3.9-1978.  (FNF) 
C   910506  Made subsidiary.  (FNF) 
C   920511  Added complete declaration section.  (WRB) 
C   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF) 
C***END PROLOGUE  DCHKW 
C     .. Scalar Arguments .. 
      DOUBLE PRECISION ERR 
      INTEGER  IERR, ITER, LENIW, LENW, LOCIW, LOCW 
      CHARACTER NAME*(*) 
C     .. Local Scalars .. 
      CHARACTER XERN1*8, XERN2*8, XERNAM*8 
C     .. External Functions .. 
      DOUBLE PRECISION D1MACH 
      EXTERNAL D1MACH 
C     .. External Subroutines .. 
      EXTERNAL XERMSG 
C***FIRST EXECUTABLE STATEMENT  DCHKW 
C 
C         Check the integer  workspace situation. 
C 
      IERR = 0 
      ITER = 0 
      ERR = D1MACH(1) 
      IF( LOCIW.GT.LENIW ) THEN 
         IERR = 1 
         ERR = D1MACH(2) 
         XERNAM = NAME 
         WRITE (XERN1, '(I8)') LOCIW 
         WRITE (XERN2, '(I8)') LENIW 
         CALL XERMSG ('SLATEC', 'DCHKW', 
     $      'In ' // XERNAM // ', INTEGER  work array too short.  ' // 
     $      'IWORK needs ' // XERN1 // '; have allocated ' // XERN2, 
     $      1, 1) 
      ENDIF 
C 
C         Check the Double Precision workspace situation. 
      IF( LOCW.GT.LENW ) THEN 
         IERR = 1 
         ERR = D1MACH(2) 
         XERNAM = NAME 
         WRITE (XERN1, '(I8)') LOCW 
         WRITE (XERN2, '(I8)') LENW 
         CALL XERMSG ('SLATEC', 'DCHKW', 
     $      'In ' // XERNAM // ', DOUBLE PRECISION work array too ' // 
     $      'short.  RWORK needs ' // XERN1 // '; have allocated ' // 
     $      XERN2, 1, 1) 
      ENDIF 
      RETURN 
C------------- LAST LINE OF DCHKW FOLLOWS ---------------------------- 
      END 
 
 
 
*DECK DSDS 
      SUBROUTINE DSDS (N, NELT, IA, JA, A, ISYM, DINV) 
C***BEGIN PROLOGUE  DSDS 
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up. 
C            Routine to compute the inverse of the diagonal of a matrix 
C            stored in the SLAP Column format. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D2E 
C***TYPE      DOUBLE PRECISION (SSDS-S, DSDS-D) 
C***KEYWORDS  DIAGONAL, SLAP SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM 
C     DOUBLE PRECISION A(NELT), DINV(N) 
C 
C     CALL DSDS( N, NELT, IA, JA, A, ISYM, DINV ) 
C 
C *Arguments: 
C N      :IN       integer . 
C         Order of the Matrix. 
C NELT   :IN       integer . 
C         Number of elements in arrays IA, JA, and A. 
C IA     :INOUT    integer  IA(NELT). 
C JA     :INOUT    integer  JA(NELT). 
C A      :INOUT    Double Precision A(NELT). 
C         These arrays should hold the matrix A in the SLAP Column 
C         format.  See "Description", below. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all non-zero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C DINV   :OUT      Double Precision DINV(N). 
C         Upon return this array holds 1./DIAG(A). 
C 
C *Description 
C       =================== S L A P Column format ================== 
C       This routine  requires that  the matrix A  be stored in  the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C 
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column): 
C 
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       With the SLAP  format  all  of  the   "inner  loops" of this 
C       routine should vectorize  on  machines with hardware support 
C       for vector   gather/scatter  operations.  Your compiler  may 
C       require a compiler directive to  convince it that  there are 
C       no  implicit  vector  dependencies.  Compiler directives for 
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are 
C       supplied with the standard SLAP distribution. 
C 
C 
C *Cautions: 
C       This routine assumes that the diagonal of A is all  non-zero 
C       and that the operation DINV = 1.0/DIAG(A) will not underflow 
C       or overflow.    This  is done so that the  loop  vectorizes. 
C       Matrices  with zero or near zero or very  large entries will 
C       have numerical difficulties  and  must  be fixed before this 
C       routine is called. 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   890404  DATE WRITTEN 
C   890404  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSDS 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT), DINV(N) 
      INTEGER  IA(NELT), JA(NELT) 
C     .. Local Scalars .. 
      INTEGER  ICOL 
C***FIRST EXECUTABLE STATEMENT  DSDS 
C 
C         Assume the Diagonal elements are the first in each column. 
C         This loop should *VECTORIZE*.  If it does not you may have 
C         to add a compiler directive.  We do not check for a zero 
C         (or near zero) diagonal element since this would interfere 
C         with vectorization.  If this makes you nervous put a check 
C         in!  It will run much slower. 
C 
      DO 10 ICOL = 1, N 
         DINV(ICOL) = 1.0D0/A(JA(ICOL)) 
 10   CONTINUE 
C 
      RETURN 
C------------- LAST LINE OF DSDS FOLLOWS ---------------------------- 
      END 
 
 
 
*DECK DSDI 
      SUBROUTINE DSDI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C***BEGIN PROLOGUE  DSDI 
C***PURPOSE  Diagonal Matrix Vector Multiply. 
C            Routine to calculate the product  X = DIAG*B, where DIAG 
C            is a diagonal matrix. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D1B4 
C***TYPE      DOUBLE PRECISION (SSDI-S, DSDI-D) 
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10) 
C     DOUBLE PRECISION B(N), X(N), A(NELT), RWORK(USER DEFINED) 
C 
C     CALL DSDI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C 
C *Arguments: 
C N      :IN       integer  
C         Order of the Matrix. 
C B      :IN       Double Precision B(N). 
C         Vector to multiply the diagonal by. 
C X      :OUT      Double Precision X(N). 
C         Result of DIAG*B. 
C NELT   :DUMMY    integer . 
C IA     :DUMMY    integer  IA(NELT). 
C JA     :DUMMY    integer  JA(NELT). 
C A      :DUMMY    Double Precision A(NELT). 
C ISYM   :DUMMY    integer . 
C         These are for compatibility with SLAP MSOLVE calling sequence. 
C RWORK  :IN       Double Precision RWORK(USER DEFINED). 
C         Work array holding the diagonal of some matrix to scale 
C         B by.  This array must be set by the user or by a call 
C         to the SLAP routine DSDS or DSD2S.  The length of RWORK 
C         must be >= IWORK(4)+N. 
C IWORK  :IN       integer  IWORK(10). 
C         IWORK(4) holds the offset into RWORK for the diagonal matrix 
C         to scale B by.  This is usually set up by the SLAP pre- 
C         conditioner setup routines DSDS or DSD2S. 
C 
C *Description: 
C         This routine is supplied with the SLAP package to perform 
C         the  MSOLVE  operation for iterative drivers that require 
C         diagonal  Scaling  (e.g., DSDCG, DSDBCG).   It  conforms 
C         to the SLAP MSOLVE CALLING CONVENTION  and hence does not 
C         require an interface routine as do some of the other pre- 
C         conditioners supplied with SLAP. 
C 
C***SEE ALSO  DSDS, DSD2S 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   871119  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSDI 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT), B(N), rwork(nelt), X(N) 
      INTEGER  IA(NELT), IWORK(10), JA(NELT) 
C     .. Local Scalars .. 
      INTEGER  I, LOCD 
C***FIRST EXECUTABLE STATEMENT  DSDI 
C 
C         Determine where the inverse of the diagonal 
C         is in the work array and then scale by it. 
C 
      LOCD = IWORK(4) - 1 
      DO 10 I = 1, N 
         X(I) = RWORK(LOCD+I)*B(I) 
 10   CONTINUE 
      RETURN 
C------------- LAST LINE OF DSDI FOLLOWS ---------------------------- 
      END 
 
 
*DECK DSMV 
      SUBROUTINE DSMV (N, X, Y, NELT, IA, JA, A, ISYM) 
C***BEGIN PROLOGUE  DSMV 
C***PURPOSE  SLAP Column Format Sparse Matrix Vector Product. 
C            Routine to calculate the sparse matrix vector product: 
C            Y = A*X. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D1B4 
C***TYPE      DOUBLE PRECISION (SSMV-S, DSMV-D) 
C***KEYWORDS  MATRIX VECTOR MULTIPLY, SLAP, SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM 
C     DOUBLE PRECISION X(N), Y(N), A(NELT) 
C 
C     CALL DSMV(N, X, Y, NELT, IA, JA, A, ISYM ) 
C 
C *Arguments: 
C N      :IN       integer . 
C         Order of the Matrix. 
C X      :IN       Double Precision X(N). 
C         The vector that should be multiplied by the matrix. 
C Y      :OUT      Double Precision Y(N). 
C         The product of the matrix and the vector. 
C NELT   :IN       integer . 
C         Number of Non-Zeros stored in A. 
C IA     :IN       integer  IA(NELT). 
C JA     :IN       integer  JA(NELT). 
C A      :IN       Double Precision A(NELT). 
C         These arrays should hold the matrix A in the SLAP Column 
C         format.  See "Description", below. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all non-zero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the upper 
C         or lower triangle of the matrix is stored. 
C 
C *Description 
C       =================== S L A P Column format ================== 
C       This routine  requires that  the matrix A  be stored in  the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C 
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column): 
C 
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       With  the SLAP  format  the "inner  loops" of  this  routine 
C       should vectorize   on machines with   hardware  support  for 
C       vector gather/scatter operations.  Your compiler may require 
C       a  compiler directive  to  convince   it that there  are  no 
C       implicit vector  dependencies.  Compiler directives  for the 
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied 
C       with the standard SLAP distribution. 
C 
C *Cautions: 
C     This   routine   assumes  that  the matrix A is stored in SLAP 
C     Column format.  It does not check  for  this (for  speed)  and 
C     evil, ugly, ornery and nasty things  will happen if the matrix 
C     data  structure  is,  in fact, not SLAP Column.  Beware of the 
C     wrong data structure!!! 
C 
C***SEE ALSO  DSMTV 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   871119  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSMV 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT), X(N), Y(N) 
      INTEGER  IA(NELT), JA(NELT) 
C     .. Local Scalars .. 
      INTEGER  I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND 
C***FIRST EXECUTABLE STATEMENT  DSMV 
C 
C         Zero out the result vector. 
C 
      DO 10 I = 1, N 
         Y(I) = 0 
 10   CONTINUE 
C 
C         Multiply by A. 
C 
CVD$R NOCONCUR 
      DO 30 ICOL = 1, N 
         IBGN = JA(ICOL) 
         IEND = JA(ICOL+1)-1 
CLLL. OPTION ASSERT (NOHAZARD) 
CDIR$ IVDEP 
CVD$ NODEPCHK 
         DO 20 I = IBGN, IEND 
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL) 
 20      CONTINUE 
 30   CONTINUE 
C 
      IF( ISYM.EQ.1 ) THEN 
C 
C         The matrix is non-symmetric.  Need to get the other half in... 
C         This loops assumes that the diagonal is the first entry in 
C         each column. 
C 
         DO 50 IROW = 1, N 
            JBGN = JA(IROW)+1 
            JEND = JA(IROW+1)-1 
            IF( JBGN.GT.JEND ) GOTO 50 
            DO 40 J = JBGN, JEND 
               Y(IROW) = Y(IROW) + A(J)*X(IA(J)) 
 40         CONTINUE 
 50      CONTINUE 
      ENDIF 
      RETURN 
C------------- LAST LINE OF DSMV FOLLOWS ---------------------------- 
      END 
 
 
*DECK DSILUS 
      SUBROUTINE DSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV, 
     $   NU, IU, JU, U, NROW, NCOL) 
C***BEGIN PROLOGUE  DSILUS 
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up. 
C            Routine to generate the incomplete LDU decomposition of a 
C            matrix.  The unit lower triangular factor L is stored by 
C            rows and the unit upper triangular factor U is stored by 
C            columns.  The inverse of the diagonal matrix D is stored. 
C            No fill in is allowed. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D2E 
C***TYPE      DOUBLE PRECISION (SSILUS-S, DSILUS-D) 
C***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION, 
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM 
C     INTEGER  NL, IL(NL), JL(NL), NU, IU(NU), JU(NU) 
C     INTEGER  NROW(N), NCOL(N) 
C     DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU) 
C 
C     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, 
C    $    DINV, NU, IU, JU, U, NROW, NCOL ) 
C 
C *Arguments: 
C N      :IN       integer  
C         Order of the Matrix. 
C NELT   :IN       integer . 
C         Number of elements in arrays IA, JA, and A. 
C IA     :IN       integer  IA(NELT). 
C JA     :IN       integer  JA(NELT). 
C A      :IN       Double Precision A(NELT). 
C         These arrays should hold the matrix A in the SLAP Column 
C         format.  See "Description", below. 
C ISYM   :IN       integer . 
C         Flag to indicate symmetric storage format. 
C         If ISYM=0, all non-zero entries of the matrix are stored. 
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored. 
C NL     :OUT      integer . 
C         Number of non-zeros in the L array. 
C IL     :OUT      integer  IL(NL). 
C JL     :OUT      integer  JL(NL). 
C L      :OUT      Double Precision L(NL). 
C         IL, JL, L  contain the unit lower triangular factor of  the 
C         incomplete decomposition  of some  matrix stored  in   SLAP 
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See 
C         "DESCRIPTION", below for more details about the SLAP format. 
C NU     :OUT      integer . 
C         Number of non-zeros in the U array. 
C IU     :OUT      integer  IU(NU). 
C JU     :OUT      integer  JU(NU). 
C U      :OUT      Double Precision     U(NU). 
C         IU, JU, U contain   the unit upper triangular factor of the 
C         incomplete  decomposition    of some matrix  stored in SLAP 
C         Column  format.   The Diagonal of ones   *IS*  stored.  See 
C         "Description", below  for  more  details  about  the   SLAP 
C         format. 
C NROW   :WORK     integer  NROW(N). 
C         NROW(I) is the number of non-zero elements in the I-th row 
C         of L. 
C NCOL   :WORK     integer  NCOL(N). 
C         NCOL(I) is the number of non-zero elements in the I-th 
C         column of U. 
C 
C *Description 
C       IL, JL, L should contain the unit  lower triangular factor of 
C       the incomplete decomposition of the A matrix  stored in SLAP 
C       Row format.  IU, JU, U should contain  the unit upper factor 
C       of the  incomplete decomposition of  the A matrix  stored in 
C       SLAP Column format This ILU factorization can be computed by 
C       the DSILUS routine. The diagonals (which are all one's) are 
C       stored. 
C 
C       =================== S L A P Column format ================== 
C 
C       This routine  requires that  the matrix A  be stored in  the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C 
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column): 
C 
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       ==================== S L A P Row format ==================== 
C 
C       This routine requires  that the matrix A  be  stored  in the 
C       SLAP  Row format.   In this format  the non-zeros are stored 
C       counting across  rows (except for the diagonal  entry, which 
C       must  appear first  in each  "row")  and  are stored  in the 
C       double precision  array A.  In other words, for each row  in 
C       the matrix  put the diagonal  entry in A.   Then put in  the 
C       other  non-zero elements  going across  the row  (except the 
C       diagonal) in order.  The JA array holds the column index for 
C       each non-zero.  The IA array holds the offsets  into the JA, 
C       A  arrays  for  the   beginning  of  each  row.    That  is, 
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW- 
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1) 
C       are  the last elements  of the  IROW-th row.   Note  that we 
C       always have  IA(N+1) = NELT+1, where N is the number of rows 
C       in the matrix  and  NELT is the  number of non-zeros  in the 
C       matrix. 
C 
C       Here is an example of the SLAP Row storage format for a  5x5 
C       Matrix (in the A and JA arrays '|' denotes the end of a row): 
C 
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53 
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C***SEE ALSO  SILUR 
C***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations, 
C                  Johns Hopkins University Press, Baltimore, Maryland, 
C                  1983. 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   890404  DATE WRITTEN 
C   890404  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   920929  Corrected format of reference.  (FNF) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSILUS 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT, NL, NU 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT), DINV(N), L(NL), U(NU) 
      INTEGER  IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU), 
     $        NCOL(N), NROW(N) 
C     .. Local Scalars .. 
      DOUBLE PRECISION TEMP 
      INTEGER  I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2, 
     $        INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC, 
     $        KR 
C***FIRST EXECUTABLE STATEMENT  DSILUS 
C 
C         Count number of elements in each row of the lower triangle. 
C 
      DO 10 I=1,N 
         NROW(I) = 0 
         NCOL(I) = 0 
 10   CONTINUE 
CVD$R NOCONCUR 
CVD$R NOVECTOR 
      DO 30 ICOL = 1, N 
         JBGN = JA(ICOL)+1 
         JEND = JA(ICOL+1)-1 
         IF( JBGN.LE.JEND ) THEN 
            DO 20 J = JBGN, JEND 
               IF( IA(J).LT.ICOL ) THEN 
                  NCOL(ICOL) = NCOL(ICOL) + 1 
               ELSE 
                  NROW(IA(J)) = NROW(IA(J)) + 1 
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1 
               ENDIF 
 20         CONTINUE 
         ENDIF 
 30   CONTINUE 
      JU(1) = 1 
      IL(1) = 1 
      DO 40 ICOL = 1, N 
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL) 
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL) 
         NROW(ICOL) = IL(ICOL) 
         NCOL(ICOL) = JU(ICOL) 
 40   CONTINUE 
C 
C         Copy the matrix A into the L and U structures. 
      DO 60 ICOL = 1, N 
         DINV(ICOL) = A(JA(ICOL)) 
         JBGN = JA(ICOL)+1 
         JEND = JA(ICOL+1)-1 
         IF( JBGN.LE.JEND ) THEN 
            DO 50 J = JBGN, JEND 
               IROW = IA(J) 
               IF( IROW.LT.ICOL ) THEN 
C         Part of the upper triangle. 
                  IU(NCOL(ICOL)) = IROW 
                  U(NCOL(ICOL)) = A(J) 
                  NCOL(ICOL) = NCOL(ICOL) + 1 
               ELSE 
C         Part of the lower triangle (stored by row). 
                  JL(NROW(IROW)) = ICOL 
                  L(NROW(IROW)) = A(J) 
                  NROW(IROW) = NROW(IROW) + 1 
                  IF( ISYM.NE.0 ) THEN 
C         Symmetric...Copy lower triangle into upper triangle as well. 
                     IU(NCOL(IROW)) = ICOL 
                     U(NCOL(IROW)) = A(J) 
                     NCOL(IROW) = NCOL(IROW) + 1 
                  ENDIF 
               ENDIF 
 50         CONTINUE 
         ENDIF 
 60   CONTINUE 
C 
C         Sort the rows of L and the columns of U. 
      DO 110 K = 2, N 
         JBGN = JU(K) 
         JEND = JU(K+1)-1 
         IF( JBGN.LT.JEND ) THEN 
            DO 80 J = JBGN, JEND-1 
               DO 70 I = J+1, JEND 
                  IF( IU(J).GT.IU(I) ) THEN 
                     ITEMP = IU(J) 
                     IU(J) = IU(I) 
                     IU(I) = ITEMP 
                     TEMP = U(J) 
                     U(J) = U(I) 
                     U(I) = TEMP 
                  ENDIF 
 70            CONTINUE 
 80         CONTINUE 
         ENDIF 
         IBGN = IL(K) 
         IEND = IL(K+1)-1 
         IF( IBGN.LT.IEND ) THEN 
            DO 100 I = IBGN, IEND-1 
               DO 90 J = I+1, IEND 
                  IF( JL(I).GT.JL(J) ) THEN 
                     JTEMP = JU(I) 
                     JU(I) = JU(J) 
                     JU(J) = JTEMP 
                     TEMP = L(I) 
                     L(I) = L(J) 
                     L(J) = TEMP 
                  ENDIF 
 90            CONTINUE 
 100        CONTINUE 
         ENDIF 
 110  CONTINUE 
C 
C         Perform the incomplete LDU decomposition. 
      DO 300 I=2,N 
C 
C           I-th row of L 
         INDX1 = IL(I) 
         INDX2 = IL(I+1) - 1 
         IF(INDX1 .GT. INDX2) GO TO 200 
         DO 190 INDX=INDX1,INDX2 
            IF(INDX .EQ. INDX1) GO TO 180 
            INDXR1 = INDX1 
            INDXR2 = INDX - 1 
            INDXC1 = JU(JL(INDX)) 
            INDXC2 = JU(JL(INDX)+1) - 1 
            IF(INDXC1 .GT. INDXC2) GO TO 180 
 160        KR = JL(INDXR1) 
 170        KC = IU(INDXC1) 
            IF(KR .GT. KC) THEN 
               INDXC1 = INDXC1 + 1 
               IF(INDXC1 .LE. INDXC2) GO TO 170 
            ELSEIF(KR .LT. KC) THEN 
               INDXR1 = INDXR1 + 1 
               IF(INDXR1 .LE. INDXR2) GO TO 160 
            ELSEIF(KR .EQ. KC) THEN 
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1) 
               INDXR1 = INDXR1 + 1 
               INDXC1 = INDXC1 + 1 
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160 
            ENDIF 
 180        L(INDX) = L(INDX)/DINV(JL(INDX)) 
 190     CONTINUE 
C 
C         I-th column of U 
 200     INDX1 = JU(I) 
         INDX2 = JU(I+1) - 1 
         IF(INDX1 .GT. INDX2) GO TO 260 
         DO 250 INDX=INDX1,INDX2 
            IF(INDX .EQ. INDX1) GO TO 240 
            INDXC1 = INDX1 
            INDXC2 = INDX - 1 
            INDXR1 = IL(IU(INDX)) 
            INDXR2 = IL(IU(INDX)+1) - 1 
            IF(INDXR1 .GT. INDXR2) GO TO 240 
 210        KR = JL(INDXR1) 
 220        KC = IU(INDXC1) 
            IF(KR .GT. KC) THEN 
               INDXC1 = INDXC1 + 1 
               IF(INDXC1 .LE. INDXC2) GO TO 220 
            ELSEIF(KR .LT. KC) THEN 
               INDXR1 = INDXR1 + 1 
               IF(INDXR1 .LE. INDXR2) GO TO 210 
            ELSEIF(KR .EQ. KC) THEN 
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1) 
               INDXR1 = INDXR1 + 1 
               INDXC1 = INDXC1 + 1 
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210 
            ENDIF 
 240        U(INDX) = U(INDX)/DINV(IU(INDX)) 
 250     CONTINUE 
C 
C         I-th diagonal element 
 260     INDXR1 = IL(I) 
         INDXR2 = IL(I+1) - 1 
         IF(INDXR1 .GT. INDXR2) GO TO 300 
         INDXC1 = JU(I) 
         INDXC2 = JU(I+1) - 1 
         IF(INDXC1 .GT. INDXC2) GO TO 300 
 270     KR = JL(INDXR1) 
 280     KC = IU(INDXC1) 
         IF(KR .GT. KC) THEN 
            INDXC1 = INDXC1 + 1 
            IF(INDXC1 .LE. INDXC2) GO TO 280 
         ELSEIF(KR .LT. KC) THEN 
            INDXR1 = INDXR1 + 1 
            IF(INDXR1 .LE. INDXR2) GO TO 270 
         ELSEIF(KR .EQ. KC) THEN 
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1) 
            INDXR1 = INDXR1 + 1 
            INDXC1 = INDXC1 + 1 
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270 
         ENDIF 
C 
 300  CONTINUE 
C 
C         Replace diagonal elements by their inverses. 
CVD$ VECTOR 
      DO 430 I=1,N 
         DINV(I) = 1.0D0/DINV(I) 
 430  CONTINUE 
C 
      RETURN 
C------------- LAST LINE OF DSILUS FOLLOWS ---------------------------- 
      END 
 
 
*DECK DSLUI 
      SUBROUTINE DSLUI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK) 
C***BEGIN PROLOGUE  DSLUI 
C***PURPOSE  SLAP MSOLVE for LDU Factorization. 
C            This routine acts as an interface between the SLAP generic 
C            MSOLVE calling convention and the routine that actually 
C                           -1 
C            computes  (LDU)  B = X. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D2E 
C***TYPE      DOUBLE PRECISION (SSLUI-S, DSLUI-D) 
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE, 
C             SLAP, SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C       It is assumed that RWORK and IWORK have initialized with 
C       the information required for DSLUI2: 
C          IWORK(1) = Starting location of IL in IWORK. 
C          IWORK(2) = Starting location of JL in IWORK. 
C          IWORK(3) = Starting location of IU in IWORK. 
C          IWORK(4) = Starting location of JU in IWORK. 
C          IWORK(5) = Starting location of L in RWORK. 
C          IWORK(6) = Starting location of DINV in RWORK. 
C          IWORK(7) = Starting location of U in RWORK. 
C       See the DESCRIPTION of DSLUI2 for details. 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  DSLUI2 
C***REVISION HISTORY  (YYMMDD) 
C   871119  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   921113  Corrected C***CATEGORY line.  (FNF) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSLUI 
C     .. Scalar Arguments .. 
      INTEGER  ISYM, N, NELT 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(NELT), B(N), rwork(nelt), X(N) 
      INTEGER  IA(NELT), IWORK(10), JA(NELT) 
C     .. Local Scalars .. 
      INTEGER  LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU 
C     .. External Subroutines .. 
      EXTERNAL DSLUI2 
C***FIRST EXECUTABLE STATEMENT  DSLUI 
C 
C         Pull out the locations of the arrays holding the ILU 
C         factorization. 
C 
      LOCIL = IWORK(1) 
      LOCJL = IWORK(2) 
      LOCIU = IWORK(3) 
      LOCJU = IWORK(4) 
      LOCL = IWORK(5) 
      LOCDIN = IWORK(6) 
      LOCU = IWORK(7) 
C 
C         Solve the system LUx = b 
      CALL DSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL), 
     $     RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) ) 
C 
      RETURN 
C------------- LAST LINE OF DSLUI FOLLOWS ---------------------------- 
      END 
 
 
 
 
*DECK XERMSG 
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL) 
C***BEGIN PROLOGUE  XERMSG 
C***PURPOSE  Process error messages for SLATEC and other libraries. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3C 
C***TYPE      ALL (XERMSG-A) 
C***KEYWORDS  ERROR MESSAGE, XERROR 
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL) 
C***DESCRIPTION 
C 
C   XERMSG processes a diagnostic message in a manner determined by the 
C   value of LEVEL and the current value of the library error control 
C   flag, KONTRL.  See subroutine XSETF for details. 
C 
C    LIBRAR   A character constant (or character variable) with the name 
C             of the library.  This will be 'SLATEC' for the SLATEC 
C             Common Math Library.  The error handling package is 
C             general enough to be used by many libraries 
C             simultaneously, so it is desirable for the routine that 
C             detects and reports an error to identify the library name 
C             as well as the routine name. 
C 
C    SUBROU   A character constant (or character variable) with the name 
C             of the routine that detected the error.  Usually it is the 
C             name of the routine that is calling XERMSG.  There are 
C             some instances where a user callable library routine calls 
C             lower level subsidiary routines where the error is 
C             detected.  In such cases it may be more informative to 
C             supply the name of the routine the user called rather than 
C             the name of the subsidiary routine that detected the 
C             error. 
C 
C    MESSG    A character constant (or character variable) with the text 
C             of the error or warning message.  In the example below, 
C             the message is a character constant that contains a 
C             generic message. 
C 
C                   CALL XERMSG ('SLATEC', 'MMPY', 
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION', 
C                  *3, 1) 
C 
C             It is possible (and is sometimes desirable) to generate a 
C             specific message--e.g., one that contains actual numeric 
C             values.  Specific numeric values can be converted into 
C             character strings using formatted WRITE statements into 
C             character variables.  This is called standard Fortran 
C             internal file I/O and is exemplified in the first three 
C             lines of the following example.  You can also catenate 
C             substrings of characters to construct the error message. 
C             Here is an example showing the use of both writing to 
C             an internal file and catenating character strings. 
C 
C                   CHARACTER*5 CHARN, CHARL 
C                   WRITE (CHARN,10) N 
C                   WRITE (CHARL,10) LDA 
C                10 FORMAT(I5) 
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN// 
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'// 
C                  *   CHARL, 3, 1) 
C 
C             There are two subtleties worth mentioning.  One is that 
C             the // for character catenation is used to construct the 
C             error message so that no single character constant is 
C             continued to the next line.  This avoids confusion as to 
C             whether there are trailing blanks at the end of the line. 
C             The second is that by catenating the parts of the message 
C             as an actual argument rather than encoding the entire 
C             message into one large character variable, we avoid 
C             having to know how long the message will be in order to 
C             declare an adequate length for that large character 
C             variable.  XERMSG calls XERPRN to print the message using 
C             multiple lines if necessary.  If the message is very long, 
C             XERPRN will break it into pieces of 72 characters (as 
C             requested by XERMSG) for printing on multiple lines. 
C             Also, XERMSG asks XERPRN to prefix each line with ' *  ' 
C             so that the total line length could be 76 characters. 
C             Note also that XERPRN scans the error message backwards 
C             to ignore trailing blanks.  Another feature is that 
C             the substring '$$' is treated as a new line sentinel 
C             by XERPRN.  If you want to construct a multiline 
C             message without having to count out multiples of 72 
C             characters, just use '$$' as a separator.  '$$' 
C             obviously must occur within 72 characters of the 
C             start of each line to have its intended effect since 
C             XERPRN is asked to wrap around at 72 characters in 
C             addition to looking for '$$'. 
C 
C    NERR     An integer  value that is chosen by the library routine's 
C             author.  It must be in the range -99 to 999 (three 
C             printable digits).  Each distinct error should have its 
C             own error number.  These error numbers should be described 
C             in the machine readable documentation for the routine. 
C             The error numbers need be unique only within each routine, 
C             so it is reasonable for each routine to start enumerating 
C             errors from 1 and proceeding to the next integer . 
C 
C    LEVEL    An integer  value in the range 0 to 2 that indicates the 
C             level (severity) of the error.  Their meanings are 
C 
C            -1  A warning message.  This is used if it is not clear 
C                that there really is an error, but the user's attention 
C                may be needed.  An attempt is made to only print this 
C                message once. 
C 
C             0  A warning message.  This is used if it is not clear 
C                that there really is an error, but the user's attention 
C                may be needed. 
C 
C             1  A recoverable error.  This is used even if the error is 
C                so serious that the routine cannot return any useful 
C                answer.  If the user has told the error package to 
C                return after recoverable errors, then XERMSG will 
C                return to the Library routine which can then return to 
C                the user's routine.  The user may also permit the error 
C                package to terminate the program upon encountering a 
C                recoverable error. 
C 
C             2  A fatal error.  XERMSG will not return to its caller 
C                after it receives a fatal error.  This level should 
C                hardly ever be used; it is much better to allow the 
C                user a chance to recover.  An example of one of the few 
C                cases in which it is permissible to declare a level 2 
C                error is a reverse communication Library routine that 
C                is likely to be called repeatedly until it integrates 
C                across some interval.  If there is a serious error in 
C                the input such that another step cannot be taken and 
C                the Library routine is called again without the input 
C                error having been corrected by the caller, the Library 
C                routine will probably be called forever with improper 
C                input.  In this case, it is reasonable to declare the 
C                error to be fatal. 
C 
C    Each of the arguments to XERMSG is input; none will be modified by 
C    XERMSG.  A routine may make multiple calls to XERMSG with warning 
C    level messages; however, after a call to XERMSG with a recoverable 
C    error, the routine should return to the user.  Do not try to call 
C    XERMSG with a second recoverable error after the first recoverable 
C    error because the error package saves the error number.  The user 
C    can retrieve this error number by calling another entry point in 
C    the error handling package and then clear the error number when 
C    recovering from the error.  Calling XERMSG in succession causes the 
C    old error number to be overwritten by the latest error number. 
C    This is considered harmless for error numbers associated with 
C    warning messages but must not be done for error numbers of serious 
C    errors.  After a call to XERMSG with a recoverable error, the user 
C    must be given a chance to call NUMXER or XERCLR to retrieve or 
C    clear the error number. 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE 
C***REVISION HISTORY  (YYMMDD) 
C   880101  DATE WRITTEN 
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988. 
C           THERE ARE TWO BASIC CHANGES. 
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO 
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES 
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS 
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE 
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER 
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY 
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE 
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76. 
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE 
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE 
C               OF LOWER CASE. 
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30. 
C           THE PRINCIPAL CHANGES ARE 
C           1.  CLARIFY COMMENTS IN THE PROLOGUES 
C           2.  RENAME XRPRNT TO XERPRN 
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES 
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE / 
C               CHARACTER FOR NEW RECORDS. 
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO 
C           CLEAN UP THE CODING. 
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN 
C           PREFIX. 
C   891013  REVISED TO CORRECT COMMENTS. 
C   891214  Prologue converted to Version 4.0 format.  (WRB) 
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but 
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added 
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and 
C           XERCTL to XERCNT.  (RWC) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XERMSG 
      CHARACTER*(*) LIBRAR, SUBROU, MESSG 
      CHARACTER*8 XLIBR, XSUBR 
      CHARACTER*72  TEMP 
      CHARACTER*20  LFIRST 
C***FIRST EXECUTABLE STATEMENT  XERMSG 
      LKNTRL = J4SAVE (2, 0, .FALSE.) 
      MAXMES = J4SAVE (4, 0, .FALSE.) 
C 
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL. 
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE 
C          SHOULD BE PRINTED. 
C 
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN 
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE, 
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2. 
C 
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR. 
     $   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN 
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // 
     $      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// 
     $      'JOB ABORT DUE TO FATAL ERROR.', 72) 
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY) 
         CALL XERHLT (' ***XERMSG -- INVALID INPUT') 
         RETURN 
      ENDIF 
C 
C       RECORD THE MESSAGE. 
C 
      I = J4SAVE (1, NERR, .TRUE.) 
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT) 
C 
C       HANDLE PRINT-ONCE WARNING MESSAGES. 
C 
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN 
C 
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG. 
C 
      XLIBR  = LIBRAR 
      XSUBR  = SUBROU 
      LFIRST = MESSG 
      LERR   = NERR 
      LLEVEL = LEVEL 
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL) 
C 
      LKNTRL = MAX(-2, MIN(2,LKNTRL)) 
      MKNTRL = ABS(LKNTRL) 
C 
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS 
C       ZERO AND THE ERROR IS NOT FATAL. 
C 
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30 
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30 
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30 
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30 
C 
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A 
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS) 
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG 
C       IS NOT ZERO. 
C 
      IF (LKNTRL .NE. 0) THEN 
         TEMP(1:21) = 'MESSAGE FROM ROUTINE ' 
         I = MIN(LEN(SUBROU), 16) 
         TEMP(22:21+I) = SUBROU(1:I) 
         TEMP(22+I:33+I) = ' IN LIBRARY ' 
         LTEMP = 33 + I 
         I = MIN(LEN(LIBRAR), 16) 
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I) 
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.' 
         LTEMP = LTEMP + I + 1 
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72) 
      ENDIF 
C 
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE 
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE 
C       FROM EACH OF THE FOLLOWING THREE OPTIONS. 
C       1.  LEVEL OF THE MESSAGE 
C              'INFORMATIVE MESSAGE' 
C              'POTENTIALLY RECOVERABLE ERROR' 
C              'FATAL ERROR' 
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE 
C              'PROG CONTINUES' 
C              'PROG ABORTED' 
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK 
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS 
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.) 
C              'TRACEBACK REQUESTED' 
C              'TRACEBACK NOT REQUESTED' 
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT 
C       EXCEED 74 CHARACTERS. 
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED. 
C 
      IF (LKNTRL .GT. 0) THEN 
C 
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL. 
C 
         IF (LEVEL .LE. 0) THEN 
            TEMP(1:20) = 'INFORMATIVE MESSAGE,' 
            LTEMP = 20 
         ELSEIF (LEVEL .EQ. 1) THEN 
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,' 
            LTEMP = 30 
         ELSE 
            TEMP(1:12) = 'FATAL ERROR,' 
            LTEMP = 12 
         ENDIF 
C 
C       THEN WHETHER THE PROGRAM WILL CONTINUE. 
C 
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR. 
     $       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN 
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,' 
            LTEMP = LTEMP + 14 
         ELSE 
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,' 
            LTEMP = LTEMP + 16 
         ENDIF 
C 
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK. 
C 
         IF (LKNTRL .GT. 0) THEN 
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED' 
            LTEMP = LTEMP + 20 
         ELSE 
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED' 
            LTEMP = LTEMP + 24 
         ENDIF 
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72) 
      ENDIF 
C 
C       NOW SEND OUT THE MESSAGE. 
C 
      CALL XERPRN (' *  ', -1, MESSG, 72) 
C 
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A 
C          TRACEBACK. 
C 
      IF (LKNTRL .GT. 0) THEN 
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR 
         DO 10 I=16,22 
            IF (TEMP(I:I) .NE. ' ') GO TO 20 
   10    CONTINUE 
C 
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72) 
         CALL FDUMP 
      ENDIF 
C 
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE. 
C 
      IF (LKNTRL .NE. 0) THEN 
         CALL XERPRN (' *  ', -1, ' ', 72) 
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72) 
         CALL XERPRN ('    ',  0, ' ', 72) 
      ENDIF 
C 
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE 
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN. 
C 
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN 
C 
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A 
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR 
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT. 
C 
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN 
         IF (LEVEL .EQ. 1) THEN 
            CALL XERPRN 
     $         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72) 
         ELSE 
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72) 
         ENDIF 
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY) 
         CALL XERHLT (' ') 
      ELSE 
         CALL XERHLT (MESSG) 
      ENDIF 
      RETURN 
      END 
 
 
 
*DECK J4SAVE 
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET) 
C***BEGIN PROLOGUE  J4SAVE 
C***SUBSIDIARY 
C***PURPOSE  Save or recall global variables needed by error 
C            handling routines. 
C***LIBRARY   SLATEC (XERROR) 
C***TYPE      INTEGER  (J4SAVE-I) 
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C     Abstract 
C        J4SAVE saves and recalls several global variables needed 
C        by the library error handling routines. 
C 
C     Description of Parameters 
C      --Input-- 
C        IWHICH - Index of item desired. 
C                = 1 Refers to current error number. 
C                = 2 Refers to current error control flag. 
C                = 3 Refers to current unit number to which error 
C                    messages are to be sent.  (0 means use standard.) 
C                = 4 Refers to the maximum number of times any 
C                     message is to be printed (as set by XERMAX). 
C                = 5 Refers to the total number of units to which 
C                     each error message is to be written. 
C                = 6 Refers to the 2nd unit for error messages 
C                = 7 Refers to the 3rd unit for error messages 
C                = 8 Refers to the 4th unit for error messages 
C                = 9 Refers to the 5th unit for error messages 
C        IVALUE - The value to be set for the IWHICH-th parameter, 
C                 if ISET is .TRUE. . 
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE 
C                 given the value, IVALUE.  If ISET=.FALSE., the 
C                 IWHICH-th parameter will be unchanged, and IVALUE 
C                 is a dummy parameter. 
C      --Output-- 
C        The (old) value of the IWHICH-th parameter will be returned 
C        in the function value, J4SAVE. 
C 
C***SEE ALSO  XERMSG 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   790801  DATE WRITTEN 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   900205  Minor modifications to prologue.  (WRB) 
C   900402  Added TYPE section.  (WRB) 
C   910411  Added KEYWORDS section.  (WRB) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  J4SAVE 
      LOGICAL ISET 
      INTEGER  IPARAM(9) 
      SAVE IPARAM 
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/ 
      DATA IPARAM(5)/1/ 
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/ 
C***FIRST EXECUTABLE STATEMENT  J4SAVE 
      J4SAVE = IPARAM(IWHICH) 
      IF (ISET) IPARAM(IWHICH) = IVALUE 
      RETURN 
      END 
 
 
*DECK XERPRN 
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP) 
C***BEGIN PROLOGUE  XERPRN 
C***SUBSIDIARY 
C***PURPOSE  Print error messages processed by XERMSG. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3C 
C***TYPE      ALL (XERPRN-A) 
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR 
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL) 
C***DESCRIPTION 
C 
C This routine sends one or more lines to each of the (up to five) 
C logical units to which error messages are to be sent.  This routine 
C is called several times by XERMSG, sometimes with a single line to 
C print and sometimes with a (potentially very long) message that may 
C wrap around into multiple lines. 
C 
C PREFIX  Input argument of type CHARACTER.  This argument contains 
C         characters to be put at the beginning of each line before 
C         the body of the message.  No more than 16 characters of 
C         PREFIX will be used. 
C 
C NPREF   Input argument of type INTEGER .  This argument is the number 
C         of characters to use from PREFIX.  If it is negative, the 
C         intrinsic function LEN is used to determine its length.  If 
C         it is zero, PREFIX is not used.  If it exceeds 16 or if 
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be 
C         used.  If NPREF is positive and the length of PREFIX is less 
C         than NPREF, a copy of PREFIX extended with blanks to length 
C         NPREF will be used. 
C 
C MESSG   Input argument of type CHARACTER.  This is the text of a 
C         message to be printed.  If it is a long message, it will be 
C         broken into pieces for printing on multiple lines.  Each line 
C         will start with the appropriate prefix and be followed by a 
C         piece of the message.  NWRAP is the number of characters per 
C         piece; that is, after each NWRAP characters, we break and 
C         start a new line.  In addition the characters '$$' embedded 
C         in MESSG are a sentinel for a new line.  The counting of 
C         characters up to NWRAP starts over for each new line.  The 
C         value of NWRAP typically used by XERMSG is 72 since many 
C         older error messages in the SLATEC Library are laid out to 
C         rely on wrap-around every 72 characters. 
C 
C NWRAP   Input argument of type INTEGER .  This gives the maximum size 
C         piece into which to break MESSG for printing on multiple 
C         lines.  An embedded '$$' ends a line, and the count restarts 
C         at the following character.  If a line break does not occur 
C         on a blank (it would split a word) that word is moved to the 
C         next line.  Values of NWRAP less than 16 will be treated as 
C         16.  Values of NWRAP greater than 132 will be treated as 132. 
C         The actual line length will be NPREF + NWRAP after NPREF has 
C         been adjusted to fall between 0 and 16 and NWRAP has been 
C         adjusted to fall between 16 and 132. 
C 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  I1MACH, XGETUA 
C***REVISION HISTORY  (YYMMDD) 
C   880621  DATE WRITTEN 
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF 
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK 
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE 
C           SLASH CHARACTER IN FORMAT STATEMENTS. 
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO 
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK 
C           LINES TO BE PRINTED. 
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF 
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH. 
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH. 
C   891214  Prologue converted to Version 4.0 format.  (WRB) 
C   900510  Added code to break messages between words.  (RWC) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XERPRN 
      CHARACTER*(*) PREFIX, MESSG 
      INTEGER  NPREF, NWRAP 
      CHARACTER*148 CBUFF 
      INTEGER  IU(5), NUNIT 
      CHARACTER*2 NEWLIN 
      PARAMETER (NEWLIN = '$$') 
C***FIRST EXECUTABLE STATEMENT  XERPRN 
      CALL XGETUA(IU,NUNIT) 
C 
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD 
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD 
C       ERROR MESSAGE UNIT. 
C 
      N = I1MACH(4) 
      DO 10 I=1,NUNIT 
         IF (IU(I) .EQ. 0) IU(I) = N 
   10 CONTINUE 
C 
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE 
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING 
C       THE REST OF THIS ROUTINE. 
C 
      IF ( NPREF .LT. 0 ) THEN 
         LPREF = LEN(PREFIX) 
      ELSE 
         LPREF = NPREF 
      ENDIF 
      LPREF = MIN(16, LPREF) 
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX 
C 
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE 
C       TIME FROM MESSG TO PRINT ON ONE LINE. 
C 
      LWRAP = MAX(16, MIN(132, NWRAP)) 
C 
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS. 
C 
      LENMSG = LEN(MESSG) 
      N = LENMSG 
      DO 20 I=1,N 
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30 
         LENMSG = LENMSG - 1 
   20 CONTINUE 
   30 CONTINUE 
C 
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE. 
C 
      IF (LENMSG .EQ. 0) THEN 
         CBUFF(LPREF+1:LPREF+1) = ' ' 
         DO 40 I=1,NUNIT 
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1) 
   40    CONTINUE 
         RETURN 
      ENDIF 
C 
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING 
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL. 
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT. 
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED. 
C 
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE 
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE 
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH 
C       OF THE SECOND ARGUMENT. 
C 
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE 
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER 
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT 
C       POSITION NEXTC. 
C 
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE 
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE 
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC, 
C                       WHICHEVER IS LESS. 
C 
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC: 
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE 
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY 
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION 
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF 
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE 
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC 
C                       SHOULD BE INCREMENTED BY 2. 
C 
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP. 
C 
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1 
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS 
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ. 
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY 
C                       AT THE END OF A LINE. 
C 
      NEXTC = 1 
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN) 
      IF (LPIECE .EQ. 0) THEN 
C 
C       THERE WAS NO NEW LINE SENTINEL FOUND. 
C 
         IDELTA = 0 
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC) 
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN 
            DO 52 I=LPIECE+1,2,-1 
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN 
                  LPIECE = I-1 
                  IDELTA = 1 
                  GOTO 54 
               ENDIF 
   52       CONTINUE 
         ENDIF 
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) 
         NEXTC = NEXTC + LPIECE + IDELTA 
      ELSEIF (LPIECE .EQ. 1) THEN 
C 
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1). 
C       DON'T PRINT A BLANK LINE. 
C 
         NEXTC = NEXTC + 2 
         GO TO 50 
      ELSEIF (LPIECE .GT. LWRAP+1) THEN 
C 
C       LPIECE SHOULD BE SET DOWN TO LWRAP. 
C 
         IDELTA = 0 
         LPIECE = LWRAP 
         DO 56 I=LPIECE+1,2,-1 
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN 
               LPIECE = I-1 
               IDELTA = 1 
               GOTO 58 
            ENDIF 
   56    CONTINUE 
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) 
         NEXTC = NEXTC + LPIECE + IDELTA 
      ELSE 
C 
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1. 
C       WE SHOULD DECREMENT LPIECE BY ONE. 
C 
         LPIECE = LPIECE - 1 
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) 
         NEXTC  = NEXTC + LPIECE + 2 
      ENDIF 
C 
C       PRINT 
C 
      DO 60 I=1,NUNIT 
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE) 
   60 CONTINUE 
C 
      IF (NEXTC .LE. LENMSG) GO TO 50 
      RETURN 
      END 
 
 
 
 
 
*DECK XERSVE 
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, 
     $   ICOUNT) 
C***BEGIN PROLOGUE  XERSVE 
C***SUBSIDIARY 
C***PURPOSE  Record that an error has occurred. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3 
C***TYPE      ALL (XERSVE-A) 
C***KEYWORDS  ERROR, XERROR 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C *Usage: 
C 
C        INTEGER   KFLAG, NERR, LEVEL, ICOUNT 
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG 
C 
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT) 
C 
C *Arguments: 
C 
C        LIBRAR :IN    is the library that the message is from. 
C        SUBROU :IN    is the subroutine that the message is from. 
C        MESSG  :IN    is the message to be saved. 
C        KFLAG  :IN    indicates the action to be performed. 
C                      when KFLAG > 0, the message in MESSG is saved. 
C                      when KFLAG=0 the tables will be dumped and 
C                      cleared. 
C                      when KFLAG < 0, the tables will be dumped and 
C                      not cleared. 
C        NERR   :IN    is the error number. 
C        LEVEL  :IN    is the error severity. 
C        ICOUNT :OUT   the number of times this message has been seen, 
C                      or zero if the table has overflowed and does not 
C                      contain this message specifically.  When KFLAG=0, 
C                      ICOUNT will not be altered. 
C 
C *Description: 
C 
C   Record that this error occurred and possibly dump and clear the 
C   tables. 
C 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  I1MACH, XGETUA 
C***REVISION HISTORY  (YYMMDD) 
C   800319  DATE WRITTEN 
C   861211  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   900413  Routine modified to remove reference to KFLAG.  (WRB) 
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling 
C           sequence, use IF-THEN-ELSE, make number of saved entries 
C           easily changeable, changed routine name from XERSAV to 
C           XERSVE.  (RWC) 
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XERSVE 
      PARAMETER (LENTAB=10) 
      INTEGER  LUN(5) 
      CHARACTER*(*) LIBRAR, SUBROU, MESSG 
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB 
      CHARACTER*20 MESTAB(LENTAB), MES 
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB) 
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG 
      DATA KOUNTX/0/, NMSG/0/ 
C***FIRST EXECUTABLE STATEMENT  XERSVE 
C 
      IF (KFLAG.LE.0) THEN 
C 
C        Dump the table. 
C 
         IF (NMSG.EQ.0) RETURN 
C 
C        Print to each unit. 
C 
         CALL XGETUA (LUN, NUNIT) 
         DO 20 KUNIT = 1,NUNIT 
            IUNIT = LUN(KUNIT) 
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4) 
C 
C           Print the table header. 
C 
            WRITE (IUNIT,9000) 
C 
C           Print body of table. 
C 
            DO 10 I = 1,NMSG 
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), 
     $            NERTAB(I),LEVTAB(I),KOUNT(I) 
   10       CONTINUE 
C 
C           Print number of other errors. 
C 
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX 
            WRITE (IUNIT,9030) 
   20    CONTINUE 
C 
C        Clear the error tables. 
C 
         IF (KFLAG.EQ.0) THEN 
            NMSG = 0 
            KOUNTX = 0 
         ENDIF 
      ELSE 
C 
C        PROCESS A MESSAGE... 
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG, 
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL. 
C 
         LIB = LIBRAR 
         SUB = SUBROU 
         MES = MESSG 
         DO 30 I = 1,NMSG 
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND. 
     $         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND. 
     $         LEVEL.EQ.LEVTAB(I)) THEN 
                  KOUNT(I) = KOUNT(I) + 1 
                  ICOUNT = KOUNT(I) 
                  RETURN 
            ENDIF 
   30    CONTINUE 
C 
         IF (NMSG.LT.LENTAB) THEN 
C 
C           Empty slot found for new message. 
C 
            NMSG = NMSG + 1 
            LIBTAB(I) = LIB 
            SUBTAB(I) = SUB 
            MESTAB(I) = MES 
            NERTAB(I) = NERR 
            LEVTAB(I) = LEVEL 
            KOUNT (I) = 1 
            ICOUNT    = 1 
         ELSE 
C 
C           Table is full. 
C 
            KOUNTX = KOUNTX+1 
            ICOUNT = 0 
         ENDIF 
      ENDIF 
      RETURN 
C 
C     Formats. 
C 
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / 
     $   ' LIBRARY    SUBROUTINE MESSAGE START             NERR', 
     $   '     LEVEL     COUNT') 
 9010 FORMAT (1X,A,3X,A,3X,A,3I10) 
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10) 
 9030 FORMAT (1X) 
      END 
 
 
 
 
*DECK XERHLT 
      SUBROUTINE XERHLT (MESSG) 
C***BEGIN PROLOGUE  XERHLT 
C***SUBSIDIARY 
C***PURPOSE  Abort program execution and print error message. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3C 
C***TYPE      ALL (XERHLT-A) 
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C     Abstract 
C        ***Note*** machine dependent routine 
C        XERHLT aborts the execution of the program. 
C        The error message causing the abort is given in the calling 
C        sequence, in case one needs it for printing on a dayfile, 
C        for example. 
C 
C     Description of Parameters 
C        MESSG is as in XERMSG. 
C 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   790801  DATE WRITTEN 
C   861211  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   900206  Routine changed from user-callable to subsidiary.  (WRB) 
C   900510  Changed calling sequence to delete length of character 
C           and changed routine name from XERABT to XERHLT.  (RWC) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XERHLT 
      CHARACTER*(*) MESSG 
C***FIRST EXECUTABLE STATEMENT  XERHLT 
      STOP 
      END  
 
 
*DECK XERCNT 
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL) 
C***BEGIN PROLOGUE  XERCNT 
C***SUBSIDIARY 
C***PURPOSE  Allow user control over handling of errors. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3C 
C***TYPE      ALL (XERCNT-A) 
C***KEYWORDS  ERROR, XERROR 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C     Abstract 
C        Allows user control over handling of individual errors. 
C        Just after each message is recorded, but before it is 
C        processed any further (i.e., before it is printed or 
C        a decision to abort is made), a call is made to XERCNT. 
C        If the user has provided his own version of XERCNT, he 
C        can then override the value of KONTROL used in processing 
C        this message by redefining its value. 
C        KONTRL may be set to any value from -2 to 2. 
C        The meanings for KONTRL are the same as in XSETF, except 
C        that the value of KONTRL changes only for this message. 
C        If KONTRL is set to a value outside the range from -2 to 2, 
C        it will be moved back into that range. 
C 
C     Description of Parameters 
C 
C      --Input-- 
C        LIBRAR - the library that the routine is in. 
C        SUBROU - the subroutine that XERMSG is being called from 
C        MESSG  - the first 20 characters of the error message. 
C        NERR   - same as in the call to XERMSG. 
C        LEVEL  - same as in the call to XERMSG. 
C        KONTRL - the current value of the control flag as set 
C                 by a call to XSETF. 
C 
C      --Output-- 
C        KONTRL - the new value of KONTRL.  If KONTRL is not 
C                 defined, it will remain at its original value. 
C                 This changed value of control affects only 
C                 the current occurrence of the current message. 
C 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   790801  DATE WRITTEN 
C   861211  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   900206  Routine changed from user-callable to subsidiary.  (WRB) 
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE 
C           names, changed routine name from XERCTL to XERCNT.  (RWC) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XERCNT 
      CHARACTER*(*) LIBRAR, SUBROU, MESSG 
C***FIRST EXECUTABLE STATEMENT  XERCNT 
      RETURN 
      END 
 
 
 
 
 
*DECK FDUMP 
      SUBROUTINE FDUMP 
C***BEGIN PROLOGUE  FDUMP 
C***PURPOSE  Symbolic dump (should be locally written). 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3 
C***TYPE      ALL (FDUMP-A) 
C***KEYWORDS  ERROR, XERMSG 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C        ***Note*** Machine Dependent Routine 
C        FDUMP is intended to be replaced by a locally written 
C        version which produces a symbolic dump.  Failing this, 
C        it should be replaced by a version which prints the 
C        subprogram nesting list.  Note that this dump must be 
C        printed on each of up to five files, as indicated by the 
C        XGETUA routine.  See XSETUA and XGETUA for details. 
C 
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee 
C 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   790801  DATE WRITTEN 
C   861211  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C***END PROLOGUE  FDUMP 
C***FIRST EXECUTABLE STATEMENT  FDUMP 
      RETURN 
      END 
 
 
 
 
*DECK QS2I1D 
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG) 
C***BEGIN PROLOGUE  QS2I1D 
C***SUBSIDIARY 
C***PURPOSE  Sort an integer  array, moving an integer  and DP array. 
C            This routine sorts the integer  array IA and makes the same 
C            interchanges in the integer  array JA and the double pre- 
C            cision array A.  The array IA may be sorted in increasing 
C            order or decreasing order.  A slightly modified QUICKSORT 
C            algorithm is used. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  N6A2A 
C***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D) 
C***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING 
C***AUTHOR  Jones, R. E., (SNLA) 
C           Kahaner, D. K., (NBS) 
C           Seager, M. K., (LLNL) seager@llnl.gov 
C           Wisniewski, J. A., (SNLA) 
C***DESCRIPTION 
C     Written by Rondall E Jones 
C     Modified by John A. Wisniewski to use the Singleton QUICKSORT 
C     algorithm. date 18 November 1976. 
C 
C     Further modified by David K. Kahaner 
C     National Bureau of Standards 
C     August, 1981 
C 
C     Even further modification made to bring the code up to the 
C     Fortran 77 level and make it more readable and to carry 
C     along one integer  array and one double precision array during 
C     the sort by 
C     Mark K. Seager 
C     Lawrence Livermore National Laboratory 
C     November, 1987 
C     This routine was adapted from the ISORT routine. 
C 
C     ABSTRACT 
C         This routine sorts an integer  array IA and makes the same 
C         interchanges in the integer  array JA and the double precision 
C         array A. 
C         The array IA may be sorted in increasing order or decreasing 
C         order.  A slightly modified quicksort algorithm is used. 
C 
C     DESCRIPTION OF PARAMETERS 
C        IA - integer  array of values to be sorted. 
C        JA - integer  array to be carried along. 
C         A - Double Precision array to be carried along. 
C         N - Number of values in integer  array IA to be sorted. 
C     KFLAG - Control parameter 
C           = 1 means sort IA in INCREASING order. 
C           =-1 means sort IA in DECREASING order. 
C 
C***SEE ALSO  DS2Y 
C***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm 
C                 for Sorting With Minimal Storage, Communications ACM 
C                 12:3 (1969), pp.185-7. 
C***ROUTINES CALLED  XERMSG 
C***REVISION HISTORY  (YYMMDD) 
C   761118  DATE WRITTEN 
C   890125  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   900805  Changed XERROR calls to calls to XERMSG.  (RWC) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   910506  Made subsidiary to DS2Y and corrected reference.  (FNF) 
C   920511  Added complete declaration section.  (WRB) 
C   920929  Corrected format of reference.  (FNF) 
C   921012  Corrected all f.p. constants to double precision.  (FNF) 
C***END PROLOGUE  QS2I1D 
CVD$R NOVECTOR 
CVD$R NOCONCUR 
C     .. Scalar Arguments .. 
      INTEGER  KFLAG, N 
C     .. Array Arguments .. 
      DOUBLE PRECISION A(N) 
      INTEGER  IA(N), JA(N) 
C     .. Local Scalars .. 
      DOUBLE PRECISION R, TA, TTA 
      INTEGER  I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN 
C     .. Local Arrays .. 
      INTEGER  IL(21), IU(21) 
C     .. External Subroutines .. 
      EXTERNAL XERMSG 
C     .. Intrinsic Functions .. 
      INTRINSIC ABS, INT 
C***FIRST EXECUTABLE STATEMENT  QS2I1D 
      NN = N 
      IF (NN.LT.1) THEN 
         CALL XERMSG ('SLATEC', 'QS2I1D', 
     $      'The number of values to be sorted was not positive.', 1, 1) 
         RETURN 
      ENDIF 
      IF( N.EQ.1 ) RETURN 
      KK = ABS(KFLAG) 
      IF ( KK.NE.1 ) THEN 
         CALL XERMSG ('SLATEC', 'QS2I1D', 
     $      'The sort control parameter, K, was not 1 or -1.', 2, 1) 
         RETURN 
      ENDIF 
C 
C     Alter array IA to get decreasing order if needed. 
C 
      IF( KFLAG.LT.1 ) THEN 
         DO 20 I=1,NN 
            IA(I) = -IA(I) 
 20      CONTINUE 
      ENDIF 
C 
C     Sort IA and carry JA and A along. 
C     And now...Just a little black magic... 
      M = 1 
      I = 1 
      J = NN 
      R = .375D0 
 210  IF( R.LE.0.5898437D0 ) THEN 
         R = R + 3.90625D-2 
      ELSE 
         R = R-.21875D0 
      ENDIF 
 225  K = I 
C 
C     Select a central element of the array and save it in location 
C     it, jt, at. 
C 
      IJ = I + INT ((J-I)*R) 
      IT = IA(IJ) 
      JT = JA(IJ) 
      TA = A(IJ) 
C 
C     If first element of array is greater than it, interchange with it. 
C 
      IF( IA(I).GT.IT ) THEN 
         IA(IJ) = IA(I) 
         IA(I)  = IT 
         IT     = IA(IJ) 
         JA(IJ) = JA(I) 
         JA(I)  = JT 
         JT     = JA(IJ) 
         A(IJ)  = A(I) 
         A(I)   = TA 
         TA     = A(IJ) 
      ENDIF 
      L=J 
C 
C     If last element of array is less than it, swap with it. 
C 
      IF( IA(J).LT.IT ) THEN 
         IA(IJ) = IA(J) 
         IA(J)  = IT 
         IT     = IA(IJ) 
         JA(IJ) = JA(J) 
         JA(J)  = JT 
         JT     = JA(IJ) 
         A(IJ)  = A(J) 
         A(J)   = TA 
         TA     = A(IJ) 
C 
C     If first element of array is greater than it, swap with it. 
C 
         IF ( IA(I).GT.IT ) THEN 
            IA(IJ) = IA(I) 
            IA(I)  = IT 
            IT     = IA(IJ) 
            JA(IJ) = JA(I) 
            JA(I)  = JT 
            JT     = JA(IJ) 
            A(IJ)  = A(I) 
            A(I)   = TA 
            TA     = A(IJ) 
         ENDIF 
      ENDIF 
C 
C     Find an element in the second half of the array which is 
C     smaller than it. 
C 
  240 L=L-1 
      IF( IA(L).GT.IT ) GO TO 240 
C 
C     Find an element in the first half of the array which is 
C     greater than it. 
C 
  245 K=K+1 
      IF( IA(K).LT.IT ) GO TO 245 
C 
C     Interchange these elements. 
C 
      IF( K.LE.L ) THEN 
         IIT   = IA(L) 
         IA(L) = IA(K) 
         IA(K) = IIT 
         JJT   = JA(L) 
         JA(L) = JA(K) 
         JA(K) = JJT 
         TTA   = A(L) 
         A(L)  = A(K) 
         A(K)  = TTA 
         GOTO 240 
      ENDIF 
C 
C     Save upper and lower subscripts of the array yet to be sorted. 
C 
      IF( L-I.GT.J-K ) THEN 
         IL(M) = I 
         IU(M) = L 
         I = K 
         M = M+1 
      ELSE 
         IL(M) = K 
         IU(M) = J 
         J = L 
         M = M+1 
      ENDIF 
      GO TO 260 
C 
C     Begin again on another portion of the unsorted array. 
C 
  255 M = M-1 
      IF( M.EQ.0 ) GO TO 300 
      I = IL(M) 
      J = IU(M) 
  260 IF( J-I.GE.1 ) GO TO 225 
      IF( I.EQ.J ) GO TO 255 
      IF( I.EQ.1 ) GO TO 210 
      I = I-1 
  265 I = I+1 
      IF( I.EQ.J ) GO TO 255 
      IT = IA(I+1) 
      JT = JA(I+1) 
      TA =  A(I+1) 
      IF( IA(I).LE.IT ) GO TO 265 
      K=I 
  270 IA(K+1) = IA(K) 
      JA(K+1) = JA(K) 
      A(K+1)  =  A(K) 
      K = K-1 
      IF( IT.LT.IA(K) ) GO TO 270 
      IA(K+1) = IT 
      JA(K+1) = JT 
      A(K+1)  = TA 
      GO TO 265 
C 
C     Clean up, if necessary. 
C 
  300 IF( KFLAG.LT.1 ) THEN 
         DO 310 I=1,NN 
            IA(I) = -IA(I) 
 310     CONTINUE 
      ENDIF 
      RETURN 
C------------- LAST LINE OF QS2I1D FOLLOWS ---------------------------- 
      END 
 
 
 
*DECK DSLUI2 
      SUBROUTINE DSLUI2 (N, B, X, IL, JL, L, DINV, IU, JU, U) 
C***BEGIN PROLOGUE  DSLUI2 
C***PURPOSE  SLAP Backsolve for LDU Factorization. 
C            Routine to solve a system of the form  L*D*U X = B, 
C            where L is a unit lower triangular matrix, D is a diagonal 
C            matrix, and U is a unit upper triangular matrix. 
C***LIBRARY   SLATEC (SLAP) 
C***CATEGORY  D2E 
C***TYPE      DOUBLE PRECISION (SSLUI2-S, DSLUI2-D) 
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE, 
C             SLAP, SPARSE 
C***AUTHOR  Greenbaum, Anne, (Courant Institute) 
C           Seager, Mark K., (LLNL) 
C             Lawrence Livermore National Laboratory 
C             PO BOX 808, L-60 
C             Livermore, CA 94550 (510) 423-3141 
C             seager@llnl.gov 
C***DESCRIPTION 
C 
C *Usage: 
C     INTEGER  N, IL(NL), JL(NL), IU(NU), JU(NU) 
C     DOUBLE PRECISION B(N), X(N), L(NL), DINV(N), U(NU) 
C 
C     CALL DSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U ) 
C 
C *Arguments: 
C N      :IN       integer  
C         Order of the Matrix. 
C B      :IN       Double Precision B(N). 
C         Right hand side. 
C X      :OUT      Double Precision X(N). 
C         Solution of L*D*U x = b. 
C IL     :IN       integer  IL(NL). 
C JL     :IN       integer  JL(NL). 
C L      :IN       Double Precision L(NL). 
C         IL, JL, L contain the unit  lower triangular factor of the 
C         incomplete decomposition of some matrix stored in SLAP Row 
C         format.  The diagonal of ones *IS* stored.  This structure 
C         can   be   set  up  by   the  DSILUS  routine.   See   the 
C         "Description", below  for more   details about   the  SLAP 
C         format.  (NL is the number of non-zeros in the L array.) 
C DINV   :IN       Double Precision DINV(N). 
C         Inverse of the diagonal matrix D. 
C IU     :IN       integer  IU(NU). 
C JU     :IN       integer  JU(NU). 
C U      :IN       Double Precision U(NU). 
C         IU, JU, U contain the unit upper triangular factor  of the 
C         incomplete decomposition  of  some  matrix stored in  SLAP 
C         Column format.   The diagonal of ones  *IS* stored.   This 
C         structure can be set up  by the DSILUS routine.  See   the 
C         "Description", below   for  more   details about  the SLAP 
C         format.  (NU is the number of non-zeros in the U array.) 
C 
C *Description: 
C       This routine is supplied with  the SLAP package as a routine 
C       to  perform  the  MSOLVE operation  in   the  SIR and   SBCG 
C       iteration routines for  the  drivers DSILUR and DSLUBC.   It 
C       must  be called  via   the  SLAP  MSOLVE  calling   sequence 
C       convention interface routine DSLUI. 
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** 
C               **** SLAP MSOLVE CALLING CONVENTION **** 
C 
C       IL, JL, L should contain the unit lower triangular factor of 
C       the incomplete decomposition of the A matrix  stored in SLAP 
C       Row format.  IU, JU, U should contain  the unit upper factor 
C       of the  incomplete decomposition of  the A matrix  stored in 
C       SLAP Column format This ILU factorization can be computed by 
C       the DSILUS routine. The diagonals (which are all one's) are 
C       stored. 
C 
C       =================== S L A P Column format ================== 
C 
C       This routine  requires that  the matrix A  be stored in  the 
C       SLAP Column format.  In this format the non-zeros are stored 
C       counting down columns (except for  the diagonal entry, which 
C       must appear first in each  "column")  and are stored  in the 
C       double precision array A.   In other words,  for each column 
C       in the matrix put the diagonal entry in  A.  Then put in the 
C       other non-zero  elements going down  the column (except  the 
C       diagonal) in order.   The  IA array holds the  row index for 
C       each non-zero.  The JA array holds the offsets  into the IA, 
C       A arrays  for  the  beginning  of each   column.   That  is, 
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the 
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), 
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. 
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the 
C       number of columns in  the matrix and NELT  is the number  of 
C       non-zeros in the matrix. 
C 
C       Here is an example of the  SLAP Column  storage format for a 
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column): 
C 
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       ==================== S L A P Row format ==================== 
C 
C       This routine requires  that the matrix A  be  stored  in the 
C       SLAP  Row format.   In this format  the non-zeros are stored 
C       counting across  rows (except for the diagonal  entry, which 
C       must  appear first  in each  "row")  and  are stored  in the 
C       double precision  array A.  In other words, for each row  in 
C       the matrix  put the diagonal  entry in A.   Then put in  the 
C       other  non-zero elements  going across  the row  (except the 
C       diagonal) in order.  The JA array holds the column index for 
C       each non-zero.  The IA array holds the offsets  into the JA, 
C       A  arrays  for  the   beginning  of  each  row.    That  is, 
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW- 
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1) 
C       are  the last elements  of the  IROW-th row.   Note  that we 
C       always have  IA(N+1) = NELT+1, where N is the number of rows 
C       in the matrix  and  NELT is the  number of non-zeros  in the 
C       matrix. 
C 
C       Here is an example of the SLAP Row storage format for a  5x5 
C       Matrix (in the A and JA arrays '|' denotes the end of a row): 
C 
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left. 
C                              1  2  3    4  5    6  7    8    9 10 11 
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53 
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12 
C       | 0  0  0 44  0| 
C       |51  0 53  0 55| 
C 
C       With  the SLAP  format  the "inner  loops" of  this  routine 
C       should vectorize   on machines with   hardware  support  for 
C       vector gather/scatter operations.  Your compiler may require 
C       a  compiler directive  to  convince   it that there  are  no 
C       implicit vector  dependencies.  Compiler directives  for the 
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied 
C       with the standard SLAP distribution. 
C 
C***SEE ALSO  DSILUS 
C***REFERENCES  (NONE) 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   871119  DATE WRITTEN 
C   881213  Previous REVISION DATE 
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS) 
C   890922  Numerous changes to prologue to make closer to SLATEC 
C           standard.  (FNF) 
C   890929  Numerous changes to reduce SP/DP differences.  (FNF) 
C   910411  Prologue converted to Version 4.0 format.  (BAB) 
C   920511  Added complete declaration section.  (WRB) 
C   921113  Corrected C***CATEGORY line.  (FNF) 
C   930701  Updated CATEGORY section.  (FNF, WRB) 
C***END PROLOGUE  DSLUI2 
C     .. Scalar Arguments .. 
      INTEGER  N 
C     .. Array Arguments .. 
      DOUBLE PRECISION B(N), DINV(N), L(*), U(*), X(N) 
      INTEGER  IL(*), IU(*), JL(*), JU(*) 
C     .. Local Scalars .. 
      INTEGER  I, ICOL, IROW, J, JBGN, JEND 
C***FIRST EXECUTABLE STATEMENT  DSLUI2 
C 
C         Solve  L*Y = B,  storing result in X, L stored by rows. 
C 
      DO 10 I = 1, N 
         X(I) = B(I) 
 10   CONTINUE 
      DO 30 IROW = 2, N 
         JBGN = IL(IROW) 
         JEND = IL(IROW+1)-1 
         IF( JBGN.LE.JEND ) THEN 
CLLL. OPTION ASSERT (NOHAZARD) 
CDIR$ IVDEP 
CVD$ ASSOC 
CVD$ NODEPCHK 
            DO 20 J = JBGN, JEND 
               X(IROW) = X(IROW) - L(J)*X(JL(J)) 
 20         CONTINUE 
         ENDIF 
 30   CONTINUE 
C 
C         Solve  D*Z = Y,  storing result in X. 
      DO 40 I=1,N 
         X(I) = X(I)*DINV(I) 
 40   CONTINUE 
C 
C         Solve  U*X = Z, U stored by columns. 
      DO 60 ICOL = N, 2, -1 
         JBGN = JU(ICOL) 
         JEND = JU(ICOL+1)-1 
         IF( JBGN.LE.JEND ) THEN 
CLLL. OPTION ASSERT (NOHAZARD) 
CDIR$ IVDEP 
CVD$ NODEPCHK 
            DO 50 J = JBGN, JEND 
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL) 
 50         CONTINUE 
         ENDIF 
 60   CONTINUE 
C 
      RETURN 
C------------- LAST LINE OF DSLUI2 FOLLOWS ---------------------------- 
      END 
 
 
 
*DECK I1MACH 
      INTEGER  FUNCTION I1MACH (I) 
C***BEGIN PROLOGUE  I1MACH 
C***PURPOSE  Return integer  machine dependent constants. 
C***LIBRARY   SLATEC 
C***CATEGORY  R1 
C***TYPE      INTEGER  (I1MACH-I) 
C***KEYWORDS  MACHINE CONSTANTS 
C***AUTHOR  Fox, P. A., (Bell Labs) 
C           Hall, A. D., (Bell Labs) 
C           Schryer, N. L., (Bell Labs) 
C***DESCRIPTION 
C 
C   I1MACH can be used to obtain machine-dependent parameters for the 
C   local machine environment.  It is a function subprogram with one 
C   (input) argument and can be referenced as follows: 
C 
C        K = I1MACH(I) 
C 
C   where I=1,...,16.  The (output) value of K above is determined by 
C   the (input) value of I.  The results for various values of I are 
C   discussed below. 
C 
C   I/O unit numbers: 
C     I1MACH( 1) = the standard input unit. 
C     I1MACH( 2) = the standard output unit. 
C     I1MACH( 3) = the standard punch unit. 
C     I1MACH( 4) = the standard error message unit. 
C 
C   Words: 
C     I1MACH( 5) = the number of bits per integer  storage unit. 
C     I1MACH( 6) = the number of characters per integer  storage unit. 
C 
C   integer s: 
C     assume integer s are represented in the S-digit, base-A form 
C 
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) 
C 
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1. 
C     I1MACH( 7) = A, the base. 
C     I1MACH( 8) = S, the number of base-A digits. 
C     I1MACH( 9) = A**S - 1, the largest magnitude. 
C 
C   Floating-Point Numbers: 
C     Assume floating-point numbers are represented in the T-digit, 
C     base-B form 
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) 
C 
C                where 0 .LE. X(I) .LT. B for I=1,...,T, 
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. 
C     I1MACH(10) = B, the base. 
C 
C   Single-Precision: 
C     I1MACH(11) = T, the number of base-B digits. 
C     I1MACH(12) = EMIN, the smallest exponent E. 
C     I1MACH(13) = EMAX, the largest exponent E. 
C 
C   Double-Precision: 
C     I1MACH(14) = T, the number of base-B digits. 
C     I1MACH(15) = EMIN, the smallest exponent E. 
C     I1MACH(16) = EMAX, the largest exponent E. 
C 
C   To alter this function for a particular environment, the desired 
C   set of DATA statements should be activated by removing the C from 
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be 
C   checked for consistency with the local operating system. 
C 
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for 
C                 a portable library, ACM Transactions on Mathematical 
C                 Software 4, 2 (June 1978), pp. 177-188. 
C***ROUTINES CALLED  (NONE) 
C***REVISION HISTORY  (YYMMDD) 
C   750101  DATE WRITTEN 
C   891012  Added VAX G-floating constants.  (WRB) 
C   891012  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   900618  Added DEC RISC constants.  (WRB) 
C   900723  Added IBM RS 6000 constants.  (WRB) 
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. 
C           (RWC) 
C   910710  Added HP 730 constants.  (SMR) 
C   911114  Added Convex IEEE constants.  (WRB) 
C   920121  Added SUN -r8 compiler option constants.  (WRB) 
C   920229  Added Touchstone Delta i860 constants.  (WRB) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C   920625  Added Convex -p8 and -pd8 compiler option constants. 
C           (BKS, WRB) 
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) 
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler 
C           options.  (DWL, RWC and WRB). 
C***END PROLOGUE  I1MACH 
C 
      INTEGER  IMACH(16),OUTPUT 
      SAVE IMACH 
      EQUIVALENCE (IMACH(4),OUTPUT) 
C 
C     MACHINE CONSTANTS FOR THE AMIGA 
C     ABSOFT COMPILER 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -126 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1022 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE APOLLO 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        129 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1025 / 
C 
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM 
C 
C     DATA IMACH( 1) /          7 / 
C     DATA IMACH( 2) /          2 / 
C     DATA IMACH( 3) /          2 / 
C     DATA IMACH( 4) /          2 / 
C     DATA IMACH( 5) /         36 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         33 / 
C     DATA IMACH( 9) / Z1FFFFFFFF / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -256 / 
C     DATA IMACH(13) /        255 / 
C     DATA IMACH(14) /         60 / 
C     DATA IMACH(15) /       -256 / 
C     DATA IMACH(16) /        255 / 
C 
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         48 / 
C     DATA IMACH( 6) /          6 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         39 / 
C     DATA IMACH( 9) / O0007777777777777 / 
C     DATA IMACH(10) /          8 / 
C     DATA IMACH(11) /         13 / 
C     DATA IMACH(12) /        -50 / 
C     DATA IMACH(13) /         76 / 
C     DATA IMACH(14) /         26 / 
C     DATA IMACH(15) /        -50 / 
C     DATA IMACH(16) /         76 / 
C 
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         48 / 
C     DATA IMACH( 6) /          6 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         39 / 
C     DATA IMACH( 9) / O0007777777777777 / 
C     DATA IMACH(10) /          8 / 
C     DATA IMACH(11) /         13 / 
C     DATA IMACH(12) /        -50 / 
C     DATA IMACH(13) /         76 / 
C     DATA IMACH(14) /         26 / 
C     DATA IMACH(15) /     -32754 / 
C     DATA IMACH(16) /      32780 / 
C 
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         64 / 
C     DATA IMACH( 6) /          8 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         63 / 
C     DATA IMACH( 9) / 9223372036854775807 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         47 / 
C     DATA IMACH(12) /      -4095 / 
C     DATA IMACH(13) /       4094 / 
C     DATA IMACH(14) /         94 / 
C     DATA IMACH(15) /      -4095 / 
C     DATA IMACH(16) /       4094 / 
C 
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /    6LOUTPUT/ 
C     DATA IMACH( 5) /         60 / 
C     DATA IMACH( 6) /         10 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         48 / 
C     DATA IMACH( 9) / 00007777777777777777B / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         47 / 
C     DATA IMACH(12) /       -929 / 
C     DATA IMACH(13) /       1070 / 
C     DATA IMACH(14) /         94 / 
C     DATA IMACH(15) /       -929 / 
C     DATA IMACH(16) /       1069 / 
C 
C     MACHINE CONSTANTS FOR THE CELERITY C1260 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          0 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / Z'7FFFFFFF' / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -126 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1022 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE CONVEX 
C     USING THE -fn COMPILER OPTION 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1023 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE CONVEX 
C     USING THE -fi COMPILER OPTION 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE CONVEX 
C     USING THE -p8 COMPILER OPTION 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         64 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         63 / 
C     DATA IMACH( 9) / 9223372036854775807 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         53 / 
C     DATA IMACH(12) /      -1023 / 
C     DATA IMACH(13) /       1023 / 
C     DATA IMACH(14) /        113 / 
C     DATA IMACH(15) /     -16383 / 
C     DATA IMACH(16) /      16383 / 
C 
C     MACHINE CONSTANTS FOR THE CONVEX 
C     USING THE -pd8 COMPILER OPTION 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         64 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         63 / 
C     DATA IMACH( 9) / 9223372036854775807 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         53 / 
C     DATA IMACH(12) /      -1023 / 
C     DATA IMACH(13) /       1023 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1023 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE CRAY 
C     USING THE 46 BIT INTEGER  COMPILER OPTION 
C 
C     DATA IMACH( 1) /        100 / 
C     DATA IMACH( 2) /        101 / 
C     DATA IMACH( 3) /        102 / 
C     DATA IMACH( 4) /        101 / 
C     DATA IMACH( 5) /         64 / 
C     DATA IMACH( 6) /          8 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         46 / 
C     DATA IMACH( 9) / 1777777777777777B / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         47 / 
C     DATA IMACH(12) /      -8189 / 
C     DATA IMACH(13) /       8190 / 
C     DATA IMACH(14) /         94 / 
C     DATA IMACH(15) /      -8099 / 
C     DATA IMACH(16) /       8190 / 
C 
C     MACHINE CONSTANTS FOR THE CRAY 
C     USING THE 64 BIT INTEGER  COMPILER OPTION 
C 
C     DATA IMACH( 1) /        100 / 
C     DATA IMACH( 2) /        101 / 
C     DATA IMACH( 3) /        102 / 
C     DATA IMACH( 4) /        101 / 
C     DATA IMACH( 5) /         64 / 
C     DATA IMACH( 6) /          8 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         63 / 
C     DATA IMACH( 9) / 777777777777777777777B / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         47 / 
C     DATA IMACH(12) /      -8189 / 
C     DATA IMACH(13) /       8190 / 
C     DATA IMACH(14) /         94 / 
C     DATA IMACH(15) /      -8099 / 
C     DATA IMACH(16) /       8190 / 
C 
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 
C 
C     DATA IMACH( 1) /         11 / 
C     DATA IMACH( 2) /         12 / 
C     DATA IMACH( 3) /          8 / 
C     DATA IMACH( 4) /         10 / 
C     DATA IMACH( 5) /         16 / 
C     DATA IMACH( 6) /          2 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         15 / 
C     DATA IMACH( 9) /      32767 / 
C     DATA IMACH(10) /         16 / 
C     DATA IMACH(11) /          6 / 
C     DATA IMACH(12) /        -64 / 
C     DATA IMACH(13) /         63 / 
C     DATA IMACH(14) /         14 / 
C     DATA IMACH(15) /        -64 / 
C     DATA IMACH(16) /         63 / 
C 
C     MACHINE CONSTANTS FOR THE DEC ALPHA 
C     USING G_FLOAT 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1023 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE DEC ALPHA 
C     USING IEEE_FLOAT 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE DEC RISC 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE DEC VAX 
C     USING D_FLOATING 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         56 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE DEC VAX 
C     USING G_FLOATING 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1023 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE ELXSI 6400 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         32 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -126 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1022 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE HARRIS 220 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          0 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         24 / 
C     DATA IMACH( 6) /          3 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         23 / 
C     DATA IMACH( 9) /    8388607 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         23 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         38 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /         43 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         36 / 
C     DATA IMACH( 6) /          6 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         35 / 
C     DATA IMACH( 9) / O377777777777 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         27 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         63 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE HP 730 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE HP 2100 
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          4 / 
C     DATA IMACH( 4) /          1 / 
C     DATA IMACH( 5) /         16 / 
C     DATA IMACH( 6) /          2 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         15 / 
C     DATA IMACH( 9) /      32767 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         23 / 
C     DATA IMACH(12) /       -128 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         39 / 
C     DATA IMACH(15) /       -128 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE HP 2100 
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          4 / 
C     DATA IMACH( 4) /          1 / 
C     DATA IMACH( 5) /         16 / 
C     DATA IMACH( 6) /          2 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         15 / 
C     DATA IMACH( 9) /      32767 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         23 / 
C     DATA IMACH(12) /       -128 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         55 / 
C     DATA IMACH(15) /       -128 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE HP 9000 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          7 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         32 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -126 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1015 / 
C     DATA IMACH(16) /       1017 / 
C 
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, 
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND 
C     THE PERKIN ELMER (INTERDATA) 7/32. 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          7 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) /  Z7FFFFFFF / 
C     DATA IMACH(10) /         16 / 
C     DATA IMACH(11) /          6 / 
C     DATA IMACH(12) /        -64 / 
C     DATA IMACH(13) /         63 / 
C     DATA IMACH(14) /         14 / 
C     DATA IMACH(15) /        -64 / 
C     DATA IMACH(16) /         63 / 
C 
C     MACHINE CONSTANTS FOR THE IBM PC 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          0 / 
C     DATA IMACH( 4) /          0 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE IBM RS 6000 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          0 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE INTEL i860 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         36 / 
C     DATA IMACH( 6) /          5 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         35 / 
C     DATA IMACH( 9) / "377777777777 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         27 / 
C     DATA IMACH(12) /       -128 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         54 / 
C     DATA IMACH(15) /       -101 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         36 / 
C     DATA IMACH( 6) /          5 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         35 / 
C     DATA IMACH( 9) / "377777777777 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         27 / 
C     DATA IMACH(12) /       -128 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         62 / 
C     DATA IMACH(15) /       -128 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING 
C     32-BIT INTEGER  ARITHMETIC. 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         56 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING 
C     16-BIT INTEGER  ARITHMETIC. 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          5 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         16 / 
C     DATA IMACH( 6) /          2 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         15 / 
C     DATA IMACH( 9) /      32767 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         56 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE SUN 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -125 / 
C     DATA IMACH(13) /        128 / 
C     DATA IMACH(14) /         53 / 
C     DATA IMACH(15) /      -1021 / 
C     DATA IMACH(16) /       1024 / 
C 
C     MACHINE CONSTANTS FOR THE SUN 
C     USING THE -r8 COMPILER OPTION 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          6 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         32 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         31 / 
C     DATA IMACH( 9) / 2147483647 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         53 / 
C     DATA IMACH(12) /      -1021 / 
C     DATA IMACH(13) /       1024 / 
C     DATA IMACH(14) /        113 / 
C     DATA IMACH(15) /     -16381 / 
C     DATA IMACH(16) /      16384 / 
C 
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER 
C 
C     DATA IMACH( 1) /          5 / 
C     DATA IMACH( 2) /          6 / 
C     DATA IMACH( 3) /          1 / 
C     DATA IMACH( 4) /          6 / 
C     DATA IMACH( 5) /         36 / 
C     DATA IMACH( 6) /          4 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         35 / 
C     DATA IMACH( 9) / O377777777777 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         27 / 
C     DATA IMACH(12) /       -128 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         60 / 
C     DATA IMACH(15) /      -1024 / 
C     DATA IMACH(16) /       1023 / 
C 
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR 
C 
C     DATA IMACH( 1) /          1 / 
C     DATA IMACH( 2) /          1 / 
C     DATA IMACH( 3) /          0 / 
C     DATA IMACH( 4) /          1 / 
C     DATA IMACH( 5) /         16 / 
C     DATA IMACH( 6) /          2 / 
C     DATA IMACH( 7) /          2 / 
C     DATA IMACH( 8) /         15 / 
C     DATA IMACH( 9) /      32767 / 
C     DATA IMACH(10) /          2 / 
C     DATA IMACH(11) /         24 / 
C     DATA IMACH(12) /       -127 / 
C     DATA IMACH(13) /        127 / 
C     DATA IMACH(14) /         56 / 
C     DATA IMACH(15) /       -127 / 
C     DATA IMACH(16) /        127 / 
C 
C***FIRST EXECUTABLE STATEMENT  I1MACH 
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10 
C 
      I1MACH = IMACH(I) 
      RETURN 
C 
   10 CONTINUE 
      WRITE (UNIT = OUTPUT, FMT = 9000) 
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS') 
C 
C     CALL FDUMP 
C 
      STOP 
      END 
 
 
 
*DECK XGETUA 
      SUBROUTINE XGETUA (IUNITA, N) 
C***BEGIN PROLOGUE  XGETUA 
C***PURPOSE  Return unit number(s) to which error messages are being 
C            sent. 
C***LIBRARY   SLATEC (XERROR) 
C***CATEGORY  R3C 
C***TYPE      ALL (XGETUA-A) 
C***KEYWORDS  ERROR, XERROR 
C***AUTHOR  Jones, R. E., (SNLA) 
C***DESCRIPTION 
C 
C     Abstract 
C        XGETUA may be called to determine the unit number or numbers 
C        to which error messages are being sent. 
C        These unit numbers may have been set by a call to XSETUN, 
C        or a call to XSETUA, or may be a default value. 
C 
C     Description of Parameters 
C      --Output-- 
C        IUNIT - an array of one to five unit numbers, depending 
C                on the value of N.  A value of zero refers to the 
C                default unit, as defined by the I1MACH machine 
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are 
C                defined by XGETUA.  The values of IUNIT(N+1),..., 
C                IUNIT(5) are not defined (for N .LT. 5) or altered 
C                in any way by XGETUA. 
C        N     - the number of units to which copies of the 
C                error messages are being sent.  N will be in the 
C                range from 1 to 5. 
C 
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC 
C                 Error-handling Package, SAND82-0800, Sandia 
C                 Laboratories, 1982. 
C***ROUTINES CALLED  J4SAVE 
C***REVISION HISTORY  (YYMMDD) 
C   790801  DATE WRITTEN 
C   861211  REVISION DATE from Version 3.2 
C   891214  Prologue converted to Version 4.0 format.  (BAB) 
C   920501  Reformatted the REFERENCES section.  (WRB) 
C***END PROLOGUE  XGETUA 
      DIMENSION IUNITA(5) 
C***FIRST EXECUTABLE STATEMENT  XGETUA 
      N = J4SAVE(5,0,.FALSE.) 
      DO 30 I=1,N 
         INDEX = I+4 
         IF (I.EQ.1) INDEX = 3 
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.) 
   30 CONTINUE 
      RETURN 
      END 
 


c	Subroutine per la risoluzione del sistema lineare 
 

 
 
 
 
 
      SUBROUTINE MATVEC(N,p,q,nelt,ia,ja,A,isym) 
      IMPLICIT  REAL*8(A-H,O-Z) 
  
      double precision A(Nelt)  
      integer  ia(Nelt),ja(Nelt)  
      integer  i,j,k  
      double precision p(n),q(n) 
      integer  l,inizio,isym  
  

c      write(6,*) "sono dentro matvec, nelt=", nelt

      do i=1,n  
      q(i)=0  
      enddo  
  
       
          
      inizio=1  
      k=1
c      write(6,*) ia(k)-1
  
      do i=1,n  
          do l=inizio,ia(k)-1  
c             write(6,*) ja(l),l  
             q(i)=q(i)+p(ja(l))*A(l)  
          enddo  
 
          inizio=ia(k)  
          k=k+1  
          enddo  
            
          return  
          end  
 
	SUBROUTINE MSOLVE( N,X,Y, nelt,ia,ja,A,isym,rwork,iwork) 
 
	IMPLICIT           REAL*8(A-H,O-Z) 
	DOUBLE PRECISION   X(N),Y(N) 
 
	double precision A(Nelt)  
      integer ia(Nelt),ja(Nelt)  
 
      integer i,j,k,m  
     
      integer l,inizio,isym  
 
 
 
C 
C     MSOLVE  solves  M*Y = X  for some symmetric pos-def matrix  M. 
C     This is a simple example for testing  SYMMLQ. 
C 
	 
 
	j=1 
	m=1 
	DO 10 I = 1, Nelt 
        if (i.lt.ia(m)) then 
		if (ja(i).eq.m) then 
						y(j)=x(j)/A(i) 
						j=j+1 
	endif 
		else 
			m=m+1 
		if (ja(i).eq.m) then 
			y(j)=x(j)/A(i) 
			j=j+1 
	endif 
	endif 
 
 
10	CONTINUE 
      RETURN 
 
 
 
C 
C     END OF MSOLVE 
      END 





c$$$	SUBROUTINE MSOLVE(N,X,Y,nelt,ia,ja,A,isym,rwork,iwork) 
c$$$ 
c$$$	IMPLICIT           REAL*8(A-H,O-Z) 
c$$$	DOUBLE PRECISION   X(N),Y(N) 
c$$$ 
c$$$	double precision A(Nelt)  
c$$$      integer  ia(Nelt),ja(Nelt)  
c$$$ 
c$$$      integer  i,j,k,m  
c$$$     
c$$$      integer  l,inizio,isym  
c$$$ 
c$$$ 
c$$$ 
c$$$C 
c$$$C     MSOLVE  solves  M*Y = X  for some symmetric pos-def matrix  M. 
c$$$C     This is a simple example for testing  SYMMLQ. 
c$$$C 
c$$$	 
c$$$      do i=1,N
c$$$         Y(i)=0.0
c$$$      enddo
c$$$
c$$$	j=1 
c$$$	m=1 
c$$$	DO 10 I = 1, Nelt 
c$$$        if (i.lt.ia(m)) then 
c$$$		if (ja(i).eq.m) then 
c$$$                   y(j)=x(j)/A(i) 
c$$$                   j=j+1 
c$$$                endif 
c$$$             else 
c$$$                m=m+1 
c$$$		if (ja(i).eq.m) then 
c$$$                   y(j)=x(j)/A(i) 
c$$$                   j=j+1 
c$$$                endif 
c$$$             endif 
c$$$             
c$$$             
c$$$ 10       CONTINUE 
c$$$          RETURN 
c$$$ 
c$$$ 
c$$$ 
c$$$C 
c$$$C     END OF MSOLVE 
c$$$      END 

c$$$	SUBROUTINE MSOLVE(N,X,Y,nelt,ia,ja,A,isym,rwork,iwork) 
c$$$ 
c$$$	IMPLICIT           REAL*8(A-H,O-Z) 
c$$$	DOUBLE PRECISION   X(N),Y(N) 
c$$$ 
c$$$	double precision A(Nelt),c2(N),c1(N),c3(N)  
c$$$      integer*4  ia(Nelt),ja(Nelt),kk,prec  
c$$$ 
c$$$      integer*4  i,j,k,m  
c$$$     
c$$$      integer*4  l,inizio,isym  
c$$$ 
c$$$ 
c$$$ 
c$$$C 
c$$$C     MSOLVE  solves  M*Y = X  for some symmetric pos-def matrix  M. 
c$$$C     This is a simple example for testing  SYMMLQ. 
c$$$C 
c$$$	 
c$$$      do i=1,N
c$$$         Y(i)=0.0
c$$$      enddo
c$$$      
c$$$      kk=1
c$$$      prec=1
c$$$      do k=1,N
c$$$         do j=1,(ia(k)-prec)
c$$$            if (k.eq.ja(kk+j-1)) c2(k)=a(kk+j-1)
c$$$            if (ja(kk+j-1).eq.(k+1)) c3(k)=a(kk+j-1)
c$$$            if (ja(kk+j-1).eq.(k-1)) c1(k)=a(kk+j-1)
c$$$         enddo
c$$$         prec=ia(k)
c$$$         kk=ia(k)
c$$$      enddo
c$$$            
c$$$      call tridag(c1,c2,c3,X,Y,N)
c$$$      
c$$$      RETURN 
c$$$C     END OF MSOLVE 
c$$$      END 
c$$$ 


*DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END






c$$$      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY) 
c$$$C 
c$$$C     COPIES A VECTOR, X, TO A VECTOR, Y. 
c$$$C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. 
c$$$C     JACK DONGARRA, LINPACK, 3/11/78. 
c$$$C 
c$$$      DOUBLE PRECISION DX(1),DY(1) 
c$$$      integer *4 I,INCX,INCY,IX,IY,M,MP1,N 
c$$$C 
c$$$      IF(N.LE.0)RETURN 
c$$$      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 
c$$$C 
c$$$C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS 
c$$$C          NOT EQUAL TO 1 
c$$$C 
c$$$      IX = 1 
c$$$      IY = 1 
c$$$      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
c$$$      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
c$$$      DO 10 I = 1,N 
c$$$        DY(IY) = DX(IX) 
c$$$        IX = IX + INCX 
c$$$        IY = IY + INCY 
c$$$   10 CONTINUE 
c$$$      RETURN 
c$$$C 
c$$$C        CODE FOR BOTH INCREMENTS EQUAL TO 1 
c$$$C 
c$$$C 
c$$$C        CLEAN-UP LOOP 
c$$$C 
c$$$   20 M = MOD(N,7) 
c$$$      IF( M .EQ. 0 ) GO TO 40 
c$$$      DO 30 I = 1,M 
c$$$        DY(I) = DX(I) 
c$$$   30 CONTINUE 
c$$$      IF( N .LT. 7 ) RETURN 
c$$$   40 MP1 = M + 1 
c$$$      DO 50 I = MP1,N,7 
c$$$        DY(I) = DX(I) 
c$$$        DY(I + 1) = DX(I + 1) 
c$$$        DY(I + 2) = DX(I + 2) 
c$$$        DY(I + 3) = DX(I + 3) 
c$$$        DY(I + 4) = DX(I + 4) 
c$$$        DY(I + 5) = DX(I + 5) 
c$$$        DY(I + 6) = DX(I + 6) 
c$$$   50 CONTINUE 
c$$$      RETURN 
c$$$      END 
c$$$ 










*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END




















 
 
c$$$ 
c$$$ 
c$$$      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY) 
c$$$C 
c$$$C     FORMS THE DOT PRODUCT OF TWO VECTORS. 
c$$$C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. 
c$$$C     JACK DONGARRA, LINPACK, 3/11/78. 
c$$$C 
c$$$      DOUBLE PRECISION DX(1),DY(1),DTEMP 
c$$$      integer *4 I,INCX,INCY,IX,IY,M,MP1,N 
c$$$C 
c$$$      DDOT = 0.0D0 
c$$$      DTEMP = 0.0D0 
c$$$      IF(N.LE.0)RETURN 
c$$$      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 
c$$$C 
c$$$C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS 
c$$$C          NOT EQUAL TO 1 
c$$$C 
c$$$      IX = 1 
c$$$      IY = 1 
c$$$      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
c$$$      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
c$$$      DO 10 I = 1,N 
c$$$        DTEMP = DTEMP + DX(IX)*DY(IY) 
c$$$        IX = IX + INCX 
c$$$        IY = IY + INCY 
c$$$   10 CONTINUE 
c$$$      DDOT = DTEMP 
c$$$      RETURN 
c$$$C 
c$$$C        CODE FOR BOTH INCREMENTS EQUAL TO 1 
c$$$C 
c$$$C 
c$$$C        CLEAN-UP LOOP 
c$$$C 
c$$$   20 M = MOD(N,5) 
c$$$      IF( M .EQ. 0 ) GO TO 40 
c$$$      DO 30 I = 1,M 
c$$$        DTEMP = DTEMP + DX(I)*DY(I) 
c$$$   30 CONTINUE 
c$$$      IF( N .LT. 5 ) GO TO 60 
c$$$   40 MP1 = M + 1 
c$$$      DO 50 I = MP1,N,5 
c$$$        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + 
c$$$     $   DX(I + 2)*DY(I + 2)+DX(I+3)*DY(I + 3) + DX(I + 4)*DY(I + 4) 
c$$$   50 CONTINUE 
c$$$   60 DDOT = DTEMP 
c$$$      RETURN 
c$$$      END 





*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END

















c$$$      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY) 
c$$$
c$$$C 
c$$$
c$$$C     CONSTANT TIMES A VECTOR PLUS A VECTOR. 
c$$$
c$$$C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. 
c$$$
c$$$C     JACK DONGARRA, LINPACK, 3/11/78. 
c$$$
c$$$C 
c$$$
c$$$      DOUBLE PRECISION DX(1),DY(1),DA 
c$$$
c$$$      INTEGER  I,INCX,INCY,IX,IY,M,MP1,N 
c$$$
c$$$C 
c$$$
c$$$      IF(N.LE.0)RETURN 
c$$$
c$$$      IF (DA .EQ. 0.0D0) RETURN 
c$$$
c$$$      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 
c$$$
c$$$C 
c$$$
c$$$C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS 
c$$$
c$$$C          NOT EQUAL TO 1 
c$$$
c$$$C 
c$$$
c$$$      IX = 1 
c$$$
c$$$      IY = 1 
c$$$
c$$$      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
c$$$
c$$$      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
c$$$
c$$$      DO 10 I = 1,N 
c$$$
c$$$        DY(IY) = DY(IY) + DA*DX(IX) 
c$$$
c$$$        IX = IX + INCX 
c$$$
c$$$        IY = IY + INCY 
c$$$
c$$$   10 CONTINUE 
c$$$
c$$$      RETURN 
c$$$
c$$$C 
c$$$
c$$$C        CODE FOR BOTH INCREMENTS EQUAL TO 1 
c$$$
c$$$C 
c$$$
c$$$C 
c$$$
c$$$C        CLEAN-UP LOOP 
c$$$
c$$$C 
c$$$
c$$$   20 M = MOD(N,4) 
c$$$
c$$$      IF( M .EQ. 0 ) GO TO 40 
c$$$
c$$$      DO 30 I = 1,M 
c$$$
c$$$        DY(I) = DY(I) + DA*DX(I) 
c$$$
c$$$   30 CONTINUE 
c$$$
c$$$      IF( N .LT. 4 ) RETURN 
c$$$
c$$$   40 MP1 = M + 1 
c$$$
c$$$      DO 50 I = MP1,N,4 
c$$$
c$$$        DY(I) = DY(I) + DA*DX(I) 
c$$$
c$$$        DY(I + 1) = DY(I + 1) + DA*DX(I + 1) 
c$$$
c$$$        DY(I + 2) = DY(I + 2) + DA*DX(I + 2) 
c$$$
c$$$        DY(I + 3) = DY(I + 3) + DA*DX(I + 3) 
c$$$
c$$$   50 CONTINUE 
c$$$
c$$$      RETURN 
c$$$
c$$$      END 
c$$$
c$$$
c$$$
c$$$



















*DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)
C***BEGIN PROLOGUE  DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNRM2
      INTEGER  NEXT
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,
     +                 ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  DNRM2
      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = I,NN,INCX
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = SQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END








c$$$
c$$$ 
c$$$
c$$$      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX) 
c$$$
c$$$      INTEGER  I, INCX, J, N, NEXT, NN 
c$$$
c$$$      DOUBLE PRECISION  DX(1),CUTLO,CUTHI,HITEST, SUM, XMAX,ZERO,ONE 
c$$$
c$$$      DATA   ZERO, ONE /0.0D0, 1.0D0/ 
c$$$
c$$$C 
c$$$
c$$$C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE 
c$$$
c$$$C     INCREMENT INCX . 
c$$$
c$$$C     IF    N .LE. 0 RETURN WITH RESULT = 0. 
c$$$
c$$$C     IF N .GE. 1 THEN INCX MUST BE .GE. 1 
c$$$
c$$$C 
c$$$
c$$$C           C.L.LAWSON, 1978 JAN 08 
c$$$
c$$$C 
c$$$
c$$$C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE 
c$$$
c$$$C     HOPEFULLY APPLICABLE TO ALL MACHINES. 
c$$$
c$$$C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES. 
c$$$
c$$$C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES. 
c$$$
c$$$C     WHERE 
c$$$
c$$$C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. 
c$$$
c$$$C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) 
c$$$
c$$$C         V   = LARGEST  NO.            (OVERFLOW  LIMIT) 
c$$$
c$$$C 
c$$$
c$$$C     BRIEF OUTLINE OF ALGORITHM.. 
c$$$
c$$$C 
c$$$
c$$$C     PHASE 1    SCANS ZERO COMPONENTS. 
c$$$
c$$$C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO 
c$$$
c$$$C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO 
c$$$
c$$$C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M 
c$$$
c$$$C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. 
c$$$
c$$$C 
c$$$
c$$$C     VALUES FOR CUTLO AND CUTHI.. 
c$$$
c$$$C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER 
c$$$
c$$$C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. 
c$$$
c$$$C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE 
c$$$
c$$$C                   UNIVAC AND DEC AT 2**(-103) 
c$$$
c$$$C                   THUS CUTLO = 2**(-51) = 4.44089E-16 
c$$$
c$$$C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. 
c$$$
c$$$C                   THUS CUTHI = 2**(63.5) = 1.30438E19 
c$$$
c$$$C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. 
c$$$
c$$$C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 
c$$$
c$$$C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 
c$$$
c$$$C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / 
c$$$
c$$$C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / 
c$$$
c$$$      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / 
c$$$
c$$$C 
c$$$
c$$$      IF(N .GT. 0) GO TO 10 
c$$$
c$$$         DNRM2  = ZERO 
c$$$
c$$$         GO TO 300 
c$$$
c$$$C 
c$$$
c$$$   10 ASSIGN 30 TO NEXT 
c$$$
c$$$      SUM = ZERO 
c$$$
c$$$      NN = N * INCX 
c$$$
c$$$C                                                 BEGIN MAIN LOOP 
c$$$
c$$$      I = 1 
c$$$
c$$$   20    GO TO NEXT,(30, 50, 70, 110) 
c$$$
c$$$   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85 
c$$$
c$$$      ASSIGN 50 TO NEXT 
c$$$
c$$$      XMAX = ZERO 
c$$$
c$$$C 
c$$$
c$$$C                        PHASE 1.  SUM IS ZERO 
c$$$
c$$$C 
c$$$
c$$$   50 IF( DX(I) .EQ. ZERO) GO TO 200 
c$$$
c$$$      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85 
c$$$
c$$$C 
c$$$
c$$$C                                PREPARE FOR PHASE 2. 
c$$$
c$$$      ASSIGN 70 TO NEXT 
c$$$
c$$$      GO TO 105 
c$$$
c$$$C 
c$$$
c$$$C                                PREPARE FOR PHASE 4. 
c$$$
c$$$C 
c$$$
c$$$  100 I = J 
c$$$
c$$$      ASSIGN 110 TO NEXT 
c$$$
c$$$      SUM = (SUM / DX(I)) / DX(I) 
c$$$
c$$$  105 XMAX = DABS(DX(I)) 
c$$$
c$$$      GO TO 115 
c$$$
c$$$C 
c$$$
c$$$C                   PHASE 2.  SUM IS SMALL. 
c$$$
c$$$C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. 
c$$$
c$$$C 
c$$$
c$$$   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75 
c$$$
c$$$C 
c$$$
c$$$C                     COMMON CODE FOR PHASES 2 AND 4. 
c$$$
c$$$C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. 
c$$$
c$$$C 
c$$$
c$$$  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115 
c$$$
c$$$         SUM = ONE + SUM * (XMAX / DX(I))**2 
c$$$
c$$$         XMAX = DABS(DX(I)) 
c$$$
c$$$         GO TO 200 
c$$$
c$$$C 
c$$$
c$$$  115 SUM = SUM + (DX(I)/XMAX)**2 
c$$$
c$$$      GO TO 200 
c$$$
c$$$C 
c$$$
c$$$C 
c$$$
c$$$C                  PREPARE FOR PHASE 3. 
c$$$
c$$$C 
c$$$
c$$$   75 SUM = (SUM * XMAX) * XMAX 
c$$$
c$$$C 
c$$$
c$$$C 
c$$$
c$$$C     FOR REAL OR D.P. SET HITEST = CUTHI/N 
c$$$
c$$$C     FOR COMPLEX      SET HITEST = CUTHI/(2*N) 
c$$$
c$$$C 
c$$$
c$$$   85 HITEST = CUTHI/FLOAT( N ) 
c$$$
c$$$C 
c$$$
c$$$C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. 
c$$$
c$$$C 
c$$$
c$$$      DO 95 J =I,NN,INCX 
c$$$
c$$$      IF(DABS(DX(J)) .GE. HITEST) GO TO 100 
c$$$
c$$$   95    SUM = SUM + DX(J)**2 
c$$$
c$$$      DNRM2 = DSQRT( SUM ) 
c$$$
c$$$      GO TO 300 
c$$$
c$$$C 
c$$$
c$$$  200 CONTINUE 
c$$$
c$$$      I = I + INCX 
c$$$
c$$$      IF ( I .LE. NN ) GO TO 20 
c$$$
c$$$C 
c$$$
c$$$C              END OF MAIN LOOP. 
c$$$
c$$$C 
c$$$
c$$$C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. 
c$$$
c$$$C 
c$$$
c$$$      DNRM2 = XMAX * DSQRT(SUM) 
c$$$
c$$$  300 CONTINUE 
c$$$
c$$$      RETURN 
c$$$
c$$$      END  
c$$$ 
 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
     

      SUBROUTINE tridag(a,b,c,r,u,n) 
      INTEGER n,NMAX 
      DOUBLE PRECISION a(n),b(n),c(n),r(n),u(n) 
      PARAMETER (NMAX=30000) 
      INTEGER j 
      DOUBLE PRECISION bet,gam(NMAX) 
      if(b(1).eq.0.)pause '*********tridag: rewrite equations***' 
      bet=b(1) 
      u(1)=r(1)/bet 
      do j=2,n 
        gam(j)=c(j-1)/bet 
        bet=b(j)-a(j)*gam(j) 
        if(bet.eq.0.) then
           write(6,*) '********inner iteration=',j
           pause '************tridag failed**********' 
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet 
      end do
      do j=n-1,1,-1 
        u(j)=u(j)-gam(j+1)*u(j+1) 
      end do
      return
      END            


c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 
 
 
      subroutine y12maf(n, z, a, snr, nn, rnr, nn1, pivot, ha,
     &     iha,aflag,iflag,b,ifail)
      implicit double precision (a-b,g,p,t-y),integer*4(c,f,h-n,r-s,z)
      double precision a(nn), pivot(n), aflag(8),b(n)
      integer*4 snr(nn), rnr(nn1)
      integer*4 iflag(10),ha(iha,11)
      aflag(1)=16
      aflag(2)=1.0e-50
c	aflag(2)=1.0e-24
      aflag(3)=1.0e16
      aflag(4)=1.0e-90
c	come era originariamente
c	aflag(4)=1.0e-24
      iflag(2)=3
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1
      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 1
      call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 1
      call y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
1     return
      end




      subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag,
     1     iflag,ifail)
c
c
c  the non-zero elements of a sparse matrix a are prepared  in order to
c  solve the system ax=b by use of sparse matrix technique/
c
c
      implicit double precision (a-b,g,p,t-y),integer*4 (c,f,h-n,r-s,z)
      double precision a(nn), aflag(8)
      integer*4 snr(nn), rnr(nn1)
      integer*4 iflag(10), ha(iha,11)
      mode=iflag(4)
      ifail=0
      if(n.lt.2)ifail=12
      if(z.le.0)ifail=13
      if(nn.lt.2*z)ifail=5
      if(nn1.lt.z)ifail=6
      if(ifail.eq.0.and.n.gt.z)ifail=14
      if(iha.lt.n)ifail=15
      if(mode.lt.0)ifail=16
      if(mode.gt.2)ifail=16
      if(ifail.ne.0) go to 22
      gt1=0.0d0
      do 10 i=1,n
      ha(i,2)=0
      ha(i,3)=0
   10 ha(i,6)=0
c
c  find the number of the non-zero elements in each row and column;move
c  the non-zero elements in the end of the arrays a and snr;find the
c  largest non-zero element in a(in absolute value).
c
      do 20 i=1,z
      t=dabs(a(i))
      l3=rnr(i)
      l4=snr(i)
      if(l4.gt.n.or.l4.lt.1)ifail=24
      if(l3.gt.n.or.l3.lt.1)ifail=25
      ha(l3,3)=ha(l3,3)+1
      ha(l4,6)=ha(l4,6)+1
      if(t.gt.gt1)gt1=t
      a(z+i)=a(i)
   20 snr(z+i)=snr(i)
      if(ifail.gt.0)go to 22
c
c  store the information of the row starts(in ha(i,1))and of the column
c  starts(in ha(i,4)).
c
      l1=1
      l2=1
      do 40 i=1,n
      l3=ha(i,3)
      l4=ha(i,6)
      if(l3.gt.0)go to 21
      ifail=17
      go to 22
   21 if(l4.gt.0)go to 23
      ifail=18
      go to 22
   23 if(mode.eq.2)go to 30
      ha(i,9)=l3
      ha(i,10)=l4
      ha(i,11)=0
      ha(l3,2)=ha(l3,2)+1
      ha(i,5)=l3
   30 ha(i,1)=l1
      ha(i,4)=l2
      l1=l1+l3
      l2=l2+l4
      ha(i,3)=0
   40 ha(i,6)=0
c
c  store the non-zero elements of matrix a(ordered in rows) in the
c  first z locations of the array a.do the same for their column numbers
c
      do 50 i=1,z
      l1=z+i
      l3=rnr(i)
      l2=ha(l3,1)+ha(l3,3)
      a(l2)=a(l1)
      snr(l2)=snr(l1)
   50 ha(l3,3)=ha(l3,3)+1
c
c  store the row numbers of the non-zero elements ordered by columns in
c  the first z locations of the array rnr. store information about row
c  ends(in ha(i,3)).
c
      l4=1
      do 70 i=1,n
      if(mode.eq.2)go to 60
      if(ha(i,2).eq.0)go to 60
      ha(i,11)=l4
      l4=l4+ha(i,2)
      ha(i,2)=ha(i,11)
   60 ha(i,3)=ha(i,1)+ha(i,3)-1
      l1=ha(i,1)
      l2=ha(i,3)
      do 70 j=l1,l2
      l3=snr(j)
      r=ha(l3,6)
      index=ha(l3,4)+r
      rnr(index)=i
      if(r.eq.0)go to 70
      if(j.eq.l1)go to 70
      if(rnr(index-1).ne.i)go to 70
      ifail=11
      go to 22
   70 ha(l3,6)=r+1
      do 90 i=1,n
      if(mode.eq.2)go to 80
      l3=ha(i,5)
      l5=ha(l3,2)
      ha(l5,8)=i
      ha(i,7)=l5
      ha(l3,2)=ha(l3,2)+1
   80 continue
   90 ha(i,6)=ha(i,4)+ha(i,6)-1
      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1
22    return
      end


      subroutine y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha, aflag,iflag
     *,ifail)
c
c  systens of linear equations are solved by use of sparse matrix tech-
c  nique and by gaussian elimination.
c
      implicit double precision (a-b,g,p,t-y),integer*4 (c,f,h-n,r-s,z)
      double precision a(nn),b(n),pivot(n),aflag(8)
c
c  information which is necessary to begin the elimination is stored.
c
      integer*4 snr(nn),rnr(nn1)
      integer*4 iflag(10),ha(iha,11)
      ifail=0
      if(iflag(1).ne.-1)ifail=2
      if(aflag(1).lt.1.0d0)aflag(1)=1.0005 d0
      if(aflag(3).lt.1.0d+5)aflag(3)=1.0d+5
      if(aflag(4).lt.0.0d0)aflag(4)=-aflag(4)
      if(iflag(2).lt.1)ifail=19
      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20
      if(iflag(5).lt.1.or.iflag(5).gt.3)ifail=21
      if(iflag(5).eq.3)ifail=22
      if(ifail.gt.0)go to 1110
      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
c
c  use the information about fill-ins if it is possible.
c
      zz=z
      nr=n*n
      if(iflag(4).ne.2)go to 100
      if(iflag(10).gt.nn)go to 50
      l1=iflag(10)
      l5=l1+1
      if(l5.le.nn)snr(l5)=0
      do 40 i=1,n
      l=n8-i
      l2=ha(l,3)+1
      l3=l2-ha(l,1)
      do 10 j=1,l3
      snr(l5-j)=snr(l2-j)
   10 a(l5-j)=a(l2-j)
      ha(l,3)=l1
      ha(l,1)=l5-l3
      l6=l1-l3
      l5=l5-ha(l,9)
      if(l5.gt.l6)go to 30
      do 20 j=l5,l6
   20 snr(j)=0
   30 continue
   40 l1=l5-1
   50 if(iflag(9).gt.nn1)go to 100
      l2=iflag(9)
      l5=l2+1
      if(l5.le.nn1)rnr(l5)=0
      do 90 i=1,n
      l=n8-i
      l1=ha(l,6)+1
      l4=l1-ha(l,4)
      do 60 j=1,l4
   60 rnr(l5-j)=rnr(l1-j)
      ha(l,4)=l5-l4
      ha(l,6)=l2
      l6=l2-l4
      l5=l5-ha(l,10)
      if(l5.gt.l6)go to 80
      do 70 j=l5,l6
   70 rnr(j)=0
   80 continue
   90 l2=l5-1
  100 r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)
      do 110 i=1,n
      pivot(i)=0.0 d0
      ha(i,2)=ha(i,1)
  110 ha(i,5)=ha(i,4)
      index=ha(n,8)
c
c  start of gaussian elimination.
c
      slut=ha(index,3)-ha(index,2)+1
      do 950 i=1,n7
      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350
      if(iflag(4).ne.2)go to 120
      rrow=ha(i,7)
      rcoll=ha(i,8)
      go to 220
  120 l4=ha(i,8)
      if(iflag(3).eq.1)go to 130
      rrow=l4
      rcoll=rrow
      rpivot=i
      go to 170
  130 r=nr
      v=0.0 d0
      index=iflag(2)
      do 160 kk=1,index
      l1=i-1+kk
      if(l1.gt.n)go to 170
      j=ha(l1,8)
      r7=ha(j,2)
      r8=ha(j,3)
      r9=r8-r7
      t=0.0 d0
      do 140 k=r7,r8
      td=dabs(a(k))
  140 if(t.lt.td)t=td
      t=t/u
      do 160 k=r7,r8
      td=dabs(a(k))
      if(td.lt.t)go to 150
      r6=snr(k)
      r3=r9*(ha(r6,6)-ha(r6,5))
      if(r3.gt.r)go to 150
      if(r3.lt.r)go to 151
      if(v.ge.td)go to 150
  151 v=td
      rrow=j
      rcoll=r6
      r=r3
      rpivot=l1
  150 continue
  160 continue
  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      ha(i,9)=r3
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
c
c  row interchanges.
c
      ha(i,8)=rcoll
  220 if(rrow.eq.i)go to 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t
      do 250 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
      r10=ha(l1,6)
  240 r=r+1
      if(rnr(r).ne.i)go to 240
      rnr(r)=rnr(r10)
  250 rnr(r10)=rrow
      rr3=ha(rrow,2)
      rr4=ha(rrow,3)
      do 270 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
  260 r=r+1
      if(rnr(r).ne.rrow)go to 260
  270 rnr(r)=i
      do 280 j=1,3
      r3=ha(rrow,j)
      ha(rrow,j)=ha(i,j)
c
c  column interchanges.
c
  280 ha(i,j)=r3
  290 if(rcoll.eq.i)go to 350
      do 310 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
      r10=ha(l1,3)
  300 r=r+1
      if(snr(r).ne.i)go to 300
      t=a(r10)
      a(r10)=a(r)
      a(r)=t
      snr(r)=snr(r10)
  310 snr(r10)=rcoll
      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)
      do 330 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
  320 r=r+1
      if(snr(r).ne.rcoll)go to 320
  330 snr(r)=i
      do 340 j=4,6
      r3=ha(rcoll,j)
      ha(rcoll,j)=ha(i,j)
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  340 ha(i,j)=r3
  350 r9=rr4-rr3
      do 360 rr=rr3,rr4
      if(snr(rr).eq.i)go to 370
  360 continue
      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=dabs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431
      do 430 j=rr3,rr4
      index=snr(j)
  430 pivot(index)=a(j)
  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0)go to 669
      r=rr1
      if(r.gt.rr2)go to 470
      do 460 l=r,rr2
      l1=snr(l)
      td=pivot(l1)
      if(td.eq.0.0d0)go to 450
      pivot(l1)=0.0 d0
      td=a(l)-td*t
      a(l)=td
      td1=dabs(td)
      if(td1.gt.aflag(7))aflag(7)=td1
c
c  too small element is created.remove it from the lists.
c
      if(td1.gt.aflag(2))go to 450
      z=z-1
      a(l)=a(rr1)
      snr(l)=snr(rr1)
      a(rr1)=a(i1)
      snr(rr1)=snr(i1)
      snr(i1)=0
      rr1=rr1+1
      i1=i1+1
      ha(r1,2)=rr1
      ha(r1,1)=i1
      r3=ha(l1,5)
      r2=r3-1
      l4=ha(l1,4)
      l5=rnr(l4)
      l6=rnr(r3)
  440 r2=r2+1
      if(rnr(r2).ne.r1)go to 440
      rnr(r2)=l6
      rnr(r3)=l5
      rnr(l4)=0
      ha(l1,5)=r3+1
      ha(l1,4)=l4+1
  450 continue
  460 continue
  470 continue
      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0d0)go to 740
      tol3=-tol2*t
      tol1=dabs(tol3)
      if(tol1.lt.aflag(2))go to 740
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2
      if(iflag(4).ne.1)go to 480
      if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
      if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
  480 if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
c
c  collection in row ordered list.
c
  500 r10=nn-lfr
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1
      do 520 jj=1,n
      l1=ha(jj,3)
      if(l1.lt.ha(jj,1))go to 510
      ha(jj,3)=snr(l1)
      snr(l1)=-jj
  510 continue
  520 continue
      l3=0
      l4=1
      do 550 jj=1,r4
      if(snr(jj).eq.0)go to 540
      l3=l3+1
      if(snr(jj).gt.0)go to 530
      l5=-snr(jj)
      snr(jj)=ha(l5,3)
      ha(l5,3)=l3
      l6=l4+ha(l5,2)-ha(l5,1)
      ha(l5,2)=l6
      ha(l5,1)=l4
      l4=l3+1
  530 a(l3)=a(jj)
      snr(l3)=snr(jj)
  540 continue
  550 continue
      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j
      if(r10.ge.r4)go to 560
      ifail=5
c
c fill-in takes place in the row ordered list.
c
      go to 1110
  560 r8=lfr-1
      rr2=r4+lfr
      if(r8.le.0)go to 579
      l3=i1-1
      do 570 ll=1,r8
      l4=r4+ll
      l5=l3+ll
      a(l4)=a(l5)
      snr(l4)=snr(l5)
  570 snr(l5)=0
  579 rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590
  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=dabs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z) iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
c
c  collection in column ordered list.
c
  630 r10=nn1-lfc
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1
      do 640 jj=i,n
      l1=ha(jj,6)
      ha(jj,6)=rnr(l1)
  640 rnr(l1)=-jj
      l3=0
      l4=1
      do 670 jj=1,r5
      if(rnr(jj).eq.0)go to 660
      l3=l3+1
      if(rnr(jj).gt.0)go to 650
      l5=-rnr(jj)
      rnr(jj)=ha(l5,6)
      ha(l5,6)=l3
      l6=l4+ha(l5,5)-ha(l5,4)
      ha(l5,5)=l6
      ha(l5,4)=l4
      l4=l3+1
  650 rnr(l3)=rnr(jj)
  660 continue
  670 continue
      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)
      if(r10.ge.r5)go to 680
      ifail=6
c
c fill-in takes place in the column ordered list.
c
      go to 1110
  680 r8=lfc-1
      cr2=r5+lfc
      if(r8.le.0)go to 699
      l3=c2-1
      do 690 l=1,r8
      l4=r5+l
      l5=l3+l
      rnr(l4)=rnr(l5)
  690 rnr(l5)=0
  699 cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710
  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
      go to 1110
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue
      if(r9.le.0)go to 882
      do 881 j=rr3,rr4
      index=snr(j)
  881 pivot(index)=0.0 d0
  882 continue
      cr3=ha(i,4)
      do 890 j=cr3,cr4
  890 rnr(j)=0
      if(r9.le.0)go to 930
      l2=ha(i,2)-1
      do 920 ll=1,r9
      r=snr(l2+ll)
      r1=ha(r,5)
      r2=ha(r,6)
      if(r2.gt.r1)go to 900
      ifail=8
      go to 1110
  900 ha(r,5)=r1+1
      r3=r1-1
  910 r3=r3+1
      if(rnr(r3).ne.i)go to 910
      rnr(r3)=rnr(r1)
  920 rnr(r1)=i
  930 aflag(5)=aflag(7)/aflag(6)
      if(aflag(5).lt.aflag(3))go to 940
      ifail=4
      go to 1110
  940 continue
c
c  preparation to begin the back substitution.
c
  950 continue
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0 d0
      td=dabs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td
      if(td.gt.grmin)go to 960
      ifail=3
      go to 1110
  960 if(iflag(4).ne.1)go to 1060
      iflag(10)=ha(n,9)
      iflag(9)=ha(n,10)
      do 990 i=1,n7
      r1=n-i
      iflag(10)=iflag(10)+ha(r1,9)
      iflag(9)=iflag(9)+ha(r1,10)
      if(iflag(3).eq.0)go to 980
      do 970 j=9,10
      r2=ha(r1,j-2)
      r6=ha(r2,j)
      ha(r2,j)=ha(r1,j)
  970 ha(r1,j)=r6
  980 continue
  990 continue
1060  continue
      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end


      subroutine y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
      implicit double precision (a-b,g,p,t-y),integer*4 (c,f,h-n,r-s,z)
      double precision a(nn), pivot(n), b(n)
      integer*4 snr(nn)
      integer*4 iflag(10),ha(iha,11)
      ifail=0
      if(iflag(1).eq.-2)go to 1000
      ifail=1
      go to 1110
1000  mode=iflag(4)
      ipiv=iflag(3)
      n8=n+1
      n7=n-1
      state=iflag(5)
c
c  solve the system with lower triangular matrix  l  (if the
c  lu-factorization is available).
c
      if(state.ne.3)go to 1051
      if(ipiv.eq.0)go to 1020
      do 1010 i=1,n7
      l1=ha(i,7)
      t=b(l1)
      b(l1)=b(i)
      b(i)=t
1010  continue
1020  continue
      do 1050 i=1,n
      rr1=ha(i,1)
      rr2=ha(i,2)-1
      if(rr1.gt.rr2)go to 1040
      do 1030 j=rr1,rr2
      l1=snr(j)
 1030 b(i)=b(i)-a(j)*b(l1)
 1040 continue
 1050 continue
c
c  solve the system with upper triagular matrix.
c
 1051 continue
      do 1090 i=1,n
      r1=n8-i
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      if(rr2.lt.rr1)   go to 1080
      do 1070 j=rr1,rr2
      r2=snr(j)
 1070 b(r1)=b(r1)-a(j)*b(r2)
 1080 continue
 1090 b(r1)=b(r1)/pivot(r1)
c
c if interchanges were used during the  elimination then a reordering in
c lution vector is made.
c
      if(ipiv.eq.0)go to 1110
      do 1100 i=1,n7
      r1=n-i
      r2=ha(r1,8)
      t=b(r2)
      b(r2)=b(r1)
	
1100   b(r1)=t
1110  return
      end
