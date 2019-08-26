      ! DECK DIR
      SUBROUTINE DIR (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,&
         ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK,&
         IWORK)
      ! ***BEGIN PROLOGUE  DIR
      ! ***PURPOSE  Preconditioned Iterative Refinement Sparse Ax = b Solver.
      !             Routine to solve a general linear system  Ax = b  using
      !             iterative refinement with a matrix splitting.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SIR-S, DIR-D)
      ! ***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      !      INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)
      !      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), DZ(N)
      !      DOUBLE PRECISION RWORK(USER DEFINED)
      !      EXTERNAL MATVEC, MSOLVE
      !
      !      CALL DIR(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, ITOL,
      !     $     TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK, IWORK)
      !
      !  *Arguments:
      !  N      :IN       Integer.
      !          Order of the Matrix.
      !  B      :IN       Double Precision B(N).
      !          Right-hand side vector.
      !  X      :INOUT    Double Precision X(N).
      !          On input X is your initial guess for solution vector.
      !          On output X is the final approximate solution.
      !  NELT   :IN       Integer.
      !          Number of Non-Zeros stored in A.
      !  IA     :IN       Integer IA(NELT).
      !  JA     :IN       Integer JA(NELT).
      !  A      :IN       Double Precision A(NELT).
      !          These arrays contain the matrix data structure for A.
      !          It could take any form.  See "Description", below,
      !          for more details.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the upper
      !          or lower triangle of the matrix is stored.
      !  MATVEC :EXT      External.
      !          Name of a routine which performs the matrix vector multiply
      !          Y = A*X given A and X.  The name of the MATVEC routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MATVEC is:
      !              CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
      !          Where N is the number of unknowns, Y is the product A*X
      !          upon return, X is an input vector, NELT is the number of
      !          non-zeros in the SLAP IA, JA, A storage for the matrix A.
      !          ISYM is a flag which, if non-zero, denotes that A is
      !          symmetric and only the lower or upper triangle is stored.
      !  MSOLVE :EXT      External.
      !          Name of a routine which solves a linear system MZ = R for
      !          Z given R with the preconditioning matrix M (M is supplied via
      !          RWORK and IWORK arrays).  The name of the MSOLVE routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MSOLVE is:
      !              CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      !          Where N is the number of unknowns, R is the right-hand side
      !          vector and Z is the solution upon return.  NELT, IA, JA, A and
      !          ISYM are defined as above.  RWORK is a double precision array
      !          that can be used to pass necessary preconditioning information
      !          and/or workspace to MSOLVE.  IWORK is an integer work array
      !          for the same purpose as RWORK.
      !  ITOL   :IN       Integer.
      !          Flag to indicate type of convergence criterion.
      !          If ITOL=1, iteration stops when the 2-norm of the residual
      !          divided by the 2-norm of the right-hand side is less than TOL.
      !          If ITOL=2, iteration stops when the 2-norm of M-inv times the
      !          residual divided by the 2-norm of M-inv times the right hand
      !          side is less than TOL, where M-inv is the inverse of the
      !          diagonal of A.
      !          ITOL=11 is often useful for checking and comparing different
      !          routines.  For this case, the user must supply the "exact"
      !          solution or a very accurate approximation (one with an error
      !          much less than TOL) through a common block,
      !              COMMON /DSLBLK/ SOLN( )
      !          If ITOL=11, iteration stops when the 2-norm of the difference
      !          between the iterative approximation and the user-supplied
      !          solution divided by the 2-norm of the user-supplied solution
      !          is less than TOL.  Note that this requires the user to set up
      !          the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
      !          The routine with this declaration should be loaded before the
      !          stop test so that the correct length is used by the loader.
      !          This procedure is not standard Fortran and may not work
      !          correctly on your system (although it has worked on every
      !          system the authors have tried).  If ITOL is not 11 then this
      !          common block is indeed standard Fortran.
      !  TOL    :INOUT    Double Precision.
      !          Convergence criterion, as described above.  (Reset if IERR=4.)
      !  ITMAX  :IN       Integer.
      !          Maximum number of iterations.
      !  ITER   :OUT      Integer.
      !          Number of iterations required to reach convergence, or
      !          ITMAX+1 if convergence criterion could not be achieved in
      !          ITMAX iterations.
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in final approximate solution, as
      !          defined by ITOL.
      !  IERR   :OUT      Integer.
      !          Return error flag.
      !            IERR = 0 => All went well.
      !            IERR = 1 => Insufficient space allocated for WORK or IWORK.
      !            IERR = 2 => Method failed to converge in ITMAX steps.
      !            IERR = 3 => Error in user input.
      !                        Check input values of N, ITOL.
      !            IERR = 4 => User error tolerance set too tight.
      !                        Reset to 500*D1MACH(3).  Iteration proceeded.
      !            IERR = 5 => Preconditioning matrix, M, is not positive
      !                        definite.  (r,z) < 0.
      !            IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.
      !  IUNIT  :IN       Integer.
      !          Unit number on which to write the error at each iteration,
      !          if this is desired for monitoring convergence.  If unit
      !          number is 0, no writing will occur.
      !  R      :WORK     Double Precision R(N).
      !  Z      :WORK     Double Precision Z(N).
      !  DZ     :WORK     Double Precision DZ(N).
      !          Double Precision arrays used for workspace.
      !  RWORK  :WORK     Double Precision RWORK(USER DEFINED).
      !          Double Precision array that can be used by  MSOLVE.
      !  IWORK  :WORK     Integer IWORK(USER DEFINED).
      !          Integer array that can be used by  MSOLVE.
      !
      !  *Description:
      !        The basic algorithm for iterative refinement (also known as
      !        iterative improvement) is:
      !
      !             n+1    n    -1       n
      !            X    = X  + M  (B - AX  ).
      !
      !            -1   -1
      !        If M =  A then this  is the  standard  iterative  refinement
      !        algorithm and the "subtraction" in the  residual calculation
      !        should be done in double precision (which it is  not in this
      !        routine).
      !        If M = DIAG(A), the diagonal of A, then iterative refinement
      !        is  known  as  Jacobi's  method.   The  SLAP  routine  DSJAC
      !        implements this iterative strategy.
      !        If M = L, the lower triangle of A, then iterative refinement
      !        is known as Gauss-Seidel.   The SLAP routine DSGS implements
      !        this iterative strategy.
      !
      !        This routine does  not care  what matrix data   structure is
      !        used for  A and M.  It simply   calls  the MATVEC and MSOLVE
      !        routines, with  the arguments as  described above.  The user
      !        could write any type of structure and the appropriate MATVEC
      !        and MSOLVE routines.  It is assumed  that A is stored in the
      !        IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
      !        stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP
      !        routines DSJAC and DSGS are examples of this procedure.
      !
      !        Two  examples  of  matrix  data structures  are the: 1) SLAP
      !        Triad  format and 2) SLAP Column format.
      !
      !        =================== S L A P Triad format ===================
      !
      !        In  this   format only the  non-zeros are  stored.  They may
      !        appear  in *ANY* order.   The user  supplies three arrays of
      !        length NELT, where  NELT  is the number  of non-zeros in the
      !        matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
      !        the  user puts   the row  and  column index   of that matrix
      !        element in the IA and JA arrays.  The  value of the non-zero
      !        matrix  element is  placed in  the corresponding location of
      !        the A  array.  This is  an extremely easy data  structure to
      !        generate.  On  the other hand it  is  not too  efficient  on
      !        vector  computers   for the  iterative  solution  of  linear
      !        systems.  Hence, SLAP  changes this input  data structure to
      !        the SLAP   Column  format for the  iteration (but   does not
      !        change it back).
      !
      !        Here is an example of the  SLAP Triad   storage format for a
      !        5x5 Matrix.  Recall that the entries may appear in any order.
      !
      !            5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
      !                               1  2  3  4  5  6  7  8  9 10 11
      !        |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
      !        |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
      !        | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
      !        | 0  0  0 44  0|
      !        |51  0 53  0 55|
      !
      !        =================== S L A P Column format ==================
      !
      !        In  this format   the non-zeros are    stored counting  down
      !        columns (except  for the diagonal  entry, which must  appear
      !        first  in each "column") and are  stored in the  double pre-
      !        cision array  A. In  other  words,  for each  column  in the
      !        matrix  first put  the diagonal entry in A.  Then put in the
      !        other non-zero  elements going  down the column  (except the
      !        diagonal)  in order.  The IA array  holds the  row index for
      !        each non-zero.  The JA array  holds the offsets into the IA,
      !        A  arrays  for  the  beginning  of  each  column.  That  is,
      !        IA(JA(ICOL)),A(JA(ICOL)) are the first elements of the ICOL-
      !        th column in IA and A, and IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)
      !        are  the last elements of the ICOL-th column.   Note that we
      !        always have JA(N+1)=NELT+1, where N is the number of columns
      !        in the matrix  and NELT  is the number  of non-zeros  in the
      !        matrix.
      !
      !        Here is an example of the  SLAP Column  storage format for a
      !        5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
      !        column):
      !
      !            5x5 Matrix      SLAP Column format for 5x5 matrix on left.
      !                               1  2  3    4  5    6  7    8    9 10 11
      !        |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
      !        |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
      !        | 0  0 33  0 35|  JA:  1  4  6    8  9   12
      !        | 0  0  0 44  0|
      !        |51  0 53  0 55|
      !
      !  *Examples:
      !        See the SLAP routines DSJAC, DSGS
      !
      !  *Cautions:
      !      This routine will attempt to write to the Fortran logical output
      !      unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
      !      this logical unit is attached to a file or terminal before calling
      !      this routine with a non-zero value for IUNIT.  This routine does
      !      not check for the validity of a non-zero IUNIT unit number.
      !
      ! ***SEE ALSO  DSJAC, DSGS
      ! ***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
      !                   Johns Hopkins University Press, Baltimore, Maryland,
      !                   1983.
      !                2. Mark K. Seager, A SLAP for the Masses, in
      !                   G. F. Carey, Ed., Parallel Supercomputing: Methods,
      !                   Algorithms and Applications, Wiley, 1989, pp.135-155.
      ! ***ROUTINES CALLED  D1MACH, ISDIR
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890921  Removed TeX from comments.  (FNF)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    891004  Added new reference.
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
      !    920407  COMMON BLOCK renamed DSLBLK.  (WRB)
      !    920511  Added complete declaration section.  (WRB)
      !    920929  Corrected format of references.  (FNF)
      !    921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
      ! ***END PROLOGUE  DIR
      !      .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), DZ(N), R(N), RWORK(*), X(N), Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
      !      .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
      !      .. Local Scalars ..
      DOUBLE PRECISION BNRM, SOLNRM, TOLMIN
      INTEGER I, K
      !      .. External Functions ..
      DOUBLE PRECISION D1MACH
      INTEGER ISDIR
      EXTERNAL D1MACH, ISDIR
      ! ***FIRST EXECUTABLE STATEMENT  DIR
      !
      !          Check some of the input data.
      !
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500*D1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
      !
      !          Calculate initial residual and pseudo-residual, and check
      !          stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I) = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      !
      IF( ISDIR(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,&
           ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK,&
           IWORK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
      IF( IERR.NE.0 ) RETURN
      !
      !          ***** iteration loop *****
      !
      DO 100 K=1,ITMAX
         ITER = K
         !
         !          Calculate new iterate x, new residual r, and new
         !          pseudo-residual z.
         DO 20 I = 1, N
            X(I) = X(I) + Z(I)
 20      CONTINUE
         CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
         DO 30 I = 1, N
            R(I) = B(I) - R(I)
 30      CONTINUE
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         !
         !          check stopping criterion.
         IF( ISDIR(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,&
              ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK,&
              IWORK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
 !
 100  CONTINUE
      !
      !          *****   end of loop  *****
      !          Stopping criterion not satisfied.
      ITER = ITMAX + 1
      IERR = 2
 !
 200  RETURN
      ! ------------- LAST LINE OF DIR FOLLOWS -------------------------------
      END
      ! DECK ISDIR
      INTEGER FUNCTION ISDIR (N, B, X, NELT, IA, JA, A, ISYM, MSOLVE,&
         ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK,&
         IWORK, BNRM, SOLNRM)
      ! ***BEGIN PROLOGUE  ISDIR
      ! ***SUBSIDIARY
      ! ***PURPOSE  Preconditioned Iterative Refinement Stop Test.
      !             This routine calculates the stop test for the iterative
      !             refinement iteration scheme.  It returns a non-zero if the
      !             error estimate (the type of which is determined by ITOL)
      !             is less than the user specified tolerance TOL.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (ISSIR-S, ISDIR-D)
      ! ***KEYWORDS  LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      !      INTEGER IERR, IUNIT, IWORK(USER DEFINED)
      !      DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N), DZ(N)
      !      DOUBLE PRECISION RWORK(USER DEFINED), BNRM, SOLNRM
      !      EXTERNAL MSOLVE
      !
      !      IF( ISDIR(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
      !     $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK, IWORK,
      !     $     BNRM, SOLNRM) .NE. 0 ) THEN ITERATION DONE
      !
      !  *Arguments:
      !  N      :IN       Integer.
      !          Order of the Matrix.
      !  B      :IN       Double Precision B(N).
      !          Right-hand side vector.
      !  X      :IN       Double Precision X(N).
      !          The current approximate solution vector.
      !  NELT   :IN       Integer.
      !          Number of Non-Zeros stored in A.
      !  IA     :IN       Integer IA(NELT).
      !  JA     :IN       Integer JA(NELT).
      !  A      :IN       Double Precision A(NELT).
      !          These arrays contain the matrix data structure for A.
      !          It could take any form.  See "C *Description" in the
      !          DIR routine.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the upper
      !          or lower triangle of the matrix is stored.
      !  MSOLVE :EXT      External.
      !          Name of a routine which solves a linear system Mz = r for
      !          z given r with the preconditioning matrix M (M is supplied via
      !          RWORK and IWORK arrays.  The name of the MSOLVE routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MSOLVE is:
      !              CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      !          Where N is the number of unknowns, R is the right-hand side
      !          vector and Z is the solution upon return.  NELT, IA, JA, A and
      !          ISYM are defined as above.  RWORK is a double precision array
      !          that can be used to pass necessary preconditioning information
      !          and/or workspace to MSOLVE.  IWORK is an integer work array
      !          for the same purpose as RWORK.
      !  ITOL   :IN       Integer.
      !          Flag to indicate type of convergence criterion.
      !          If ITOL=1, iteration stops when the 2-norm of the residual
      !          divided by the 2-norm of the right-hand side is less than TOL.
      !          If ITOL=2, iteration stops when the 2-norm of M-inv times the
      !          residual divided by the 2-norm of M-inv times the right hand
      !          side is less than TOL, where M-inv is the inverse of the
      !          diagonal of A.
      !          ITOL=11 is often useful for checking and comparing different
      !          routines.  For this case, the user must supply the "exact"
      !          solution or a very accurate approximation (one with an error
      !          much less than TOL) through a common block,
      !              COMMON /DSLBLK/ SOLN( )
      !          If ITOL=11, iteration stops when the 2-norm of the difference
      !          between the iterative approximation and the user-supplied
      !          solution divided by the 2-norm of the user-supplied solution
      !          is less than TOL.  Note that this requires the user to set up
      !          the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
      !          The routine with this declaration should be loaded before the
      !          stop test so that the correct length is used by the loader.
      !          This procedure is not standard Fortran and may not work
      !          correctly on your system (although it has worked on every
      !          system the authors have tried).  If ITOL is not 11 then this
      !          common block is indeed standard Fortran.
      !  TOL    :IN       Double Precision.
      !          Convergence criterion, as described above.
      !  ITMAX  :IN       Integer.
      !          Maximum number of iterations.
      !  ITER   :IN       Integer.
      !          Current iteration count.  (Must be zero on first call.)
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in the X(N) approximate solution, as
      !          defined by ITOL.
      !  IERR   :OUT      Integer.
      !          Error flag.  IERR is set to 3 if ITOL is not one of the
      !          acceptable values, see above.
      !  IUNIT  :IN       Integer.
      !          Unit number on which to write the error at each iteration,
      !          if this is desired for monitoring convergence.  If unit
      !          number is 0, no writing will occur.
      !  R      :IN       Double Precision R(N).
      !          The residual R = B-AX.
      !  Z      :WORK     Double Precision Z(N).
      !          Workspace used to hold the pseudo-residual M z = r.
      !  DZ     :WORK     Double Precision DZ(N).
      !          Workspace used to hold temporary vector(s).
      !  RWORK  :WORK     Double Precision RWORK(USER DEFINED).
      !          Double Precision array that can be used by  MSOLVE.
      !  IWORK  :WORK     Integer IWORK(USER DEFINED).
      !          Integer array that can be used by MSOLVE.
      !  BNRM   :INOUT    Double Precision.
      !          Norm of the right hand side.  Type of norm depends on ITOL.
      !          Calculated only on the first call.
      !  SOLNRM :INOUT    Double Precision.
      !          2-Norm of the true solution, SOLN.  Only computed and used
      !          if ITOL = 11.
      !
      !  *Function Return Values:
      !        0 : Error estimate (determined by ITOL) is *NOT* less than the
      !            specified tolerance, TOL.  The iteration must continue.
      !        1 : Error estimate (determined by ITOL) is less than the
      !            specified tolerance, TOL.  The iteration can be considered
      !            complete.
      !
      !  *Cautions:
      !      This routine will attempt to write to the Fortran logical output
      !      unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
      !      this logical unit is attached to a file or terminal before calling
      !      this routine with a non-zero value for IUNIT.  This routine does
      !      not check for the validity of a non-zero IUNIT unit number.
      !
      ! ***SEE ALSO  DIR, DSJAC, DSGS
      ! ***ROUTINES CALLED  D1MACH, DNRM2
      ! ***COMMON BLOCKS    DSLBLK
      ! ***REVISION HISTORY  (YYMMDD)
      !    871119  DATE WRITTEN
      !    880320  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    891003  Removed C***REFER TO line, per MKS.
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
      !    910506  Made subsidiary to DIR.  (FNF)
      !    920407  COMMON BLOCK renamed DSLBLK.  (WRB)
      !    920511  Added complete declaration section.  (WRB)
      !    921026  Changed 1.0E10 to D1MACH(2) and corrected E to D in
      !            output format.  (FNF)
      ! ***END PROLOGUE  ISDIR
      !      .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, SOLNRM, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), DZ(N), R(N), RWORK(*), X(N), Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
      !      .. Subroutine Arguments ..
      EXTERNAL MSOLVE
      !      .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
      !      .. Local Scalars ..
      INTEGER I
      !      .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
      !      .. Common blocks ..
      COMMON /DSLBLK/ SOLN
      ! ***FIRST EXECUTABLE STATEMENT  ISDIR
      ISDIR = 0
      IF( ITOL.EQ.1 ) THEN
         !          err = ||Residual||/||RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)
         ERR = DNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
         !                   -1              -1
         !          err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM2(N, DZ, 1)
         ENDIF
         ERR = DNRM2(N, Z, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
         !          err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF( ITER.EQ.0 ) SOLNRM = DNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM2(N, DZ, 1)/SOLNRM
      ELSE
         !
         !          If we get here ITOL is not one of the acceptable values.
         ERR = D1MACH(2)
         IERR = 3
      ENDIF
      !
      IF( IUNIT.NE.0 ) THEN
         WRITE(IUNIT,1000) ITER,ERR
      ENDIF
      !
      IF( ERR.LE.TOL ) ISDIR = 1
      !
      RETURN
 1000 FORMAT(5X,'ITER = ',I4,' Error Estimate = ',D16.7)
      ! ------------- LAST LINE OF ISDIR FOLLOWS -----------------------------
      END
