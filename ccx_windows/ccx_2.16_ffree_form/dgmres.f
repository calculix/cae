      !
      !    SLATEC: public domain
      !
      ! DECK DGMRES
      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,&
         ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,&
         IGWK, LIGW, RWORK, IWORK)
      ! ***BEGIN PROLOGUE  DGMRES
      ! ***PURPOSE  Preconditioned GMRES iterative sparse Ax=b solver.
      !             This routine uses the generalized minimum residual
      !             (GMRES) method with preconditioning to solve
      !             non-symmetric linear systems of the form: Ax = b.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SGMRES-S, DGMRES-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !
      !  *Usage:
      !       INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      !       INTEGER   ITER, IERR, IUNIT, LRGW, IGWK(LIGW), LIGW
      !       INTEGER   IWORK(USER DEFINED)
      !       DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)
      !       DOUBLE PRECISION RGWK(LRGW), RWORK(USER DEFINED)
      !       EXTERNAL  MATVEC, MSOLVE
      !
      !       CALL DGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
      !      $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX,
      !      $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)
      !
      !  *Arguments:
      !  N      :IN       Integer.
      !          Order of the Matrix.
      !  B      :IN       Double Precision B(N).
      !          Right-hand side vector.
      !  X      :INOUT    Double Precision X(N).
      !          On input X is your initial guess for the solution vector.
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
      !              CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
      !          where N is the number of unknowns, Y is the product A*X
      !          upon return, X is an input vector, and NELT is the number of
      !          non-zeros in the SLAP IA, JA, A storage for the matrix A.
      !          ISYM is a flag which, if non-zero, denotes that A is
      !          symmetric and only the lower or upper triangle is stored.
      !  MSOLVE :EXT      External.
      !          Name of the routine which solves a linear system Mz = r for
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
      !          Flag to indicate the type of convergence criterion used.
      !          ITOL=0  Means the  iteration stops when the test described
      !                  below on  the  residual RL  is satisfied.  This is
      !                  the  "Natural Stopping Criteria" for this routine.
      !                  Other values  of   ITOL  cause  extra,   otherwise
      !                  unnecessary, computation per iteration and     are
      !                  therefore  much less  efficient.  See  ISDGMR (the
      !                  stop test routine) for more information.
      !          ITOL=1  Means   the  iteration stops   when the first test
      !                  described below on  the residual RL  is satisfied,
      !                  and there  is either right  or  no preconditioning
      !                  being used.
      !          ITOL=2  Implies     that   the  user    is   using    left
      !                  preconditioning, and the second stopping criterion
      !                  below is used.
      !          ITOL=3  Means the  iteration stops   when  the  third test
      !                  described below on Minv*Residual is satisfied, and
      !                  there is either left  or no  preconditioning being
      !                  used.
      !          ITOL=11 is    often  useful  for   checking  and comparing
      !                  different routines.  For this case, the  user must
      !                  supply  the  "exact" solution or  a  very accurate
      !                  approximation (one with  an  error much less  than
      !                  TOL) through a common block,
      !                      COMMON /DSLBLK/ SOLN( )
      !                  If ITOL=11, iteration stops when the 2-norm of the
      !                  difference between the iterative approximation and
      !                  the user-supplied solution  divided by the  2-norm
      !                  of the  user-supplied solution  is  less than TOL.
      !                  Note that this requires  the  user to  set up  the
      !                  "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
      !                  routine.  The routine with this declaration should
      !                  be loaded before the stop test so that the correct
      !                  length is used by  the loader.  This procedure  is
      !                  not standard Fortran and may not work correctly on
      !                  your   system (although  it  has  worked  on every
      !                  system the authors have tried).  If ITOL is not 11
      !                  then this common block is indeed standard Fortran.
      !  TOL    :INOUT    Double Precision.
      !          Convergence criterion, as described below.  If TOL is set
      !          to zero on input, then a default value of 500*(the smallest
      !          positive magnitude, machine epsilon) is used.
      !  ITMAX  :DUMMY    Integer.
      !          Maximum number of iterations in most SLAP routines.  In
      !          this routine this does not make sense.  The maximum number
      !          of iterations here is given by ITMAX = MAXL*(NRMAX+1).
      !          See IGWK for definitions of MAXL and NRMAX.
      !  ITER   :OUT      Integer.
      !          Number of iterations required to reach convergence, or
      !          ITMAX if convergence criterion could not be achieved in
      !          ITMAX iterations.
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in final approximate solution, as
      !          defined by ITOL.  Letting norm() denote the Euclidean
      !          norm, ERR is defined as follows..
      !
      !          If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
      !                                for right or no preconditioning, and
      !                          ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
      !                                 norm(SB*(M-inverse)*B),
      !                                for left preconditioning.
      !          If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
      !                                since right or no preconditioning
      !                                being used.
      !          If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
      !                                 norm(SB*(M-inverse)*B),
      !                                since left preconditioning is being
      !                                used.
      !          If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
      !                                i=1,n
      !          If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
      !  IERR   :OUT      Integer.
      !          Return error flag.
      !                IERR = 0 => All went well.
      !                IERR = 1 => Insufficient storage allocated for
      !                            RGWK or IGWK.
      !                IERR = 2 => Routine DGMRES failed to reduce the norm
      !                            of the current residual on its last call,
      !                            and so the iteration has stalled.  In
      !                            this case, X equals the last computed
      !                            approximation.  The user must either
      !                            increase MAXL, or choose a different
      !                            initial guess.
      !                IERR =-1 => Insufficient length for RGWK array.
      !                            IGWK(6) contains the required minimum
      !                            length of the RGWK array.
      !                IERR =-2 => Illegal value of ITOL, or ITOL and JPRE
      !                            values are inconsistent.
      !          For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
      !          left-hand-side of the relevant stopping test defined
      !          below associated with the residual for the current
      !          approximation X(L).
      !  IUNIT  :IN       Integer.
      !          Unit number on which to write the error at each iteration,
      !          if this is desired for monitoring convergence.  If unit
      !          number is 0, no writing will occur.
      !  SB     :IN       Double Precision SB(N).
      !          Array of length N containing scale factors for the right
      !          hand side vector B.  If JSCAL.eq.0 (see below), SB need
      !          not be supplied.
      !  SX     :IN       Double Precision SX(N).
      !          Array of length N containing scale factors for the solution
      !          vector X.  If JSCAL.eq.0 (see below), SX need not be
      !          supplied.  SB and SX can be the same array in the calling
      !          program if desired.
      !  RGWK   :INOUT    Double Precision RGWK(LRGW).
      !          Double Precision array used for workspace by DGMRES.
      !          On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.
      !  LRGW   :IN       Integer.
      !          Length of the double precision workspace, RGWK.
      !          LRGW >= 1 + N*(MAXL+6) + MAXL*(MAXL+3).
      !          See below for definition of MAXL.
      !          For the default values, RGWK has size at least 131 + 16*N.
      !  IGWK   :INOUT    Integer IGWK(LIGW).
      !          The following IGWK parameters should be set by the user
      !          before calling this routine.
      !          IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in
      !             which X - X0 is to be found (where, X0 is the initial
      !             guess).  The default value of MAXL is 10.
      !          IGWK(2) = KMP.  Maximum number of previous Krylov basis
      !             vectors to which each new basis vector is made orthogonal.
      !             The default value of KMP is MAXL.
      !          IGWK(3) = JSCAL.  Flag indicating whether the scaling
      !             arrays SB and SX are to be used.
      !             JSCAL = 0 => SB and SX are not used and the algorithm
      !                will perform as if all SB(I) = 1 and SX(I) = 1.
      !             JSCAL = 1 =>  Only SX is used, and the algorithm
      !                performs as if all SB(I) = 1.
      !             JSCAL = 2 =>  Only SB is used, and the algorithm
      !                performs as if all SX(I) = 1.
      !             JSCAL = 3 =>  Both SB and SX are used.
      !          IGWK(4) = JPRE.  Flag indicating whether preconditioning
      !             is being used.
      !             JPRE = 0  =>  There is no preconditioning.
      !             JPRE > 0  =>  There is preconditioning on the right
      !                only, and the solver will call routine MSOLVE.
      !             JPRE < 0  =>  There is preconditioning on the left
      !                only, and the solver will call routine MSOLVE.
      !          IGWK(5) = NRMAX.  Maximum number of restarts of the
      !             Krylov iteration.  The default value of NRMAX = 10.
      !             if IWORK(5) = -1,  then no restarts are performed (in
      !             this case, NRMAX is set to zero internally).
      !          The following IWORK parameters are diagnostic information
      !          made available to the user after this routine completes.
      !          IGWK(6) = MLWK.  Required minimum length of RGWK array.
      !          IGWK(7) = NMS.  The total number of calls to MSOLVE.
      !  LIGW   :IN       Integer.
      !          Length of the integer workspace, IGWK.  LIGW >= 20.
      !  RWORK  :WORK     Double Precision RWORK(USER DEFINED).
      !          Double Precision array that can be used for workspace in
      !          MSOLVE.
      !  IWORK  :WORK     Integer IWORK(USER DEFINED).
      !          Integer array that can be used for workspace in MSOLVE.
      !
      !  *Description:
      !        DGMRES solves a linear system A*X = B rewritten in the form:
      !
      !         (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
      !
      !        with right preconditioning, or
      !
      !         (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
      !
      !        with left preconditioning, where A is an N-by-N double precision
      !        matrix, X and B are N-vectors,  SB and SX  are diagonal scaling
      !        matrices,   and M is  a preconditioning    matrix.   It uses
      !        preconditioned  Krylov   subpace  methods  based     on  the
      !        generalized minimum residual  method (GMRES).   This routine
      !        optionally performs  either  the  full     orthogonalization
      !        version of the  GMRES  algorithm or an incomplete variant of
      !        it.  Both versions use restarting of the linear iteration by
      !        default, although the user can disable this feature.
      !
      !        The GMRES  algorithm generates a sequence  of approximations
      !        X(L) to the  true solution of the above  linear system.  The
      !        convergence criteria for stopping the  iteration is based on
      !        the size  of the  scaled norm of  the residual  R(L)  =  B -
      !        A*X(L).  The actual stopping test is either:
      !
      !                norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
      !
      !        for right preconditioning, or
      !
      !                norm(SB*(M-inverse)*(B-A*X(L))) .le.
      !                        TOL*norm(SB*(M-inverse)*B),
      !
      !        for left preconditioning, where norm() denotes the Euclidean
      !        norm, and TOL is  a positive scalar less  than one  input by
      !        the user.  If TOL equals zero  when DGMRES is called, then a
      !        default  value  of 500*(the   smallest  positive  magnitude,
      !        machine epsilon) is used.  If the  scaling arrays SB  and SX
      !        are used, then  ideally they  should be chosen  so  that the
      !        vectors SX*X(or SX*M*X) and  SB*B have all their  components
      !        approximately equal  to  one in  magnitude.  If one wants to
      !        use the same scaling in X  and B, then  SB and SX can be the
      !        same array in the calling program.
      !
      !        The following is a list of the other routines and their
      !        functions used by DGMRES:
      !        DPIGMR  Contains the main iteration loop for GMRES.
      !        DORTH   Orthogonalizes a new vector against older basis vectors.
      !        DHEQR   Computes a QR decomposition of a Hessenberg matrix.
      !        DHELS   Solves a Hessenberg least-squares system, using QR
      !                factors.
      !        DRLCAL  Computes the scaled residual RL.
      !        DXLCAL  Computes the solution XL.
      !        ISDGMR  User-replaceable stopping routine.
      !
      !        This routine does  not care  what matrix data   structure is
      !        used for  A and M.  It simply   calls  the MATVEC and MSOLVE
      !        routines, with  the arguments as  described above.  The user
      !        could write any type of structure and the appropriate MATVEC
      !        and MSOLVE routines.  It is assumed  that A is stored in the
      !        IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
      !        stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
      !        routines DSDCG and DSICCG are examples of this procedure.
      !
      !        Two  examples  of  matrix  data structures  are the: 1) SLAP
      !        Triad  format and 2) SLAP Column format.
      !
      !        =================== S L A P Triad format ===================
      !        This routine requires that the  matrix A be   stored in  the
      !        SLAP  Triad format.  In  this format only the non-zeros  are
      !        stored.  They may appear in  *ANY* order.  The user supplies
      !        three arrays of  length NELT, where  NELT is  the number  of
      !        non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
      !        each non-zero the user puts the row and column index of that
      !        matrix element  in the IA and  JA arrays.  The  value of the
      !        non-zero   matrix  element is  placed  in  the corresponding
      !        location of the A array.   This is  an  extremely  easy data
      !        structure to generate.  On  the  other hand it   is  not too
      !        efficient on vector computers for  the iterative solution of
      !        linear systems.  Hence,   SLAP changes   this  input    data
      !        structure to the SLAP Column format  for  the iteration (but
      !        does not change it back).
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
      !        This routine  requires that  the matrix A  be stored in  the
      !        SLAP Column format.  In this format the non-zeros are stored
      !        counting down columns (except for  the diagonal entry, which
      !        must appear first in each  "column")  and are stored  in the
      !        double precision array A.   In other words,  for each column
      !        in the matrix put the diagonal entry in  A.  Then put in the
      !        other non-zero  elements going down  the column (except  the
      !        diagonal) in order.   The  IA array holds the  row index for
      !        each non-zero.  The JA array holds the offsets  into the IA,
      !        A arrays  for  the  beginning  of each   column.   That  is,
      !        IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
      !        ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
      !        A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
      !        Note that we always have  JA(N+1) = NELT+1,  where N is  the
      !        number of columns in  the matrix and NELT  is the number  of
      !        non-zeros in the matrix.
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
      !  *Cautions:
      !      This routine will attempt to write to the Fortran logical output
      !      unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
      !      this logical unit is attached to a file or terminal before calling
      !      this routine with a non-zero value for IUNIT.  This routine does
      !      not check for the validity of a non-zero IUNIT unit number.
      !
      ! ***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage
      !                   Matrix Methods in Stiff ODE Systems, Lawrence Liver-
      !                   more National Laboratory Report UCRL-95088, Rev. 1,
      !                   Livermore, California, June 1987.
      !                2. Mark K. Seager, A SLAP for the Masses, in
      !                   G. F. Carey, Ed., Parallel Supercomputing: Methods,
      !                   Algorithms and Applications, Wiley, 1989, pp.135-155.
      ! ***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DPIGMR
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    891004  Added new reference.
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910506  Corrected errors in C***ROUTINES CALLED list.  (FNF)
      !    920407  COMMON BLOCK renamed DSLBLK.  (WRB)
      !    920511  Added complete declaration section.  (WRB)
      !    920929  Corrected format of references.  (FNF)
      !    921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
      !    921026  Added check for valid value of ITOL.  (FNF)
      ! ***END PROLOGUE  DGMRES
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LIGW, LRGW, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RGWK(LRGW), RWORK(*), SB(N),&
                       SX(N), X(N)
      INTEGER IA(NELT), IGWK(LIGW), IWORK(*), JA(NELT)
      !      .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
      !      .. Local Scalars ..
      DOUBLE PRECISION BNRM, RHOL, SUM
      INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV,&
              LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS
      !      .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
      !      .. External Subroutines ..
      EXTERNAL DCOPY, DPIGMR
      !      .. Intrinsic Functions ..
      INTRINSIC SQRT
      ! ***FIRST EXECUTABLE STATEMENT  DGMRES
      IERR = 0
      !    ------------------------------------------------------------------
      !          Load method parameters with user values or defaults.
      !    ------------------------------------------------------------------
      MAXL = IGWK(1)
      IF (MAXL .EQ. 0) MAXL = 10
      IF (MAXL .GT. N) MAXL = N
      KMP = IGWK(2)
      IF (KMP .EQ. 0) KMP = MAXL
      IF (KMP .GT. MAXL) KMP = MAXL
      JSCAL = IGWK(3)
      JPRE = IGWK(4)
      !          Check for valid value of ITOL.
      IF( (ITOL.LT.0) .OR. ((ITOL.GT.3).AND.(ITOL.NE.11)) ) GOTO 650
      !          Check for consistent values of ITOL and JPRE.
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650
      NRMAX = IGWK(5)
      IF( NRMAX.EQ.0 ) NRMAX = 10
      !          If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
      IF( NRMAX.EQ.-1 ) NRMAX = 0
      !          If input value of TOL is zero, set it to its default value.
      IF( TOL.EQ.0.0D0 ) TOL = 500*D1MACH(3)
      !
      !          Initialize counters.
      ITER = 0
      NMS = 0
      NRSTS = 0
      !    ------------------------------------------------------------------
      !          Form work array segment pointers.
      !    ------------------------------------------------------------------
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + N*MAXLP1
      LHES = LR + N + 1
      LQ = LHES + MAXL*MAXLP1
      LDL = LQ + 2*MAXL
      LW = LDL + N
      LXL = LW + N
      LZ = LXL + N
      !
      !          Load IGWK(6) with required minimum length of the RGWK array.
      IGWK(6) = LZ + N - 1
      IF( LZ+N-1.GT.LRGW ) GOTO 640
      !    ------------------------------------------------------------------
      !          Calculate scaled-preconditioned norm of RHS vector b.
      !    ------------------------------------------------------------------
      IF (JPRE .LT. 0) THEN
         CALL MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM,&
              RWORK, IWORK)
         NMS = NMS + 1
      ELSE
         CALL DCOPY(N, B, 1, RGWK(LR), 1)
      ENDIF
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN
         SUM = 0
         DO 10 I = 1,N
            SUM = SUM + (RGWK(LR-1+I)*SB(I))**2
 10      CONTINUE
         BNRM = SQRT(SUM)
      ELSE
         BNRM = DNRM2(N,RGWK(LR),1)
      ENDIF
      !    ------------------------------------------------------------------
      !          Calculate initial residual.
      !    ------------------------------------------------------------------
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
      DO 50 I = 1,N
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
 !    ------------------------------------------------------------------
 !          If performing restarting, then load the residual into the
 !          correct location in the RGWK array.
 !    ------------------------------------------------------------------
 100  CONTINUE
      IF( NRSTS.GT.NRMAX ) GOTO 610
      IF( NRSTS.GT.0 ) THEN
         !          Copy the current residual to a different location in the RGWK
         !          array.
         CALL DCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)
      ENDIF
      !    ------------------------------------------------------------------
      !          Use the DPIGMR algorithm to solve the linear system A*Z = R.
      !    ------------------------------------------------------------------
      CALL DPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,&
             NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),&
             RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),&
             RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,&
             TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      ITER = ITER + LGMR
      NMS = NMS + NMSL
      !
      !          Increment X by the current approximate solution Z of A*Z = R.
      !
      LZM1 = LZ - 1
      DO 110 I = 1,N
         X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
      IF( IFLAG.EQ.0 ) GOTO 600
      IF( IFLAG.EQ.1 ) THEN
         NRSTS = NRSTS + 1
         GOTO 100
      ENDIF
      IF( IFLAG.EQ.2 ) GOTO 620
 !    ------------------------------------------------------------------
 !          All returns are made through this section.
 !    ------------------------------------------------------------------
 !          The iteration has converged.
 !
 600  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 0
      RETURN
 !
 !          Max number((NRMAX+1)*MAXL) of linear iterations performed.
 610  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 1
      RETURN
 !
 !          GMRES failed to reduce last residual in MAXL iterations.
 !          The iteration has stalled.
 620  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 2
      RETURN
 !          Error return.  Insufficient length for RGWK array.
 640  CONTINUE
      ERR = TOL
      IERR = -1
      RETURN
 !          Error return.  Inconsistent ITOL and JPRE values.
 650  CONTINUE
      ERR = TOL
      IERR = -2
      RETURN
      ! ------------- LAST LINE OF DGMRES FOLLOWS ----------------------------
      END
      ! DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)
      ! ***BEGIN PROLOGUE  DNRM2
      ! ***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
      ! ***LIBRARY   SLATEC (BLAS)
      ! ***CATEGORY  D1A3B
      ! ***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
      ! ***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
      !              LINEAR ALGEBRA, UNITARY, VECTOR
      ! ***AUTHOR  Lawson, C. L., (JPL)
      !            Hanson, R. J., (SNLA)
      !            Kincaid, D. R., (U. of Texas)
      !            Krogh, F. T., (JPL)
      ! ***DESCRIPTION
      !
      !                 B L A S  Subprogram
      !     Description of parameters
      !
      !      --Input--
      !         N  number of elements in input vector(s)
      !        DX  double precision vector with N elements
      !      INCX  storage spacing between elements of DX
      !
      !      --Output--
      !     DNRM2  double precision result (zero if N .LE. 0)
      !
      !      Euclidean norm of the N-vector stored in DX with storage
      !      increment INCX.
      !      If N .LE. 0, return with result = 0.
      !      If N .GE. 1, then INCX must be .GE. 1
      !
      !      Four phase method using two built-in constants that are
      !      hopefully applicable to all machines.
      !          CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
      !          CUTHI = minimum of  SQRT(V)      over all known machines.
      !      where
      !          EPS = smallest no. such that EPS + 1. .GT. 1.
      !          U   = smallest positive no.   (underflow limit)
      !          V   = largest  no.            (overflow  limit)
      !
      !      Brief outline of algorithm.
      !
      !      Phase 1 scans zero components.
      !      move to phase 2 when a component is nonzero and .LE. CUTLO
      !      move to phase 3 when a component is .GT. CUTLO
      !      move to phase 4 when a component is .GE. CUTHI/M
      !      where M = N for X() real and M = 2*N for complex.
      !
      !      Values for CUTLO and CUTHI.
      !      From the environmental parameters listed in the IMSL converter
      !      document the limiting values are as follows:
      !      CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
      !                    Univac and DEC at 2**(-103)
      !                    Thus CUTLO = 2**(-51) = 4.44089E-16
      !      CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
      !                    Thus CUTHI = 2**(63.5) = 1.30438E19
      !      CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
      !                    Thus CUTLO = 2**(-33.5) = 8.23181D-11
      !      CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
      !      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
      !      DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
      !
      ! ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
      !                  Krogh, Basic linear algebra subprograms for Fortran
      !                  usage, Algorithm No. 539, Transactions on Mathematical
      !                  Software 5, 3 (September 1979), pp. 308-323.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    791001  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DNRM2
      INTEGER NEXT
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,&
                       ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
      !
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
      ! ***FIRST EXECUTABLE STATEMENT  DNRM2
      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
   !
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
      !
      !                                                  BEGIN MAIN LOOP
      !
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
   !
   !                         PHASE 1.  SUM IS ZERO
   !
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      !
      !                                 PREPARE FOR PHASE 2.
      !
      ASSIGN 70 TO NEXT
      GO TO 105
  !
  !                                 PREPARE FOR PHASE 4.
  !
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
   !
   !                    PHASE 2.  SUM IS SMALL.
   !                              SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
   !
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
  !
  !                      COMMON CODE FOR PHASES 2 AND 4.
  !                      IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
  !
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
  !
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
   !
   !                   PREPARE FOR PHASE 3.
   !
   75 SUM = (SUM * XMAX) * XMAX
   !
   !      FOR REAL OR D.P. SET HITEST = CUTHI/N
   !      FOR COMPLEX      SET HITEST = CUTHI/(2*N)
   !
   85 HITEST = CUTHI / N
      !
      !                    PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
      !
      DO 95 J = I,NN,INCX
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = SQRT(SUM)
      GO TO 300
  !
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
      !
      !               END OF MAIN LOOP.
      !
      !               COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
      !
      DNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
      ! DECK DPIGMR
      SUBROUTINE DPIGMR (N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS,&
         JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK,&
         DL, RHOL, NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A,&
         ISYM, IUNIT, IFLAG, ERR)
      ! ***BEGIN PROLOGUE  DPIGMR
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !          This routine solves the linear system A * Z = R0 using a
      !          scaled preconditioned version of the generalized minimum
      !          residual method.  An initial guess of Z = 0 is assumed.
      !
      !  *Usage:
      !       INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR
      !       INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT)
      !       INTEGER ISYM, IUNIT, IFLAG
      !       DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1),
      !      $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED),
      !      $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N),
      !      $                 TOL, A(NELT), ERR
      !       EXTERNAL MATVEC, MSOLVE
      !
      !       CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP,
      !      $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
      !      $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
      !      $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          The order of the matrix A, and the lengths
      !          of the vectors SR, SZ, R0 and Z.
      !  R0     :IN       Double Precision R0(N)
      !          R0 = the right hand side of the system A*Z = R0.
      !          R0 is also used as workspace when computing
      !          the final approximation.
      !          (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.)
      !  SR     :IN       Double Precision SR(N)
      !          SR is a vector of length N containing the non-zero
      !          elements of the diagonal scaling matrix for R0.
      !  SZ     :IN       Double Precision SZ(N)
      !          SZ is a vector of length N containing the non-zero
      !          elements of the diagonal scaling matrix for Z.
      !  JSCAL  :IN       Integer
      !          A flag indicating whether arrays SR and SZ are used.
      !          JSCAL=0 means SR and SZ are not used and the
      !                  algorithm will perform as if all
      !                  SR(i) = 1 and SZ(i) = 1.
      !          JSCAL=1 means only SZ is used, and the algorithm
      !                  performs as if all SR(i) = 1.
      !          JSCAL=2 means only SR is used, and the algorithm
      !                  performs as if all SZ(i) = 1.
      !          JSCAL=3 means both SR and SZ are used.
      !  MAXL   :IN       Integer
      !          The maximum allowable order of the matrix H.
      !  MAXLP1 :IN       Integer
      !          MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
      !  KMP    :IN       Integer
      !          The number of previous vectors the new vector VNEW
      !          must be made orthogonal to.  (KMP .le. MAXL)
      !  NRSTS  :IN       Integer
      !          Counter for the number of restarts on the current
      !          call to DGMRES.  If NRSTS .gt. 0, then the residual
      !          R0 is already scaled, and so scaling of it is
      !          not necessary.
      !  JPRE   :IN       Integer
      !          Preconditioner type flag.
      !  MATVEC :EXT      External.
      !          Name of a routine which performs the matrix vector multiply
      !          Y = A*X given A and X.  The name of the MATVEC routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MATVEC is:
      !              CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
      !          where N is the number of unknowns, Y is the product A*X
      !          upon return, X is an input vector, and NELT is the number of
      !          non-zeros in the SLAP IA, JA, A storage for the matrix A.
      !          ISYM is a flag which, if non-zero, denotes that A is
      !          symmetric and only the lower or upper triangle is stored.
      !  MSOLVE :EXT      External.
      !          Name of the routine which solves a linear system Mz = r for
      !          z given r with the preconditioning matrix M (M is supplied via
      !          RPAR and IPAR arrays.  The name of the MSOLVE routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MSOLVE is:
      !              CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
      !          Where N is the number of unknowns, R is the right-hand side
      !          vector and Z is the solution upon return.  NELT, IA, JA, A and
      !          ISYM are defined as below.  RPAR is a double precision array
      !          that can be used to pass necessary preconditioning information
      !          and/or workspace to MSOLVE.  IPAR is an integer work array
      !          for the same purpose as RPAR.
      !  NMSL   :OUT      Integer
      !          The number of calls to MSOLVE.
      !  Z      :OUT      Double Precision Z(N)
      !          The final computed approximation to the solution
      !          of the system A*Z = R0.
      !  V      :OUT      Double Precision V(N,MAXLP1)
      !          The N by (LGMR+1) array containing the LGMR
      !          orthogonal vectors V(*,1) to V(*,LGMR).
      !  HES    :OUT      Double Precision HES(MAXLP1,MAXL)
      !          The upper triangular factor of the QR decomposition
      !          of the (LGMR+1) by LGMR upper Hessenberg matrix whose
      !          entries are the scaled inner-products of A*V(*,I)
      !          and V(*,K).
      !  Q      :OUT      Double Precision Q(2*MAXL)
      !          A double precision array of length 2*MAXL containing the
      !          components of the Givens rotations used in the QR
      !          decomposition of HES.  It is loaded in DHEQR and used in
      !          DHELS.
      !  LGMR   :OUT      Integer
      !          The number of iterations performed and
      !          the current order of the upper Hessenberg
      !          matrix HES.
      !  RPAR   :IN       Double Precision RPAR(USER DEFINED)
      !          Double Precision workspace passed directly to the MSOLVE
      !          routine.
      !  IPAR   :IN       Integer IPAR(USER DEFINED)
      !          Integer workspace passed directly to the MSOLVE routine.
      !  WK     :IN       Double Precision WK(N)
      !          A double precision work array of length N used by routines
      !          MATVEC and MSOLVE.
      !  DL     :INOUT    Double Precision DL(N)
      !          On input, a double precision work array of length N used for
      !          calculation of the residual norm RHO when the method is
      !          incomplete (KMP.lt.MAXL), and/or when using restarting.
      !          On output, the scaled residual vector RL.  It is only loaded
      !          when performing restarts of the Krylov iteration.
      !  RHOL   :OUT      Double Precision
      !          A double precision scalar containing the norm of the final
      !          residual.
      !  NRMAX  :IN       Integer
      !          The maximum number of restarts of the Krylov iteration.
      !          NRMAX .gt. 0 means restarting is active, while
      !          NRMAX = 0 means restarting is not being used.
      !  B      :IN       Double Precision B(N)
      !          The right hand side of the linear system A*X = b.
      !  BNRM   :IN       Double Precision
      !          The scaled norm of b.
      !  X      :IN       Double Precision X(N)
      !          The current approximate solution as of the last
      !          restart.
      !  XL     :IN       Double Precision XL(N)
      !          An array of length N used to hold the approximate
      !          solution X(L) when ITOL=11.
      !  ITOL   :IN       Integer
      !          A flag to indicate the type of convergence criterion
      !          used.  See the driver for its description.
      !  TOL    :IN       Double Precision
      !          The tolerance on residuals R0-A*Z in scaled norm.
      !  NELT   :IN       Integer
      !          The length of arrays IA, JA and A.
      !  IA     :IN       Integer IA(NELT)
      !          An integer array of length NELT containing matrix data.
      !          It is passed directly to the MATVEC and MSOLVE routines.
      !  JA     :IN       Integer JA(NELT)
      !          An integer array of length NELT containing matrix data.
      !          It is passed directly to the MATVEC and MSOLVE routines.
      !  A      :IN       Double Precision A(NELT)
      !          A double precision array of length NELT containing matrix
      !          data. It is passed directly to the MATVEC and MSOLVE routines.
      !  ISYM   :IN       Integer
      !          A flag to indicate symmetric matrix storage.
      !          If ISYM=0, all non-zero entries of the matrix are
      !          stored.  If ISYM=1, the matrix is symmetric and
      !          only the upper or lower triangular part is stored.
      !  IUNIT  :IN       Integer
      !          The i/o unit number for writing intermediate residual
      !          norm values.
      !  IFLAG  :OUT      Integer
      !          An integer error flag..
      !          0 means convergence in LGMR iterations, LGMR.le.MAXL.
      !          1 means the convergence test did not pass in MAXL
      !            iterations, but the residual norm is .lt. norm(R0),
      !            and so Z is computed.
      !          2 means the convergence test did not pass in MAXL
      !            iterations, residual .ge. norm(R0), and Z = 0.
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in final approximate solution, as
      !          defined by ITOL.
      !
      !  *Cautions:
      !      This routine will attempt to write to the Fortran logical output
      !      unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
      !      this logical unit is attached to a file or terminal before calling
      !      this routine with a non-zero value for IUNIT.  This routine does
      !      not check for the validity of a non-zero IUNIT unit number.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL,
      !                     DSCAL, ISDGMR
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DPIGMR
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, RHOL, TOL
      INTEGER IFLAG, ISYM, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, MAXL,&
              MAXLP1, N, NELT, NMSL, NRMAX, NRSTS
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(*), DL(*), HES(MAXLP1,*), Q(*), R0(*),&
                       RPAR(*), SR(*), SZ(*), V(N,*), WK(*), X(*),&
                       XL(*), Z(*)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
      !      .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
      !      .. Local Scalars ..
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1
      !      .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISDGMR
      EXTERNAL DNRM2, ISDGMR
      !      .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS, DHEQR, DORTH, DRLCAL, DSCAL
      !      .. Intrinsic Functions ..
      INTRINSIC ABS
      ! ***FIRST EXECUTABLE STATEMENT  DPIGMR
      !
      !          Zero out the Z array.
      !
      DO 5 I = 1,N
         Z(I) = 0
 5    CONTINUE
      !
      IFLAG = 0
      LGMR = 0
      NMSL = 0
      !          Load ITMAX, the maximum number of iterations.
      ITMAX =(NRMAX+1)*MAXL
      !    -------------------------------------------------------------------
      !          The initial residual is the vector R0.
      !          Apply left precon. if JPRE < 0 and this is not a restart.
      !          Apply scaling to R0 if JSCAL = 2 or 3.
      !    -------------------------------------------------------------------
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
         CALL DCOPY(N, R0, 1, WK, 1)
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      R0NRM = DNRM2(N, V, 1)
      ITER = NRSTS*MAXL
      !
      !          Call stopping routine ISDGMR.
      !
      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,&
          NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK,&
          RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL,&
          KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,&
          HES, JPRE) .NE. 0) RETURN
      TEM = 1.0D0/R0NRM
      CALL DSCAL(N, TEM, V(1,1), 1)
      !
      !          Zero out the HES array.
      !
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0
 40      CONTINUE
 50   CONTINUE
      !    -------------------------------------------------------------------
      !          Main loop to compute the vectors V(*,2) to V(*,MAXL).
      !          The running product PROD is needed for the convergence test.
      !    -------------------------------------------------------------------
      PROD = 1
      DO 90 LL = 1,MAXL
         LGMR = LL
        !    -------------------------------------------------------------------
        !         Unscale  the  current V(LL)  and store  in WK.  Call routine
        !         MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the
        !         preconditioner matrix.  Save the answer in Z.   Call routine
        !         MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
        !         matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
        !         routine DORTH  to  orthogonalize the    new vector VNEW   =
        !         V(*,LL+1).  Call routine DHEQR to update the factors of HES.
        !    -------------------------------------------------------------------
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
           DO 60 I = 1,N
              WK(I) = V(I,LL)/SZ(I)
 60        CONTINUE
        ELSE
           CALL DCOPY(N, V(1,LL), 1, WK, 1)
        ENDIF
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ELSE
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ENDIF
        IF (JPRE .LT. 0) THEN
           CALL DCOPY(N, V(1,LL+1), 1, WK, 1)
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
        !    -------------------------------------------------------------------
        !          Update RHO, the estimate of the norm of the residual R0-A*ZL.
        !          If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
        !          necessarily orthogonal for LL > KMP.  The vector DL must then
        !          be computed, and its norm used in the calculation of RHO.
        !    -------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*R0NRM)
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL DCOPY(N, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
           DLNRM = DNRM2(N, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
        RHOL = RHO
        !    -------------------------------------------------------------------
        !          Test for convergence.  If passed, compute approximation ZL.
        !          If failed and LL < MAXL, then continue iterating.
        !    -------------------------------------------------------------------
        ITER = NRSTS*MAXL + LGMR
        IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,&
            NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK,&
            RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL,&
            KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,&
            HES, JPRE) .NE. 0) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
        !    -------------------------------------------------------------------
        !          Rescale so that the norm of V(1,LL+1) is one.
        !    -------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL(N, TEM, V(1,LL+1), 1)
 90   CONTINUE
 100  CONTINUE
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
      !
      !          Load approximate solution with zero.
      !
      DO 130 I = 1,N
         Z(I) = 0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
      !
      !          Tolerance not met, but residual norm reduced.
      !
      IF (NRMAX .GT. 0) THEN
         !
         !         If performing restarting (NRMAX > 0)  calculate the residual
         !         vector RL and  store it in the DL  array.  If the incomplete
         !         version is being used (KMP < MAXL) then DL has  already been
         !         calculated up to a scaling factor.   Use DRLCAL to calculate
         !         the scaled residual vector.
         !
         CALL DRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD,&
              R0NRM)
      ENDIF
 !    -------------------------------------------------------------------
 !          Compute the approximation ZL to the solution.  Since the
 !          vector Z was used as workspace, and the initial guess
 !          of the linear iteration is zero, Z must be reset to zero.
 !    -------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0
 210  CONTINUE
      R0(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0
 220  CONTINUE
      DO 230 I = 1,LL
         CALL DAXPY(N, R0(I), V(1,I), 1, Z, 1)
 230  CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, Z, 1, WK, 1)
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      RETURN
      ! ------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------
      END
      ! DECK DRLCAL
      SUBROUTINE DRLCAL (N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,&
         R0NRM)
      ! ***BEGIN PROLOGUE  DRLCAL
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SRLCAL-S, DRLCAL-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !          This routine calculates the scaled residual RL from the
      !          V(I)'s.
      !  *Usage:
      !       INTEGER N, KMP, LL, MAXL
      !       DOUBLE PRECISION V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM
      !
      !       CALL DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM)
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          The order of the matrix A, and the lengths
      !          of the vectors SR, SZ, R0 and Z.
      !  KMP    :IN       Integer
      !          The number of previous V vectors the new vector VNEW
      !          must be made orthogonal to. (KMP .le. MAXL)
      !  LL     :IN       Integer
      !          The current dimension of the Krylov subspace.
      !  MAXL   :IN       Integer
      !          The maximum dimension of the Krylov subspace.
      !  V      :IN       Double Precision V(N,LL)
      !          The N x LL array containing the orthogonal vectors
      !          V(*,1) to V(*,LL).
      !  Q      :IN       Double Precision Q(2*MAXL)
      !          A double precision array of length 2*MAXL containing the
      !          components of the Givens rotations used in the QR
      !          decomposition of HES.  It is loaded in DHEQR and used in
      !          DHELS.
      !  RL     :OUT      Double Precision RL(N)
      !          The residual vector RL.  This is either SB*(B-A*XL) if
      !          not preconditioning or preconditioning on the right,
      !          or SB*(M-inverse)*(B-A*XL) if preconditioning on the
      !          left.
      !  SNORMW :IN       Double Precision
      !          Scale factor.
      !  PROD   :IN       Double Precision
      !          The product s1*s2*...*sl = the product of the sines of the
      !          Givens rotations used in the QR factorization of
      !          the Hessenberg matrix HES.
      !  R0NRM  :IN       Double Precision
      !          The scaled norm of initial residual R0.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  DCOPY, DSCAL
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DRLCAL
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      DOUBLE PRECISION PROD, R0NRM, SNORMW
      INTEGER KMP, LL, MAXL, N
      !      .. Array Arguments ..
      DOUBLE PRECISION Q(*), RL(N), V(N,*)
      !      .. Local Scalars ..
      DOUBLE PRECISION C, S, TEM
      INTEGER I, I2, IP1, K, LLM1, LLP1
      !      .. External Subroutines ..
      EXTERNAL DCOPY, DSCAL
      ! ***FIRST EXECUTABLE STATEMENT  DRLCAL
      IF (KMP .EQ. MAXL) THEN
         !
         !          calculate RL.  Start by copying V(*,1) into RL.
         !
         CALL DCOPY(N, V(1,1), 1, RL, 1)
         LLM1 = LL - 1
         DO 20 I = 1,LLM1
            IP1 = I + 1
            I2 = I*2
            S = Q(I2)
            C = Q(I2-1)
            DO 10 K = 1,N
               RL(K) = S*RL(K) + C*V(K,IP1)
 10         CONTINUE
 20      CONTINUE
         S = Q(2*LL)
         C = Q(2*LL-1)/SNORMW
         LLP1 = LL + 1
         DO 30 K = 1,N
            RL(K) = S*RL(K) + C*V(K,LLP1)
 30      CONTINUE
      ENDIF
      !
      !          When KMP < MAXL, RL vector already partially calculated.
      !          Scale RL by R0NRM*PROD to obtain the residual RL.
      !
      TEM = R0NRM*PROD
      CALL DSCAL(N, TEM, RL, 1)
      RETURN
      ! ------------- LAST LINE OF DRLCAL FOLLOWS ----------------------------
      END
      ! DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
      use omp_lib
      ! ***BEGIN PROLOGUE  DAXPY
      ! ***PURPOSE  Compute a constant times a vector plus a vector.
      ! ***LIBRARY   SLATEC (BLAS)
      ! ***CATEGORY  D1A7
      ! ***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
      ! ***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
      ! ***AUTHOR  Lawson, C. L., (JPL)
      !            Hanson, R. J., (SNLA)
      !            Kincaid, D. R., (U. of Texas)
      !            Krogh, F. T., (JPL)
      ! ***DESCRIPTION
      !
      !                 B L A S  Subprogram
      !     Description of Parameters
      !
      !      --Input--
      !         N  number of elements in input vector(s)
      !        DA  double precision scalar multiplier
      !        DX  double precision vector with N elements
      !      INCX  storage spacing between elements of DX
      !        DY  double precision vector with N elements
      !      INCY  storage spacing between elements of DY
      !
      !      --Output--
      !        DY  double precision result (unchanged if N .LE. 0)
      !
      !      Overwrite double precision DY with double precision DA*DX + DY.
      !      For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
      !        DY(LY+I*INCY),
      !      where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
      !      defined in a similar way using INCY.
      !
      ! ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
      !                  Krogh, Basic linear algebra subprograms for Fortran
      !                  usage, Algorithm No. 539, Transactions on Mathematical
      !                  Software 5, 3 (September 1979), pp. 308-323.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    791001  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
      ! ***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
    !
    !      Code for unequal or nonpositive increments.
    !
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
   !
   !      Code for both increments equal to 1.
   !
   !      Clean-up loop so remaining vector length is a multiple of 4.
   !
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
   !
   !      Code for equal, positive, non-unit increments.
   !
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
      ! DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
      use omp_lib
      ! ***BEGIN PROLOGUE  DCOPY
      ! ***PURPOSE  Copy a vector.
      ! ***LIBRARY   SLATEC (BLAS)
      ! ***CATEGORY  D1A5
      ! ***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
      ! ***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
      ! ***AUTHOR  Lawson, C. L., (JPL)
      !            Hanson, R. J., (SNLA)
      !            Kincaid, D. R., (U. of Texas)
      !            Krogh, F. T., (JPL)
      ! ***DESCRIPTION
      !
      !                 B L A S  Subprogram
      !     Description of Parameters
      !
      !      --Input--
      !         N  number of elements in input vector(s)
      !        DX  double precision vector with N elements
      !      INCX  storage spacing between elements of DX
      !        DY  double precision vector with N elements
      !      INCY  storage spacing between elements of DY
      !
      !      --Output--
      !        DY  copy of vector DX (unchanged if N .LE. 0)
      !
      !      Copy double precision DX to double precision DY.
      !      For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
      !      where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
      !      defined in a similar way using INCY.
      !
      ! ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
      !                  Krogh, Basic linear algebra subprograms for Fortran
      !                  usage, Algorithm No. 539, Transactions on Mathematical
      !                  Software 5, 3 (September 1979), pp. 308-323.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    791001  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
      ! ***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
    !
    !      Code for unequal or nonpositive increments.
    !
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
   !
   !      Code for both increments equal to 1.
   !
   !      Clean-up loop so remaining vector length is a multiple of 7.
   !
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
   !
   !      Code for equal, positive, non-unit increments.
   !
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
      ! DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
      ! ***BEGIN PROLOGUE  DHELS
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SHELS-S, DHELS-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !         This routine is extracted from the LINPACK routine SGESL with
      !         changes due to the fact that A is an upper Hessenberg matrix.
      !
      !         DHELS solves the least squares problem:
      !
      !                    MIN(B-A*X,B-A*X)
      !
      !         using the factors computed by DHEQR.
      !
      !  *Usage:
      !       INTEGER LDA, N
      !       DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)
      !
      !       CALL DHELS(A, LDA, N, Q, B)
      !
      !  *Arguments:
      !  A       :IN       Double Precision A(LDA,N)
      !           The output from DHEQR which contains the upper
      !           triangular factor R in the QR decomposition of A.
      !  LDA     :IN       Integer
      !           The leading dimension of the array A.
      !  N       :IN       Integer
      !           A is originally an (N+1) by N matrix.
      !  Q       :IN       Double Precision Q(2*N)
      !           The coefficients of the N Givens rotations
      !           used in the QR factorization of A.
      !  B       :INOUT    Double Precision B(N+1)
      !           On input, B is the right hand side vector.
      !           On output, B is the solution vector X.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  DAXPY
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DHELS
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      INTEGER LDA, N
      !      .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
      !      .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
      !      .. External Subroutines ..
      EXTERNAL DAXPY
      ! ***FIRST EXECUTABLE STATEMENT  DHELS
      !
      !          Minimize(B-A*X,B-A*X).  First form Q*B.
      !
      DO 20 K = 1, N
         KP1 = K + 1
         IQ = 2*(K-1) + 1
         C = Q(IQ)
         S = Q(IQ+1)
         T1 = B(K)
         T2 = B(KP1)
         B(K) = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
      !
      !          Now solve  R*X = Q*B.
      !
      DO 40 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
      ! ------------- LAST LINE OF DHELS FOLLOWS ----------------------------
      END
      ! DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
      ! ***BEGIN PROLOGUE  DHEQR
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SHEQR-S, DHEQR-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !         This   routine  performs  a QR   decomposition  of an  upper
      !         Hessenberg matrix A using Givens  rotations.  There  are two
      !         options  available: 1)  Performing  a fresh decomposition 2)
      !         updating the QR factors by adding a row and  a column to the
      !         matrix A.
      !
      !  *Usage:
      !       INTEGER LDA, N, INFO, IJOB
      !       DOUBLE PRECISION A(LDA,N), Q(2*N)
      !
      !       CALL DHEQR(A, LDA, N, Q, INFO, IJOB)
      !
      !  *Arguments:
      !  A      :INOUT    Double Precision A(LDA,N)
      !          On input, the matrix to be decomposed.
      !          On output, the upper triangular matrix R.
      !          The factorization can be written Q*A = R, where
      !          Q is a product of Givens rotations and R is upper
      !          triangular.
      !  LDA    :IN       Integer
      !          The leading dimension of the array A.
      !  N      :IN       Integer
      !          A is an (N+1) by N Hessenberg matrix.
      !  Q      :OUT      Double Precision Q(2*N)
      !          The factors c and s of each Givens rotation used
      !          in decomposing A.
      !  INFO   :OUT      Integer
      !          = 0  normal value.
      !          = K  if  A(K,K) .eq. 0.0 .  This is not an error
      !            condition for this subroutine, but it does
      !            indicate that DHELS will divide by zero
      !            if called.
      !  IJOB   :IN       Integer
      !          = 1     means that a fresh decomposition of the
      !                  matrix A is desired.
      !          .ge. 2  means that the current decomposition of A
      !                  will be updated by the addition of a row
      !                  and a column.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DHEQR
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      INTEGER IJOB, INFO, LDA, N
      !      .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), Q(*)
      !      .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      !      .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
      ! ***FIRST EXECUTABLE STATEMENT  DHEQR
      IF (IJOB .GT. 1) GO TO 70
      !    -------------------------------------------------------------------
      !          A new factorization is desired.
      !    -------------------------------------------------------------------
      !          QR decomposition without pivoting.
      !
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
         !
         !            Compute K-th column of R.
         !            First, multiply the K-th column of A by the previous
         !            K-1 Givens rotations.
         !
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
 !
 !          Compute Givens components C and S.
 !
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0D0 ) THEN
            C = 1
            S = 0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0D0/SQRT(1.0D0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0D0/SQRT(1.0D0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0D0 ) INFO = K
 60   CONTINUE
      RETURN
 !    -------------------------------------------------------------------
 !          The old factorization of a will be updated.  A row and a
 !          column has been added to the matrix A.  N by N-1 is now
 !          the old size of the matrix.
 !    -------------------------------------------------------------------
 70   CONTINUE
      NM1 = N - 1
      !    -------------------------------------------------------------------
      !          Multiply the new column by the N previous Givens rotations.
      !    -------------------------------------------------------------------
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
      !    -------------------------------------------------------------------
      !          Complete update of decomposition by forming last Givens
      !          rotation, and multiplying it times the column
      !          vector(A(N,N),A(NP1,N)).
      !    -------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0D0 ) THEN
         C = 1
         S = 0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0D0/SQRT(1.0D0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0D0/SQRT(1.0D0+T*T)
         S = -C*T
      ENDIF
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      ! ------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
      END
      ! DECK ISDGMR
      INTEGER FUNCTION ISDGMR (N, B, X, XL, NELT, IA, JA, A, ISYM,&
         MSOLVE, NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,&
         RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL,&
         MAXLP1, V, Q, SNORMW, PROD, R0NRM, HES, JPRE)
      ! ***BEGIN PROLOGUE  ISDGMR
      ! ***SUBSIDIARY
      ! ***PURPOSE  Generalized Minimum Residual Stop Test.
      !             This routine calculates the stop test for the Generalized
      !             Minimum RESidual (GMRES) iteration scheme.  It returns a
      !             non-zero if the error estimate (the type of which is
      !             determined by ITOL) is less than the user specified
      !             tolerance TOL.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (ISSGMR-S, ISDGMR-D)
      ! ***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !
      !  *Usage:
      !       INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL
      !       INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL
      !       INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE
      !       DOUBLE PRECISION B(N), X(N), XL(MAXL), A(NELT), TOL, ERR,
      !      $                 R(N), Z(N), DZ(N), RWORK(USER DEFINED),
      !      $                 RNRM, BNRM, SB(N), SX(N), V(N,MAXLP1),
      !      $                 Q(2*MAXL), SNORMW, PROD, R0NRM,
      !      $                 HES(MAXLP1,MAXL)
      !       EXTERNAL MSOLVE
      !
      !       IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
      !      $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
      !      $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
      !      $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
      !      $     HES, JPRE) .NE. 0) THEN ITERATION DONE
      !
      !  *Arguments:
      !  N      :IN       Integer.
      !          Order of the Matrix.
      !  B      :IN       Double Precision B(N).
      !          Right-hand-side vector.
      !  X      :IN       Double Precision X(N).
      !          Approximate solution vector as of the last restart.
      !  XL     :OUT      Double Precision XL(N)
      !          An array of length N used to hold the approximate
      !          solution as of the current iteration.  Only computed by
      !          this routine when ITOL=11.
      !  NELT   :IN       Integer.
      !          Number of Non-Zeros stored in A.
      !  IA     :IN       Integer IA(NELT).
      !  JA     :IN       Integer JA(NELT).
      !  A      :IN       Double Precision A(NELT).
      !          These arrays contain the matrix data structure for A.
      !          It could take any form.  See "Description", in the DGMRES,
      !          DSLUGM and DSDGMR routines for more details.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the upper
      !          or lower triangle of the matrix is stored.
      !  MSOLVE :EXT      External.
      !          Name of a routine which solves a linear system Mz = r for  z
      !          given r with the preconditioning matrix M (M is supplied via
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
      !  NMSL   :INOUT    Integer.
      !          A counter for the number of calls to MSOLVE.
      !  ITOL   :IN       Integer.
      !          Flag to indicate the type of convergence criterion used.
      !          ITOL=0  Means the  iteration stops when the test described
      !                  below on  the  residual RL  is satisfied.  This is
      !                  the  "Natural Stopping Criteria" for this routine.
      !                  Other values  of   ITOL  cause  extra,   otherwise
      !                  unnecessary, computation per iteration and     are
      !                  therefore much less efficient.
      !          ITOL=1  Means   the  iteration stops   when the first test
      !                  described below on  the residual RL  is satisfied,
      !                  and there  is either right  or  no preconditioning
      !                  being used.
      !          ITOL=2  Implies     that   the  user    is   using    left
      !                  preconditioning, and the second stopping criterion
      !                  below is used.
      !          ITOL=3  Means the  iteration stops   when  the  third test
      !                  described below on Minv*Residual is satisfied, and
      !                  there is either left  or no  preconditioning begin
      !                  used.
      !          ITOL=11 is    often  useful  for   checking  and comparing
      !                  different routines.  For this case, the  user must
      !                  supply  the  "exact" solution or  a  very accurate
      !                  approximation (one with  an  error much less  than
      !                  TOL) through a common block,
      !                      COMMON /DSLBLK/ SOLN( )
      !                  If ITOL=11, iteration stops when the 2-norm of the
      !                  difference between the iterative approximation and
      !                  the user-supplied solution  divided by the  2-norm
      !                  of the  user-supplied solution  is  less than TOL.
      !                  Note that this requires  the  user to  set up  the
      !                  "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
      !                  routine.  The routine with this declaration should
      !                  be loaded before the stop test so that the correct
      !                  length is used by  the loader.  This procedure  is
      !                  not standard Fortran and may not work correctly on
      !                  your   system (although  it  has  worked  on every
      !                  system the authors have tried).  If ITOL is not 11
      !                  then this common block is indeed standard Fortran.
      !  TOL    :IN       Double Precision.
      !          Convergence criterion, as described above.
      !  ITMAX  :IN       Integer.
      !          Maximum number of iterations.
      !  ITER   :IN       Integer.
      !          The iteration for which to check for convergence.
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in final approximate solution, as
      !          defined by ITOL.  Letting norm() denote the Euclidean
      !          norm, ERR is defined as follows..
      !
      !          If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
      !                                for right or no preconditioning, and
      !                          ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
      !                                 norm(SB*(M-inverse)*B),
      !                                for left preconditioning.
      !          If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
      !                                since right or no preconditioning
      !                                being used.
      !          If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
      !                                 norm(SB*(M-inverse)*B),
      !                                since left preconditioning is being
      !                                used.
      !          If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
      !                                i=1,n
      !          If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
      !  IUNIT  :IN       Integer.
      !          Unit number on which to write the error at each iteration,
      !          if this is desired for monitoring convergence.  If unit
      !          number is 0, no writing will occur.
      !  R      :INOUT    Double Precision R(N).
      !          Work array used in calling routine.  It contains
      !          information necessary to compute the residual RL = B-A*XL.
      !  Z      :WORK     Double Precision Z(N).
      !          Workspace used to hold the pseudo-residual M z = r.
      !  DZ     :WORK     Double Precision DZ(N).
      !          Workspace used to hold temporary vector(s).
      !  RWORK  :WORK     Double Precision RWORK(USER DEFINED).
      !          Double Precision array that can be used by MSOLVE.
      !  IWORK  :WORK     Integer IWORK(USER DEFINED).
      !          Integer array that can be used by MSOLVE.
      !  RNRM   :IN       Double Precision.
      !          Norm of the current residual.  Type of norm depends on ITOL.
      !  BNRM   :IN       Double Precision.
      !          Norm of the right hand side.  Type of norm depends on ITOL.
      !  SB     :IN       Double Precision SB(N).
      !          Scaling vector for B.
      !  SX     :IN       Double Precision SX(N).
      !          Scaling vector for X.
      !  JSCAL  :IN       Integer.
      !          Flag indicating if scaling arrays SB and SX are being
      !          used in the calling routine DPIGMR.
      !          JSCAL=0 means SB and SX are not used and the
      !                  algorithm will perform as if all
      !                  SB(i) = 1 and SX(i) = 1.
      !          JSCAL=1 means only SX is used, and the algorithm
      !                  performs as if all SB(i) = 1.
      !          JSCAL=2 means only SB is used, and the algorithm
      !                  performs as if all SX(i) = 1.
      !          JSCAL=3 means both SB and SX are used.
      !  KMP    :IN       Integer
      !          The number of previous vectors the new vector VNEW
      !          must be made orthogonal to.  (KMP .le. MAXL)
      !  LGMR   :IN       Integer
      !          The number of GMRES iterations performed on the current call
      !          to DPIGMR (i.e., # iterations since the last restart) and
      !          the current order of the upper Hessenberg
      !          matrix HES.
      !  MAXL   :IN       Integer
      !          The maximum allowable order of the matrix H.
      !  MAXLP1 :IN       Integer
      !          MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
      !  V      :IN       Double Precision V(N,MAXLP1)
      !          The N by (LGMR+1) array containing the LGMR
      !          orthogonal vectors V(*,1) to V(*,LGMR).
      !  Q      :IN       Double Precision Q(2*MAXL)
      !          A double precision array of length 2*MAXL containing the
      !          components of the Givens rotations used in the QR
      !          decomposition of HES.
      !  SNORMW :IN       Double Precision
      !          A scalar containing the scaled norm of VNEW before it
      !          is renormalized in DPIGMR.
      !  PROD   :IN       Double Precision
      !          The product s1*s2*...*sl = the product of the sines of the
      !          Givens rotations used in the QR factorization of the
      !          Hessenberg matrix HES.
      !  R0NRM  :IN       Double Precision
      !          The scaled norm of initial residual R0.
      !  HES    :IN       Double Precision HES(MAXLP1,MAXL)
      !          The upper triangular factor of the QR decomposition
      !          of the (LGMR+1) by LGMR upper Hessenberg matrix whose
      !          entries are the scaled inner-products of A*V(*,I)
      !          and V(*,K).
      !  JPRE   :IN       Integer
      !          Preconditioner type flag.
      !          (See description of IGWK(4) in DGMRES.)
      !
      !  *Description
      !        When using the GMRES solver,  the preferred value  for ITOL
      !        is 0.  This is due to the fact that when ITOL=0 the norm of
      !        the residual required in the stopping test is  obtained for
      !        free, since this value is already  calculated  in the GMRES
      !        algorithm.   The  variable  RNRM contains the   appropriate
      !        norm, which is equal to norm(SB*(RL - A*XL))  when right or
      !        no   preconditioning is  being  performed,   and equal   to
      !        norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
      !        Here, norm() is the Euclidean norm.  Nonzero values of ITOL
      !        require  additional work  to  calculate the  actual  scaled
      !        residual  or its scaled/preconditioned  form,  and/or   the
      !        approximate solution XL.  Hence, these values of  ITOL will
      !        not be as efficient as ITOL=0.
      !
      !  *Cautions:
      !      This routine will attempt to write to the Fortran logical output
      !      unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
      !      this logical unit is attached to a file or terminal before calling
      !      this routine with a non-zero value for IUNIT.  This routine does
      !      not check for the validity of a non-zero IUNIT unit number.
      !
      !      This routine does not verify that ITOL has a valid value.
      !      The calling routine should make such a test before calling
      !      ISDGMR, as is done in DGMRES.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DRLCAL, DSCAL, DXLCAL
      ! ***COMMON BLOCKS    DSLBLK
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Corrected conversion errors, etc.  (FNF)
      !    910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920407  COMMON BLOCK renamed DSLBLK.  (WRB)
      !    920511  Added complete declaration section.  (WRB)
      !    921026  Corrected D to E in output format.  (FNF)
      !    921113  Corrected C***CATEGORY line.  (FNF)
      ! ***END PROLOGUE  ISDGMR
      !      .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, PROD, R0NRM, RNRM, SNORMW, TOL
      INTEGER ISYM, ITER, ITMAX, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR,&
              MAXL, MAXLP1, N, NELT, NMSL
      !      .. Array Arguments ..
      DOUBLE PRECISION A(*), B(*), DZ(*), HES(MAXLP1, MAXL), Q(*), R(*),&
                       RWORK(*), SB(*), SX(*), V(N,*), X(*), XL(*), Z(*)
      INTEGER IA(*), IWORK(*), JA(*)
      !      .. Subroutine Arguments ..
      EXTERNAL MSOLVE
      !      .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
      !      .. Local Scalars ..
      DOUBLE PRECISION DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM
      INTEGER I, IELMAX
      !      .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
      !      .. External Subroutines ..
      EXTERNAL DCOPY, DRLCAL, DSCAL, DXLCAL
      !      .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SQRT
      !      .. Common blocks ..
      COMMON /DSLBLK/ SOLN
      !      .. Save statement ..
      SAVE SOLNRM
      ! ***FIRST EXECUTABLE STATEMENT  ISDGMR
      ISDGMR = 0
      IF ( ITOL.EQ.0 ) THEN
         !
         !        Use input from DPIGMR to determine if stop conditions are met.
         !
         ERR = RNRM/BNRM
      ENDIF
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN
         !
         !        Use DRLCAL to calculate the scaled residual vector.
         !        Store answer in R.
         !
         IF ( LGMR.NE.0 ) CALL DRLCAL(N, KMP, LGMR, MAXL, V, Q, R,&
                                      SNORMW, PROD, R0NRM)
         IF ( ITOL.LE.2 ) THEN
            !          err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = DNRM2(N, R, 1)/BNRM
            !
            !          Unscale R by R0NRM*PROD when KMP < MAXL.
            !
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
         ELSEIF ( ITOL.EQ.3 ) THEN
            !          err = Max |(Minv*Residual)(i)/x(i)|
            !          When JPRE .lt. 0, R already contains Minv*Residual.
            IF ( JPRE.GT.0 ) THEN
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK,&
                    IWORK)
               NMSL = NMSL + 1
            ENDIF
            !
            !          Unscale R by R0NRM*PROD when KMP < MAXL.
            !
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
            !
            FUZZ = D1MACH(1)
            IELMAX = 1
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ)
            DO 25 I = 2, N
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ)
               IF( RAT.GT.RATMAX ) THEN
                  IELMAX = I
                  RATMAX = RAT
               ENDIF
 25         CONTINUE
            ERR = RATMAX
            IF( RATMAX.LE.TOL ) ISDGMR = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
         !
         !        Use DXLCAL to calculate the approximate solution XL.
         !
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL DXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,&
                 DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK,&
                 NELT, IA, JA, A, ISYM)
         ELSEIF ( ITER.EQ.0 ) THEN
            !          Copy X to XL to check if initial guess is good enough.
            CALL DCOPY(N, X, 1, XL, 1)
         ELSE
            !          Return since this is the first call to DPIGMR on a restart.
            RETURN
         ENDIF
         !
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN
            !          err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF ( ITER.EQ.0 ) SOLNRM = DNRM2(N, SOLN, 1)
            DO 30 I = 1, N
               DZ(I) = XL(I) - SOLN(I)
 30         CONTINUE
            ERR = DNRM2(N, DZ, 1)/SOLNRM
         ELSE
            IF (ITER .EQ. 0) THEN
               SOLNRM = 0
               DO 40 I = 1,N
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2
 40            CONTINUE
               SOLNRM = SQRT(SOLNRM)
            ENDIF
            DXNRM = 0
            DO 50 I = 1,N
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2
 50         CONTINUE
            DXNRM = SQRT(DXNRM)
            !          err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM
         ENDIF
      ENDIF
      !
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISDGMR = 1
      !
      RETURN
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ',&
           'N, ITOL = ',I5, I5,&
           /' ITER','   Natural Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7)
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5,&
           ' |R(IELMAX)/X(IELMAX)| = ',D12.5)
      ! ------------- LAST LINE OF ISDGMR FOLLOWS ----------------------------
      END
      ! DECK DORTH
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
      ! ***BEGIN PROLOGUE  DORTH
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SORTH-S, DORTH-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !         This routine  orthogonalizes  the  vector  VNEW  against the
      !         previous KMP  vectors in the   V array.  It uses  a modified
      !         Gram-Schmidt   orthogonalization procedure with  conditional
      !         reorthogonalization.
      !
      !  *Usage:
      !       INTEGER N, LL, LDHES, KMP
      !       DOUBLE PRECISION VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
      !
      !       CALL DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
      !
      !  *Arguments:
      !  VNEW   :INOUT    Double Precision VNEW(N)
      !          On input, the vector of length N containing a scaled
      !          product of the Jacobian and the vector V(*,LL).
      !          On output, the new vector orthogonal to V(*,i0) to V(*,LL),
      !          where i0 = max(1, LL-KMP+1).
      !  V      :IN       Double Precision V(N,LL)
      !          The N x LL array containing the previous LL
      !          orthogonal vectors V(*,1) to V(*,LL).
      !  HES    :INOUT    Double Precision HES(LDHES,LL)
      !          On input, an LL x LL upper Hessenberg matrix containing,
      !          in HES(I,K), K.lt.LL, the scaled inner products of
      !          A*V(*,K) and V(*,i).
      !          On return, column LL of HES is filled in with
      !          the scaled inner products of A*V(*,LL) and V(*,i).
      !  N      :IN       Integer
      !          The order of the matrix A, and the length of VNEW.
      !  LL     :IN       Integer
      !          The current order of the matrix HES.
      !  LDHES  :IN       Integer
      !          The leading dimension of the HES array.
      !  KMP    :IN       Integer
      !          The number of previous vectors the new vector VNEW
      !          must be made orthogonal to (KMP .le. MAXL).
      !  SNORMW :OUT      DOUBLE PRECISION
      !          Scalar containing the l-2 norm of VNEW.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  DAXPY, DDOT, DNRM2
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DORTH
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      DOUBLE PRECISION SNORMW
      INTEGER KMP, LDHES, LL, N
      !      .. Array Arguments ..
      DOUBLE PRECISION HES(LDHES,*), V(N,*), VNEW(*)
      !      .. Local Scalars ..
      DOUBLE PRECISION ARG, SUMDSQ, TEM, VNRM
      INTEGER I, I0
      !      .. External Functions ..
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL DDOT, DNRM2
      !      .. External Subroutines ..
      EXTERNAL DAXPY
      !      .. Intrinsic Functions ..
      INTRINSIC MAX, SQRT
      ! ***FIRST EXECUTABLE STATEMENT  DORTH
      !
      !          Get norm of unaltered VNEW for later use.
      !
      VNRM = DNRM2(N, VNEW, 1)
      !    -------------------------------------------------------------------
      !          Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
      !          Scaled inner products give new column of HES.
      !          Projections of earlier vectors are subtracted from VNEW.
      !    -------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
         HES(I,LL) = DDOT(N, V(1,I), 1, VNEW, 1)
         TEM = -HES(I,LL)
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
      !    -------------------------------------------------------------------
      !          Compute SNORMW = norm of VNEW.  If VNEW is small compared
      !          to its input value (in norm), then reorthogonalize VNEW to
      !          V(*,1) through V(*,LL).  Correct if relative correction
      !          exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
      !          the dot products involved.
      !    -------------------------------------------------------------------
      SNORMW = DNRM2(N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0
      DO 30 I = I0,LL
         TEM = -DDOT(N, V(1,I), 1, VNEW, 1)
         IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
      !
      RETURN
      ! ------------- LAST LINE OF DORTH FOLLOWS ----------------------------
      END
      ! DECK DXLCAL
      SUBROUTINE DXLCAL (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,&
         WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A,&
         ISYM)
      ! ***BEGIN PROLOGUE  DXLCAL
      ! ***SUBSIDIARY
      ! ***PURPOSE  Internal routine for DGMRES.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SXLCAL-S, DXLCAL-D)
      ! ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
      !            Hindmarsh, Alan, (LLNL), alanh@llnl.gov
      !            Seager, Mark K., (LLNL), seager@llnl.gov
      !              Lawrence Livermore National Laboratory
      !              PO Box 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      ! ***DESCRIPTION
      !         This  routine computes the solution  XL,  the current DGMRES
      !         iterate, given the  V(I)'s and  the  QR factorization of the
      !         Hessenberg  matrix HES.   This routine  is  only called when
      !         ITOL=11.
      !
      !  *Usage:
      !       INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED)
      !       INTEGER NELT, IA(NELT), JA(NELT), ISYM
      !       DOUBLE PRECISION X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL),
      !      $                 V(N,MAXLP1), R0NRM, WK(N), SZ(N),
      !      $                 RPAR(USER DEFINED), A(NELT)
      !       EXTERNAL MSOLVE
      !
      !       CALL DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
      !      $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
      !      $     NELT, IA, JA, A, ISYM)
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          The order of the matrix A, and the lengths
      !          of the vectors SR, SZ, R0 and Z.
      !  LGMR   :IN       Integer
      !          The number of iterations performed and
      !          the current order of the upper Hessenberg
      !          matrix HES.
      !  X      :IN       Double Precision X(N)
      !          The current approximate solution as of the last restart.
      !  XL     :OUT      Double Precision XL(N)
      !          An array of length N used to hold the approximate
      !          solution X(L).
      !          Warning: XL and ZL are the same array in the calling routine.
      !  ZL     :IN       Double Precision ZL(N)
      !          An array of length N used to hold the approximate
      !          solution Z(L).
      !  HES    :IN       Double Precision HES(MAXLP1,MAXL)
      !          The upper triangular factor of the QR decomposition
      !          of the (LGMR+1) by LGMR upper Hessenberg matrix whose
      !          entries are the scaled inner-products of A*V(*,i) and V(*,k).
      !  MAXLP1 :IN       Integer
      !          MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
      !          MAXL is the maximum allowable order of the matrix HES.
      !  Q      :IN       Double Precision Q(2*MAXL)
      !          A double precision array of length 2*MAXL containing the
      !          components of the Givens rotations used in the QR
      !          decomposition of HES.  It is loaded in DHEQR.
      !  V      :IN       Double Precision V(N,MAXLP1)
      !          The N by(LGMR+1) array containing the LGMR
      !          orthogonal vectors V(*,1) to V(*,LGMR).
      !  R0NRM  :IN       Double Precision
      !          The scaled norm of the initial residual for the
      !          current call to DPIGMR.
      !  WK     :IN       Double Precision WK(N)
      !          A double precision work array of length N.
      !  SZ     :IN       Double Precision SZ(N)
      !          A vector of length N containing the non-zero
      !          elements of the diagonal scaling matrix for Z.
      !  JSCAL  :IN       Integer
      !          A flag indicating whether arrays SR and SZ are used.
      !          JSCAL=0 means SR and SZ are not used and the
      !                  algorithm will perform as if all
      !                  SR(i) = 1 and SZ(i) = 1.
      !          JSCAL=1 means only SZ is used, and the algorithm
      !                  performs as if all SR(i) = 1.
      !          JSCAL=2 means only SR is used, and the algorithm
      !                  performs as if all SZ(i) = 1.
      !          JSCAL=3 means both SR and SZ are used.
      !  JPRE   :IN       Integer
      !          The preconditioner type flag.
      !  MSOLVE :EXT      External.
      !          Name of the routine which solves a linear system Mz = r for
      !          z given r with the preconditioning matrix M (M is supplied via
      !          RPAR and IPAR arrays.  The name of the MSOLVE routine must
      !          be declared external in the calling program.  The calling
      !          sequence to MSOLVE is:
      !              CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
      !          Where N is the number of unknowns, R is the right-hand side
      !          vector and Z is the solution upon return.  NELT, IA, JA, A and
      !          ISYM are defined as below.  RPAR is a double precision array
      !          that can be used to pass necessary preconditioning information
      !          and/or workspace to MSOLVE.  IPAR is an integer work array
      !          for the same purpose as RPAR.
      !  NMSL   :IN       Integer
      !          The number of calls to MSOLVE.
      !  RPAR   :IN       Double Precision RPAR(USER DEFINED)
      !          Double Precision workspace passed directly to the MSOLVE
      !          routine.
      !  IPAR   :IN       Integer IPAR(USER DEFINED)
      !          Integer workspace passed directly to the MSOLVE routine.
      !  NELT   :IN       Integer
      !          The length of arrays IA, JA and A.
      !  IA     :IN       Integer IA(NELT)
      !          An integer array of length NELT containing matrix data.
      !          It is passed directly to the MATVEC and MSOLVE routines.
      !  JA     :IN       Integer JA(NELT)
      !          An integer array of length NELT containing matrix data.
      !          It is passed directly to the MATVEC and MSOLVE routines.
      !  A      :IN       Double Precision A(NELT)
      !          A double precision array of length NELT containing matrix
      !          data.
      !          It is passed directly to the MATVEC and MSOLVE routines.
      !  ISYM   :IN       Integer
      !          A flag to indicate symmetric matrix storage.
      !          If ISYM=0, all non-zero entries of the matrix are
      !          stored.  If ISYM=1, the matrix is symmetric and
      !          only the upper or lower triangular part is stored.
      !
      ! ***SEE ALSO  DGMRES
      ! ***ROUTINES CALLED  DAXPY, DCOPY, DHELS
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
      !    910506  Made subsidiary to DGMRES.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      ! ***END PROLOGUE  DXLCAL
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Scalar Arguments ..
      DOUBLE PRECISION R0NRM
      INTEGER ISYM, JPRE, JSCAL, LGMR, MAXLP1, N, NELT, NMSL
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), HES(MAXLP1,*), Q(*), RPAR(*), SZ(*),&
                       V(N,*), WK(N), X(N), XL(N), ZL(N)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
      !      .. Subroutine Arguments ..
      EXTERNAL MSOLVE
      !      .. Local Scalars ..
      INTEGER I, K, LL, LLP1
      !      .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS
      ! ***FIRST EXECUTABLE STATEMENT  DXLCAL
      LL = LGMR
      LLP1 = LL + 1
      DO 10 K = 1,LLP1
         WK(K) = 0
 10   CONTINUE
      WK(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, WK)
      DO 20 K = 1,N
         ZL(K) = 0
 20   CONTINUE
      DO 30 I = 1,LL
         CALL DAXPY(N, WK(I), V(1,I), 1, ZL, 1)
 30   CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 40 K = 1,N
            ZL(K) = ZL(K)/SZ(K)
 40      CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, ZL, 1, WK, 1)
         CALL MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      !          calculate XL from X and ZL.
      DO 50 K = 1,N
         XL(K) = X(K) + ZL(K)
 50   CONTINUE
      RETURN
      ! ------------- LAST LINE OF DXLCAL FOLLOWS ----------------------------
      END
      ! DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
      use omp_lib
      ! ***BEGIN PROLOGUE  DDOT
      ! ***PURPOSE  Compute the inner product of two vectors.
      ! ***LIBRARY   SLATEC (BLAS)
      ! ***CATEGORY  D1A4
      ! ***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
      ! ***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
      ! ***AUTHOR  Lawson, C. L., (JPL)
      !            Hanson, R. J., (SNLA)
      !            Kincaid, D. R., (U. of Texas)
      !            Krogh, F. T., (JPL)
      ! ***DESCRIPTION
      !
      !                 B L A S  Subprogram
      !     Description of Parameters
      !
      !      --Input--
      !         N  number of elements in input vector(s)
      !        DX  double precision vector with N elements
      !      INCX  storage spacing between elements of DX
      !        DY  double precision vector with N elements
      !      INCY  storage spacing between elements of DY
      !
      !      --Output--
      !      DDOT  double precision dot product (zero if N .LE. 0)
      !
      !      Returns the dot product of double precision DX and DY.
      !      DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
      !      where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
      !      defined in a similar way using INCY.
      !
      ! ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
      !                  Krogh, Basic linear algebra subprograms for Fortran
      !                  usage, Algorithm No. 539, Transactions on Mathematical
      !                  Software 5, 3 (September 1979), pp. 308-323.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    791001  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
      ! ***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
    !
    !      Code for unequal or nonpositive increments.
    !
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
   !
   !      Code for both increments equal to 1.
   !
   !      Clean-up loop so remaining vector length is a multiple of 5.
   !
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +&
                    DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
   !
   !      Code for equal, positive, non-unit increments.
   !
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
