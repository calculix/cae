      ! DECK DSLUGM
      SUBROUTINE DSLUGM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,&
         TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
      ! ***BEGIN PROLOGUE  DSLUGM
      ! ***PURPOSE  Incomplete LU GMRES iterative sparse Ax=b solver.
      !             This routine uses the generalized minimum residual
      !             (GMRES) method with incomplete LU factorization for
      !             preconditioning to solve possibly non-symmetric linear
      !             systems of the form: Ax = b.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2A4, D2B4
      ! ***TYPE      DOUBLE PRECISION (SSLUGM-S, DSLUGM-D)
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
      !       INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
      !       INTEGER   ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      !       DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      !
      !       CALL DSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
      !      $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
      !      $     RWORK, LENW, IWORK, LENIW)
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
      !          These arrays should hold the matrix A in either the SLAP
      !          Triad format or the SLAP Column format.  See "Description",
      !          below.  If the SLAP Triad format is chosen it is changed
      !          internally to the SLAP Column format.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the upper
      !          or lower triangle of the matrix is stored.
      !  NSAVE  :IN       Integer.
      !          Number of direction vectors to save and orthogonalize against.
      !          Must be greater than 1.
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
      !  TOL    :INOUT    Double Precision.
      !          Convergence criterion, as described below.  If TOL is set
      !          to zero on input, then a default value of 500*(the smallest
      !          positive magnitude, machine epsilon) is used.
      !  ITMAX  :IN       Integer.
      !          Maximum number of iterations.  This routine uses the default
      !          of NRMAX = ITMAX/NSAVE to determine the when each restart
      !          should occur.  See the description of NRMAX and MAXL in
      !          DGMRES for a full and frightfully interesting discussion of
      !          this topic.
      !  ITER   :OUT      Integer.
      !          Number of iterations required to reach convergence, or
      !          ITMAX+1 if convergence criterion could not be achieved in
      !          ITMAX iterations.
      !  ERR    :OUT      Double Precision.
      !          Error estimate of error in final approximate solution, as
      !          defined by ITOL.  Letting norm() denote the Euclidean
      !          norm, ERR is defined as follows...
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
      !                IERR = 2 => Routine DPIGMR failed to reduce the norm
      !                            of the current residual on its last call,
      !                            and so the iteration has stalled.  In
      !                            this case, X equals the last computed
      !                            approximation.  The user must either
      !                            increase MAXL, or choose a different
      !                            initial guess.
      !                IERR =-1 => Insufficient length for RGWK array.
      !                            IGWK(6) contains the required minimum
      !                            length of the RGWK array.
      !                IERR =-2 => Inconsistent ITOL and JPRE values.
      !          For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
      !          left-hand-side of the relevant stopping test defined
      !          below associated with the residual for the current
      !          approximation X(L).
      !  IUNIT  :IN       Integer.
      !          Unit number on which to write the error at each iteration,
      !          if this is desired for monitoring convergence.  If unit
      !          number is 0, no writing will occur.
      !  RWORK  :WORK    Double Precision RWORK(LENW).
      !          Double Precision array of size LENW.
      !  LENW   :IN       Integer.
      !          Length of the double precision workspace, RWORK.
      !          LENW >= 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NL+NU.
      !          Here NL is the number of non-zeros in the lower triangle of
      !          the matrix (including the diagonal) and NU is the number of
      !          non-zeros in the upper triangle of the matrix (including the
      !          diagonal).
      !          For the recommended values,  RWORK  has size at least
      !          131 + 17*N + NL + NU.
      !  IWORK  :INOUT    Integer IWORK(LENIW).
      !          Used to hold pointers into the RWORK array.
      !          Upon return the following locations of IWORK hold information
      !          which may be of use to the user:
      !          IWORK(9)  Amount of Integer workspace actually used.
      !          IWORK(10) Amount of Double Precision workspace actually used.
      !  LENIW  :IN       Integer.
      !          Length of the integer workspace, IWORK.
      !          LENIW >= NL+NU+4*N+32.
      !
      !  *Description:
      !        DSLUGM solves a linear system A*X = B rewritten in the form:
      !
      !         (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
      !
      !        with right preconditioning, or
      !
      !         (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
      !
      !        with left preconditioning, where A is an n-by-n double precision
      !        matrix, X and B are N-vectors, SB and SX are  diagonal scaling
      !        matrices, and M is the Incomplete LU factorization of A.  It
      !        uses preconditioned  Krylov subpace   methods  based on  the
      !        generalized minimum residual  method (GMRES).   This routine
      !        is a  driver  routine  which  assumes a SLAP   matrix   data
      !        structure   and  sets  up  the  necessary  information to do
      !        diagonal  preconditioning  and calls the main GMRES  routine
      !        DGMRES for the   solution   of the linear   system.   DGMRES
      !        optionally   performs  either  the full    orthogonalization
      !        version of the  GMRES algorithm or  an incomplete variant of
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
      !        the user.  If TOL equals zero  when DSLUGM is called, then a
      !        default  value  of 500*(the   smallest  positive  magnitude,
      !        machine epsilon) is used.  If the  scaling arrays SB  and SX
      !        are used, then  ideally they  should be chosen  so  that the
      !        vectors SX*X(or SX*M*X) and  SB*B have all their  components
      !        approximately equal  to  one in  magnitude.  If one wants to
      !        use the same scaling in X  and B, then  SB and SX can be the
      !        same array in the calling program.
      !
      !        The following is a list of the other routines and their
      !        functions used by GMRES:
      !        DGMRES  Contains the matrix structure independent driver
      !                routine for GMRES.
      !        DPIGMR  Contains the main iteration loop for GMRES.
      !        DORTH   Orthogonalizes a new vector against older basis vectors.
      !        DHEQR   Computes a QR decomposition of a Hessenberg matrix.
      !        DHELS   Solves a Hessenberg least-squares system, using QR
      !                factors.
      !        RLCALC  Computes the scaled residual RL.
      !        XLCALC  Computes the solution XL.
      !        ISDGMR  User-replaceable stopping routine.
      !
      !        The Sparse Linear Algebra Package (SLAP) utilizes two matrix
      !        data structures: 1) the  SLAP Triad  format or  2)  the SLAP
      !        Column format.  The user can hand this routine either of the
      !        of these data structures and SLAP  will figure out  which on
      !        is being used and act accordingly.
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
      !  *Side Effects:
      !        The SLAP Triad format (IA, JA, A) is modified internally to be
      !        the SLAP Column format.  See above.
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
      ! ***ROUTINES CALLED  DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    920407  COMMON BLOCK renamed DSLBLK.  (WRB)
      !    920511  Added complete declaration section.  (WRB)
      !    920929  Corrected format of references.  (FNF)
      !    921019  Corrected NEL to NL.  (FNF)
      ! ***END PROLOGUE  DSLUGM
      !          The following is for optimized compilation on LLNL/LTSS Crays.
      ! LLL. OPTIMIZE
      !      .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
      !      .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,&
              NELT, NSAVE
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
      !      .. Local Scalars ..
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCIGW, LOCIL, LOCIU, LOCIW,&
              LOCJL, LOCJU, LOCL, LOCNC, LOCNR, LOCRGW, LOCU, LOCW,&
              MYITOL, NL, NU
      !      .. External Subroutines ..
      EXTERNAL DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV
      ! ***FIRST EXECUTABLE STATEMENT  DSLUGM
      !
      IERR = 0
      ERR  = 0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      !
      !          Change the SLAP input matrix IA, JA, A to SLAP-Column format.
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
      !
      !          Count number of Non-Zero elements preconditioner ILU matrix.
      !          Then set up the work arrays.  We assume MAXL=KMP=NSAVE.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
         !          Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            ! VD$ NOVECTOR
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
      !
      LOCIGW = LOCIB
      LOCIL = LOCIGW + 20
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
      !
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCRGW = LOCU + NU
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
      !
      !          Check the workspace allocations.
      CALL DCHKW( 'DSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
      !
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
      !
      !          Compute the Incomplete LU decomposition.
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),&
           IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),&
           IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
      !
      !          Perform the Incomplete LU Preconditioned Generalized Minimum
      !          Residual iteration algorithm.  The following DGMRES
      !          defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
      !          JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
      !
      CALL DGMRES( N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI,&
           MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,&
           RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,&
           RWORK, IWORK )
      !
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
      ! ------------- LAST LINE OF DSLUGM FOLLOWS ----------------------------
      END
      ! DECK DSLUI
      SUBROUTINE DSLUI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      ! ***BEGIN PROLOGUE  DSLUI
      ! ***PURPOSE  SLAP MSOLVE for LDU Factorization.
      !             This routine acts as an interface between the SLAP generic
      !             MSOLVE calling convention and the routine that actually
      !                            -1
      !             computes  (LDU)  B = X.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2E
      ! ***TYPE      DOUBLE PRECISION (SSLUI-S, DSLUI-D)
      ! ***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
      !              SLAP, SPARSE
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !        It is assumed that RWORK and IWORK have initialized with
      !        the information required for DSLUI2:
      !           IWORK(1) = Starting location of IL in IWORK.
      !           IWORK(2) = Starting location of JL in IWORK.
      !           IWORK(3) = Starting location of IU in IWORK.
      !           IWORK(4) = Starting location of JU in IWORK.
      !           IWORK(5) = Starting location of L in RWORK.
      !           IWORK(6) = Starting location of DINV in RWORK.
      !           IWORK(7) = Starting location of U in RWORK.
      !        See the DESCRIPTION of DSLUI2 for details.
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  DSLUI2
      ! ***REVISION HISTORY  (YYMMDD)
      !    871119  DATE WRITTEN
      !    881213  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    920511  Added complete declaration section.  (WRB)
      !    921113  Corrected C***CATEGORY line.  (FNF)
      !    930701  Updated CATEGORY section.  (FNF, WRB)
      ! ***END PROLOGUE  DSLUI
      !      .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(10), JA(NELT)
      !      .. Local Scalars ..
      INTEGER LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU
      !      .. External Subroutines ..
      EXTERNAL DSLUI2
      ! ***FIRST EXECUTABLE STATEMENT  DSLUI
      !
      !          Pull out the locations of the arrays holding the ILU
      !          factorization.
      !
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
      !
      !          Solve the system LUx = b
      CALL DSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL),&
           RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) )
      !
      RETURN
      ! ------------- LAST LINE OF DSLUI FOLLOWS ----------------------------
      END
      ! DECK DSLUI2
      SUBROUTINE DSLUI2 (N, B, X, IL, JL, L, DINV, IU, JU, U)
      use omp_lib
      ! ***BEGIN PROLOGUE  DSLUI2
      ! ***PURPOSE  SLAP Backsolve for LDU Factorization.
      !             Routine to solve a system of the form  L*D*U X = B,
      !             where L is a unit lower triangular matrix, D is a diagonal
      !             matrix, and U is a unit upper triangular matrix.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2E
      ! ***TYPE      DOUBLE PRECISION (SSLUI2-S, DSLUI2-D)
      ! ***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
      !              SLAP, SPARSE
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER N, IL(NL), JL(NL), IU(NU), JU(NU)
      !      DOUBLE PRECISION B(N), X(N), L(NL), DINV(N), U(NU)
      !
      !      CALL DSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U )
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          Order of the Matrix.
      !  B      :IN       Double Precision B(N).
      !          Right hand side.
      !  X      :OUT      Double Precision X(N).
      !          Solution of L*D*U x = b.
      !  IL     :IN       Integer IL(NL).
      !  JL     :IN       Integer JL(NL).
      !  L      :IN       Double Precision L(NL).
      !          IL, JL, L contain the unit  lower triangular factor of the
      !          incomplete decomposition of some matrix stored in SLAP Row
      !          format.  The diagonal of ones *IS* stored.  This structure
      !          can   be   set  up  by   the  DSILUS  routine.   See   the
      !          "Description", below  for more   details about   the  SLAP
      !          format.  (NL is the number of non-zeros in the L array.)
      !  DINV   :IN       Double Precision DINV(N).
      !          Inverse of the diagonal matrix D.
      !  IU     :IN       Integer IU(NU).
      !  JU     :IN       Integer JU(NU).
      !  U      :IN       Double Precision U(NU).
      !          IU, JU, U contain the unit upper triangular factor  of the
      !          incomplete decomposition  of  some  matrix stored in  SLAP
      !          Column format.   The diagonal of ones  *IS* stored.   This
      !          structure can be set up  by the DSILUS routine.  See   the
      !          "Description", below   for  more   details about  the SLAP
      !          format.  (NU is the number of non-zeros in the U array.)
      !
      !  *Description:
      !        This routine is supplied with  the SLAP package as a routine
      !        to  perform  the  MSOLVE operation  in   the  SIR and   SBCG
      !        iteration routines for  the  drivers DSILUR and DSLUBC.   It
      !        must  be called  via   the  SLAP  MSOLVE  calling   sequence
      !        convention interface routine DSLUI.
      !          **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
      !                **** SLAP MSOLVE CALLING CONVENTION ****
      !
      !        IL, JL, L should contain the unit lower triangular factor of
      !        the incomplete decomposition of the A matrix  stored in SLAP
      !        Row format.  IU, JU, U should contain  the unit upper factor
      !        of the  incomplete decomposition of  the A matrix  stored in
      !        SLAP Column format This ILU factorization can be computed by
      !        the DSILUS routine. The diagonals (which are all one's) are
      !        stored.
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
      !        ==================== S L A P Row format ====================
      !
      !        This routine requires  that the matrix A  be  stored  in the
      !        SLAP  Row format.   In this format  the non-zeros are stored
      !        counting across  rows (except for the diagonal  entry, which
      !        must  appear first  in each  "row")  and  are stored  in the
      !        double precision  array A.  In other words, for each row  in
      !        the matrix  put the diagonal  entry in A.   Then put in  the
      !        other  non-zero elements  going across  the row  (except the
      !        diagonal) in order.  The JA array holds the column index for
      !        each non-zero.  The IA array holds the offsets  into the JA,
      !        A  arrays  for  the   beginning  of  each  row.    That  is,
      !        JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
      !        th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
      !        are  the last elements  of the  IROW-th row.   Note  that we
      !        always have  IA(N+1) = NELT+1, where N is the number of rows
      !        in the matrix  and  NELT is the  number of non-zeros  in the
      !        matrix.
      !
      !        Here is an example of the SLAP Row storage format for a  5x5
      !        Matrix (in the A and JA arrays '|' denotes the end of a row):
      !
      !            5x5 Matrix         SLAP Row format for 5x5 matrix on left.
      !                               1  2  3    4  5    6  7    8    9 10 11
      !        |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
      !        |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
      !        | 0  0 33  0 35|  IA:  1  4  6    8  9   12
      !        | 0  0  0 44  0|
      !        |51  0 53  0 55|
      !
      !        With  the SLAP  format  the "inner  loops" of  this  routine
      !        should vectorize   on machines with   hardware  support  for
      !        vector gather/scatter operations.  Your compiler may require
      !        a  compiler directive  to  convince   it that there  are  no
      !        implicit vector  dependencies.  Compiler directives  for the
      !        Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
      !        with the standard SLAP distribution.
      !
      ! ***SEE ALSO  DSILUS
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    871119  DATE WRITTEN
      !    881213  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    920511  Added complete declaration section.  (WRB)
      !    921113  Corrected C***CATEGORY line.  (FNF)
      !    930701  Updated CATEGORY section.  (FNF, WRB)
      ! ***END PROLOGUE  DSLUI2
      !      .. Scalar Arguments ..
      INTEGER N
      !      .. Array Arguments ..
      DOUBLE PRECISION B(N), DINV(N), L(*), U(*), X(N)
      INTEGER IL(*), IU(*), JL(*), JU(*)
      !      .. Local Scalars ..
      INTEGER I, ICOL, IROW, J, JBGN, JEND
      ! ***FIRST EXECUTABLE STATEMENT  DSLUI2
      !
      !          Solve  L*Y = B,  storing result in X, L stored by rows.
      !
      ! $omp parallel default(none)
      ! $omp& shared(x,b,n)
      ! $omp& private(i)
      ! $omp do
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
      ! $omp end do
      ! $omp end parallel
      !
      DO 30 IROW = 2, N
         JBGN = IL(IROW)
         JEND = IL(IROW+1)-1
         IF( JBGN.LE.JEND ) THEN
            ! LLL. OPTION ASSERT (NOHAZARD)
            ! DIR$ IVDEP
            ! VD$ ASSOC
            ! VD$ NODEPCHK
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW) - L(J)*X(JL(J))
 20         CONTINUE
         ENDIF
 30   CONTINUE
      !
      !
      !          Solve  D*Z = Y,  storing result in X.
      ! $omp parallel default(none)
      ! $omp& shared(n,x,dinv)
      ! $omp& private(i)
      ! $omp do
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
      ! $omp end do
      ! $omp end parallel
      !
      !          Solve  U*X = Z, U stored by columns.
      DO 60 ICOL = N, 2, -1
         JBGN = JU(ICOL)
         JEND = JU(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            ! LLL. OPTION ASSERT (NOHAZARD)
            ! DIR$ IVDEP
            ! VD$ NODEPCHK
            DO 50 J = JBGN, JEND
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
      !
      RETURN
      ! ------------- LAST LINE OF DSLUI2 FOLLOWS ----------------------------
      END
      ! DECK DSMV
      SUBROUTINE DSMV (N, X, Y, NELT, IA, JA, A, ISYM)
      use omp_lib
      ! ***BEGIN PROLOGUE  DSMV
      ! ***PURPOSE  SLAP Column Format Sparse Matrix Vector Product.
      !             Routine to calculate the sparse matrix vector product:
      !             Y = A*X.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D1B4
      ! ***TYPE      DOUBLE PRECISION (SSMV-S, DSMV-D)
      ! ***KEYWORDS  MATRIX VECTOR MULTIPLY, SLAP, SPARSE
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM
      !      DOUBLE PRECISION X(N), Y(N), A(NELT)
      !
      !      CALL DSMV(N, X, Y, NELT, IA, JA, A, ISYM )
      !
      !  *Arguments:
      !  N      :IN       Integer.
      !          Order of the Matrix.
      !  X      :IN       Double Precision X(N).
      !          The vector that should be multiplied by the matrix.
      !  Y      :OUT      Double Precision Y(N).
      !          The product of the matrix and the vector.
      !  NELT   :IN       Integer.
      !          Number of Non-Zeros stored in A.
      !  IA     :IN       Integer IA(NELT).
      !  JA     :IN       Integer JA(NELT).
      !  A      :IN       Double Precision A(NELT).
      !          These arrays should hold the matrix A in the SLAP Column
      !          format.  See "Description", below.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the upper
      !          or lower triangle of the matrix is stored.
      !
      !  *Description
      !        =================== S L A P Column format ==================
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
      !        With  the SLAP  format  the "inner  loops" of  this  routine
      !        should vectorize   on machines with   hardware  support  for
      !        vector gather/scatter operations.  Your compiler may require
      !        a  compiler directive  to  convince   it that there  are  no
      !        implicit vector  dependencies.  Compiler directives  for the
      !        Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
      !        with the standard SLAP distribution.
      !
      !  *Cautions:
      !      This   routine   assumes  that  the matrix A is stored in SLAP
      !      Column format.  It does not check  for  this (for  speed)  and
      !      evil, ugly, ornery and nasty things  will happen if the matrix
      !      data  structure  is,  in fact, not SLAP Column.  Beware of the
      !      wrong data structure!!!
      !
      ! ***SEE ALSO  DSMTV
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    871119  DATE WRITTEN
      !    881213  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    920511  Added complete declaration section.  (WRB)
      !    930701  Updated CATEGORY section.  (FNF, WRB)
      ! ***END PROLOGUE  DSMV
      !      .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), X(N), Y(N)
      INTEGER IA(NELT), JA(NELT)
      !      .. Local Scalars ..
      INTEGER I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND
      ! ***FIRST EXECUTABLE STATEMENT  DSMV
      if(isym.ne.1) then
         !
         !          Zero out the result vector.
         !
         ! $omp parallel default(none)
         ! $omp& shared(n,y,ja,ia,a,x)
         ! $omp& private(i,icol,ibgn,iend)
         ! $omp do
         DO 10 I = 1, N
            Y(I) = 0
 10      CONTINUE
         ! $omp end do
         !
         !          Multiply by A.
         !
         ! $omp do
         DO 30 ICOL = 1, N
            IBGN = JA(ICOL)
            IEND = JA(ICOL+1)-1
            DO 20 I = IBGN, IEND
               ! $omp atomic
               Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20         CONTINUE
 30      CONTINUE
      ! $omp end do
      ! $omp end parallel
      !
      else
         !
         !          The matrix is symmetric.  Need to get the other half in...
         !          This loops assumes that the diagonal is the first entry in
         !          each column.
         !
         !
         !          Zero out the result vector.
         !
         ! $omp parallel default(none)
         ! $omp& shared(n,y,ja,ia,a,x)
         ! $omp& private(i,j,icol,irow,ibgn,iend,jbgn,jend)
         ! $omp do
         DO 100 I = 1, N
            Y(I) = 0
 100     CONTINUE
         ! $omp end do
         !
         !          Multiply by A.
         !
         ! $omp do
         DO 300 ICOL = 1, N
            IBGN = JA(ICOL)
            IEND = JA(ICOL+1)-1
            DO 200 I = IBGN, IEND
               ! $omp atomic
               Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 200        CONTINUE
 300     CONTINUE
         ! $omp end do
         !
         ! $omp do
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ! $omp end do
      ! $omp end parallel
      ENDIF
      RETURN
      ! ------------- LAST LINE OF DSMV FOLLOWS ----------------------------
      END
      ! DECK DCHKW
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR)
      ! ***BEGIN PROLOGUE  DCHKW
      ! ***SUBSIDIARY
      ! ***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
      !             This routine checks the work array lengths and interfaces
      !             to the SLATEC error handler if a problem is found.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  R2
      ! ***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D)
      ! ***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING
      ! ***AUTHOR  Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      CHARACTER*(*) NAME
      !      INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
      !      DOUBLE PRECISION ERR
      !
      !      CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      !
      !  *Arguments:
      !  NAME   :IN       Character*(*).
      !          Name of the calling routine.  This is used in the output
      !          message, if an error is detected.
      !  LOCIW  :IN       Integer.
      !          Location of the first free element in the integer workspace
      !          array.
      !  LENIW  :IN       Integer.
      !          Length of the integer workspace array.
      !  LOCW   :IN       Integer.
      !          Location of the first free element in the double precision
      !          workspace array.
      !  LENRW  :IN       Integer.
      !          Length of the double precision workspace array.
      !  IERR   :OUT      Integer.
      !          Return error flag.
      !                IERR = 0 => All went well.
      !                IERR = 1 => Insufficient storage allocated for
      !                            WORK or IWORK.
      !  ITER   :OUT      Integer.
      !          Set to zero on return.
      !  ERR    :OUT      Double Precision.
      !          Set to the smallest positive magnitude if all went well.
      !          Set to a very large number if an error is detected.
      !
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  D1MACH, XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    880225  DATE WRITTEN
      !    881213  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI
      !            X3.9-1978.  (FNF)
      !    910506  Made subsidiary.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      !    921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)
      ! ***END PROLOGUE  DCHKW
      !      .. Scalar Arguments ..
      DOUBLE PRECISION ERR
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW
      CHARACTER NAME*(*)
      !      .. Local Scalars ..
      CHARACTER XERN1*8, XERN2*8, XERNAM*8
      !      .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
      !      .. External Subroutines ..
      EXTERNAL XERMSG
      ! ***FIRST EXECUTABLE STATEMENT  DCHKW
      !
      !          Check the Integer workspace situation.
      !
      IERR = 0
      ITER = 0
      ERR = D1MACH(1)
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCIW
         WRITE (XERN2, '(I8)') LENIW
         CALL XERMSG ('SLATEC', 'DCHKW',&
            'In ' // XERNAM // ', INTEGER work array too short.  ' //&
            'IWORK needs ' // XERN1 // '; have allocated ' // XERN2,&
            1, 1)
      ENDIF
      !
      !          Check the Double Precision workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCW
         WRITE (XERN2, '(I8)') LENW
         CALL XERMSG ('SLATEC', 'DCHKW',&
            'In ' // XERNAM // ', DOUBLE PRECISION work array too ' //&
            'short.  RWORK needs ' // XERN1 // '; have allocated ' //&
            XERN2, 1, 1)
      ENDIF
      RETURN
      ! ------------- LAST LINE OF DCHKW FOLLOWS ----------------------------
      END
      ! DECK DS2Y
      SUBROUTINE DS2Y (N, NELT, IA, JA, A, ISYM)
      ! ***BEGIN PROLOGUE  DS2Y
      ! ***PURPOSE  SLAP Triad to SLAP Column Format Converter.
      !             Routine to convert from the SLAP Triad to SLAP Column
      !             format.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D1B9
      ! ***TYPE      DOUBLE PRECISION (SS2Y-S, DS2Y-D)
      ! ***KEYWORDS  LINEAR SYSTEM, SLAP SPARSE
      ! ***AUTHOR  Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      !      DOUBLE PRECISION A(NELT)
      !
      !      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          Order of the Matrix.
      !  NELT   :IN       Integer.
      !          Number of non-zeros stored in A.
      !  IA     :INOUT    Integer IA(NELT).
      !  JA     :INOUT    Integer JA(NELT).
      !  A      :INOUT    Double Precision A(NELT).
      !          These arrays should hold the matrix A in either the SLAP
      !          Triad format or the SLAP Column format.  See "Description",
      !          below.  If the SLAP Triad format is used, this format is
      !          translated to the SLAP Column format by this routine.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the lower
      !          triangle of the matrix is stored.
      !
      !  *Description:
      !        The Sparse Linear Algebra Package (SLAP) utilizes two matrix
      !        data structures: 1) the  SLAP Triad  format or  2)  the SLAP
      !        Column format.  The user can hand this routine either of the
      !        of these data structures.  If the SLAP Triad format is give
      !        as input then this routine transforms it into SLAP Column
      !        format.  The way this routine tells which format is given as
      !        input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
      !        have the SLAP Column format.  If that equality does not hold
      !        then it is assumed that the IA, JA, A arrays contain the
      !        SLAP Triad format.
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
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  QS2I1D
      ! ***REVISION HISTORY  (YYMMDD)
      !    871119  DATE WRITTEN
      !    881213  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      !    930701  Updated CATEGORY section.  (FNF, WRB)
      ! ***END PROLOGUE  DS2Y
      !      .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT)
      INTEGER IA(NELT), JA(NELT)
      !      .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, ITEMP, J
      !      .. External Subroutines ..
      EXTERNAL QS2I1D
      ! ***FIRST EXECUTABLE STATEMENT  DS2Y
      !
      !          Check to see if the (IA,JA,A) arrays are in SLAP Column
      !          format.  If it's not then transform from SLAP Triad.
      !
      IF( JA(N+1).EQ.NELT+1 ) RETURN
      !
      !          Sort into ascending order by COLUMN (on the ja array).
      !          This will line up the columns.
      !
      CALL QS2I1D( JA, IA, A, NELT, 1 )
      !
      !          Loop over each column to see where the column indices change
      !          in the column index array ja.  This marks the beginning of the
      !          next column.
      !
      ! VD$R NOVECTOR
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
      !
      !          Mark the n+2 element so that future calls to a SLAP routine
      !          utilizing the YSMP-Column storage format will be able to tell.
      !
      JA(N+2) = 0
      !
      !          Now loop through the IA array making sure that the diagonal
      !          matrix element appears first in the column.  Then sort the
      !          rest of the column in ascending order.
      !
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
               !
               !               Swap the diagonal element with the first element in the
               !               column.
               !
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
      ! ------------- LAST LINE OF DS2Y FOLLOWS ----------------------------
      END
      ! DECK DSILUS
      SUBROUTINE DSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV,&
         NU, IU, JU, U, NROW, NCOL)
      ! ***BEGIN PROLOGUE  DSILUS
      ! ***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
      !             Routine to generate the incomplete LDU decomposition of a
      !             matrix.  The unit lower triangular factor L is stored by
      !             rows and the unit upper triangular factor U is stored by
      !             columns.  The inverse of the diagonal matrix D is stored.
      !             No fill in is allowed.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  D2E
      ! ***TYPE      DOUBLE PRECISION (SSILUS-S, DSILUS-D)
      ! ***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION,
      !              NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
      ! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
      !            Seager, Mark K., (LLNL)
      !              Lawrence Livermore National Laboratory
      !              PO BOX 808, L-60
      !              Livermore, CA 94550 (510) 423-3141
      !              seager@llnl.gov
      ! ***DESCRIPTION
      !
      !  *Usage:
      !      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      !      INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)
      !      INTEGER NROW(N), NCOL(N)
      !      DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
      !
      !      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,
      !     $    DINV, NU, IU, JU, U, NROW, NCOL )
      !
      !  *Arguments:
      !  N      :IN       Integer
      !          Order of the Matrix.
      !  NELT   :IN       Integer.
      !          Number of elements in arrays IA, JA, and A.
      !  IA     :IN       Integer IA(NELT).
      !  JA     :IN       Integer JA(NELT).
      !  A      :IN       Double Precision A(NELT).
      !          These arrays should hold the matrix A in the SLAP Column
      !          format.  See "Description", below.
      !  ISYM   :IN       Integer.
      !          Flag to indicate symmetric storage format.
      !          If ISYM=0, all non-zero entries of the matrix are stored.
      !          If ISYM=1, the matrix is symmetric, and only the lower
      !          triangle of the matrix is stored.
      !  NL     :OUT      Integer.
      !          Number of non-zeros in the L array.
      !  IL     :OUT      Integer IL(NL).
      !  JL     :OUT      Integer JL(NL).
      !  L      :OUT      Double Precision L(NL).
      !          IL, JL, L  contain the unit lower triangular factor of  the
      !          incomplete decomposition  of some  matrix stored  in   SLAP
      !          Row format.     The   Diagonal  of ones  *IS*  stored.  See
      !          "DESCRIPTION", below for more details about the SLAP format.
      !  NU     :OUT      Integer.
      !          Number of non-zeros in the U array.
      !  IU     :OUT      Integer IU(NU).
      !  JU     :OUT      Integer JU(NU).
      !  U      :OUT      Double Precision     U(NU).
      !          IU, JU, U contain   the unit upper triangular factor of the
      !          incomplete  decomposition    of some matrix  stored in SLAP
      !          Column  format.   The Diagonal of ones   *IS*  stored.  See
      !          "Description", below  for  more  details  about  the   SLAP
      !          format.
      !  NROW   :WORK     Integer NROW(N).
      !          NROW(I) is the number of non-zero elements in the I-th row
      !          of L.
      !  NCOL   :WORK     Integer NCOL(N).
      !          NCOL(I) is the number of non-zero elements in the I-th
      !          column of U.
      !
      !  *Description
      !        IL, JL, L should contain the unit  lower triangular factor of
      !        the incomplete decomposition of the A matrix  stored in SLAP
      !        Row format.  IU, JU, U should contain  the unit upper factor
      !        of the  incomplete decomposition of  the A matrix  stored in
      !        SLAP Column format This ILU factorization can be computed by
      !        the DSILUS routine. The diagonals (which are all one's) are
      !        stored.
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
      !        ==================== S L A P Row format ====================
      !
      !        This routine requires  that the matrix A  be  stored  in the
      !        SLAP  Row format.   In this format  the non-zeros are stored
      !        counting across  rows (except for the diagonal  entry, which
      !        must  appear first  in each  "row")  and  are stored  in the
      !        double precision  array A.  In other words, for each row  in
      !        the matrix  put the diagonal  entry in A.   Then put in  the
      !        other  non-zero elements  going across  the row  (except the
      !        diagonal) in order.  The JA array holds the column index for
      !        each non-zero.  The IA array holds the offsets  into the JA,
      !        A  arrays  for  the   beginning  of  each  row.    That  is,
      !        JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
      !        th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
      !        are  the last elements  of the  IROW-th row.   Note  that we
      !        always have  IA(N+1) = NELT+1, where N is the number of rows
      !        in the matrix  and  NELT is the  number of non-zeros  in the
      !        matrix.
      !
      !        Here is an example of the SLAP Row storage format for a  5x5
      !        Matrix (in the A and JA arrays '|' denotes the end of a row):
      !
      !            5x5 Matrix         SLAP Row format for 5x5 matrix on left.
      !                               1  2  3    4  5    6  7    8    9 10 11
      !        |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
      !        |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
      !        | 0  0 33  0 35|  IA:  1  4  6    8  9   12
      !        | 0  0  0 44  0|
      !        |51  0 53  0 55|
      !
      ! ***SEE ALSO  SILUR
      ! ***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
      !                   Johns Hopkins University Press, Baltimore, Maryland,
      !                   1983.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    890404  DATE WRITTEN
      !    890404  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    920511  Added complete declaration section.  (WRB)
      !    920929  Corrected format of reference.  (FNF)
      !    930701  Updated CATEGORY section.  (FNF, WRB)
      ! ***END PROLOGUE  DSILUS
      !      .. Scalar Arguments ..
      INTEGER ISYM, N, NELT, NL, NU
      !      .. Array Arguments ..
      DOUBLE PRECISION A(NELT), DINV(N), L(NL), U(NU)
      INTEGER IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU),&
              NCOL(N), NROW(N)
      !      .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2,&
              INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC,&
              KR
      ! ***FIRST EXECUTABLE STATEMENT  DSILUS
      !
      !          Count number of elements in each row of the lower triangle.
      !
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
      ! VD$R NOCONCUR
      ! VD$R NOVECTOR
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
      !
      !          Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
                  !          Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  !          Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
                     !          Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
      !
      !          Sort the rows of L and the columns of U.
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
      !
      !          Perform the incomplete LDU decomposition.
      DO 300 I=2,N
         !
         !            I-th row of L
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
 !
 !          I-th column of U
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
 !
 !          I-th diagonal element
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
 !
 300  CONTINUE
      !
      !          Replace diagonal elements by their inverses.
      ! VD$ VECTOR
      DO 430 I=1,N
         DINV(I) = 1.0D0/DINV(I)
 430  CONTINUE
      !
      RETURN
      ! ------------- LAST LINE OF DSILUS FOLLOWS ----------------------------
      END
      ! DECK QS2I1D
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)
      ! ***BEGIN PROLOGUE  QS2I1D
      ! ***SUBSIDIARY
      ! ***PURPOSE  Sort an integer array, moving an integer and DP array.
      !             This routine sorts the integer array IA and makes the same
      !             interchanges in the integer array JA and the double pre-
      !             cision array A.  The array IA may be sorted in increasing
      !             order or decreasing order.  A slightly modified QUICKSORT
      !             algorithm is used.
      ! ***LIBRARY   SLATEC (SLAP)
      ! ***CATEGORY  N6A2A
      ! ***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
      ! ***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
      ! ***AUTHOR  Jones, R. E., (SNLA)
      !            Kahaner, D. K., (NBS)
      !            Seager, M. K., (LLNL) seager@llnl.gov
      !            Wisniewski, J. A., (SNLA)
      ! ***DESCRIPTION
      !      Written by Rondall E Jones
      !      Modified by John A. Wisniewski to use the Singleton QUICKSORT
      !      algorithm. date 18 November 1976.
      !
      !      Further modified by David K. Kahaner
      !      National Bureau of Standards
      !      August, 1981
      !
      !      Even further modification made to bring the code up to the
      !      Fortran 77 level and make it more readable and to carry
      !      along one integer array and one double precision array during
      !      the sort by
      !      Mark K. Seager
      !      Lawrence Livermore National Laboratory
      !      November, 1987
      !      This routine was adapted from the ISORT routine.
      !
      !      ABSTRACT
      !          This routine sorts an integer array IA and makes the same
      !          interchanges in the integer array JA and the double precision
      !          array A.
      !          The array IA may be sorted in increasing order or decreasing
      !          order.  A slightly modified quicksort algorithm is used.
      !
      !      DESCRIPTION OF PARAMETERS
      !         IA - Integer array of values to be sorted.
      !         JA - Integer array to be carried along.
      !          A - Double Precision array to be carried along.
      !          N - Number of values in integer array IA to be sorted.
      !      KFLAG - Control parameter
      !            = 1 means sort IA in INCREASING order.
      !            =-1 means sort IA in DECREASING order.
      !
      ! ***SEE ALSO  DS2Y
      ! ***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
      !                  for Sorting With Minimal Storage, Communications ACM
      !                  12:3 (1969), pp.185-7.
      ! ***ROUTINES CALLED  XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    761118  DATE WRITTEN
      !    890125  Previous REVISION DATE
      !    890915  Made changes requested at July 1989 CML Meeting.  (MKS)
      !    890922  Numerous changes to prologue to make closer to SLATEC
      !            standard.  (FNF)
      !    890929  Numerous changes to reduce SP/DP differences.  (FNF)
      !    900805  Changed XERROR calls to calls to XERMSG.  (RWC)
      !    910411  Prologue converted to Version 4.0 format.  (BAB)
      !    910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
      !    920511  Added complete declaration section.  (WRB)
      !    920929  Corrected format of reference.  (FNF)
      !    921012  Corrected all f.p. constants to double precision.  (FNF)
      ! ***END PROLOGUE  QS2I1D
      ! VD$R NOVECTOR
      ! VD$R NOCONCUR
      !      .. Scalar Arguments ..
      INTEGER KFLAG, N
      !      .. Array Arguments ..
      DOUBLE PRECISION A(N)
      INTEGER IA(N), JA(N)
      !      .. Local Scalars ..
      DOUBLE PRECISION R, TA, TTA
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
      !      .. Local Arrays ..
      INTEGER IL(21), IU(21)
      !      .. External Subroutines ..
      EXTERNAL XERMSG
      !      .. Intrinsic Functions ..
      INTRINSIC ABS, INT
      ! ***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N
      IF (NN.LT.1) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',&
            'The number of values to be sorted was not positive.', 1, 1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',&
            'The sort control parameter, K, was not 1 or -1.', 2, 1)
         RETURN
      ENDIF
      !
      !      Alter array IA to get decreasing order if needed.
      !
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
      !
      !      Sort IA and carry JA and A along.
      !      And now...Just a little black magic...
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
      !
      !      Select a central element of the array and save it in location
      !      it, jt, at.
      !
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
      !
      !      If first element of array is greater than it, interchange with it.
      !
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
      !
      !      If last element of array is less than it, swap with it.
      !
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
         !
         !      If first element of array is greater than it, swap with it.
         !
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
  !
  !      Find an element in the second half of the array which is
  !      smaller than it.
  !
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
  !
  !      Find an element in the first half of the array which is
  !      greater than it.
  !
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
      !
      !      Interchange these elements.
      !
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
      !
      !      Save upper and lower subscripts of the array yet to be sorted.
      !
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
  !
  !      Begin again on another portion of the unsorted array.
  !
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
  !
  !      Clean up, if necessary.
  !
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      ! ------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END
