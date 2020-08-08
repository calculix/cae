*DECK DSLUGM
      SUBROUTINE DSLUGM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,
     +   TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C***BEGIN PROLOGUE  DSLUGM
C***PURPOSE  Incomplete LU GMRES iterative sparse Ax=b solver.
C            This routine uses the generalized minimum residual
C            (GMRES) method with incomplete LU factorization for
C            preconditioning to solve possibly non-symmetric linear
C            systems of the form: Ax = b.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SSLUGM-S, DSLUGM-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
C      INTEGER   ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
C      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
C
C      CALL DSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
C     $     RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description",
C         below.  If the SLAP Triad format is chosen it is changed
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C NSAVE  :IN       Integer.
C         Number of direction vectors to save and orthogonalize against.
C         Must be greater than 1.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISDGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /DSLBLK/ SOLN( )
C                 If ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Double Precision.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.  This routine uses the default
C         of NRMAX = ITMAX/NSAVE to determine the when each restart
C         should occur.  See the description of NRMAX and MAXL in
C         DGMRES for a full and frightfully interesting discussion of
C         this topic.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX+1 if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows...
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for
C                           RGWK or IGWK.
C               IERR = 2 => Routine DPIGMR failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Inconsistent ITOL and JPRE values.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C RWORK  :WORK    Double Precision RWORK(LENW).
C         Double Precision array of size LENW.
C LENW   :IN       Integer.
C         Length of the double precision workspace, RWORK.
C         LENW >= 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NL+NU.
C         Here NL is the number of non-zeros in the lower triangle of
C         the matrix (including the diagonal) and NU is the number of
C         non-zeros in the upper triangle of the matrix (including the
C         diagonal).
C         For the recommended values,  RWORK  has size at least
C         131 + 17*N + NL + NU.
C IWORK  :INOUT    Integer IWORK(LENIW).
C         Used to hold pointers into the RWORK array.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Double Precision workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.
C         LENIW >= NL+NU+4*N+32.
C
C *Description:
C       DSLUGM solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where A is an n-by-n double precision
C       matrix, X and B are N-vectors, SB and SX are  diagonal scaling
C       matrices, and M is the Incomplete LU factorization of A.  It
C       uses preconditioned  Krylov subpace   methods  based on  the
C       generalized minimum residual  method (GMRES).   This routine
C       is a  driver  routine  which  assumes a SLAP   matrix   data
C       structure   and  sets  up  the  necessary  information to do
C       diagonal  preconditioning  and calls the main GMRES  routine
C       DGMRES for the   solution   of the linear   system.   DGMRES
C       optionally   performs  either  the full    orthogonalization
C       version of the  GMRES algorithm or  an incomplete variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le.
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the Euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when DSLUGM is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by GMRES:
C       DGMRES  Contains the matrix structure independent driver
C               routine for GMRES.
C       DPIGMR  Contains the main iteration loop for GMRES.
C       DORTH   Orthogonalizes a new vector against older basis vectors.
C       DHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       DHELS   Solves a Hessenberg least-squares system, using QR
C               factors.
C       RLCALC  Computes the scaled residual RL.
C       XLCALC  Computes the solution XL.
C       ISDGMR  User-replaceable stopping routine.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
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
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage
C                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-
C                  more National Laboratory Report UCRL-95088, Rev. 1,
C                  Livermore, California, June 1987.
C***ROUTINES CALLED  DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of references.  (FNF)
C   921019  Corrected NEL to NL.  (FNF)
C***END PROLOGUE  DSLUGM
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,
     +        NELT, NSAVE
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
C     .. Local Scalars ..
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCIGW, LOCIL, LOCIU, LOCIW,
     +        LOCJL, LOCJU, LOCL, LOCNC, LOCNR, LOCRGW, LOCU, LOCW,
     +        MYITOL, NL, NU
C     .. External Subroutines ..
      EXTERNAL DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV
C***FIRST EXECUTABLE STATEMENT  DSLUGM
C
      IERR = 0
      ERR  = 0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.  We assume MAXL=KMP=NSAVE.
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
      LOCIGW = LOCIB
      LOCIL = LOCIGW + 20
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCRGW = LOCU + NU
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
C
C         Check the workspace allocations.
      CALL DCHKW( 'DSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
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
C         Perform the Incomplete LU Preconditioned Generalized Minimum
C         Residual iteration algorithm.  The following DGMRES
C         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
C         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
C
      CALL DGMRES( N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI,
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,
     $     RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,
     $     RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF DSLUGM FOLLOWS ----------------------------
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
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(10), JA(NELT)
C     .. Local Scalars ..
      INTEGER LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU
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
*DECK DSLUI2
      SUBROUTINE DSLUI2 (N, B, X, IL, JL, L, DINV, IU, JU, U)
      use omp_lib
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
C     INTEGER N, IL(NL), JL(NL), IU(NU), JU(NU)
C     DOUBLE PRECISION B(N), X(N), L(NL), DINV(N), U(NU)
C
C     CALL DSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right hand side.
C X      :OUT      Double Precision X(N).
C         Solution of L*D*U x = b.
C IL     :IN       Integer IL(NL).
C JL     :IN       Integer JL(NL).
C L      :IN       Double Precision L(NL).
C         IL, JL, L contain the unit  lower triangular factor of the
C         incomplete decomposition of some matrix stored in SLAP Row
C         format.  The diagonal of ones *IS* stored.  This structure
C         can   be   set  up  by   the  DSILUS  routine.   See   the
C         "Description", below  for more   details about   the  SLAP
C         format.  (NL is the number of non-zeros in the L array.)
C DINV   :IN       Double Precision DINV(N).
C         Inverse of the diagonal matrix D.
C IU     :IN       Integer IU(NU).
C JU     :IN       Integer JU(NU).
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
      INTEGER N
C     .. Array Arguments ..
      DOUBLE PRECISION B(N), DINV(N), L(*), U(*), X(N)
      INTEGER IL(*), IU(*), JL(*), JU(*)
C     .. Local Scalars ..
      INTEGER I, ICOL, IROW, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  DSLUI2
C
C         Solve  L*Y = B,  storing result in X, L stored by rows.
C
c$omp parallel default(none)
c$omp& shared(x,b,n)
c$omp& private(i)
c$omp do
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
c$omp end do
c$omp end parallel
c
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
c
C
C         Solve  D*Z = Y,  storing result in X.
c$omp parallel default(none)
c$omp& shared(n,x,dinv)
c$omp& private(i)
c$omp do
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
c$omp end do
c$omp end parallel
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
*DECK DSMV
      SUBROUTINE DSMV (N, X, Y, NELT, IA, JA, A, ISYM)
      use omp_lib
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
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM
C     DOUBLE PRECISION X(N), Y(N), A(NELT)
C
C     CALL DSMV(N, X, Y, NELT, IA, JA, A, ISYM )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C X      :IN       Double Precision X(N).
C         The vector that should be multiplied by the matrix.
C Y      :OUT      Double Precision Y(N).
C         The product of the matrix and the vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
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
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), X(N), Y(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  DSMV
      if(isym.ne.1) then
C
C         Zero out the result vector.
C
c$omp parallel default(none)
c$omp& shared(n,y,ja,ia,a,x)
c$omp& private(i,icol,ibgn,iend)
c$omp do
         DO 10 I = 1, N
            Y(I) = 0
 10      CONTINUE
c$omp end do
C
C         Multiply by A.
C
c$omp do
         DO 30 ICOL = 1, N
            IBGN = JA(ICOL)
            IEND = JA(ICOL+1)-1
            DO 20 I = IBGN, IEND
c$omp atomic
               Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20         CONTINUE
 30      CONTINUE
c$omp end do
c$omp end parallel
C
      else
C
C         The matrix is symmetric.  Need to get the other half in...
C         This loops assumes that the diagonal is the first entry in
C         each column.
C
C
C         Zero out the result vector.
C
c$omp parallel default(none)
c$omp& shared(n,y,ja,ia,a,x)
c$omp& private(i,j,icol,irow,ibgn,iend,jbgn,jend)
c$omp do
         DO 100 I = 1, N
            Y(I) = 0
 100     CONTINUE
c$omp end do
C
C         Multiply by A.
C
c$omp do
         DO 300 ICOL = 1, N
            IBGN = JA(ICOL)
            IEND = JA(ICOL+1)-1
            DO 200 I = IBGN, IEND
c$omp atomic
               Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 200        CONTINUE
 300     CONTINUE
c$omp end do
C
c$omp do
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
c$omp end do
c$omp end parallel
      ENDIF
      RETURN
C------------- LAST LINE OF DSMV FOLLOWS ----------------------------
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
C     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
C     DOUBLE PRECISION ERR
C
C     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
C
C *Arguments:
C NAME   :IN       Character*(*).
C         Name of the calling routine.  This is used in the output
C         message, if an error is detected.
C LOCIW  :IN       Integer.
C         Location of the first free element in the integer workspace
C         array.
C LENIW  :IN       Integer.
C         Length of the integer workspace array.
C LOCW   :IN       Integer.
C         Location of the first free element in the double precision
C         workspace array.
C LENRW  :IN       Integer.
C         Length of the double precision workspace array.
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for
C                           WORK or IWORK.
C ITER   :OUT      Integer.
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
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW
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
C         Check the Integer workspace situation.
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
     $      'In ' // XERNAM // ', INTEGER work array too short.  ' //
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
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     DOUBLE PRECISION A(NELT)
C
C     CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description",
C         below.  If the SLAP Triad format is used, this format is
C         translated to the SLAP Column format by this routine.
C ISYM   :IN       Integer.
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
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, ITEMP, J
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
*DECK DSILUS
      SUBROUTINE DSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV,
     +   NU, IU, JU, U, NROW, NCOL)
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
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)
C     INTEGER NROW(N), NCOL(N)
C     DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
C
C     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,
C    $    DINV, NU, IU, JU, U, NROW, NCOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NL     :OUT      Integer.
C         Number of non-zeros in the L array.
C IL     :OUT      Integer IL(NL).
C JL     :OUT      Integer JL(NL).
C L      :OUT      Double Precision L(NL).
C         IL, JL, L  contain the unit lower triangular factor of  the
C         incomplete decomposition  of some  matrix stored  in   SLAP
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See
C         "DESCRIPTION", below for more details about the SLAP format.
C NU     :OUT      Integer.
C         Number of non-zeros in the U array.
C IU     :OUT      Integer IU(NU).
C JU     :OUT      Integer JU(NU).
C U      :OUT      Double Precision     U(NU).
C         IU, JU, U contain   the unit upper triangular factor of the
C         incomplete  decomposition    of some matrix  stored in SLAP
C         Column  format.   The Diagonal of ones   *IS*  stored.  See
C         "Description", below  for  more  details  about  the   SLAP
C         format.
C NROW   :WORK     Integer NROW(N).
C         NROW(I) is the number of non-zero elements in the I-th row
C         of L.
C NCOL   :WORK     Integer NCOL(N).
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
      INTEGER ISYM, N, NELT, NL, NU
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), DINV(N), L(NL), U(NU)
      INTEGER IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU),
     +        NCOL(N), NROW(N)
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2,
     +        INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC,
     +        KR
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
*DECK QS2I1D
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)
C***BEGIN PROLOGUE  QS2I1D
C***SUBSIDIARY
C***PURPOSE  Sort an integer array, moving an integer and DP array.
C            This routine sorts the integer array IA and makes the same
C            interchanges in the integer array JA and the double pre-
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
C     along one integer array and one double precision array during
C     the sort by
C     Mark K. Seager
C     Lawrence Livermore National Laboratory
C     November, 1987
C     This routine was adapted from the ISORT routine.
C
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         A - Double Precision array to be carried along.
C         N - Number of values in integer array IA to be sorted.
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
      INTEGER KFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(N)
      INTEGER IA(N), JA(N)
C     .. Local Scalars ..
      DOUBLE PRECISION R, TA, TTA
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
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
