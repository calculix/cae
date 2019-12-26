      SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
      !
      !   -- LAPACK routine (version 2.0) --
      !      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
      !      Courant Institute, Argonne National Lab, and Rice University
      !      March 31, 1993
      !
      !      .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
      !      ..
      !      .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AP( * )
      !      ..
      !
      !   Purpose
      !   =======
      !
      !   DSPTRF computes the factorization of a real symmetric matrix A stored
      !   in packed format using the Bunch-Kaufman diagonal pivoting method:
      !
      !      A = U*D*U**T  or  A = L*D*L**T
      !
      !   where U (or L) is a product of permutation and unit upper (lower)
      !   triangular matrices, and D is symmetric and block diagonal with
      !   1-by-1 and 2-by-2 diagonal blocks.
      !
      !   Arguments
      !   =========
      !
      !   UPLO    (input) CHARACTER*1
      !           = 'U':  Upper triangle of A is stored;
      !           = 'L':  Lower triangle of A is stored.
      !
      !   N       (input) INTEGER
      !           The order of the matrix A.  N >= 0.
      !
      !   AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
      !           On entry, the upper or lower triangle of the symmetric matrix
      !           A, packed columnwise in a linear array.  The j-th column of A
      !           is stored in the array AP as follows:
      !           if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
      !           if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
      !
      !           On exit, the block diagonal matrix D and the multipliers used
      !           to obtain the factor U or L, stored as a packed triangular
      !           matrix overwriting A (see below for further details).
      !
      !   IPIV    (output) INTEGER array, dimension (N)
      !           Details of the interchanges and the block structure of D.
      !           If IPIV(k) > 0, then rows and columns k and IPIV(k) were
      !           interchanged and D(k,k) is a 1-by-1 diagonal block.
      !           If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
      !           columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
      !           is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
      !           IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
      !           interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
      !
      !   INFO    (output) INTEGER
      !           = 0: successful exit
      !           < 0: if INFO = -i, the i-th argument had an illegal value
      !           > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
      !                has been completed, but the block diagonal matrix D is
      !                exactly singular, and division by zero will occur if it
      !                is used to solve a system of equations.
      !
      !   Further Details
      !   ===============
      !
      !   If UPLO = 'U', then A = U*D*U', where
      !      U = P(n)*U(n)* ... *P(k)U(k)* ...,
      !   i.e., U is a product of terms P(k)*U(k), where k decreases from n to
      !   1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
      !   and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
      !   defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
      !   that if the diagonal block D(k) is of order s (s = 1 or 2), then
      !
      !              (   I    v    0   )   k-s
      !      U(k) =  (   0    I    0   )   s
      !              (   0    0    I   )   n-k
      !                 k-s   s   n-k
      !
      !   If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
      !   If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
      !   and A(k,k), and v overwrites A(1:k-2,k-1:k).
      !
      !   If UPLO = 'L', then A = L*D*L', where
      !      L = P(1)*L(1)* ... *P(k)*L(k)* ...,
      !   i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
      !   n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
      !   and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
      !   defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
      !   that if the diagonal block D(k) is of order s (s = 1 or 2), then
      !
      !              (   I    0     0   )  k-1
      !      L(k) =  (   0    I     0   )  s
      !              (   0    v     I   )  n-k-s+1
      !                 k-1   s  n-k-s+1
      !
      !   If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
      !   If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
      !   and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
      !
      !   =====================================================================
      !
      !      .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      !      ..
      !      .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            IMAX, J, JMAX, K, KC, KK, KNC, KP, KPC, KSTEP,&
                         KX, NPP
      DOUBLE PRECISION   ABSAKK, ALPHA, C, COLMAX, R1, R2, ROWMAX, S, T
      !      ..
      !      .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
      !      ..
      !      .. External Subroutines ..
      EXTERNAL           DLAEV2, DROT, DSCAL, DSPR, DSWAP, XERBLA
      !      ..
      !      .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
      !      ..
      !      .. Executable Statements ..
      !
      !      Test the input parameters.
      !
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPTRF', -INFO )
         RETURN
      END IF
      !
      !      Initialize ALPHA for use in choosing pivot block size.
      !
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
      !
      IF( UPPER ) THEN
         !
         !         Factorize A as U*D*U' using the upper triangle of A
         !
         !         K is the main loop index, decreasing from N to 1 in steps of
         !         1 or 2
         !
         K = N
         KC = ( N-1 )*N / 2 + 1
   10    CONTINUE
         KNC = KC
         !
         !         If K < 1, exit from loop
         !
         IF( K.LT.1 )&
            GO TO 70
         KSTEP = 1
         !
         !         Determine rows and columns to be interchanged and whether
         !         a 1-by-1 or 2-by-2 pivot block will be used
         !
         ABSAKK = ABS( AP( KC+K-1 ) )
         !
         !         IMAX is the row-index of the largest off-diagonal element in
         !         column K, and COLMAX is its absolute value
         !
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, AP( KC ), 1 )
            COLMAX = ABS( AP( KC+IMAX-1 ) )
         ELSE
            COLMAX = ZERO
         END IF
         !
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            !
            !            Column K is zero: set INFO and continue
            !
            IF( INFO.EQ.0 )&
               INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
               !
               !               no interchange, use 1-by-1 pivot block
               !
               KP = K
            ELSE
               !
               !               JMAX is the column-index of the largest off-diagonal
               !               element in row IMAX, and ROWMAX is its absolute value
               !
               ROWMAX = ZERO
               JMAX = IMAX
               KX = IMAX*( IMAX+1 ) / 2 + IMAX
               DO 20 J = IMAX + 1, K
                  IF( ABS( AP( KX ) ).GT.ROWMAX ) THEN
                     ROWMAX = ABS( AP( KX ) )
                     JMAX = J
                  END IF
                  KX = KX + J
   20          CONTINUE
               KPC = ( IMAX-1 )*IMAX / 2 + 1
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, AP( KPC ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( AP( KPC+JMAX-1 ) ) )
               END IF
               !
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                  !
                  !                  no interchange, use 1-by-1 pivot block
                  !
                  KP = K
               ELSE IF( ABS( AP( KPC+IMAX-1 ) ).GE.ALPHA*ROWMAX ) THEN
                  !
                  !                  interchange rows and columns K and IMAX, use 1-by-1
                  !                  pivot block
                  !
                  KP = IMAX
               ELSE
                  !
                  !                  interchange rows and columns K-1 and IMAX, use 2-by-2
                  !                  pivot block
                  !
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
            !
            KK = K - KSTEP + 1
            IF( KSTEP.EQ.2 )&
               KNC = KNC - K + 1
            IF( KP.NE.KK ) THEN
               !
               !               Interchange rows and columns KK and KP in the leading
               !               submatrix A(1:k,1:k)
               !
               CALL DSWAP( KP-1, AP( KNC ), 1, AP( KPC ), 1 )
               KX = KPC + KP - 1
               DO 30 J = KP + 1, KK - 1
                  KX = KX + J - 1
                  T = AP( KNC+J-1 )
                  AP( KNC+J-1 ) = AP( KX )
                  AP( KX ) = T
   30          CONTINUE
               T = AP( KNC+KK-1 )
               AP( KNC+KK-1 ) = AP( KPC+KP-1 )
               AP( KPC+KP-1 ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = AP( KC+K-2 )
                  AP( KC+K-2 ) = AP( KC+KP-1 )
                  AP( KC+KP-1 ) = T
               END IF
            END IF
            !
            !            Update the leading submatrix
            !
            IF( KSTEP.EQ.1 ) THEN
               !
               !               1-by-1 pivot block D(k): column k now holds
               !
               !               W(k) = U(k)*D(k)
               !
               !               where U(k) is the k-th column of U
               !
               !               Perform a rank-1 update of A(1:k-1,1:k-1) as
               !
               !               A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
               !
               R1 = ONE / AP( KC+K-1 )
               CALL DSPR( UPLO, K-1, -R1, AP( KC ), 1, AP )
               !
               !               Store U(k) in column k
               !
               CALL DSCAL( K-1, R1, AP( KC ), 1 )
            ELSE
               !
               !               2-by-2 pivot block D(k): columns k and k-1 now hold
               !
               !               ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
               !
               !               where U(k) and U(k-1) are the k-th and (k-1)-th columns
               !               of U
               !
               !               Perform a rank-2 update of A(1:k-2,1:k-2) as
               !
               !               A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
               !                  = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
               !
               !               Convert this to two rank-1 updates by using the eigen-
               !               decomposition of D(k)
               !
               CALL DLAEV2( AP( KC-1 ), AP( KC+K-2 ), AP( KC+K-1 ), R1,&
                            R2, C, S )
               R1 = ONE / R1
               R2 = ONE / R2
               CALL DROT( K-2, AP( KNC ), 1, AP( KC ), 1, C, S )
               CALL DSPR( UPLO, K-2, -R1, AP( KNC ), 1, AP )
               CALL DSPR( UPLO, K-2, -R2, AP( KC ), 1, AP )
               !
               !               Store U(k) and U(k-1) in columns k and k-1
               !
               CALL DSCAL( K-2, R1, AP( KNC ), 1 )
               CALL DSCAL( K-2, R2, AP( KC ), 1 )
               CALL DROT( K-2, AP( KNC ), 1, AP( KC ), 1, C, -S )
            END IF
         END IF
         !
         !         Store details of the interchanges in IPIV
         !
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
         !
         !         Decrease K and return to the start of the main loop
         !
         K = K - KSTEP
         KC = KNC - K
         GO TO 10
      !
      ELSE
         !
         !         Factorize A as L*D*L' using the lower triangle of A
         !
         !         K is the main loop index, increasing from 1 to N in steps of
         !         1 or 2
         !
         K = 1
         KC = 1
         NPP = N*( N+1 ) / 2
   40    CONTINUE
         KNC = KC
         !
         !         If K > N, exit from loop
         !
         IF( K.GT.N )&
            GO TO 70
         KSTEP = 1
         !
         !         Determine rows and columns to be interchanged and whether
         !         a 1-by-1 or 2-by-2 pivot block will be used
         !
         ABSAKK = ABS( AP( KC ) )
         !
         !         IMAX is the row-index of the largest off-diagonal element in
         !         column K, and COLMAX is its absolute value
         !
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, AP( KC+1 ), 1 )
            COLMAX = ABS( AP( KC+IMAX-K ) )
         ELSE
            COLMAX = ZERO
         END IF
         !
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            !
            !            Column K is zero: set INFO and continue
            !
            IF( INFO.EQ.0 )&
               INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
               !
               !               no interchange, use 1-by-1 pivot block
               !
               KP = K
            ELSE
               !
               !               JMAX is the column-index of the largest off-diagonal
               !               element in row IMAX, and ROWMAX is its absolute value
               !
               ROWMAX = ZERO
               KX = KC + IMAX - K
               DO 50 J = K, IMAX - 1
                  IF( ABS( AP( KX ) ).GT.ROWMAX ) THEN
                     ROWMAX = ABS( AP( KX ) )
                     JMAX = J
                  END IF
                  KX = KX + N - J
   50          CONTINUE
               KPC = NPP - ( N-IMAX+1 )*( N-IMAX+2 ) / 2 + 1
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, AP( KPC+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( AP( KPC+JMAX-IMAX ) ) )
               END IF
               !
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                  !
                  !                  no interchange, use 1-by-1 pivot block
                  !
                  KP = K
               ELSE IF( ABS( AP( KPC ) ).GE.ALPHA*ROWMAX ) THEN
                  !
                  !                  interchange rows and columns K and IMAX, use 1-by-1
                  !                  pivot block
                  !
                  KP = IMAX
               ELSE
                  !
                  !                  interchange rows and columns K+1 and IMAX, use 2-by-2
                  !                  pivot block
                  !
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
            !
            KK = K + KSTEP - 1
            IF( KSTEP.EQ.2 )&
               KNC = KNC + N - K + 1
            IF( KP.NE.KK ) THEN
               !
               !               Interchange rows and columns KK and KP in the trailing
               !               submatrix A(k:n,k:n)
               !
               IF( KP.LT.N )&
                  CALL DSWAP( N-KP, AP( KNC+KP-KK+1 ), 1, AP( KPC+1 ),&
                              1 )
               KX = KNC + KP - KK
               DO 60 J = KK + 1, KP - 1
                  KX = KX + N - J + 1
                  T = AP( KNC+J-KK )
                  AP( KNC+J-KK ) = AP( KX )
                  AP( KX ) = T
   60          CONTINUE
               T = AP( KNC )
               AP( KNC ) = AP( KPC )
               AP( KPC ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = AP( KC+1 )
                  AP( KC+1 ) = AP( KC+KP-K )
                  AP( KC+KP-K ) = T
               END IF
            END IF
            !
            !            Update the trailing submatrix
            !
            IF( KSTEP.EQ.1 ) THEN
               !
               !               1-by-1 pivot block D(k): column k now holds
               !
               !               W(k) = L(k)*D(k)
               !
               !               where L(k) is the k-th column of L
               !
               IF( K.LT.N ) THEN
                  !
                  !                  Perform a rank-1 update of A(k+1:n,k+1:n) as
                  !
                  !                  A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
                  !
                  R1 = ONE / AP( KC )
                  CALL DSPR( UPLO, N-K, -R1, AP( KC+1 ), 1,&
                             AP( KC+N-K+1 ) )
                  !
                  !                  Store L(k) in column K
                  !
                  CALL DSCAL( N-K, R1, AP( KC+1 ), 1 )
               END IF
            ELSE
               !
               !               2-by-2 pivot block D(k): columns K and K+1 now hold
               !
               !               ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
               !
               !               where L(k) and L(k+1) are the k-th and (k+1)-th columns
               !               of L
               !
               IF( K.LT.N-1 ) THEN
                  !
                  !                  Perform a rank-2 update of A(k+2:n,k+2:n) as
                  !
                  !                  A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
                  !                     = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
                  !
                  !                  Convert this to two rank-1 updates by using the eigen-
                  !                  decomposition of D(k)
                  !
                  CALL DLAEV2( AP( KC ), AP( KC+1 ), AP( KNC ), R1, R2,&
                               C, S )
                  R1 = ONE / R1
                  R2 = ONE / R2
                  CALL DROT( N-K-1, AP( KC+2 ), 1, AP( KNC+1 ), 1, C,&
                             S )
                  CALL DSPR( UPLO, N-K-1, -R1, AP( KC+2 ), 1,&
                             AP( KNC+N-K ) )
                  CALL DSPR( UPLO, N-K-1, -R2, AP( KNC+1 ), 1,&
                             AP( KNC+N-K ) )
                  !
                  !                  Store L(k) and L(k+1) in columns k and k+1
                  !
                  CALL DSCAL( N-K-1, R1, AP( KC+2 ), 1 )
                  CALL DSCAL( N-K-1, R2, AP( KNC+1 ), 1 )
                  CALL DROT( N-K-1, AP( KC+2 ), 1, AP( KNC+1 ), 1, C,&
                             -S )
               END IF
            END IF
         END IF
         !
         !         Store details of the interchanges in IPIV
         !
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
         !
         !         Increase K and return to the start of the main loop
         !
         K = K + KSTEP
         KC = KNC + N - K + 2
         GO TO 40
      !
      END IF
   !
   70 CONTINUE
      RETURN
      !
      !      End of DSPTRF
      !
      END
      SUBROUTINE DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
      !      .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
      CHARACTER*1        UPLO
      !      .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), X( * )
      !      ..
      !
      !   Purpose
      !   =======
      !
      !   DSPR    performs the symmetric rank 1 operation
      !
      !      A := alpha*x*x' + A,
      !
      !   where alpha is a real scalar, x is an n element vector and A is an
      !   n by n symmetric matrix, supplied in packed form.
      !
      !   Parameters
      !   ==========
      !
      !   UPLO   - CHARACTER*1.
      !            On entry, UPLO specifies whether the upper or lower
      !            triangular part of the matrix A is supplied in the packed
      !            array AP as follows:
      !
      !               UPLO = 'U' or 'u'   The upper triangular part of A is
      !                                   supplied in AP.
      !
      !               UPLO = 'L' or 'l'   The lower triangular part of A is
      !                                   supplied in AP.
      !
      !            Unchanged on exit.
      !
      !   N      - INTEGER.
      !            On entry, N specifies the order of the matrix A.
      !            N must be at least zero.
      !            Unchanged on exit.
      !
      !   ALPHA  - DOUBLE PRECISION.
      !            On entry, ALPHA specifies the scalar alpha.
      !            Unchanged on exit.
      !
      !   X      - DOUBLE PRECISION array of dimension at least
      !            ( 1 + ( n - 1 )*abs( INCX ) ).
      !            Before entry, the incremented array X must contain the n
      !            element vector x.
      !            Unchanged on exit.
      !
      !   INCX   - INTEGER.
      !            On entry, INCX specifies the increment for the elements of
      !            X. INCX must not be zero.
      !            Unchanged on exit.
      !
      !   AP     - DOUBLE PRECISION array of DIMENSION at least
      !            ( ( n*( n + 1 ) )/2 ).
      !            Before entry with  UPLO = 'U' or 'u', the array AP must
      !            contain the upper triangular part of the symmetric matrix
      !            packed sequentially, column by column, so that AP( 1 )
      !            contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
      !            and a( 2, 2 ) respectively, and so on. On exit, the array
      !            AP is overwritten by the upper triangular part of the
      !            updated matrix.
      !            Before entry with UPLO = 'L' or 'l', the array AP must
      !            contain the lower triangular part of the symmetric matrix
      !            packed sequentially, column by column, so that AP( 1 )
      !            contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
      !            and a( 3, 1 ) respectively, and so on. On exit, the array
      !            AP is overwritten by the lower triangular part of the
      !            updated matrix.
      !
      !
      !   Level 2 Blas routine.
      !
      !   -- Written on 22-October-1986.
      !      Jack Dongarra, Argonne National Lab.
      !      Jeremy Du Croz, Nag Central Office.
      !      Sven Hammarling, Nag Central Office.
      !      Richard Hanson, Sandia National Labs.
      !
      !
      !      .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      !      .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
      !      .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
      !      .. External Subroutines ..
      EXTERNAL           XERBLA
      !      ..
      !      .. Executable Statements ..
      !
      !      Test the input parameters.
      !
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.&
               .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSPR  ', INFO )
         RETURN
      END IF
      !
      !      Quick return if possible.
      !
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )&
         RETURN
      !
      !      Set the start point in X if the increment is not unity.
      !
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
      !
      !      Start the operations. In this version the elements of the array AP
      !      are accessed sequentially with one pass through AP.
      !
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
         !
         !         Form  A  when upper triangle is stored in AP.
         !
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  K    = KK
                  DO 10, I = 1, J
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K       = K       + 1
   10             CONTINUE
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, K = KK, KK + J - 1
                     AP( K ) = AP( K ) + X( IX )*TEMP
                     IX      = IX      + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
         !
         !         Form  A  when lower triangle is stored in AP.
         !
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  K    = KK
                  DO 50, I = J, N
                     AP( K ) = AP( K ) + X( I )*TEMP
                     K       = K       + 1
   50             CONTINUE
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, K = KK, KK + N - J
                     AP( K ) = AP( K ) + X( IX )*TEMP
                     IX      = IX      + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
      !
      RETURN
      !
      !      End of DSPR  .
      !
      END
