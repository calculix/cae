      SUBROUTINE HYBSVD(NA, NU, NV, NZ, NB, M, N, A, W, MATU, U, MATV,&
       V, Z, B, IRHS, IERR, RV1)
      INTEGER NA, NU, NV, NZ, M, N, IRHS, IERR, MIN0
      REAL*8 A(NA,1), W(1), U(NU,1), V(NV,1), Z(NZ,1), B(NB,IRHS)
      REAL*8 RV1(1)
      LOGICAL MATU, MATV
      !
      !      THIS ROUTINE IS A MODIFICATION OF THE GOLUB-REINSCH PROCEDURE (1)
      !                                                            T
      !      FOR COMPUTING THE SINGULAR VALUE DECOMPOSITION A = UWV  OF A
      !      REAL M BY N RECTANGULAR MATRIX. U IS M BY MIN(M,N) CONTAINING
      !      THE LEFT SINGULAR VECTORS, W IS A MIN(M,N) BY MIN(M,N) DIAGONAL
      !      MATRIX CONTAINING THE SINGULAR VALUES, AND V IS N BY MIN(M,N)
      !      CONTAINING THE RIGHT SINGULAR VECTORS.
      !
      !      THE ALGORITHM IMPLEMENTED IN THIS
      !      ROUTINE HAS A HYBRID NATURE.  WHEN M IS APPROXIMATELY EQUAL TO N,
      !      THE GOLUB-REINSCH ALGORITHM IS USED, BUT WHEN EITHER OF THE RATIOS
      !      M/N OR N/M IS GREATER THAN ABOUT 2,
      !      A MODIFIED VERSION OF THE GOLUB-REINSCH
      !      ALGORITHM IS USED.  THIS MODIFIED ALGORITHM FIRST TRANSFORMS A
      !                                                                 T
      !      INTO UPPER TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATIONS L
      !      AND THEN USES THE GOLUB-REINSCH ALGORITHM TO FIND THE SINGULAR
      !      VALUE DECOMPOSITION OF THE RESULTING UPPER TRIANGULAR MATRIX R.
      !      WHEN U IS NEEDED EXPLICITLY IN THE CASE M.GE.N (OR V IN THE CASE
      !      M.LT.N), AN EXTRA ARRAY Z (OF SIZE AT LEAST
      !      MIN(M,N)**2) IS NEEDED, BUT OTHERWISE Z IS NOT REFERENCED
      !      AND NO EXTRA STORAGE IS REQUIRED.  THIS HYBRID METHOD
      !      SHOULD BE MORE EFFICIENT THAN THE GOLUB-REINSCH ALGORITHM WHEN
      !      M/N OR N/M IS LARGE.  FOR DETAILS, SEE (2).
      !
      !      WHEN M .GE. N,
      !      HYBSVD CAN ALSO BE USED TO COMPUTE THE MINIMAL LENGTH LEAST
      !      SQUARES SOLUTION TO THE OVERDETERMINED LINEAR SYSTEM A*X=B.
      !      IF M .LT. N (I.E. FOR UNDERDETERMINED SYSTEMS), THE RHS B
      !      IS NOT PROCESSED.
      !
      !      NOTICE THAT THE SINGULAR VALUE DECOMPOSITION OF A MATRIX
      !      IS UNIQUE ONLY UP TO THE SIGN OF THE CORRESPONDING COLUMNS
      !      OF U AND V.
      !
      !      THIS ROUTINE HAS BEEN CHECKED BY THE PFORT VERIFIER (3) FOR
      !      ADHERENCE TO A LARGE, CAREFULLY DEFINED, PORTABLE SUBSET OF
      !      AMERICAN NATIONAL STANDARD FORTRAN CALLED PFORT.
      !
      !      REFERENCES:
      !
      !      (1) GOLUB,G.H. AND REINSCH,C. (1970) 'SINGULAR VALUE
      !          DECOMPOSITION AND LEAST SQUARES SOLUTIONS,'
      !          NUMER. MATH. 14,403-420, 1970.
      !
      !      (2) CHAN,T.F. (1982) 'AN IMPROVED ALGORITHM FOR COMPUTING
      !          THE SINGULAR VALUE DECOMPOSITION,' ACM TOMS, VOL.8,
      !          NO. 1, MARCH, 1982.
      !
      !      (3) RYDER,B.G. (1974) 'THE PFORT VERIFIER,' SOFTWARE -
      !          PRACTICE AND EXPERIENCE, VOL.4, 359-377, 1974.
      !
      !      ON INPUT:
      !
      !         NA MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY PARAMETER A AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT.  NOTE THAT NA MUST BE AT LEAST
      !           AS LARGE AS M.
      !
      !         NU MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY U AS DECLARED IN THE CALLING PROGRAM DIMENSION
      !           STATEMENT. NU MUST BE AT LEAST AS LARGE AS M.
      !
      !         NV MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY PARAMETER V AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT. NV MUST BE AT LEAST AS LARGE AS N.
      !
      !         NZ MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY PARAMETER Z AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT.  NOTE THAT NZ MUST BE AT LEAST
      !           AS LARGE AS MIN(M,N).
      !
      !         NB MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY PARAMETER B AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT. NB MUST BE AT LEAST AS LARGE AS M.
      !
      !         M IS THE NUMBER OF ROWS OF A (AND U).
      !
      !         N IS THE NUMBER OF COLUMNS OF A (AND NUMBER OF ROWS OF V).
      !
      !         A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED.
      !
      !         B CONTAINS THE IRHS RIGHT-HAND-SIDES OF THE OVERDETERMINED
      !          LINEAR SYSTEM A*X=B. IF IRHS .GT. 0 AND M .GE. N,
      !          THEN ON OUTPUT, THE FIRST N COMPONENTS OF THESE IRHS COLUMNS
      !                        T
      !          WILL CONTAIN U B. THUS, TO COMPUTE THE MINIMAL LENGTH LEAST
      !                                                +
      !          SQUARES SOLUTION, ONE MUST COMPUTE V*W  TIMES THE COLUMNS OF
      !                    +                        +
      !          B, WHERE W  IS A DIAGONAL MATRIX, W (I)=0 IF W(I) IS
      !          NEGLIGIBLE, OTHERWISE IS 1/W(I). IF IRHS=0 OR M.LT.N,
      !          B IS NOT REFERENCED.
      !
      !         IRHS IS THE NUMBER OF RIGHT-HAND-SIDES OF THE OVERDETERMINED
      !          SYSTEM A*X=B. IRHS SHOULD BE SET TO ZERO IF ONLY THE SINGULAR
      !          VALUE DECOMPOSITION OF A IS DESIRED.
      !
      !         MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
      !           DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
      !
      !         MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
      !           DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
      !
      !         WHEN HYBSVD IS USED TO COMPUTE THE MINIMAL LENGTH LEAST
      !         SQUARES SOLUTION TO AN OVERDETERMINED SYSTEM, MATU SHOULD
      !         BE SET TO .FALSE. , AND MATV SHOULD BE SET TO .TRUE.  .
      !
      !      ON OUTPUT:
      !
      !         A IS UNALTERED (UNLESS OVERWRITTEN BY U OR V).
      !
      !         W CONTAINS THE (NON-NEGATIVE) SINGULAR VALUES OF A (THE
      !           DIAGONAL ELEMENTS OF W).  THEY ARE SORTED IN DESCENDING
      !           ORDER.  IF AN ERROR EXIT IS MADE, THE SINGULAR VALUES
      !           SHOULD BE CORRECT AND SORTED FOR INDICES IERR+1,...,MIN(M,N).
      !
      !         U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
      !           DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  IF MATU IS
      !           FALSE, THEN U IS EITHER USED AS A TEMPORARY STORAGE (IF
      !           M .GE. N) OR NOT REFERENCED (IF M .LT. N).
      !           U MAY COINCIDE WITH A IN THE CALLING SEQUENCE.
      !           IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
      !           TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT.
      !
      !         V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
      !           MATV HAS BEEN SET TO .TRUE.  IF MATV IS
      !           FALSE, THEN V IS EITHER USED AS A TEMPORARY STORAGE (IF
      !           M .LT. N) OR NOT REFERENCED (IF M .GE. N).
      !           IF M .GE. N, V MAY ALSO COINCIDE WITH A.  IF AN ERROR
      !           EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO INDICES OF
      !           CORRECT SINGULAR VALUES SHOULD BE CORRECT.
      !
      !         Z CONTAINS THE MATRIX X IN THE SINGULAR VALUE DECOMPOSITION
      !                   T
      !           OF R=XSY,  IF THE MODIFIED ALGORITHM IS USED. IF THE
      !           GOLUB-REINSCH PROCEDURE IS USED, THEN IT IS NOT REFERENCED.
      !           IF MATU HAS BEEN SET TO .FALSE. IN THE CASE M.GE.N (OR
      !           MATV SET TO .FALSE. IN THE CASE M.LT.N), THEN Z IS NOT
      !           REFERENCED AND NO EXTRA STORAGE IS REQUIRED.
      !
      !         IERR IS SET TO
      !           ZERO       FOR NORMAL RETURN,
      !           K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
      !                      DETERMINED AFTER 30 ITERATIONS.
      !           -1         IF IRHS .LT. 0 .
      !           -2         IF M .LT. 1 .OR. N .LT. 1
      !           -3         IF NA .LT. M .OR. NU .LT. M .OR. NB .LT. M.
      !           -4         IF NV .LT. N .
      !           -5         IF NZ .LT. MIN(M,N).
      !
      !         RV1 IS A TEMPORARY STORAGE ARRAY OF LENGTH AT LEAST MIN(M,N).
      !
      !      PROGRAMMED BY : TONY CHAN
      !                      BOX 2158, YALE STATION,
      !                      COMPUTER SCIENCE DEPT, YALE UNIV.,
      !                      NEW HAVEN, CT 06520.
      !      LAST MODIFIED : JANUARY, 1982.
      !
      !      HYBSVD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
      !        INTERNAL  GRSVD, MGNSVD, SRELPR
      !        FORTRAN   MIN0,DABS,DSQRT,DFLOAT,DSIGN,DMAX1
      !        BLAS      DSWAP
      !
      !      -----------------------------------------------------------------
      !      ERROR CHECK.
      !
      IERR = 0
      IF (IRHS.GE.0) GO TO 10
      IERR = -1
      RETURN
   10 IF (M.GE.1 .AND. N.GE.1) GO TO 20
      IERR = -2
      RETURN
   20 IF (NA.GE.M .AND. NU.GE.M .AND. NB.GE.M) GO TO 30
      IERR = -3
      RETURN
   30 IF (NV.GE.N) GO TO 40
      IERR = -4
      RETURN
   40 IF (NZ.GE.MIN0(M,N)) GO TO 50
      IERR = -5
      RETURN
   50 CONTINUE
      !
      !      FIRST COPIES A INTO EITHER U OR V ACCORDING TO WHETHER
      !      M .GE. N OR M .LT. N, AND THEN CALLS SUBROUTINE MGNSVD
      !      WHICH ASSUMES THAT NUMBER OF ROWS .GE. NUMBER OF COLUMNS.
      !
      IF (M.LT.N) GO TO 80
      !
      !        M .GE. N  CASE.
      !
      DO 70 I=1,M
        DO 60 J=1,N
          U(I,J) = A(I,J)
   60   CONTINUE
   70 CONTINUE
      !
      CALL MGNSVD(NU, NV, NZ, NB, M, N, W, MATU, U, MATV, V, Z, B,&
       IRHS, IERR, RV1)
      RETURN
   !
   80 CONTINUE
      !                               T
      !        M .LT. N CASE. COPIES A  INTO V.
      !
      DO 100 I=1,M
        DO 90 J=1,N
          V(J,I) = A(I,J)
   90   CONTINUE
  100 CONTINUE
      CALL MGNSVD(NV, NU, NZ, NB, N, M, W, MATV, V, MATU, U, Z, B, 0,&
       IERR, RV1)
      RETURN
      END
      !                                                                        MGN   10
      SUBROUTINE MGNSVD(NU, NV, NZ, NB, M, N, W, MATU, U, MATV, V, Z,&
       B, IRHS, IERR, RV1)
      !
      !      THE DESCRIPTION OF SUBROUTINE MGNSVD IS ALMOST IDENTICAL
      !      TO THAT FOR SUBROUTINE HYBSVD ABOVE, WITH THE EXCEPTION
      !      THAT MGNSVD ASSUMES M .GE. N.
      !      IT ALSO ASSUMES THAT A COPY OF THE MATRIX A IS IN THE ARRAY U.
      !
      INTEGER NU, NV, NZ, M, N, IRHS, IERR, IP1, I, J, K, IM1, IBACK
      REAL*8 W(1), U(NU,1), V(NV,1), Z(NZ,1), B(NB,IRHS), RV1(1)
      REAL*8 XOVRPT, C, R, G, SCALE, DSIGN, DABS, DSQRT, F, S, H
      REAL*8 DFLOAT
      LOGICAL MATU, MATV
      !
      !      SET VALUE FOR C. THE VALUE FOR C DEPENDS ON THE RELATIVE
      !      EFFICIENCY OF FLOATING POINT MULTIPLICATIONS, FLOATING POINT
      !      ADDITIONS AND TWO-DIMENSIONAL ARRAY INDEXINGS ON THE
      !      COMPUTER WHERE THIS SUBROUTINE IS TO BE RUN.  C SHOULD
      !      USUALLY BE BETWEEN 2 AND 4.  FOR DETAILS ON CHOOSING C, SEE
      !      (2).  THE ALGORITHM IS NOT SENSITIVE TO THE VALUE OF C
      !      ACTUALLY USED AS LONG AS C IS BETWEEN 2 AND 4.
      !
      C = 4.d0
      !
      !      DETERMINE CROSS-OVER POINT
      !
      IF (MATU .AND. MATV) XOVRPT = (C+7.d0/3.d0)/C
      IF (MATU .AND. .NOT.MATV) XOVRPT = (C+7.d0/3.d0)/C
      IF (.NOT.MATU .AND. MATV) XOVRPT = 5.d0/3.d0
      IF (.NOT.MATU .AND. .NOT.MATV) XOVRPT = 5.d0/3.d0
      !
      !      DETERMINE WHETHER TO USE GOLUB-REINSCH OR THE MODIFIED
      !      ALGORITHM.
      !
      R = DFLOAT(M)/DFLOAT(N)
      IF (R.GE.XOVRPT) GO TO 10
      !
      !      USE GOLUB-REINSCH PROCEDURE
      !
      CALL GRSVD(NU, NV, NB, M, N, W, MATU, U, MATV, V, B, IRHS, IERR,&
       RV1)
      GO TO 330
   !
   !      USE MODIFIED ALGORITHM
   !
   10 CONTINUE
      !
      !      TRIANGULARIZE U BY HOUSEHOLDER TRANSFORMATIONS, USING
      !      W AND RV1 AS TEMPORARY STORAGE.
      !
      DO 110 I=1,N
        G = 0.d0
        S = 0.d0
        SCALE = 0.d0
        !
        !          PERFORM SCALING OF COLUMNS TO AVOID UNNECSSARY OVERFLOW
        !          OR UNDERFLOW
        !
        DO 20 K=I,M
          SCALE = SCALE + DABS(U(K,I))
   20   CONTINUE
        IF (SCALE.EQ.0.d0) GO TO 110
        DO 30 K=I,M
          U(K,I) = U(K,I)/SCALE
          S = S + U(K,I)*U(K,I)
   30   CONTINUE
        !
        !          THE VECTOR E OF THE HOUSEHOLDER TRANSFORMATION I + EE'/H
        !          WILL BE STORED IN COLUMN I OF U. THE TRANSFORMED ELEMENT
        !          U(I,I) WILL BE STORED IN W(I) AND THE SCALAR H IN
        !          RV1(I).
        !
        F = U(I,I)
        G = -DSIGN(DSQRT(S),F)
        H = F*G - S
        U(I,I) = F - G
        RV1(I) = H
        W(I) = SCALE*G
        !
        IF (I.EQ.N) GO TO 70
        !
        !          APPLY TRANSFORMATIONS TO REMAINING COLUMNS OF A
        !
        IP1 = I + 1
        DO 60 J=IP1,N
          S = 0.d0
          DO 40 K=I,M
            S = S + U(K,I)*U(K,J)
   40     CONTINUE
          F = S/H
          DO 50 K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
   50     CONTINUE
   60   CONTINUE
   !
   !          APPLY TRANSFORMATIONS TO COLUMNS OF B IF IRHS .GT. 0
   !
   70   IF (IRHS.EQ.0) GO TO 110
        DO 100 J=1,IRHS
          S = 0.d0
          DO 80 K=I,M
            S = S + U(K,I)*B(K,J)
   80     CONTINUE
          F = S/H
          DO 90 K=I,M
            B(K,J) = B(K,J) + F*U(K,I)
   90     CONTINUE
  100   CONTINUE
  110 CONTINUE
      !
      !      COPY R INTO Z IF MATU = .TRUE.
      !
      IF (.NOT.MATU) GO TO 290
      DO 130 I=1,N
        DO 120 J=I,N
          Z(J,I) = 0.d0
          Z(I,J) = U(I,J)
  120   CONTINUE
        Z(I,I) = W(I)
  130 CONTINUE
      !
      !      ACCUMULATE HOUSEHOLDER TRANSFORMATIONS IN U
      !
      DO 240 IBACK=1,N
        I = N - IBACK + 1
        IP1 = I + 1
        G = W(I)
        H = RV1(I)
        IF (I.EQ.N) GO TO 150
        !
        DO 140 J=IP1,N
          U(I,J) = 0.d0
  140   CONTINUE
  !
  150   IF (H.EQ.0.d0) GO TO 210
        IF (I.EQ.N) GO TO 190
        !
        DO 180 J=IP1,N
          S = 0.d0
          DO 160 K=IP1,M
            S = S + U(K,I)*U(K,J)
  160     CONTINUE
          F = S/H
          DO 170 K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
  170     CONTINUE
  180   CONTINUE
  !
  190   S = U(I,I)/H
        DO 200 J=I,M
          U(J,I) = U(J,I)*S
  200   CONTINUE
        GO TO 230
  !
  210   DO 220 J=I,M
          U(J,I) = 0.d0
  220   CONTINUE
  230   U(I,I) = U(I,I) + 1.d0
  240 CONTINUE
      !
      !      COMPUTE SVD OF R (WHICH IS STORED IN Z)
      !
      CALL GRSVD(NZ, NV, NB, N, N, W, MATU, Z, MATV, V, B, IRHS, IERR,&
       RV1)
      !
      !                                       T
      !      FORM L*X TO OBTAIN U (WHERE R=XWY ). X IS RETURNED IN Z
      !      BY GRSVD. THE MATRIX MULTIPLY IS DONE ONE ROW AT A TIME,
      !      USING RV1 AS SCRATCH SPACE.
      !
      DO 280 I=1,M
        DO 260 J=1,N
          S = 0.d0
          DO 250 K=1,N
            S = S + U(I,K)*Z(K,J)
  250     CONTINUE
          RV1(J) = S
  260   CONTINUE
        DO 270 J=1,N
          U(I,J) = RV1(J)
  270   CONTINUE
  280 CONTINUE
      GO TO 330
  !
  !      FORM R IN U BY ZEROING THE LOWER TRIANGULAR PART OF R IN U
  !
  290 IF (N.EQ.1) GO TO 320
      DO 310 I=2,N
        IM1 = I - 1
        DO 300 J=1,IM1
          U(I,J) = 0.d0
  300   CONTINUE
        U(I,I) = W(I)
  310 CONTINUE
  320 U(1,1) = W(1)
      !
      CALL GRSVD(NU, NV, NB, N, N, W, MATU, U, MATV, V, B, IRHS, IERR,&
       RV1)
  330 CONTINUE
      IERRP1 = IERR + 1
      IF (IERR.LT.0 .OR. N.LE.1 .OR. IERRP1.EQ.N) RETURN
      !
      !      SORT SINGULAR VALUES AND EXCHANGE COLUMNS OF U AND V ACCORDINGLY.
      !      SELECTION SORT MINIMIZES SWAPPING OF U AND V.
      !
      NM1 = N - 1
      DO 360 I=IERRP1,NM1
        ! ...    FIND INDEX OF MAXIMUM SINGULAR VALUE
        ID = I
        IP1 = I + 1
        DO 340 J=IP1,N
          IF (W(J).GT.W(ID)) ID = J
  340   CONTINUE
        IF (ID.EQ.I) GO TO 360
        ! ...    SWAP SINGULAR VALUES AND VECTORS
        T = W(I)
        W(I) = W(ID)
        W(ID) = T
        IF (MATV) CALL DSWAP(N, V(1,I), 1, V(1,ID), 1)
        IF (MATU) CALL DSWAP(M, U(1,I), 1, U(1,ID), 1)
        IF (IRHS.LT.1) GO TO 360
        DO 350 KRHS=1,IRHS
          T = B(I,KRHS)
          B(I,KRHS) = B(ID,KRHS)
          B(ID,KRHS) = T
  350   CONTINUE
  360 CONTINUE
      RETURN
      !      ************** LAST CARD OF HYBSVD *****************
      END
      SUBROUTINE GRSVD(NU, NV, NB, M, N, W, MATU, U, MATV, V, B, IRHS,&
       IERR, RV1)
      !
      INTEGER I, J, K, L, M, N, II, I1, KK, K1, LL, L1, MN, NU, NV, NB,&
       ITS, IERR, IRHS
      REAL*8 W(1), U(NU,1), V(NV,1), B(NB,IRHS), RV1(1)
      REAL*8 C, F, G, H, S, X, Y, Z, EPS, SCALE, SRELPR, DUMMY
      REAL*8 DSQRT, DMAX1, DABS, DSIGN
      LOGICAL MATU, MATV
      !
      !      THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
      !      NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
      !      HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
      !
      !      THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
      !           T
      !      A=USV  OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
      !      BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
      !      GRSVD ASSUMES THAT A COPY OF THE MATRIX A IS IN THE ARRAY U. IT
      !      ALSO ASSUMES M .GE. N.  IF M .LT. N, THEN COMPUTE THE SINGULAR
      !                              T       T    T             T
      !      VALUE DECOMPOSITION OF A .  IF A =UWV  , THEN A=VWU  .
      !
      !      GRSVD CAN ALSO BE USED TO COMPUTE THE MINIMAL LENGTH LEAST SQUARES
      !      SOLUTION TO THE OVERDETERMINED LINEAR SYSTEM A*X=B.
      !
      !      ON INPUT-
      !
      !         NU MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
      !           ARRAY PARAMETERS U AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT.  NOTE THAT NU MUST BE AT LEAST
      !           AS LARGE AS M,
      !
      !         NV MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
      !           ARRAY PARAMETER V AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT.  NV MUST BE AT LEAST AS LARGE AS N,
      !
      !         NB MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
      !           ARRAY PARAMETERS B AS DECLARED IN THE CALLING PROGRAM
      !           DIMENSION STATEMENT.  NOTE THAT NB MUST BE AT LEAST
      !           AS LARGE AS M,
      !
      !         M IS THE NUMBER OF ROWS OF A (AND U),
      !
      !         N IS THE NUMBER OF COLUMNS OF A (AND U) AND THE ORDER OF V,
      !
      !         A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED,
      !
      !         B CONTAINS THE IRHS RIGHT-HAND-SIDES OF THE OVERDETERMINED
      !           LINEAR SYSTEM A*X=B.  IF IRHS .GT. 0,  THEN ON OUTPUT,
      !           THE FIRST N COMPONENTS OF THESE IRHS COLUMNS OF B
      !                         T
      !           WILL CONTAIN U B.  THUS, TO COMPUTE THE MINIMAL LENGTH LEAST
      !                                                 +
      !           SQUARES SOLUTION, ONE MUST COMPUTE V*W  TIMES THE COLUMNS OF
      !                     +                        +
      !           B, WHERE W  IS A DIAGONAL MATRIX, W (I)=0 IF W(I) IS
      !           NEGLIGIBLE, OTHERWISE IS 1/W(I).  IF IRHS=0, B MAY COINCIDE
      !           WITH A OR U AND WILL NOT BE REFERENCED,
      !
      !         IRHS IS THE NUMBER OF RIGHT-HAND-SIDES OF THE OVERDETERMINED
      !           SYSTEM A*X=B.  IRHS SHOULD BE SET TO ZERO IF ONLY THE SINGULA
      !           VALUE DECOMPOSITION OF A IS DESIRED,
      !
      !         MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
      !           DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE,
      !
      !         MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
      !           DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
      !
      !      ON OUTPUT-
      !
      !         W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
      !           DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
      !           ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
      !           FOR INDICES IERR+1,IERR+2,...,N,
      !
      !         U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
      !           DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  OTHERWISE
      !           U IS USED AS A TEMPORARY ARRAY.
      !           IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
      !           TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
      !
      !         V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
      !           MATV HAS BEEN SET TO .TRUE.  OTHERWISE V IS NOT REFERENCED.
      !           IF AN ERROR EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO
      !           INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
      !
      !         IERR IS SET TO
      !           ZERO       FOR NORMAL RETURN,
      !           K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
      !                      DETERMINED AFTER 30 ITERATIONS,
      !           -1         IF IRHS .LT. 0 ,
      !           -2         IF M .LT. N ,
      !           -3         IF NU .LT. M .OR. NB .LT. M,
      !           -4         IF NV .LT. N .
      !
      !         RV1 IS A TEMPORARY STORAGE ARRAY.
      !
      !         THIS SUBROUTINE HAS BEEN CHECKED BY THE PFORT VERIFIER
      !         (RYDER, B.G. 'THE PFORT VERIFIER', SOFTWARE - PRACTICE AND
      !         EXPERIENCE, VOL.4, 359-377, 1974) FOR ADHERENCE TO A LARGE,
      !         CAREFULLY DEFINED, PORTABLE SUBSET OF AMERICAN NATIONAL STANDAR
      !         FORTRAN CALLED PFORT.
      !
      !         ORIGINAL VERSION OF THIS CODE IS SUBROUTINE SVD IN RELEASE 2 OF
      !         EISPACK.
      !
      !         MODIFIED BY TONY F. CHAN,
      !                     COMP. SCI. DEPT, YALE UNIV.,
      !                     BOX 2158, YALE STATION,
      !                     CT 06520
      !         LAST MODIFIED : JANUARY, 1982.
      !
      !      ------------------------------------------------------------------
      !
      !      ********** SRELPR IS A MACHINE-DEPENDENT FUNCTION SPECIFYING
      !                 THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
      !
      !                 **********
      !
      IERR = 0
      IF (IRHS.GE.0) GO TO 10
      IERR = -1
      RETURN
   10 IF (M.GE.N) GO TO 20
      IERR = -2
      RETURN
   20 IF (NU.GE.M .AND. NB.GE.M) GO TO 30
      IERR = -3
      RETURN
   30 IF (NV.GE.N) GO TO 40
      IERR = -4
      RETURN
   40 CONTINUE
      !
      !      ********** HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM **********
      G = 0.d0
      SCALE = 0.d0
      X = 0.d0
      !
      DO 260 I=1,N
        L = I + 1
        RV1(I) = SCALE*G
        G = 0.d0
        S = 0.d0
        SCALE = 0.d0
        !
        !      COMPUTE LEFT TRANSFORMATIONS THAT ZERO THE SUBDIAGONAL ELEMENTS
        !      OF THE I-TH COLUMN.
        !
        DO 50 K=I,M
          SCALE = SCALE + DABS(U(K,I))
   50   CONTINUE
        !
        IF (SCALE.EQ.0.d0) GO TO 160
        !
        DO 60 K=I,M
          U(K,I) = U(K,I)/SCALE
          S = S + U(K,I)**2
   60   CONTINUE
        !
        F = U(I,I)
        G = -DSIGN(DSQRT(S),F)
        H = F*G - S
        U(I,I) = F - G
        IF (I.EQ.N) GO TO 100
        !
        !      APPLY LEFT TRANSFORMATIONS TO REMAINING COLUMNS OF A.
        !
        DO 90 J=L,N
          S = 0.d0
          !
          DO 70 K=I,M
            S = S + U(K,I)*U(K,J)
   70     CONTINUE
          !
          F = S/H
          !
          DO 80 K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
   80     CONTINUE
   90   CONTINUE
  !
  !       APPLY LEFT TRANSFORMATIONS TO THE COLUMNS OF B IF IRHS .GT. 0
  !
  100   IF (IRHS.EQ.0) GO TO 140
        DO 130 J=1,IRHS
          S = 0.d0
          DO 110 K=I,M
            S = S + U(K,I)*B(K,J)
  110     CONTINUE
          F = S/H
          DO 120 K=I,M
            B(K,J) = B(K,J) + F*U(K,I)
  120     CONTINUE
  130   CONTINUE
  !
  !      COMPUTE RIGHT TRANSFORMATIONS.
  !
  140   DO 150 K=I,M
          U(K,I) = SCALE*U(K,I)
  150   CONTINUE
  !
  160   W(I) = SCALE*G
        G = 0.d0
        S = 0.d0
        SCALE = 0.d0
        IF (I.GT.M .OR. I.EQ.N) GO TO 250
        !
        DO 170 K=L,N
          SCALE = SCALE + DABS(U(I,K))
  170   CONTINUE
        !
        IF (SCALE.EQ.0.d0) GO TO 250
        !
        DO 180 K=L,N
          U(I,K) = U(I,K)/SCALE
          S = S + U(I,K)**2
  180   CONTINUE
        !
        F = U(I,L)
        G = -DSIGN(DSQRT(S),F)
        H = F*G - S
        U(I,L) = F - G
        !
        DO 190 K=L,N
          RV1(K) = U(I,K)/H
  190   CONTINUE
        !
        IF (I.EQ.M) GO TO 230
        !
        DO 220 J=L,M
          S = 0.d0
          !
          DO 200 K=L,N
            S = S + U(J,K)*U(I,K)
  200     CONTINUE
          !
          DO 210 K=L,N
            U(J,K) = U(J,K) + S*RV1(K)
  210     CONTINUE
  220   CONTINUE
  !
  230   DO 240 K=L,N
          U(I,K) = SCALE*U(I,K)
  240   CONTINUE
  !
  250   X = DMAX1(X,DABS(W(I))+DABS(RV1(I)))
  260 CONTINUE
      !      ********** ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS **********
      IF (.NOT.MATV) GO TO 350
      !      ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 340 II=1,N
        I = N + 1 - II
        IF (I.EQ.N) GO TO 330
        IF (G.EQ.0.d0) GO TO 310
        !
        DO 270 J=L,N
          !      ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
          V(J,I) = (U(I,J)/U(I,L))/G
  270   CONTINUE
        !
        DO 300 J=L,N
          S = 0.d0
          !
          DO 280 K=L,N
            S = S + U(I,K)*V(K,J)
  280     CONTINUE
          !
          DO 290 K=L,N
            V(K,J) = V(K,J) + S*V(K,I)
  290     CONTINUE
  300   CONTINUE
  !
  310   DO 320 J=L,N
          V(I,J) = 0.d0
          V(J,I) = 0.d0
  320   CONTINUE
  !
  330   V(I,I) = 1.d0
        G = RV1(I)
        L = I
  340 CONTINUE
  !      ********** ACCUMULATION OF LEFT-HAND TRANSFORMATIONS **********
  350 IF (.NOT.MATU) GO TO 470
      !      **********FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- **********
      MN = N
      IF (M.LT.N) MN = M
      !
      DO 460 II=1,MN
        I = MN + 1 - II
        L = I + 1
        G = W(I)
        IF (I.EQ.N) GO TO 370
        !
        DO 360 J=L,N
          U(I,J) = 0.d0
  360   CONTINUE
  !
  370   IF (G.EQ.0.d0) GO TO 430
        IF (I.EQ.MN) GO TO 410
        !
        DO 400 J=L,N
          S = 0.d0
          !
          DO 380 K=L,M
            S = S + U(K,I)*U(K,J)
  380     CONTINUE
          !      ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
          F = (S/U(I,I))/G
          !
          DO 390 K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
  390     CONTINUE
  400   CONTINUE
  !
  410   DO 420 J=I,M
          U(J,I) = U(J,I)/G
  420   CONTINUE
        !
        GO TO 450
  !
  430   DO 440 J=I,M
          U(J,I) = 0.d0
  440   CONTINUE
  !
  450   U(I,I) = U(I,I) + 1.d0
  460 CONTINUE
  !      ********** DIAGONALIZATION OF THE BIDIAGONAL FORM **********
  470 EPS = SRELPR(DUMMY)*X
      !      ********** FOR K=N STEP -1 UNTIL 1 DO -- **********
      DO 650 KK=1,N
        K1 = N - KK
        K = K1 + 1
        ITS = 0
  !      ********** TEST FOR SPLITTING.
  !                 FOR L=K STEP -1 UNTIL 1 DO -- **********
  480   DO 490 LL=1,K
          L1 = K - LL
          L = L1 + 1
          IF (DABS(RV1(L)).LE.EPS) GO TO 550
          !      ********** RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
          !                 THROUGH THE BOTTOM OF THE LOOP **********
          IF (DABS(W(L1)).LE.EPS) GO TO 500
  490   CONTINUE
  !      ********** CANCELLATION OF RV1(L) IF L GREATER THAN 1 **********
  500   C = 0.d0
        S = 1.d0
        !
        DO 540 I=L,K
          F = S*RV1(I)
          RV1(I) = C*RV1(I)
          IF (DABS(F).LE.EPS) GO TO 550
          G = W(I)
          H = DSQRT(F*F+G*G)
          W(I) = H
          C = G/H
          S = -F/H
          !
          !      APPLY LEFT TRANSFORMATIONS TO B IF IRHS .GT. 0
          !
          IF (IRHS.EQ.0) GO TO 520
          DO 510 J=1,IRHS
            Y = B(L1,J)
            Z = B(I,J)
            B(L1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
  510     CONTINUE
  520     CONTINUE
          !
          IF (.NOT.MATU) GO TO 540
          !
          DO 530 J=1,M
            Y = U(J,L1)
            Z = U(J,I)
            U(J,L1) = Y*C + Z*S
            U(J,I) = -Y*S + Z*C
  530     CONTINUE
  !
  540   CONTINUE
  !      ********** TEST FOR CONVERGENCE **********
  550   Z = W(K)
        IF (L.EQ.K) GO TO 630
        !      ********** SHIFT FROM BOTTOM 2 BY 2 MINOR **********
        IF (ITS.EQ.30) GO TO 660
        ITS = ITS + 1
        X = W(L)
        Y = W(K1)
        G = RV1(K1)
        H = RV1(K)
        F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.d0*H*Y)
        G = DSQRT(F*F+1.0)
        F = ((X-Z)*(X+Z)+H*(Y/(F+DSIGN(G,F))-H))/X
        !      ********** NEXT QR TRANSFORMATION **********
        C = 1.0
        S = 1.0
        !
        DO 620 I1=L,K1
          I = I1 + 1
          G = RV1(I)
          Y = W(I)
          H = S*G
          G = C*G
          Z = DSQRT(F*F+H*H)
          RV1(I1) = Z
          C = F/Z
          S = H/Z
          F = X*C + G*S
          G = -X*S + G*C
          H = Y*S
          Y = Y*C
          IF (.NOT.MATV) GO TO 570
          !
          DO 560 J=1,N
            X = V(J,I1)
            Z = V(J,I)
            V(J,I1) = X*C + Z*S
            V(J,I) = -X*S + Z*C
  560     CONTINUE
  !
  570     Z = DSQRT(F*F+H*H)
          W(I1) = Z
          !      ********** ROTATION CAN BE ARBITRARY IF Z IS ZERO **********
          IF (Z.EQ.0.d0) GO TO 580
          C = F/Z
          S = H/Z
  580     F = C*G + S*Y
          X = -S*G + C*Y
          !
          !      APPLY LEFT TRANSFORMATIONS TO B IF IRHS .GT. 0
          !
          IF (IRHS.EQ.0) GO TO 600
          DO 590 J=1,IRHS
            Y = B(I1,J)
            Z = B(I,J)
            B(I1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
  590     CONTINUE
  600     CONTINUE
          !
          IF (.NOT.MATU) GO TO 620
          !
          DO 610 J=1,M
            Y = U(J,I1)
            Z = U(J,I)
            U(J,I1) = Y*C + Z*S
            U(J,I) = -Y*S + Z*C
  610     CONTINUE
  !
  620   CONTINUE
        !
        RV1(L) = 0.d0
        RV1(K) = F
        W(K) = X
        GO TO 480
  !      ********** CONVERGENCE **********
  630   IF (Z.GE.0.d0) GO TO 650
        !      ********** W(K) IS MADE NON-NEGATIVE **********
        W(K) = -Z
        IF (.NOT.MATV) GO TO 650
        !
        DO 640 J=1,N
          V(J,K) = -V(J,K)
  640   CONTINUE
  !
  650 CONTINUE
      !
      GO TO 670
  !      ********** SET ERROR -- NO CONVERGENCE TO A
  !                 SINGULAR VALUE AFTER 30 ITERATIONS **********
  660 IERR = K
  670 RETURN
      !      ********** LAST CARD OF GRSVD **********
      END
      !       SUBROUTINE DSWAP(N, SX, INCX, SY, INCY)                           SSW   10
      ! C
      ! C     INTERCHANGES TWO VECTORS.
      ! C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
      ! C     JACK DONGARRA, LINPACK, 3/11/78.
      ! C
      !       REAL*8 SX(1), SY(1), STEMP
      !       INTEGER I, INCX, INCY, IX, IY, M, MP1, N
      ! C
      !       IF (N.LE.0) RETURN
      !       IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      ! C
      ! C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
      ! C         TO 1
      ! C
      !       IX = 1
      !       IY = 1
      !       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      !       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      !       DO 10 I=1,N
      !         STEMP = SX(IX)
      !         SX(IX) = SY(IY)
      !         SY(IY) = STEMP
      !         IX = IX + INCX
      !         IY = IY + INCY
      !    10 CONTINUE
      !       RETURN
      ! C
      ! C       CODE FOR BOTH INCREMENTS EQUAL TO 1
      ! C
      ! C
      ! C       CLEAN-UP LOOP
      ! C
      !    20 M = MOD(N,3)
      !       IF (M.EQ.0) GO TO 40
      !       DO 30 I=1,M
      !         STEMP = SX(I)
      !         SX(I) = SY(I)
      !         SY(I) = STEMP
      !    30 CONTINUE
      !       IF (N.LT.3) RETURN
      !    40 MP1 = M + 1
      !       DO 50 I=MP1,N,3
      !         STEMP = SX(I)
      !         SX(I) = SY(I)
      !         SY(I) = STEMP
      !         STEMP = SX(I+1)
      !         SX(I+1) = SY(I+1)
      !         SY(I+1) = STEMP
      !         STEMP = SX(I+2)
      !         SX(I+2) = SY(I+2)
      !         SY(I+2) = STEMP
      !    50 CONTINUE
      !       RETURN
      !       END
      REAL*8 FUNCTION SRELPR(DUMMY)
      REAL*8 DUMMY
      !
      !      SRELPR COMPUTES THE RELATIVE PRECISION OF THE FLOATING POINT
      !      ARITHMETIC OF THE MACHINE.
      !
      !      IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
      !      THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
      !      ASSUME THE COMPUTER HAS
      !
      !         B = BASE OF ARITHMETIC
      !         T = NUMBER OF BASE  B  DIGITS
      !
      !      THEN
      !
      !         SRELPR = B**(1-T)
      !
      REAL*8 S
      !
      SRELPR = 1.d0
   10 SRELPR = SRELPR/2.d0
      S = 1.d0 + SRELPR
      IF (S.GT.1.d0) GO TO 10
      SRELPR = 2.d0*SRELPR
      RETURN
      END
