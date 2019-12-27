      !
      !  Module 496 in TOMS
      !  Based upon the LZ-algorithm
      !  (see L.C. Kaufman, ACM TOMS 1 (1975) pp. 271-281)
      !
      SUBROUTINE DLZHES(N, A, NA, B, NB, X, NX, WANTX)
      !  THIS SUBROUTINE REDUCES THE COMPLEX MATRIX A TO UPPER
      !  HESSENBERG FORM AND REDUCES THE COMPLEX MATRIX B TO
      !  TRIANGULAR FORM
      !  INPUT PARAMETERS..
      !  N   THE ORDER OF THE A AND B MATRICES
      !  A   A COMPLEX MATRIX
      !  NA  THE ROW DIMENSION OF THE A MATRIX
      !  B   A COMPLEX MATRIX
      !  NB  THE ROW DIMENSION OF THE B MATRIX
      !  NX  THE ROW DIMENSION OF THE X MATRIX
      !  WANTX A LOGICAL VARIABLE WHICH IS SET TO .TRUE. IF
      !        THE EIGENVECTORS ARE WANTED. OTHERWISE IT SHOULD
      !      BE SET TO .FALSE.
      !  OUTPUT PARAMETERS..
      !  A  ON OUTPUT A IS AN UPPER HESSENBERG MATRIX, THE
      !     ORIGINAL MATRIX HAS BEEN DESTROYED
      !  B  AN UPPER TRIANGULAR MATRIX, THE ORIGINAL MATRIX
      !     HAS BEEN DESTROYED
      !  X  CONTAINS THE TRANSFORMATIONS NEEDED TO COMPUTE
      !     THE EIGENVECTORS OF THE ORIGINAL SYSTEM
      implicit none
      !
      LOGICAL WANTX
      !
      integer k,nm1,nm2,jm2,jp1,ip1,imj,j,ii,im1,i,nx,nb,na,n
      !
      REAL*8 C, D
      !
      COMPLEX*16 Y, W, Z, A(NA,N), B(NB,N), X(NX,N)
      !
      NM1 = N - 1
      !  REDUCE B TO TRIANGULAR FORM USING ELEMENTARY
      !  TRANSFORMATIONS
      DO 80 I=1,NM1
        D = 0.D0
        IP1 = I + 1
        DO 10 K=IP1,N
          Y = B(K,I)
          C = DABS(DREAL(Y)) + DABS(DIMAG(Y))
          IF (C.LE.D) GO TO 10
          D = C
          II = K
   10   CONTINUE
        IF (D.EQ.0.D0) GO TO 80
        Y = B(I,I)
        IF (D.LE.DABS(DREAL(Y))+DABS(DIMAG(Y))) GO TO 40
        !  MUST INTERCHANGE
        DO 20 J=1,N
          Y = A(I,J)
          A(I,J) = A(II,J)
          A(II,J) = Y
   20   CONTINUE
        DO 30 J=I,N
          Y = B(I,J)
          B(I,J) = B(II,J)
          B(II,J) = Y
   30   CONTINUE
   40   DO 70 J=IP1,N
          Y = B(J,I)/B(I,I)
          IF (DREAL(Y).EQ.0.D0 .AND. DIMAG(Y).EQ.0.D0) GO TO 70
          DO 50 K=1,N
            A(J,K) = A(J,K) - Y*A(I,K)
   50     CONTINUE
          DO 60 K=IP1,N
            B(J,K) = B(J,K) - Y*B(I,K)
   60     CONTINUE
   70   CONTINUE
        B(IP1,I) = (0.D0,0.D0)
   80 CONTINUE
      !  INITIALIZE X
      IF (.NOT.WANTX) GO TO 110
      DO 100 I=1,N
        DO 90 J=1,N
          X(I,J) = (0.D0,0.D0)
   90   CONTINUE
        X(I,I) = (1.D0,0.D0)
  100 CONTINUE
  !  REDUCE A TO UPPER HESSENBERG FORM
  110 NM2 = N - 2
      IF (NM2.LT.1) GO TO 270
      DO 260 J=1,NM2
        JM2 = NM1 - J
        JP1 = J + 1
        DO 250 II=1,JM2
          I = N + 1 - II
          IM1 = I - 1
          IMJ = I - J
          W = A(I,J)
          Z = A(IM1,J)
          IF (DABS(DREAL(W))+DABS(DIMAG(W)).LE.DABS(DREAL(Z))&
           +DABS(DIMAG(Z))) GO TO 140
          !  MUST INTERCHANGE ROWS
          DO 120 K=J,N
            Y = A(I,K)
            A(I,K) = A(IM1,K)
            A(IM1,K) = Y
  120     CONTINUE
          DO 130 K=IM1,N
            Y = B(I,K)
            B(I,K) = B(IM1,K)
            B(IM1,K) = Y
  130     CONTINUE
  140     Z = A(I,J)
          IF (DREAL(Z).EQ.0.D0 .AND. DIMAG(Z).EQ.0.D0) GO TO 170
          Y = Z/A(IM1,J)
          DO 150 K=JP1,N
            A(I,K) = A(I,K) - Y*A(IM1,K)
  150     CONTINUE
          DO 160 K=IM1,N
            B(I,K) = B(I,K) - Y*B(IM1,K)
  160     CONTINUE
  !  TRANSFORMATION FROM THE RIGHT
  170     W = B(I,IM1)
          Z = B(I,I)
          IF (DABS(DREAL(W))+DABS(DIMAG(W)).LE.DABS(DREAL(Z))&
           +DABS(DIMAG(Z))) GO TO 210
          !  MUST INTERCHANGE COLUMNS
          DO 180 K=1,I
            Y = B(K,I)
            B(K,I) = B(K,IM1)
            B(K,IM1) = Y
  180     CONTINUE
          DO 190 K=1,N
            Y = A(K,I)
            A(K,I) = A(K,IM1)
            A(K,IM1) = Y
  190     CONTINUE
          IF (.NOT.WANTX) GO TO 210
          DO 200 K=IMJ,N
            Y = X(K,I)
            X(K,I) = X(K,IM1)
            X(K,IM1) = Y
  200     CONTINUE
  210     Z = B(I,IM1)
          IF (DREAL(Z).EQ.0.D0 .AND. DIMAG(Z).EQ.0.D0) GO TO 250
          Y = Z/B(I,I)
          DO 220 K=1,IM1
            B(K,IM1) = B(K,IM1) - Y*B(K,I)
  220     CONTINUE
          B(I,IM1) = (0.D0,0.D0)
          DO 230 K=1,N
            A(K,IM1) = A(K,IM1) - Y*A(K,I)
  230     CONTINUE
          IF (.NOT.WANTX) GO TO 250
          DO 240 K=IMJ,N
            X(K,IM1) = X(K,IM1) - Y*X(K,I)
  240     CONTINUE
  250   CONTINUE
        A(JP1+1,J) = (0.D0,0.D0)
  260 CONTINUE
  270 RETURN
      END
      SUBROUTINE DLZIT(N, A, NA, B, NB, X, NX, WANTX, ITER, EIGA,&
       EIGB)
      !  THIS SUBROUTINE SOLVES THE GENERALIZED EIGENVALUE PROBLEM
      !               A X  = LAMBDA B X
      !  WHERE A IS A COMPLEX UPPER HESSENBERG MATRIX OF
      !  ORDER N AND B IS A COMPLEX UPPER TRIANGULAR MATRIX OF ORDER N
      !  INPUT PARAMETERS
      !  N      ORDER OF A AND B
      !  A      AN N X N UPPER HESSENBERG COMPLEX MATRIX
      !  NA     THE ROW DIMENSION OF THE A MATRIX
      !  B      AN N X N UPPER TRIANGULAR COMPLEX MATRIX
      !  NB     THE ROW DIMENSION OF THE B MATRIX
      !  X      CONTAINS TRANSFORMATIONS TO OBTAIN EIGENVECTORS OF
      !         ORIGINAL SYSTEM. IF EIGENVECTORS ARE REQUESTED AND QZHES
      !         IS NOT CALLED, X SHOULD BE SET TO THE IDENTITY MATRIX
      !  NX     THE ROW DIMENSION OF THE X MATRIX
      !  WANTX  LOGICAL VARIABLE WHICH SHOULD BE SET TO .TRUE.
      !         IF EIGENVECTORS ARE WANTED. OTHERWISE IT
      !         SHOULD BE SET TO .FALSE.
      !  OUTPUT PARAMETERS
      !  X      THE ITH COLUMN CONTAINS THE ITH EIGENVECTOR
      !         IF EIGENVECTORS ARE REQUESTED
      !  ITER   AN INTEGER ARRAY OF LENGTH N WHOSE ITH ENTRY
      !         CONTAINS THE NUMBER OF ITERATIONS NEEDED TO FIND
      !         THE ITH EIGENVALUE. FOR ANY I IF ITER(I) =-1 THEN
      !         AFTER 30 ITERATIONS THERE HAS NOT BEEN A SUFFICIENT
      !         DECREASE IN THE LAST SUBDIAGONAL ELEMENT OF A
      !         TO CONTINUE ITERATING.
      !  EIGA   A COMPLEX ARRAY OF LENGTH N CONTAINING THE DIAGONAL OF A
      !  EIGB   A COMPLEX ARRAY OF LENGTH N CONTAINING THE DIAGONAL OF B
      !  THE ITH EIGENVALUE CAN BE FOUND BY DIVIDING EIGA(I) BY
      !  EIGB(I). WATCH OUT FOR EIGB(I) BEING ZERO
      COMPLEX*16 A(NA,N), B(NB,N), EIGA(N), EIGB(N)
      COMPLEX*16 S, W, Y, Z, CDSQRT
      COMPLEX*16 X(NX,N)
      INTEGER ITER(N)
      COMPLEX*16 ANNM1, ALFM, BETM, D, SL, DEN, NUM, ANM1M1
      REAL*8 EPSA, EPSB, SS, R, ANORM, BNORM, ANI, BNI, C
      REAL*8 D0, D1, D2, E0, E1
      LOGICAL WANTX
      NN = N
      !  COMPUTE THE MACHINE PRECISION TIMES THE NORM OF A AND B
      ANORM = 0.D0
      BNORM = 0.D0
      DO 30 I=1,N
        ANI = 0.D0
        IF (I.EQ.1) GO TO 10
        Y = A(I,I-1)
        ANI = ANI + DABS(DREAL(Y)) + DABS(DIMAG(Y))
   10   BNI = 0.D0
        DO 20 J=I,N
          ANI = ANI + DABS(DREAL(A(I,J))) + DABS(DIMAG(A(I,J)))
          BNI = BNI + DABS(DREAL(B(I,J))) + DABS(DIMAG(B(I,J)))
   20   CONTINUE
        IF (ANI.GT.ANORM) ANORM = ANI
        IF (BNI.GT.BNORM) BNORM = BNI
   30 CONTINUE
      IF (ANORM.EQ.0.D0) ANORM = 1.D0
      IF (BNORM.EQ.0.D0) BNORM = 1.D0
      EPSB = BNORM
      EPSA = ANORM
   40 EPSA = EPSA/2.D0
      EPSB = EPSB/2.D0
      C = ANORM + EPSA
      IF (C.GT.ANORM) GO TO 40
      IF (N.LE.1) GO TO 320
   50 ITS = 0
      NM1 = NN - 1
   !  CHECK FOR NEGLIGIBLE SUBDIAGONAL ELEMENTS
   60 D2 = DABS(DREAL(A(NN,NN))) + DABS(DIMAG(A(NN,NN)))
      DO 70 LB=2,NN
        L = NN + 2 - LB
        SS = D2
        Y = A(L-1,L-1)
        D2 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
        SS = SS + D2
        Y = A(L,L-1)
        R = SS + DABS(DREAL(Y)) + DABS(DIMAG(Y))
        IF (R.EQ.SS) GO TO 80
   70 CONTINUE
      L = 1
   80 IF (L.EQ.NN) GO TO 320
      IF (ITS.LT.30) GO TO 90
      ITER(NN) = -1
      IF (DABS(DREAL(A(NN,NM1)))+DABS(DIMAG(A(NN,NM1))).GT.0.8D0*&
       DABS(DREAL(ANNM1))+DABS(DIMAG(ANNM1))) RETURN
   90 IF (ITS.EQ.10 .OR. ITS.EQ.20) GO TO 110
      !  COMPUTE SHIFT AS EIGENVALUE OF LOWER 2 BY 2
      ANNM1 = A(NN,NM1)
      ANM1M1 = A(NM1,NM1)
      S = A(NN,NN)*B(NM1,NM1) - ANNM1*B(NM1,NN)
      W = ANNM1*B(NN,NN)*(A(NM1,NN)*B(NM1,NM1)-B(NM1,NN)*ANM1M1)
      Y = (ANM1M1*B(NN,NN)-S)/2.
      Z = CDSQRT(Y*Y+W)
      IF (DREAL(Z).EQ.0.D0 .AND. DIMAG(Z).EQ.0.D0) GO TO 100
      D0 = DREAL(Y/Z)
      IF (D0.LT.0.D0) Z = -Z
  100 DEN = (Y+Z)*B(NM1,NM1)*B(NN,NN)
      IF (DREAL(DEN).EQ.0.D0 .AND. DIMAG(DEN).EQ.0.D0) DEN =&
       CMPLX(EPSA,0.D0)
      NUM = (Y+Z)*S - W
      GO TO 120
  !  AD-HOC SHIFT
  110 Y = A(NM1,NN-2)
      NUM = CMPLX(DABS(DREAL(ANNM1))+DABS(DIMAG(ANNM1)),DABS(DREAL(Y))&
       +DABS(DIMAG(Y)))
      DEN = (1.D0,0.D0)
  !  CHECK FOR 2 CONSECUTIVE SMALL SUBDIAGONAL ELEMENTS
  120 IF (NN.EQ.L+1) GO TO 140
      D2 = DABS(DREAL(A(NM1,NM1))) + DABS(DIMAG(A(NM1,NM1)))
      E1 = DABS(DREAL(ANNM1)) + DABS(DIMAG(ANNM1))
      D1 = DABS(DREAL(A(NN,NN))) + DABS(DIMAG(A(NN,NN)))
      NL = NN - (L+1)
      DO 130 MB=1,NL
        M = NN - MB
        E0 = E1
        Y = A(M,M-1)
        E1 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
        D0 = D1
        D1 = D2
        Y = A(M-1,M-1)
        D2 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
        Y = A(M,M)*DEN - B(M,M)*NUM
        D0 = (D0+D1+D2)*(DABS(DREAL(Y))+DABS(DIMAG(Y)))
        E0 = E0*E1*(DABS(DREAL(DEN))+DABS(DIMAG(DEN))) + D0
        IF (E0.EQ.D0) GO TO 150
  130 CONTINUE
  140 M = L
  150 CONTINUE
      ITS = ITS + 1
      W = A(M,M)*DEN - B(M,M)*NUM
      Z = A(M+1,M)*DEN
      D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
      D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
      !  FIND L AND M AND SET A=LAM AND B=LBM
      !      NP1 = N + 1
      LOR1 = L
      NNORN = NN
      IF (.NOT.WANTX) GO TO 160
      LOR1 = 1
      NNORN = N
  160 DO 310 I=M,NM1
        J = I + 1
        !  FIND ROW TRANSFORMATIONS TO RESTORE A TO
        !  UPPER HESSENBERG FORM. APPLY TRANSFORMATIONS
        !  TO A AND B
        IF (I.EQ.M) GO TO 170
        W = A(I,I-1)
        Z = A(J,I-1)
        D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
        D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
        IF (D1.EQ.0.D0) GO TO 60
  170   IF (D2.GT.D1) GO TO 190
        !  MUST INTERCHANGE ROWS
        DO 180 K=I,NNORN
          Y = A(I,K)
          A(I,K) = A(J,K)
          A(J,K) = Y
          Y = B(I,K)
          B(I,K) = B(J,K)
          B(J,K) = Y
  180   CONTINUE
        IF (I.GT.M) A(I,I-1) = A(J,I-1)
        IF (D2.EQ.0.D0) GO TO 220
        !  THE SCALING OF W AND Z IS DESIGNED TO AVOID A DIVISION BY ZERO
        !  WHEN THE DENOMINATOR IS SMALL
        Y = CMPLX(DREAL(W)/D1,DIMAG(W)/D1)/CMPLX(DREAL(Z)/D1,DIMAG(Z)/&
         D1)
        GO TO 200
  190   Y = CMPLX(DREAL(Z)/D2,DIMAG(Z)/D2)/CMPLX(DREAL(W)/D2,DIMAG(W)/&
         D2)
  200   DO 210 K=I,NNORN
          A(J,K) = A(J,K) - Y*A(I,K)
          B(J,K) = B(J,K) - Y*B(I,K)
  210   CONTINUE
  220   IF (I.GT.M) A(J,I-1) = (0.D0,0.D0)
        !  PERFORM TRANSFORMATIONS FROM RIGHT TO RESTORE B TO
        !    TRIANGLULAR FORM
        !  APPLY TRANSFORMATIONS TO A
        Z = B(J,I)
        W = B(J,J)
        D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
        D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
        IF (D1.EQ.0.D0) GO TO 60
        IF (D2.GT.D1) GO TO 270
        !  MUST INTERCHANGE COLUMNS
        DO 230 K=LOR1,J
          Y = A(K,J)
          A(K,J) = A(K,I)
          A(K,I) = Y
          Y = B(K,J)
          B(K,J) = B(K,I)
          B(K,I) = Y
  230   CONTINUE
        IF (I.EQ.NM1) GO TO 240
        Y = A(J+1,J)
        A(J+1,J) = A(J+1,I)
        A(J+1,I) = Y
  240   IF (.NOT.WANTX) GO TO 260
        DO 250 K=1,N
          Y = X(K,J)
          X(K,J) = X(K,I)
          X(K,I) = Y
  250   CONTINUE
  260   IF (D2.EQ.0.D0) GO TO 310
        Z = CMPLX(DREAL(W)/D1,DIMAG(W)/D1)/CMPLX(DREAL(Z)/D1,DIMAG(Z)/&
         D1)
        GO TO 280
  270   Z = CMPLX(DREAL(Z)/D2,DIMAG(Z)/D2)/CMPLX(DREAL(W)/D2,DIMAG(W)/&
         D2)
  280   DO 290 K=LOR1,J
          A(K,I) = A(K,I) - Z*A(K,J)
          B(K,I) = B(K,I) - Z*B(K,J)
  290   CONTINUE
        B(J,I) = (0.D0,0.D0)
        IF (I.LT.NM1) A(I+2,I) = A(I+2,I) - Z*A(I+2,J)
        IF (.NOT.WANTX) GO TO 310
        DO 300 K=1,N
          X(K,I) = X(K,I) - Z*X(K,J)
  300   CONTINUE
  310 CONTINUE
      GO TO 60
  320 CONTINUE
      EIGA(NN) = A(NN,NN)
      EIGB(NN) = B(NN,NN)
      IF (NN.EQ.1) GO TO 330
      ITER(NN) = ITS
      NN = NM1
      IF (NN.GT.1) GO TO 50
      ITER(1) = 0
      GO TO 320
  !  FIND EIGENVECTORS USING B FOR INTERMEDIATE STORAGE
  330 IF (.NOT.WANTX) RETURN
      M = N
  340 CONTINUE
      ALFM = A(M,M)
      BETM = B(M,M)
      B(M,M) = (1.D0,0.D0)
      L = M - 1
      IF (L.EQ.0) GO TO 370
  350 CONTINUE
      L1 = L + 1
      SL = (0.D0,0.D0)
      DO 360 J=L1,M
        SL = SL + (BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
  360 CONTINUE
      Y = BETM*A(L,L) - ALFM*B(L,L)
      IF (DREAL(Y).EQ.0.D0 .AND. DIMAG(Y).EQ.0.D0) Y =&
       CMPLX((EPSA+EPSB)/2.D0,0.D0)
      B(L,M) = -SL/Y
      L = L - 1
  370 IF (L.GT.0) GO TO 350
      M = M - 1
      IF (M.GT.0) GO TO 340
      !  TRANSFORM TO ORIGINAL COORDINATE SYSTEM
      M = N
  380 CONTINUE
      DO 400 I=1,N
        S = (0.D0,0.D0)
        DO 390 J=1,M
          S = S + X(I,J)*B(J,M)
  390   CONTINUE
        X(I,M) = S
  400 CONTINUE
      M = M - 1
      IF (M.GT.0) GO TO 380
      !  NORMALIZE SO THAT LARGEST COMPONENT = 1.
      M = N
  410 CONTINUE
      SS = 0.D0
      DO 420 I=1,N
        R = DABS(DREAL(X(I,M))) + DABS(DIMAG(X(I,M)))
        IF (R.LT.SS) GO TO 420
        SS = R
        D = X(I,M)
  420 CONTINUE
      IF (SS.EQ.0.D0) GO TO 440
      DO 430 I=1,N
        X(I,M) = X(I,M)/D
  430 CONTINUE
  440 M = M - 1
      IF (M.GT.0) GO TO 410
      RETURN
      END
