      !      ALGORITHM 584, COLLECTED ALGORITHMS FROM ACM.
      !      ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.8, NO. 2,
      !      JUN., 1982, P. 210. Author: D.P. Laurie
      !      PROGRAM KONYN(OUTPUT,TAPE6=OUTPUT)
      SUBROUTINE CUBTRI(F, T, EPS, MCALLS, ANS, ERR, NCALLS, W, NW,&
       IDATA, RDATA, IER)
      !
      !      changed on June 6th 2011 in order to allow parallelization:
      !      the common statements were commented with cparallel, the data
      !      statements were converted.
      !
      !        ADAPTIVE CUBATURE OVER A TRIANGLE
      !
      !        PARAMETERS
      !           F     - USER SUPPLIED EXTERNAL FUNCTION OF THE FORM
      !                   F(X,Y,IDATA,RDATA)
      !                   WHERE X AND Y ARE THE CARTESIAN COORDINATES OF A
      !                   POINT IN THE PLANE, AND IDATA AND RDATA ARE INTEGER
      !                   AND REAL*8 VECTORS IN WHICH DATA MAY BE PASSED.
      !           T     - ARRAY OF DIMENSION (2,3) WHERE T(1,J) AND T(2,J)
      !                   ARE THE X AND Y COORDINATES OF THE J-TH VERTEX OF
      !                   THE GIVEN TRIANGLE (INPUT)
      !           EPS   - REQUIRED TOLERANCE (INPUT).  IF THE COMPUTED
      !                   INTEGRAL IS BETWEEN-1 AND 1, AN ABSOLUTE ERROR
      !                   TEST IS USED, ELSE A RELATIVE ERROR TEST IS USED.
      !           MCALLS- MAXIMUM PERMITTED NUMBER OF CALLS TO F (INPUT)
      !           ANS   - ESTIMATE FOR THE VALUE OF THE INTEGRAL OF F OVER
      !                   THE GIVEN TRIANGLE (OUTPUT)
      !           ERR   - ESTIMATED ABSOLUTE ERROR IN ANS (OUTPUT)
      !           NCALLS- ACTUAL NUMBER OF CALLS TO F (OUTPUT).  THIS
      !                   PARAMETER MUST BE INITIALIZED TO 0 ON THE FIRST
      !                   CALL TO CUBTRI FOR A GIVEN INTEGRAL (INPUT)
      !           W     - WORK SPACE.  MAY NOT BE DESTROYED BETWEEN CALLS TO
      !                   CUBTRI IF RESTARTING IS INTENDED
      !           NW    - LENGTH OF WORK SPACE (INPUT).
      !                   IF NW .GE. 3*(19+3*MCALLS)/38, TERMINATION DUE TO
      !                   FULL WORK SPACE WILL NOT OCCUR.
      !           IER   - TERMINATION INDICATOR (OUTPUT)
      !                   IER=0   NORMAL TERMINATION, TOLERANCE SATISFIED
      !                   IER=1   MAXIMUM NUMBER OF CALLS REACHED
      !                   IER=2   WORK SPACE FULL
      !                   IER=3   FURTHER SUBDIVISION OF TRIANGLES IMPOSSIBLE
      !                   IER=4   NO FURTHER IMPROVEMENT IN ACCURACY IS
      !                         POSSIBLE DUE TO ROUNDING ERRORS IN FUNCTION
      !                         VALUES
      !                   IER=5   NO FURTHER IMPROVEMENT IN ACCURACY IS
      !                         POSSIBLE BECAUSE SUBDIVISION DOES NOT
      !                         CHANGE THE ESTIMATED INTEGRAL. MACHINE
      !                         ACCURACY HAS PROBABLY BEEN REACHED BUT
      !                         THE ERROR ESTIMATE IS NOT SHARP ENOUGH.
      !
      !        CUBTRI IS DESIGNED TO BE CALLED REPEATEDLY WITHOUT WASTING
      !        EARLIER WORK.  THE PARAMETER NCALLS IS USED TO INDICATE TO
      !        CUBTRI AT WHAT POINT TO RESTART, AND MUST BE RE-INITIALIZED
      !        TO 0 WHEN A NEW INTEGRAL IS TO BE COMPUTED.  AT LEAST ONE OF
      !        THE PARAMETERS EPS, MCALLS AND NW MUST BE CHANGED BETWEEN
      !        CALLS TO CUBTRI, ACCORDING TO THE RETURNED VALUE OF IER. NONE
      !        OF THE OTHER PARAMETERS MAY BE CHANGED IF RESTARTING IS DONE.
      !        IF IER=3 IS ENCOUNTERED, THERE PROBABLY IS A SINGULARITY
      !        SOMEWHERE IN THE REGION.  THE ERROR MESSAGE INDICATES THAT
      !        FURTHER SUBDIVISION IS IMPOSSIBLE BECAUSE THE VERTICES OF THE
      !        SMALLER TRIANGLES PRODUCED WILL BEGIN TO COALESCE TO THE
      !        PRECISION OF THE COMPUTER.  THIS SITUATION CAN USUALLY BE
      !        RELIEVED BY SPECIFYING THE REGION IN SUCH A WAY THAT THE
      !        SINGULARITY IS LOCATED AT THE THIRD VERTEX OF THE TRIANGLE.
      !        IF IER=4 IS ENCOUNTERED, THE VALUE OF THE INTEGRAL CANNOT BE
      !        IMPROVED ANY FURTHER. THE ONLY EXCEPTION TO THIS OCCURS WHEN A
      !        FUNCTION WITH HIGHLY IRREGULAR BEHAVIOUR IS INTEGRATED (E.G.
      !        FUNCTIONS WITH JUMP DISCONTINUITIES OR VERY HIGHLY OSCILLATORY
      !        FUNCTIONS). IN SUCH A CASE THE USER CAN DISABLE THE ROUNDING
      !        ERROR TEST BY REMOVING THE IF STATEMENT SHORTLY AFTER STATEMENT
      !        NUMBER 70.
      !
      implicit none
      EXTERNAL F,rnderr
      INTEGER IDATA(1), IER, MCALLS, NCALLS, NW,jkp,i,j,k,l,maxc,maxk,&
        mw,nfe
      REAL*8 ALFA, ANS, ANSKP, AREA, EPS, ERR, ERRMAX, H, Q1, Q2, R1,R2,&
       RDATA(1), D(2,4), S(4), T(2,3), VEC(2,3), W(6,NW), X(2),zero,&
       point5,one,rnderr
      !        ACTUAL DIMENSION OF W IS (6,NW/6)
      !
      REAL*8 TANS, TERR, DZERO
      ! parallel      COMMON /CUBSTA/ TANS, TERR
      !        THIS COMMON IS REQUIRED TO PRESERVE TANS AND TERR BETWEEN CALLS
      !        AND TO SAVE VARIABLES IN FUNCTION RNDERR
      nfe=19
      s=(/1.d0,1.d0,1.d0,-1.d0/)
      d=reshape((/0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,0.0d0,1.0d0,1.0d0/),&
        (/2,4/))
      zero=0.d0
      one=1.d0
      dzero=0.d0
      point5=.5d0
      ! parallel      DATA NFE /19/, S(1), S(2), S(3), S(4) /3*1E0,-1E0/, D(1,1),
      ! parallel     * D(2,1) /0.0,0.0/, D(1,2), D(2,2) /0.0,1.0/, D(1,3), D(2,3)
      ! parallel     * /1.0,0.0/, D(1,4), D(2,4) /1.0,1.0/
      !        NFE IS THE NUMBER OF FUNCTION EVALUATIONS PER CALL TO CUBRUL.
      ! parallel      DATA ZERO /0.E0/, ONE /1.E0/, DZERO /0.D0/, POINT5 /.5E0/
      !
      !       CALCULATE DIRECTION VECTORS, AREA AND MAXIMUM NUMBER
      !       OF SUBDIVISIONS THAT MAY BE PERFORMED
      DO 20 I=1,2
        VEC(I,3) = T(I,3)
        DO 10 J=1,2
          VEC(I,J) = T(I,J) - T(I,3)
   10   CONTINUE
   20 CONTINUE
      MAXC = (MCALLS/NFE+3)/4
      IER = 1
      MAXK = MIN0(MAXC,(NW/6+2)/3)
      IF (MAXC.GT.MAXK) IER = 2
      AREA = ABS(VEC(1,1)*VEC(2,2)-VEC(1,2)*VEC(2,1))*POINT5
      K = (NCALLS/NFE+3)/4
      MW = 3*(K-1) + 1
      IF (NCALLS.GT.0) GO TO 30
      !
      !        TEST FOR TRIVIAL CASES
      TANS = DZERO
      TERR = DZERO
      IF (AREA.EQ.ZERO) GO TO 90
      IF (MCALLS.LT.NFE) GO TO 100
      IF (NW.LT.6) GO TO 110
      !
      !        INITIALIZE DATA LIST
      K = 1
      MW = 1
      W(1,1) = ZERO
      W(2,1) = ZERO
      W(3,1) = ONE
      CALL CUBRUL(F, VEC, W(1,1), IDATA, RDATA)
      TANS = W(5,1)
      TERR = W(6,1)
      NCALLS = NFE
   !
   !        TEST TERMINATION CRITERIA
   30 ANS = TANS
      ERR = TERR
      IF (ERR.LT.DMAX1(ONE,ABS(ANS))*EPS) GO TO 90
      IF (K.EQ.MAXK) GO TO 120
      !
      !        FIND TRIANGLE WITH LARGEST ERROR
      ERRMAX = ZERO
      DO 40 I=1,MW
        IF (W(6,I).LE.ERRMAX) GO TO 40
        ERRMAX = W(6,I)
        J = I
   40 CONTINUE
      !
      !        SUBDIVIDE TRIANGLE INTO FOUR SUBTRIANGLES AND UPDATE DATA LIST
      DO 50 I=1,2
        X(I) = W(I,J)
   50 CONTINUE
      H = W(3,J)*POINT5
      IF (RNDERR(X(1),H,X(1),H).NE.ZERO) GO TO 130
      IF (RNDERR(X(2),H,X(2),H).NE.ZERO) GO TO 130
      ANSKP = (TANS)
      TANS = TANS - (W(5,J))
      TERR = TERR - (W(6,J))
      R1 = W(4,J)
      R2 = W(5,J)
      JKP = J
      Q1 = ZERO
      Q2 = ZERO
      DO 70 I=1,4
        DO 60 L=1,2
          W(L,J) = X(L) + H*D(L,I)
   60   CONTINUE
        W(3,J) = H*S(I)
        CALL CUBRUL(F, VEC, W(1,J), IDATA, RDATA)
        Q2 = Q2 + W(5,J)
        Q1 = Q1 + W(4,J)
        J = MW + I
   70 CONTINUE
      ALFA = 1E15
      IF (Q2.NE.R2) ALFA = ABS((Q1-R1)/(Q2-R2)-ONE)
      J = JKP
      DO 80 I=1,4
        W(6,J) = W(6,J)/ALFA
        TANS = TANS + W(5,J)
        TERR = TERR + W(6,J)
        J = MW + I
   80 CONTINUE
      MW = MW + 3
      NCALLS = NCALLS + 4*NFE
      K = K + 1
      !
      !        IF ANSWER IS UNCHANGED, IT CANNOT BE IMPROVED
      IF (ANSKP.EQ.(TANS)) GO TO 150
      !
      !        REMOVE THIS IF STATEMENT TO DISABLE ROUNDING ERROR TEST
      IF (K.GT.3 .AND. ABS(Q2-R2).GT.ABS(Q1-R1)) GO TO 140
      GO TO 30
   !
   !        EXITS FROM SUBROUTINE
   90 IER = 0
      GO TO 120
  100 IER = 1
      GO TO 120
  110 IER = 2
  120 ANS = TANS
      ERR = TERR
      RETURN
  130 IER = 3
      GO TO 120
  140 IER = 4
      GO TO 120
  150 IER = 5
      GO TO 120
      END
      real*8 FUNCTION RNDERR(X, A, Y, B)
      !        THIS FUNCTION COMPUTES THE ROUNDING ERROR COMMITTED WHEN THE
      !        SUM X+A IS FORMED.  IN THE CALLING PROGRAM, Y MUST BE THE SAME
      !        AS X AND B MUST BE THE SAME AS A.  THEY ARE DECLARED AS
      !        DISTINCT VARIABLES IN THIS FUNCTION, AND THE INTERMEDIATE
      !        VARIABLES S AND T ARE PUT INTO COMMON, IN ORDER TO DEFEND
      !        AGAINST THE WELL-MEANING ACTIONS OF SOME OFFICIOUS OPTIMIZING
      !        FORTRAN COMPILERS.
      implicit none
      real*8 x,a,y,b,s,t
      ! parallel      COMMON /CUBATB/ S, T
      S = X + A
      T = S - Y
      RNDERR = T - B
      RETURN
      END
      SUBROUTINE CUBRUL(F, VEC, P, IDATA, RDATA)
      !
      !        BASIC CUBATURE RULE PAIR OVER A TRIANGLE
      !
      !        PARAMETERS
      !          F  - EXTERNAL FUNCTION - SEE COMMENTS TO CUBTRI
      !          VEC- MATRIX OF BASE VECTORS AND ORIGIN (INPUT)
      !          P  - TRIANGLE DESCRIPTION VECTOR OF DIMENSION 6
      !                P(1) - TRANSFORMED X COORDINATE OF ORIGIN VERTEX(INPUT)
      !                P(2) - TRANSFORMED Y COORDINATE OF ORIGIN VERTEX(INPUT)
      !                P(3) - DISTANCE OF OTHER VERTICES IN THE DIRECTIONS
      !                      OF THE BASE VECTORS (INPUT)
      !                P(4) - LESS ACCURATE ESTIMATED INTEGRAL (OUTPUT)
      !                P(5) - MORE ACCURATE ESTIMATED INTEGRAL (OUTPUT)
      !                P(6) - ABS(P(5)-P(4))   (OUTPUT)
      !
      !        CUBRUL EVALUATES A LINEAR COMBINATION OF BASIC INTEGRATION
      !        RULES HAVING D3 SYMMETRY.  THE AREAL*8 COORDINATES PERTAINING TO
      !        THE J-TH RULE ARE STORED IN W(I,J),I=1,2,3.  THE CORRESPONDING
      !        WEIGHTS ARE W(4,J) AND W(5,J), WITH W(5,J) BELONGING TO THE
      !        MORE ACCURATE FORMULA.  IF W(1,J).EQ.W(2,J), THE INTEGRATION
      !        POINT IS THE CENTROID, ELSE IF W(2,J).EQ.W(3,J), THE EVALUATION
      !        POINTS ARE ON THE MEDIANS.  IN BOTH CASES ADVANTAGE IS TAKEN OF
      !        SYMMETRY TO AVOID REPEATING FUNCTION EVALUATIONS.
      !
      !        THE FOLLOWING REAL*8 VARIABLES ARE USED TO AVOID
      !        UNNECESSARY ROUNDING ERRORS IN FLOATING POINT ADDITION.
      !        THEY MAY BE DECLARED SINGLE PRECISION IF REAL*8 IS
      !        NOT AVAILABLE AND FULL ACCURACY IS NOT NEEDED.
      implicit none
      REAL*8 A1, A2, S, SN, DZERO, DONE, DTHREE, DSIX,f,&
        point5,x,y
      REAL*8 AREA, ORIGIN(2), P(6), RDATA(1), TVEC(2,3), VEC(2,3),W(5,6)
      INTEGER IDATA(1),nquad,i,j,k
      !
      !        W CONTAINS POINTS AND WEIGHTS OF THE INTEGRATION FORMULAE
      !        NQUAD - NUMBER OF BASIC RULES USED
      !
      !        THIS PARTICULAR RULE IS THE 19 POINT EXTENSION (DEGREE 8) OF
      !        THE FAMILIAR 7 POINT RULE (DEGREE 5).
      !
      !      SIGMA=SQRT(7)
      !      PHI=SQRT(15)
      !      W(1,1),W(2,1),W(3,1) = 1/3
      !      W(4,1) = 9/40
      !      W(5,1) = 7137/62720 - 45*SIGMA/1568
      !      W(1,2) = 3/7 + 2*PHI/21
      !      W(2,2),W(3,2) = 2/7 - PHI/21
      !      W(4,2) = 31/80 - PHI/400
      !      W(5,2) = - 9301697/4695040 - 13517313*PHI/23475200
      !             + 764885*SIGMA/939008 + 198763*PHI*SIGMA/939008
      !      W(*,3) = W(*,2) WITH PHI REPLACED BY -PHI
      !      W(1,5) = 4/9 + PHI/9 + SIGMA/9 - SIGMA*PHI/45
      !      W(2,5),W(3,5) = 5/18 - PHI/18 - SIGMA/18 + SIGMA*PHI/90
      !      W(4,5) = 0
      !      W(5,5) = 102791225/59157504 + 23876225*PHI/59157504
      !             - 34500875*SIGMA/59157504 - 9914825*PHI*SIGMA/59157504
      !      W(*,4) = W(*,5) WITH PHI REPLACED BY -PHI
      !      W(1,6) = 4/9 + SIGMA/9
      !      W(2,6) = W(2,4)
      !      W(3,6) = W(2,5)
      !      W(4,6) = 0
      !      W(5,6) = 11075/8064 - 125*SIGMA/288
      !
      nquad=6
      w=reshape((/.3333333333333333333333333d0,&
       .3333333333333333333333333d0,.3333333333333333333333333d0,.225d0,&
       .3786109120031468330830822d-1,.7974269853530873223980253d0,&
       .1012865073234563388009874d0,.1012865073234563388009874d0,&
       .3778175416344814577870518d0,&
       .1128612762395489164329420d0,.5971587178976982045911758d-1,&
       .4701420641051150897704412d0,.4701420641051150897704412d0,&
       .3971824583655185422129482d0,&
       .2350720567323520126663380d0,.5357953464498992646629509d0,&
       .2321023267750503676685246d0,.2321023267750503676685246d0,&
       0.d0,.3488144389708976891842461d0,&
       .9410382782311208665596304d0,&
       .2948086088443956672018481d-1,.2948086088443956672018481d-1,&
       0.d0,.4033280212549620569433320d-1,.7384168123405100656112906d0,&
       .2321023267750503676685246d0,.2948086088443956672018481d-1,&
       0.d0,.2250583347313904927138324d0/),(/5,6/))
      ! parallel
      ! parallel      DATA NQUAD /6/, W(1,1), W(2,1), W(3,1)/3*.33333333333333333333333
      ! parallel     * 33E0/, W(4,1), W(5,1) /.225E0,.3786109120031468330830822E-1/,
      ! parallel     * W(1,2), W(2,2), W(3,2) /.7974269853530873223980253E0,2*
      ! parallel     * .1012865073234563388009874E0/, W(4,2), W(5,2)
      ! parallel     * /.3778175416344814577870518E0,.1128612762395489164329420E0/,
      ! parallel     * W(1,3), W(2,3), W(3,3) /.5971587178976982045911758E-1,2*
      ! parallel     * .4701420641051150897704412E0/, W(4,3), W(5,3)
      ! parallel     * /.3971824583655185422129482E0,.2350720567323520126663380E0/
      ! parallel      DATA W(1,4), W(2,4), W(3,4) /.5357953464498992646629509E0,2*
      ! parallel     * .2321023267750503676685246E0/, W(4,4), W(5,4)
      ! parallel     * /0.E0,.3488144389708976891842461E0/, W(1,5), W(2,5), W(3,5)
      ! parallel     * /.9410382782311208665596304E0,2*.2948086088443956672018481E-1/,
      ! parallel     * W(4,5), W(5,5) /0.E0,.4033280212549620569433320E-1/, W(1,6),
      ! parallel     * W(2,6), W(3,6) /.7384168123405100656112906E0,
      ! parallel     * .2321023267750503676685246E0,.2948086088443956672018481E-1/,
      ! parallel     * W(4,6), W(5,6) /0.E0,.2250583347313904927138324E0/
      ! parallel
      !
      dzero=0.d0
      done=1.d0
      dthree=3.d0
      dsix=6.d0
      point5=.5d0
      ! parallel      DATA DZERO /0.D0/, DONE /1.D0/, DTHREE /3.D0/, DSIX /6.D0/,
      ! parallel     * POINT5 /.5E0/
      !
      !        SCALE BASE VECTORS AND OBTAIN AREA
      DO 20 I=1,2
        ORIGIN(I) = VEC(I,3) + P(1)*VEC(I,1) + P(2)*VEC(I,2)
        DO 10 J=1,2
          TVEC(I,J) = P(3)*VEC(I,J)
   10   CONTINUE
   20 CONTINUE
      AREA = POINT5*ABS(TVEC(1,1)*TVEC(2,2)-TVEC(1,2)*TVEC(2,1))
      A1 = DZERO
      A2 = DZERO
      !
      !        COMPUTE ESTIMATES FOR INTEGRAL AND ERROR
      DO 40 K=1,NQUAD
        X = ORIGIN(1) + W(1,K)*TVEC(1,1) + W(2,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(1,K)*TVEC(2,1) + W(2,K)*TVEC(2,2)
        S = (F(X,Y,IDATA,RDATA))
        SN = DONE
        IF (W(1,K).EQ.W(2,K)) GO TO 30
        X = ORIGIN(1) + W(2,K)*TVEC(1,1) + W(1,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(2,K)*TVEC(2,1) + W(1,K)*TVEC(2,2)
        S = S + (F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(2,K)*TVEC(1,1) + W(3,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(2,K)*TVEC(2,1) + W(3,K)*TVEC(2,2)
        S = S + (F(X,Y,IDATA,RDATA))
        SN = DTHREE
        IF (W(2,K).EQ.W(3,K)) GO TO 30
        X = ORIGIN(1) + W(1,K)*TVEC(1,1) + W(3,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(1,K)*TVEC(2,1) + W(3,K)*TVEC(2,2)
        S = S + (F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(3,K)*TVEC(1,1) + W(1,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(3,K)*TVEC(2,1) + W(1,K)*TVEC(2,2)
        S = S + (F(X,Y,IDATA,RDATA))
        X = ORIGIN(1) + W(3,K)*TVEC(1,1) + W(2,K)*TVEC(1,2)
        Y = ORIGIN(2) + W(3,K)*TVEC(2,1) + W(2,K)*TVEC(2,2)
        S = S + (F(X,Y,IDATA,RDATA))
        SN = DSIX
   30   S = S/SN
        A1 = A1 + W(4,K)*S
        A2 = A2 + W(5,K)*S
   40 CONTINUE
      P(4) = (A1)*AREA
      P(5) = (A2)*AREA
      P(6) = ABS(P(5)-P(4))
      RETURN
      END
