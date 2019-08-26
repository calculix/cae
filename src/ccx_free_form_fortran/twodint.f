 !
 !     CalculiX - A 3-dimensional finite element program
 !     Copyright (C) 1998-2018 Guido Dhondt
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
 !     published by the Free Software Foundation(version 2);
 !
 !
 !     This program is distributed in the hope that it will be useful,
 !     but WITHOUT ANY WARRANTY; without even the implied warranty of
 !     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 !     GNU General Public License for more details.
 !
 !     You should have received a copy of the GNU General Public License
 !     along with this program; if not, write to the Free Software
 !     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 !
 !
 !   1.TASK          INTERPOLATION OF A TWO DIMENSIONAL FUNCTION DEFINED POINT BY POINT
 !   *********       THE X COORDINATES ARE USER SPECIFIED.
 !                   THE INTERPOLATION TYPE CAN BE INDEPENDANTLY CHOSEN IN THE TWO DIRECTIONS
 !                   EITHER CONSTANT, LINEAR OR DOUBLE QUADRATIC.
 !                   BEYOND THE FIELD OF INTERPOLATION AN EXTRAOLATION IS CARRIED OUT.
 !                   FOR ALL FOUR EXTRAPOLATION DIRECTIONS DIFFERENT EXTRAPOLATION METHOD
 !                  (C ONSTANT,LINEAR,QUADRATIC) CAN BE CHOSEN, WHICH ORDER MUST NOT BE HIGHER
 !                   THAN THE IONTERPOLATION ORDER
 !
 !   2.UP-AUFRUF     CALL TWODINT(T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
 !   ***********                 T = MATRIX OF THE SAMPLE POINTS FORMATED AS FOLLOW
 !                                   T(1,1) = NX + NY * 0.001
 !                                          NX = NUMBER OF LINES T
 !                                          NY = NUMBER OF COLUMNS T
 !                                   T(1,2) ... T(1,NY)
 !                                         VECTOR OF THE Y COORDINATES OF THE T MATRIX
 !                                   T(2,1) ... T(NX,1)
 !                                          VECTOR OF THE X COORDINATES OF THE T MATRIX
 !                                   REST  OF T-MATRIX:
 !                                           POINT(X,Y) OF THE T MATRIX
 !
 !                               LSP  = COLUMN STEPOF  T
 !                               IART = TYPE OF INTERPOLATION
 !                                      IART = INTX * 10 + INTY
 !                                      INTX INTERPOLATION TYPE IN X-DIRECTION
 !                                      INTY INTERPOLATION TYPE IN Y-DIRECTION
 !                               XA = VECTOR OF THE X COORDINATES OF THE VALUE TO BE INTERPOLATED
 !                               YA = VECTOR OF THE Y COORDINATES OF THE VALUE TO BE INTERPOLATED
 !                               ZA = VECTOR OF THE INTERPOLATED VALUES
 !                               NA = ACTUAL LENGTH OF THE 3 PREVIOUS VECTORS
 !                               IEXP = TWO ELEMENT VECTOR CONTRAINING THE TYPE OF EXTRAPOLATION
 !                                       CHOSEN BEYOND THE INTERPOLATION DOMAIN
 !                                      IEXP(1): EXTRAPOLATION IN X-DIRECTION
 !                                      IEXP(1) = IEXPX1 * 10 + IEXPXN
 !                                      IEXPX1: EXTRAPOLATION BENEATH THE FIRST POINT
 !                                      IEXPXN: EXTRAPOLATION BEYOND THE LAST POINT
 !                                      IEXP(2): EXTRAPOLATION IN Y-DIRECTION
 !                                      IEXP(2) = IEXPY1 * 10 + IEXPYN
 !                                      SAME METHOD AS FOR IEXP(1):
 !                               IER = ERROR CODE
 !                                      IER = 0: NORMAL PROCEEDING
 !                                      IER = -1: ERROR  INPUTDATA
 !
 !                            REMARK: CHOICE OF THE INTER- EXTRAPOLATION TYPE  IART AND IEXP -
 !                          --------             ASSIGNEMENT OF  INTX,INTY,IEXPX1,
 !                                        IEXPXN,IEXPY1,IEXPYN:
 !                                         = 0 :   CONSTANT
 !                                         = 1 :   LINEAR
 !                                         = 2 :   DOUBLE QUADRATIC FROM
 !                                                 THE  SECOND UNTIL PENULTIMATE
 !                                                 INTERVAL IN THE INTERPOLATION MATRIX T,OTHERWISE QUADRATIC
 !
 !  3.RESTRICTIONS   THE SAMPLING POINT VECTORS (X UND Y COORDINATES
 !  ***************  OF THE MATRICX T MUST BE  STRICTLY MONOTONIC INCREASING SORTED
 !                   THE PARAMETER FOR THE TYPE OF EXTRAPOLATION
 !                   MUST NOT BE GREATER THAN THE ONE FOR TH EINTERPOLATION TYPE
 !                   OTHERWISE THE  VALUE IS AUTOMATICALLY ADAPTATED
 !                   IF THE NUMBER OF THE SAMPLING POINTS FOR THE REQUIRED TYPE OF INTERPOLATION IS TOO SMALL,
 !                   THE DEGREE OF INTERPOLATION WILL BE ACCORDINGLY ADAPTATED
 !
 !   4.USED UP'S     ONEDINT (ONE DIMENSIONAL INTERPOLATION ANALOG TO THIS PROGRAMM)
 !

      SUBROUTINE TWODINT (T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
      implicit none
      INTEGER IEXP(2),IYU,IYO,IXU,IXO,IDX,IDY,LL,INPY,IEXPX1,IEXPXN,&
        IEXPY1,IEXPYN,LX,LY,INPX,IART,LSP,IER,NX,NY,L,NA,one
      REAL*8 T(LSP,1),XA(1),YA(1),ZA(1)
      REAL*8 Z1(4),Z2(4)
      !       ENTRY ZWEINT (T,LSP,IART,XA,YA,ZA,NA,IEXP,IER)
      IER = 0
      one=1
      NX = T(1,1)
      NY = (T(1,1)-NX)*1000 + 0.1d0
      !
      !  TESTING INPUT
      ! --------------
      IF ((NX-2).lt.0) then
         go to 900
      elseif((nx-2).eq.0) then
         go to 30
      else
         go to 10
      endif
   10 DO 20 L = 3,NX
   20 IF ((T(L,1)-T(L-1,1)) .LE. 0) GO TO 900
   30 IF ((NY-2).lt.0) then
         go to 900
      elseif((ny-2).eq.0) then
         go to 60
      else
         go to 40
      endif
   40 DO 50 L = 3,NY
   50 IF ((T(1,L)-T(1,L-1)) .LE. 0) GO TO 900
   60 IF (NA .LE. 0) GO TO 900
      !
      !  DEFINING THE CONTROL VALUES
      ! ---------------------------
      INPX = IART/10
      INPY = IART - INPX*10 + 0.1d0
      IEXPX1 = IEXP(1)/10
      IEXPXN = IEXP(1) - IEXPX1*10
      IEXPY1 = IEXP(2)/10
      IEXPYN = IEXP(2) - IEXPY1*10
      IF (NX-2 .LT. INPX) INPX = NX - 2
      IF (NY-2 .LT. INPY) INPY = NY - 2
      IF (IEXPX1 .GT. INPX) IEXPX1 = INPX
      IF (IEXPXN .GT. INPX) IEXPXN = INPX
      IF (IEXPY1 .GT. INPY) IEXPY1 = INPY
      IF (IEXPYN .GT. INPY) IEXPYN = INPY
      !
      !  SUCCESSIVE PROCESSING THE INTERPOLATION EXIGENCES
      ! -------------------------------------------------------
      DO 400 L = 1,NA
      LX = 2
  !
  !  SETTING REFERENCE POINTS (LX,LY)
  ! ---------------------------------
  200 IF (XA(L) .LT. T(LX,1)) GO TO 220
      LX = LX + 1
      IF ((LX-NX).le.0) then
         go to 200
      else
         go to 210
      endif
  210 LX = NX
  220 DO 230 LY = 2,NY
  230 IF (YA(L) .LT. T(1,LY)) GO TO 235
      LY = NY
  235 IYU = LY - INPY
      IYO = LY + INPY - 1
      IF (IYU .GE. 2) GO TO 240
      IYU = 2
      IYO = IYU + INPY
  240 IF (IYO .GT. NY) IYO = NY
      IXU = LX - INPX
      IXO = LX + INPX - 1
      IF (IXU .GE. 2) GO TO 245
      IXU = 2
      IXO = IXU + INPX
  245 IF (IXO .GT. NX) IXO = NX
      IDX = IXO - IXU + 1
      IF (IXU .LT. IXO) GO TO 270
      IF (IYU .LT. IYO) GO TO 250
      !
      !  CONSTANT INTERPOLATION
      ! ------------------------
      IF (LX .GT. 2 .AND. XA(L) .LT. T(NX,1)) LX = LX - 1
      IF (LY .GT. 2 .AND. YA(L) .LT. T(1,NY)) LY = LY - 1
      ZA(L) = T(LX,LY)
      GO TO 400
  !
  !  LINEAR AND  QUADRATIC INTERPOLATION USING ONEDINT (ONEDIMENSIONAL)
  ! ---------------------------------------------------------------------
  !
  !  INTERPOLATION ONLY IN Y-DIRECTION
  !
  250 IDY = 0
      DO 260 LL = IYU,IYO
      IDY = IDY + 1
      Z1(IDY) = T(1,LL)
  260 Z2(IDY) = T(LX,LL)
      GO TO 300
  !
  !  INTERPOLATION ONLY IN X-DIRECTION
  !
  270 IF (IYU .LT. IYO) GO TO 280
      CALL ONEDINT(T(IXU,1),T(IXU,LY),IDX,XA(L),ZA(L),one,INPX,IEXP(1),&
                 IER)
      IF (IER.eq.0) then
         go to 400
      else
         go to 900
      endif
  !
  !  1.INTERPOLATION STEP IN X-DIRECTION
  !
  280 IDY = 0
      DO 290 LL = IYU,IYO
      IDY = IDY + 1
      Z1(IDY) = T(1,LL)
      CALL ONEDINT (T(IXU,1),T(IXU,LL),IDX,XA(L),Z2(IDY),one,INPX,&
          IEXP(1),IER)
      IF (IER.eq.0) then
         go to 290
      else
         go to 900
      endif
  290 CONTINUE
  !
  !  1.OR 2.INTERPOLATION STEP IN Y-DIRECTION
  !
  300 CALL ONEDINT (Z1,Z2,IDY,YA(L),ZA(L),one,INPY,IEXP(2),IER)
      IF (IER.eq.0) then
         go to 400
      else
         go to 900
      endif
  !
  !  RETURN BY NORMAL PROCEEDING
  ! --------------------------------
  400 CONTINUE
      IER = 0
      RETURN
  !
  !  ERROR RETURN
  ! -------------
  900 IER = -1
      RETURN
      END
