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
      !   1. TASK         INTERPOLATION OF A FUNCTION DEFINED POINT BY POINT
      !   *********       THE X COORDINATES ARE USER SPECIFIED.
      !                   THE INTERPOLATION PROCESS CAN BE EITHER CONSTANT,LINEAR
      !                   OR EVEN dOUBLE QUADRATIC WITH EXTRAPOLATION USING THE
      !                   POLYNOM HIGHEST ORDER
      !                   thE DOUBLE QUADRATIC INTERPOLATION IS A 3RD ORDER METHOD
      !                   BY WHICH 2 PARABOLS ENCOMPASSING EACH 3 AND 4 sAMPLING POINTS
      !                   ARE DEFINED.
      !                   THE SOLUTION IS A LINEAR COMBINATION OF THE CONCERNED
      !                   PARABOLS VALUES DEPENDING ON THE DEFINITION OF THE ACTUAL
      !                   SAMPLING POINT INTERVAL
      !
      !
      !   2.INPUT     CALL ONEDINT(XE,YE,NE,XA,YA,NA,IART,IEXP,IER)
      !   ***********                XE = ABSCISSE VECTOR OF THE SAMPLING POINTS
      !                              YE = ORDINATE VECTOR OF THE SAMPLING POINTS
      !                              NE = LENGHT OF THE SAMPLING POINT VECTOR
      !                              XA = ASCISSE VECTOR OF THE INTERPOLATION POINT(INPUT)
      !                              YA = ORDINATE VECTOR OF THE INTERPOLATION POINT(OUTPUT)
      !                              NA = LENGTH OF THE INTERPOLATION VECTOR                         c                             IART = tYPE OF INTERPOLATION
      !                                     =0: CONSTANT
      !                                     =1: LINEAR
      !                                     =2: DOUBLE QUADRATIC
      !                              IEXP = TYPE OF EXTRAPOLATION
      !                                     IEXP = 10*IEX1 + IEXN
      !                                     IEX1 EXTRAPOLATIONS BEYOND THE
      !                                          1. SAMPLING POINT IN THE VECTOR
      !                                     IEXN EXTRAPOLATION BEYOND THE
      !                                          LAST SAMPLING POINT IN THE VECTOR
      !                                     SELECTION OF THE EXTRAPOLATION TYPE AS
      !                                     FOR IART.
      !                              IER  = ERROR CODE
      !                                     = 0: NORMAL PROCEEDING
      !                                     =-1:PROBLEM IN TH EGIVEN VALUES
      !                                          PROGRAMM STOPS.
      !
      !   3.RESTRICTION    ABSCISSE VECTOR XE MUST BE STRICTLY MONOTONIC INCREASING SORTED
      !   ***************  AUTOMATIC CONTROL INSIDE TEH SUBROUTINE:
      !                    NE = 0: ERROR INTERRUPTION
      !                    NE = 1: ONLY CONSTANT INTER- EXTRAPOLATION
      !                    NE = 2: MAXIMAL LINEAR INTER- EXTRAPOLATION
      !                    NE = 3: MAXIMAL QUADRATIC INTER- EXTRAPOLATIO
      !                    THE PARAMETER FOR THE TYPE OF EXTRAPOLATION
      !                    MUST NOT BE GREATER THAN THE ONE FOR TH EINTERPOLATION TYPE
      !                    OTHERWISE THE  VALUE IS AUTOMATICALLY ADAPTATED
      !
      SUBROUTINE ONEDINT(XE,YE,NE,XA,YA,NA,IART,IEXP,IER)
      implicit none
      INTEGER NE,NA,NA1,NE1,IG,IER,IA,IART,IE2,I,IEXP,IE1,L
      REAL*8 XE(NE),YE(NE),XA(NA),YA(NA),ZW1,ZW2,XO,YO,RAB,XD,YD,&
        XZ,YZ,XU,YU,EQ,EQD,X
      !
      !  INTERPOLATION FUNCTION
      !  ------------------------
      EQ(X) = YU + YU * (X-XU) / XU +&
              ((YZ-YU)/(XZ-XU) - YU/XU) * (X-XU) * X / XZ
      EQD(X) = YZ * X / XZ +&
              (YD / XD - YZ / XZ) * X * (X - XZ) / (XD - XZ)
      !
      !  INPUT/DATA TEST,INTERPOLATION DIVERGENCE,EXTRAPOLATION LIMIT
      ! ----------------------------------------------------------------
      NA1 = NA - 1
      IF (NA .LE. 0) GO TO 900
      NE1 = NE - 1
      IF (NE1.lt.0) then
         go to 900
      elseif(ne1.eq.0) then
         go to 22
      else
         go to 18
      endif
   18 DO 20 L = 1,NE1
   20 IF ((XE(L+1)-XE(L)) .LE. 0) GO TO 900
   22 IE1 = IEXP / 10
      IE2 = IEXP - 10*IE1
      IA = IART
      IF (NE1 .LT. IA)   IA = NE1
      IF (IA .LT. IE1)   IE1 = IA
      IF (IA .LT. IE2)   IE2 = IA
      !
      !  SUCCESSIVE PROCESSING THE INTERPOLATION EXIGENCES
      ! -------------------------------------------------------
      !
      !      ZUR ERHOEHUNG DER NUMERISCHEN GENAUIGKEIT WIRD EINE
      !      TRANSLATION VON (XO,YO) IN (0,0) DURCHGEFUEHRT. DIES
      !      BEWIRKT AUSSERDEM  EINE BESCHLEUNIGUNG DES VERFAHRENS.
      !
      DO 100 I = 1,NA
      DO 24 L = 1,NE
      IF (XA(I) .LT. XE(L)) GO TO 30
   24 CONTINUE
      L = NE
      IF ((IE2 - 1).lt.0) then
         go to 50
      elseif((ie2-1).eq.0) then
         go to 35
      else
         go to 70
      endif
   30 IF (L .GT. 1) GO TO 40
      IF ((IE1 - 1).lt.0) then
         go to 50
      elseif((ie1-1).eq.0) then
         go to 25
      else
         go to 70
      endif
   40 IF ((IA-1).lt.0) then
         go to 45
      elseif((ia-1).eq.0) then
         go to 60
      else
         go to 70
      endif
   !
   !  CONSTANT INTERPOLATION
   !  -----------------------
   45 L = L - 1
   50 YA(I) = YE(L)
      GO TO 100
   !
   !  LINEAR EXTRAPOLATION
   !  ------------------------------
   25 IF (IA .EQ. 1) GO TO 60
      XO = XE(2)
      XU = XE(1) - XO
      YO = YE(2)
      YU = YE(1) - YO
      XZ = XE(3) - XO
      YZ = YE(3) - YO
      GO TO 38
   35 IF (IA .EQ. 1) GO TO 60
      XO = XE(NE1)
      XZ = XE(NE1-1) - XO
      XU = XE(NE) - XO
      YO = YE(NE1)
      YZ = YE(NE1-1) - YO
      YU = YE(NE) - YO
   !
   !  LINEAR EXTRAPOLATION WITH QUADRATIC INTERPOLATION
   !  -----------------------------------------------------
   38 RAB = YU / XU + XU * ((YZ-YU) / (XZ-XU) - YU/XU) / XZ
      YA(I) = YU + YO + (XA(I) -XU-XO)*RAB
      GO TO 100
  !
  !  LINEAR INTERPOLATION
  !  ---------------------
  60  IG = L - 1
      IF (IG .LT. 1) IG = 1
      YA(I) = YE(IG) + (XA(I)-XE(IG))*(YE(IG+1)-YE(IG))&
              / (XE(IG+1)-XE(IG))
      GO TO 100
   70 IF (L .GT. 2) GO TO 80
      XO = XE(2)
      XU = XE(1) - XO
      YO = YE(2)
      YU = YE(1) - YO
      XZ = XE(3) - XO
      YZ = YE(3) - YO
      GO TO 85
   80 IF (L .LT. NE) GO TO 90
      XO = XE(NE1)
      XU = XE(NE1-1) - XO
      XZ = XE(NE) - XO
      YO = YE(NE1)
      YU = YE(NE1-1) - YO
      YZ = YE(NE) - YO
   85 YA(I) = EQ(XA(I)-XO) + YO
      GO TO 100
   !
   !  DOUBLE QUADRATIC INTERPOLATION
   !  ----------------------------------
   90 XO = XE(L-1)
      XU = XE(L-2) - XO
      XZ = XE(L) - XO
      XD = XE(L+1) - XO
      YO = YE(L-1)
      YU = YE(L-2) - YO
      YZ = YE(L) - YO
      YD = YE(L+1) - YO
      ZW1 = EQ(XA(I)-XO)
      ZW2 = EQD(XA(I)-XO)
      YA(I) = ZW1 + (ZW2 - ZW1) * (XA(I) - XO)/XZ + YO
  100 CONTINUE
      !
      !  RETURN BY NORMAL PROCEEDING
      !  -------------------------------
      IER = 0
      RETURN
  !
  !  ERROR RETURN
  !  ------------
  900 IER = -1
      RETURN
      END
