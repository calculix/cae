      !
      !     CalculiX - A 3-dimensional finite element program
      !              Copyright (C) 1998-2019 Guido Dhondt
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
      ! -----MATRIX-VECTOR MULTIPLY FOR REAL SPARSE NONSYMMETRIC MATRICES---------
      !
      SUBROUTINE OPNONSYMt(n,p,W,U,ad,au,jq,irow)
      !
      implicit none
      !
      ! -----------------------------------------------------------------------
      integer  IROW(*),JQ(*),n,l,i,j,llast
      real*8   U(*),W(*),Au(*),AD(*),p(*)
      ! -----------------------------------------------------------------------
      ! >    SPARSE MATRIX-VECTOR MULTIPLY FOR LANCZS  U = A^T*W
      ! >    SEE USPEC SUBROUTINE FOR DESCRIPTION OF THE ARRAYS THAT DEFINE
      ! >    THE MATRIX
      ! >    the vector p is not needed but is kept for compatibility reasons
      ! >    with the calling program
      ! -----------------------------------------------------------------------
      !
      !      COMPUTE THE DIAGONAL TERMS
      DO 10 I = 1,N
         U(I) = AD(I)*W(I)+U(I)
 10   CONTINUE
      !
      !      COMPUTE BY COLUMN
      LLAST = 0
      DO 30 J = 1,N
         !
         DO 20 L = JQ(J),JQ(J+1)-1
            I = IROW(L)
            !
            U(j) = U(j) + Au(L)*W(i)
 !
 20      CONTINUE
 !
 30   CONTINUE
      !
      RETURN
      END




