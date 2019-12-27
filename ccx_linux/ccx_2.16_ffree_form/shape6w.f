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
      subroutine shape6w(xi,et,ze,xl,xsj,shp,iflag)
      !
      !     shape functions and derivatives for a 6-node linear
      !     isoparametric wedge element. 0<=xi,et<=1,xi+et<=1,-1<=ze<=1.
      !
      !     iflag=1: calculate only the value of the shape functions
      !     iflag=2: calculate the value of the shape functions and
      !              the Jacobian determinant
      !     iflag=3: calculate the value of the shape functions, the
      !              value of their derivatives w.r.t. the global
      !              coordinates and the Jacobian determinant
      !
      !
      !    Copyright (c) 2003 WB
      !
      !    Written January 2003 on the basis of the Guido's shape function files
      !
      implicit none
      !
      integer i,j,k,iflag
      !
      real*8 shp(4,6),xs(3,3),xsi(3,3),xl(3,6),sh(3)
      !
      real*8 xi,et,ze,xsj,a
      !
      intent(in) xi,et,ze,xl,iflag
      !
      intent(out) shp,xsj
      !
      !     shape functions and their glocal derivatives
      !
      a=1.d0-xi-et
      !
      !     shape functions
      !
      shp(4, 1)=0.5d0*a *(1.d0-ze)
      shp(4, 2)=0.5d0*xi*(1.d0-ze)
      shp(4, 3)=0.5d0*et*(1.d0-ze)
      shp(4, 4)=0.5d0*a *(1.d0+ze)
      shp(4, 5)=0.5d0*xi*(1.d0+ze) 
      shp(4, 6)=0.5d0*et*(1.d0+ze) 
      !
      if(iflag.eq.1) return
      !
      !     local derivatives of the shape functions: xi-derivative
      !
      shp(1, 1)=-0.5d0*(1.d0-ze)
      shp(1, 2)= 0.5d0*(1.d0-ze)
      shp(1, 3)= 0.d0
      shp(1, 4)=-0.5d0*(1.d0+ze)
      shp(1, 5)= 0.5d0*(1.d0+ze)
      shp(1, 6)= 0.d0
      !
      !     local derivatives of the shape functions: eta-derivative
      !
      shp(2, 1)=-0.5d0*(1.d0-ze)
      shp(2, 2)= 0.d0
      shp(2, 3)= 0.5d0*(1.d0-ze)
      shp(2, 4)=-0.5d0*(1.d0+ze)
      shp(2, 5)= 0.d0
      shp(2, 6)= 0.5d0*(1.d0+ze)

      !
      !     local derivatives of the shape functions: zeta-derivative
      !
      shp(3, 1)=-0.5d0*a
      shp(3, 2)=-0.5d0*xi
      shp(3, 3)=-0.5d0*et
      shp(3, 4)= 0.5d0*a
      shp(3, 5)= 0.5d0*xi
      shp(3, 6)= 0.5d0*et
      !
      !
      !     computation of the local derivative of the global coordinates
      !     (xs)
      !
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,6
            xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
          enddo
        enddo
      enddo
      !
      !     computation of the jacobian determinant
      !
      xsj=xs(1,1)*(xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))&
         -xs(1,2)*(xs(2,1)*xs(3,3)-xs(2,3)*xs(3,1))&
         +xs(1,3)*(xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))
      !
      if(iflag.eq.2) return
      !
      !     computation of the global derivative of the local coordinates
      !     (xsi) (inversion of xs)
      !
      xsi(1,1)=(xs(2,2)*xs(3,3)-xs(3,2)*xs(2,3))/xsj
      xsi(1,2)=(xs(1,3)*xs(3,2)-xs(1,2)*xs(3,3))/xsj
      xsi(1,3)=(xs(1,2)*xs(2,3)-xs(2,2)*xs(1,3))/xsj
      xsi(2,1)=(xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3))/xsj
      xsi(2,2)=(xs(1,1)*xs(3,3)-xs(3,1)*xs(1,3))/xsj
      xsi(2,3)=(xs(1,3)*xs(2,1)-xs(1,1)*xs(2,3))/xsj
      xsi(3,1)=(xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2))/xsj
      xsi(3,2)=(xs(1,2)*xs(3,1)-xs(1,1)*xs(3,2))/xsj
      xsi(3,3)=(xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2))/xsj
      !
      !     computation of the global derivatives of the shape functions
      !
      do k=1,6
        do j=1,3
          sh(j)=shp(1,k)*xsi(1,j)+shp(2,k)*xsi(2,j)+shp(3,k)*xsi(3,j)
        enddo
        do j=1,3
          shp(j,k)=sh(j)
        enddo
      enddo
      !
      return
      end
