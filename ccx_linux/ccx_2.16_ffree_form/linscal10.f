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
      subroutine linscal10(scal,konl,scall,idim,shp)
      !
      !     calculates a linear approximation to the quadratic interpolation
      !     of the temperatures in a C3D10 element. A
      !     quadratic interpolation of the temperatures leads to quadratic
      !     thermal stresses, which cannot be handled by the elements
      !     displacement functions (which lead to linear stresses). Thus,
      !     the temperatures are approximated by a linear function.
      !
      implicit none
      !
      integer konl(*),idim
      !
      real*8 scal(0:idim,*),scall,shp(4,*)
      !
      scall=&
          (shp(4,1)+(shp(4,5)+shp(4,7)+shp(4,8))/2.d0)*scal(0,konl(1))&
          +(shp(4,2)+(shp(4,5)+shp(4,6)+shp(4,9))/2.d0)*scal(0,konl(2))&
          +(shp(4,3)+(shp(4,6)+shp(4,7)+shp(4,10))/2.d0)*scal(0,konl(3))&
          +(shp(4,4)+(shp(4,8)+shp(4,9)+shp(4,10))/2.d0)*scal(0,konl(4))
      !
      return
      end
