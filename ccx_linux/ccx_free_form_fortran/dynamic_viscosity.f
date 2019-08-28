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
!     This subroutine enables to calculate the dynamic viscosity of air
!     using Sutherland's formula for gas
!     dvi=dvi0*dsqrt(T/273.d0)*(1.d0+113/273.15)/(1+113/T)
!
!     author: Yannick Muller
!
      subroutine dynamic_viscosity (kgas,T,dvi)
      !
      implicit none
      !
      integer kgas
      !
      real*8 T,dvi
      kgas=kgas
      !
      dvi=0.00001711d0*dsqrt(T/273.15d0)*(1d0+113d0/273.15d0)&
           /(1.d0+113.d0/T)
      !
      return
      !
      end
!
