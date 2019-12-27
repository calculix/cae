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
!     Calculating the residual in CFD-calculations
!
      subroutine calcresvfluid2(n,x,xmax,nef)
      !
      implicit none
      !
      integer i,n,nef
      !
      real*8 xmax,x(*),vel
      !
      xmax=0.d0
      !
      do i=1,n
         vel=dsqrt(x(i)**2+x(nef+i)**2+x(2*nef+i)**2)
         if(vel.gt.xmax) xmax=vel
      enddo
      !
      return
      end
