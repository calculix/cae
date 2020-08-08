!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2020 Guido Dhondt
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
!     Calculating the residual in CFD-calculations: 
!     special case: dp-equation (change in pressure)
!
      subroutine calcrespfluid(n,b,au,x,res)
!
      implicit none
!
      integer i,n,nflnei
!
      real*8 xmax,xmin,b(*),au(*),x(*),res
!
      xmax=0.d0
      xmin=0.d0
      res=0.d0
!
      do i=1,n
         res=res+(b(i)/au(i))**2
         xmax=max(xmax,x(i))
         xmin=min(xmin,x(i))
      enddo
      res=dsqrt(res/n)/(xmax-xmin)
!
      return
      end
