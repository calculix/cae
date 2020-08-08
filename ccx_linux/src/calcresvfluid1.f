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
!     Calculating the residual in CFD-calculations
!
      subroutine calcresvfluid1(n,a,b,au,ia,ja,x,res)
!
      implicit none
!
      integer i,j,n,ia(*),ja(*)
!
      real*8 a(*),b(*),au(*),x(*),res,resi,vel
!
      res=0.d0
!
      do i=1,n
         resi=a(ja(i)+1)*x(ia(ja(i)+1))
         do j=ja(i)+2,ja(i+1)
            resi=resi+a(j)*x(ia(j))
         enddo
         res=res+((resi-b(i))/au(i))**2
      enddo
!
      return
      end
