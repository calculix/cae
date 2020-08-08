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
      subroutine smoothshock(adb,aub,adl,addiv,sol,aux,icol,irow,jq,
     &  neq,nzl,sa)
!
!     smoothing the finite element solution
!
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer icol(*),irow(*),jq(*),neq,nzl,i,j,k
!
      real*8 adb(*),aub(*),adl(*),sol(*),aux(*),p,sa(*),addiv(*)
!
!     multiplying M-ML with the solution
!
c      call op(neq,p,sol,aux,adb,aub,icol,irow,nzl)
      call op(neq,sol,aux,adb,aub,jq,irow)
!
!     smoothing the solution
!
      do i=1,neq
         sol(i)=sol(i)+sa(i)*aux(i)*adl(i)
      enddo
!
      return
      end
