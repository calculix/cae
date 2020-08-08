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
      subroutine correctvel(adv,nef,volume,gradpcel,vel,nefa,nefb)
!
!     correction of the velocity at the element centers due to the
!     pressure change (balance of mass)
!
!     the solution is stored in field bv.
!
      implicit none
!
      integer i,k,nef,nefa,nefb
!
      real*8 adv(*),volume(*),gradpcel(3,*),vel(nef,0:7)
!
!
!
      do i=nefa,nefb
         do k=1,3
            vel(i,k)=vel(i,k)-volume(i)*gradpcel(k,i)/adv(i)
         enddo
      enddo
!  
      return
      end
