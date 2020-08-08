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
      subroutine correctvel_simplec(adv,nef,volume,gradpel,vel,
     &  ipnei,auv)
!
!     correction of the velocity at the element centers due to the
!     pressure change (balance of mass)
!
!     the solution is stored in field bv.
!
      implicit none
!
      integer i,k,nef,ipnei(*),indexf
!
      real*8 adv(*),volume(*),gradpel(3,*),vel(nef,0:7),auv(*),a1
!
      do i=1,nef
!
         a1=adv(i)
         do indexf=ipnei(i)+1,ipnei(i+1)
            a1=a1+auv(indexf)
         enddo
!
         do k=1,3
            vel(i,k)=vel(i,k)-volume(i)*gradpel(k,i)/a1
         enddo
      enddo
!  
      return
      end
