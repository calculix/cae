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
      subroutine resultsp(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,
     &  mi)
!
!     calculates the pressure correction (STEP 2) in the nodes
!
      implicit none
!
      integer ipompc(*),nodempc(3,*),nmpc,nk,nactdoh(0:4,*),i,ist,
     &  node,ndir,index,mi(*)
!
      real*8 coefmpc(*),sol(*),v(0:mi(2),*),fixed_disp
!
!     extracting the pressure correction from the solution
!
      do i=1,nk
         if(nactdoh(4,i).gt.0) then
            v(4,i)=sol(nactdoh(4,i))
c            write(*,*) 'dpressureee ',i,v(4,i)
         else
            v(4,i)=0.d0
         endif
      enddo
c      write(*,*) 'sol307',v(4,307)
!     
!     inserting the mpc information: it is assumed that the
!     temperature MPC's also apply to the pressure
!     
c      do i=1,nmpc
c         ist=ipompc(i)
c         node=nodempc(1,ist)
c         ndir=nodempc(2,ist)
c         if(ndir.ne.0) cycle
c         index=nodempc(3,ist)
c         fixed_disp=0.d0
c         if(index.ne.0) then
c            do
c               fixed_disp=fixed_disp-coefmpc(index)*
c     &              v(4,nodempc(1,index))
c               index=nodempc(3,index)
c               if(index.eq.0) exit
c            enddo
c         endif
c         fixed_disp=fixed_disp/coefmpc(ist)
c         v(4,node)=fixed_disp
c      enddo
!
      return
      end
