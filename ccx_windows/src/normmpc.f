!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine normmpc(nmpc,ipompc,nodempc,coefmpc,inomat)
!     
!     normalizing the coefficients of MPC's for fluid applications
!     (CFD)
!     
      implicit none
!     
      integer i,node,index,nmpc,inomat(*),nodempc(3,*),ipompc(*)
!     
      real*8 coefmpc(*),size
!     
!     normalizing the MPC-coefficients
!     
      do i=1,nmpc
        index=ipompc(i)
!     
!     check whether fluid node
!     
        node=nodempc(1,index)
        if(inomat(node).eq.0) cycle
!     
!     calculating sum of the square of the MPC coefficients
!     
        size=coefmpc(index)**2
        do
          index=nodempc(3,index)
          if(index.eq.0) exit
          size=size+coefmpc(index)**2
        enddo
!
        size=dsqrt(size)
!     
!     normalizing all terms of the MPC
!     
        index=ipompc(i)
        do
          coefmpc(index)=coefmpc(index)/size
c          write(*,*) 'normmpc',i,coefmpc(index)
          index=nodempc(3,index)
          if(index.eq.0) exit
        enddo
      enddo
!     
      return
      end
      
