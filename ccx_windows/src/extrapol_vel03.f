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
      subroutine extrapol_vel03(ipbounv3,ibounv3,xbounv3,xbounact,vold,
     &     mi,nbpv3a,nbpv3b)
!
!     application of boundary conditions to the nodes
!     velocity in z-direction
!     
      implicit none
!
      integer ipbounv3(2,*),node,nbpv3a,nbpv3b,mi(*),i,j,ibounv3(*)
!
      real*8 vold(0:mi(2),*),xbounact(*),xbounv3(*)
!
!
!
      do i=nbpv3a,nbpv3b
        node=ipbounv3(1,i)
        vold(3,node)=0.d0
        do j=ipbounv3(2,i)+1,ipbounv3(2,i+1)
          vold(3,node)=vold(3,node)+xbounact(ibounv3(j))
        enddo
        vold(3,node)=vold(3,node)*xbounv3(i)
      enddo
!
      return
      end
