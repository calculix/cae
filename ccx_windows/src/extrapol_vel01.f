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
      subroutine extrapol_vel01(ipbounv1,ibounv1,xbounv1,xbounact,vold,
     &     mi,nbpv1a,nbpv1b)
!
!     application of boundary conditions to the nodes
!     velocity in x-direction
!     
      implicit none
!
      integer ipbounv1(2,*),node,nbpv1a,nbpv1b,mi(*),i,j,ibounv1(*)
!
      real*8 vold(0:mi(2),*),xbounact(*),xbounv1(*)
!
!
!
      do i=nbpv1a,nbpv1b
        node=ipbounv1(1,i)
        vold(1,node)=0.d0
        do j=ipbounv1(2,i)+1,ipbounv1(2,i+1)
          vold(1,node)=vold(1,node)+xbounact(ibounv1(j))
        enddo
        vold(1,node)=vold(1,node)*xbounv1(i)
      enddo
!
      return
      end
