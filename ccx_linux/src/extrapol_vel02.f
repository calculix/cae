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
      subroutine extrapol_vel02(ipbounv2,ibounv2,xbounv2,xbounact,vold,
     &     mi,nbpv2a,nbpv2b)
!
!     application of boundary conditions to the nodes
!     velocity in y-direction
!     
      implicit none
!
      integer ipbounv2(2,*),node,nbpv2a,nbpv2b,mi(*),i,j,ibounv2(*)
!
      real*8 vold(0:mi(2),*),xbounact(*),xbounv2(*)
!
!
!
      do i=nbpv2a,nbpv2b
        node=ipbounv2(1,i)
        vold(2,node)=0.d0
        do j=ipbounv2(2,i)+1,ipbounv2(2,i+1)
          vold(2,node)=vold(2,node)+xbounact(ibounv2(j))
        enddo
        vold(2,node)=vold(2,node)*xbounv2(i)
      enddo
!
      return
      end
