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
      subroutine extrapol_tel01(ipbount,ibount,xbount,xbounact,vold,
     &     mi,nbpta,nbptb)
!
!     application of boundary conditions to the nodes
!     pressure
!     
      implicit none
!
      integer ipbount(2,*),node,nbpta,nbptb,mi(*),i,j,ibount(*)
!
      real*8 vold(0:mi(2),*),xbounact(*),xbount(*)
!
!
!
      do i=nbpta,nbptb
        node=ipbount(1,i)
        vold(0,node)=0.d0
        do j=ipbount(2,i)+1,ipbount(2,i+1)
          vold(0,node)=vold(0,node)+xbounact(ibount(j))
        enddo
        vold(0,node)=vold(0,node)*xbount(i)
      enddo
!
      return
      end
