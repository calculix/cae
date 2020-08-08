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
      subroutine extrapol_dpel01(ipbounp,vold,
     &     mi,nbppa,nbppb)
!
!     application of boundary conditions to the nodes
!     pressure
!     
      implicit none
!
      integer ipbounp(2,*),node,nbppa,nbppb,mi(*),i
!
      real*8 vold(0:mi(2),*)
!
!
!
      do i=nbppa,nbppb
        node=ipbounp(1,i)
        vold(4,node)=0.d0
      enddo
!
      return
      end
