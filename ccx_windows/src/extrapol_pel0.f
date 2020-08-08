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
      subroutine extrapol_pel0(vel,nef,nkfa,
     &     nkfb,vold,ikf,ratio,konl,mi)
!
!     interpolation of the element center values of the pressure to
!     the nodes
!     
      implicit none
!
      integer i,nef,nkfa,nkfb,ikf(*),konl(4,*),mi(*)
!
      real*8 vel(nef,0:7),ratio(4,*),vold(0:mi(2),*)
!
!
!
      do i=nkfa,nkfb
        vold(4,ikf(i))=ratio(1,i)*vel(konl(1,i),4)+
     &       ratio(2,i)*vel(konl(2,i),4)+
     &       ratio(3,i)*vel(konl(3,i),4)+
     &       ratio(4,i)*vel(konl(4,i),4)
      enddo
!
      return
      end
