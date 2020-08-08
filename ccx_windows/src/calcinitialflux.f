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
      subroutine calcinitialflux(area,vfa,xxna,
     &  ipnei,nef,neifa,lakonf,flux)
!
!     correction of v due to the balance of mass
!     the correction is in normal direction to the face
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,j,indexf,ipnei(*),ifa,nef,neifa(*),numfaces
!
      real*8 area(*),vfa(0:7,*),xxna(3,*),flux(*)
!
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            flux(indexf)=vfa(5,ifa)*
     &               (vfa(1,ifa)*xxna(1,indexf)+
     &                vfa(2,ifa)*xxna(2,indexf)+
     &                vfa(3,ifa)*xxna(3,indexf))
         enddo
      enddo
!  
      return
      end
