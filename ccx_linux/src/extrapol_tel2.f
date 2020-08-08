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
      subroutine extrapol_tel2(ipnei,neifa,vfa,area,xxna,volume,gradtel,
     &  nefa,nefb,ncfd)
!
!     calculate the gradient of the temperature at the center of
!     the elements
!
      implicit none
!
      integer ipnei(*),neifa(*),nefa,nefb,ifa,i,l,indexf,ncfd
!
      real*8 vfa(0:7,*),area(*),xxna(3,*),volume(*),gradtel(3,*)
!
!
!
      do i=nefa,nefb
!
!        initialization
!     
         do l=1,ncfd
            gradtel(l,i)=0.d0
         enddo
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do l=1,ncfd
               gradtel(l,i)=gradtel(l,i)+
     &              vfa(0,ifa)*xxna(l,indexf)
            enddo
         enddo
!     
!        dividing by the volume of the element
!     
         do l=1,ncfd
            gradtel(l,i)=gradtel(l,i)/volume(i)
         enddo
      enddo
!            
      return
      end
