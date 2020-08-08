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
      subroutine extrapol_vel2(ipnei,neifa,vfa,area,xxna,volume,gradvel,
     &  nefa,nefb,ncfd)
!
!     calculate the gradient of the velocities at the center of
!     the elements
!
      implicit none
!
      integer ipnei(*),neifa(*),nefa,nefb,i,k,l,ifa,indexf,ncfd
!
      real*8 vfa(0:7,*),area(*),xxna(3,*),volume(*),gradvel(3,3,*)
!
!
!
      do i=nefa,nefb
!
!           initialization
!
         do k=1,ncfd
            do l=1,ncfd
               gradvel(k,l,i)=0.d0
            enddo
         enddo
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do k=1,ncfd
               do l=1,ncfd
                  gradvel(k,l,i)=gradvel(k,l,i)+
     &                 vfa(k,ifa)*xxna(l,indexf)
               enddo
            enddo
         enddo
!     
!     dividing by the volume of the element
!     
         do k=1,ncfd
            do l=1,ncfd
               gradvel(k,l,i)=gradvel(k,l,i)/volume(i)
            enddo
         enddo
      enddo
c      do i=1,120
c         write(*,*) 'extrapol_vel2',i,gradvel(1,1,i),gradvel(1,2,i)
c      enddo
!
      return
      end
