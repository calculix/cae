!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine extrapol_kel2(ipnei,neifa,vfap,area,xxn,volume,gradkel,&
        nefa,nefb)
      !
      !     calculate the gradient of the temperature at the center of
      !     the elements
      !
      implicit none
      !
      integer ipnei(*),neifa(*),nefa,nefb,ifa,i,l,indexf
      !
      real*8 vfap(0:7,*),area(*),xxn(3,*),volume(*),gradkel(3,*)
      !
      intent(in) ipnei,neifa,vfap,area,xxn,volume,&
        nefa,nefb
      !
      intent(inout) gradkel
      !
      do i=nefa,nefb
         !
         !        initialization
         !
         do l=1,3
            gradkel(l,i)=0.d0
         enddo
         !
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do l=1,3
               gradkel(l,i)=gradkel(l,i)+&
                    vfap(6,ifa)*area(ifa)*xxn(l,indexf)
            enddo
         enddo
         !
         !        dividing by the volume of the element
         !
         do l=1,3
            gradkel(l,i)=gradkel(l,i)/volume(i)
         enddo
      enddo
      !
      return
      end
