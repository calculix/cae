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
!     completing hel:
!
!     at the start of the subroutine: rhs of conservation of momentum
!                                     without pressure contribution
!     at the end of the subroutine: neighboring velocity terms subtracted
!
      subroutine complete_hel_cyclic_blk(vel,hel,auv6,c,ipnei,neiel,&
               neifa,ifatie,nef)
      !
      implicit none
      !
      integer i,j,k,indexf,ipnei(*),neiel(*),iel,nef,&
        ifa,neifa(*),ifatie(*)
      !
      real*8 hel(3,*),vel(nef,0:7),auv6(6,*),c(3,3)
      !
      do i=1,nef
         indexf=ipnei(i)
         do j=1,6
            indexf=indexf+1
            iel=neiel(indexf)
            ifa=neifa(indexf)
            if(iel.eq.0) cycle
            !
            if(ifatie(ifa).eq.0) then
               !
               !              no cyclic symmetry face
               !
               do k=1,3
                  hel(k,i)=hel(k,i)-auv6(j,i)*vel(iel,k)
               enddo
            elseif(ifatie(ifa).gt.0) then
               do k=1,3
                  hel(k,i)=hel(k,i)-auv6(j,i)*(c(k,1)*vel(iel,1)+&
                        c(k,2)*vel(iel,2)+c(k,3)*vel(iel,3))
               enddo
            else
               do k=1,3
                  hel(k,i)=hel(k,i)-auv6(j,i)*(c(1,k)*vel(iel,1)+&
                        c(2,k)*vel(iel,2)+c(3,k)*vel(iel,3))
               enddo
            endif
         enddo
      enddo
      !
      return
      end
