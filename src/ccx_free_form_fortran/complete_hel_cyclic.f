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
      subroutine complete_hel_cyclic(nef,bv,hel,adv,auv,jq,irow,&
        ipnei,neiel,ifatie,c,lakonf,neifa,nzs)
      !
      implicit none
      !
      character*8 lakonf(*)
      !
      integer irow(*),nef,nzs,j,k,l,jdof1,jq(*),ifa,neifa(*),&
        indexf,ipnei(*),neiel(1),ifatie(*),iel,i
      !
      real*8 hel(3,*),bv(nef,3),auv(*),adv(*),c(3,3)
      !
      !     off-diagonal terms
      !
      ! $omp parallel default(none)
      ! $omp& shared(nef,ipnei,neiel,neifa,ifatie,hel,auv,bv,c)
      ! $omp& private(i,indexf,iel,ifa,k)
      ! $omp do
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            !
            iel=neiel(indexf)
            if(iel.eq.0) cycle
            ifa=neifa(indexf)
            !
            if(ifatie(ifa).eq.0) then
               do k=1,3
                  hel(k,i)=hel(k,i)-auv(indexf)*bv(iel,k)
               enddo
            elseif(ifatie(ifa).gt.0) then
                  !                if(i.gt.iel) then
                  do k=1,3
                     hel(k,i)=hel(k,i)-auv(indexf)*&
                          (c(k,1)*bv(iel,1)+c(k,2)*bv(iel,2)&
                          +c(k,3)*bv(iel,3))
                  enddo
            !                else
            !                   do k=1,3
            !                      hel(k,i)=hel(k,i)-auv(indexf)*
            !      &                    (c(1,k)*bv(iel,1)+c(2,k)*bv(iel,2)
            !      &                    +c(3,k)*bv(iel,3))
            !                   enddo
            !                endif
            else
                  !                if(i.gt.iel) then
                  do k=1,3
                     hel(k,i)=hel(k,i)-auv(indexf)*&
                          (c(1,k)*bv(iel,1)+c(2,k)*bv(iel,2)&
                          +c(3,k)*bv(iel,3))
                  enddo
            !                else
            !                   do k=1,3
            !                      hel(k,i)=hel(k,i)-auv(indexf)*
            !      &                    (c(k,1)*bv(iel,1)+c(k,2)*bv(iel,2)
            !      &                    +c(k,3)*bv(iel,3))
            !                   enddo
            !                endif
            endif
         enddo
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
