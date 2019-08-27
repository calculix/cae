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
      subroutine complete_hel(nef,bv,hel,adv,auv,ipnei,neiel,nzs)
      !
      implicit none
      !
      integer neiel(*),nef,nzs,j,k,l,jdof1,ipnei(*),indexf,i,iel
      !
      real*8 hel(3,*),bv(nef,3),auv(*),adv(*)
      !
      !     off-diagonal terms
      !
      ! $omp parallel default(none)
      ! $omp& shared(nef,ipnei,neiel,hel,auv,bv)
      ! $omp& private(i,indexf,iel,k)
      ! $omp do
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            !
            !           neighboring element
            !
            iel=neiel(indexf)
            if(iel.ne.0) then
               do k=1,3
                  hel(k,i)=hel(k,i)-auv(indexf)*bv(iel,k)
               enddo
            endif
         enddo
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
