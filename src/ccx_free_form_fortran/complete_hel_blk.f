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
      subroutine complete_hel_blk(vel,hel,auv6,ipnei,neiel,nef,&
              nactdohinv)
      !
      implicit none
      !
      integer i,j,k,indexf,ipnei(*),neiel(*),iel,nef,nactdohinv(*)
      !
      real*8 hel(3,*),vel(nef,0:7),auv6(6,*)
      !
      do i=1,nef
         indexf=ipnei(i)
         do j=1,6
            indexf=indexf+1
            iel=neiel(indexf)
            if(iel.eq.0) cycle
            do k=1,3
               hel(k,i)=hel(k,i)-auv6(j,i)*vel(iel,k)
            !                write(*,*) 'complete_hel_blk ',nactdohinv(i),
            !      &                nactdohinv(iel),k,
            !      &                -auv6(j,i)
            enddo
         enddo
      enddo
      !
      !       do j=1,nef
      !          write(*,*) 'complete_hel ',j,(hel(k,j),k=1,3)
      !       enddo
      !
      return
      end
