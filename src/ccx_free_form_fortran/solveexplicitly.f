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
      subroutine solveexplicitly(nef,vel,bv,auv,ipnei,neiel,nflnei)
      !
      !     explicitly solving the momentum equations
      !     the solution is stored in the rhs
      !
      implicit none
      !
      integer i,j,indexf,ipnei(*),iel,neiel(*),nflnei,nef
      !
      real*8 vel(nef,0:7),bv(nef,3),auv(*)
      !
      ! $omp parallel default(none)
      ! $omp& shared(nef,vel,bv,neiel,nflnei,auv,ipnei)
      ! $omp& private(i,j,indexf,iel)
      ! $omp do
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            iel=neiel(indexf)
            do j=1,3
               bv(i,j)=bv(i,j)-auv(indexf)*vel(iel,j)
            enddo
         enddo
         do j=1,3
            bv(i,j)=bv(i,j)/auv(nflnei+i)
         enddo
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
