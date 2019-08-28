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
      subroutine norm(vel,velnorm,nef)
      !
      !     calculation of the norm of all field components at
      !     the element centers
      !
      implicit none
      !
      integer i,j,nef
      !
      real*8 velnorm(0:4),vel(nef,0:7)
      !
      ! $omp parallel default(none)
      ! $omp& shared(velnorm,vel,nef)
      ! $omp& private(i,j)
      ! $omp do
      do i=1,nef
         do j=0,4
            velnorm(j)=velnorm(j)+vel(i,j)**2
         enddo
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
