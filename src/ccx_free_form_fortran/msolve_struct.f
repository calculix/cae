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
!     matrix preconditioning: used in dgmres.f
!
      subroutine msolve_struct(n,r,z,nelt,ia,ja,a,isym,rwork,iwork)
      !
      implicit none
      !
      integer n,nelt,ia(*),ja(*),isym,iwork(*),i,nd
      !
      real*8 r(*),z(*),a(*),rwork(*)
      !
      ! $omp parallel default(none)
      ! $omp& shared(n,z,r,rwork)
      ! $omp& private(i)
      ! $omp do
      do i=1,n
         z(i)=r(i)*rwork(i)
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
