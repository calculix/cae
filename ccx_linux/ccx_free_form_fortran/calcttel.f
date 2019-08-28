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
      subroutine calcttel(nef,vel,shcon,nshcon,ielmatf,&
        ntmat_,mi,physcon,ttel)
      !
      !     calculation of material properties at elements centers and
      !     face centers (compressible fluids)
      !
      implicit none
      !
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),&
        nshcon(2,*)
      !
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),&
        cp,physcon(*),&
        ttel(*)
      !
      intent(in) nef,shcon,nshcon,ielmatf,ntmat_,mi,physcon
      !
      intent(inout) vel,ttel
      !
      ! $omp parallel default(none)
      ! $omp& shared(nef,vel,ielmatf,ntmat_,shcon,nshcon,physcon,ttel)
      ! $omp& private(i,t1l,imat,cp)
      !
      !     element (cell) values
      !
      ! $omp do
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmatf(1,i)
         !
         !     heat capacity at constant volume
         !
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,&
              physcon)
         !
         ttel(i)=t1l+(vel(i,1)*vel(i,1)+&
              vel(i,2)*vel(i,2)+&
              vel(i,3)*vel(i,3))/(2.d0*cp)
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
      
