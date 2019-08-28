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
      subroutine calcttfaext(nfaext,vfa,shcon,nshcon,ielmatf,&
        ntmat_,mi,physcon,ttfa,ifaext,ielfa)
      !
      !     calculation of the total temperature at the external faces
      !     based on the primary variables:
      !
      !     Tt=T+v**2/(2*cp)
      !
      implicit none
      !
      integer nfaext,i,j,imat,ntmat_,mi(*),ielmatf(mi(3),*),&
        nshcon(2,*),ifaext(*),ielfa(4,*)
      !
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),&
        cp,physcon(*),ttfa(*)
      !
      intent(in) nfaext,shcon,nshcon,ielmatf,ntmat_,mi,physcon,ielfa
      !
      intent(inout) vfa,ttfa
      !
      ! $omp parallel default(none)
      ! $omp& shared(nfaext,vfa,ielmatf,ntmat_,shcon,nshcon,physcon,ttfa,ielfa,
      ! $omp&        ifaext)
      ! $omp& private(i,j,t1l,imat,cp)
      !
      !     external face values
      !
      ! $omp do
      do j=1,nfaext
         i=ifaext(j)
         t1l=vfa(0,i)
         imat=ielmatf(1,ielfa(1,i))
         !
         !     heat capacity at constant volume
         !
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,&
              physcon)
         !
         ttfa(i)=t1l+(vfa(1,i)*vfa(1,i)+&
              vfa(2,i)*vfa(2,i)+&
              vfa(3,i)*vfa(3,i))/(2.d0*cp)
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
      
