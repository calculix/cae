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
      subroutine correctrhofacomp(nface,vfa,shcon,ielmat,ntmat_,&
        mi,ielfa,ipnei,vel,nef,flux,gradpel,gradtel,xxj,betam,&
        xlet)
      !
      !     calculation of the density at the face centers
      !     (compressible fluids)
      !
      implicit none
      !
      integer nface,i,j,imat,ntmat_,mi(*),ipnei(*),nef,iel1,iel2,&
        ielmat(mi(3),*),ielfa(4,*),indexf
      !
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),vel(nef,0:7),flux(*),&
        r,gradpel(3,*),gradtel(3,*),xxj(3,*),gamma,betam,phic,vud,&
        vcd,xlet(*)
      !
      ! $omp parallel default(none)
      ! $omp& shared(nface,vfa,ielmat,ielfa,shcon,ipnei,flux,gradpel,gradtel,
      ! $omp&        xxj,betam,xlet,vel)
      ! $omp& private(i,j,t1l,imat,iel1,iel2,indexf,vcd,vud,r,gamma,phic)
      ! $omp do
      do i=1,nface
         t1l=vfa(0,i)
         !
         !        take the material of the first adjacent element
         !
         imat=ielmat(1,ielfa(1,i))
         r=shcon(3,1,imat)
         iel2=ielfa(2,i)
         !
         !        faces with only one neighbor need not be treated
         !
         if(iel2.le.0) then
            vfa(5,i)=vfa(5,i)+0.2d0*vfa(4,i)/(r*t1l)
            cycle
         endif
         !
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
         !
         vcd=vel(iel2,5)-vel(iel1,5)
         if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,5))) vcd=0.d0
         !
         if(flux(indexf).ge.0.d0) then
            vud=2.d0*xlet(indexf)*&
                 ((gradpel(1,iel1)-vel(iel1,5)*r*gradtel(1,iel1))&
                   *xxj(1,indexf)+&
                  (gradpel(2,iel1)-vel(iel1,5)*r*gradtel(2,iel1))&
                   *xxj(2,indexf)+&
                  (gradpel(3,iel1)-vel(iel1,5)*r*gradtel(3,iel1))&
                   *xxj(3,indexf))/(r*vel(iel1,0))
         else
            vud=2.d0*xlet(indexf)*&
                 ((gradpel(1,iel2)-vel(iel2,5)*r*gradtel(1,iel2))&
                   *xxj(1,indexf)+&
                  (gradpel(2,iel2)-vel(iel2,5)*r*gradtel(2,iel2))&
                   *xxj(2,indexf)+&
                  (gradpel(3,iel2)-vel(iel2,5)*r*gradtel(3,iel2))&
                   *xxj(3,indexf))/(r*vel(iel2,0))
         endif
         !
         if(dabs(vud).lt.1.d-20) then
            gamma=0.d0
         else
            !
            phic=1.d0-vcd/vud
            !
            if(phic.ge.1.d0) then
               gamma=0.d0
            elseif(phic.le.0.d0) then
               gamma=0.d0
            elseif(betam.le.phic) then
               gamma=1.d0
            else
               gamma=phic/betam
            endif
         endif
         !
         !        original value plus 20 % of a
         !        mixture of central difference and upwind difference
         !        of the correction
         !
         !             vfa(5,i)=vfa(5,i)+0.2d0*
         !      &              vfa(4,i)/(r*t1l)
         if(flux(indexf).ge.0.d0) then
            vfa(5,i)=vfa(5,i)+0.2d0*&
                    (gamma*vfa(4,i)+(1.d0-gamma)*vel(iel1,4))/(r*t1l)
         else
            vfa(5,i)=vfa(5,i)+0.2d0*&
                    (gamma*vfa(4,i)+(1.d0-gamma)*vel(iel2,4))/(r*t1l)
         endif
      !
      enddo
      ! $omp end do
      ! $omp end parallel
      !
      return
      end
