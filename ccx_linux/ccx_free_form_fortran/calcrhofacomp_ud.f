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
      subroutine calcrhofacomp_ud(vfa,shcon,ielmat,ntmat_,&
        mi,ielfa,ipnei,vel,nef,flux,nfacea,nfaceb,xxi,xle,&
        gradpel,gradtel,neij)
      !
      !     calculation of the density at the face centers
      !     (compressible fluids)
      !
      !     facial temperature and pressure is only used for external
      !     faces
      !
      implicit none
      !
      integer i,j,imat,ntmat_,mi(*),ipnei(*),nef,iel1,iel2,&
        ielmat(mi(3),*),ielfa(4,*),indexf,nfacea,nfaceb,neij(*)
      !
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),vel(nef,0:7),flux(*),&
        r,dd,xxv(3),xxi(3,*),qp(3),p,t,xle(*),gradpel(3,*),&
        gradtel(3,*)
      !
      intent(in) shcon,ielmat,ntmat_,mi,ielfa,ipnei,vel,nef,flux,&
        nfacea,nfaceb
      !
      intent(inout) vfa
      !
      do i=nfacea,nfaceb
         !
         iel2=ielfa(2,i)
         !
         !        faces with only one neighbor need not be treated
         !        unless outlet
         !
         !          if((iel2.le.0).and.(ielfa(3,i).ge.0)) then
         if(iel2.le.0) then
            t1l=vfa(0,i)
            !
            !           take the material of the first adjacent element
            !
            imat=ielmat(1,ielfa(1,i))
            r=shcon(3,1,imat)
            !
            !           specific gas constant
            !
            vfa(5,i)=vfa(4,i)/(r*t1l)
            cycle
         endif
         !
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
         !
         if(flux(indexf).ge.0.d0) then
            !
            !           outflow && (neighbor || outlet)
            !
            vfa(5,i)=vel(iel1,5)
         elseif(iel2.gt.0) then
            !
            !           inflow && neighbor
            !
            vfa(5,i)=vel(iel2,5)
         endif
      !
      enddo
      !
      return
      end
