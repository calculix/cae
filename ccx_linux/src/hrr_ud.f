!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine hrr_ud(vfa,shcon,ielmat,ntmat_,
     &  mi,ielfa,ipnei,vel,nef,flux,nfacea,nfaceb,xxi,xle,
     &  gradpel,gradtel,neij)
!
!     calculation of the density at the face centers
!     (compressible fluids)
!
!     facial temperature and pressure is only used for external
!     faces
!
      implicit none
!
      integer i,j,imat,ntmat_,mi(*),ipnei(*),nef,iel1,iel2,
     &  ielmat(mi(3),*),ielfa(4,*),indexf,nfacea,nfaceb,neij(*)
!
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),vel(nef,0:7),flux(*),
     &  r,dd,xxv(3),xxi(3,*),qp(3),p,t,xle(*),gradpel(3,*),
     &  gradtel(3,*)
!
!
!     
      do i=nfacea,nfaceb
!
         iel2=ielfa(2,i)
!     
!        take the material of the first adjacent element
!     
         imat=ielmat(1,ielfa(1,i))
         r=shcon(3,1,imat)
!
!        no neighbor
!
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
            vfa(5,i)=vfa(4,i)/(r*vfa(0,i))
            cycle
         endif
!
!        neighbor
!
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
!           outflow 
!
c            vfa(5,i)=vel(iel1,5)
            vfa(5,i)=vel(iel1,4)/(r*vfa(0,i))
         elseif(iel2.gt.0) then
!
!           inflow && neighbor
!
c            vfa(5,i)=vel(iel2,5)
            vfa(5,i)=vel(iel2,4)/(r*vfa(0,i))
         endif
!
      enddo
!            
      return
      end
