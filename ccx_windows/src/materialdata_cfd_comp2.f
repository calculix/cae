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
      subroutine materialdata_cfd_comp2(shcon,nshcon,ielmatf,
     &  ntmat_,mi,vfa,cocon,ncocon,physcon,cvfa,ithermal,
     &  umfa,ielfa,hcfa,nfacea,nfaceb)
!
!     calculation of material properties at elements centers and
!     face centers (compressible fluids)
!
      implicit none
!
      integer i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal(*),
     &  nshcon(2,*),ncocon(2,*),ielfa(4,*),nfacea,nfaceb
!
      real*8 t1l,shcon(0:3,ntmat_,*),vfa(0:7,*),
     &  cp,cocon(0:6,ntmat_,*),physcon(*),cvfa(*),umfa(*),
     &  hcfa(*)
!
!
!
!     facial values
!
      do i=nfacea,nfaceb
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmatf(1,ielfa(1,i))
!
!        density
!
c         vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*t1l)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvfa(i)=cp-shcon(3,1,imat)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umfa(i),t1l,ntmat_,
     &            ithermal)
!
!        heat conduction
!
         call materialdata_cond(imat,ntmat_,t1l,cocon,ncocon,hcfa(i))
      enddo
!            
      return
      end
