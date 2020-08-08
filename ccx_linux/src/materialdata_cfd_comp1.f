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
      subroutine materialdata_cfd_comp1(nef,vel,shcon,nshcon,ielmatf,
     &  ntmat_,mi,cvel,physcon,ithermal,umel,nefa,nefb)
!
!     calculation of material properties at element centers
!     (compressible fluids)
!
      implicit none
!
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal(*),
     &  nshcon(2,*),nefa,nefb
!
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),cvel(*),
     &  cp,physcon(*),umel(*)
!
!
!
!     element (cell) values
!
      do i=nefa,nefb
         t1l=vel(i,0)
         imat=ielmatf(1,i)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvel(i)=cp-shcon(3,1,imat)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umel(i),t1l,ntmat_,
     &            ithermal)
      enddo
!            
      return
      end
