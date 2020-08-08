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
      subroutine calcrhofa(nface,vfa,rhcon,nrhcon,ielmat,ntmat_,
     &  ithermal,mi,ielfa)
!
!     calculation of the density at the face centers
!     (incompressible fluids)
!
      implicit none
!
      integer nface,i,nrhcon(*),imat,ithermal(*),ntmat_,mi(*),
     &  ielmat(mi(3),*),ielfa(4,*)
!
      real*8 t1l,vfa(0:7,*),rhcon(0:1,ntmat_,*) 
!     
      do i=1,nface
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmat(1,ielfa(1,i))
         call materialdata_rho(rhcon,nrhcon,imat,vfa(5,i),t1l,ntmat_,
     &            ithermal)
      enddo
!            
      return
      end
