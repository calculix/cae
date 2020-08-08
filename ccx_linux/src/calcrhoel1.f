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
      subroutine calcrhoel1(nef,vel,rhcon,nrhcon,ielmat,ntmat_,
     &  ithermal,mi,nefa,nefb)
!
!     calculation of rho in the element centers (incompressible
!     fluids)
!
      implicit none
!
      integer nef,i,nrhcon(*),imat,ithermal(*),ntmat_,mi(*),
     &  ielmat(mi(3),*),nefa,nefb
!
      real*8 t1l,vel(nef,0:7),rhcon(0:1,ntmat_,*)
!
!
!     
      do i=nefa,nefb
         t1l=vel(i,0)
         imat=ielmat(1,i)
         call materialdata_rho(rhcon,nrhcon,imat,vel(i,5),t1l,ntmat_,
     &            ithermal)
      enddo
!            
      return
      end
