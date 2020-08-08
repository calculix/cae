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
      subroutine calcumel(nef,vel,shcon,nshcon,ielmat,ntmat_,
     &  ithermal,mi,umel)
!
!     calculation of the dynamic viscosity in the element centers
!
      implicit none
!
      integer nef,i,nshcon(*),imat,ithermal(*),ntmat_,mi(*),
     &  ielmat(mi(3),*)
!
      real*8 t1l,vel(nef,0:7),dvi,shcon(0:3,ntmat_,*),umel(*)
!     
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmat(1,i)
         call materialdata_dvi(shcon,nshcon,imat,dvi,t1l,ntmat_,
     &            ithermal)
         umel(i)=dvi
      enddo
!            
      return
      end
