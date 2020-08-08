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
      subroutine calcstressheatflux(sti,umel,gradvel,qfx,hcel,gradtel,
     &  nef,isti,iqfx,mi)
!
!     calculation of the stress and the heat flux for output
!     in CFD calculations
!
      implicit none
!
      integer nef,isti,iqfx,mi(*),i
!
      real*8 sti(6,mi(1),*),qfx(3,mi(1),*),umel(*),hcel(*),
     &  gradvel(3,3,*),gradtel(3,*)
!
      do i=1,nef
         if(isti.gt.0) then
            sti(1,1,i)=2.d0*umel(i)*gradvel(1,1,i)
            sti(2,1,i)=2.d0*umel(i)*gradvel(2,2,i)
            sti(3,1,i)=2.d0*umel(i)*gradvel(3,3,i)
            sti(4,1,i)=umel(i)*(gradvel(1,2,i)+gradvel(2,1,i))
            sti(5,1,i)=umel(i)*(gradvel(1,3,i)+gradvel(3,1,i))
            sti(6,1,i)=umel(i)*(gradvel(2,3,i)+gradvel(3,2,i))
         endif
         if(iqfx.gt.0) then
            qfx(1,1,i)=-hcel(i)*gradtel(1,i)
            qfx(2,1,i)=-hcel(i)*gradtel(2,i)
            qfx(3,1,i)=-hcel(i)*gradtel(3,i)
         endif
      enddo
!
      return
      end
