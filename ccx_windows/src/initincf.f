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
      subroutine initincf(nface,hmin,vfa,umfa,cvfa,hcfa,ithermal,
     &        dtimef,compressible)
!
!     calculates a guess for tincf based on the minimum of:
!       hmin/(local velocity)
!       density*hmin**2/(2*dynamic_viscosity)
!
!       hmin is the smallest edge length of the mesh
!
      implicit none
!
      integer nface,i,ithermal(*),compressible
!
      real*8 vfa(0:7,*),umax,hmin,umfa(*),dtimef,cvfa(*),hcfa(*)
!
      dtimef=1.d30
      do i=1,nface
!
!        convection
!
         umax=dsqrt(vfa(1,i)*vfa(1,i)+
     &              vfa(2,i)*vfa(2,i)+
     &              vfa(3,i)*vfa(3,i))
         if(umax.gt.1.d-30) dtimef=min(dtimef,hmin/umax)
c         write(*,*) 'calcguesstincf conv ',dtimef
c!
!        viscous diffusion
!
         if((umfa(i).gt.0.d0).and.(vfa(5,i).gt.0.d0)) then
            dtimef=min(dtimef,vfa(5,i)*hmin*hmin/
     &           (2.d0*umfa(i)))
         endif
c         write(*,*) 'calcguesstincf diff ',dtimef
!
!        thermal diffusion
!
         if(ithermal(1).gt.0) then
            if((hcfa(i).gt.0.d0).and.(cvfa(i).gt.0.d0).and.
     &         (vfa(5,i).gt.0.d0)) then
               dtimef=min(dtimef,vfa(5,i)*cvfa(i)*hmin*hmin/
     &              (2.d0*hcfa(i)))
            endif
         endif
c         write(*,*) 'calcguesstincf ther ',dtimef
      enddo
!     
      return
      end
