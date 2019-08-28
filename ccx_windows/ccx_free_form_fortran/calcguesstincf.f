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
      subroutine calcguesstincf(nface,dmin,vfa,umfa,cvfa,hcfa,ithermal,&
              tincfguess,compressible)
      !
      !     calculates a guess for tincf based on the minimum of:
      !       dmin/(local velocity)
      !       density*dmin**2/(2*dynamic_viscosity)
      !
      !       dmin is the smallest edge length of the mesh
      !
      implicit none
      !
      integer nface,i,ithermal,compressible
      !
      real*8 vfa(0:7,*),umax,dmin,umfa(*),tincfguess,cvfa(*),hcfa(*)
      !
      tincfguess=1.d30
      do i=1,nface
         !
         !        convection
         !
         umax=dsqrt(vfa(1,i)*vfa(1,i)+&
                    vfa(2,i)*vfa(2,i)+&
                    vfa(3,i)*vfa(3,i))
         if(umax.gt.1.d-30) tincfguess=min(tincfguess,dmin/umax)
         !          write(*,*) 'calcguesstincf ',umax,tincfguess,dmin
         !
         !        viscous diffusion
         !
         if((umfa(i).gt.0.d0).and.(vfa(5,i).gt.0.d0)) then
            tincfguess=min(tincfguess,vfa(5,i)*dmin*dmin/&
                 (2.d0*umfa(i)))
         endif
         !          write(*,*) 'calcguesstincf ',tincfguess,vfa(5,i),dmin,umfa(i)
         !
         !        thermal diffusion
         !
         if(ithermal.gt.0) then
            if((hcfa(i).gt.0.d0).and.(cvfa(i).gt.0.d0).and.&
               (vfa(5,i).gt.0.d0)) then
               tincfguess=min(tincfguess,vfa(5,i)*cvfa(i)*dmin*dmin/&
                    (2.d0*hcfa(i)))
            endif
         !             write(*,*) 'calcguesstincf ',tincfguess,vfa(5,i),cvfa(i),
         !      &            dmin,hcfa(i)
         endif
      enddo
      !       tincfguess=tincfguess/5.d0
      !       write(*,*) 'calcguesstincf ',tincfguess
      !
      return
      end
