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
      subroutine extrapol_oel1(ielfa,xrlfa,vfap,vel,ifabou,
     &  nef,umfa,constant,dy,vfa,inlet,nfacea,nfaceb)
!
!     extrapolation of turbulent kinetic energy values to the faces
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  ipointer,inlet(*)
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),xl1,dy(*),
     &  umfa(*),constant,vfa(0:7,*)
!
!
!     
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements: interpolation
!
            vfap(7,i)=xl1*vel(iel1,7)+xrlfa(2,i)*vel(iel2,7)
!
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; no zero gradient
!
!           iel2=0 is not possible: if iel2=0, there are no
!           boundary conditions on the face, hence it is an
!           exit, which means zero gradient and 
!           ielfa(3,i) <= 0
!     
            ipointer=-iel2
!     
            if(ifabou(ipointer+5).gt.0) then
!     
!              wall: turbulent dissipation rate known
!     
               vfap(7,i)=dy(ifabou(ipointer+5))*umfa(i)/vfa(5,i)
            elseif((inlet(i).eq.1).or.
     &             (ifabou(ipointer+5).lt.0)) then
!     
!              inlet or sliding conditions: kinetic turbulent energy known
!     
               vfap(7,i)=constant*umfa(i)
            else
!
!              extrapolation
!
               vfap(7,i)=xl1*vel(iel1,7)+xrlfa(3,i)*vel(ielfa(3,i),7)
            endif
         else
!
!           boundary face; zero gradient in i-direction
!
            vfap(7,i)=vel(iel1,7)
         endif
      enddo
!            
      return
      end
