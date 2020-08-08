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
      subroutine extrapol_kel4(ielfa,vfa,vfap,gradkfa,rf,ifabou,
     &  ipnei,vel,xxi,xle,nef,inlet,nfacea,nfaceb)
!
!     extrapolation of turbulent kinetic energy values to the faces (taking the
!     skewness of the elements into account)
!
!     for turbulent kinetic energy calculations the external faces have
!     as boundary condition either a specified turbulent kinetic energy or
!     a specified flux. If the user did not apply any of these,
!     a zero specified flux is implicitly assumed
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  indexf,ipnei(*),ipointer,inlet(*)
!
      real*8 vfap(0:7,*),vel(nef,0:7),vfa(0:7,*),rf(3,*),gradkfa(3,*),
     &  xle(*),xxi(3,*)
!
!
!
!     Moukalled et al. p 279
!
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!     
!        interpolation
!     
            vfa(6,i)=vfap(6,i)+gradkfa(1,i)*rf(1,i)
     &           +gradkfa(2,i)*rf(2,i)
     &           +gradkfa(3,i)*rf(3,i)
         elseif(ielfa(3,i).gt.0) then
!     
!           no implicit zero gradient
!     
            ipointer=-iel2
!     
            if(ifabou(ipointer+5).gt.0) then
!     
!              wall: kinetic turbulent energy known
!     
               vfa(6,i)=vfap(6,i)
            elseif((inlet(i).eq.1).or.
     &              (ifabou(ipointer+5).lt.0)) then
!     
!              inlet or sliding conditions: kinetic turbulent energy known
!     
               vfa(6,i)=vfap(6,i)
            else
!     
!              turbulent kinetic energy is not given
!     
               vfa(6,i)=vfap(6,i)+gradkfa(1,i)*rf(1,i)+
     &              gradkfa(2,i)*rf(2,i)+
     &              gradkfa(3,i)*rf(3,i)
            endif
         else
!     
!           zero gradient in i-direction
!     
            vfa(6,i)=vel(iel1,6)
         endif
      enddo
!     
      return
      end
