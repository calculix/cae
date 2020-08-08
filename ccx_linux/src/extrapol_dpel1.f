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
      subroutine extrapol_dpel1(ielfa,xrlfa,vfap,vel,ifabou,
     &  nef,nfacea,nfaceb)
!
!     extrapolating the pressure correction from the element
!     centers to the faces
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  ibou
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),xl1
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
            vfap(4,i)=xl1*vel(iel1,4)+xrlfa(2,i)*vel(iel2,4)
!
         elseif(ielfa(3,i).ne.0) then
!
!           boundary face; more than one layer
!            
            ibou=0
            if(iel2.lt.0) then
               if(ifabou(-iel2+4).gt.0) then
                  ibou=ifabou(-iel2+4)
               endif
            endif
!
            if(ibou.gt.0) then
!
!              pressure boundary condition
!
               vfap(4,i)=0.d0
            else
!
!              extrapolation
!
               vfap(4,i)=xl1*vel(iel1,4)
     &              +xrlfa(3,i)*vel(abs(ielfa(3,i)),4)
           endif
         else
!
!           boundary face; one layer
!
            vfap(4,i)=vel(iel1,4)
         endif
      enddo
!            
      return
      end
