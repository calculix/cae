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
      subroutine extrapol_vel1(ielfa,xrlfa,icyclic,ifatie,vfap,vel,
     &  c,ifabou,xbounact,ipnei,xxn,nef,nfacea,nfaceb,ncfd)
!
!     inter/extrapolation of v at the center of the elements
!     to the center of the faces
!
      implicit none
!
      integer ielfa(4,*),icyclic,ifatie(*),ifabou(*),ipnei(*),nfacea,
     &  nfaceb,i,j,ipointer,iel1,iel2,iel3,indexf,nef,ncfd
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),c(3,3),xbounact(*),
     &  xxn(3,*),xl1,xl2,dd
!
!
!
!     initialization of the facial velocities
!
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements: interpolation
!
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do j=1,ncfd
                  vfap(j,i)=xl1*vel(iel1,j)+xl2*vel(iel2,j)
               enddo
            elseif(ifatie(i).gt.0) then
               do j=1,ncfd
                  vfap(j,i)=xl1*vel(iel1,j)+xl2*
     &              (c(j,1)*vel(iel2,1)
     &              +c(j,2)*vel(iel2,2)
     &              +c(j,3)*vel(iel2,3))
               enddo
            else
               do j=1,ncfd
                  vfap(j,i)=xl1*vel(iel1,j)+xl2*
     &              (c(1,j)*vel(iel2,1)
     &              +c(2,j)*vel(iel2,2)
     &              +c(3,j)*vel(iel2,3))
               enddo
            endif
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; no zero gradient
!
            iel3=ielfa(3,i)
            ipointer=-iel2
!
!           global x-direction
!     
            if(ifabou(ipointer+1).gt.0) then
!     
!              v_1 given
!     
               vfap(1,i)=xbounact(ifabou(ipointer+1))
            else
!     
!              extrapolation
!     
               vfap(1,i)=xl1*vel(iel1,1)+xrlfa(3,i)*vel(iel3,1)
            endif
!     
!           global y-direction
!     
            if(ifabou(ipointer+2).gt.0) then
!     
!              v_2 given
!     
               vfap(2,i)=xbounact(ifabou(ipointer+2))
            else
!     
!              extrapolation
!     
               vfap(2,i)=xl1*vel(iel1,2)+xrlfa(3,i)*vel(iel3,2)
            endif
!     
!           global z-direction
!     
            if(ifabou(ipointer+3).gt.0) then
!     
!              v_3 given
!     
               vfap(3,i)=xbounact(ifabou(ipointer+3))
            else
!     
!              extrapolation
!     
               vfap(3,i)=xl1*vel(iel1,3)+xrlfa(3,i)*vel(iel3,3)
            endif
!     
!           correction for sliding boundary conditions        
!     
            if(ifabou(ipointer+5).lt.0) then
               indexf=ipnei(iel1)+ielfa(4,i)
               dd=vfap(1,i)*xxn(1,indexf)+
     &              vfap(2,i)*xxn(2,indexf)+
     &              vfap(3,i)*xxn(3,indexf)
               do j=1,ncfd
                  vfap(j,i)=vfap(j,i)-dd*xxn(j,indexf)
               enddo
            endif
!     
         else
!     
!           boundary face; zero gradient
!     
            do j=1,ncfd
               vfap(j,i)=vel(iel1,j)
            enddo
         endif
      enddo
!
      return
      end
