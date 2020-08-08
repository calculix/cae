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
      subroutine extrapol_vel4(ielfa,vfa,vfap,gradvfa,rf,ifabou,
     &  ipnei,xxn,vel,xxi,xle,nef,nfacea,nfaceb,ncfd)
!
!     inter/extrapolation of v at the center of the elements
!     to the center of the faces: taking the skewness of the 
!     elements into account (using the gradient at the center
!     of the faces)
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),ipnei(*),nef,nfacea,nfaceb,i,j,
     &  indexf,ipointer,iel1,iel2,ncfd
!
      real*8 vfa(0:7,*),vfap(0:7,*),gradvfa(3,3,*),rf(3,*),xxn(3,*),
     &  vel(nef,0:7),xxi(3,*),xle(*),dd
!
!
!
!        Moukalled et al. p 279
!
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!     
!              face between two elements
!
            do j=1,ncfd
               vfa(j,i)=vfap(j,i)+gradvfa(j,1,i)*rf(1,i)+
     &              gradvfa(j,2,i)*rf(2,i)+
     &              gradvfa(j,3,i)*rf(3,i)
            enddo
         elseif(ielfa(3,i).gt.0) then
!
!              boundary face; no zero gradient
!
            ipointer=-iel2
!     
!           x-direction
!     
            if(ifabou(ipointer+1).gt.0) then
               vfa(1,i)=vfap(1,i)
            else
               vfa(1,i)=vfap(1,i)+gradvfa(1,1,i)*rf(1,i)+
     &              gradvfa(1,2,i)*rf(2,i)+
     &              gradvfa(1,3,i)*rf(3,i)
            endif
!     
!           y-direction
!     
            if(ifabou(ipointer+2).gt.0) then
               vfa(2,i)=vfap(2,i)
            else
               vfa(2,i)=vfap(2,i)+gradvfa(2,1,i)*rf(1,i)+
     &              gradvfa(2,2,i)*rf(2,i)+
     &              gradvfa(2,3,i)*rf(3,i)
            endif
!     
!           z-direction
!     
            if(ifabou(ipointer+3).gt.0) then
               vfa(3,i)=vfap(3,i)
            else
               vfa(3,i)=vfap(3,i)+gradvfa(3,1,i)*rf(1,i)+
     &              gradvfa(3,2,i)*rf(2,i)+
     &              gradvfa(3,3,i)*rf(3,i)
            endif
!     
!           correction for sliding boundary conditions        
!     
            if(ifabou(ipointer+5).lt.0) then
               indexf=ipnei(iel1)+ielfa(4,i)
               dd=vfa(1,i)*xxn(1,indexf)+
     &              vfa(2,i)*xxn(2,indexf)+
     &              vfa(3,i)*xxn(3,indexf)
               do j=1,ncfd
                  vfa(j,i)=vfa(j,i)-dd*xxn(j,indexf)
               enddo
            endif
         else
!     
!           boundary face; zero gradient
!     
c            indexf=ipnei(iel1)+ielfa(4,i)
            do j=1,ncfd
               vfa(j,i)=vel(iel1,j)
c     &              +(gradvfa(j,1,i)*xxi(1,indexf)+
c     &              gradvfa(j,2,i)*xxi(2,indexf)+
c     &              gradvfa(j,3,i)*xxi(3,indexf))*xle(indexf)
            enddo
         endif
      enddo
!
      return
      end
