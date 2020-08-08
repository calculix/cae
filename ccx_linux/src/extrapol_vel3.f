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
      subroutine extrapol_vel3(ielfa,xrlfa,icyclic,ifatie,gradvfa,
     &  gradvel,c,ipnei,xxi,nfacea,nfaceb,ncfd)
!
!     interpolate/extrapolate the velocity gradient from the
!     center of the elements to the center of the faces
!
      implicit none
!
      integer ielfa(4,*),icyclic,ifatie(*),ipnei(*),nfacea,nfaceb,
     &  i,k,l,indexf,iel1,iel2,ncfd
!
      real*8 xrlfa(3,*),gradvfa(3,3,*),gradvel(3,3,*),c(3,3),xxi(3,*),
     &  xl1,xl2,gradnor
!
!
! 
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements
!
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do k=1,ncfd
                  do l=1,ncfd
                     gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+
     &                    xl2*gradvel(k,l,iel2)
                  enddo
               enddo
            elseif(ifatie(i).gt.0) then
               do k=1,ncfd
                  do l=1,ncfd
                     gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+xl2*
     &                   (c(k,1)*gradvel(1,1,iel2)*c(l,1)+
     &                    c(k,1)*gradvel(1,2,iel2)*c(l,2)+
     &                    c(k,1)*gradvel(1,3,iel2)*c(l,3)+
     &                    c(k,2)*gradvel(2,1,iel2)*c(l,1)+
     &                    c(k,2)*gradvel(2,2,iel2)*c(l,2)+
     &                    c(k,2)*gradvel(2,3,iel2)*c(l,3)+
     &                    c(k,3)*gradvel(3,1,iel2)*c(l,1)+
     &                    c(k,3)*gradvel(3,2,iel2)*c(l,2)+
     &                    c(k,3)*gradvel(3,3,iel2)*c(l,3))
                  enddo
               enddo
            else
               do k=1,ncfd
                  do l=1,ncfd
                     gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+xl2*
     &                   (c(1,k)*gradvel(1,1,iel2)*c(1,l)+
     &                    c(1,k)*gradvel(1,2,iel2)*c(2,l)+
     &                    c(1,k)*gradvel(1,3,iel2)*c(3,l)+
     &                    c(2,k)*gradvel(2,1,iel2)*c(1,l)+
     &                    c(2,k)*gradvel(2,2,iel2)*c(2,l)+
     &                    c(2,k)*gradvel(2,3,iel2)*c(3,l)+
     &                    c(3,k)*gradvel(3,1,iel2)*c(1,l)+
     &                    c(3,k)*gradvel(3,2,iel2)*c(2,l)+
     &                    c(3,k)*gradvel(3,3,iel2)*c(3,l))
                  enddo
               enddo
            endif
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; no zero gradient
!
            do k=1,ncfd
               do l=1,ncfd
                  gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+
     &                           xrlfa(3,i)*gradvel(k,l,ielfa(3,i))
               enddo
            enddo
         else
!
!           boundary face; zero gradient in i-direction
!
            indexf=ipnei(iel1)+ielfa(4,i)
            do k=1,ncfd
               gradnor=gradvel(k,1,iel1)*xxi(1,indexf)
     &                +gradvel(k,2,iel1)*xxi(2,indexf)
     &                +gradvel(k,3,iel1)*xxi(3,indexf)
               do l=1,ncfd
                  gradvfa(k,l,i)=gradvel(k,l,iel1)
     &                          -gradnor*xxi(l,indexf)
               enddo
            enddo
         endif
      enddo
c      do i=1,120
c         write(*,*) 'extrapol_vel3',i,gradvfa(1,1,i),gradvfa(1,2,i)
c      enddo
!
      return
      end
