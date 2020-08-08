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
      subroutine extrapol_tel3(ielfa,xrlfa,icyclic,ifatie,gradtfa,
     &  gradtel,c,ipnei,xxi,ifabou,xxn,xload,nfacea,nfaceb,ncfd)
!
!     interpolate/extrapolate the temperature gradient from the
!     center of the elements to the center of the faces
!           
      implicit none
!
      integer ielfa(4,*),icyclic,ifatie(*),ipnei(*),nfacea,nfaceb,
     &  iel1,iel2,i,l,indexf,ifabou(*),ncfd
!
      real*8 xrlfa(3,*),gradtfa(3,*),gradtel(3,*),c(3,3),xxi(3,*),
     &  gradnor,xl1,xl2,q,xxn(3,*),xload(2,*)
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
               do l=1,ncfd
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &                 xl2*gradtel(l,iel2)
               enddo
            elseif(ifatie(i).gt.0) then
               do l=1,ncfd
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+xl2*
     &                  (gradtel(1,iel2)*c(l,1)+
     &                   gradtel(2,iel2)*c(l,2)+
     &                   gradtel(3,iel2)*c(l,3))
               enddo
            else
               do l=1,ncfd
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+xl2*
     &                  (gradtel(1,iel2)*c(1,l)+
     &                   gradtel(2,iel2)*c(2,l)+
     &                   gradtel(3,iel2)*c(3,l))
               enddo
            endif
         elseif(ielfa(3,i).gt.0) then
!
!           the above condition implies iel2!=0
!           if iel2 were zero, no b.c. would apply and
!           ielfa(3,i) would be zero or negative
!     
!           boundary face; no zero gradient
!
            if(ifabou(-iel2).gt.0) then
!
!              temperature given: extrapolate gradient
!
               do l=1,ncfd
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &                 xrlfa(3,i)*gradtel(l,abs(ielfa(3,i)))
               enddo
            else
!
!              facial temperature gradient given, may be nonzero
!
               indexf=ipnei(iel1)+ielfa(4,i)
!
               if(ifabou(-iel2+6).eq.0) then
                  q=0.d0
               else
                  q=xload(1,ifabou(-iel2+6))
               endif
               gradnor=gradtel(1,iel1)*xxn(1,indexf)
     &                +gradtel(2,iel1)*xxn(2,indexf)
     &                +gradtel(3,iel1)*xxn(3,indexf)-q
               do l=1,ncfd
                  gradtfa(l,i)=gradtel(l,iel1)
     &                        -gradnor*xxn(l,indexf)
               enddo
            endif
         else
!     
!           boundary face; zero gradient
!   
            indexf=ipnei(iel1)+ielfa(4,i)
            gradnor=gradtel(1,iel1)*xxi(1,indexf)+
     &              gradtel(2,iel1)*xxi(2,indexf)+
     &              gradtel(3,iel1)*xxi(3,indexf)
            do l=1,ncfd
                  gradtfa(l,i)=gradtel(l,iel1)
     &                        -gradnor*xxi(l,indexf)
            enddo
         endif
      enddo
!            
      return
      end
