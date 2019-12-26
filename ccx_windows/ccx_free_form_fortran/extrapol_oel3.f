!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
      subroutine extrapol_oel3(ielfa,xrlfa,icyclic,ifatie,gradofa,&
        gradoel,c,ipnei,xxi,nfacea,nfaceb,ncfd)
      !
      !     interpolate/extrapolate the turbulent frequency gradient
      !     from the center of the elements to the center of the faces
      !
      implicit none
      !
      integer ielfa(4,*),icyclic,ifatie(*),ipnei(*),nfacea,nfaceb,&
        iel1,iel2,i,l,indexf,ncfd
      !
      real*8 xrlfa(3,*),gradofa(3,*),gradoel(3,*),c(3,3),xxi(3,*),&
        gradnor,xl1,xl2
      !
      intent(in) ielfa,xrlfa,icyclic,ifatie,&
        gradoel,c,ipnei,xxi,nfacea,nfaceb,ncfd
      !
      intent(inout) gradofa
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
                  gradofa(l,i)=xl1*gradoel(l,iel1)+&
                       xl2*gradoel(l,iel2)
               enddo
            elseif(ifatie(i).gt.0) then
               do l=1,ncfd
                  gradofa(l,i)=xl1*gradoel(l,iel1)+xl2*&
                        (gradoel(1,iel2)*c(l,1)+&
                         gradoel(2,iel2)*c(l,2)+&
                         gradoel(3,iel2)*c(l,3))
               enddo
            else
               do l=1,ncfd
                  gradofa(l,i)=xl1*gradoel(l,iel1)+xl2*&
                        (gradoel(1,iel2)*c(1,l)+&
                         gradoel(2,iel2)*c(2,l)+&
                         gradoel(3,iel2)*c(3,l))
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
            do l=1,ncfd
               gradofa(l,i)=xl1*gradoel(l,iel1)+&
                    xrlfa(3,i)*gradoel(l,abs(ielfa(3,i)))
            enddo
         else
            !
            !           boundary face; zero gradient
            !
            indexf=ipnei(iel1)+ielfa(4,i)
            gradnor=gradoel(1,iel1)*xxi(1,indexf)+&
                    gradoel(2,iel1)*xxi(2,indexf)+&
                    gradoel(3,iel1)*xxi(3,indexf)
            do l=1,ncfd
                  gradofa(l,i)=gradoel(l,iel1)&
                              -gradnor*xxi(l,indexf)
            enddo
         endif
      enddo
      !
      return
      end
