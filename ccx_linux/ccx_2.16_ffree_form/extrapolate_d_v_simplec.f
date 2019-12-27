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
      subroutine extrapolate_d_v_simplec(ielfa,xrlfa,adv,advfa,&
                        hfa,icyclic,c,ifatie,vel,nef,volume,&
                        auv,ipnei,nfacea,nfaceb,ncfd)
      !
      !     inter/extrapolation of volume/adv at the center of the elements
      !     to the center of the faces
      !
      !     inter/extrapolation of v* at the center of the elements
      !     to the center of the faces;
      !
      implicit none
      !
      integer ielfa(4,*),iel1,iel2,iel3,i,j,icyclic,ifatie(*),&
        nef,ipnei(*),indexf,nfacea,nfaceb,ncfd
      !
      real*8 xrlfa(3,*),xl1,xl2,advfa(*),adv(*),vel(nef,0:7),hfa(3,*),&
           c(3,3),volume(*),auv(*),a1,a2,a3
      !
      intent(in) ielfa,xrlfa,adv,icyclic,c,ifatie,vel,nef,volume,&
                        auv,ipnei,nfacea,nfaceb,ncfd
      !
      intent(inout) advfa,hfa
      !
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            !
            !           internal face
            !
            xl2=xrlfa(2,i)
            !
            a1=adv(iel1)
            do indexf=ipnei(iel1)+1,ipnei(iel1+1)
               a1=a1+auv(indexf)
            enddo
            !
            a2=adv(iel2)
            do indexf=ipnei(iel2)+1,ipnei(iel2+1)
               a2=a2+auv(indexf)
            enddo
            !
            advfa(i)=xl1*volume(iel1)/a1&
                    +xl2*volume(iel2)/a2
            !
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do j=1,ncfd
                  hfa(j,i)=(xl1*vel(iel1,j)&
                       +xl2*vel(iel2,j))
               enddo
            elseif(ifatie(i).gt.0) then
               do j=1,ncfd
                  hfa(j,i)=(xl1*vel(iel1,j)&
                       +xl2*(c(j,1)*vel(iel2,1)+&
                             c(j,2)*vel(iel2,2)+&
                             c(j,3)*vel(iel2,3)))
               enddo
            else
               do j=1,ncfd
                  hfa(j,i)=(xl1*vel(iel1,j)&
                       +xl2*(c(1,j)*vel(iel2,1)+&
                             c(2,j)*vel(iel2,2)+&
                             c(3,j)*vel(iel2,3)))
               enddo
            endif
         elseif(ielfa(3,i).ne.0) then
            !
            !           external face; linear extrapolation
            !
            a1=adv(iel1)
            do indexf=ipnei(iel1)+1,ipnei(iel1+1)
               a1=a1+auv(indexf)
            enddo
            !
            iel3=abs(ielfa(3,i))
            a3=adv(iel3)
            do indexf=ipnei(iel3)+1,ipnei(iel3+1)
               a3=a3+auv(indexf)
            enddo
            !
            advfa(i)=xl1*volume(iel1)/a1&
                    +xrlfa(3,i)*volume(iel3)/a3
            !
            do j=1,ncfd
               hfa(j,i)=(xl1*vel(iel1,j)+xrlfa(3,i)*vel(iel3,j))
            enddo
         else
            !
            !           external face: constant extrapolation (only one adjacent
            !           element layer)
            !
            a1=adv(iel1)
            do indexf=ipnei(iel1)+1,ipnei(iel1+1)
               a1=a1+auv(indexf)
            enddo
            !
            advfa(i)=volume(iel1)/a1
            !
            do j=1,ncfd
               hfa(j,i)=vel(iel1,j)
            enddo
         endif
      enddo
      !
      return
      end
