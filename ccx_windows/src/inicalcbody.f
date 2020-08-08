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
      subroutine inicalcbody(nef,body,ipobody,ibody,xbody,coel,vel,
     &  lakon,nactdohinv,icent)
!
!     calculation of the actual body force in each element. The body
!     force is the sum of gravity and centrifugal/Coriolis forces
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,j,nef,index,ipobody(2,*),ibody(3,*),nactdohinv(*),icent
!
      real*8 om,body(0:3,*),p1(3),p2(3),xbody(7,*),omcor,q(3),coel(3,*),
     &  vel(nef,0:7),const,corio(3)
!
      do i=1,nef
         om=0.d0
         do j=1,3
            body(j,i)=0.d0
         enddo
!
         index=nactdohinv(i)
!
         do
            j=ipobody(1,index)
            if(j.eq.0) exit
            if(ibody(1,j).eq.1) then
               om=xbody(1,j)
               p1(1)=xbody(2,j)
               p1(2)=xbody(3,j)
               p1(3)=xbody(4,j)
               p2(1)=xbody(5,j)
               p2(2)=xbody(6,j)
               p2(3)=xbody(7,j)
!
               icent=1
!     
!     assigning gravity forces
!     
            elseif(ibody(1,j).eq.2) then
               body(1,i)=body(1,i)+xbody(1,j)*xbody(2,j)
               body(2,i)=body(2,i)+xbody(1,j)*xbody(3,j)
               body(3,i)=body(3,i)+xbody(1,j)*xbody(4,j)
            endif
            index=ipobody(2,index)
            if(index.eq.0) exit
         enddo
!
!        adding the centrifugal/Coriolis force (if any) to the
!        gravity loading
!
         if(om.gt.0.d0) then
            omcor=2.d0*dsqrt(om)
            do j=1,3
               q(j)=coel(j,i)-p1(j)
            enddo
            const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!     
!     Coriolis forces
!     
            corio(1)=vel(i,2)*p2(3)-vel(i,3)*p2(2)
            corio(2)=vel(i,3)*p2(1)-vel(i,1)*p2(3)
            corio(3)=vel(i,1)*p2(2)-vel(i,2)*p2(1)
!     
!     inclusion of the centrifugal force into the body force
!     
            do j=1,3
               body(j,i)=body(j,i)+(q(j)-const*p2(j))*om+
     &              corio(j)*omcor
            enddo
         endif
!
      enddo
!     
      return
      end
