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
      subroutine gencycsymelemcfd(cs,islav,nslav,imast,nmast,inomat,
     &  nk,co,ne,ipkon,lakon,kon,nkon,ielmat,mi,vold,ielslav,ielmast,
     &  inoslav,inomast,iponoel,inoel)
!
!     creates an additional layer of elements on the slave and the
!     master surface of a cyclic symmetric structure
!
      implicit none
!
      character*8 lakon(*)
!
      integer imast(*),i,j,k,nmast,islav(*),nslav,
     &  nk,ne,ipkon(*),kon(*),nkon,
     &  mi(*),index,nodes,nodem,nope,
     &  ielmat(mi(3),*),indexe,inoslav(*),inomast(*),ielslav(*),
     &  ielmast(*),iponoel(*),inoel(3,*),nodeslav,nodemast,ielemslav,
     &  ielemmast,inomat(*)
!
      real*8 cs(17,*),phibasis,phi,xa(3),xn(3),dd,cphi,sphi,dphi,
     &  c(3,3),xp(3),xq(3),al,co(3,*),vold(0:mi(2),*)
!
!     mcs=1 is assumed
!     phi is the angle of the cyclic basis sector
!
      phibasis=8.d0*datan(1.d0)/cs(1,1)
!
!     xn is an normed vector along the axis
!     xa is a point on the axis
!
      do i=1,3
         xa(i)=cs(5+i,1)
         xn(i)=cs(8+i,1)-xa(i)
      enddo
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!
!     new nodes on the master side
!
!     phi is the thickness of the layer in radians
!
      phi=phibasis
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      c(1,1)=cphi+dphi*xn(1)*xn(1)
      c(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      c(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      c(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      c(2,2)=cphi+dphi*xn(2)*xn(2)
      c(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      c(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      c(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      c(3,3)=cphi+dphi*xn(3)*xn(3)
!
!     new nodes and elements on the master side
!
      do i=1,nmast
!
!        corresponding slave node
!
         nodeslav=inoslav(imast(i))
!
!        loop over all elements to which nodeslav belongs
!
         index=iponoel(nodeslav)
         if(index.eq.0) cycle
!
         do
            ielemslav=inoel(1,index)
            if(ielmast(ielemslav).eq.0) then
               read(lakon(ielemslav)(4:4),'(i1)') nope
               indexe=ipkon(ielemslav)
               do k=1,nope
                  nodes=kon(indexe+k)
                  if(inomast(nodes).eq.0) then
!
!                    generate new node on the master side
!
                     nk=nk+1
                     inomast(nodes)=nk
!
!                    calculating the coordinates of the node
!
                     do j=1,3
                        xp(j)=co(j,nodes)
                     enddo
!
                     al=(xp(1)-xa(1))*xn(1)+
     &                    (xp(2)-xa(2))*xn(2)+
     &                    (xp(3)-xa(3))*xn(3)
!     
!                    xq is the orthogonal projection of xp on the axis
!     
                     do j=1,3
                        xq(j)=xa(j)+al*xn(j)
                     enddo
!     
                     do j=1,3
                        co(j,nk)=xq(j)+c(j,1)*(xp(1)-xq(1))+
     &                       c(j,2)*(xp(2)-xq(2))+
     &                       c(j,3)*(xp(3)-xq(3))
                     enddo
!
!                    initial conditions (global rectangular system)
!                     
                     do j=0,4
                        vold(j,nk)=vold(j,nodes)
                     enddo
!
!                    assign material to node
!
                     inomat(nk)=inomat(nodes)
!
                  endif
               enddo
!
!              create a new element on the master side
!
               ne=ne+1
               ielmast(ielemslav)=ne
               ipkon(ne)=nkon
               lakon(ne)=lakon(ielemslav)
               ielmat(1,ne)=ielmat(1,ielemslav)
               do j=1,nope
                  kon(nkon+j)=inomast(kon(indexe+j))
               enddo
               nkon=nkon+nope
            endif
!
            index=inoel(3,index)
            if(index.eq.0) exit
         enddo
      enddo
!
!     new nodes on the slave side
!
      phi=-phibasis
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      c(1,1)=cphi+dphi*xn(1)*xn(1)
      c(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      c(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      c(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      c(2,2)=cphi+dphi*xn(2)*xn(2)
      c(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      c(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      c(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      c(3,3)=cphi+dphi*xn(3)*xn(3)
!
!     new nodes and elements on the slave side
!
      do i=1,nslav
!
!        corresponding slave node
!
         nodemast=inomast(islav(i))
!
!        loop over all elements to which nodeslav belongs
!
         index=iponoel(nodemast)
         if(index.eq.0) cycle
!
         do
            ielemmast=inoel(1,index)
            if(ielslav(ielemmast).eq.0) then
               read(lakon(ielemmast)(4:4),'(i1)') nope
               indexe=ipkon(ielemmast)
               do k=1,nope
                  nodem=kon(indexe+k)
                  if(inoslav(nodem).eq.0) then
!
!                    generate new node on the master side
!
                     nk=nk+1
                     inoslav(nodem)=nk
!
!                    calculating the coordinates of the node
!
                     do j=1,3
                        xp(j)=co(j,nodem)
                     enddo
!
                     al=(xp(1)-xa(1))*xn(1)+
     &                    (xp(2)-xa(2))*xn(2)+
     &                    (xp(3)-xa(3))*xn(3)
!     
!                    xq is the orthogonal projection of xp on the axis
!     
                     do j=1,3
                        xq(j)=xa(j)+al*xn(j)
                     enddo
!     
                     do j=1,3
                        co(j,nk)=xq(j)+c(j,1)*(xp(1)-xq(1))+
     &                       c(j,2)*(xp(2)-xq(2))+
     &                       c(j,3)*(xp(3)-xq(3))
                     enddo
!
!                    initial conditions (global rectangular system)
!                     
                     do j=0,4
                        vold(j,nk)=vold(j,nodem)
                     enddo
!
!                    assign material to node
!
                     inomat(nk)=inomat(nodem)
!
                  endif
               enddo
!
!              create a new element on the slave side
!
               ne=ne+1
               ielslav(ielemmast)=ne
               ipkon(ne)=nkon
               lakon(ne)=lakon(ielemmast)
               ielmat(1,ne)=ielmat(1,ielemmast)
               do j=1,nope
                  kon(nkon+j)=inoslav(kon(indexe+j))
               enddo
               nkon=nkon+nope
            endif
!
            index=inoel(3,index)
            if(index.eq.0) exit
         enddo
      enddo
!
      return
      end

