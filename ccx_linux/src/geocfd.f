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
      subroutine geocfd(nef,ipkonf,konf,lakonf,co,coel,cofa,nface,
     &  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa,
     &  volume,neifa,xxj,cosb,hmin,ifatie,cs,tieset,icyclic,c,neij,
     &  physcon,isolidsurf,nsolidsurf,dy,xxni,xxnj,xxicn,nflnei,
     &  iturbulent,rf,yy,vel,velo,veloo,xxna,ale,alet,h)
!
!     calculating geometric variables of the cells and their faces
!
      implicit none
!
      character*8 lakonf(*)
      character*81 tieset(3,*)
!
      integer nef,ipkonf(*),konf(*),nface,ielfa(4,*),ipnei(*),neiel(*),
     &  ifaceq(8,6),i,j,k,indexe,kflag,index1,index2,j1,nope,
     &  nodes(4),iel1,iel2,iel3,iface,indexf,neifa(*),nf(5),ifacet(7,4),
     &  ifacew(8,5),ied4(2,6),ied6(2,9),ied8(2,12),ifatie(*),
     &  ics,itie,neighface,ifirst_occurrence,icyclic,neij(*),jface,
     &  isolidsurf(*),nsolidsurf,iopp8(4,6),iopp6(3,5),iopp4(1,4),
     &  node,nelem,nflnei,iturbulent,nx(nsolidsurf),ny(nsolidsurf),
     &  nz(nsolidsurf),kneigh,neigh(1),nfa(nsolidsurf)
!
      real*8 co(3,*),coel(3,*),cofa(3,*),area(*),xxn(3,*),xxi(3,*),
     &  xle(*),xlen(*),xlet(*),xrlfa(3,*),cosa(*),xsj2(3),xi,et,
     &  shp2(7,4),xs2(3,7),xl2(3,8),xl13,volume(*),dxsj2,xl(3,8),
     &  xxj(3,*),cosb(*),hmin,cs(17,*),xn(3),theta,pi,dc,ds,dd,
     &  c(3,3),diff(3),p(3),q(3),a(3),physcon(*),dy(*),xs(3,3),
     &  aa,bb,cc,dist,xxni(3,*),xxnj(3,*),xxicn(3,*),rf(3,*),x13(3),
     &  yy(*),x(nsolidsurf),y(nsolidsurf),z(nsolidsurf),xo(nsolidsurf),
     &  yo(nsolidsurf),zo(nsolidsurf),xp,yp,zp,vel(nef,0:7),h(*),
     &  velo(nef,0:7),veloo(nef,0:7),areal,xxna(3,*),ale(*),alet(*)
!
!     nodes belonging to the cell faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data nf /3,3,4,4,4/
      data ied4 /1,2,2,3,3,1,1,4,2,4,3,4/
      data ied6 /1,2,2,3,3,1,4,5,5,6,6,4,1,4,2,5,3,6/
      data ied8 /1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,1,1,5,2,6,3,7,4,8/
      data iopp4 /4,3,1,2/
      data iopp6 /4,5,6,1,3,2,3,6,0,1,4,0,2,5,0/
      data iopp8 /5,6,7,8,4,3,2,1,3,4,8,7,4,1,5,8,1,2,6,5,2,3,7,6/
!
      ifirst_occurrence=1
      icyclic=0
!
!     coordinates of the center of the cells
!
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         indexe=ipkonf(i)
         if(lakonf(i)(4:4).eq.'8') then
            nope=8
         else if(lakonf(i)(4:4).eq.'6') then
            nope=6
         else
            nope=4
         endif
         do j=1,3
            do k=1,nope
               coel(j,i)=coel(j,i)+co(j,konf(indexe+k))
            enddo
            coel(j,i)=coel(j,i)/nope
         enddo
      enddo
!
      kflag=2
!
!     loop over all faces
!
      do i=1,nface
!
!        check for cyclic symmetry
!
         if(ifatie(i).ne.0) then
            ics=abs(ifatie(i))
            itie=int(cs(17,ics))
            if(tieset(1,itie)(81:81).eq.'P') then
               if(ifirst_occurrence.eq.1) then
                  do k=1,3
                     diff(k)=-cs(5+k,ics)
                  enddo
                  ifirst_occurrence=0
               endif
            elseif(tieset(1,itie)(81:81).eq.'Z') then
               if(ifirst_occurrence.eq.1) then
                  icyclic=1
                  pi=4.d0*datan(1.d0)
!
!                 normal along the cyclic symmetry axis such that
!                 the slave surface is rotated clockwise through the
!                 body into the master surface while looking in 
!                 the direction of xn
!
                  do k=1,3
                     a(k)=cs(5+k,ics)
                     xn(k)=cs(8+k,ics)-a(k)
                  enddo
                  dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
                  do k=1,3
                     xn(k)=xn(k)/dd
                  enddo
!
!                 angle from the master to the slave surface
!
                  theta=-2.d0*pi/cs(1,ics)
!
!                 rotation matrix rotating a vector in the master
!                 surface into a vector in the slave surface
!
                  dc=dcos(theta)
                  ds=dsin(theta)
!
!                 C-matrix from Guido Dhondt, The Finite Element
!                 Method for Three-Dimensional Thermomechanical
!                 Applications p 158
!     
                  c(1,1)=dc+(1.d0-dc)*xn(1)*xn(1)
                  c(1,2)=   (1.d0-dc)*xn(1)*xn(2)-ds*xn(3)
                  c(1,3)=   (1.d0-dc)*xn(1)*xn(3)+ds*xn(2)
                  c(2,1)=   (1.d0-dc)*xn(2)*xn(1)+ds*xn(3)
                  c(2,2)=dc+(1.d0-dc)*xn(2)*xn(2)
                  c(2,3)=   (1.d0-dc)*xn(2)*xn(3)-ds*xn(1)
                  c(3,1)=   (1.d0-dc)*xn(3)*xn(1)-ds*xn(2)
                  c(3,2)=   (1.d0-dc)*xn(3)*xn(2)+ds*xn(1)
                  c(3,3)=dc+(1.d0-dc)*xn(3)*xn(3)
                  ifirst_occurrence=0
               endif
            else
               write(*,*) '*ERROR in geocfd'
               write(*,*) '       kind of cyclic symmetry'
               write(*,*) '       not known'
               stop
            endif
         endif
!
         iel1=ielfa(1,i)
         indexe=ipkonf(iel1)
         j1=ielfa(4,i)
         if(lakonf(iel1)(4:4).eq.'8') then
!
!           hexahedral element
!
!           coordinates of the face centers
!
            do j=1,4
               nodes(j)=konf(indexe+ifaceq(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/4.d0
            enddo
!
            xi=0.d0
            et=0.d0
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            area(i)=4.d0*dxsj2
!
            neighface=6
!
         else if(lakonf(iel1)(4:4).eq.'6') then
!
!           wedge element
!
!           coordinates of the face centers
!
            do j=1,nf(j1)
               nodes(j)=konf(indexe+ifacew(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/nf(j1)
            enddo
!
            xi=0.d0
            et=0.d0
            if(nf(j1).eq.3) then
               call shape3tri(xi,et,xl2,xsj2,xs2,shp2,kflag)
            else
               call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
            endif
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            if(nf(j1).eq.3) then
               area(i)=dxsj2/2.d0
            else
               area(i)=4.d0*dxsj2
            endif
!
            neighface=5
!
         else
!
!           tetrahedral element
!
!           coordinates of the face centers
!
            do j=1,3
               nodes(j)=konf(indexe+ifacet(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/3.d0
            enddo
!
            xi=0.d0
            et=0.d0
            call shape3tri(xi,et,xl2,xsj2,xs2,shp2,kflag)
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            area(i)=dxsj2/2.d0
!
            neighface=4
!
         endif
!
!        normal and xi-vector on face viewed from cell 1
!
         index1=ipnei(iel1)+j1
         do k=1,3
            xxn(k,index1)=xsj2(k)/dxsj2
            xxi(k,index1)=cofa(k,i)-coel(k,iel1)
            rf(k,i)=xxi(k,index1)
         enddo
!     
!     distance from face center to the center of cell 1
!
         xle(index1)=dsqrt(xxi(1,index1)**2+xxi(2,index1)**2+
     &        xxi(3,index1)**2)
         do k=1,3
            xxi(k,index1)=xxi(k,index1)/xle(index1)
         enddo
!     
!     angle between the normal and the xi-vector
!     
         cosa(index1)=xxn(1,index1)*xxi(1,index1)+
     &        xxn(2,index1)*xxi(2,index1)+
     &        xxn(3,index1)*xxi(3,index1)
!     
         iel2=ielfa(2,i)
!     
!     check whether there is an adjacent cell
!     
         if(iel2.ne.0) then
            index2=ipnei(iel2)+neij(index1)
!     
!     normal and xi-vector on face viewed from cell 2
!     
            if(ifatie(i).eq.0) then
!
!              genuine neighbor
!
               do k=1,3
                  xxi(k,index2)=cofa(k,i)-coel(k,iel2)
                  xxn(k,index2)=-xxn(k,index1)
               enddo
!     
               xle(index2)=dsqrt(xxi(1,index2)**2+xxi(2,index2)**2+
     &              xxi(3,index2)**2)
               do k=1,3
                  xxi(k,index2)=xxi(k,index2)/xle(index2)
               enddo
!     
!     angle between the normal and the xi-vector: xxn.xxi
!     
               cosa(index2)=xxn(1,index2)*xxi(1,index2)+
     &              xxn(2,index2)*xxi(2,index2)+
     &              xxn(3,index2)*xxi(3,index2)
!     
!     distance from the face center to the center of the
!     adjacent cell
!     
               xlen(index1)=xle(index2)
               xlen(index2)=xle(index1)
!     
               do k=1,3
                  xxj(k,index2)=coel(k,iel1)-coel(k,iel2)
               enddo
!     
!     distance between the cell center and the center of the
!     adjacent cell
!     
               xlet(index1)=dsqrt(xxj(1,index2)**2+xxj(2,index2)**2
     &              +xxj(3,index2)**2)
               xlet(index2)=xlet(index1)
!     
!     xxj is the unit vector connecting neighboring cell centers
!     
               do k=1,3
                  xxj(k,index2)=xxj(k,index2)/xlet(index2)
                  xxj(k,index1)=-xxj(k,index2)
               enddo
!     
!     xxn.xxj
!     
               cosb(index2)=xxn(1,index1)*xxj(1,index1)+
     &              xxn(2,index1)*xxj(2,index1)+
     &              xxn(3,index1)*xxj(3,index1)
               cosb(index1)=cosb(index2)
!     
c               xrlfa(1,i)=xle(index2)/(xle(index1)+xle(index2))
c               xrlfa(2,i)=xle(index1)/(xle(index1)+xle(index2))
            else
!
!              cyclic symmetry face: some quantities are
!              calculated on the cyclic symmetric face
!
               xlen(index2)=xle(index1)
!
!              rotational cyclic symmetry
!   
!              vector from the axis to the center of the
!              cyclic symmetric cell and orthogonal to the
!              axis
!            
               if(tieset(1,itie)(81:81).eq.'Z') then
                  do k=1,3
                     p(k)=coel(k,iel2)-a(k)
                  enddo
                  dd=p(1)*xn(1)+p(2)*xn(2)+p(3)*xn(3)
                  do k=1,3
                     p(k)=p(k)-dd*xn(k)
                  enddo
!
                  if(ifatie(i).gt.0) then
!
!                    vector rotated in the direction of the slave surface
!                    (iel2 is adjacent to the master surface)
!     
                     do k=1,3
                        q(k)=c(k,1)*p(1)+c(k,2)*p(2)+c(k,3)*p(3)
                     enddo
                  else
!
!                    vector rotated in the direction of the master surface
!                    (iel2 is adjacent to the slave surface)
!     
                     do k=1,3
                        q(k)=c(1,k)*p(1)+c(2,k)*p(2)+c(3,k)*p(3)
                     enddo
                  endif
!
!                 vector connecting the center of the cyclic
!                 symmetry cell with the center of its rotated
!                 ghost cell
!
                  do k=1,3
                     diff(k)=q(k)-p(k)
                  enddo
               endif
!
               do k=1,3
                  xxj(k,index1)=coel(k,iel2)-coel(k,iel1)
     &                         +diff(k)
               enddo
!     
!     distance between the cell center and the center of the
!     adjacent cell
!     
               xlet(index1)=dsqrt(xxj(1,index1)**2+xxj(2,index1)**2
     &              +xxj(3,index1)**2)
!     
!     xxj is the unit vector connecting neighboring cell centers
!     
               do k=1,3
                  xxj(k,index1)=xxj(k,index1)/xlet(index1)
               enddo
!     
!     xxn.xxj
!     
               cosb(index1)=xxn(1,index1)*xxj(1,index1)+
     &              xxn(2,index1)*xxj(2,index1)+
     &              xxn(3,index1)*xxj(3,index1)
            endif
!
!           calculating rf = the shortest vector between the
!                            center of the face and the line connecting
!                            the centers of the adjacent elements. The
!                            corresponding point on this line is called p
!           calculating xrlfa = the ratio of the segments created by
!                               the point p to the total distance between
!                               the centers of the adjacent elements of
!                               the face
!
            dd=rf(1,i)*xxj(1,index1)+
     &         rf(2,i)*xxj(2,index1)+
     &         rf(3,i)*xxj(3,index1)
!
            do k=1,3
               rf(k,i)=rf(k,i)-dd*xxj(k,index1)
            enddo
!
            xrlfa(2,i)=dd/xlet(index1)
            xrlfa(1,i)=1.d0-xrlfa(2,i)
         else
!     
!     xxi and xxj coincide
!     
            do k=1,3
               xxj(k,index1)=xxi(k,index1)
            enddo
            cosb(index1)=cosa(index1)
!     
!     external face: determining the cell next to the
!     adjacent cell
!     
            iel3=ielfa(3,i)
            if(iel3.eq.0) cycle
c            xl13=dsqrt((coel(1,iel1)-coel(1,iel3))**2+
c     &           (coel(2,iel1)-coel(2,iel3))**2+
c     &           (coel(3,iel1)-coel(3,iel3))**2)
c            xrlfa(1,i)=(xl13+xle(index1))/xl13
c            xrlfa(3,i)=1.d0-xrlfa(1,i)
!
!           unit vector pointing from the center in element iel3
!           to the center in element iel1
!
            do k=1,3
               x13(k)=coel(k,iel1)-coel(k,iel3)
            enddo
!
            xl13=dsqrt(x13(1)*x13(1)+x13(2)*x13(2)+x13(3)*x13(3))
!
            do k=1,3
               x13(k)=x13(k)/xl13
            enddo
!
            dd=rf(1,i)*x13(1)+rf(2,i)*x13(2)+rf(3,i)*x13(3)
!
            do k=1,3
               rf(k,i)=rf(k,i)-dd*x13(k)
            enddo
!
            xrlfa(3,i)=-dd/xl13
            xrlfa(1,i)=1.d0-xrlfa(3,i)
!
         endif
      enddo
!
!     for cyclic symmetric faces xrlfa has not been filled yet
!
c      if(ifirst_occurrence.eq.0) then
c         do i=1,nface
c            if(ifatie(i).ne.0) then
c               index1=ipnei(ielfa(1,i))+ielfa(4,i)
c               xrlfa(1,i)=xlen(index1)/(xle(index1)+xlen(index1))
c               xrlfa(2,i)=1.d0-xrlfa(1,i)
c            endif
c         enddo
c      endif
!
!     calculation of the volume of the elements
!
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         indexf=ipnei(i)
         volume(i)=0.d0
         do j=1,ipnei(i+1)-ipnei(i)
            iface=neifa(indexf+j)
            volume(i)=volume(i)+
     &            area(iface)*cofa(1,iface)*xxn(1,indexf+j)
         enddo
      enddo
!     
!     calculation of the minimum length within the cells
!
c      hmin=1.d30
c      do i=1,nef
c         indexe=ipkonf(i)
c         read(lakonf(i)(4:4),'(i1)') nope
c         do j=1,nope
c            do k=1,3
c               xl(k,j)=co(k,konf(indexe+j))
c            enddo
c         enddo
c         if(nope.eq.4) then
c            do j=1,6
c               hmin=min(hmin,(xl(1,ied4(1,j))-xl(1,ied4(2,j)))**2+
c     &                       (xl(2,ied4(1,j))-xl(2,ied4(2,j)))**2+
c     &                       (xl(3,ied4(1,j))-xl(3,ied4(2,j)))**2)
c            enddo
c         elseif(nope.eq.6) then
c            do j=1,9
c               hmin=min(hmin,(xl(1,ied6(1,j))-xl(1,ied6(2,j)))**2+
c     &                       (xl(2,ied6(1,j))-xl(2,ied6(2,j)))**2+
c     &                       (xl(3,ied6(1,j))-xl(3,ied6(2,j)))**2)
c            enddo
c         else
c            do j=1,12
c               hmin=min(hmin,(xl(1,ied8(1,j))-xl(1,ied8(2,j)))**2+
c     &                       (xl(2,ied8(1,j))-xl(2,ied8(2,j)))**2+
c     &                       (xl(3,ied8(1,j))-xl(3,ied8(2,j)))**2)
c            enddo
c         endif
c      enddo
c      hmin=dsqrt(hmin)
c      write(*,*) 'hmin first ',hmin
!     
!     calculation of the minimum height within the cells
!
      hmin=1.d30
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         if(lakonf(i)(4:4).eq.'8') then
            nope=8
         else if(lakonf(i)(4:4).eq.'6') then
            nope=6
         else
            nope=4
         endif
         indexf=ipnei(i)
c         write(*,*) 'geocfd ',i,volume(i)
         do j=1,ipnei(i+1)-ipnei(i)
            indexf=indexf+1
            iface=neifa(indexf)
            if(nope.eq.4) then
               h(indexf)=volume(i)*3.d0/area(iface)
            elseif(nope.eq.6) then
               if(j.le.2) then
                  h(indexf)=volume(i)/area(iface)
               else
                  h(indexf)=volume(i)*2.d0/area(iface)
               endif
            else
               h(indexf)=volume(i)/area(iface)
c               write(*,*) 'geocfd ',indexf,area(iface),h(indexf)
            endif
            hmin=min(h(indexf),hmin)
         enddo
      enddo
c      write(*,*) 'hmin second ',hmin
!
!     calculate the distance to the nearest node for solid surface
!     faces
!
      if(iturbulent.gt.0) then
!
         if(dabs(physcon(5)).le.0.d0) then
            write(*,*) '*ERROR in geocfd: velocity at infinity'
            write(*,*) '       is nonpositive;'
            write(*,*) '       wrong *VALUES AT INFINITY'  
            call exit(201)
         endif
!
         if(dabs(physcon(7)).le.0.d0) then
            write(*,*) '*ERROR in geocfd: density at infinity'
            write(*,*) '       is nonpositive;'
            write(*,*) '       wrong *VALUES AT INFINITY'  
            call exit(201)
         endif
!
         if(dabs(physcon(8)).le.0.d0) then
            write(*,*) '*ERROR in geocfd: length of the '
            write(*,*) '       computational domain is nonpositive;'
            write(*,*) '       wrong *VALUES AT INFINITY'  
            call exit(201)
         endif
!
         do i=1,nsolidsurf
            iface=isolidsurf(i)
            nelem=int(iface/10)
            indexe=ipkonf(nelem)
            jface=iface-nelem*10
!
!           xl contains the coordinates of the nodes belonging
!           to the face
!
            if(lakonf(nelem)(4:4).eq.'8') then
               do j=1,4
                  node=konf(indexe+ifaceq(j,jface))
                  do k=1,3
                     xl(k,j)=co(k,node)
                  enddo
               enddo
            elseif(lakonf(nelem)(4:4).eq.'6') then
               if(jface.le.2) then
                  do j=1,3
                     node=konf(indexe+ifacew(j,jface))
                     do k=1,3
                        xl(k,j)=co(k,node)
                     enddo
                  enddo
               else
                  do j=1,4
                     node=konf(indexe+ifacew(j,jface))
                     do k=1,3
                        xl(k,j)=co(k,node)
                     enddo
                  enddo
               endif
            else
               do j=1,3
                  node=konf(indexe+ifacet(j,jface))
                  do k=1,3
                     xl(k,j)=co(k,node)
                  enddo
               enddo
            endif
!
!           determine the face number in field ielfa
!
            nfa(i)=neifa(ipnei(nelem)+jface)
!
!           determine the plane through the face (exact for
!           3-node face, approximate for a 4-node face)
!
            if((lakonf(nelem)(4:4).eq.'8').or.
     &         ((lakonf(nelem)(4:4).eq.'6').and.(jface.gt.2))) then
!
!              computation of the local derivative of the global coordinates
!             (xs)
!
               do j=1,3
                  xs(j,1)=-xl(j,1)+xl(j,2)+xl(j,3)-xl(j,4)
                  xs(j,2)=-xl(j,1)-xl(j,2)+xl(j,3)+xl(j,4)
               enddo
!
!              computation of the jacobian vector for xi,et=0
!
               aa=xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2)
               bb=xs(1,2)*xs(3,1)-xs(3,2)*xs(1,1)
               cc=xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
               dd=dsqrt(aa*aa+bb*bb+cc*cc)
               aa=aa/dd
               bb=bb/dd
               cc=cc/dd
               dd=-(aa*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))
     &             +bb*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))
     &             +cc*(xl(3,1)+xl(3,2)+xl(3,3)+xl(3,4)))/4.d0
            else
!
!              computation of the local derivative of the global coordinates
!             (xs)
!
               do j=1,3
                  xs(j,1)=-xl(j,1)+xl(j,2)
                  xs(j,2)=-xl(j,1)+xl(j,3)
               enddo
!     
!              computation of the jacobian vector (unique for triangle)
!     
               aa=xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2)
               bb=xs(1,2)*xs(3,1)-xs(3,2)*xs(1,1)
               cc=xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
               dd=dsqrt(aa*aa+bb*bb+cc*cc)
               aa=aa/dd
               bb=bb/dd
               cc=cc/dd
               dd=-(aa*xl(1,1)+bb*xl(2,1)+cc*xl(3,1))
            endif
!
!           determine the shortest distance within the element
!           from the solid surface face
!
            dist=1.d30
            if(lakonf(nelem)(4:4).eq.'8') then
               do j=1,4
                  node=konf(indexe+iopp8(j,jface))
                  dist=min(dist,-(aa*co(1,node)+bb*co(2,node)
     &                           +cc*co(3,node)+dd))
               enddo
            elseif(lakonf(nelem)(4:4).eq.'6') then
               if(jface.le.2) then
                  do j=1,3
                     node=konf(indexe+iopp6(j,jface))
                     dist=min(dist,-(aa*co(1,node)+bb*co(2,node)
     &                    +cc*co(3,node)+dd))
                  enddo
               else
                  do j=1,2
                     node=konf(indexe+iopp6(j,jface))
                     dist=min(dist,-(aa*co(1,node)+bb*co(2,node)
     &                    +cc*co(3,node)+dd))
                  enddo
               endif
            else
               node=konf(indexe+iopp4(1,jface))
               dist=min(dist,-(aa*co(1,node)+bb*co(2,node)
     &                        +cc*co(3,node)+dd))
            endif
!
!           60.d0/(0.075d0*delta(y)**2)
!
            dy(i)=800.d0/(dist*dist)
         enddo
      endif
!
!     calculate for each element center the shortest distance to a solid
!     surface (only for the BSL and SST turbulence model)
!
      if(iturbulent.gt.2) then
!
         do i=1,nsolidsurf
            iface=nfa(i)
            x(i)=cofa(1,iface)
            y(i)=cofa(2,iface)
            z(i)=cofa(3,iface)
            xo(i)=x(i)
            yo(i)=y(i)
            zo(i)=z(i)
            nx(i)=i
            ny(i)=i
            nz(i)=i
         enddo
!
         kflag=2
         call dsort(x,nx,nsolidsurf,kflag)
         call dsort(y,ny,nsolidsurf,kflag)
         call dsort(z,nz,nsolidsurf,kflag)
!
         kneigh=1
         do i=1,nef
            xp=coel(1,i)
            yp=coel(2,i)
            zp=coel(3,i)
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,nsolidsurf,
     &            neigh,kneigh)
            yy(i)=dsqrt((xp-xo(neigh(1)))**2+
     &                  (yp-yo(neigh(1)))**2+
     &                  (zp-zo(neigh(1)))**2)
         enddo
      endif
!
!     auxiliary fields
!
      do i=1,nflnei
         areal=area(neifa(i))
         do k=1,3
            xxni(k,i)=xxn(k,i)-xxi(k,i)
            xxnj(k,i)=(xxn(k,i)-xxj(k,i))*areal
            xxicn(k,i)=xxn(k,i)-xxj(k,i)/cosb(i)
            xxna(k,i)=xxn(k,i)*areal
         enddo
         ale(i)=areal/xle(i)
         alet(i)=areal/xlet(i)
         cosa(i)=ale(i)/cosa(i)
      enddo
c
c     initial conditions
c
c         do i=1,nef
c            vel(i,0)=1.d0+coel(2,i)*(1.d0-coel(2,i))/2.d0
c            vel(i,1)=coel(2,i)
c            vel(i,2)=0.d0
c            vel(i,3)=0.d0
c            vel(i,4)=1.d0
cc            do j=0,4
cc               velo(i,j)=vel(i,j)
cc               veloo(i,j)=vel(i,j)
cc            enddo
c         enddo
c      write(*,*) 'geocfd neifa,neiel'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),neifa(i),neiel(i)
c      enddo
c      write(*,*) 'geocfd xle,xlen,xlet'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),xle(i),xlen(i),xlet(i)
c      enddo
c      write(*,*) 'geocfd xxn'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),(xxn(j,i),j=1,3)
c      enddo
c      write(*,*) 'geocfd xxi'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),(xxi(j,i),j=1,3)
c      enddo
c      write(*,*) 'geocfd xxj'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),(xxj(j,i),j=1,3)
c      enddo
c      write(*,*) 'geocfd cosa,cosb'
c      do i=1,6*nef
c         write(*,*) (i-1)/6+1,i-6*((i-1)/6),cosa(i),cosb(i)
c      enddo
c      do i=1,nef
c         write(*,*) 'coef ',i,(coel(j,i),j=1,3)
c      enddo
c      do i=1,nface
c         write(*,*) 'cofa ',i,(cofa(j,i),j=1,3)
c      enddo
c      do i=1,nface
c         write(*,*) 'rf ',i,(rf(j,i),j=1,3)
c      enddo
!
      return
      end
