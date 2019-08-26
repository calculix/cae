!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine normalsonsurface_se(ipkon,kon,lakon,extnor,co,nk,&
            ipoface,nodface,nactdof,mi,nodedesiinv,iregion)
      !
      !     calculating the normal direction onto the external surface;
      !     the design variables are moved in this direction.
      !
      !     if the design variables constitute a region, i.e. they
      !     are constitute a set of faces, the normals at the boundary
      !     of this set is determined based on the faces belonging to
      !     this set only (so faces external to this set do not
      !     contribute to the normal at the boundary)
      !
      !     if it is not a region, the normal in a node is the mean
      !     of the normal of all external faces to which this node
      !     belongs, no matter how many other design variables
      !     belong to these faces (if any)
      !
      implicit none
      !
      character*8 lakon(*)
      !
      integer j,nelemm,jfacem,indexe,ipkon(*),kon(*),nopem,node,&
        ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),&
        konl(26),ipoface(*),nodface(5,*),mi(*),nodedesiinv(*),&
        nactdof(0:mi(2),*),nopesurf(9),nnodes,iregion,nope,&
        nopedesi,l,m,iflag,k,nk
      !
      real*8 extnor(3,*),xsj2(3),shp2(7,9),xs2(3,2),xi,et,dd,&
        xquad(2,9),xtri(2,7),xl2m(3,9),co(3,*)
      !
      intent(in) ipkon,kon,lakon,co,nk,&
            ipoface,nodface,nactdof,mi,nodedesiinv,iregion
      !
      intent(inout) extnor
      !
      !     nodes per face for hex elements
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
                  5,6,7,8,13,14,15,16,&
                  1,2,6,5,9,18,13,17,&
                  2,3,7,6,10,19,14,18,&
                  3,4,8,7,11,20,15,19,&
                  4,1,5,8,12,17,16,20/
      !
      !     nodes per face for tet elements
      !
      data ifacet /1,3,2,7,6,5,&
                   1,2,4,5,9,8,&
                   2,3,4,6,10,9,&
                   1,4,3,8,10,7/
      !
      !     nodes per face for linear wedge elements
      !
      data ifacew1 /1,3,2,0,&
                   4,5,6,0,&
                   1,2,5,4,&
                   2,3,6,5,&
                   3,1,4,6/
      !
      !     nodes per face for quadratic wedge elements
      !
      data ifacew2 /1,3,2,9,8,7,0,0,&
                   4,5,6,10,11,12,0,0,&
                   1,2,5,4,7,14,10,13,&
                   2,3,6,5,8,15,11,14,&
                   3,1,4,6,9,13,12,15/
      !
      !     new added data for the local coodinates for nodes
      !
      data xquad /-1.d0,-1.d0,&
                   1.d0,-1.d0,&
                   1.d0,1.d0,&
                  -1.d0,1.d0,&
                   0.d0,-1.d0,&
                   1.d0,0.d0,&
                   0.d0,1.d0,&
                  -1.d0,0.d0,&
                   0.d0,0.d0/
      !
      data xtri /0.d0,0.d0,&
                 1.d0,0.d0,&
                 0.d0,1.d0,&
                 .5d0,0.d0,&
                 .5d0,.5d0,&
                 0.d0,.5d0,&
                 0.333333333333333d0,0.333333333333333d0/
      !
      data iflag /2/
      !
      !     calculation of the normal to the external surface;
      !     each external face contributes its normal
      !     to each node belonging to the face providing that
      !     1) the node is not a design variable OR
      !     2) the node is a design variable and belongs to a design face
      !     A design face is an external face for which more than nopedesi
      !     nodes are design variables
      !
      !     The appropriate normal component is set to zero for fixed dofs
      !
      do j=1,nk
         !
         if(ipoface(j).eq.0) cycle
         indexe=ipoface(j)
         !
         do
            !
            nelemm=nodface(3,indexe)
            jfacem=nodface(4,indexe)
            !
            !           nopem: # of nodes in the surface
            !           nope: # of nodes in the element
            !
            if(lakon(nelemm)(4:4).eq.'8') then
               nopem=4
               nope=8
               nopedesi=3
            elseif(lakon(nelemm)(4:5).eq.'20') then
               nopem=8
               nope=20
               nopedesi=5
            elseif(lakon(nelemm)(4:5).eq.'10') then
               nopem=6
               nope=10
               !                nopedesi=3
               nopedesi=4
            elseif(lakon(nelemm)(4:4).eq.'4') then
               nopem=3
               nope=4
               nopedesi=3
            !
            !     treatment of wedge faces
            !
            elseif(lakon(nelemm)(4:4).eq.'6') then
               nope=6
               if(jfacem.le.2) then
                  nopem=3
               else
                  nopem=4
               endif
               nopedesi=3
            elseif(lakon(nelemm)(4:5).eq.'15') then
               nope=15
               if(jfacem.le.2) then
                  nopem=6
                  nopedesi=4
               else
                  nopem=8
                  nopedesi=5
               endif
            !                nopedesi=3
            endif
            if(iregion.eq.0) nopedesi=0
            !
            !     actual position of the nodes belonging to the
            !     master surface
            !
            do k=1,nope
               konl(k)=kon(ipkon(nelemm)+k)
            enddo
            !
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopem
                  nopesurf(m)=konl(ifaceq(m,jfacem))
                  do k=1,3
                     xl2m(k,m)=co(k,nopesurf(m))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4))&
                    then
               do m=1,nopem
                  nopesurf(m)=konl(ifacet(m,jfacem))
                  do k=1,3
                     xl2m(k,m)=co(k,nopesurf(m))
                  enddo
               enddo
            elseif(nope.eq.15) then
               do m=1,nopem
                  nopesurf(m)=konl(ifacew2(m,jfacem))
                  do k=1,3
                     xl2m(k,m)=co(k,nopesurf(m))
                  enddo
               enddo
            else
               do m=1,nopem
                  nopesurf(m)=konl(ifacew1(m,jfacem))
                  do k=1,3
                     xl2m(k,m)=co(k,nopesurf(m))
                  enddo
               enddo
            endif
            !
            !     sum up how many designvariables are on that surface
            !
            nnodes=0
            do m=1,nopem
               if(nodedesiinv(nopesurf(m)).eq.1) then
                  nnodes=nnodes+1
               endif
            enddo         
            !
            !     calculate the normal vector in the nodes belonging to the master surface
            !
            if(nopem.eq.8) then
               do m=1,nopem
                  xi=xquad(1,m)
                  et=xquad(2,m)
                  call shape8q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                  dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)&
                       + xsj2(3)*xsj2(3))
                  xsj2(1)=xsj2(1)/dd
                  xsj2(2)=xsj2(2)/dd
                  xsj2(3)=xsj2(3)/dd
                  !
                  if(nope.eq.20) then
                     node=konl(ifaceq(m,jfacem))
                  elseif(nope.eq.15) then
                     node=konl(ifacew2(m,jfacem))
                  endif
                  if((nodedesiinv(node).eq.0).or.&
                     ((nodedesiinv(node).eq.1).and.&
                      (nnodes.ge.nopedesi))) then
                     extnor(1,node)=extnor(1,node)&
                          +xsj2(1)
                     extnor(2,node)=extnor(2,node)&
                          +xsj2(2)
                     extnor(3,node)=extnor(3,node)&
                          +xsj2(3)
                  endif
               enddo
            elseif(nopem.eq.4) then
               do m=1,nopem
                  xi=xquad(1,m)
                  et=xquad(2,m)
                  call shape4q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                  dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)&
                       + xsj2(3)*xsj2(3))
                  xsj2(1)=xsj2(1)/dd
                  xsj2(2)=xsj2(2)/dd
                  xsj2(3)=xsj2(3)/dd
                  !
                  if(nope.eq.8) then
                     node=konl(ifaceq(m,jfacem))
                  elseif(nope.eq.6) then
                     node=konl(ifacew1(m,jfacem))
                  endif
                  if((nodedesiinv(node).eq.0).or.&
                     ((nodedesiinv(node).eq.1).and.&
                      (nnodes.ge.nopedesi))) then
                     extnor(1,node)=extnor(1,node)&
                          +xsj2(1)
                     extnor(2,node)=extnor(2,node)&
                          +xsj2(2)
                     extnor(3,node)=extnor(3,node)&
                          +xsj2(3)
                  endif
               enddo
            elseif(nopem.eq.6) then
               do m=1,nopem
                  xi=xtri(1,m)
                  et=xtri(2,m)
                  call shape6tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                  dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)&
                       + xsj2(3)*xsj2(3))
                  xsj2(1)=xsj2(1)/dd
                  xsj2(2)=xsj2(2)/dd
                  xsj2(3)=xsj2(3)/dd
                  !
                  if(nope.eq.10) then
                     node=konl(ifacet(m,jfacem))
                  elseif(nope.eq.15) then
                     node=konl(ifacew2(m,jfacem))
                  endif
                  if((nodedesiinv(node).eq.0).or.&
                     ((nodedesiinv(node).eq.1).and.&
                      (nnodes.ge.nopedesi))) then
                     extnor(1,node)=extnor(1,node)&
                          +xsj2(1)
                     extnor(2,node)=extnor(2,node)&
                          +xsj2(2)
                     extnor(3,node)=extnor(3,node)&
                          +xsj2(3)
                  endif
               enddo
            else
               do m=1,nopem
                  xi=xtri(1,m)
                  et=xtri(2,m)
                  call shape3tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                  dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)&
                       + xsj2(3)*xsj2(3))
                  xsj2(1)=xsj2(1)/dd
                  xsj2(2)=xsj2(2)/dd
                  xsj2(3)=xsj2(3)/dd
                  !
                  if(nope.eq.6) then
                     node=konl(ifacew1(m,jfacem))
                  elseif(nope.eq.4) then
                     node=konl(ifacet(m,jfacem))
                  endif
                  if((nodedesiinv(node).eq.0).or.&
                     ((nodedesiinv(node).eq.1).and.&
                      (nnodes.ge.nopedesi))) then
                     extnor(1,node)=extnor(1,node)&
                          +xsj2(1)
                     extnor(2,node)=extnor(2,node)&
                          +xsj2(2)
                     extnor(3,node)=extnor(3,node)&
                          +xsj2(3)
                  endif
               enddo
            endif
            !
            indexe=nodface(5,indexe)
            if(indexe.eq.0) exit
         !
         enddo
      enddo
      !
      !     not considering the directions with fixed displacements
      !
      do l=1,nk
         do m=1,3
            if(nactdof(m,l).le.0) then
               extnor(m,l)=0.d0
            endif
         enddo
      enddo 
      !
      !     normalizing the normals
      !
      do l=1,nk
         dd=dsqrt(extnor(1,l)**2+extnor(2,l)**2+&
              extnor(3,l)**2)
         if(dd.gt.0.d0) then
            do m=1,3
               extnor(m,l)=extnor(m,l)/dd
            enddo
         endif
      enddo    
      !
      return
      end
