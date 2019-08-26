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
      subroutine projectnodes(nktet_,ipoed,iedgmid,iexternedg,&
        iedgext,cotet,nktet,iedg,iquad,integerglob,doubleglob,&
        iexternfa,ifacext,itreated,ialsete,nexternel,ilist,&
        isharp,ipofa,ifac,iedgextfa,ifacexted,nexternedg)
      !
      !     1. projects nodes lying on external edges of the unrefined mesh
      !     on those edges
      !
      !     2. generates middles nodes if necessary
      !
      implicit none
      !
      integer nktet_,ipoed(*),index,i,j,k,iedgmid(*),iexternedg(*),&
        imasted,iedgext(3,*),nterms,node1,node2,iedg(3,*),iquad,&
        nktet,integerglob(*),nktri,netet,ne,nkon,nfaces,nfield,&
        nselect,iselect(1),imastset,imastfa,m,istartset(1),&
        iendset(1),ialsete(*),ifacext(6,*),nelem,iexternfa(*),&
        iface,konl(20),indexf,itreated(*),ii,node,loopa,&
        nexternel,kontypindex,ipkonindex,konindex,isharp(*),&
        ilist(*),nlist,ifac(4,*),iedgextfa(2,*),ifacexted(3,*),&
        id,indexe,isol,ipofa(*),nexternedg,iflag
      !
      real*8 pneigh(3,9),cotet(3,*),pnode(3),ratio(9),dist,xi,et,&
        doubleglob(*),value,distmin,pnodeproj(3),xl(3,3),xsj(3),&
        xs(3,7),shp(7,3),xn1(3),xn2(3),dd
      !
      intent(in) nktet_,ipoed,iexternedg,nexternedg,&
        iedgext,iedg,iquad,integerglob,doubleglob,&
        iexternfa,ifacext,ialsete,nexternel,&
        ipofa,ifac,iedgextfa,ifacexted
      !
      intent(inout) iedgmid,cotet,nktet,isharp,itreated,ilist
      !
      !       loopa=2
      ! !
      ! !     initializing the fields for the projection on the faces
      ! !     (only if the output is a quadratic mesh, i.e. iquad=1)
      ! !
      ! c      if(iquad.eq.1) then
      !          nktri=integerglob(1)
      !          netet=integerglob(2)
      !          ne=integerglob(3)
      !          nkon=integerglob(4)
      !          nfaces=integerglob(5)
      !          nfield=1
      !          nselect=1
      !          iselect(1)=1
      !          istartset(1)=1
      !          iendset(1)=nexternel
      !          imastset=1
      ! !
      ! !        changing the element type to C3D10 (kontyp 6) if applicable
      ! !
      ! !        the next three indices point to a position in integerglob
      ! !        immediately before the storage of the appropriate field
      ! !
      !          kontypindex=7*netet+5
      !          ipkonindex=ne+7*netet+5
      !          konindex=2*ne+7*netet+5
      ! !
      !          do i=1,netet
      ! c            write(*,*) 'projectnodes ',i,
      ! c     &           integerglob(konindex+integerglob(ipkonindex+i)+5)
      !             if(integerglob(konindex+integerglob(ipkonindex+i)+5).ne.0)
      !      &             integerglob(kontypindex+i)=6
      !          enddo
      ! c      endif
      !
      !     loop over all edges: projection on edges from the unrefined mesh
      !
      loop1: do i=1,nktet_
         index=ipoed(i)
         !
         do
            if(index.eq.0) cycle loop1
            ! !
            ! !           end nodes belonging to the edge
            ! !
            !             node1=iedg(1,index)
            !             node2=iedg(2,index)
            !
            if(iexternedg(index).gt.0) then
               !
               !              edge is a subset of an external edge
               !              of the unrefined mesh
               !
               !              end nodes belonging to the edge
               !
               node1=iedg(1,index)
               node2=iedg(2,index)
               !
               !              recovering the nodes of the external edge of the
               !              unrefined mesh
               !
               imasted=iexternedg(index)
               !
               !              an external edge from the unrefined mesh used
               !              as parent edge is labeled as "sharp"
               !
               isharp(imasted)=1
               !
               if(iedgext(2,imasted).ne.0) then
                  nterms=3
                  do j=1,3
                     do k=1,3
                        pneigh(k,j)=cotet(k,iedgext(j,imasted))
                     enddo
                  enddo
               else
                  nterms=2
               endif
               !
               !              projection is only needed if the master edge is
               !              quadratic
               !
               if(nterms.eq.3) then
                  !
                  !                 attach the first end node of the edge
                  !
                  if(itreated(node1).eq.0) then
                     do k=1,3
                        pnode(k)=cotet(k,node1)
                     enddo
                     call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                     do k=1,3
                        cotet(k,node1)=pnode(k)
                     enddo
                     itreated(node1)=1
                  !                      write(*,*) 'projectnodes',node1,itreated(node1)
                  endif
                  !
                  !                 attach the other end node of the edge
                  !
                  if(itreated(node2).eq.0) then
                     do k=1,3
                        pnode(k)=cotet(k,node2)
                     enddo
                     call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                     do k=1,3
                        cotet(k,node2)=pnode(k)
                     enddo
                     itreated(node2)=1
                  !                      write(*,*) 'projectnodes',node2,itreated(node2)
                  endif
               else
                  itreated(node1)=1
                  itreated(node2)=1
               !                      write(*,*) 'projectnodes',node1,itreated(node1)
               !                      write(*,*) 'projectnodes',node2,itreated(node2)
               endif
               !
               !              create a middle node (if necessary) and attach it
               !              (if necessary)
               !
               if(iquad.eq.1) then
                  nktet=nktet+1
                  iedgmid(index)=nktet
                  !
                  if(nterms.eq.3) then
                     do k=1,3
                        pnode(k)=(cotet(k,node1)+cotet(k,node2))/2.d0
                     enddo
                     !
                     !                    projection for quadratic master edges
                     !
                     call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                     !
                     do k=1,3
                        cotet(k,nktet)=pnode(k)
                     enddo
                  else
                     do k=1,3
                        cotet(k,nktet)=&
                            (cotet(k,node1)+cotet(k,node2))/2.d0
                     enddo
                  endif
               endif
            endif
            !
            index=iedg(3,index)
         enddo
      enddo loop1
      !
      !     check whether the sharp external edges are really sharp
      !
      do i=1,nexternedg
         if(isharp(i).eq.0) cycle
         !
         !        first neighboring face
         !
         imastfa=iedgextfa(1,i)
         do j=1,3
            do k=1,3
               xl(k,j)=cotet(k,ifacext(j,imastfa))
            enddo
         enddo
         iflag=2
         call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
         dd=dsqrt(xsj(1)*xsj(1)+xsj(2)*xsj(2)+xsj(3)*xsj(3))
         do j=1,3
            xn1(j)=xsj(j)/dd
         enddo
         !
         !        second neighboring face
         !
         imastfa=iedgextfa(2,i)
         do j=1,3
            do k=1,3
               xl(k,j)=cotet(k,ifacext(j,imastfa))
            enddo
         enddo
         iflag=2
         call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
         dd=dsqrt(xsj(1)*xsj(1)+xsj(2)*xsj(2)+xsj(3)*xsj(3))
         do j=1,3
            xn2(j)=xsj(j)/dd
         enddo
         !
         !        if the normals are nearly parallel, the edge is no sharp edge
         !        "nearly parallel" means that the angle between the vectors
         !        is smaller than 0.0464 degrees.
         !
         if(dabs(xn1(1)*xn2(1)+xn1(2)*xn2(2)+xn1(3)*xn2(3)-1.d0)&
             .lt.1.d-10) isharp(i)=0
      enddo
      !
      !     projection of the external nodes not treated yet onto the
      !     faces of the unrefined mesh
      !
      loop2: do i=1,nktet_
         !
         !        loop over all faces of the refined mesh
         !
         index=ipofa(i)
         do
            if(index.eq.0) cycle loop2
            !
            !           if no external face: loop
            !
            if(iexternfa(index).le.0) then
               index=ifac(4,index)
               cycle
            endif
            !
            !           external face; treat the vertex nodes of the face
            !
            do ii=1,3
               node=ifac(ii,index)
               if(itreated(node).ne.0) cycle
               !                write(*,*) 'projectnodes face projection node ',node
               !
               !              parent face
               !
               imastfa=iexternfa(index)
               ilist(1)=imastfa
               nlist=1
               do k=1,3
                  pnode(k)=cotet(k,node)
               enddo
               !
               !              start the loop looking for the correct face;
               !              starting with the parent face
               !
               do
                  !                   write(*,*) 'imastfa ',imastfa
                  if(ifacext(4,imastfa).ne.0) then
                     nterms=6
                     do j=1,6
                        do k=1,3
                           pneigh(k,j)=cotet(k,ifacext(j,imastfa))
                        enddo
                     enddo
                  else
                     nterms=3
                     do j=1,3
                        do k=1,3
                           pneigh(k,j)=cotet(k,ifacext(j,imastfa))
                        enddo
                     enddo
                  endif
                  !
                  do k=1,3
                     pnodeproj(k)=pnode(k)
                  enddo
                  call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,&
                       xi,et)
                  !                   write(*,*) 'xi et ',xi,et
                  !
                  !                 check whether this face is the correct one;
                  !                 if dabs(xi)=1 or dabs(et)=1 or xi+et=0 this
                  !                 may not be the case
                  !
                  !                 the solution is found (isol=1) unless proved
                  !                 otherwise
                  !
                  isol=1
                  !
                  if(dabs(et).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 1-2 unless sharp
                     !
                     imasted=ifacexted(1,imastfa)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     endif
                  endif
                  !
                  if(dabs(xi+et-1.d0).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 2-3 unless sharp
                     !
                     imasted=ifacexted(2,imastfa)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     endif
                  endif
                  !
                  if(dabs(xi).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 3-1 unless sharp
                     !
                     imasted=ifacexted(3,imastfa)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     endif
                  endif
                  !
                  !                 if solution is found: copy projected coordinates
                  !                 else continue with a neighbor
                  !
                  if(isol.eq.1) then
                     do k=1,3
                        cotet(k,node)=pnodeproj(k)
                     enddo
                     itreated(node)=2
                     !                      write(*,*) 'projectnodes',node,itreated(node)
                     exit
                  else
                     !
                     !                 update list; exit if an element in the list is
                     !                 revisited
                     !
                     !                      write(*,*) 'revisit needed!'
                     call nident(ilist,imastfa,nlist,id)
                     if(id.gt.0) then
                        if(ilist(id).eq.imastfa) then
                           write(*,*) '*WARNING in projectnodes:'
                           write(*,*) '         face in list revisited'
                           write(*,*) '         no projection for node',&
                                 node
                           itreated(node)=2
                           !                      write(*,*) 'projectnodes',node,itreated(node)
                           exit
                        endif
                     endif
                     nlist=nlist+1
                     do k=nlist,id+2,-1
                        ilist(k)=ilist(k-1)
                     enddo
                     ilist(id+1)=imastfa
                  !
                  endif
               enddo
            enddo
            !
            !           if no middle nodes are needed: proceed with the
            !           next face
            !
            if(iquad.ne.1) then
               index=ifac(4,index)
               cycle
            endif
            !
            !           quadratic: proceed with the middle nodes
            !           check each edge of the face whether it has already
            !           been treated (i.e. a middle node was generated,
            !           cf. field iedgmid)
            !
            do ii=1,3
               if(ii.eq.1) then
                  node1=ifac(1,index)
                  node2=ifac(2,index)
               elseif(ii.eq.2) then
                  node1=ifac(2,index)
                  node2=ifac(3,index)
               else
                  node1=ifac(1,index)
                  node2=ifac(3,index)
               endif
               !
               !              determine the number of the edge
               !
               if(node1.gt.node2) then
                  node=node2
                  node2=node1
                  node1=node
               endif
               indexe=ipoed(node1)
               do
                  if(iedg(2,indexe).eq.node2) exit
                  indexe=iedg(3,indexe)
               enddo
               !
               !              check whether a middle node was generated
               !              if so, nothing to do -> check the next edge
               !
               if(iedgmid(indexe).ne.0) cycle
               !
               !              generate a middle node
               !
               nktet=nktet+1
               iedgmid(indexe)=nktet
               do k=1,3
                  pnode(k)=(cotet(k,node1)+cotet(k,node2))/2.d0
                  cotet(k,nktet)=pnode(k)
               enddo
               node=nktet
               !                if(node.eq.731) write(*,*) 'projectnodes 1',node
               !
               !              parent face
               !
               imastfa=iexternfa(index)
               ilist(1)=imastfa
               !                if(node.eq.731) write(*,*) 'projectnodes 1a',imastfa
               nlist=1
               !
               !              start the loop looking for the correct face;
               !              starting with the parent face
               !
               do
                  if(ifacext(4,imastfa).ne.0) then
                     nterms=6
                     do j=1,6
                        do k=1,3
                           pneigh(k,j)=cotet(k,ifacext(j,imastfa))
                        enddo
                     enddo
                  !                      if(node.eq.731) write(*,*)
                  !      &                  'projectnodes 2',(ifacext(j,imastfa),j=1,6)
                  else
                     nterms=3
                     do j=1,3
                        do k=1,3
                           pneigh(k,j)=cotet(k,ifacext(j,imastfa))
                        enddo
                     enddo
                  endif
                  !
                  do k=1,3
                     pnodeproj(k)=pnode(k)
                  enddo
                  !                      if(node.eq.731) write(*,*)
                  !      &                  'projectnodes 3',(pnodeproj(k),k=1,3)
                  call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,&
                       xi,et)
                  !                      if(node.eq.731) write(*,*)
                  !      &                  'projectnodes 4',(pnodeproj(k),k=1,3)
                  !                      if(node.eq.731) write(*,*)
                  !      &                  'projectnodes 5',xi,et
                  !
                  !                 check whether this face is the correct one;
                  !                 if dabs(xi)=1 or dabs(et)=1 or xi+et=0 this
                  !                 may not be the case
                  !
                  !                 the solution is found (isol=1) unless proved
                  !                 otherwise
                  !
                  isol=1
                  !
                  if(dabs(et).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 1-2 unless sharp
                     !
                     imasted=ifacexted(1,imastfa)
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 6a',imasted,isharp(imasted)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 6',imastfa
                     endif
                  endif
                  !
                  if(dabs(xi+et-1.d0).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 2-3 unless sharp
                     !
                     imasted=ifacexted(2,imastfa)
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 7a',imasted,isharp(imasted)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 7',imastfa
                     endif
                  endif
                  !
                  if(dabs(xi).lt.1.d-10) then
                     !
                     !                    take neighboring face across edge 3-1 unless sharp
                     !
                     imasted=ifacexted(3,imastfa)
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 8a',imasted,isharp(imasted)
                     if(isharp(imasted).eq.0) then
                        if(iedgextfa(1,imasted).eq.imastfa) then
                           imastfa=iedgextfa(2,imasted)
                        else
                           imastfa=iedgextfa(1,imasted)
                        endif
                        isol=0
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 8',imastfa
                     endif
                  endif
                  !
                  !                 if solution is found: copy projected coordinates
                  !                 else continue with a neighbor
                  !
                  if(isol.eq.1) then
                     do k=1,3
                        cotet(k,node)=pnodeproj(k)
                     enddo
                     !                      if(node.eq.731) write(*,*)
                     !      &                  'projectnodes 9',(pnodeproj(k),k=1,3)
                     exit
                  else
                     !
                     !                 update list; exit if an element in the list is
                     !                 revisited
                     !
                     call nident(ilist,imastfa,nlist,id)
                     if(id.gt.0) then
                        if(ilist(id).eq.imastfa) then
                           write(*,*) '*WARNING in projectnodes:'
                           write(*,*) '         face in list revisited'
                           write(*,*) '         no projection for node',&
                                 node
                           exit
                        endif
                     endif
                     nlist=nlist+1
                     do k=nlist,id+2,-1
                        ilist(k)=ilist(k-1)
                     enddo
                     ilist(id+1)=imastfa
                  !
                  endif
               enddo
            !
            enddo
            !
            !           treat the next face
            !
            index=ifac(4,index)
         enddo
      enddo loop2
               
      ! !
      ! !     loop over all edges: projection on surfaces of the unrefined mesh
      ! !
      !       loop2: do i=1,nktet_
      !          index=ipoed(i)
      ! !
      !          do
      !             if(index.eq.0) cycle loop2
      ! !
      ! !           end nodes belonging to the edge
      ! !
      !             node1=iedg(1,index)
      !             node2=iedg(2,index)
      ! !
      !             if(iexternedg(index).lt.0) then
      ! !
      ! !              external edge not belonging to any external edge of
      ! !              the unrefined mesh: projection onto the external
      ! !              surface
      ! !
      !                do ii=1,3
      ! !
      ! !                 pnode(1..3) contains the coordinates before projection
      ! !
      !                   if(ii.eq.1) then
      ! c                     if(itreated(node1).eq.1) cycle
      !                      if(itreated(node1).ne.0) cycle
      !                      node=node1
      !                      do k=1,3
      !                         pnode(k)=cotet(k,node1)
      !                      enddo
      !                      itreated(node1)=2
      !                   elseif(ii.eq.2) then
      ! c                     if(itreated(node2).eq.1) cycle
      !                      if(itreated(node2).ne.0) cycle
      !                      node=node2
      !                      do k=1,3
      !                         pnode(k)=cotet(k,node2)
      !                      enddo
      !                      itreated(node2)=2
      !                   else
      !                      if(iquad.ne.1) cycle
      !                      nktet=nktet+1
      !                      iedgmid(index)=nktet
      !                      do k=1,3
      !                         pnode(k)=
      !      &                       (cotet(k,node1)+cotet(k,node2))/2.d0
      !                         cotet(k,nktet)=pnode(k)
      !                      enddo
      ! c                     if(nktet.eq.3366) then
      ! c                        write(*,*) 'projectnodes 3366 ',(pnode(k),k=1,3)
      ! c                     endif
      !                      node=nktet
      !                   endif
      ! !
      ! !                 determining the element in the unrefined mesh
      ! !                 to which the node belongs
      ! !
      !                   call basis(doubleglob(1),doubleglob(netet+1),
      !      &                 doubleglob(2*netet+1),
      !      &                 doubleglob(3*netet+1),doubleglob(4*netet+1),
      !      &                 doubleglob(5*netet+1),integerglob(6),
      !      &                 integerglob(netet+6),
      !      &                 integerglob(2*netet+6),doubleglob(6*netet+1),
      !      &                 integerglob(3*netet+6),nktri,netet,
      !      &                 doubleglob(4*nfaces+6*netet+1),nfield,
      !      &                 doubleglob(nktri+4*nfaces+6*netet+1),
      !      &                 integerglob(7*netet+6),integerglob(ne+7*netet+6),
      !      &                 integerglob(2*ne+7*netet+6),
      !      &                 integerglob(nkon+2*ne+7*netet+6),
      !      &                 pnode(1),pnode(2),pnode(3),value,ratio,iselect,
      !      &                 nselect,istartset,iendset,ialsete,imastset,
      !      &                 integerglob(nkon+2*ne+8*netet+6),nterms,konl,
      !      &                 nelem,loopa)
      ! c                  if(nktet.eq.3366) then
      ! c                     write(*,*) 'projectnodes 3366 nelem',nelem
      ! c                  endif
      ! c                  if(nktet.eq.2559) then
      ! c                     write(*,*) 'projectnodes 2559 nelem',nelem
      ! c                  endif
      ! !
      !                   indexf=3*netet+5+4*(nelem-1)
      ! !
      !                   distmin=1.d30
      ! !
      !                   do m=1,4
      !                      iface=abs(integerglob(indexf+m))
      !                      if(iexternfaor(iface).eq.0) cycle
      ! !
      ! !                    index into field ifacext
      ! !
      !                      imastfa=iexternfaor(iface)
      !                      if(ifacext(4,imastfa).ne.0) then
      !                         nterms=6
      !                         do j=1,6
      !                            do k=1,3
      !                               pneigh(k,j)=cotet(k,ifacext(j,imastfa))
      !                            enddo
      !                         enddo
      !                      else
      !                         nterms=3
      !                         do j=1,3
      !                            do k=1,3
      !                               pneigh(k,j)=cotet(k,ifacext(j,imastfa))
      !                            enddo
      !                         enddo
      !                      endif
      ! !
      !                      do k=1,3
      !                         pnodeproj(k)=pnode(k)
      !                      enddo
      !                      call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,
      !      &                           xi,et)
      !                      if(dist.lt.distmin) then
      !                      if(nktet.eq.3366) then
      ! c                   write(*,*) 'projectnodes 3366 ',(pnodeproj(k),k=1,3)
      ! c                   write(*,*) (ifacext(j,imastfa),j=1,6)
      !                      endif
      !                         do k=1,3
      !                            cotet(k,node)=pnodeproj(k)
      !                         enddo
      !                         distmin=dist
      !                      endif
      ! !
      !                   enddo
      !                enddo
      !             endif
      ! !
      !             index=iedg(3,index)
      !          enddo
      !       enddo loop2
      !
      !     loop over all edges: volumetric edges
      !
      loop3: do i=1,nktet_
         index=ipoed(i)
         !
         do
            if(index.eq.0) cycle loop3
            !
            !           end nodes belonging to the edge
            !
            node1=iedg(1,index)
            node2=iedg(2,index)
            !
            if(iexternedg(index).eq.0) then
               !
               !              internal edge: middle node (if necessary) in the
               !              middle of the edge
               !
               if(iquad.eq.1) then
                  nktet=nktet+1
                  iedgmid(index)=nktet
                  do k=1,3
                     cotet(k,nktet)=&
                          (cotet(k,node1)+cotet(k,node2))/2.d0
                  enddo
               endif
            endif
            !
            index=iedg(3,index)
         enddo
      enddo loop3
      !
      return
      end
