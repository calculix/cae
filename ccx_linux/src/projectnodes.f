!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine projectnodes(nktet_,ipoed,iedgmid,iexternedg,iedgext,
     &     cotet,nktet,iedg,iquad,iexternfa,ifacext,itreated,ilist,
     &     isharp,ipofa,ifac,iedgextfa,ifacexted,jfix,co,idimsh)
!     
!     1. projects nodes lying on external edges of the unrefined mesh
!     on those edges
!     
!     2. generates middles nodes if necessary
!     
      implicit none
!     
      integer nktet_,ipoed(*),index,i,j,k,iedgmid(*),iexternedg(*),
     &     imasted,iedgext(3,*),nterms,node1,node2,iedg(3,*),iquad,
     &     nktet,imastfa,ifacext(6,*),iexternfa(*),itreated(*),ii,node,
     &     isharp(*),ilist(*),nlist,ifac(4,*),iedgextfa(2,*),
     &     ifacexted(3,*),id,indexe,isol,ipofa(*),jfix(*),idimsh(*)
!     
      real*8 pneigh(3,9),cotet(3,*),pnode(3),ratio(9),dist,xi,et,
     &     pnodeproj(3),co(3,*)
!     
!     
!
!     set all nodes on "non-treated"
!      
      do i=1,nktet
        itreated(i)=0
      enddo
!     
!     loop over all edges: projection on sharp edges from the unrefined mesh 
!     
      loop1: do i=1,nktet_
        index=ipoed(i)
!     
        do
          if(index.eq.0) cycle loop1
!     
          if(iexternedg(index).gt.0) then
!     
!     recovering master edge
!     
            imasted=iexternedg(index)
!
!     if the master edge is not sharp: loop
!
            if(isharp(imasted).ne.1) then
              index=iedg(3,index)
              cycle
            endif
!     
!     edge is a subset of an external sharp edge
!     of the unrefined mesh
!     
!     end nodes belonging to the edge
!     
            node1=iedg(1,index)
            node2=iedg(2,index)
!     
!     check whether the external edge has a middle node
!     
            if(iedgext(2,imasted).ne.0) then
              nterms=3
              do j=1,3
                do k=1,3
                  pneigh(k,j)=cotet(k,iedgext(j,imasted))
                enddo
              enddo
              if(idimsh(node1).ne.2) itreated(node1)=1
              if(idimsh(node2).ne.2) itreated(node2)=1
            else
              nterms=2
            endif
!     
!     projection is only needed if the master edge is
!     quadratic
!     
            if(nterms.eq.3) then
!     
!     attach the first end node of the edge
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
              endif
!     
!     attach the other end node of the edge
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
              endif
            else
              itreated(node1)=1
              itreated(node2)=1
            endif
!     
!     create a middle node (if necessary) and attach it
!     (if necessary) 
!     
            if(iquad.eq.1) then
              if((jfix(node1).eq.1).and.(jfix(node2).eq.1)) then
                iedgmid(index)=iedgext(2,iexternedg(index))
              else
                nktet=nktet+1
                iedgmid(index)=nktet
!     
                if(nterms.eq.3) then
                  do k=1,3
                    pnode(k)=(cotet(k,node1)+cotet(k,node2))/2.d0
                  enddo
!     
!     projection for quadratic master edges
!     
                  call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
!     
                  do k=1,3
                    cotet(k,nktet)=pnode(k)
                  enddo
                else
                  do k=1,3
                    cotet(k,nktet)=
     &                   (cotet(k,node1)+cotet(k,node2))/2.d0
                  enddo
                endif
              endif
            endif
          endif
!     
          index=iedg(3,index)
        enddo
      enddo loop1
!     
!     projection of the external nodes not treated yet onto the
!     faces of the unrefined mesh
!     
      loop2: do i=1,nktet_
!     
!     loop over all faces of the refined mesh
!     
        index=ipofa(i)
        do
          if(index.eq.0) cycle loop2
!     
!     if no external face: loop
!     
          if(iexternfa(index).le.0) then
            index=ifac(4,index)
            cycle
          endif
!     
!     external face; treat the vertex nodes of the face
!     
          do ii=1,3
            node=ifac(ii,index)
            if(itreated(node).ne.0) cycle
!     
!     parent face
!     
            imastfa=iexternfa(index)
            ilist(1)=imastfa
            nlist=1
            do k=1,3
              pnode(k)=cotet(k,node)
            enddo
!     
!     start the loop looking for the correct face; 
!     starting with the parent face 
!     
            do
              if(ifacext(4,imastfa).ne.0) then
                nterms=6
                do j=1,6
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              else
                nterms=3
                do j=1,3
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              endif
!     
              do k=1,3
                pnodeproj(k)=pnode(k)
              enddo
              call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,
     &             xi,et)
!     
!     check whether this face is the correct one;
!     if dabs(xi)=1 or dabs(et)=1 or xi+et=0 this
!     may not be the case
!     
!     the solution is found (isol=1) unless proved
!     otherwise
!     
              isol=1
!     
              if(dabs(et).lt.1.d-10) then
!     
!     take neighboring face across edge 1-2 unless sharp
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
!     take neighboring face across edge 2-3 unless sharp
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
!     take neighboring face across edge 3-1 unless sharp
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
!     if solution is found: copy projected coordinates
!     else continue with a neighbor
!     
              if(isol.eq.1) then
                do k=1,3
                  cotet(k,node)=pnodeproj(k)
                enddo
                itreated(node)=2
                exit
              else
!     
!     update list; exit if an element in the list is
!     revisited
!     
                call nident(ilist,imastfa,nlist,id)
                if(id.gt.0) then
                  if(ilist(id).eq.imastfa) then
                    write(*,*) '*WARNING in projectnodes:'
                    write(*,*) '         face in list revisited'
                    write(*,*) '         no projection for node',
     &                   node
                    itreated(node)=2
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
!     if no middle nodes are needed: proceed with the
!     next face
!     
          if(iquad.ne.1) then
            index=ifac(4,index)
            cycle
          endif
!     
!     quadratic: proceed with the middle nodes
!     check each edge of the face whether it has already
!     been treated (i.e. a middle node was generated,
!     cf. field iedgmid)
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
!     determine the number of the edge
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
!     check whether a middle node was generated
!     if so, nothing to do -> check the next edge
!     
            if(iedgmid(indexe).ne.0) cycle
!     
!     generate a middle node
!
            if((jfix(node1).eq.1).and.(jfix(node2).eq.1)) then
              iedgmid(indexe)=iedgext(2,iexternedg(indexe))
              cycle
            else
              nktet=nktet+1
              iedgmid(indexe)=nktet
              do k=1,3
                pnode(k)=(cotet(k,node1)+cotet(k,node2))/2.d0
                cotet(k,nktet)=pnode(k)
              enddo
              node=nktet
            endif
!     
!     parent face
!     
            imastfa=iexternfa(index)
            ilist(1)=imastfa
            nlist=1
!     
!     start the loop looking for the correct face; 
!     starting with the parent face 
!     
            do
              if(ifacext(4,imastfa).ne.0) then
                nterms=6
                do j=1,6
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              else
                nterms=3
                do j=1,3
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              endif
!     
              do k=1,3
                pnodeproj(k)=pnode(k)
              enddo
              call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,
     &             xi,et)
!     
!     check whether this face is the correct one;
!     if dabs(xi)=1 or dabs(et)=1 or xi+et=0 this
!     may not be the case
!     
!     the solution is found (isol=1) unless proved
!     otherwise
!     
              isol=1
!     
              if(dabs(et).lt.1.d-10) then
!     
!     take neighboring face across edge 1-2 unless sharp
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
!     take neighboring face across edge 2-3 unless sharp
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
!     take neighboring face across edge 3-1 unless sharp
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
!     if solution is found: copy projected coordinates
!     else continue with a neighbor
!     
              if(isol.eq.1) then
                do k=1,3
                  cotet(k,node)=pnodeproj(k)
                enddo
                exit
              else
!     
!     update list; exit if an element in the list is
!     revisited
!     
                call nident(ilist,imastfa,nlist,id)
                if(id.gt.0) then
                  if(ilist(id).eq.imastfa) then
                    write(*,*) '*WARNING in projectnodes:'
                    write(*,*) '         face in list revisited'
                    write(*,*) '         no projection for node',
     &                   node
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
!     treat the next face
!     
          index=ifac(4,index)
        enddo
      enddo loop2
!     
!     loop over all edges: volumetric edges (only if middles nodes are to be created)
!     
      if(iquad.eq.1) then
        loop3: do i=1,nktet_
          index=ipoed(i)
!     
          do
            if(index.eq.0) cycle loop3
!     
            if(iexternedg(index).eq.0) then
!     
!     end nodes belonging to the edge
!     
              node1=iedg(1,index)
              node2=iedg(2,index)
!     
!     internal edge: middle node (if necessary) in the 
!     middle of the edge
!     
              nktet=nktet+1
              iedgmid(index)=nktet
              do k=1,3
                cotet(k,nktet)=
     &               (cotet(k,node1)+cotet(k,node2))/2.d0
              enddo
            endif
!     
            index=iedg(3,index)
          enddo
        enddo loop3
      endif
!     
      return
      end
