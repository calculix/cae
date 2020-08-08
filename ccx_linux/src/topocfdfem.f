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
      subroutine topocfdfem(nelemface,sideface,nface,ipoface,nodface,
     &  ne,ipkon,kon,lakon,ikboun,ilboun,xboun,nboun,nk,isolidsurf,
     &  nsolidsurf,ifreestream,nfreestream,neighsolidsurf,iponoel,inoel,
     &  inoelfree,nef,co,ipompc,nodempc,ikmpc,ilmpc,nmpc,set,istartset,
     &  iendset,ialset,nset,iturbulent,inomat,ielmat)
!
!     preliminary calculations for cfd applicatons:
!     - determining the external faces of the mesh and storing
!       them in fields nelemface and sideface
!     - determining the nodes belonging to solid surfaces and
!       storing them in isolidsurf (in ascending order)
!     - determining the nodes belonging to freestream surfaces
!       and storing them in ifreestream (in ascending order)
!     - determining the fluid elements belonging to a given node
!       and storing them in fields iponoel and inoel
!
      implicit none
!
      logical solidboun,solidmpc,freestream
!
      character*1 sideface(*)
      character*8 lakon(*)
      character*81 set(*),noset
!
      integer nelemface(*),nface,ipoface(*),nodface(5,*),nodes(4),
     &  ne,ipkon(*),kon(*),indexe,ifaceq(8,6),ifacet(7,4),index,
     &  ifacew(8,5),ithree,ifour,iaux,kflag,nnodes,ikboun(*),
     &  ilboun(*),nboun,isolidsurf(*),nsolidsurf,ifreestream(*),
     &  nfreestream,id,nk,node,idof,i,j,k,l,m,neighsolidsurf(*),
     &  iponoel(*),noden,idn,nope,nodemin,ifree,nef,indexold,
     &  inoel(3,*),ifreenew,inoelfree,ikmpc(*),nmpc,indexi,
     &  nodempc(3,*),ipompc(*),ilmpc(*),kmax,impc,idofi,idi,
     &  iturbulent,istartset(*),iendset(*),ialset(*),nset,inomat(*),
     &  ielmat(*)
!
      real*8 xboun(*),dist,distmin,co(3,*)
!
!     nodes belonging to the element faces
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
!
      kflag=1
      ithree=3
      ifour=4
!
!     determining the external element faces of the fluid mesh 
!     the faces are catalogued by the three lowes nodes numbers
!     in ascending order. ipoface(i) points to a face for which
!     node i is the lowest node and nodface(1,ipoface(i)) and
!     nodface(2,ipoface(i)) are the next lower ones. 
!     nodface(3,ipoface(i)) contains the element number,
!     nodface(4,ipoface(i)) the face number and nodface(5,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     An external element face is one which belongs to one element
!     only
!
      ifree=1
      do i=1,6*nef-1
         nodface(5,i)=i+1
      enddo
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'F') cycle
         indexe=ipkon(i)
         if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
            do j=1,6
               do k=1,4
                  nodes(k)=kon(indexe+ifaceq(k,j))
               enddo
               call isortii(nodes,iaux,ifour,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) then
            do j=1,4
               do k=1,3
                  nodes(k)=kon(indexe+ifacet(k,j))
               enddo
               call isortii(nodes,iaux,ithree,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         else
            do j=1,5
               if(j.le.2) then
                  do k=1,3
                     nodes(k)=kon(indexe+ifacew(k,j))
                  enddo
                  call isortii(nodes,iaux,ithree,kflag)
               else
                  do k=1,4
                     nodes(k)=kon(indexe+ifacew(k,j))
                  enddo
                  call isortii(nodes,iaux,ifour,kflag)
               endif
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         endif
      enddo
!
!     storing the external faces in nelemface and sideface
!     catalogueing the external nodes in isolidsurf and ifreestream
!
!     only the nodes which
!      - belong to external faces AND
!      - in which all velocity components are set
!        by SPC or MPC boundary conditions
!     are considered as solid surface nodes
!
!     all other external face nodes are freestream nodes
!
      nface=0
      nsolidsurf=0
      nfreestream=0
!
      do m=1,nk
         index=ipoface(m)
         do
            if(index.eq.0) exit
            nface=nface+1      
            i=nodface(3,index)
            j=nodface(4,index)
!
            nelemface(nface)=i
            write(sideface(nface)(1:1),'(i1)') j
!
c            indexe=ipkon(i)
c            if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
c               if(lakon(i)(4:4).eq.'2') then
c                  nnodes=8
c               else
c                  nnodes=4
c               endif
c               do k=1,nnodes
c                  node=kon(indexe+ifaceq(k,j))
c                  solidboun=.true.
c                  solidmpc=.true.
c                  loop1: do l=1,3
c                     idof=8*(node-1)+l
c                     call nident(ikboun,idof,nboun,id)
c                     if(id.le.0) then
c                        solidboun=.false.
cc                        exit
c                     elseif(ikboun(id).ne.idof) then
c                        solidboun=.false.
cc                        exit
c                     elseif(dabs(xboun(ilboun(id))).gt.1.d-20) then
c                        solidboun=.false.
cc                        exit
c                     endif
c!
c!                    if the degree of freedom was not fixed by a SPC
c!                    check whether it is fixed by a MPC
c!
c                     if(.not.solidboun) then
c                        call nident(ikmpc,idof,nmpc,id)
c                        if(id.le.0) then
c                           solidmpc=.false.
c                           exit
c                        elseif(ikmpc(id).ne.idof) then
c                           solidmpc=.false.
c                           exit
c                        else
c                           impc=ilmpc(id)
c                           indexi=nodempc(3,ipompc(impc))
c                           do
c                              if(indexi.eq.0) exit
c                              idofi=8*(nodempc(1,indexi)-1)+
c     &                                 nodempc(2,indexi)
c                              call nident(ikboun,idofi,nboun,idi)
c                              if(idi.le.0) then
c                                 solidmpc=.false.
c                                 exit loop1
c                              elseif(ikboun(idi).ne.idofi) then
c                                 solidmpc=.false.
c                                 exit loop1
c                              endif
c                              indexi=nodempc(3,indexi)
c                            enddo
c                         endif
c                      endif
c                  enddo loop1
c!
c                  if((solidboun).or.(solidmpc)) then
c                     call nident(isolidsurf,node,nsolidsurf,id)
c                     if(id.gt.0) then
c                        if(isolidsurf(id).eq.node) cycle
c                     endif
c                     nsolidsurf=nsolidsurf+1
c                     do l=nsolidsurf,id+2,-1
c                        isolidsurf(l)=isolidsurf(l-1)
c                     enddo
c                     isolidsurf(id+1)=node
c                  else
c                     call nident(ifreestream,node,nfreestream,id)
c                     if(id.gt.0) then
c                        if(ifreestream(id).eq.node) cycle
c                     endif
c                     nfreestream=nfreestream+1
c                     do l=nfreestream,id+2,-1
c                        ifreestream(l)=ifreestream(l-1)
c                     enddo
c                     ifreestream(id+1)=node
c                  endif
c               enddo
c            elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) 
c     &        then
c               if(lakon(i)(4:4).eq.'4') then
c                  nnodes=3
c               else
c                  nnodes=6
c               endif
c               do k=1,nnodes
c                  node=kon(indexe+ifacet(k,j))
c                  solidboun=.true.
c                  solidmpc=.true.
c                  loop2: do l=1,3
c                     idof=8*(node-1)+l
c                     call nident(ikboun,idof,nboun,id)
c                     if(id.le.0) then
c                        solidboun=.false.
cc                        exit
c                     elseif(ikboun(id).ne.idof) then
c                        solidboun=.false.
cc                        exit
c                     elseif(dabs(xboun(ilboun(id))).gt.1.d-20) then
c                        solidboun=.false.
cc                        exit
c                     endif
c!
c!                    if the degree of freedom was not fixed by a SPC
c!                    check whether it is fixed by a MPC
c!
c                     if(.not.solidboun) then
c                        call nident(ikmpc,idof,nmpc,id)
c                        if(id.le.0) then
c                           solidmpc=.false.
c                           exit
c                        elseif(ikmpc(id).ne.idof) then
c                           solidmpc=.false.
c                           exit
c                        else
c                           impc=ilmpc(id)
c                           indexi=nodempc(3,ipompc(impc))
c                           do
c                              if(indexi.eq.0) exit
c                              idofi=8*(nodempc(1,indexi)-1)+
c     &                                 nodempc(2,indexi)
c                              call nident(ikboun,idofi,nboun,idi)
c                              if(idi.le.0) then
c                                 solidmpc=.false.
c                                 exit loop2
c                              elseif(ikboun(idi).ne.idofi) then
c                                 solidmpc=.false.
c                                 exit loop2
c                              endif
c                              indexi=nodempc(3,indexi)
c                            enddo
c                         endif
c                      endif
c                  enddo loop2
c!
c                  if((solidboun).or.(solidmpc)) then
c                     call nident(isolidsurf,node,nsolidsurf,id)
c                     if(id.gt.0) then
c                        if(isolidsurf(id).eq.node) cycle
c                     endif
c                     nsolidsurf=nsolidsurf+1
c                     do l=nsolidsurf,id+2,-1
c                        isolidsurf(l)=isolidsurf(l-1)
c                     enddo
c                     isolidsurf(id+1)=node
c                  else
c                     call nident(ifreestream,node,nfreestream,id)
c                     if(id.gt.0) then
c                        if(ifreestream(id).eq.node) cycle
c                     endif
c                     nfreestream=nfreestream+1
c                     do l=nfreestream,id+2,-1
c                        ifreestream(l)=ifreestream(l-1)
c                     enddo
c                     ifreestream(id+1)=node
c                  endif
c               enddo
c            else
c               if(lakon(i)(4:4).eq.'6') then
c                  if(j.le.2) then
c                     nnodes=3
c                  else
c                     nnodes=4
c                  endif
c               else
c                  if(j.le.2) then
c                     nnodes=6
c                  else
c                     nnodes=8
c                  endif
c               endif
c               do k=1,nnodes
c                  node=kon(indexe+ifacew(k,j))
c                  solidboun=.true.
c                  solidmpc=.true.
c                  loop3: do l=1,3
c                     idof=8*(node-1)+l
c                     call nident(ikboun,idof,nboun,id)
c                     if(id.le.0) then
c                        solidboun=.false.
c                        exit
c                     elseif(ikboun(id).ne.idof) then
c                        solidboun=.false.
cc                        exit
c                     elseif(dabs(xboun(ilboun(id))).gt.1.d-20) then
c                        solidboun=.false.
cc                        exit
c                     endif
c!
c!                    if the degree of freedom was not fixed by a SPC
c!                    check whether it is fixed by a MPC
c!
c                     if(.not.solidboun) then
c                        call nident(ikmpc,idof,nmpc,id)
c                        if(id.le.0) then
c                           solidmpc=.false.
c                           exit
c                        elseif(ikmpc(id).ne.idof) then
c                           solidmpc=.false.
c                           exit
c                        else
c                           impc=ilmpc(id)
c                           indexi=nodempc(3,ipompc(impc))
c                           do
c                              if(indexi.eq.0) exit
c                              idofi=8*(nodempc(1,indexi)-1)+
c     &                                 nodempc(2,indexi)
c                              call nident(ikboun,idofi,nboun,idi)
c                              if(idi.le.0) then
c                                 solidmpc=.false.
c                                 exit loop3
c                              elseif(ikboun(idi).ne.idofi) then
c                                 solidmpc=.false.
c                                 exit loop3
c                              endif
c                              indexi=nodempc(3,indexi)
c                            enddo
c                         endif
c                      endif
c                  enddo loop3
c!
c                  if((solidboun).or.(solidmpc)) then
c                     call nident(isolidsurf,node,nsolidsurf,id)
c                     if(id.gt.0) then
c                        if(isolidsurf(id).eq.node) cycle
c                     endif
c                     nsolidsurf=nsolidsurf+1
c                     do l=nsolidsurf,id+2,-1
c                        isolidsurf(l)=isolidsurf(l-1)
c                     enddo
c                     isolidsurf(id+1)=node
c                  else
c                     call nident(ifreestream,node,nfreestream,id)
c                     if(id.gt.0) then
c                        if(ifreestream(id).eq.node) cycle
c                     endif
c                     nfreestream=nfreestream+1
c                     do l=nfreestream,id+2,-1
c                        ifreestream(l)=ifreestream(l-1)
c                     enddo
c                     ifreestream(id+1)=node
c                  endif
c               enddo
c            endif
            index=nodface(5,index)
         enddo
      enddo
!
!     storing the nodes of the solid surfaces
!
      noset(1:13)='SOLIDSURFACEN'
      do i=1,nset
         if(set(i)(1:13).eq.noset(1:13)) exit
      enddo
      if((i.gt.nset).and.(iturbulent.gt.0)) then
         write(*,*) '*WARNING in precfd: node set SOLID SURFACE '
         write(*,*) '         has not been defined. This set may'
         write(*,*) '         be needed in a turbulent calculation'
      elseif(i.le.nset) then
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               nsolidsurf=nsolidsurf+1
               isolidsurf(nsolidsurf)=ialset(j)
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  nsolidsurf=nsolidsurf+1
                  isolidsurf(nsolidsurf)=k
               enddo
            endif
         enddo
         call isortii(isolidsurf,iaux,nsolidsurf,kflag)
      endif
!
!     storing the nodes of freestream surfaces
!
      noset(1:18)='FREESTREAMSURFACEN'
      do i=1,nset
         if(set(i)(1:18).eq.noset(1:18)) exit
      enddo
      if((i.gt.nset).and.(iturbulent.gt.0)) then
         write(*,*) '*WARNING in precfd: node set FREESTREAM SURFACE '
         write(*,*) '         has not been defined. This set may'
         write(*,*) '         be needed in a turbulent calculation'
      elseif(i.le.nset) then
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               nfreestream=nfreestream+1
               ifreestream(nfreestream)=ialset(j)
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  nfreestream=nfreestream+1
                  ifreestream(nfreestream)=k
               enddo
            endif
         enddo
         call isortii(ifreestream,iaux,nfreestream,kflag)
      endif
!
!     all nodes belonging to MPC's are removed from the
!     ifreestream stack
!
c      do i=1,nmpc
c         index=ipompc(i)
c         do
c            if(index.eq.0) exit
c            node=nodempc(1,index)
c            call nident(ifreestream,node,nfreestream,id)
c            if(id.gt.0) then
c               if(ifreestream(id).eq.node) then
c                  nfreestream=nfreestream-1
c                  do j=id,nfreestream
c                     ifreestream(j)=ifreestream(j+1)
c                  enddo
c               endif
c            endif
c            index=nodempc(3,index)
c         enddo
c      enddo
!
!     reject external faces which contain end nodes which belong
!     neither to solid surfaces nor to freestream faces
!     (basically containing nodes which are part of MPC's)
!
c      do m=1,nface
c         i=nelemface(m)
c         indexe=ipkon(i)
c         read(sideface(m)(1:1),'(i1)') j
c         if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
c            do k=1,4
c               freestream=.true.
c               node=kon(indexe+ifaceq(k,j))
c               call nident(ifreestream,node,nfreestream,id)
c               if(id.le.0) then
c                  freestream=.false.
c               elseif(ifreestream(id).ne.node) then
c                  freestream=.false.
c               endif
c               if(.not.freestream) then
c                  call nident(isolidsurf,node,nsolidsurf,id)
c                  if(id.le.0) then
c                     nelemface(m)=0
c                     exit
c                  elseif(isolidsurf(id).ne.node) then
c                     nelemface(m)=0
c                     exit
c                  endif
c               endif
c            enddo
c         elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) then
c            do k=1,3
c               node=kon(indexe+ifacet(k,j))
c               call nident(ifreestream,node,nfreestream,id)
c               if(id.le.0) then
c                  freestream=.false.
c               elseif(ifreestream(id).ne.node) then
c                  freestream=.false.
c               endif
c               if(.not.freestream) then
c                  call nident(isolidsurf,node,nsolidsurf,id)
c                  if(id.le.0) then
c                     nelemface(m)=0
c                     exit
c                  elseif(isolidsurf(id).ne.node) then
c                     nelemface(m)=0
c                     exit
c                  endif
c               endif
c            enddo
c         else
c            if(j.le.2) then
c               kmax=3
c            else
c               kmax=4
c            endif
c            do k=1,kmax
c               node=kon(indexe+ifacew(k,j))
c               call nident(ifreestream,node,nfreestream,id)
c               if(id.le.0) then
c                  freestream=.false.
c               elseif(ifreestream(id).ne.node) then
c                  freestream=.false.
c               endif
c               if(.not.freestream) then
c                  call nident(isolidsurf,node,nsolidsurf,id)
c                  if(id.le.0) then
c                     nelemface(m)=0
c                     exit
c                  elseif(isolidsurf(id).ne.node) then
c                     nelemface(m)=0
c                     exit
c                  endif
c               endif
c            enddo
c         endif
c      enddo
c!
c!     remove the faces of zero elements
c!
c      i=0
c      do m=1,nface
c         if(nelemface(m).ne.0) then
c            i=i+1
c            nelemface(i)=nelemface(m)
c            sideface(i)=sideface(m)
c         endif
c      enddo
c      nface=i
!
!     storing the in-stream neighbors of the solid surface external
!       nodes in neighsolidsurf
!
      do m=1,nface
         i=nelemface(m)
         read(sideface(m)(1:1),'(i1)') j
         indexe=ipkon(i)
!
         if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
            if(lakon(i)(4:4).eq.'2') then
               nnodes=8
               nope=20
            else
               nnodes=4
               nope=8
            endif
            do k=1,nnodes
               node=kon(indexe+ifaceq(k,j))
!
!              node must belong to solid surface
!
               call nident(isolidsurf,node,nsolidsurf,id)
               if(id.le.0) then
                  cycle
               elseif(isolidsurf(id).ne.node) then
                  cycle
               endif
!
!              check whether neighbor was already found
!
               if(neighsolidsurf(id).ne.0) cycle
!               
               distmin=1.d30
               nodemin=0
!
               do l=1,nope
                  noden=kon(indexe+l)
!
!                 node must not belong to solid surface
!
                  call nident(isolidsurf,noden,nsolidsurf,idn)
                  if(idn.gt.0) then
                     if(isolidsurf(idn).eq.noden) cycle
                  endif
                  dist=dsqrt((co(1,node)-co(1,noden))**2+
     &                       (co(2,node)-co(2,noden))**2+
     &                       (co(3,node)-co(3,noden))**2)
                  if(dist.lt.distmin) then
                     distmin=dist
                     nodemin=noden
                  endif
               enddo
               if(nodemin.ne.0) then
                  neighsolidsurf(id)=nodemin
               endif
            enddo
         elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) 
     &           then
            if(lakon(i)(4:4).eq.'4') then
               nnodes=3
               nope=4
            else
               nnodes=6
               nope=10
            endif
            do k=1,nnodes
               node=kon(indexe+ifacet(k,j))
!
!              node must belong to solid surface
!
               call nident(isolidsurf,node,nsolidsurf,id)
               if(id.le.0) then
                  cycle
               elseif(isolidsurf(id).ne.node) then
                  cycle
               endif
!
!              check whether neighbor was already found
!
               if(neighsolidsurf(id).ne.0) cycle
!               
               distmin=1.d30
               nodemin=0
!
               do l=1,nope
                  noden=kon(indexe+l)
!
!                 node must not belong to solid surface
!
                  call nident(isolidsurf,noden,nsolidsurf,idn)
                  if(idn.gt.0) then
                     if(isolidsurf(idn).eq.noden) cycle
                  endif
                  dist=dsqrt((co(1,node)-co(1,noden))**2+
     &                       (co(2,node)-co(2,noden))**2+
     &                       (co(3,node)-co(3,noden))**2)
                  if(dist.lt.distmin) then
                     distmin=dist
                     nodemin=noden
                  endif
               enddo
               if(nodemin.ne.0) then
                  neighsolidsurf(id)=nodemin
               endif
            enddo
         else
            if(lakon(i)(4:4).eq.'6') then
               nope=6
               if(j.le.2) then
                  nnodes=3
               else
                  nnodes=4
               endif
            else
               nope=15
               if(j.le.2) then
                  nnodes=6
               else
                  nnodes=8
               endif
            endif
            do k=1,nnodes
               node=kon(indexe+ifacew(k,j))
 !
!              node must belong to solid surface
!
               call nident(isolidsurf,node,nsolidsurf,id)
               if(id.le.0) then
                  cycle
               elseif(isolidsurf(id).ne.node) then
                  cycle
               endif
!
!              check whether neighbor was already found
!
               if(neighsolidsurf(id).ne.0) cycle
!               
               distmin=1.d30
               nodemin=0
!
               do l=1,nope
                  noden=kon(indexe+l)
!
!                 node must not belong to solid surface
!
                  call nident(isolidsurf,noden,nsolidsurf,idn)
                  if(idn.gt.0) then
                     if(isolidsurf(idn).eq.noden) cycle
                  endif
                  dist=dsqrt((co(1,node)-co(1,noden))**2+
     &                       (co(2,node)-co(2,noden))**2+
     &                       (co(3,node)-co(3,noden))**2)
                  if(dist.lt.distmin) then
                     distmin=dist
                     nodemin=noden
                  endif
               enddo
               if(nodemin.ne.0) then
                  neighsolidsurf(id)=nodemin
               endif
           enddo
         endif
      enddo 
!
!     determining the fluid elements belonging to edge nodes of
!     the elements
!
      inoelfree=1
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'F') cycle
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            nope=15
         endif
         indexe=ipkon(i)
         do j=1,nope
            node=kon(indexe+j)
            inoel(1,inoelfree)=i
            inoel(2,inoelfree)=j
            inoel(3,inoelfree)=iponoel(node)
            iponoel(node)=inoelfree
            inoelfree=inoelfree+1
         enddo
      enddo
!
!     sorting nelemface
!
      kflag=2
      call isortic(nelemface,sideface,nface,kflag)
!
!     filling inomat: asigns a material to fluid nodes. 
!     (a fluid nodes is not assumed to be part of two
!      different fluids)
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'F') cycle
!
         indexe=ipkon(i)
         read(lakon(i)(4:4),'(i1)')nope
!
         do j=1,nope
            inomat(kon(indexe+j))=ielmat(i)
         enddo
      enddo
!
c      write(*,*) 'nfreestream ',nfreestream
c      do i=1,nfreestream
c         write(*,*) 'nfreestream ',i,ifreestream(i)
c      enddo
c      write(*,*) 'nsolidsurf ',nsolidsurf
c      do i=1,nsolidsurf
c         write(*,*) 'nsolidsurf ',i,isolidsurf(i),neighsolidsurf(i)
c      enddo
c      write(*,*) 'external faces'
c      do i=1,nface
c         write(*,*) nelemface(i),sideface(i)
c      enddo
!
      return
      end
