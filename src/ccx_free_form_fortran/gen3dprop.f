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
      subroutine gen3dprop(prop,ielprop,iponoel,inoel,iponoelmax,kon,&
        ipkon,lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,&
        nmpc_,mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,&
        nk_,co,nmethod,iperturb)
      !
      !     connects nodes of 1-D and 2-D elements which are used in fluid
      !     property definitions to the nodes of their expanded counterparts
      !
      implicit none
      !
      character*8 lakon(*)
      character*20 labmpc(*)
      !
      integer iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),ne,&
        iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,&
        ikmpc(*),ilmpc(*),rig(*),ntrans,inotr(2,*),i,node,ielprop(*),&
        index,ielem,j,indexe,indexk,idir,nk,nk_,&
        newnode,idof,id,mpcfreenew,k,nam,nmethod,iperturb,ii
      !
      real*8 xnor(*),coefmpc(*),trab(7,*),co(3,*),prop(*)
      !
      do i=1,ne
         !          if((lakon(i).ne.'DLIPIMAF').and.(lakon(i).ne.'DLIPIWCF')) cycle
         if((lakon(i).ne.'DLIPIMAF').and.(lakon(i).ne.'DLIPIWCF')&
              .and.(lakon(i)(1:5).ne.'DLABF')&
              .and.(lakon(i)(1:6).ne.'DGAPFF')&
              .and.(lakon(i)(1:5).ne.'DORFL')&
              .and.(lakon(i)(1:6).ne.'DGAPIF')) cycle
         do ii=1,6
            node=int(prop(ielprop(i)+int((ii+2.5d0)/3.d0)))
            if(node.gt.iponoelmax) cycle
            index=iponoel(node)
            if(index.eq.0) cycle
            ielem=inoel(1,index)
            j=inoel(2,index)
            indexe=ipkon(ielem)
            indexk=iponor(2,indexe+j)
            idir=ii-3*(int((ii+2.5d0)/3.d0)-1)
            !             write(*,*) 'gen3dprop,node,idir',node,idir
            !
            if(rig(node).ne.0) cycle
            !
            !     2d element shell element: generate MPC's
            !
            if(lakon(ielem)(7:7).eq.'L') then
               newnode=knor(indexk+1)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=knor(indexk+3)
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
            elseif(lakon(ielem)(7:7).eq.'B') then
               !
               !                       1d beam element: generate MPC's
               !
               newnode=knor(indexk+1)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase memmpc_'
                     call exit(201)
                  endif
                  do k=2,4
                     nodempc(1,mpcfree)=knor(indexk+k)
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*)&
                             '*ERROR in gen3dprop: increase memmpc_'
                        call exit(201)
                     endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dprop: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
            else
               !
               !                       2d plane stress, plane strain or axisymmetric
               !                       element: SPC
               !
               newnode=knor(indexk+2)
               idof=8*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if(((id.le.0).or.(ikmpc(id).ne.idof)).and.&
                    (idir.ne.3)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*)&
                          '*ERROR in gen3dmpc: increase nmpc_'
                     call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dmpc: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-1.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*)&
                          '*ERROR in gen3dmpc: increase memmpc_'
                     call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
            endif
         enddo
      enddo
      !
      return
      end


