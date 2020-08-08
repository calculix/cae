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
      subroutine getnodesinitetmesh(ne,lakon,ipkon,kon,istartset,
     &     iendset,ialset,set,nset,filab,inodestet,nnodestet,
     &     nboun,nodeboun,nforc,nodeforc,impctet)
!     
!     reading the initial tet mesh which will be refined 
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 set(*),elset
      character*87 filab(*)
!     
      integer ipkon(*),kon(*),istartset(*),iendset(*),ialset(*),
     &     inodestet(*),nnodestet,i,j,k,m,n,node,ne,nset,indexe,id,
     &     nboun,nodeboun(*),nforc,nodeforc(2,*),impctet(*)
!     
!     identify the set name/s if they exist
!     
      read(filab(48)(27:87),'(a61)') elset(1:61)
!     
      do i=62,81
        elset(i:i)=' '
      enddo
!     
      do i=1,nset
        if(set(i).eq.elset)exit
      enddo   
!     
      if(i.le.nset) then     
        do j=istartset(i),iendset(i)
          if(ialset(j).gt.0) then
            k=ialset(j)
!
!     the elements belonging to the unrefined mesh have been deactivated,
!     i.e. the transformation ipkon(k)=-2-ipkon(k) was performed in
!     modelchanges.f
!
            indexe=-2-ipkon(k)
            if(indexe.lt.0) cycle
            if(lakon(k)(1:4).eq.'C3D4') then
              do n=1,4
                node=kon(indexe+n)
                call nident(inodestet,node,nnodestet,id)
                if(id.gt.0) then
                  if(inodestet(id).eq.node) cycle
                endif
                nnodestet=nnodestet+1
                do m=nnodestet,id+2,-1
                  inodestet(m)=inodestet(m-1)
                enddo
                inodestet(id+1)=node
              enddo
            elseif(lakon(k)(1:5).eq.'C3D10') then
              do n=1,10
                node=kon(indexe+n)
                call nident(inodestet,node,nnodestet,id)
                if(id.gt.0) then
                  if(inodestet(id).eq.node) cycle
                endif
                nnodestet=nnodestet+1
                do m=nnodestet,id+2,-1
                  inodestet(m)=inodestet(m-1)
                enddo
                inodestet(id+1)=node
              enddo
            endif
          endif
        enddo
      else
        do j=1,ne
          k=j
!
!     the elements belonging to the unrefined mesh have been deactivated,
!     i.e. the transformation ipkon(k)=-2-ipkon(k) was performed in
!     modelchanges.f
!
          indexe=-2-ipkon(k)
          if(indexe.lt.0) cycle
          if(lakon(k)(1:4).eq.'C3D4') then
            do n=1,4
              node=kon(indexe+n)
              call nident(inodestet,node,nnodestet,id)
c     call nident(nodeboun,node,nboun,idboun)
              if(id.gt.0) then
                if(inodestet(id).eq.node) cycle
              endif
              nnodestet=nnodestet+1
              do m=nnodestet,id+2,-1
                inodestet(m)=inodestet(m-1)
              enddo
              inodestet(id+1)=node
            enddo
          elseif(lakon(k)(1:5).eq.'C3D10') then
            do n=1,10
              node=kon(indexe+n)
              call nident(inodestet,node,nnodestet,id)
              if(id.gt.0) then
                if(inodestet(id).eq.node) cycle
              endif
              nnodestet=nnodestet+1
              do m=nnodestet,id+2,-1
                inodestet(m)=inodestet(m-1)
              enddo
              inodestet(id+1)=node
            enddo
          endif
        enddo      
      endif
!     
!     node in the old mesh in which a SPC is defined       
!     
      do i=1,nboun
        node=nodeboun(i)
        call nident(inodestet,node,nnodestet,id)
        if(id.gt.0) then
          if(inodestet(id).eq.node) then
            impctet(id)=1
          endif
        endif
      enddo
!     
!     node in which a point force is defined
!     
      do i=1,nforc
        node=nodeforc(1,i)
        call nident(inodestet,node,nnodestet,id)
        if(id.gt.0) then
          if(inodestet(id).eq.node) then
            impctet(id)=1
          endif
        endif
      enddo
!     
!     keeping the nodes in which a SPC or point force was defined
!     
      j=0
      do i=1,nnodestet
        if(impctet(i).eq.1) then
          j=j+1
          inodestet(j)=inodestet(i)
        endif
      enddo
      nnodestet=j
!     
      return
      end

