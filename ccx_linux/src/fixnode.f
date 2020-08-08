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
      subroutine fixnode(nobject,nk,set,nset,istartset,iendset,
     &     ialset,iobject,nodedesiinv,dgdxglob,objectset)                       
!     
!     determination of the fixed nodes for the sensitivity analysis     
!     
      implicit none
!     
      character*81 objectset(4,*),set(*)
!     
      integer nk,nobject,nset,istartset(*),iendset(*),ialset(*),
     &     iobject,nodedesiinv(*),i,j,node
!     
      real*8 dgdxglob(2,nk,nobject)
!     
!     determining the set of fixed nodes
!     
      do i=1,nset
        if(objectset(3,iobject).eq.set(i)) exit
      enddo
!     
      if(i.le.nset) then
!     
        do j=istartset(i),iendset(i)
          if(ialset(j).gt.0) then
            node=ialset(j)
            if(nodedesiinv(node).eq.1) then
              dgdxglob(1,node,iobject)=1.0d0
              dgdxglob(2,node,iobject)=1.0d0
            endif
          else
            node=ialset(j-2)
            do
              node=node-ialset(j)
              if(node.ge.ialset(j-1)) exit
              if(nodedesiinv(node).eq.1) then
                dgdxglob(1,node,iobject)=1.0d0
                dgdxglob(2,node,iobject)=1.0d0
              endif
            enddo
          endif
        enddo
      endif
!     
      return        
      end
