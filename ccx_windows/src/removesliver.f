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
      subroutine removesliver(netet_,kontet,iexternnode,iedtet,
     &     iexternedg,quality)
!     
!     removes the slivers on the surface of a tetrahedral mesh
!     
      implicit none
!
      integer i,j,m,isliver,index,icount,kontet(4,*),iexternedg(*),
     &     iexternnode(*),nodes(4),netet_,iedtet(6,*)
!
      real*8 tol,quality(*)
!
      isliver=0
!
      do m=1,2
        loop1: do i=1,netet_
!
!     if a tet has exactly 5 external edges AND has a quality
!     (aspect ratio) bigger than the desired value THEN this element
!     is removed from the mesh
!
          if(kontet(1,i).eq.0) cycle
!     
          do j=1,4
            nodes(j)=kontet(j,i)
            if(iexternnode(nodes(j)).eq.0) cycle loop1
          enddo
!     
          tol=10.d0
!     
          if(quality(i).gt.tol) then
!     
            icount=0
            do j=1,6
              index=iedtet(j,i)
              if(iexternedg(index).ne.0) icount=icount+1
            enddo
!     
            if(((m.eq.1).and.(icount.eq.5)).or.
     &           ((m.eq.2).and.(icount.eq.4))) then
              write(*,*) 'removesliver sliver element found ',i
              isliver=isliver+1
              kontet(1,i)=0
            endif
!     
          endif
!     
        enddo loop1
      enddo
!
      write(*,*) 'Total number of sliver elements: ',isliver
      write(*,*)
!
      return
      end
