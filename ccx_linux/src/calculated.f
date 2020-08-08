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
      subroutine calculated(nktet,d,dmin,ipoed,iedg,cotet)
!
!     determine the length of all edges in the actual mesh 
!
      implicit none
!
      integer i,nktet,ipoed(*),iedg(3,*),index,n1,n2
!
      real*8 d(*),dmin,cotet(3,*)
!
!
!
!     determine the size of all edges
!
      dmin=1.d30
!
      loop: do i=1,nktet
         index=ipoed(i)
         do
            if(index.eq.0) cycle loop
!
            n1=iedg(1,index)
            n2=iedg(2,index)
!
            d(index)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &                     (cotet(2,n1)-cotet(2,n2))**2+
     &                     (cotet(3,n1)-cotet(3,n2))**2)
!
            if(d(index).lt.dmin) dmin=d(index)
!
            index=iedg(3,index)
         enddo
      enddo loop
c!
c!     calculating the desired edge size through interpolation
c!
c!     initializing fields
c!
c      nktri=integerglob(1)
c      netet=integerglob(2)
c      ne=integerglob(3)
c      nkon=integerglob(4)
c      nfaces=integerglob(5)
c      nfield=1
c      nselect=1
c      iselect(1)=1
c      imastset=0
c!
c      do i=1,nktet
c         if(ipoeln(i).eq.0) cycle
c!
c!        perform the interpolation for the internal node
c!
c         do j=1,3
c            coords(j)=cotet(j,i)
c         enddo
c         call basis(doubleglob(1),doubleglob(netet+1),
c     &        doubleglob(2*netet+1),
c     &        doubleglob(3*netet+1),doubleglob(4*netet+1),
c     &        doubleglob(5*netet+1),integerglob(6),
c     &        integerglob(netet+6),
c     &        integerglob(2*netet+6),doubleglob(6*netet+1),
c     &        integerglob(3*netet+6),nktri,netet,
c     &        doubleglob(4*nfaces+6*netet+1),nfield,
c     &        doubleglob(nktri+4*nfaces+6*netet+1),
c     &        integerglob(7*netet+6),integerglob(ne+7*netet+6),
c     &        integerglob(2*ne+7*netet+6),
c     &        integerglob(nkon+2*ne+7*netet+6),
c     &        coords(1),coords(2),coords(3),h(i),ratio,iselect,
c     &        nselect,istartset,iendset,ialset,imastset,
c     &        integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem)
c!
c      enddo
!
      return
      end

