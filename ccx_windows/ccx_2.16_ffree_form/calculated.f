!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
      intent(in) nktet,ipoed,iedg,cotet
      !
      intent(inout) d,dmin
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
            d(index)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+&
                           (cotet(2,n1)-cotet(2,n2))**2+&
                           (cotet(3,n1)-cotet(3,n2))**2)
            !
            if(d(index).lt.dmin) dmin=d(index)
            !
            index=iedg(3,index)
         enddo
      enddo loop
      ! !
      ! !     calculating the desired edge size through interpolation
      ! !
      ! !     initializing fields
      ! !
      !       nktri=integerglob(1)
      !       netet=integerglob(2)
      !       ne=integerglob(3)
      !       nkon=integerglob(4)
      !       nfaces=integerglob(5)
      !       nfield=1
      !       nselect=1
      !       iselect(1)=1
      !       imastset=0
      ! !
      !       do i=1,nktet
      !          if(ipoeln(i).eq.0) cycle
      ! !
      ! !        perform the interpolation for the internal node
      ! !
      !          do j=1,3
      !             coords(j)=cotet(j,i)
      !          enddo
      !          call basis(doubleglob(1),doubleglob(netet+1),
      !      &        doubleglob(2*netet+1),
      !      &        doubleglob(3*netet+1),doubleglob(4*netet+1),
      !      &        doubleglob(5*netet+1),integerglob(6),
      !      &        integerglob(netet+6),
      !      &        integerglob(2*netet+6),doubleglob(6*netet+1),
      !      &        integerglob(3*netet+6),nktri,netet,
      !      &        doubleglob(4*nfaces+6*netet+1),nfield,
      !      &        doubleglob(nktri+4*nfaces+6*netet+1),
      !      &        integerglob(7*netet+6),integerglob(ne+7*netet+6),
      !      &        integerglob(2*ne+7*netet+6),
      !      &        integerglob(nkon+2*ne+7*netet+6),
      !      &        coords(1),coords(2),coords(3),h(i),ratio,iselect,
      !      &        nselect,istartset,iendset,ialset,imastset,
      !      &        integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem)
      ! !
      !       enddo
      !
      return
      end

