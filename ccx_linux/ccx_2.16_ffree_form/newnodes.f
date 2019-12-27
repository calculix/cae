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
      subroutine newnodes(nktet_,ipoed,n,iedg,h,d,r,&
        conewnodes,cotet,ibasenewnodes,ipoeled,ieled,doubleglob,&
        integerglob,nnewnodes,iedgnewnodes,hnewnodes,iexternedg)
      !
      !     determines the edges to be divided and the number of sub-intervals
      !     for each such edge based on the h-field
      !
      implicit none
      !
      integer indexn,i,nktet_,ipoed(*),n(*),iedg(3,*),&
        n1,n2,k,ibasenewnodes(*),ipoeled(*),index,j,ieled(2,*),&
        nnewnodes,nx(nnewnodes),ny(nnewnodes),nz(nnewnodes),kflag,&
        neigh(2),integerglob(*),nktri,netet,ne,nkon,nfaces,nfield,&
        nselect,iselect(1),imastset,istartset(1),iendset(1),ialset(1),&
        nterms,konl(20),kneigh,iedgnewnodes(*),nelem,loopa,&
        iexternedg(*)
      !
      real*8 alphasum,h(*),h1,d(*),r(*),scale,hscale,rscale,&
        conewnodes(3,*),cotet(3,*),x(nnewnodes),y(nnewnodes),&
        z(nnewnodes),xo(nnewnodes),yo(nnewnodes),zo(nnewnodes),&
        doubleglob(*),coords(3),ratio(20),dnew,hnewnodes(*)
      !
      indexn=0
      !
      loop: do i=1,nktet_
         index=ipoed(i)
         do
            if(index.eq.0) cycle loop
            !
            !           no subdivision or external edge: cycle
            !
            if(n(index).eq.0) then
               index=iedg(3,index)
               cycle
            endif
            !
            !           scale subdivision lengths:
            !
            alphasum=0.d0
            n1=iedg(1,index)
            n2=iedg(2,index)
            h1=h(n1)
            !
            !           the subintervals are numbered 0,...,n(index)
            !           the length of the subintervals is progressively
            !           increased by r(index) (may be negative)
            !
            do k=0,n(index)
               alphasum=alphasum+h1+(k+1)*r(index)
            enddo
            scale=d(index)/alphasum
            !
            !           loop over all new positions: numbered from
            !           1,...,n(index)
            !
            alphasum=0.d0
            hscale=scale*h1
            rscale=scale*r(index)
            do j=1,n(index)
               alphasum=alphasum+hscale+j*rscale
               do k=1,3
                  conewnodes(k,indexn+j)=cotet(k,n1)+alphasum*&
                       (cotet(k,n2)-cotet(k,n1))/d(index)
               enddo
               !                hnewnodes(indexn+j)=(h(n1)+h(n2))/2.d0
               !                hnewnodes(indexn+j)=h(n1)+alphasum*
               !      &              (h(n2)-h(n1))/d(index)
               !                write(*,*) 'newnodes indexn j ',indexn,j
               !                hhhh(indexn+j)=h(n1)+alphasum*
               !                hhhh=h(n1)+alphasum*
               !      &              (h(n2)-h(n1))/d(index)
               !
               !              determine an element to which the new node belongs
               !
               ibasenewnodes(indexn+j)=ieled(1,ipoeled(index))
               iedgnewnodes(indexn+j)=index
            enddo
            indexn=indexn+n(index)
            index=iedg(3,index)
         enddo
      enddo loop
      !
      !     initializing fields
      !
      nktri=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=1
      nselect=1
      iselect(1)=1
      imastset=0
      !
      !     interpolating the desired element size h for the new
      !     nodes
      !
      !       kflag=2
      !       kneigh=2
      !       loopa=2
      !
      !     preparing a field with the new nodes for routine near3d
      !
      !       do i=1,nnewnodes
      !          x(i)=conewnodes(1,i)
      !          y(i)=conewnodes(2,i)
      !          z(i)=conewnodes(3,i)
      !          xo(i)=x(i)
      !          yo(i)=y(i)
      !          zo(i)=z(i)
      !          nx(i)=i
      !          ny(i)=i
      !          nz(i)=i
      !       enddo
      ! !
      !       if(nnewnodes.gt.0) then
      !          call dsort(x,nx,nnewnodes,kflag)
      !          call dsort(y,ny,nnewnodes,kflag)
      !          call dsort(z,nz,nnewnodes,kflag)
      !       endif
      !
      !     check for each new nodes the distance to the closest new
      !     node
      !
      do i=1,nnewnodes
         !          if(ibasenewnodes(i).eq.0) cycle
         !
         !        looking for the closest new node
         !
         !          call near3d(xo,yo,zo,x,y,z,nx,ny,nz,conewnodes(1,i),
         !      &                  conewnodes(2,i),conewnodes(3,i),nnewnodes,
         !      &                  neigh,kneigh)
         !          if(ibasenewnodes(neigh(2)).eq.0) cycle
         !
         !          dnew=dsqrt((conewnodes(1,i)-xo(neigh(2)))**2+
         !      &              (conewnodes(2,i)-yo(neigh(2)))**2+
         !      &              (conewnodes(3,i)-zo(neigh(2)))**2)
         !
         !        perform the interpolation of h for the internal node
         !
         !          coords(1)=xo(i)
         !          coords(2)=yo(i)
         !          coords(3)=zo(i)
         coords(1)=conewnodes(1,i)
         coords(2)=conewnodes(2,i)
         coords(3)=conewnodes(3,i)
            !          if(iexternedg(iedgnewnodes(i)).ne.0) then
            loopa=8
         !          else
         !             loopa=2
         !          endif
         !          write(*,*) 'newnodes ',i,iexternedg(iedgnewnodes(i))
         !          hhhh=hnewnodes(i)
         call basis(doubleglob(1),doubleglob(netet+1),&
                 doubleglob(2*netet+1),&
                 doubleglob(3*netet+1),doubleglob(4*netet+1),&
                 doubleglob(5*netet+1),integerglob(6),&
                 integerglob(netet+6),&
                 integerglob(2*netet+6),doubleglob(6*netet+1),&
                 integerglob(3*netet+6),nktri,netet,&
                 doubleglob(4*nfaces+6*netet+1),nfield,&
                 doubleglob(nktri+4*nfaces+6*netet+1),&
                 integerglob(7*netet+6),integerglob(ne+7*netet+6),&
                 integerglob(2*ne+7*netet+6),&
                 integerglob(nkon+2*ne+7*netet+6),&
                 coords(1),coords(2),coords(3),hnewnodes(i),&
                 ratio,iselect,nselect,&
                 istartset,iendset,ialset,imastset,&
                 integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,&
                 loopa)
      !          write(*,*) 'newnodes',hnewnodes(i),
      !      &           dabs((hhhh-hnewnodes(i))/hnewnodes(i))*100.d0
      !          hnewnodes(i)=hhhh(i)
      !
      !             if(hnewnodes(i).gt.dnew) ibasenewnodes(neigh(2))=0
      enddo
      ! !
      ! !        removing the deleted new nodes
      ! !
      !          m=0
      !          do i=1,nnewnodes
      !             if(ibasenewnodes(i).eq.0) cycle
      !             m=m+1
      !             do k=1,3
      !                conewnodes(k,m)=conewnodes(k,i)
      !             enddo
      !             ibasenewnodes(m)=ibasenewnodes(i)
      !             iedgnewnodes(m)=iedgnewnodes(i)
      !             hnewnodes(m)=hnewnodes(i)
      !          enddo
      !          nnewnodes=m
      !
      return
      end
