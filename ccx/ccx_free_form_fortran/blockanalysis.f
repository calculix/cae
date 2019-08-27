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
      subroutine blockanalysis(set,nset,istartset,iendset,ialset,&
        nblk,ipkon,kon,ielfa,nodface,neiel,neij,neifa,&
        ipoface,ipnei,konf,istartblk,iendblk,nactdoh,&
        nblket,nblkze,neielsize,ielblk,nk,nactdohinv)
      !
      !     renumbering the blocks in a CFD-analysis
      !
      implicit none
      !
      logical iexternal(6)
      !
      character*81 set(*)
      !
      integer nset,istartset(*),iendset(*),ialset(*),&
        nblk,i,j,k,neighbor(8,8,8),iel,ii,indexe,ifreenei,ifour,&
        iaux,iset,j2,index,indexold,ifree,nodes(4),ipkon(*),kon(*),&
        kflag,ielfa(4,*),nodface(5,*),neiel(*),neij(*),neifa(*),&
        ipoface(*),ipnei(*),ifaceq(8,6),ifreenei2,iel2,ifa,&
        kon2(8),kon3(8),kon4(8),kon5(8),kon6(8),nactdohinv(*),&
        kon7(8),kon8(8),ielsav,ifasav,jopposite8(6),&
        konf(*),nef,indexet,indexeet,indexze,indexeze,ifaf,&
        indexe2,istartblk(*),iendblk(*),nactdoh(*),nblket(*),&
        nblkze(*),nblketl,nblkzel,neielsize,ielblk(*),nk
      !
      !     nodes belonging to the element faces
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
                  5,6,7,8,13,14,15,16,&
                  1,2,6,5,9,18,13,17,&
                  2,3,7,6,10,19,14,18,&
                  3,4,8,7,11,20,15,19,&
                  4,1,5,8,12,17,16,20/
      data jopposite8 /2,1,5,6,3,4/
      data kon2 /2,3,4,1,6,7,8,5/
      data kon3 /3,4,1,2,7,8,5,6/
      data kon4 /4,1,2,3,8,5,6,7/
      data kon5 /5,8,7,6,1,4,3,2/
      data kon6 /6,5,8,7,2,1,4,3/
      data kon7 /7,6,5,8,3,2,1,4/
      data kon8 /8,7,6,5,4,3,2,1/
      !
      nblk=0
      do iset=1,nset
         if(set(iset)(1:10).eq.'FLUIDBLOCK') then
            nblk=nblk+1
         endif
      enddo
      !
      if(nblk.eq.0) return
      !
      !     initialize neighbor
      !
      do i=1,8
         do j=1,8
            do k=1,8
               neighbor(i,j,k)=0
            enddo
         enddo
      enddo
      neighbor(5,1,2)=4
      neighbor(2,1,5)=4
      neighbor(1,2,6)=3
      neighbor(6,2,1)=3
      neighbor(2,6,5)=7
      neighbor(5,6,2)=7
      neighbor(6,5,1)=8
      neighbor(1,5,6)=8
      neighbor(2,1,4)=5
      neighbor(4,1,2)=5
      neighbor(1,4,3)=8
      neighbor(3,4,1)=8
      neighbor(4,3,2)=7
      neighbor(2,3,4)=7
      neighbor(3,2,1)=6
      neighbor(1,2,3)=6
      neighbor(5,6,7)=2
      neighbor(7,6,5)=2
      neighbor(6,7,8)=3
      neighbor(8,7,6)=3
      neighbor(7,8,5)=4
      neighbor(5,8,7)=4
      neighbor(8,5,6)=1
      neighbor(6,5,8)=1
      neighbor(2,3,7)=4
      neighbor(7,3,2)=4
      neighbor(3,7,6)=8
      neighbor(6,7,3)=8
      neighbor(7,6,2)=5
      neighbor(2,6,7)=5
      neighbor(6,2,3)=1
      neighbor(3,2,6)=1
      neighbor(8,7,3)=6
      neighbor(3,7,8)=6
      neighbor(7,3,4)=2
      neighbor(4,3,7)=2
      neighbor(3,4,8)=1
      neighbor(8,4,3)=1
      neighbor(4,8,7)=5
      neighbor(7,8,4)=5
      neighbor(5,8,4)=7
      neighbor(4,8,5)=7
      neighbor(8,4,1)=3
      neighbor(1,4,8)=3
      neighbor(4,1,5)=2
      neighbor(5,1,4)=2
      neighbor(1,5,8)=6
      neighbor(8,5,1)=6
      !
      ifour=4
      kflag=1
      nef=0
      !
      !     loop over all blocks
      !
      nblk=0
      do iset=1,nset
         if(set(iset)(1:10).ne.'FLUIDBLOCK') cycle
         nblk=nblk+1
         iel=ialset(istartset(iset))
         !
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
         ifreenei=0
         do i=1,nk
            ipoface(i)=0
         enddo
         do i=1,6*neielsize
            neiel(i)=0
         enddo
         !
         do ii=istartset(iset),iendset(iset)
            i=ialset(ii)
            indexe=ipkon(i)
            !
            !           hex element
            !
            ipnei(i)=ifreenei
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
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifreenei=ifreenei+1
                     neiel(ifreenei)=0
                     neifa(ifreenei)=ifree
                     ielfa(1,ifree)=i
                     ielfa(4,ifree)=j
                     ifree=ifree+1
                     exit
                  endif
                  !
                  !                 removing a surface which has already
                  !                 been catalogued
                  !
                  if((nodface(1,index).eq.nodes(2)).and.&
                     (nodface(2,index).eq.nodes(3))) then
                     ifreenei=ifreenei+1
                     !
                     !                    completing the facial info in neifa
                     !
                     neifa(ifreenei)=index
                     ielfa(2,index)=i
                     !
                     !                    the neighboring elements to the face are i and iel2
                     !                    with local face number for the face of j and j2
                     !
                     iel2=ielfa(1,index)
                     j2=ielfa(4,index)
                     !
                     !                    completing the neighboring info for (element i,side j)
                     !
                     neiel(ifreenei)=iel2
                     neij(ifreenei)=j2
                     !
                     !                    completing the neighboring info for (element iel2,side j2)
                     !
                     ifreenei2=ipnei(iel2)+j2
                     neiel(ifreenei2)=i
                     neij(ifreenei2)=j
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         enddo
         !
         !        looking for a corner element of the block
         !
         index=ipnei(iel)
         ifa=1
         !
         !        marching in one direction up to the wall
         !
         do
            if(neiel(index+ifa).eq.0) then
               if((ifa.eq.1).or.(ifa.eq.2)) then
                  ifa=3
               elseif((ifa.eq.3).or.(ifa.eq.5)) then
                  ifa=4
               else
                  ifa=1
               endif
               exit
            endif
            iel=neiel(index+ifa)
            ifa=neij(index+ifa)
            ifa=jopposite8(ifa)
            index=ipnei(iel)
         enddo
         !
         !        marching in the next direction up to the wall
         !
         do
            if(neiel(index+ifa).eq.0) then
               if((ifa.eq.1).or.(ifa.eq.2)) then
                  ifa=3
               elseif((ifa.eq.3).or.(ifa.eq.5)) then
                  ifa=4
               else
                  ifa=1
               endif
               exit
            endif
            iel=neiel(index+ifa)
            ifa=neij(index+ifa)
            ifa=jopposite8(ifa)
            index=ipnei(iel)
         enddo
         !
         !        marching in the third direction up to the wall
         !
         ielsav=iel
         ifasav=ifa
         !
         do
            if(neiel(index+ifa).eq.0) then
               exit
            endif
            iel=neiel(index+ifa)
            ifa=neij(index+ifa)
            ifa=jopposite8(ifa)
            index=ipnei(iel)
         enddo
         !
         !        check whether the element is a corner element and
         !        modifying the topology in such
         !        a way that the local axes are pointing inwards into the
         !        block
         !
         index=ipnei(iel)
         do j=1,6
            if(neiel(index+j).eq.0) then
               iexternal(j)=.true.
            else
               iexternal(j)=.false.
            endif
         enddo
         indexe=ipkon(iel)
         if(iexternal(1).and.iexternal(6).and.iexternal(3)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+j)
            enddo
         elseif(iexternal(1).and.iexternal(3).and.iexternal(4)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon2(j))
            enddo
         elseif(iexternal(1).and.iexternal(4).and.iexternal(5)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon3(j))
            enddo
         elseif(iexternal(1).and.iexternal(5).and.iexternal(6)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon4(j))
            enddo
         elseif(iexternal(2).and.iexternal(3).and.iexternal(6)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon5(j))
            enddo
         elseif(iexternal(2).and.iexternal(3).and.iexternal(6)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon5(j))
            enddo
         elseif(iexternal(2).and.iexternal(3).and.iexternal(4)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon6(j))
            enddo
         elseif(iexternal(2).and.iexternal(4).and.iexternal(5)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon7(j))
            enddo
         elseif(iexternal(2).and.iexternal(5).and.iexternal(6)) then
            do j=1,8
               konf(indexe+j)=kon(indexe+kon8(j))
            enddo
         else
            !
            !        if no corner: look in the other direction
            !
            if((ifasav.eq.1).or.(ifasav.eq.2)) then
               ifa=4
            elseif((ifasav.eq.3).or.(ifasav.eq.5)) then
               ifa=1
            else
               ifa=3
            endif
            iel=ielsav
            !
            do
               if(neiel(index+ifa).eq.0) then
                  exit
               endif
               iel=neiel(index+ifa)
               ifa=neij(index+ifa)
               ifa=jopposite8(ifa)
               index=ipnei(iel)
            enddo
            !
            !           iel is now a corner element
            !
            !           modifying the topology of the corner element in such
            !           a way that the local axes are pointing inwards into the
            !           block
            !
            index=ipnei(iel)
            do j=1,6
               if(neiel(index+j).eq.0) then
                  iexternal(j)=.true.
               else
                  iexternal(j)=.false.
               endif
            enddo
            indexe=ipkon(iel)
            if(iexternal(1).and.iexternal(6).and.iexternal(3)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+j)
               enddo
            elseif(iexternal(1).and.iexternal(3).and.iexternal(4)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon2(j))
               enddo
            elseif(iexternal(1).and.iexternal(4).and.iexternal(5)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon3(j))
               enddo
            elseif(iexternal(1).and.iexternal(5).and.iexternal(6)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon4(j))
               enddo
            elseif(iexternal(2).and.iexternal(3).and.iexternal(6)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon5(j))
               enddo
            elseif(iexternal(2).and.iexternal(3).and.iexternal(6)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon5(j))
               enddo
            elseif(iexternal(2).and.iexternal(3).and.iexternal(4)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon6(j))
               enddo
            elseif(iexternal(2).and.iexternal(4).and.iexternal(5)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon7(j))
               enddo
            elseif(iexternal(2).and.iexternal(5).and.iexternal(6)) then
               do j=1,8
                  konf(indexe+j)=kon(indexe+kon8(j))
               enddo
            else
               write(*,*) '*ERROR in blockanalysis'
               write(*,*) '       no block structure'
               write(*,*) '       no corner element found'
               call exit(201)
            endif
         endif
         !
         istartblk(nblk)=nef+1
         index=ipnei(iel)
         nblketl=0
         nblkzel=0
         !
         !        renumbering the elements
         !
         do
            indexze=index
            indexeze=indexe
            !
            !           check whether all layers have the same number of
            !           elements
            !
            if(nblkzel.ne.0) then
               if(nblkze(nblk).eq.0) then
                  nblkze(nblk)=nblkzel
               elseif(nblkze(nblk).ne.nblkzel) then
                  write(*,*) '*ERROR in blockanalysis'
                  write(*,*) '       block ',set(iset)
                  write(*,*) '       is no block'
                  call exit(201)
               endif
               nblkzel=0
            endif
            !
            do
               indexet=index
               indexeet=indexe
               !
               if(nblketl.ne.0) then
                  if(nblket(nblk).eq.0) then
                     nblket(nblk)=nblketl
                  elseif(nblket(nblk).ne.nblketl) then
                     write(*,*) '*ERROR in blockanalysis'
                     write(*,*) '       block ',set(iset)
                     write(*,*) '       is no block'
                     call exit(201)
                  endif
                  nblketl=0
               endif
               !
               nblketl=nblketl+1
               nblkzel=nblkzel+1
               !
               nef=nef+1
               nactdoh(iel)=nef
               nactdohinv(nef)=iel
               ielblk(nef)=nblk
               ifaf=4
               call identifyface(konf(indexe+1),kon(indexe+1),ifaf,ifa)
               !
               do
                  iel=neiel(index+ifa)
                  if(iel.eq.0) exit
                  ifa=neij(index+ifa)
                  !
                  nblketl=nblketl+1
                  nblkzel=nblkzel+1
                  !
                  nef=nef+1
                  nactdoh(iel)=nef
                  nactdohinv(nef)=iel
                  ielblk(nef)=nblk
                  ifa=jopposite8(ifa)
                  !
                  indexe2=ipkon(iel)
                  call adaptconnectivity(konf(indexe+1),ifaf,&
                       kon(indexe2+1),konf(indexe2+1),neighbor)
                  indexe=indexe2
                  index=ipnei(iel)
               enddo
               !
               ifaf=5
               call identifyface(konf(indexeet+1),kon(indexeet+1),&
                                 ifaf,ifa)
               iel=neiel(indexet+ifa)
               if(iel.eq.0) exit
               !
               indexe2=ipkon(iel)
               call adaptconnectivity(konf(indexeet+1),ifaf,&
                    kon(indexe2+1),konf(indexe2+1),neighbor)
               indexe=indexe2
               index=ipnei(iel)
            enddo
            !
            ifaf=2
            call identifyface(konf(indexeze+1),kon(indexeze+1),&
                              ifaf,ifa)
            iel=neiel(indexze+ifa)
            if(iel.eq.0) exit
            !
            indexe2=ipkon(iel)
            call adaptconnectivity(konf(indexeze+1),ifaf,&
                 kon(indexe2+1),konf(indexe2+1),neighbor)
            indexe=indexe2
            index=ipnei(iel)
         enddo
         iendblk(nblk)=nef
      !
      enddo
      !
      !       open(91,file='kon',status='unknown')
      !       do i=1,nef
      !          iel=nactdohinv(i)
      !          indexe=ipkon(iel)
      !          write(91,100) i
      !  100     format(' -1',i10,'    1    0WATER')
      !          write(91,101) (konf(indexe+j),j=1,8)
      !  101     format(' -2',8i10)
      !       enddo
      !       close(91)
      !
      return
      end
