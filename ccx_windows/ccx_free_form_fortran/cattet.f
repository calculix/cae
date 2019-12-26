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
      subroutine cattet(kontet,netet_,ifac,ne,ipkon,kon,ifatet,ifreetet,&
        bc,itetfa,ifreefa,planfa,ipofa,cotet,cg,ipoeln,ieln,ifreeln,&
        lakon,kontetor,iquad)
      !
      !     catalogueing the tetrahedral elements of the mesh
      !
      implicit none
      !
      character*8 lakon(*)
      !
      integer kontet(4,*),netet_,i,ifac(4,*),ne,ipkon(*),kon(*),index,j,&
        nodes(4),ifatet(4,*),ifreetet,itetfa(2,*),ifreefa,ipofa(*),&
        ipoeln(*),ieln(2,*),node,ifreeln,kontetor(6,*),iquad
      !
      real*8 bc(4,*),planfa(4,*),cotet(3,*),cg(3,*)
      !
      !     determine the first tetrahedral element or first
      !     unused element
      !
      do i=1,ne
         if((lakon(i)(1:4).eq.'C3D4').or.&
            (lakon(i)(1:5).eq.'C3D10').or.&
            (lakon(i)(1:1).eq.char(0))) exit
      enddo
      ifreetet=i
      !
      do 
         do j=i+1,netet_
            !
            !           element number supersedes largest one
            !
            if(j.gt.ne) then
               kontet(4,i)=j
               i=j
               exit
            endif
            !
            !           tetrahedral element
            !
            if((lakon(j)(1:4).eq.'C3D4').or.&
               (lakon(j)(1:5).eq.'C3D10')) then
               kontet(4,i)=j
               i=j
               exit
            endif
            !
            !           unused element
            !
            if(lakon(j)(1:1).eq.char(0)) then
               kontet(4,i)=j
               i=j
               exit
            endif
         enddo
         if(j.eq.netet_) exit
      enddo
      kontet(4,netet_)=0
      !
      !     initialization of ipofa and ifac
      !
      do i=1,4*netet_
         ifac(4,i)=i+1
      enddo
      ifac(4,4*netet_)=0
      !
      !     initialization of ieln
      !
      do i=1,4*netet_
         ieln(2,i)=i+1
      enddo
      ieln(2,4*netet_)=0
      !
      !     adding the tetrahedral elements one by one
      !
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if((lakon(i)(1:4).ne.'C3D4').and.&
            (lakon(i)(1:5).ne.'C3D10')) cycle
         index=ipkon(i)
         do j=1,4
            nodes(j)=kon(index+j)
         enddo
         !
         !        if C3D10: store the middle nodes
         !        If there is at least one C3D10 element iquad is set to 1
         !        which means that the refined mesh will be fully quadratic
         !
         if(lakon(i)(4:4).eq.'1') then
            iquad=1
            do j=1,6
               kontetor(j,ifreetet)=kon(index+4+j)
            enddo
         endif
         !
         call generatetet(kontet,ifatet,ifreetet,bc,ifac,itetfa,&
                          ifreefa,planfa,ipofa,nodes,cotet,cg)
      enddo
      !
      !     generating the element per node relationship
      !
      do j=1,netet_
         if(kontet(1,j).eq.0) cycle
         do i=1,4
            node=kontet(i,j)
            index=ifreeln
            ieln(1,index)=j
            ifreeln=ieln(2,index)
            if(ifreeln.eq.0) then
               write(*,*) '*ERROR in cattet: increase the'
               write(*,*) '       dimension of ieln'
               call exit(201)
            endif
            ieln(2,index)=ipoeln(node)
            ipoeln(node)=index
         enddo
      enddo
      !
      return
      end

