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
      subroutine midexternalfaces(iexternfa,nexternfa,ifacext,&
         ifreefa,itetfa,ifac,kontet,kontetor,ialsetexternel,&
         nexternel,iedgextfa,ifacexted,ipoed,iedg,iexternedg)
      !
      !     stores the nodes belonging to the external faces of the
      !     unrefined mesh
      !
      implicit none
      !
      integer nexternfa,iexternfa(*),ifacext(6,*),i,j,ifreefa,&
        iel,itetfa(2,*),ifac(4,*),kontetor(6,*),nofa(6,4),&
        kontet(4,*),node,k,ialsetexternel(*),nexternel,id,&
        iedgextfa(2,*),ifacexted(3,*),n1,n2,n3,indexe,iedg(3,*),&
        iexternedg(*),ipoed(*),iedge
      !
      intent(in) ifreefa,itetfa,ifac,kontet,kontetor,ipoed,iedg,&
        iexternedg
      !
      intent(inout) nexternfa,iexternfa,ifacext,ialsetexternel,&
        nexternel,iedgextfa,ifacexted
      !
      !     local node numbers of the nodes belonging to the 4 faces
      !     of a tetrahedron. From the middle nodes 4 is subtracted,
      !     since they are stored in position 1..6 in field kontetor
      !
      data nofa /2,3,4,2,6,5,&
                 3,1,4,3,4,6,&
                 1,2,4,1,5,4,&
                 1,3,2,3,2,1/
      !
      nexternfa=0
      nexternel=0
      !
      !     loop over all faces. For the unrefined mesh the faces are stored
      !     in ifac in a consecutive order.
      !     A face i is external if iexternfa(i)!=0; at exit the value
      !     of iexternfa(i) is the number of the external face in field
      !     ifacext
      !
      do i=1,ifreefa-1
         if(iexternfa(i).ne.0) then
            nexternfa=nexternfa+1
            iexternfa(i)=nexternfa
            !
            !           recovering the element containing the face
            !
            iel=itetfa(1,i)
            !
            !           the local number of the face is the local number of the
            !           vertex node not belonging to the face
            !
            loop: do j=1,4
               node=kontet(j,iel)
               do k=1,3
                  if(ifac(k,i).eq.node) cycle loop
               enddo
               exit
            enddo loop
            !
            do k=1,3
               ifacext(k,nexternfa)=kontet(nofa(k,j),iel)
               ifacext(k+3,nexternfa)=kontetor(nofa(k+3,j),iel)
            enddo
            !
            !           catalogueing the elements adjacent to an external face =
            !           external elements
            !
            !           these elements are needed for the projection on the
            !           external surface for newly generated nodes in the
            !           refined mesh
            !
            !             call nident(ialsetexternel,iel,nexternel,id)
            !             if(id.gt.0) then
            !                if(ialsetexternel(id).eq.iel) cycle
            !             endif
            !             nexternel=nexternel+1
            !             do j=nexternel,id+2,-1
            !                ialsetexternel(j)=ialsetexternel(j-1)
            !             enddo
            !             ialsetexternel(id+1)=iel
            !
            !           catalogueing the (external) edges per external face in
            !           ifacexted and the external faces per external edge in
            !           iedgextfa
            !
            !           first (external) edge of external face
            !
            n1=ifacext(1,nexternfa)
            n2=ifacext(2,nexternfa)
            if(n2.lt.n1) then
               n3=n1
               n1=n2
               n2=n3
            endif
            !             write(*,*) 'midexternafaces1 ',n1,n2,nexternfa
            indexe=ipoed(n1)
            do
               if(iedg(2,indexe).eq.n2) then
                  iedge=iexternedg(indexe)
                  ifacexted(1,nexternfa)=iedge
                  if(iedgextfa(1,iedge).eq.0) then
                     iedgextfa(1,iedge)=nexternfa
                  else
                     iedgextfa(2,iedge)=nexternfa
                  endif
                  !                   write(*,*) iedge,iedgextfa(1,iedge),iedgextfa(2,iedge)
                  exit
               else
                  indexe=iedg(3,indexe)
               endif
            enddo
            !
            !           second (external) edge of external face
            !
            n1=ifacext(2,nexternfa)
            n2=ifacext(3,nexternfa)
            if(n2.lt.n1) then
               n3=n1
               n1=n2
               n2=n3
            endif
            !             write(*,*) 'midexternafaces2 ',n1,n2,nexternfa
            indexe=ipoed(n1)
            do
               if(iedg(2,indexe).eq.n2) then
                  iedge=iexternedg(indexe)
                  ifacexted(2,nexternfa)=iedge
                  if(iedgextfa(1,iedge).eq.0) then
                     iedgextfa(1,iedge)=nexternfa
                  else
                     iedgextfa(2,iedge)=nexternfa
                  endif
                  !                   write(*,*) iedge,iedgextfa(1,iedge),iedgextfa(2,iedge)
                  exit
               else
                  indexe=iedg(3,indexe)
               endif
            enddo
            !
            !           third (external) edge of external face
            !
            n1=ifacext(3,nexternfa)
            n2=ifacext(1,nexternfa)
            if(n2.lt.n1) then
               n3=n1
               n1=n2
               n2=n3
            endif
            !             write(*,*) 'midexternafaces3 ',n1,n2,nexternfa
            indexe=ipoed(n1)
            do
               if(iedg(2,indexe).eq.n2) then
                  iedge=iexternedg(indexe)
                  ifacexted(3,nexternfa)=iedge
                  if(iedgextfa(1,iedge).eq.0) then
                     iedgextfa(1,iedge)=nexternfa
                  else
                     iedgextfa(2,iedge)=nexternfa
                  endif
                  !                   write(*,*) iedge,iedgextfa(1,iedge),iedgextfa(2,iedge)
                  exit
               else
                  indexe=iedg(3,indexe)
               endif
            enddo
         !
         endif
      enddo
      !
      !       do i=1,54
      !          write(*,*) 'iedgext/fa',
      !      &            (iedgextfa(j,i),j=1,2)
      !       enddo
      !       do i=1,36
      !          write(*,*) 'ifacext/ed',(ifacext(j,i),j=1,6),
      !      &              (ifacexted(j,i),j=1,3)
      !       enddo
      !       do i=1,69
      !          write(*,*) 'iedg',(iedg(j,i),j=1,2),iexternedg(i)
      !       enddo
      !       do i=1,80
      !          write(*,*) 'ifac',(ifac(j,i),j=1,3),iexternfa(i)
      !       enddo
      !
      return
      end
