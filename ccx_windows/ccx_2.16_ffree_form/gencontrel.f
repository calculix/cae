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
!
! >    \brief Calculating the normals and tangential vectors in the nodes of the slave surface (slavnor, slavtan)
! >
! >
! > @param   [in]     tieset     name and dependent surface of tie set
! > @param   [in]     ntie       number of contraints
! > @param   [in]     itietri    (1,i)pointer to node where trangulation starts for i (2,i) pointer to end
! > @param   [in]     ipkon      pointer into field kon
! > @param   [in]     kon        Field containing the connectivity of the elements in succesive order
! > @param   [in]     lakon      element label
! > @param   [in]     set        (i)name of set_i
! > @param   [in]     cg         field containing centers of gravity
! > @param   [in]     straight   (1:4 5:8 9:13,i)coeffs of plane equation for edges of triagle_i (13:16,i) coeffs of plane containing triagle
! > @param   [in]     koncont    (1:3,i) nodes of triagle_i (4,i) element face
! > @param   [in]     co         field containing the coordinates of all nodes
! > @param   [in]     vold       field containing the displacements
! > @param   [in]     nset       number of sets
! > @param   [in]     islavsurf   islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
! > @param   [in]     itiefac     pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
! > @param   [in]     islavnode   fields containing nodes of slace surfaces
! > @param   [in]     nslavnode   (i) for contraint i pointer into field islavnode
! > @param   [out]    slavnor     normal vektors in the nods of slave surface
! > @param   [out]    slavtan     tangetial vektors in the nodes of slave surface
! > @param   [in]     mi              (1) max # of integration points per element (2) max degree of freedom per element
! >
      subroutine gencontrel(tieset,ntie,itietri,ipkon,kon,&
           lakon,set,cg,straight,&
           koncont,co,vold,nset,&
           islavsurf,itiefac,&
           islavnode,nslavnode,slavnor,slavtan,mi)
      !
      !     Calculating the normals and tangential vectors in the nodes of
      !     the slave surface (slavnor, slavtan)
      !     Determining the coefficients of the dual shape functions on
      !     the slave surface
      !
      !     Author: Li, Yang; Rakotonanahary, Samoela; Sitzmann, Saskia
      !
      implicit none
      !
      logical checkbiorthogonality
      !
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
      !
      integer ntie,nset,ifree,kmax(3),ncont,&
           itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,&
           iflag,kneigh,i,j,l,islav,&
           itri,kflag,ipos,ijk,&
           index1,&
           nface,nope,m1,km1,km2,km3,number,&
           islavsurf(2,*),islavnode(*),nslavnode(ntie+1),&
           itiefac(2,*),ifaces,nelems,jfaces,mi(*),&
           mint2d,m,nopes,konl(20),id,indexnode(8),&
           line,&
           nintpoint,&
           ipiv(4),ifac,getiface,nodesf,ifs,&
           flagtan,numintpoints,numcontact,numslavint,numslavtan
      !
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),&
           et,xi,xl2s(3,8),xsj2(3),&
           shp2(7,8),t1(6),t2(6),xlnode(3),&
           xs2(3,2),slavnor(3,*),slavtan(6,*), xquad(2,8), xtri(2,6),dd,&
           al2,xn(3),xnabs(3),e(3,3),&
           err,xs2m(3,2),xsj2m(3)
      !
      include "gauss.f"
      !
      data iflag /2/
      ijk=0
      !
      !     new added data for the local coodinates for nodes
      !
      data xquad /-1, -1,&
           1, -1,&
           1, 1,&
           -1, 1,&
           0, -1,&
           1, 0,&
           0, 1,&
           -1, 0/
      !
      data xtri /0, 0,&
           1, 0,&
           0, 1,&
           0.5, 0,&
           0.5, 0.5,&
           0, 0.5/
      !
      data e /1.d0 , 0.d0 , 0.d0,&
           0.d0 , 1.d0 , 0.d0,&
           0.d0 , 0.d0 , 1.d0/
      !
      checkbiorthogonality=.false.
      flagtan=7
      !
      numcontact=40
      numslavtan=50
      numslavint=20
      numintpoints=30
      open(numcontact,file='contact.fbd',status='unknown')
      open(numslavtan,file='slavtan.fbd',status='unknown')
      open(numslavint,file='slavintmortar.out',status='unknown')
      open(numintpoints,file='intpoints.out',status='unknown')
      !
      !     maximum number of neighboring master triangles for a slave node
      !
      kflag=2
      ifree = 0     
      err=1.d-6
      do i=1,ntie  
         do l=nslavnode(i)+1,nslavnode(i+1)
            do m=1,3
               slavnor(m,l)=0.0
            enddo
         enddo
      enddo
      !
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         kneigh=1     
         slavset=tieset(2,i)
         ipos=index(slavset,' ')
         if(slavset(ipos-1:ipos-1).eq.'S') cycle
         !
         !     determining the slave set
         !
         do j=1,nset
            !      write(*,*) set(j)(1:ipos-2)
            if(set(j)(1:ipos-2).eq.slavset(1:ipos-2)) exit
         enddo
         if(j.gt.nset) then
            write(*,*) '*ERROR in gencontrel: contact slave set',&
                 slavset
            write(*,*) '       does not exist'
            call exit(201)
         endif
         islav=j
         
         do l = itiefac(1,i), itiefac(2,i)
            ifaces = islavsurf(1,l)
            nelems = int(ifaces/10)
            jfaces = ifaces - nelems*10
            call getnumberofnodes(nelems,jfaces,lakon,nope,&
                 nopes,mint2d)
            !
            !     actual position of the nodes belonging to the
            !     slave surface
            !
            do j=1,nope
               konl(j)=kon(ipkon(nelems)+j)
            enddo
            !
            do m=1,nopes
               ifac=getiface(m,jfaces,nope)
               do j=1,3
                  xl2s(j,m)=co(j,konl(ifac))+&
                       vold(j,konl(ifac))
               enddo
            enddo          
            !     calculate the normal vector in the nodes belonging to the slave surface
            !
            do m = 1, nopes
               if(nopes.eq.4 .or. nopes.eq.8)then
                  xi = xquad(1,m)
                  et = xquad(2,m)
               else
                  xi = xtri(1,m)
                  et = xtri(2,m)
               endif
               if(nopes.eq.8)then
                  call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4)then
                  call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.6)then
                  call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               else
                  call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
               endif   
               dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)&
                    + xsj2(3)*xsj2(3))
               xsj2(1) = xsj2(1)/dd
               xsj2(2) = xsj2(2)/dd
               xsj2(3) = xsj2(3)/dd                 
               ifac=getiface(m,jfaces,nope)
               node = konl(ifac)
               call nident(islavnode(nslavnode(i)+1), node,&
                    nslavnode(i+1)-nslavnode(i), id)
               index1=nslavnode(i)+id
               indexnode(m)=index1
               slavnor(1,index1) = slavnor(1,index1)&
                    +xsj2(1)
               slavnor(2,index1) = slavnor(2,index1)&
                    +xsj2(2)
               slavnor(3,index1) = slavnor(3,index1)&
                    +xsj2(3)
            enddo
 !
 105        format(4(1x,e15.8))
         enddo
         !
         !     FIRST SLAVE SURFACE LOOP DONE
         !
         !     normalizing the normals
         !
         do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
            dd=dsqrt(slavnor(1,l)**2+slavnor(2,l)**2+&
                 slavnor(3,l)**2)
            do m=1,3
               slavnor(m,l)=slavnor(m,l)/dd
            enddo
            !
            !     determining the tangential directions
            !
            do m=1,3
               xn(m)=slavnor(m,l)
               xnabs(m)=dabs(xn(m))
               xlnode(m)=co(m,node)+vold(m,node)
            enddo
            number=3
            kmax(1)=1
            kmax(2)=2
            kmax(3)=3
            kflag=2
            !
            !     sorting the components of the normal
            !
            call dsort(xnabs,kmax,number,kflag)
            
            !      elseif(flagtan==7)then
            if(1.d0 - dabs(xn(1)).lt.1.5231d-6) then       
               !
               !     calculating the local directions on master surface
               !
               slavtan(1,l)=-xn(3)*xn(1)
               slavtan(2,l)=-xn(3)*xn(2)
               slavtan(3,l)=1.d0-xn(3)*xn(3)
            else
               slavtan(1,l)=1.d0-xn(1)*xn(1)
               slavtan(2,l)=-xn(1)*xn(2)
               slavtan(3,l)=-xn(1)*xn(3)
            endif
            dd=dsqrt(slavtan(1,l)**2+ slavtan(2,l)**2&
                 +slavtan(3,l)**2)
            do m=1,3
               slavtan(m,l)=slavtan(m,l)/dd
            enddo 
            slavtan(4,l)=-(xn(2)*slavtan(3,l)-xn(3)*slavtan(2,l))
            slavtan(5,l)=-(xn(3)*slavtan(1,l)-xn(1)*slavtan(3,l))
            slavtan(6,l)=-(xn(1)*slavtan(2,l)-xn(2)*slavtan(1,l)) 
 !      endif
 !
 100        format('PNT ',i10,'P',3(1x,e15.8))
 101        format('LINE ',i10,'L',i10,'P ',i10,'P')
         enddo
      !
      enddo
      !
      close(50)
      return
      end
      
