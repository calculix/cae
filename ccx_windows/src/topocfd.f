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
      subroutine topocfd(ne,ipkon,kon,lakon,ipnei,neifa,neiel,ipoface,
     &     nodface,ielfa,nflnei,nface,ifaext,nfaext,
     &     isolidsurf,nsolidsurf,set,nset,istartset,iendset,ialset,
     &     vel,vold,mi,neij,nef,nactdoh,ipkonf,lakonf,ielmatf,ielmat,
     &     ielorienf,ielorien,norien,cs,mcs,tieset,x,y,z,xo,yo,zo,nx,
     &     ny,nz,co,ifatie,velo,veloo,initial)
!     
!     gathering topological information for computational fluid 
!     dynamics calculations
!     
      implicit none
!     
      character*8 lakon(*),lakonf(*)
      character*81 set(*),noset,tieset(3,*),slavset,mastset
!     
      integer ne,ipkon(*),ipnei(*),ipoface(*),nodface(5,*),neifa(*),
     &     ielfa(4,*),nflnei,nface,i,j,k,index,indexe,neiel(*),ithree,
     &     nfaext,ifaext(*),isolidsurf(*),nsolidsurf,indexf,nneigh,
     &     nset,istartset(*),iendset(*),ialset(*),iaux,kflag,ifour,iel2,
     &     ifaceq(8,6),ifacet(7,4),ifacew(8,5),kon(*),nodes(4),iel1,j2,
     &     indexold,ifree,ifreenew,ifreenei,mi(*),neij(*),ifreenei2,nef,
     &     nactdoh(*),ipkonf(*),ielmatf(mi(3),*),ielmat(mi(3),*),nf(5),
     &     nope,ielorien(mi(3),*),ielorienf(mi(3),*),norien,
     &     jopposite6(5),itie,nx(*),ny(*),nz(*),noden(1),nelemm,nelems,
     &     n,mcs,l,jfacem,jfaces,islav,imast,ifaces,ifacem,ifatie(*),
     &     nodeinface,nodeoutface,nopes,jop,initial,jopposite8(6)
!     
      real*8 vel(nef,0:7),vold(0:mi(2),*),coords(3),cs(17,*),x(*),
     &     y(*),z(*),xo(*),yo(*),zo(*),co(3,*),a(3),b(3),xn(3),p(3),
     &     q(3),c(3,3),dot,dc,ds,dd,theta,pi,velo(nef,0:7),
     &     veloo(nef,0:7)
!     
!     nodes belonging to the element faces
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &     1,2,4,5,9,8,12,
     &     2,3,4,6,10,9,13,
     &     1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data jopposite6 /2,1,0,0,0/
      data jopposite8 /2,1,5,6,3,4/
      data nf /3,3,4,4,4/
!     
!     transfer element information from ipkon, lakon, ielmat
!     and ielorien to ipkonf, lakonf, ielmatf and ielorienf. 
!     The latter fields were created
!     for CFD applications in which the elements are numerated in
!     ascending order (without gaps)
!     
      do i=1,ne
c     write(*,*) 'topocfd ',i,nactdoh(i)
        if(nactdoh(i).ne.0) then
          ipkonf(nactdoh(i))=ipkon(i)
          lakonf(nactdoh(i))=lakon(i)
          do j=1,mi(3)
            ielmatf(j,nactdoh(i))=ielmat(j,i)
            if(norien.gt.0) ielorienf(j,nactdoh(i))=ielorien(j,i)
          enddo
        endif
      enddo
!     
      kflag=1
      ithree=3
      ifour=4
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
!     
      do i=1,6*nef
        nodface(5,i)=i+1
      enddo
      do i=1,nef
        indexe=ipkonf(i)
        if(lakonf(i)(4:4).eq.'8') then
!     
!     hex element
!     
          ipnei(i)=ifreenei
          do j=1,6
            do k=1,4
              nodes(k)=kon(indexe+ifaceq(k,j))
            enddo
            call insertsorti(nodes,ifour)
c     call isortii(nodes,iaux,ifour,kflag)
c     write(*,*) 'topocfd ',i,(nodes(k),k=1,4)
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifreenei=ifreenei+1
                neifa(ifreenei)=ifree
                ielfa(1,ifree)=i
                ielfa(4,ifree)=j
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                ifreenei=ifreenei+1
!     
!     completing the facial info in neifa
!     
                neifa(ifreenei)=index
                ielfa(2,index)=i
!     
!     the neighboring elements to the face are i and iel2
!     with local face number for the face of j and j2
!     
                iel2=ielfa(1,index)
                j2=ielfa(4,index)
!     
!     completing the neighboring info for (element i,side j)
!     
                neiel(ifreenei)=iel2
                neij(ifreenei)=j2
!     
!     completing the neighboring info for (element iel2,side j2)
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
        else if(lakonf(i)(4:4).eq.'6') then
!     
!     wedge element
!     
          ipnei(i)=ifreenei
          do j=1,5
            do k=1,nf(j)
              nodes(k)=kon(indexe+ifacew(k,j))
            enddo
            call insertsorti(nodes,nf(j))
c     call isortii(nodes,iaux,nf(j),kflag)
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifreenei=ifreenei+1
                neifa(ifreenei)=ifree
                ielfa(1,ifree)=i
                ielfa(4,ifree)=j
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                ifreenei=ifreenei+1
!     
!     completing the facial info in neifa
!     
                neifa(ifreenei)=index
                ielfa(2,index)=i
!     
!     the neighboring elements to the face are i and iel2
!     with local face number for the face of j and j2
!     
                iel2=ielfa(1,index)
                j2=ielfa(4,index)
!     
!     completing the neighboring info for (element i,side j)
!     
                neiel(ifreenei)=iel2
                neij(ifreenei)=j2
!     
!     completing the neighboring info for (element iel2,side j2)
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
        else
!     
!     tet element
!     
          ipnei(i)=ifreenei
          do j=1,4
            do k=1,3
              nodes(k)=kon(indexe+ifacet(k,j))
            enddo
            call insertsorti(nodes,ithree)
c     call isortii(nodes,iaux,ithree,kflag)
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifreenei=ifreenei+1
                neifa(ifreenei)=ifree
                ielfa(1,ifree)=i
                ielfa(4,ifree)=j
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                ifreenei=ifreenei+1
!     
!     completing the facial info in neifa
!     
                neifa(ifreenei)=index
                ielfa(2,index)=i
!     
!     the neighboring elements to the face are i and iel2
!     with local face number for the face of j and j2
!     
                iel2=ielfa(1,index)
                j2=ielfa(4,index)
!     
!     completing the neighboring info for (element i,side j)
!     
                neiel(ifreenei)=iel2
                neij(ifreenei)=j2
!     
!     completing the neighboring info for (element iel2,side j2)
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
        endif
      enddo
!     
      nflnei=ifreenei
      ipnei(nef+1)=ifreenei
      nface=ifree-1
!     
!     check for cyclic symmetry (translational or rotational)
!     
      do i=1,mcs
        itie=int(cs(17,i))
        if((tieset(1,itie)(81:81).ne.'P').and.
     &       (tieset(1,itie)(81:81).ne.'Z')) cycle
!     
        slavset=tieset(2,itie)
        mastset=tieset(3,itie)
!     
        do j=1,nset
          if(set(j).eq.slavset) exit
        enddo
        islav=j
!     
        do j=1,nset
          if(set(j).eq.mastset) exit
        enddo
        imast=j
!     
!     catalogueing the center of gravity of the master faces
!     
        n=0
        do j=istartset(imast),iendset(imast)
          ifacem=ialset(j)
!     
          nelemm=int(ifacem/10)
          jfacem=ifacem-nelemm*10
!     
!     taking the renumbering of the elements into account
!     
          nelemm=nactdoh(nelemm)
!     
          indexe=ipkonf(nelemm)
!     
!     determination of the face center of the master face
!     
          if(lakonf(nelemm)(4:4).eq.'8') then
            do l=1,3
              coords(l)=0.d0
            enddo
            do k=1,4
              nodes(k)=kon(indexe+ifaceq(k,jfacem))
              do l=1,3
                coords(l)=coords(l)+co(l,nodes(k))
              enddo
            enddo
            do l=1,3
              coords(l)=coords(l)/4.d0
            enddo
          elseif(lakonf(nelemm)(4:4).eq.'6') then
            do l=1,3
              coords(l)=0.d0
            enddo
            do k=1,nf(jfacem)
              nodes(k)=kon(indexe+ifacew(k,jfacem))
              do l=1,3
                coords(l)=coords(l)+co(l,nodes(k))
              enddo
            enddo
            do l=1,3
              coords(l)=coords(l)/nf(jfacem)
            enddo
          else
            do l=1,3
              coords(l)=0.d0
            enddo
            do k=1,3
              nodes(k)=kon(indexe+ifaceq(k,jfacem))
              do l=1,3
                coords(l)=coords(l)+co(l,nodes(k))
              enddo
            enddo
            do l=1,3
              coords(l)=coords(l)/3.d0
            enddo
          endif
!     
          if(j.eq.istartset(imast)) then
            if(tieset(1,itie)(81:81).eq.'Z') then
!     
!     check the correct oprientation of the axis vector
!     
              if(lakonf(nelemm)(4:4).eq.'8') then
                nope=8
                nopes=4
              elseif(lakonf(nelemm)(4:4).eq.'6') then
                nope=6
                nopes=nf(jfacem)
              else
                nope=4
                nopes=3
              endif
!     
!     two nodes on the axis
!     
              do k=1,3
                a(k)=cs(5+k,i)
                b(k)=cs(8+k,i)
              enddo
!     
!     looking for a master node not on the axis
!     
              do k=1,nopes
                if((co(1,nodes(k))-a(1))**2+
     &               (co(2,nodes(k))-a(2))**2+
     &               (co(3,nodes(k))-a(3))**2.gt.1.d-20) exit
              enddo
              nodeinface=nodes(k)
!     
!     looking for a node belonging to the same element
!     but not in the master face
!     
              loop:do k=1,nope
              nodeoutface=kon(indexe+k)
              do l=1,nopes
                if(nodes(l).eq.nodeoutface) cycle loop
              enddo
              exit
            enddo loop
!     
!     normalized vector along the axis
!     
            dd=dsqrt((b(1)-a(1))**2+
     &           (b(2)-a(2))**2+
     &           (b(3)-a(3))**2)
            do k=1,3
              xn(k)=(b(k)-a(k))/dd
            enddo
!     
!     vector connecting a on the axis with 
!     nodeinface and nodeoutface
!     
            do k=1,3
              p(k)=co(k,nodeinface)-a(k)
              q(k)=co(k,nodeoutface)-a(k)
            enddo
!     
!     n.(pxq) must be negative for the axis vector
!     to have the right orientation (= you turn from the
!     slave surface to the master surface through the
!     body in clockwise direction while looking in 
!     the direction of xn)
!     
            dot=xn(1)*(p(2)*q(3)-q(2)*p(3))
     &           +xn(2)*(p(3)*q(1)-q(3)*p(1))
     &           +xn(3)*(p(1)*q(2)-q(1)*p(2))
!     
!     revert the direction of the axis if not correct
!     
            if(dot.gt.0.d0) then
              do k=1,3
                b(k)=2.d0*a(k)-b(k)
                cs(8+k,i)=b(k)
                xn(k)=-xn(k)
              enddo
            endif
!     
!     angle from the master to the slave surface
!     
            pi=4.d0*datan(1.d0)
            theta=-2.d0*pi/cs(1,i)
!     
!     rotation matrix rotating a vector in the master
!     surface into a vector in the slave surface
!     
            dc=dcos(theta)
            ds=dsin(theta)
!     
!     C-matrix from Guido Dhondt, The Finite Element
!     Method for Three-Dimensional Thermomechanical
!     Applications p 158
!     
            c(1,1)=dc+(1.d0-dc)*xn(1)*xn(1)
            c(1,2)=   (1.d0-dc)*xn(1)*xn(2)-ds*xn(3)
            c(1,3)=   (1.d0-dc)*xn(1)*xn(3)+ds*xn(2)
            c(2,1)=   (1.d0-dc)*xn(2)*xn(1)+ds*xn(3)
            c(2,2)=dc+(1.d0-dc)*xn(2)*xn(2)
            c(2,3)=   (1.d0-dc)*xn(2)*xn(3)-ds*xn(1)
            c(3,1)=   (1.d0-dc)*xn(3)*xn(1)-ds*xn(2)
            c(3,2)=   (1.d0-dc)*xn(3)*xn(2)+ds*xn(1)
            c(3,3)=dc+(1.d0-dc)*xn(3)*xn(3)
          endif
        endif
!     
!     storing the face centers (translated from the master
!     face into the slave face) in the vectors x,y and z
!     
        n=n+1
        if(tieset(1,itie)(81:81).eq.'P') then
          x(n)=coords(1)-cs(6,i)
          y(n)=coords(2)-cs(7,i)
          z(n)=coords(3)-cs(8,i)
        else
!     
!     vector orthogonal to axis (in the master surface)
!     
          do k=1,3
            p(k)=coords(k)-a(k)
          enddo
          dd=p(1)*xn(1)+p(2)*xn(2)+p(3)*xn(3)
          do k=1,3
            p(k)=p(k)-dd*xn(k)
          enddo
!     
!     vector rotated into the slave surface
!     
          do k=1,3
            q(k)=c(k,1)*p(1)+c(k,2)*p(2)+c(k,3)*p(3)
          enddo
!     
!     rotated node in the slave surface
!     
          x(n)=coords(1)+q(1)-p(1)
          y(n)=coords(2)+q(2)-p(2)
          z(n)=coords(3)+q(3)-p(3)
        endif
      enddo
!     
!     preparing the fields for near3d
!     
      do j=1,n
        nx(j)=j
        ny(j)=j
        nz(j)=j
        xo(j)=x(j)
        yo(j)=y(j)
        zo(j)=z(j)
      enddo
      kflag=2
      call dsort(x,nx,n,kflag)
      call dsort(y,ny,n,kflag)
      call dsort(z,nz,n,kflag)
!     
!     loop over all slave faces
!     
      do j=istartset(islav),iendset(islav)
        ifaces=ialset(j)
!     
        nelems=int(ifaces/10)
        jfaces=ifaces-nelems*10
!     
!     taking the renumbering of the elements into account
!     
        nelems=nactdoh(nelems)
!     
        indexe=ipkonf(nelems)
!     
!     determination of the face center of the slave face
!     
        if(lakonf(nelems)(4:4).eq.'8') then
          do l=1,3
            coords(l)=0.d0
          enddo
          do k=1,4
            nodes(k)=kon(indexe+ifaceq(k,jfaces))
            do l=1,3
              coords(l)=coords(l)+co(l,nodes(k))
            enddo
          enddo
          do l=1,3
            coords(l)=coords(l)/4.d0
          enddo
        elseif(lakonf(nelems)(4:4).eq.'6') then
          do l=1,3
            coords(l)=0.d0
          enddo
          do k=1,nf(jfaces)
            nodes(k)=kon(indexe+ifacew(k,jfaces))
            do l=1,3
              coords(l)=coords(l)+co(l,nodes(k))
            enddo
          enddo
          do l=1,3
            coords(l)=coords(l)/nf(jfaces)
          enddo
        else
          do l=1,3
            coords(l)=0.d0
          enddo
          do k=1,3
            nodes(k)=kon(indexe+ifaceq(k,jfaces))
            do l=1,3
              coords(l)=coords(l)+co(l,nodes(k))
            enddo
          enddo
          do l=1,3
            coords(l)=coords(l)/3.d0
          enddo
        endif
!     
        nneigh=1
        call near3d(xo,yo,zo,x,y,z,nx,ny,nz,coords(1),
     &       coords(2),coords(3),n,noden,nneigh)
!     
        ifacem=ialset(istartset(imast)+noden(1)-1)
!     
        nelemm=int(ifacem/10)
        jfacem=ifacem-nelemm*10
!     
!     taking the renumbering of the elements into account
!     
        nelemm=nactdoh(nelemm)
!     
        ielfa(2,neifa(ipnei(nelems)+jfaces))=nelemm
        ielfa(2,neifa(ipnei(nelemm)+jfacem))=nelems
!     
        neiel(ipnei(nelems)+jfaces)=nelemm
        neiel(ipnei(nelemm)+jfacem)=nelems
!     
        neij(ipnei(nelems)+jfaces)=jfacem
        neij(ipnei(nelemm)+jfacem)=jfaces
!     
!     field indicating whether a face (slave: +, master: -)
!     is part of a cyclic symmetry condition i
!     
        ifatie(neifa(ipnei(nelems)+jfaces))=i
        ifatie(neifa(ipnei(nelemm)+jfacem))=-i
      enddo
      enddo
!     
!     catalogueing external faces
!     
      nfaext=0
!     
      do i=1,nface
        if(ielfa(2,i).ne.0) cycle
        nfaext=nfaext+1
        ifaext(nfaext)=i
        iel1=ielfa(1,i)
        indexf=ipnei(iel1)
        if(lakonf(iel1)(4:4).eq.'8') then
!     
!     hexes
!     
          j=ielfa(4,i)
          j=jopposite8(j)
        elseif(lakonf(iel1)(4:4).eq.'6') then
!     
!     wedges
!     
          j=ielfa(4,i)
          jop=jopposite6(j)
          if(jop.eq.0) then
            jop=j+1
            if(jop.gt.5) jop=3
            if(neiel(indexf+jop).eq.0) then
              jop=j-1
              if(jop.lt.3) jop=5
              if(neiel(indexf+jop).eq.0) cycle
            endif
            j=jop
          endif
        else
          cycle
        endif
        ielfa(3,i)=neiel(indexf+j)
      enddo
!     
c     do i=1,nef
c     write(*,*)'topocfd neiel',i,ipnei(i),(neiel(ipnei(i)+j),j=1,6)
c     write(*,*)'topocfd neij',i,ipnei(i),(neij(ipnei(i)+j),j=1,6)
c     write(*,*)'topocfd neifa',i,ipnei(i),(neifa(ipnei(i)+j),j=1,6)
c     enddo
c     do k=1,nface
c     write(*,*) 'topocfd ielfa',k,(ielfa(j,k),j=1,4),ifatie(k)
c     enddo
c     do k=1,nfaext
c     write(*,*) 'topocfd ifaext',k,ifaext(k)
c     enddo
c     write(*,*)
!     
!     faces belonging to solid surfaces
!     
      nsolidsurf=0
      noset(1:13)='SOLIDSURFACET'
      do i=1,nset
        if(set(i)(1:13).eq.noset(1:13)) exit
      enddo
      if(i.gt.nset) then
        write(*,*) '*WARNING in topocfd: facial surface SOLID SURFACE '
        write(*,*) '         has not been defined.'
        write(*,*)
      else
        do j=istartset(i),iendset(i)
          nsolidsurf=nsolidsurf+1
          isolidsurf(nsolidsurf)=ialset(j)
        enddo
        call isortii(isolidsurf,iaux,nsolidsurf,kflag)
      endif
c     do i=1,nsolidsurf
c     write(*,*) 'topocfd solidsurf ',i,isolidsurf(i)
c     enddo
!     
!     initial conditions: element values
!     (unless some value is nonzero which points to a restart
!     or an ongoing calculation)
!     
      initial=1
      loop1: do i=1,nef
      do j=0,4
        if(vel(i,j).ne.0.d0) then
          initial=0
          exit loop1
        endif
      enddo
      enddo loop1
      if(initial.eq.1) then
        do i=1,nef
          indexe=ipkonf(i)
          if(lakonf(i)(4:4).eq.'8') then
            nope=8
          elseif(lakonf(i)(4:4).eq.'6') then
            nope=6
          else
            nope=4
          endif
          do j=0,4
            do k=1,nope
              vel(i,j)=vel(i,j)+vold(j,kon(indexe+k))
            enddo
            vel(i,j)=vel(i,j)/nope
          enddo
        enddo
      endif
c     do i=1,nef
c     write(*,*) 'topocfd vel ',i,(vel(j,i),j=1,3)
c     enddo
!     
      return
      end
