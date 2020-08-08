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
      subroutine couplings(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nboun,nk,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,ikboun,ikmpc,ilmpc,co,labmpc,istat,n,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,norien,orname,orab,irstrt,ipkon,
     &  kon,lakon,istep,ics,dcs,nk_,nboun_,nodeboun,ndirboun,
     &  typeboun,ilboun,xboun,ier)
!
!     reading the input deck: *COUPLING in combination with
!     *KINEMATIC or *DISTRIBUTING
!
      implicit none
!
      character*1 inpc(*),surfkind,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*80 orname(*),orientation
      character*81 set(*),surfset
      character*132 textpart(16),name
!
      integer istartset(*),iendset(*),ialset(*),norien,irstrt(*),nface,
     &  iorientation,iface,jface,nset,nboun,istat,n,i,j,k,ibounstart,
     &  ibounend,key,nk,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikboun(*),ikmpc(*),ilmpc(*),ipos,m,node,iline,ipol,inl,nope,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),jsurf,irefnode,indexe,nopes,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),ipkon(*),
     &  kon(*),istep,nelem,ics(2,*),nodef(8),nboun_,nodeboun(*),
     &  npt,mint,index1,mpcfreeold,id,idof,iflag,inhomnode,nk_,kk,
     &  irotnode(3),l,index2,indexold(3),indexnew(3),idirold,idirmax,
     &  idir,indexmax,nodeold,node1,nodemax,ndirboun(*),ilboun(*),
     &  irotnode_kin,idupnode,ier
!
      real*8 coefmpc(*),co(3,*),orab(7,*),dcs(*),areanodal(8),xl2(3,8),
     &  shp2(7,8),xsj2(3),xsj,xi,et,weight,xs2(3,2),area,a(3,3),cgx(3),
     &  aa(3),pi(3),c1,c4,c9,c10,amax,xcoef,coef(3),xboun(*),stdev
!
      include "gauss.f"
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
!     flag for shape functions
!
      data iflag /2/
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *COUPLING: *COUPLING'
         write(*,*)'       should be placed before all step definitions'
         ier=1
         return
      endif
!
      label='                    '
      orientation='
     &                           '
      do i=1,81
         surfset(i:i)=' '
      enddo
!
      name(1:1)=' '
      do i=2,n
         if(textpart(i)(1:8).eq.'REFNODE=') then
            read(textpart(i)(9:18),'(i10)',iostat=istat) irefnode
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*COUPLING%",ier)
               return
            endif
            if(irefnode.gt.nk) then
               write(*,*) '*ERROR reading *COUPLING: ref node',irefnode
               write(*,*) '       has not been defined'
               ier=1
               return
            endif
         else if(textpart(i)(1:8).eq.'SURFACE=') then
            surfset(1:80)=textpart(i)(9:88)
!     
            ipos=index(surfset,' ')
            surfkind='S'
            surfset(ipos:ipos)=surfkind
            do j=1,nset
               if(set(j).eq.surfset) exit
            enddo
            if(j.gt.nset) then
               surfkind='T'
               surfset(ipos:ipos)=surfkind
               do j=1,nset
                  if(set(j).eq.surfset) exit
               enddo
               if(j.gt.nset) then
                  write(*,*) '*ERROR reading *COUPLING:'
                  write(*,*) '       surface ',surfset
                  write(*,*) '       has not yet been defined.' 
                  ier=1
                  return
               endif
            endif
            jsurf=j
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:15).eq.'CONSTRAINTNAME=') then
            name(1:117)=textpart(i)(16:132)
         else
            write(*,*) 
     &           '*WARNING reading *COUPLING: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*COUPLING%")
         endif
      enddo
!
      if(name(1:1).eq.' ') then
         write(*,*)
     &        '*ERROR reading *COUPLING: no CONTRAINT NAME given'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*COUPLING%",ier)
         return
      endif
!
      if(orientation.eq.'                    ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &       '*ERROR reading *COUPLING: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*COUPLING%",ier)
            return
         endif
         iorientation=i
      endif
!
!     next keyword should be *KINEMATIC or *DISTRIBUTING
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if(textpart(1)(2:10).eq.'KINEMATIC') then
!
!        catalogueing the nodes
!
         npt=0
!
         do j=istartset(jsurf),iendset(jsurf)
            if(ialset(j).gt.0) then
               if(surfkind.eq.'T') then
!     
!                 facial surface
!     
                  iface=ialset(j)
                  nelem=int(iface/10)
                  jface=iface-nelem*10
                  indexe=ipkon(nelem)
!     
                  if(lakon(nelem)(4:5).eq.'20') then
                     nopes=8
                  elseif(lakon(nelem)(4:4).eq.'2') then
                     nopes=9
                  elseif(lakon(nelem)(4:4).eq.'8') then
                     nopes=4
                  elseif(lakon(nelem)(4:5).eq.'10') then
                     nopes=6
                  elseif(lakon(nelem)(4:4).eq.'4') then
                     nopes=3
                  endif
!     
                  if(lakon(nelem)(4:4).eq.'6') then
                     if(jface.le.2) then
                        nopes=3
                     else
                        nopes=4
                     endif
                  endif
                  if(lakon(nelem)(4:5).eq.'15') then
                     if(jface.le.2) then
                        nopes=6
                     else
                        nopes=8
                     endif
                  endif   
               else
!
!                 nodal surface
!
                  nopes=1
               endif
!     
               do m=1,nopes
                  if(surfkind.eq.'T') then
                     if((lakon(nelem)(4:4).eq.'2').or.
     &                    (lakon(nelem)(4:4).eq.'8')) then
                        node=kon(indexe+ifaceq(m,jface))
                     elseif((lakon(nelem)(4:4).eq.'4').or.
     &                       (lakon(nelem)(4:5).eq.'10')) then
                        node=kon(indexe+ifacet(m,jface))
                     elseif(lakon(nelem)(4:4).eq.'6') then
                        node=kon(indexe+ifacew1(m,jface))
                     elseif(lakon(nelem)(4:5).eq.'15') then
                        node=kon(indexe+ifacew2(m,jface))
                     endif
                  else
                     node =ialset(j)
                  endif
!
                  call nident2(ics,node,npt,id)
                  if(id.gt.0) then
                     if(ics(1,id).eq.node) then
                        cycle
                     endif
                  endif
!
!                 updating ics
!
                  npt=npt+1
                  do l=npt,id+2,-1
                     ics(1,l)=ics(1,l-1)
                  enddo
                  ics(1,id+1)=node
               enddo
            else
!
!           if a negative value occurs the surface has to be
!           nodal
!
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  node=k
                  call nident2(ics,node,npt,id)
                  if(id.gt.0) then
                     if(ics(1,id).eq.node) then
                        cycle
                     endif
                  endif
!
!              updating ics
!
                  npt=npt+1
                  do l=npt,id+2,-1
                     ics(1,l)=ics(1,l-1)
                  enddo
                  ics(1,id+1)=node
               enddo
            endif
         enddo
!
!        generating a rotational node and connecting the
!        rotational dofs of the reference node with the
!        translational dofs of the rotational node
!
!        generating a rotational reference node
!
         nk=nk+1
         if(nk.gt.nk_) then
            write(*,*) 
     &           '*ERROR reading *KINEMATIC: increase nk_'
            ier=1
            return
         endif
         irotnode_kin=nk
         do l=1,3
            co(l,nk)=co(l,irefnode)
         enddo
!
!        generating connecting MPCs between the rotational
!        dofs of irefnode and the translational dofs of
!        irotnode_kin
!
         do k=1,3
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) 
     &              '*ERROR reading *KINEMATIC: increase nmpc_'
               ier=1
               return
            endif
!     
!           the internal dofs for rotation are 4, 5 and 6
!     
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)='ROTTRACOUPLING      '
            idof=8*(irefnode-1)+k+3
            call nident(ikmpc,idof,nmpc-1,id)
            do l=nmpc,id+2,-1
               ikmpc(l)=ikmpc(l-1)
               ilmpc(l)=ilmpc(l-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            nodempc(1,mpcfree)=irefnode
            nodempc(2,mpcfree)=k+4
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
!     
            nodempc(1,mpcfree)=irotnode_kin
            nodempc(2,mpcfree)=k
            coefmpc(mpcfree)=-1.d0
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            nodempc(3,mpcfreeold)=0
         enddo
!
         if(iorientation.gt.0) then
!
!           duplicating the nodes
!           generating rigid body MPC's for all dofs in the new nodes
!
            do m=1,npt
               nk=nk+1
               if(nk.gt.nk_) then
                  write(*,*) 
     &                 '*ERROR reading *KINEMATIC: increase nk_'
                  ier=1
                  return
               endif
               ics(2,m)=nk
               do k=1,3
                  co(k,nk)=co(k,ics(1,m))
               enddo
               node=nk
               ibounstart=1
               ibounend=3
               call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &              irotnode_kin,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &              nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,node,
     &              typeboun,co,ibounstart,ibounend)
            enddo
         endif
!
!        reading the degrees of freedom
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
!     
            read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            endif
            if(ibounstart.lt.1) then
               write(*,*) '*ERROR reading *KINEMATIC'
               write(*,*) '       starting degree of freedom cannot'
               write(*,*) '       be less than 1'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            endif
!     
            if(textpart(2)(1:1).eq.' ') then
               ibounend=ibounstart
            else
               read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*BOUNDARY%",ier)
                  return
               endif
            endif
            if(ibounend.gt.3) then
               write(*,*) '*ERROR reading *KINEMATIC'
               write(*,*) '       final degree of freedom cannot'
               write(*,*) '       exceed 3'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            elseif(ibounend.lt.ibounstart) then
               write(*,*) '*ERROR reading *KINEMATIC'
               write(*,*) '       initial degree of freedom cannot'
               write(*,*) '       exceed final degree of freedom'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            endif
!
!           generating the MPCs
!
            if(iorientation.eq.0) then
!
!              generating rigid body MPC's for the appropriate dofs
!
               do j=1,npt
                  node=ics(1,j)
                  call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &                 irotnode_kin,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &                 nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &                 node,typeboun,co,ibounstart,ibounend)
               enddo
            else
!
!              connecting the original nodes with the duplicated nodes
!              for the appropriate dofs
!
               do j=1,npt
                  node=ics(1,j)
                  idupnode=ics(2,j)
                  call mpcadd(node,ibounstart,ibounend,nboun,ipompc,
     &                 nodempc,coefmpc,nmpc,nmpc_,mpcfree,orab,ikboun,
     &                 ikmpc,ilmpc,co,labmpc,label,idupnode,
     &                 iorientation)
               enddo
            endif
         enddo
      elseif(textpart(1)(2:13).eq.'DISTRIBUTING') then
         if(surfkind.eq.'S') then
            write(*,*) '*ERROR reading *DISTRIBUTING'
            write(*,*) '       a nodal surface is not allowed'
            write(*,*) '       please use a facial surface on'
            write(*,*) '       the *COUPLING card'
            ier=1
            return
         endif
!
         npt=0
         area=0.d0
!
!        catalogueing the nodes belonging to the surface (ics(1,*))
!        catalogueing the area (dcs(*))
!
         do k=istartset(jsurf),iendset(jsurf)
!     
!           facial surface
!               
            iface=ialset(k)
            nelem=int(iface/10)
            jface=iface-nelem*10
            indexe=ipkon(nelem)
!
!        nodes: #nodes in the face
!        the nodes are stored in nodef(*)
!
            if(lakon(nelem)(4:4).eq.'2') then
               nopes=8
               nface=6
            elseif(lakon(nelem)(3:4).eq.'D8') then
               nopes=4
               nface=6
            elseif(lakon(nelem)(4:5).eq.'10') then
               nopes=6
               nface=4
               nope=10
            elseif(lakon(nelem)(4:4).eq.'4') then
               nopes=3
               nface=4
               nope=4
            elseif(lakon(nelem)(4:5).eq.'15') then
               if(jface.le.2) then
                  nopes=6
               else
                  nopes=8
               endif
               nface=5
               nope=15
            elseif(lakon(nelem)(3:4).eq.'D6') then
               if(jface.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
               nface=5
               nope=6
            else
               cycle
            endif
!     
!     determining the nodes of the face
!     
            if(nface.eq.4) then
               do i=1,nopes
                  nodef(i)=kon(indexe+ifacet(i,jface))
               enddo
            elseif(nface.eq.5) then
               if(nope.eq.6) then
                  do i=1,nopes
                     nodef(i)=kon(indexe+ifacew1(i,jface))
                  enddo
               elseif(nope.eq.15) then
                  do i=1,nopes
                     nodef(i)=kon(indexe+ifacew2(i,jface))
                  enddo
               endif
            elseif(nface.eq.6) then
               do i=1,nopes
                  nodef(i)=kon(indexe+ifaceq(i,jface))
               enddo
            endif
!
!        loop over the nodes belonging to the face   
!        ics(1,*): pretension node
!        dcs(*): area corresponding to pretension node   
!         
            do i=1,nopes
               node=nodef(i)
               call nident2(ics,node,npt,id)
               if(id.gt.0) then
                  if(ics(1,id).eq.node) then
                     cycle
                  endif
               endif
!
!              updating ics
!
               npt=npt+1
               do j=npt,id+2,-1
                  ics(1,j)=ics(1,j-1)
                  dcs(j)=dcs(j-1)
               enddo
               ics(1,id+1)=node
               dcs(id+1)=0.d0
            enddo
!
!        calculating the area of the face and its contributions
!        to the facial nodes
!
!        number of integration points
!         
            if(lakon(nelem)(3:5).eq.'D8R') then
               mint=1
            elseif(lakon(nelem)(3:4).eq.'D8') then
               mint=4
            elseif(lakon(nelem)(4:6).eq.'20R') then
               mint=4
            elseif(lakon(nelem)(4:4).eq.'2') then
               mint=9
            elseif(lakon(nelem)(4:5).eq.'10') then
               mint=3
            elseif(lakon(nelem)(4:4).eq.'4') then
               mint=1
            elseif(lakon(nelem)(3:4).eq.'D6') then
               mint=1
            elseif(lakon(nelem)(4:5).eq.'15') then
               if(jface.le.2) then
                  mint=3
               else
                  mint=4
               endif
            endif
!     
            do i=1,nopes
               areanodal(i)=0.d0
               do j=1,3
                  xl2(j,i)=co(j,nodef(i))
               enddo
            enddo
!     
            do m=1,mint
               if((lakon(nelem)(3:5).eq.'D8R').or.
     &              ((lakon(nelem)(3:4).eq.'D6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,m)
                  et=gauss2d1(2,m)
                  weight=weight2d1(m)
               elseif((lakon(nelem)(3:4).eq.'D8').or.
     &                 (lakon(nelem)(4:6).eq.'20R').or.
     &                 ((lakon(nelem)(4:5).eq.'15').and.
     &                 (nopes.eq.8))) then
                  xi=gauss2d2(1,m)
                  et=gauss2d2(2,m)
                  weight=weight2d2(m)
               elseif(lakon(nelem)(4:4).eq.'2') then
                  xi=gauss2d3(1,m)
                  et=gauss2d3(2,m)
                  weight=weight2d3(m)
               elseif((lakon(nelem)(4:5).eq.'10').or.
     &                 ((lakon(nelem)(4:5).eq.'15').and.
     &                 (nopes.eq.6))) then
                  xi=gauss2d5(1,m)
                  et=gauss2d5(2,m)
                  weight=weight2d5(m)
               elseif((lakon(nelem)(4:4).eq.'4').or.
     &                 ((lakon(nelem)(3:4).eq.'D6').and.
     &                 (nopes.eq.3))) then
                  xi=gauss2d4(1,m)
                  et=gauss2d4(2,m)
                  weight=weight2d4(m)
               endif
!     
               if(nopes.eq.8) then
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.3) then
                  call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               endif
!     
!     calculating the total area and nodal area
!     
               xsj=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
               area=area+xsj
               do i=1,nopes
                  areanodal(i)=areanodal(i)+xsj*shp2(4,i)
               enddo
!     
            enddo
!     
!     inserting the nodal area into field dcs
!     
            do i=1,nopes
               node=nodef(i)
               call nident2(ics,node,npt,id)
               dcs(id)=dcs(id)+areanodal(i)
            enddo
!     
         enddo
!
!        create for each node in the surface a new node.
!        the ratio of the displacements between both nodes
!        is governed by the area for which each node is
!        representative
!
!        initializing the location of the center of gravity
!
         do k=1,3
            cgx(k)=0.d0
         enddo
!
         do m=1,npt
            nk=nk+1
            if(nk.gt.nk_) then
               write(*,*) 
     &              '*ERROR reading *DISTRIBUTING: increase nk_'
               ier=1
               return
            endif
            node=ics(1,m)
            ics(1,m)=nk
            do k=1,3
               co(k,nk)=co(k,node)
            enddo
!     
!           generate the connecting MPC's
!     
            do k=1,3
!
!              contribution to the location of the center of gravity
!
               cgx(k)=cgx(k)+dcs(m)*co(k,node)
c               cgx(k)=cgx(k)+co(k,node)
!     
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR reading *DISTRIBUTING: increase nmpc_'
                  ier=1
                  return
               endif
!     
!              MPC: u(old node)= 
!              (mean area)/(area_m) * u(new node) (in all directions)
!
!              check whether MPC was already used
!
               ipompc(nmpc)=mpcfree
               labmpc(nmpc)='                    '
               idof=8*(node-1)+k
               call nident(ikmpc,idof,nmpc-1,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     write(*,*) '*ERROR reading *COUPLING: dof',k
                     write(*,*) '       in node ',node
                     write(*,*) '       was already used'
                     ier=1
                     return
                  endif
               endif
               do l=nmpc,id+2,-1
                  ikmpc(l)=ikmpc(l-1)
                  ilmpc(l)=ilmpc(l-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
!     
!              generating the terms of the MPC
!     
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=k
               coefmpc(mpcfree)=1.d0
               mpcfree=nodempc(3,mpcfree)
!     
               nodempc(1,mpcfree)=nk
               nodempc(2,mpcfree)=k
               coefmpc(mpcfree)=-area/(npt*dcs(m))
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
               nodempc(3,mpcfreeold)=0
            enddo
         enddo
!
         do k=1,3
            cgx(k)=cgx(k)/area
c            cgx(k)=cgx(k)/npt
         enddo
!
!        calculating a standard deviation; this quantity will
!        serve as a limit for checking the closeness of individual
!        nodes to the center of gravity
!
         stdev=0.d0
         do m=1,npt
            node=ics(1,m)
            do k=1,3
c               stdev=stdev+(dcs(m)*co(k,node)-cgx(k))**2
               stdev=stdev+(co(k,node)-cgx(k))**2
            enddo
         enddo
         stdev=stdev/npt
!
!        generating the translational MPC's (default)
!
         do m=1,3
            node=ics(1,1)
            idof=8*(node-1)+m
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR reading *COUPLING: increase nmpc_'
               ier=1
               return
            endif
            labmpc(nmpc)='                    '
            ipompc(nmpc)=mpcfree
!     
!           updating ikmpc and ilmpc
!     
            do j=nmpc,id+2,-1
               ikmpc(j)=ikmpc(j-1)
               ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!
            do j=1,npt
               node=ics(1,j)
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=m
               coefmpc(mpcfree)=1.d0
               mpcfree=nodempc(3,mpcfree)
            enddo
            nodempc(1,mpcfree)=irefnode
            nodempc(2,mpcfree)=m
            coefmpc(mpcfree)=-npt
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            nodempc(3,mpcfreeold)=0
         enddo
!
!        generating the rotational MPC's
!
         irotnode(1)=0
!
!        reading the degrees of freedom
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
!     
            read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            endif
            if(ibounstart.lt.1) then
               write(*,*) '*ERROR reading *KINEMATIC'
               write(*,*) '       starting degree of freedom cannot'
               write(*,*) '       be less than 1'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*KINEMATIC%",ier)
               return
            endif
!     
            if(textpart(2)(1:1).eq.' ') then
               ibounend=ibounstart
            else
               read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*BOUNDARY%",ier)
                  return
               endif
            endif
!
            ibounstart=max(4,ibounstart)
            if(ibounend.gt.6) then
               write(*,*) '*ERROR reading *DISTRIBUTING'
               write(*,*) '       final degree of freedom cannot'
               write(*,*) '       exceed 6'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTING%",ier)
               return
            elseif(ibounend.lt.ibounstart) then
               cycle
            endif
!
            if((ibounend.gt.3).and.(irotnode(1).eq.0)) then
!
!              check the orientation definition: the local reference
!              system is not allowed to be cylindrical
!
               if(iorientation.ne.0) then
                  if(orab(7,iorientation).lt.0.d0) then
                     write(*,*) '*ERROR reading *DISTRIBUTING'
                     write(*,*) '       a cylindrical local coordinate'
                     write(*,*) '       system is not allowed'
                     ier=1
                     return
                  endif
!
                  call transformatrix(orab(1,iorientation),
     &                 co(1,irefnode),a)
               endif
!
!              generating a rotational reference node
!
               do k=1,3
                  nk=nk+1
                  if(nk.gt.nk_) then
                     write(*,*) 
     &                    '*ERROR reading *DISTRIBUTING: increase nk_'
                     ier=1
                     return
                  endif
                  irotnode(k)=nk
                  do l=1,3
                     co(l,nk)=co(l,irefnode)
                  enddo
               enddo
!
!              generating connecting MPCs between the rotational
!              dofs of irefnode and the translational dofs of
!              irotnode
!
               do k=1,3
!     
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                 '*ERROR reading *DISTRIBUTING: increase nmpc_'
                     ier=1
                     return
                  endif
!     
!                 the internal dofs for rotation are 4, 5 and 7
!     
                  ipompc(nmpc)=mpcfree
                  labmpc(nmpc)='ROTTRACOUPLING      '
                  idof=8*(irefnode-1)+k+3
                  call nident(ikmpc,idof,nmpc-1,id)
                  do l=nmpc,id+2,-1
                     ikmpc(l)=ikmpc(l-1)
                     ilmpc(l)=ilmpc(l-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
!     
                  nodempc(1,mpcfree)=irefnode
                  nodempc(2,mpcfree)=k+4
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
!
                  nodempc(1,mpcfree)=irotnode(k)
                  nodempc(2,mpcfree)=1
                  coefmpc(mpcfree)=-1.d0
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  nodempc(3,mpcfreeold)=0
               enddo
!
!              generating a node for the inhomogeneous term
!
               nk=nk+1
               if(nk.gt.nk_) then
                  write(*,*) 
     &               '*ERROR reading *DISTRIBUTING: increase nk_'
                  ier=1
                  return
               endif
               inhomnode=nk
!
!              fictituous center of gravity (for compatibility
!              reasons with usermpc.f)
!
c               do k=1,3
c                  cgx(k)=co(k,irefnode)
c               enddo
            endif
!
!           generating rotational MPC's
!                  
            do kk=ibounstart,ibounend
               k=kk-3
!
!              axis of rotation
!
               if(iorientation.eq.0) then
                  do j=1,3
                     aa(j)=0.d0
                  enddo
                  aa(k)=1.d0
               else
                  do j=1,3
                     aa(j)=a(j,k)
                  enddo
               endif
c               write(*,*) 'couplings aa ',(aa(j),j=1,3)
!
!              storing the axis of rotation as coordinates of the
!              rotation nodes
!
               do j=1,3
                  co(j,irotnode(k))=aa(j)
               enddo
!
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR reading *DISTRIBUTING: increase nmpc_'
                  ier=1
                  return
               endif
!     
               ipompc(nmpc)=mpcfree
c               labmpc(nmpc)='MEANROTDISTRIB      '
               labmpc(nmpc)='                    '
!
!              defining the terms of the MPC
!
c               write(*,*) 'couplings, npt ',npt
               do m=1,npt
                  do j=1,3
                     nodempc(1,mpcfree)=ics(1,m)
                     nodempc(2,mpcfree)=0
                     coefmpc(mpcfree)=0.d0
                     mpcfree=nodempc(3,mpcfree)
                  enddo
c                  write(*,*) 'couplings, co',m,(co(j,ics(1,m)),j=1,3)
               enddo
               nodempc(1,mpcfree)=irotnode(k)
               nodempc(2,mpcfree)=1
               coefmpc(mpcfree)=0.d0
               mpcfree=nodempc(3,mpcfree)
               nodempc(1,mpcfree)=inhomnode
               nodempc(2,mpcfree)=k
c changed on 2. november 2018
               coefmpc(mpcfree)=1.d0
c               coefmpc(mpcfree)=0.d0
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
               nodempc(3,mpcfreeold)=0
!
!              storing the inhomogeneous value as a boundary
!              condition in the inhomogeneous node
!
               idof=8*(inhomnode-1)+k
               call nident(ikboun,idof,nboun,id)
               nboun=nboun+1
               if(nboun.gt.nboun_) then
                  write(*,*) '*ERROR reading *COUPLING: increase nboun_'
                  ier=1
                  return
               endif
               nodeboun(nboun)=inhomnode
               ndirboun(nboun)=k
               xboun(nboun)=0.d0
               typeboun(nboun)='U'
               do l=nboun,id+2,-1
                  ikboun(l)=ikboun(l-1)
                  ilboun(l)=ilboun(l-1)
               enddo
               ikboun(id+1)=idof
               ilboun(id+1)=nboun
c               write(*,*) 'couplings cgx',(cgx(l),l=1,3)
c               write(*,*) 'couplings stdev',stdev
!
!              calculating the coefficients of the rotational MPC
!
!              loop over all nodes belonging to the mean rotation MPC
!
               index2=ipompc(nmpc)
               do
                  node=nodempc(1,index2)
                  if(node.eq.irotnode(k)) exit
!     
!     relative positions
!     
                  do j=1,3
                     pi(j)=co(j,node)-cgx(j)
                  enddo
c                  write(*,*) 'couplings pi',(pi(j),j=1,3)
!     
!              projection on a plane orthogonal to the rotation vector
!
                  c1=pi(1)*aa(1)+pi(2)*aa(2)+pi(3)*aa(3)
                  do j=1,3
                     pi(j)=pi(j)-c1*aa(j)
                  enddo
c                  write(*,*) 'couplings c1 pi',c1,(pi(j),j=1,3)
!     
                  c1=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
c                  if(c1.lt.1.d-20) then
                  if(c1.lt.stdev*1.d-10) then
                     write(*,*) 
     &                    '*WARNING reading *DISTRIBUTING: node ',node
                     write(*,*) '         is very close to the '
                     write(*,*) '         rotation axis through the'
                     write(*,*) '         center of gravity of'
                     write(*,*) '         the nodal cloud in a'
                     write(*,*) '         mean rotation MPC.'
                     write(*,*) '         This node is not taken'
                     write(*,*) '         into account in the MPC'
                     index2=nodempc(3,nodempc(3,nodempc(3,index2)))
                     cycle
                  endif
!     
                  do j=1,3
                     if(j.eq.1) then
                        c4=aa(2)*pi(3)-aa(3)*pi(2)
                     elseif(j.eq.2) then
                        c4=aa(3)*pi(1)-aa(1)*pi(3)
                     else
                        c4=aa(1)*pi(2)-aa(2)*pi(1)
                     endif
                     c9=c4/c1
c                     write(*,*) 'couplings c4,c9',j,c4,c9
!     
                     index1=ipompc(nmpc)
                     do
                        node1=nodempc(1,index1)
                        if(node1.eq.irotnode(k)) exit
                        if(node1.eq.node) then
                           c10=c9*(1.d0-1.d0/real(npt))
                        else
                           c10=-c9/real(npt)
                        endif
                        if(j.eq.1) then
                           coefmpc(index1)=coefmpc(index1)+c10
c                           write(*,*) 'couplings c10',
c     &                        j,c10,coefmpc(index1)
                        elseif(j.eq.2) then
                           coefmpc(nodempc(3,index1))=
     &                          coefmpc(nodempc(3,index1))+c10
c                           write(*,*) 'couplings c10',
c     &                        j,c10,coefmpc(nodempc(3,index1))
                        else
                           coefmpc(nodempc(3,nodempc(3,index1)))=
     &                         coefmpc(nodempc(3,nodempc(3,index1)))+c10
c                           write(*,*) 'couplings c10',
c     &                       j,c10,coefmpc(nodempc(3,nodempc(3,index1)))
                        endif
                        index1=
     &                       nodempc(3,nodempc(3,nodempc(3,index1)))
                     enddo
                  enddo
                  index2=nodempc(3,nodempc(3,nodempc(3,index2)))
               enddo
               coefmpc(index2)=-npt
!     
!     assigning the degrees of freedom
!
               j=0
               index2=ipompc(nmpc)
               do
                  j=j+1
                  if(j.gt.3) j=1
                  nodempc(2,index2)=j
c                  write(*,*) 'couplings a',(coefmpc(index2))
                  index2=nodempc(3,index2)
                  if(nodempc(1,index2).eq.irotnode(k)) exit
               enddo
!
!              look for biggest coefficient of all but the last
!              regular node. The general form of the MPC is:
!              regular nodes, rotationl dof and inhomogeneous dof
!
               indexmax=0
               index2=ipompc(nmpc)
               amax=1.d-5
!
!              coefficients of the first terms in the translational
!              dofs
!
               coef(1)=coefmpc(index2)
               coef(2)=coefmpc(nodempc(3,index2))
               coef(3)=coefmpc(nodempc(3,nodempc(3,index2)))
!
!              loop over all regular nodes
!
               do i=1,3*npt
                  if(dabs(coefmpc(index2)).gt.amax) then
                     idir=nodempc(2,index2)
!
!                    in the translational mpcs the coefficients of
!                    all npt terms are equal to 1
!                    in the rotational mpcs the coefficient of the
!                    dependent term should be different from the
!                    coefficient of the node corresponding to the
!                    translational dependent term with the same dof
!
c                     if(dabs(coefmpc(index2)-coef(idir)).lt.1.d-10) then
                     if(dabs(coefmpc(index2)-coef(idir)).lt.
     &                       dabs(coefmpc(index2))/1000.)then
                        index2=nodempc(3,index2)
                        cycle
                     endif
!     
                     node=nodempc(1,index2)
                     idof=8*(node-1)+idir
!     
!                    check whether the node was used in another MPC
!
                     call nident(ikmpc,idof,nmpc-1,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           index2=nodempc(3,index2)
                           cycle
                        endif
                     endif
!
!                    check whether the node was used in a SPC
!
                     call nident(ikboun,idof,nboun,id)
                     if(id.gt.0) then
                        if(ikboun(id).eq.idof) then
                           index2=nodempc(3,index2)
                           cycle
                        endif
                     endif
!     
                     amax=dabs(coefmpc(index2))
                     indexmax=index2
                  endif
!
!                 after each node has been treated, a check is performed
!                 whether indexmax is already nonzero. Taking too many
!                 nodes into account may lead to dependent equations if
!                 all rotational dofs are constrained.
!
                  if((i/3)*3.eq.i) then
                     if(indexmax.ne.0) exit
                  endif
!     
                  index2=nodempc(3,index2)
               enddo
!     
               if(indexmax.eq.0) then
                  write(*,*) 
     &                 '*WARNING reading *DISTRIBUTING: no MPC is '
                  write(*,*) '         generated for the mean'
                  write(*,*) '         rotation in node ',irefnode
                  write(*,*) '         and direction ',kk
!
!                 removing the MPC
!
                  mpcfreeold=mpcfree
                  index2=ipompc(nmpc)
                  ipompc(nmpc)=0
                  mpcfree=index2
!     
!                 removing the MPC from fields nodempc and coefmpc
!     
                  do
                     nodempc(1,index2)=0
                     nodempc(2,index2)=0
                     coefmpc(index2)=0
                     if(nodempc(3,index2).ne.0) then
                        index2=nodempc(3,index2)
                     else
                        nodempc(3,index2)=mpcfreeold
                        exit
                     endif
                  enddo
!     
!                 decrementing nmpc
!     
                  nmpc=nmpc-1
!
                  cycle
               endif
!     
               nodemax=nodempc(1,indexmax)
               idirmax=nodempc(2,indexmax)
               index2=ipompc(nmpc)
!     
               nodeold=nodempc(1,index2)
               idirold=nodempc(2,index2)
!
!              exchange the node information in the MPC
!
               if(nodemax.ne.nodeold) then
!     
                  indexold(1)=index2
                  indexold(2)=nodempc(3,index2)
                  indexold(3)=nodempc(3,indexold(2))
!     
                  index2=nodempc(3,indexold(3))
                  do
                     if(nodempc(1,index2).eq.irotnode(k)) exit
                     if(nodempc(1,index2).eq.nodemax) then
                        if(nodempc(2,index2).eq.1) then
                           indexnew(1)=index2
                        elseif(nodempc(2,index2).eq.2) then
                           indexnew(2)=index2
                        else
                           indexnew(3)=index2
                        endif
                     endif
                     index2=nodempc(3,index2)
                  enddo
!     
                  do j=1,3
                     node=nodempc(1,indexold(j))
                     idir=nodempc(2,indexold(j))
                     xcoef=coefmpc(indexold(j))
                     nodempc(1,indexold(j))=nodempc(1,indexnew(j))
                     nodempc(2,indexold(j))=nodempc(2,indexnew(j))
                     coefmpc(indexold(j))=coefmpc(indexnew(j))
                     nodempc(1,indexnew(j))=node
                     nodempc(2,indexnew(j))=idir
                     coefmpc(indexnew(j))=xcoef
                  enddo
               endif
!
!              exchange the direction information in the MPC
!
               index2=ipompc(nmpc)
               if(idirmax.ne.1) then
                  indexold(1)=index2
                  if(idirmax.eq.2) then
                     indexnew(1)=nodempc(3,index2)
                  else
                     indexnew(1)=nodempc(3,nodempc(3,index2))
                  endif
!     
                  do j=1,1
                     idir=nodempc(2,indexold(j))
                     xcoef=coefmpc(indexold(j))
                     nodempc(2,indexold(j))=nodempc(2,indexnew(j))
                     coefmpc(indexold(j))=coefmpc(indexnew(j))
                     nodempc(2,indexnew(j))=idir
                     coefmpc(indexnew(j))=xcoef
                  enddo
               endif
!     
!              determining the dependent dof of the MPC
!     
               idof=8*(nodemax-1)+idirmax
               call nident(ikmpc,idof,nmpc-1,id)
               do l=nmpc,id+2,-1
                  ikmpc(l)=ikmpc(l-1)
                  ilmpc(l)=ilmpc(l-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
!
            enddo
         enddo
!
!        generating the translational MPC's (default)
!
c         do m=1,3
c            node=ics(1,1)
c            idof=8*(node-1)+m
c            call nident(ikmpc,idof,nmpc,id)
c!     
c            nmpc=nmpc+1
c            if(nmpc.gt.nmpc_) then
c               write(*,*) '*ERROR reading *COUPLING: increase nmpc_'
c               ier=1
c               return
c            endif
c            labmpc(nmpc)='                    '
c            ipompc(nmpc)=mpcfree
c!     
c!           updating ikmpc and ilmpc
c!     
c            do j=nmpc,id+2,-1
c               ikmpc(j)=ikmpc(j-1)
c               ilmpc(j)=ilmpc(j-1)
c            enddo
c            ikmpc(id+1)=idof
c            ilmpc(id+1)=nmpc
c!
c            do j=1,npt
c               node=ics(1,j)
c               nodempc(1,mpcfree)=node
c               nodempc(2,mpcfree)=m
c               coefmpc(mpcfree)=1.d0
c               mpcfree=nodempc(3,mpcfree)
c            enddo
c            nodempc(1,mpcfree)=irefnode
c            nodempc(2,mpcfree)=m
c            coefmpc(mpcfree)=-npt
c            mpcfreeold=mpcfree
c            mpcfree=nodempc(3,mpcfree)
c            nodempc(3,mpcfreeold)=0
c         enddo
      else
         write(*,*)
     &        '*ERROR reading *COUPLING: the line following'
         write(*,*) '       *COUPLING must contain the'
         write(*,*) '       *KINEMATIC keyword or the'
         write(*,*) '       *DISTRIBUTING keyword'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*COUPLING%",ier)
         return
      endif
!
      return
      end

