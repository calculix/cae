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
! >
! >  \brief subroutine to catalogue SPC'S and MPC'S for master and slave node
! >     needed for mortar contact
! > Author: Saskia Sitzmann
! > e-mail: saskia.sitzmann@fau.de
! >
! > @param [in] lakon              (i) label for element i
! > @param [in] ipkon              pointer into field kon...
! > @param [in] kon               .. for element i storing the connectivity list of elem. in succ. order
! > @param [in] ntie              number of ties
! > @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface
! > @param [in] nset              number of sets
! > @param [in] set              (i) name of set i
! > @param [in] itiefac               pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
! > @param [in] islavsurf       islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i
! > @param [in] islavnode       field storing the nodes of the slave surface
! > @param [in] imastnode       field storing the nodes of the master surfaces
! > @param [in] nslavnode       (i)pointer into field isalvnode for contact tie i
! > @param [in] nmastnode       (i)pointer into field imastnode for contact tie i
! > @param [in] nslavs              number of slave nodes
! > @param [in] iponoels         (i) pointer to inoels
! > @param [in] inoels           (3,i) element number, local node number and pointer to another entry
! > @param [in] nk              number of nodes
! > @param [in] nboun            number of SPCs
! > @param [in] ndirboun              (i) direction of SPC i
! > @param [in] nodeboun         (i) node of SPC i
! > @param [in] xboun            (i) value of SPC i
! > @param [in] nmpc              number of mpcs
! > @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
! > @param [in] nodempc          nodes and directions of MPCs
! > @param [in] coefmpc          coefficients of MPCs
! > @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
! > @param [in] ilboun           SPC numbers for sorted dofs
! > @param [in] ikmpc               sorted dofs idof=8*(node-1)+dir for MPCs
! > @param [in] ilmpc              SPC numbers for sorted dofs
! > @param [in] nboun2           number of transformed SPCs
! > @param [in] ndirboun2       (i) direction of transformed SPC i
! > @param [in] nodeboun2         (i) node of transformed SPC i
! > @param [in] xboun2            (i) value of transformed SPC i
! > @param [in] nmpc2              number of transformed mpcs
! > @param [in] ipompc2           (i) pointer to nodempc and coeffmpc for transformed MPC i
! > @param [in] nodempc2          nodes and directions of transformed MPCs
! > @param [in] coefmpc2          coefficients of transformed MPCs
! > @param [in] ikboun2           sorted dofs idof=8*(node-1)+dir for transformed SPCs
! > @param [in] ilboun2           transformed SPC numbers for sorted dofs
! > @param [in] ikmpc2               sorted dofs idof=8*(node-1)+dir for transformed MPCs
! > @param [in] ilmpc2              transformed SPC numbers for sorted dofs
! > @param [out] nslavspc       (2*i) pointer to islavspc...
! > @param [out] islavspc         ... which stores SPCs for slave node i
! > @param [out] nsspc            number of SPC for slave nodes
! > @param [out] nslavmpc       (2*i) pointer to islavmpc...
! > @param [out] islavmpc       ... which stores MPCs for slave node i
! > @param [out] nsmpc              number of MPC for slave nodes
! > @param [out] nslavspc2       (2*i) pointer to islavspc2...
! > @param [out] islavspc2       ... which stores transformed SPCs for slave node i
! > @param [out] nsspc2          number of transformed SPC for slave nodes
! > @param [out] nslavmpc2       (2*i) pointer to islavmpc2...
! > @param [out] islavmpc2       ... which stores transformed MPCs for slave node i
! > @param [out] nsmpc2              number of transformed MPC for slave nodes
! > @param [out] nmastspc       (2*i) pointer to imastspc...
! > @param [out] imastspc        ... which stores SPCs for master node i
! > @param [out] nmspc           number of SPC for master nodes
! > @param [out] nmastmpc       (2*i) pointer to imastmpc...
! > @param [out] imastmpc       ... which stores MPCs for master node i
! > @param [out] nmmpc              number of MPC for master nodes
! >
      subroutine conttiemortar(lakon,ipkon,kon,ntie,tieset,nset,set,&
           itiefac,islavsurf,islavnode,&
           imastnode,nslavnode,nmastnode,nslavs,&
           iponoels,inoels,&
           nk,&
           nboun,ndirboun,nodeboun,xboun,&
           nmpc,ipompc,nodempc,coefmpc,&
           ikboun,ilboun,ikmpc,ilmpc,&
           nboun2,ndirboun2,nodeboun2,xboun2,&
           nmpc2,ipompc2,nodempc2,coefmpc2,&
           ikboun2,ilboun2,ikmpc2,ilmpc2,&
           nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,&
           nslavspc2,islavspc2,nsspc2,nslavmpc2,islavmpc2,nsmpc2,&
           nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,&
           nmastmpc2,imastmpc2,nmmpc2)
      !      >
      !      >     subroutine to catalogue SPC'S and MPC'S for master and slave node
      !      >     needed for mortar contact
      !      >
      !      >     author: Sitzmann,Saskia
      !      >
      implicit none
      !
      logical nodeslavsurf,debug
      !
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,mastset,set(*)
      !
      integer ntie,i,j,k,l,nset,&
           ifaces,nelems,jfaces,ifacem,nelemm,nslavs,&
           jfacem,indexe,nopes,nopem,ipkon(*),kon(*),id,&
           node,&
           itiefac(2,*),islavsurf(2,*),islavnode(*),imastnode(*),&
           nslavnode(ntie+1),nmastnode(ntie+1),ifacecount,islav,imast,&
           ipos,index1,iponoels(*),inoels(2,*),ifreenoelold,&
           numbern,numberf,iface,kflag,nk,nmmpc2,&
           index,&
           nboun,ndirboun(*),nodeboun(*),&
           nmpc,ipompc(*),nodempc(3,*),dof,&
           nboun2,ndirboun2(*),nodeboun2(*),&
           nmpc2,ipompc2(*),nodempc2(3,*),&
           ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),&
           ikboun2(*),ilboun2(*),ikmpc2(*),ilmpc2(*),&
           nslavspc(2,*),islavspc(2,*),nsspc,nslavmpc(2,*),&
           islavmpc(2,*),&
           nsmpc,nslavspc2(2,*),islavspc2(2,*),&
           nsspc2,nslavmpc2(2,*),islavmpc2(2,*),&
           nsmpc2,nmastspc(2,*),imastspc(2,*),nmspc,nmastmpc(2,*),&
           imastmpc(2,*),nmmpc,isspc,imspc,ismpc,immpc,&
           nodem,nodes,ist,nmastmpc2(2,*),imastmpc2(2,*)
      !
      real*8  xboun(*),coefmpc(*),xboun2(*),coefmpc2(*)
      !
      debug=.false.
      !
      isspc=0
      ismpc=0
      i=0
      do i=1,ntie
         do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
            !     check for SPCs
            nslavspc(1,l)=isspc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikboun, dof,&
                    nboun, id)
               if(id>0)then
                  if(dof.eq.ikboun(id))then
                     isspc=isspc+1
                     islavspc(1,isspc)=ilboun(id)
                  endif
               endif
            enddo
            nslavspc(2,l)=isspc
            !     check for MPCs
            nslavmpc(1,l)=ismpc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikmpc, dof,&
                    nmpc, id)
               if(id>0)then
                  if(dof.eq.ikmpc(id))then
                     ismpc=ismpc+1
                     islavmpc(1,ismpc)=ipompc(ilmpc(id))
                  endif
               endif
            enddo
            nslavmpc(2,l)=ismpc
         enddo
      enddo
      nsspc=isspc
      nsmpc=ismpc
      !
      isspc=0
      ismpc=0
      i=0
      do i=1,ntie
         do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
            !     check for SPCs
            nslavspc2(1,l)=isspc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikboun2, dof,&
                    nboun2, id)
               if(id>0)then
                  if(dof.eq.ikboun2(id))then
                     isspc=isspc+1
                     islavspc2(1,isspc)=ilboun2(id)
                  endif
               endif
            enddo
            nslavspc2(2,l)=isspc
            !     check for MPCs
            nslavmpc2(1,l)=ismpc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikmpc2, dof,&
                    nmpc2, id)
               if(id>0)then
                  if(dof.eq.ikmpc2(id))then
                     ismpc=ismpc+1
                     islavmpc2(1,ismpc)=ipompc2(ilmpc2(id))
                  endif
               endif
            enddo
            nslavmpc2(2,l)=ismpc
         enddo
      enddo
      nsspc2=isspc
      nsmpc2=ismpc
      !
      imspc=0
      immpc=0
      do i=1,ntie
         do l=nmastnode(i)+1,nmastnode(i+1)
            node=imastnode(l)
            !     check for SPCs
            nmastspc(1,l)=imspc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikboun, dof,&
                    nboun, id)
               if(id>0)then
                  if(dof.eq.ikboun(id))then
                     imspc=imspc+1
                     imastspc(1,imspc)=ilboun(id)
                  endif
               endif
            enddo
            nmastspc(2,l)=imspc
            !     check for MPCs
            nmastmpc(1,l)=immpc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikmpc, dof,&
                    nmpc, id)
               if(id>0)then
                  if(dof.eq.ikmpc(id))then
                     immpc=immpc+1
                     imastmpc(1,immpc)=ipompc(ilmpc(id))
                  endif
               endif
            enddo
            nmastmpc(2,l)=immpc
         enddo
      enddo
      nmspc=imspc
      nmmpc=immpc

      immpc=0
      do i=1,ntie
         do l=nmastnode(i)+1,nmastnode(i+1)
            node=imastnode(l)
            !     check for MPCs
            nmastmpc2(1,l)=immpc
            do k=1,3
               dof=8*(node-1)+k
               call nident(ikmpc2, dof,&
                    nmpc2, id)
               if(id>0)then
                  if(dof.eq.ikmpc2(id))then
                     immpc=immpc+1
                     imastmpc2(1,immpc)=ipompc2(ilmpc2(id))
                  endif
               endif
            enddo
            nmastmpc2(2,l)=immpc
         enddo
      enddo
      nmmpc2=immpc
      !
      if(debug)then
         do i=1,ntie
            do l=nslavnode(i)+1,nslavnode(i+1)
               node=islavnode(l)
               write(*,*)'***node',node,'***'
               write(*,*) 'nodes-spc',node, nslavspc(1,l),nslavspc(2,l)
               do j=nslavspc(1,l)+1,nslavspc(2,l)
                  ist=islavspc(1,j)
                  write(*,*)' spc', ist, nodeboun(ist),ndirboun(ist)
               enddo
               write(*,*) 'nodes-mpc',node, nslavmpc(1,l),nslavmpc(2,l)
               do j=nslavmpc(1,l)+1,nslavmpc(2,l)
                  ist=islavmpc(1,j)
                  index=nodempc(3,ist)
                  write(*,*)' mpc',ist,nodempc(2,ist),nodempc(1,index)
               enddo
               write(*,*)'nodes-mpc2',node,nslavmpc2(1,l),nslavmpc2(2,l)
               do j=nslavmpc2(1,l)+1,nslavmpc2(2,l)
                  ist=islavmpc2(1,j)
                  index=nodempc2(3,ist)
                  write(*,*)' mpc',ist,nodempc2(2,ist),nodempc2(1,index)
               enddo
            enddo
         enddo
         do i=1,ntie
            do l=nmastnode(i)+1,nmastnode(i+1)
               node=imastnode(l)
               write(*,*)'***node',node,'***'
               write(*,*) 'nodem-spc',node, nmastspc(1,l),nmastspc(2,l)
               do j=nmastspc(1,l)+1,nmastspc(2,l)
                  ist=imastspc(1,j)
                  write(*,*)'spc', ist, nodeboun(ist),ndirboun(ist)
               enddo
               write(*,*) 'nodem-mpc',node, nmastmpc(1,l),nmastmpc(2,l)
               do j=nmastmpc(1,l)+1,nmastmpc(2,l)
                  ist=imastmpc(1,j)
                  index=nodempc(3,ist)
                  write(*,*)'mpc', ist, nodempc(2,ist), nodempc(1,index)
               enddo
               write(*,*) 'nodem-mpc2',node, nmastmpc2(1,l),&
                nmastmpc2(2,l)
               do j=nmastmpc2(1,l)+1,nmastmpc2(2,l)
                  ist=imastmpc2(1,j)
                  index=nodempc2(3,ist)
                  write(*,*)'mpc', ist, nodempc2(2,ist),&
                   nodempc2(1,index)
               enddo
            enddo
         enddo
      endif
      !
      return
      end
      
