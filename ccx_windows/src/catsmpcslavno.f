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
!     Author: Saskia Sitzmann
!     
!     [in] nboun2           number of transformed SPCs
!     [in] nmpc2		number of transformed mpcs
!     [in] ipompc2           (i) pointer to nodempc and coeffmpc for transformed MPC i
!     [in] nodempc2          nodes and directions of transformed MPCs
!     [in] ikboun2           sorted dofs idof=8*(node-1)+dir for transformed SPCs
!     [in] ilboun2           transformed SPC numbers for sorted dofs
!     [in] ikmpc2 		sorted dofs idof=8*(node-1)+dir for transformed MPCs
!     [in] ilmpc2		transformed SPC numbers for sorted dofs
!     [out] nslavspc	(2*i) pointer to islavspc...
!     [out] islavspc         ... which stores SPCs for slave node i
!     [out] nsspc            number of SPC for slave nodes
!     [out] nslavmpc	(2*i) pointer to islavmpc...
!     [out] islavmpc	... which stores MPCs for slave node i
!     [out] nsmpc		number of MPC for slave nodes
!     [out] nslavspc2	(2*i) pointer to islavspc2...
!     [out] islavspc2       ... which stores transformed SPCs for slave node i
!     [out] nsspc2          number of transformed SPC for slave nodes
!     [out] nslavmpc2	(2*i) pointer to islavmpc2...
!     [out] islavmpc2	... which stores transformed MPCs for slave node i
!     [out] nsmpc2		number of transformed MPC for slave nodes 
!     [out] nmastspc	(2*i) pointer to imastspc...
!     [out] imastspc        ... which stores SPCs for master node i
!     [out] nmspc           number of SPC for master nodes
!     [out] nmastmpc	(2*i) pointer to imastmpc...
!     [out] imastmpc	... which stores MPCs for master node i
!     [out] nmmpc		number of MPC for master nodes
!     
      subroutine catsmpcslavno(ntie,islavnode,imastnode,
     &     nslavnode,nmastnode,nboun,ndirboun,nodeboun,nmpc,ipompc,
     &     nodempc,ikboun,ilboun,ikmpc,ilmpc,nboun2,nmpc2,ipompc2,
     &     nodempc2,ikboun2,ilboun2,ikmpc2,ilmpc2,
     &     nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
     &     nslavspc2,islavspc2,nsspc2,nslavmpc2,islavmpc2,nsmpc2,
     &     nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
     &     nmastmpc2,imastmpc2,nmmpc2)
!     
!     subroutine to catalogue SPC'S and MPC'S for master and slave node    
!     needed for mortar contact
!     
!     correlate the transformed and untransformed spc's and mpc's to the
!     slave and master nodes
!     
!     author: Sitzmann,Saskia
!     
      implicit none
!     
      logical debug
!     
      integer ntie,i,j,k,l,
     &     id,node,islavnode(*),imastnode(*),nslavnode(ntie+1),
     &     nmastnode(ntie+1),nmmpc2,index,nboun,ndirboun(*),nodeboun(*),
     &     nmpc,ipompc(*),nodempc(3,*),dof,nboun2,nmpc2,ipompc2(*),
     &     nodempc2(3,*),ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &     ikboun2(*),ilboun2(*),ikmpc2(*),ilmpc2(*),nslavspc(2,*),
     &     islavspc(2,*),nsspc,nslavmpc(2,*),islavmpc(2,*),nsmpc,
     &     nslavspc2(2,*),islavspc2(2,*),nsspc2,nslavmpc2(2,*),
     &     islavmpc2(2,*),nsmpc2,nmastspc(2,*),imastspc(2,*),nmspc,
     &     nmastmpc(2,*),imastmpc(2,*),nmmpc,isspc,imspc,ismpc,immpc,
     &     ist,nmastmpc2(2,*),imastmpc2(2,*)
!     
      debug=.false.
!     
!     slave surfaces
!     
      isspc=0
      ismpc=0
      i=0
      do i=1,ntie
        do l=nslavnode(i)+1,nslavnode(i+1)
          node=islavnode(l)
!     
!     check for SPCs
!     
          nslavspc(1,l)=isspc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikboun,dof,nboun,id)
            if(id>0)then
              if(dof.eq.ikboun(id))then
                isspc=isspc+1
                islavspc(1,isspc)=ilboun(id)
              endif
            endif
          enddo
          nslavspc(2,l)=isspc
!     
!     check for MPCs
!     
          nslavmpc(1,l)=ismpc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikmpc,dof,nmpc,id)
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
!     
!     check for SPCs (2)
!     
          nslavspc2(1,l)=isspc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikboun2,dof,nboun2,id)
            if(id>0)then
              if(dof.eq.ikboun2(id))then
                isspc=isspc+1
                islavspc2(1,isspc)=ilboun2(id)
              endif
            endif
          enddo
          nslavspc2(2,l)=isspc
!     
!     check for MPCs (2)
!     
          nslavmpc2(1,l)=ismpc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikmpc2,dof,nmpc2,id)
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
!     master surfaces
!     
      imspc=0
      immpc=0
      do i=1,ntie
        do l=nmastnode(i)+1,nmastnode(i+1)
          node=imastnode(l)
!     
!     check for SPCs
!     
          nmastspc(1,l)=imspc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikboun,dof,nboun,id)
            if(id>0)then
              if(dof.eq.ikboun(id))then
                imspc=imspc+1
                imastspc(1,imspc)=ilboun(id)
              endif
            endif
          enddo
          nmastspc(2,l)=imspc
!     
!     check for MPCs
!     
          nmastmpc(1,l)=immpc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikmpc,dof,nmpc,id)
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
!     
      immpc=0
      do i=1,ntie
        do l=nmastnode(i)+1,nmastnode(i+1)
          node=imastnode(l)
!     
!     check for MPCs (2)
!     
          nmastmpc2(1,l)=immpc
          do k=1,3
            dof=8*(node-1)+k
            call nident(ikmpc2,dof,nmpc2,id)
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
            write(*,*) 'nodes-spc',node,nslavspc(1,l),nslavspc(2,l)
            do j=nslavspc(1,l)+1,nslavspc(2,l)
              ist=islavspc(1,j)
              write(*,*)' spc',ist,nodeboun(ist),ndirboun(ist)
            enddo
            write(*,*) 'nodes-mpc',node,nslavmpc(1,l),nslavmpc(2,l)
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
            write(*,*) 'nodem-spc',node,nmastspc(1,l),nmastspc(2,l)
            do j=nmastspc(1,l)+1,nmastspc(2,l)
              ist=imastspc(1,j)
              write(*,*)'spc',ist,nodeboun(ist),ndirboun(ist)
            enddo
            write(*,*) 'nodem-mpc',node,nmastmpc(1,l),nmastmpc(2,l)
            do j=nmastmpc(1,l)+1,nmastmpc(2,l)
              ist=imastmpc(1,j)
              index=nodempc(3,ist)
              write(*,*)'mpc',ist,nodempc(2,ist),nodempc(1,index)
            enddo
            write(*,*) 'nodem-mpc2',node,nmastmpc2(1,l),
     &           nmastmpc2(2,l)
            do j=nmastmpc2(1,l)+1,nmastmpc2(2,l)
              ist=imastmpc2(1,j)
              index=nodempc2(3,ist)
              write(*,*)'mpc',ist,nodempc2(2,ist),
     &             nodempc2(1,index)
            enddo
          enddo
        enddo
      endif
!     
      return
      end
      
