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
      subroutine applyboun(ifaext,nfaext,ielfa,ikboun,ilboun,
     &  nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
     &  ifabou,nfabou,nface,nodeboun,ndirboun,ikmpc,ilmpc,labmpc,nmpc,
     &  nactdohinv,compressible,iatleastonepressurebc,ipkonf,kon,konf,
     &  inlet)
!
!     stores pointers to ifabou in ielfa(2,*) at those locations
!     which are zero (external faces)
!     stores pointers to the boundary conditions in field ifabou
!
      implicit none
!
      character*1 typeboun(*)
      character*20 sideload(*),labmpc(*)
!
      integer nfabou,ifaext(*),nfaext,ielem,ielfa(4,*),ifa,jface,
     &  j,idof,ikboun(*),nboun,id,ilboun(*),iboun,nelemload(2,*),
     &  nload,isolidsurf(*),nsolidsurf,ifabou(*),i,nface,indexb,
     &  nodeboun(*),ndirboun(*),jsum,ig,ikmpc(*),ilmpc(*),nmpc,mpc,
     &  nactdohinv(*),compressible,iatleastonepressurebc,iface,
     &  ipkonf(*),kon(*),konf(*),indexe,inlet(*)
!
      nfabou=1
      iatleastonepressurebc=0
!
      do i=1,nfaext
!
!        number of the face in field ielfa
!
         ifa=ifaext(i)
!
!        adjacent element number (global number)
!
         ielem=nactdohinv(ielfa(1,ifa))
!
!        face label used to apply the SPC
!
         iface=ielfa(4,ifa)
!
         jface=10*ielem+iface
!
!        SPC's: loop over the degrees of freedom
!
         jsum=0
         do j=0,4
            idof=-(8*(jface-1)+j)
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.idof) then
                  iboun=ilboun(id)
                  if(typeboun(iboun).ne.'F') cycle
                  if(ielfa(2,ifa).eq.0) then
                     ielfa(2,ifa)=-nfabou
c                     write(*,*) 'applyboun bc',ielfa(2,ifa)
                     nfabou=nfabou+7
                  endif
!
!                 if all velocity components are known no pressure
!                 should be defined (only for incompressible fluids)
!
                  if(compressible.eq.0) then
                     if(j.eq.4) then
                        if(jsum.eq.6) then
                           write(*,*) '*WARNING in applyboun: a pressure
     & SPC is being applied to'
                           write(*,*) '         face ',iface,
     &                      'of element ',ielem,'for which all'
                           write(*,*) '         velocities are known (by
     & SPCs or MPCs). The pressure'
                           write(*,*) '         SPC is discarded'
                           write(*,*)
                           exit
                        endif
                        iatleastonepressurebc=1
                     endif
                  endif
                  jsum=jsum+j
c                  write(*,*) 'applyboun jsum ',ifa,jsum
!
                  ifabou(-ielfa(2,ifa)+j)=iboun
               endif
            endif
!
!           MPC's: loop over the degrees of freedom
!
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  mpc=ilmpc(id)
                  if(labmpc(mpc)(1:5).ne.'FLUID') cycle
                  if(ielfa(2,ifa).eq.0) then
                     ielfa(2,ifa)=-nfabou
                     nfabou=nfabou+7
                  else if(ifabou(-ielfa(2,ifa)+j).ne.0) then
                     write(*,*) '*ERROR in applyboun: MPC is applied'
                     write(*,*) '       to degree of freedom ',j
                     write(*,*) '       in face ',iface
                     write(*,*) '       of element ',ielem,'.'
                     write(*,*) '       To this degree of freedom '
                     write(*,*) '       another SPC or MPC has already'
                     write(*,*) '       been applied'
                     call exit(201)
                  endif
!
!                 if all velocity components are known no pressure
!                 should be defined (only for incompressible fluids)
!
                  if(compressible.eq.0) then
                     if(j.eq.4) then
                        if(jsum.eq.6) then
                           write(*,*) '*WARNING in applyboun: a pressure
     & MPC is being applied to'
                           write(*,*) '         face ',iface,
     &                       'of element ',ielem,'for which all'
                           write(*,*) '         velocities are known (by
     & SPCs or MPCs). The pressure'
                           write(*,*) '         MPC is discarded'
                           exit
                        endif
                     endif
                  endif
                  jsum=jsum+j
!
                  ifabou(-ielfa(2,ifa)+j)=-mpc
               endif
            endif
         enddo
!
!        heat flux
!
         call nident2(nelemload,ielem,nload,id)
!
         do
            if(id.gt.0) then
               if(nelemload(1,id).eq.ielem) then
                  if(sideload(id)(1:1).eq.'S') then
                     read(sideload(id)(2:2),'(i1)') ig
                     if(ig.eq.iface) then
                        if(ielfa(2,ifa).eq.0) then
                           ielfa(2,ifa)=-nfabou
                           nfabou=nfabou+7
                        endif
                        ifabou(-ielfa(2,ifa)+6)=id
                     endif
                  endif
                  id=id-1
                  cycle
               else
                  exit
               endif
            else
               exit
            endif
         enddo
!
!        sliding conditions
!
         call nident2(nelemload,ielem,nload,id)
!
         do
            if(id.gt.0) then
c               write(*,*) 'applyboun ',nelemload(1,id),sideload(id)
               if(nelemload(1,id).eq.ielem) then
                  if(sideload(id)(1:1).eq.'M') then
                     read(sideload(id)(2:2),'(i1)') ig
                     if(ig.eq.iface) then
c                        write(*,*) 'store '
                        if(ielfa(2,ifa).eq.0) then
                           ielfa(2,ifa)=-nfabou
                           nfabou=nfabou+7
                        endif
c                        ifabou(-ielfa(2,ifa)+5)=2
                        ifabou(-ielfa(2,ifa)+5)=-1
                     endif
                  endif
                  id=id-1
                  cycle
               else
                  exit
               endif
            else
               exit
            endif
         enddo
!
!        wall
!
         call nident(isolidsurf,jface,nsolidsurf,id)
         if(id.gt.0) then
            if(isolidsurf(id).eq.jface) then
               if(ielfa(2,ifa).eq.0) then
                  ielfa(2,ifa)=-nfabou
                  nfabou=nfabou+7
               endif
               indexb=-ielfa(2,ifa)
               if((ifabou(indexb+1).eq.0).or.
     &            (ifabou(indexb+2).eq.0).or.
     &            (ifabou(indexb+3).eq.0)) then
                  write(*,*) '*ERROR in applyboun: face',iface
                  write(*,*) '       of element ',ielem,'is defined'
                  write(*,*) '       as solid surface but not all'
                  write(*,*) '       velocity components are defined'
                  write(*,*) '       as boundary conditions'
                  call exit(201)
               endif
c               ifabou(indexb+5)=1
               ifabou(indexb+5)=id
            endif
         endif
!
!        checking for:
!           - absent boundary conditions or
!           - absent velocity boundary conditions without 
!             sliding conditions
!        -> zero velocity gradient and zero temperature gradient
!           
!        tagged by negative ielfa(3,ifa) (more than one layer) or
!               by zero ielfa(3,ifa) (one layer)
!
         if(ielfa(2,ifa).eq.0) then
!
!           no boundary conditions
!
            ielfa(3,ifa)=-ielfa(3,ifa)
         elseif(ifabou(-ielfa(2,ifa)+5).ne.-1) then
c            write(*,*) 'applyboun absent ',ifa,ifabou(-ielfa(2,ifa)+5),
c     &          jsum
!
!           no sliding conditions
!            
            if((jsum.eq.0).or.(jsum.eq.4)) then
!
!              no velocity boundary conditions
!
               ielfa(3,ifa)=-ielfa(3,ifa)
            endif
!
!           check for inlet conditions (all velocity components given
!           and no wall)
!
            indexb=-ielfa(2,ifa)
            if((ifabou(indexb+1).ne.0).and.
     &           (ifabou(indexb+2).ne.0).and.
     &           (ifabou(indexb+3).ne.0).and.
     &           (ifabou(indexb+5).le.0)) then
               inlet(ifa)=1
            endif
         endif
      enddo
!
!     dimension of field ifabou containing the pointers to the
!     boundary conditions
!
      nfabou=nfabou-1
!
c      write(*,*)
c      do i=1,nface
c         write(*,*) 'applyboun ielfa ',i,(ielfa(j,i),j=1,4)
c      enddo
c      do i=1,nfabou
c         write(*,*) 'applyboun ifabou',i,ifabou(i)
c      enddo
c      do i=1,nboun
c         write(*,*) 'applyboun nodeboun',i,nodeboun(i),ndirboun(i)
c      enddo
c      write(*,*)
!     
      return
      end
