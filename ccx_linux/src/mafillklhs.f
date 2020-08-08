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
      subroutine mafillklhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
     &  nactdok,icolk,jqk,irowk,neqk,nzlk,
     &  ikmpc,ilmpc,ikboun,ilboun,nzsk,adbk,aubk,ipvar,var)
!
!     filling the lhs turbulence matrix in spare matrix format (sm)
!
!     it is assumed that the temperature MPC's also apply to the
!     turbulence. Temperature MPC's are not allowed to contain
!     other variables than the temperature.
!
      implicit none
!
      character*8 lakon(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  icolk(*),jqk(*),ikmpc(*),nzsk,nmethod,
     &  ilmpc(*),ikboun(*),ilboun(*),nactdok(*),konl(20),irowk(*),
     &  ipkon(*),ipvar(*)
!
      integer nk,ne,nboun,nmpc,neqk,nzlk,i,j,jj,
     &  ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &  mpc1,mpc2,index1,index2,node1,node2,
     &  indexe,nope,i0
!
      real*8 co(3,*),xboun(*),coefmpc(*),sm(78,78),adbk(*),aubk(*),
     &  var(*)
!
      real*8 value
!
      i0=0
!
!     determining nzlk
!
      nzlk=0
      do i=neqk,1,-1
         if(icolk(i).gt.0) then
            nzlk=i
            exit
         endif
      enddo
!
      do i=1,neqk
         adbk(i)=0.d0
      enddo
      do i=1,nzsk
         aubk(i)=0.d0
      enddo
!
!     loop over all fluid elements
!
      do i=1,ne
!
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'F') cycle
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'2') then
           nope=20
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
        else
           cycle
        endif
!
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
!
!       the temperature element routine = the turbulence element
!       routine
!
        call e_c3d_tlhs(co,nk,konl,lakon(i),sm,i,ipvar,var)
!
        do jj=1,nope
!
          node1=kon(indexe+jj)
          jdof1=nactdok(node1)
!
          do ll=jj,nope
!
            node2=kon(indexe+ll)
            jdof2=nactdok(node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.gt.0).and.(jdof2.gt.0)) then
                  call add_sm_fl(aubk,adbk,jqk,irowk,jdof1,jdof2,
     &                 sm(jj,ll),jj,ll)
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!
!              idof1: genuine DOF
!              idof2: nominal DOF of the SPC/MPC
!
               if(jdof1.le.0) then
                  idof1=jdof2
                  idof2=jdof1
               else
                  idof1=jdof1
                  idof2=jdof2
               endif
               if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
!
!                    regular DOF / MPC
!
                     id=(-idof2+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdok(nodempc(1,index))
                        if(idof2.gt.0) then
                              value=-coefmpc(index)*sm(jj,ll)/
     &                               coefmpc(ist)
                              if(idof1.eq.idof2) value=2.d0*value
                              call add_sm_fl(aubk,adbk,jqk,irowk,
     &                             idof1,idof2,value,i0,i0)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
            else
               idof1=jdof1
               idof2=jdof2
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
!
!                    MPC id1 / MPC id1
!
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdok(nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdok(nodempc(1,index2))
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_fl(aubk,adbk,jqk,
     &                             irowk,idof1,idof2,value,i0,i0)
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
!
!                    MPC id1 / MPC id2
!
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdok(nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdok(nodempc(1,index2))
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                                 if(idof1.eq.idof2) value=2.d0*value
                                 call add_sm_fl(aubk,adbk,jqk,
     &                             irowk,idof1,idof2,value,i0,i0)
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
        enddo
      enddo
!
      return
      end
