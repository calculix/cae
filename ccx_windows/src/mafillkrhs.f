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
      subroutine mafillkrhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
     &     nface,nactdok,neqk,nmethod,ikmpc,ilmpc,
     &     ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,vcon,nzsk,
     &     dtime,matname,mi,ncmat_,shcon,nshcon,theta1,
     &     bk,bt,vcontu,isolidsurf,nsolidsurf,ifreestream,nfreestream,
     &     xsolidsurf,yy,compressible,turbulent,ithermal,ipvar,var,
     &     ipvarf,varf,nea,neb,dt,ck,ct)
!     
!     filling the rhs b of the turbulence equations (step 5)
!     
!     it is assumed that the temperature MPC's also apply to the
!     turbulence. The temperature MPC's are not allowed to contain
!     any other variables but temperatures
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakon(*)
      character*80 matname(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nelemface(*),ikmpc(*),ilmpc(*),ikboun(*),compressible,
     &     ilboun(*),nactdok(*),konl(20),nrhcon(*),mi(*),
     &     ipkon(*),nshcon(*),ifreestream(*),nfreestream,isolidsurf(*),
     &     nsolidsurf,turbulent,ithermal(*),ipvar(*),ipvarf(*),
     &     ielmat(mi(3),*)
!     
      integer nk,ne,nboun,nmpc,nface,neqk,nmethod,nzsk,i,j,k,jj,
     &     id,ist,index,jdof1,idof1,node1,kflag,ntmat_,indexe,nope,
     &     i0,ncmat_,nea,neb
!     
      real*8 co(3,*),xboun(*),coefmpc(*),bk(*),ck(*),ct(*),
     &     vold(0:mi(2),*),var(*),varf(*),dt(*),ggk(78),ggt(78),
     &     vcon(0:4,*),ffk(78),rhcon(0:1,ntmat_,*),yy(*),
     &     shcon(0:3,ntmat_,*),theta1,bt(*),fft(78),vcontu(2,*),
     &     xsolidsurf(*)
!     
      real*8 dtime
!     
      kflag=2
      i0=0
!     
      do i=1,neqk
        bk(i)=0.d0
        bt(i)=0.d0
        ck(i)=0.d0
        ct(i)=0.d0
      enddo
!     
      do i=nea,neb
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
        call e_c3d_krhs(co,nk,kon(indexe+1),lakon(i),ffk,fft,i,nmethod,
     &       rhcon,
     &       nrhcon,ielmat,ntmat_,vold,vcon,dtime,matname,mi(1),
     &       shcon,nshcon,vcontu,compressible,yy,nelemface,sideface,
     &       nface,turbulent,ithermal,ipvar,var,ipvarf,varf,dt,ggk,ggt)
!     
        do jj=1,nope
!     
          j=jj
          k=jj-3*(j-1)
!     
          node1=kon(indexe+j)
          jdof1=nactdok(node1)
!     
!     inclusion of ffk and fft
!     
          if(jdof1.le.0) then
            if(nmpc.ne.0) then
              if(jdof1.ne.2*(jdof1/2)) then
                id=(-jdof1+1)/2
                ist=ipompc(id)
                index=nodempc(3,ist)
                if(index.eq.0) cycle
                do
                  jdof1=nactdok(nodempc(1,index))
                  if(jdof1.gt.0) then
                    bk(jdof1)=bk(jdof1)
     &                   -coefmpc(index)*ffk(jj)
     &                   /coefmpc(ist)
                    bt(jdof1)=bt(jdof1)
     &                   -coefmpc(index)*fft(jj)
     &                   /coefmpc(ist)
                    ck(jdof1)=ck(jdof1)
     &                   -coefmpc(index)*ggk(jj)
     &                   /coefmpc(ist)
                    ct(jdof1)=ct(jdof1)
     &                   -coefmpc(index)*ggt(jj)
     &                   /coefmpc(ist)
                  endif
                  index=nodempc(3,index)
                  if(index.eq.0) exit
                enddo
              endif
            endif
            cycle
          endif
          bk(jdof1)=bk(jdof1)+ffk(jj)
          bt(jdof1)=bt(jdof1)+fft(jj)
          ck(jdof1)=ck(jdof1)+ggk(jj)
          ct(jdof1)=ct(jdof1)+ggt(jj)
!     
        enddo
      enddo
!     
      return
      end
