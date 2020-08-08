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
      subroutine mafillv2rhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,ipompc,nodempc,coefmpc,nmpc,b,nactdoh,
     &     icolv,jqv,irowv,neqv,nzlv,nmethod,ikmpc,ilmpc,ikboun,
     &     ilboun,vold,nzsv,dt,v,theta2,iexplicit,nea,neb,mi,dtimef,
     &     ipvar,var,ipvarf,varf)
!     
!     filling the rhs b of the velocity equations (step 3)
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     icolv(*),jqv(*),ikmpc(*),ilmpc(*),ikboun(*),ilboun(*),
     &     nactdoh(0:4,*),konl(20),irowv(*),ipkon(*),nea,neb,mi(*),
     &     ipvar(*),ipvarf(*)
!     
      integer nk,ne,nboun,nmpc,neqv,nzlv,nmethod,nzsv,i,j,k,jj,
     &     id,ist,index,jdof1,iexplicit,node1,kflag,indexe,nope,i0
!     
      real*8 co(3,*),xboun(*),coefmpc(*),b(*),v(0:mi(2),*),theta2,
     &     vold(0:mi(2),*),ff(78),dtimef,var(*),varf(*),dt(*)
!     
      kflag=2
      i0=0
!     
      do i=1,neqv
        b(i)=0.d0
      enddo
!     
      do i=nea,neb
!     
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'F') cycle
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'4') then
          nope=4
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
        call e_c3d_v2rhs(co,nk,kon(indexe+1),lakon(i),
     &       ff,i,nmethod,vold,v,dtimef,theta2,iexplicit,mi,
     &       ipvar,var,ipvarf,varf,dt)
!     
        do jj=1,3*nope
!     
          j=(jj-1)/3+1
          k=jj-3*(j-1)
!     
          node1=kon(indexe+j)
          jdof1=nactdoh(k,node1)
!     
!     distributed forces
!     
          if(jdof1.le.0) then
            if(nmpc.ne.0) then
              if(jdof1.ne.2*(jdof1/2)) then
                id=(-jdof1+1)/2
                ist=ipompc(id)
                index=nodempc(3,ist)
                if(index.eq.0) cycle
                do
                  jdof1=nactdoh(nodempc(2,index),
     &                 nodempc(1,index))
                  if(jdof1.gt.0) then
                    b(jdof1)=b(jdof1)
     &                   -coefmpc(index)*ff(jj)
     &                   /coefmpc(ist)
                  endif
                  index=nodempc(3,index)
                  if(index.eq.0) exit
                enddo
              endif
            endif
            cycle
          endif
          b(jdof1)=b(jdof1)+ff(jj)
!     
        enddo
      enddo
!     
      return
      end
