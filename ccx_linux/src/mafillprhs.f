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
      subroutine mafillprhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
     &     nface,b,nactdoh,icolp,jqp,irowp,neqp,nzlp,nmethod,ikmpc,
     &     ilmpc,ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,vcon,
     &     nzsp,dt,matname,mi,ncmat_,shcon,nshcon,v,theta1,
     &     iexplicit,physcon,nea,neb,dtimef,ipvar,var,ipvarf,varf)
!     
!     filling the rhs b of the pressure equations (step 2)
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakon(*)
      character*80 matname(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nelemface(*),icolp(*),jqp(*),ikmpc(*),ilmpc(*),ikboun(*),
     &     ilboun(*),nactdoh(0:4,*),konl(20),irowp(*),nrhcon(*),
     &     mi(*),ielmat(mi(3),*),
     &     ipkon(*),nshcon(*),iexplicit,nea,neb,ipvar(*),ipvarf(*),
     &     nk,ne,nboun,nmpc,nface,neqp,nzlp,nmethod,nzsp,i,j,k,l,jj,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,kflag,ntmat_,indexe,nope,
     &     i0,ncmat_,idof3
!     
      real*8 co(3,*),xboun(*),coefmpc(*),b(*),v(0:mi(2),*),
     &     vold(0:mi(2),*),dt(*),
     &     vcon(0:4,*),ff(78),sm(78,78),rhcon(0:1,ntmat_,*),
     &     shcon(0:3,ntmat_,*),theta1,physcon(*),var(*),varf(*)
!     
      real*8 value,dtimef
!     
      kflag=2
      i0=0
!     
      do i=1,neqp
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
        call e_c3d_prhs(co,nk,kon(indexe+1),lakon(i),sm,ff,i,nmethod,
     &       rhcon,
     &       nrhcon,ielmat,ntmat_,v,vold,vcon,nelemface,sideface,
     &       nface,dtimef,matname,mi(1),shcon,nshcon,theta1,physcon,
     &       iexplicit,ipvar,var,ipvarf,varf,dt)
!     
        do jj=1,nope
!     
          j=jj
          k=jj-3*(j-1)
!     
          node1=kon(indexe+j)
          jdof1=nactdoh(4,node1)
!     
!     inclusion of ff
!     
          if(jdof1.le.0) then
            if(nmpc.ne.0) then
              if(jdof1.ne.2*(jdof1/2)) then
                id=(-jdof1+1)/2
                ist=ipompc(id)
                index=nodempc(3,ist)
                if(index.eq.0) cycle
                do
                  jdof1=nactdoh(4,nodempc(1,index))
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
      
