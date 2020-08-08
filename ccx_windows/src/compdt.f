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
      subroutine compdt(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,
     &  iponoel,inoel,dtimef,iexplicit,ielmat,physcon,dh,cocon,
     &  ncocon,ithermal,mi,ipkon,kon,lakon,ne,v,co,turbulent,vcontu,
     &  vcon)
!
!     - determine the time step for each node (stored in field dt
!       and the minimum value across all nodes (dtimef)
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer nk,i,j,k,iponoel(*),inoel(3,*),index,nelem,ithermal(*),
     &  mi(*),
     &  nshcon(*),nrhcon(*),ntmat_,ielmat(mi(3),*),imat,ncocon(2,*),
     &  ipkon(*),kon(*),ne,nope,indexe,iflag,iexplicit,turbulent
!
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),vold(0:mi(2),*),temp,vel,dtcon,dtmed,
     &  physcon(*),dh(*),cocon(0:6,ntmat_,*),dtthd,cond,voldl(3,20),
     &  xl(3,20),vertex6(3,6),vertex8(3,8),xi,et,ze,xsj,shp(4,20),
     &  h,v(0:mi(2),*),co(3,*),dd,vcontu(2,*),dttud,vcon(0:4,*)
!
      data vertex6 /0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     &              0.d0,1.d0,0.d0,0.d0,0.d0,1.d0,
     &              1.d0,0.d0,1.d0,0.d0,1.d0,1.d0/
      data vertex8 /-1.d0,-1.d0,-1.d0,1.d0,-1.d0,-1.d0,
     &              1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,
     &              -1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &              1.d0,1.d0,1.d0,-1.d0,1.d0,1.d0/
      data iflag /3/
c!
c!     determining the element height in flow direction
c!
c      if(iexplicit.eq.1) then
c         do i=1,ne
c            indexe=ipkon(i)
c            if(indexe.lt.0) cycle
c            lakonl(1:8)=lakon(i)(1:8)
c!     
c!     number of nodes in the element
c!     
c            if(lakonl(4:4).eq.'2') then
c               nope=20
c            elseif(lakonl(4:4).eq.'8') then
c               nope=8
c            elseif(lakonl(4:5).eq.'10') then
c               nope=10
c            elseif(lakonl(4:4).eq.'4') then
c               nope=4
c            elseif(lakonl(4:5).eq.'15') then
c               nope=15
c            elseif(lakonl(4:4).eq.'6') then
c               nope=6
c            else
c               cycle
c            endif
c!     
c!     velocity at the nodes
c!     
c            do j=1,nope
c               do k=1,3
c                  voldl(k,j)=vold(k,kon(indexe+j))
c                  xl(k,j)=co(k,kon(indexe+j))
c               enddo
c            enddo
c!     
c!     element height
c!     
c            h=0.d0
c            do j=1,nope
c               if(nope.eq.20) then
c                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.8) then
c                  xi=vertex8(1,j)
c                  et=vertex8(2,j)
c                  ze=vertex8(3,j)
c                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.10) then
c                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.4) then
c                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.15) then
c                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.6) then
c                  xi=vertex6(1,j)
c                  et=vertex6(2,j)
c                  ze=vertex6(3,j)
c                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
c               endif
c!     
c               dd=dsqrt(voldl(1,j)*voldl(1,j)+
c     &              voldl(2,j)*voldl(2,j)+voldl(3,j)*voldl(3,j))
c               if(dd.lt.1.d-10) then
c                  cycle
c               else
c                  h=h+dabs(shp(1,j)*voldl(1,j)+shp(2,j)*voldl(2,j)+
c     &              shp(3,j)*voldl(3,j))/dd
c               endif
c            enddo
c!     
cc            if(h.gt.0.d0) h=2.d0/h
c            if(h.gt.0.d0) h=nope/h
c!     
c!        height at the nodes of the elements is replaced by the
c!        element height of the latter is smaller
c!
c            do j=1,nope
c               if(dtl(kon(indexe+j)).gt.h) dtl(kon(indexe+j))=h
c            enddo
c         enddo
c      endif
!     
!     determining the time increment dt for each node.
!
!     edge nodes (fields iponoel and inoel are determined in precfd.f)
!
      dtimef=1.d30
!
      do i=1,nk
         index=iponoel(i)
         if(index.le.0) cycle
!
!        look at an element belonging to the edge node
!
         nelem=inoel(1,index)
!     
!        determining the time increment
!
         imat=ielmat(1,nelem)
         temp=vold(0,i)
!
         vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
!
         if(iexplicit.eq.1) then
!
!           gas
!
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            r=shcon(3,1,imat)
            rho=vcon(4,i)
!
!           convective time step (dtcon)
!
            dtcon=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
         else
!
!           liquid
!
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
!
!           convective time step (dtcon)
!
            if(vel.lt.1.d-10) vel=1.d-10
            dtcon=dh(i)/vel
         endif
!
!        mechanical diffusion time step (dtmed)
!     
         call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
         if(dvi.lt.1.d-10) dvi=1.d-10
!
         dtmed=dh(i)*dh(i)*rho/(2.d0*dvi)
!     
         dt(i)=dtcon*dtmed/(dtcon+dtmed)
!     
!        thermal diffusion time step (dtthd)
!     
         if(ithermal(1).gt.1) then
            call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &           cond)
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            if(cond.lt.1.d-10) cond=1.d-10
            dtthd=dh(i)*dh(i)*rho*cp/(2.d0*cond)
            dt(i)=(dt(i)*dtthd)/(dt(i)+dtthd)
         endif
!     
!        turbulent diffusion time step (dttud)
!     
         if(turbulent.ne.0) then
            dttud=1.d0/(1.d0+0.1656d0*dabs(vcontu(2,i)))
            dt(i)=(dt(i)*dttud)/(dt(i)+dttud)
         endif
!
!        safety factor for compressible fluids
!
         if(iexplicit.eq.1) dt(i)=dt(i)/1.25d0
!
         if(dt(i).lt.dtimef) dtimef=dt(i)
!
      enddo
!
!     increased damping for incompressible fluids
!     the use of an internal (for the damping) and an external
!     (for the time derivative) time step stems from Zienkiewicz,
!     Taylor and Nithiarasu, The Finite Element Method for Fluid
!     Dynamics, 6th edition, p94 bottom and p 95 top.
!
      if(iexplicit.eq.1) then
         do i=1,nk
            dt(i)=dtimef
         enddo
       endif
!
      return
      end
