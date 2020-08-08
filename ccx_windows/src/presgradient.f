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
      subroutine presgradient(iponoel,inoel,sa,sav,nk,dt,shockcoef,
     &  dtimef,ipkon,kon,lakon,vold,mi,compressible,nmethod,dtl,
     &  isolidsurf,nsolidsurf,co,euler)
!
!     determining measure for the pressure gradient
!
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponoel(*),inoel(3,*),nk,i,j,index,indexe,nope,
     &  ipkon(*),kon(*),node,ielem,mi(*),compressible,nmethod,
     &  isolidsurf(*),nsolidsurf,id,isum,euler
!
      real*8 sa(*),sav(*),dt(*),shockcoef,dtimef,ca,sum,xmaxsum,pa,
     &  vold(0:mi(2),*),dtl(*),co(3,*),dd,v(3),p(3),cosang,
     &  xmaxshear,sumabs
!
c      if(euler.eq.1) then
!
!        nonviscous (euler): pressure switch is 1 if local
!                            maximum or minimum occurs, independent
!                            of its size
!
         do i=1,nk
            if(iponoel(i).le.0) cycle
            sum=0.d0
            sumabs=0.d0
            pa=vold(4,i)
            index=iponoel(i)
            do
               ielem=inoel(1,index)
               if(ipkon(ielem).lt.0) cycle
               if(lakon(ielem)(1:1).ne.'F') cycle
               if(lakon(ielem)(4:4).eq.'8') then
                  nope=8
               elseif(lakon(ielem)(4:4).eq.'4') then
                  nope=4
               elseif(lakon(ielem)(4:4).eq.'6') then
                  nope=6
               endif
               indexe=ipkon(ielem)
               do j=1,nope
                  node=kon(indexe+j)
                  sum=sum+pa-vold(4,node)
                  sumabs=sumabs+dabs(pa-vold(4,node))
               enddo
               index=inoel(3,index)
               if(index.eq.0) exit
            enddo
            if(sumabs.lt.1.d-10) then
               sum=0.d0
               sumabs=1.d0
            endif
            sa(i)=dabs(sum)/(sumabs*dt(i))
         enddo
c       endif
c!     
c!     viscous: pressure switch is calculated based on 
c!                 actual second derivative
c!
c         xmaxsum=0.d0
c         do i=1,nk
c            if(iponoel(i).le.0) cycle
c            sum=0.d0
c            pa=vold(4,i)
c            index=iponoel(i)
c            dd=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
c            if(dd.lt.1.d-10) then
c               v(1)=1.d0
c               v(2)=0.d0
c               v(3)=0.d0
c            else
c               v(1)=vold(1,i)/dd
c               v(2)=vold(2,i)/dd
c               v(3)=vold(3,i)/dd
c            endif
c            do
c               ielem=inoel(1,index)
c               if(ipkon(ielem).lt.0) cycle
c               if(lakon(ielem)(1:1).ne.'F') cycle
c               if(lakon(ielem)(4:4).eq.'2') then
c                  nope=20
c               elseif(lakon(ielem)(4:4).eq.'8') then
c                  nope=8
c               elseif(lakon(ielem)(4:4).eq.'4') then
c                  nope=4
c               elseif(lakon(ielem)(4:5).eq.'10') then
c                  nope=10
c               elseif(lakon(ielem)(4:4).eq.'6') then
c                  nope=6
c               elseif(lakon(ielem)(4:5).eq.'15') then
c                  nope=15
c               endif
c               indexe=ipkon(ielem)
c               do j=1,nope
c                  node=kon(indexe+j)
c                  if(node.eq.i) cycle
c                  p(1)=co(1,node)-co(1,i)
c                  p(2)=co(2,node)-co(2,i)
c                  p(3)=co(3,node)-co(3,i)
c                  dd=dsqrt(p(1)**2+
c     &                 p(2)**2+
c     &                 p(3)**2)
c                  cosang=dabs((p(1)*v(1)+p(2)*v(2)+p(3)*v(3)))/dd
c                  call nident(isolidsurf,node,nsolidsurf,id)
c                  if(id.gt.0) then
c                     if(isolidsurf(id).eq.node) then
c                        sum=sum+2.d0*cosang*(pa-vold(4,node))/dd**2
c                        cycle
c                     endif
c                  endif
c                  sum=sum+cosang*(pa-vold(4,node))/dd**2
c               enddo
c               index=inoel(3,index)
c               if(index.eq.0) exit
c            enddo
c            sa(i)=dabs(sum)
c            if(xmaxsum.lt.sa(i)) xmaxsum=sa(i)
c         enddo
c!     
c         if(xmaxsum.lt.1.e-30) xmaxsum=1.
c!     
c!     a lower exponent in the next line (but >0) creates
c!     more smoothing
c!     
c         if(euler.eq.1) then
c!
c!           nonviscous: basic smoothing of 0.1 augmented by
c!           a term dependent on the second pressure derivative
c!
c            do i=1,nk
c               sa(i)=(0.1d0+0.9d0*(sa(i)/xmaxsum)**(1.d0/2.d0))/dt(i)
c            enddo
c         else
c!
c!           viscous: only the second order derivative of the pressure
c!           is taken
c!
c            do i=1,nk
c               sa(i)=(sa(i)/xmaxsum)**(1.d0/2.d0)/dt(i)
c            enddo
c         endif
cc      endif
!     
      if(nmethod.eq.4) then
!
!        transient compressible
!
         ca=shockcoef*dtimef
         do i=1,nk
            sa(i)=ca*sa(i)
            sav(3*i-2)=sa(i)
            sav(3*i-1)=sa(i)
            sav(3*i)=sa(i)
         enddo
      else
!
!        steady state compressible
!
        ca=shockcoef*dtimef
         do i=1,nk
            sa(i)=ca*sa(i)
            sav(3*i-2)=sa(i)
            sav(3*i-1)=sa(i)
            sav(3*i)=sa(i)
c            vold(0,i)=sa(i)
         enddo
      endif
!
      return
      end
