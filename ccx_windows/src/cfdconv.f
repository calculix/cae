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
      subroutine cfdconv(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,iturbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor,vconini,
     &  dtimef,del,sum,sumx,sumxx,sumy,sumxy,nstart,shockcoef)
!
!     calculates the change in solution
!
      implicit none
!
      integer convergence,compressible,
     &  nrhcon(*),ntmat_,nactdoh(0:4,*),iit,iturbulent,mi(*),
     &  nshcon(*),ielmat(mi(3),*),nk,ithermal(*),i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun,nstart,
     &  nstartold
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),vcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:6),dummy,press,
     &  vconmax(0:6),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  vcontu(2,*),vtu(2,*),co(3,*),factor,vconini(0:4,*),
     &  dtimef,del(0:6,*),sum,sumx,sumxx,sumy(0:6),sumxy(0:6),
     &  b(0:6),bmax,denominator,shockcoef
!     
      do j=0,6
         vmax(j)=0.d0
         vconmax(j)=0.d0
      enddo
!
      if(compressible.eq.1) then
         do i=1,nk
            if(inomat(i).eq.0) cycle
            do j=0,4
               vmax(j)=vmax(j)+(vcon(j,i)-vconini(j,i))**2
               vconmax(j)=vconmax(j)+vconini(j,i)**2
               vconini(j,i)=vcon(j,i)
            enddo
         enddo
      else
         do i=1,nk
            if(inomat(i).eq.0) cycle
            do j=0,3
               vmax(j)=vmax(j)+(vcon(j,i)-vconini(j,i))**2
               vconmax(j)=vconmax(j)+vconini(j,i)**2
               vconini(j,i)=vcon(j,i)
            enddo
!
!           for incompressible fluids the pressure is stored
!           in vold(4,*), the initial pressure in 
!           vconini(4,*)
!
            do j=4,4
               vmax(j)=vmax(j)+(vold(j,i)-vconini(j,i))**2
               vconmax(j)=vconmax(j)+vconini(j,i)**2
               vconini(j,i)=vold(j,i)
            enddo
         enddo
      endif
      if(iturbulent.ne.0) then
         do i=1,nk
            if(inomat(i).eq.0) cycle
c            write(*,*) 'cfdconv.f ',i,vtu(1,i),vtu(2,i)
            do j=1,2
               vmax(4+j)=vmax(4+j)+vtu(j,i)**2
               vconmax(4+j)=vconmax(4+j)+vcontu(j,i)**2
            enddo
         enddo
      endif
!     
!     for steady state calculations: check convergence
!     
      convergence=0
!
      if(iturbulent.eq.0) then
!
!        laminar
!
         do i=0,4
            vmax(i)=dsqrt(vmax(i))
            vconmax(i)=dsqrt(vconmax(i))
            if(vconmax(i).lt.1.d-10) then
               vconmax(i)=1.d0
               vmax(i)=0.d0
            endif
         enddo
         if(nmethod.eq.1) then
            if(((vmax(0).lt.1.d-8*vconmax(0)).or.
     &           (vconmax(0).lt.1.d-10)).and.
     &           ((vmax(1).lt.1.d-8*vconmax(1)).or.
     &           (vconmax(1).lt.1.d-10)).and.
     &           ((vmax(2).lt.1.d-8*vconmax(2)).or.
     &           (vconmax(2).lt.1.d-10)).and.
     &           ((vmax(3).lt.1.d-8*vconmax(3)).or.
     &           (vconmax(3).lt.1.d-10)).and.
     &           ((vmax(4).lt.1.d-8*vconmax(4)).or.
     &           (vconmax(4).lt.1.d-10)).and.
     &           (iit.gt.1)) convergence=1
         endif
         if(iit.gt.1)
     &   write(12,'(i7,15(1x,e10.3))') iit,vmax(0)/vconmax(0),
     &     vmax(1)/vconmax(1),vmax(2)/vconmax(2),
     &     vmax(3)/vconmax(3),vmax(4)/vconmax(4),
     &     dtimef
!
!        linear regression of the last half of the iterations
!        for iterations exceeding 1000. If the regression is nearly
!        constant (zero slope) no improvement is deemed possible by further
!        iterations: convergence is assumed. This applies to
!        calculations with a nonzero shock coefficient 
!
         if((compressible.eq.1).and.(shockcoef>0)) then
            if(iit.gt.500) then
               sum=sum+1
               sumx=sumx+iit
               sumxx=sumxx+iit*iit
               do i=0,4
                  if(vmax(i).lt.1.d-10*vconmax(i)) then
                     del(i,iit)=-10.d0
                  else
                     del(i,iit)=dlog10(vmax(i)/vconmax(i))
                  endif
                  sumy(i)=sumy(i)+del(i,iit)
                  sumxy(i)=sumxy(i)+iit*del(i,iit)
               enddo
            endif
!     
!     check for zero slope
!     
            if(iit.gt.1000) then
               nstartold=nstart
               nstart=int(iit/2.d0)
               do j=nstartold,nstart-1
                  sum=sum-1
                  sumx=sumx-j
                  sumxx=sumxx-j*j
                  do i=0,4
                     sumy(i)=sumy(i)-del(i,j)
                     sumxy(i)=sumxy(i)-j*del(i,j)
                  enddo
               enddo
               bmax=0.d0
               do i=0,4
                  denominator=sumx*sumx-sum*sumxx
                  if(dabs(denominator).gt.1.d-30) then
                     b(i)=(sumx*sumy(i)-sum*sumxy(i))/denominator
                  else
                     b(i)=0.d0
                  endif
                  if(dabs(b(i)).gt.bmax) bmax=dabs(b(i))
               enddo
c               write(*,*) 'cfdconv ',iit,bmax
c               if(bmax.lt.1.d-4) convergence=1
               if(bmax.lt.1.d-6) convergence=1
            endif
         endif
      else
!
!        turbulent
!
         do i=0,6
            vmax(i)=dsqrt(vmax(i))
            vconmax(i)=dsqrt(vconmax(i))
            if(vconmax(i).lt.1.d-10) then
               vconmax(i)=1.d0
               vmax(i)=0.d0
            endif
         enddo
         if(ithermal(1).eq.0) vconmax(0)=1.d0
         if(nmethod.eq.1) then
            if(((vmax(0).lt.1.d-8*vconmax(0)).or.
     &           (vconmax(0).lt.1.d-10)).and.
     &           ((vmax(1).lt.1.d-8*vconmax(1)).or.
     &           (vconmax(1).lt.1.d-10)).and.
     &           ((vmax(2).lt.1.d-8*vconmax(2)).or.
     &           (vconmax(2).lt.1.d-10)).and.
     &           ((vmax(3).lt.1.d-8*vconmax(3)).or.
     &           (vconmax(3).lt.1.d-10)).and.
     &           ((vmax(4).lt.1.d-8*vconmax(4)).or.
     &           (vconmax(4).lt.1.d-10)).and.
     &           ((vmax(5).lt.1.d-8*vconmax(5)).or.
     &           (vconmax(5).lt.1.d-10)).and.
     &           ((vmax(6).lt.1.d-8*vconmax(6)).or.
     &           (vconmax(6).lt.1.d-10)).and.
     &           (iit.gt.1)) convergence=1
         endif
         if(iit.gt.1)
     &   write(12,'(i7,15(1x,e10.3))') iit,vmax(0)/vconmax(0),
     &     vmax(1)/vconmax(1),vmax(2)/vconmax(2),
     &     vmax(3)/vconmax(3),vmax(4)/vconmax(4),
     &     vmax(5)/vconmax(5),vmax(6)/vconmax(6),dtimef
      endif
!
!     checking change exceeding 0.1 %
!
c      factor=min(1.d0,1.01d0*factor)
c      if(dabs(vconmax(0)).gt.1.d-3) then
c         factor=min(factor,vconmax(0)/vmax(0)*0.001)
c      endif
c      if(dabs(vconmax(1)).gt.1.d-3) then
c         factor=min(factor,vconmax(1)/vmax(1)*0.001)
c      endif
c      if(dabs(vconmax(2)).gt.1.d-3) then
c         factor=min(factor,vconmax(2)/vmax(2)*0.001)
c      endif
c      if(dabs(vconmax(3)).gt.1.d-3) then
c         factor=min(factor,vconmax(3)/vmax(3)*0.001)
c      endif
c      if(dabs(vconmax(4)).gt.1.d-3) then
c         factor=min(factor,vconmax(4)/vmax(4)*0.001)
c      endif
c      if(iturbulent.ne.0) then
c         if(dabs(vconmax(5)).gt.1.d-3) then
c            factor=min(factor,vconmax(5)/vmax(5)*0.001)
c         endif
c         if(dabs(vconmax(6)).gt.1.d-3) then
c            factor=min(factor,vconmax(6)/vmax(6)*0.001)
c         endif
c      endif
!     
      return
      end
      
