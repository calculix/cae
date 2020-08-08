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
      subroutine e_c3d_v2rhs(co,nk,konl,lakonl,
     &     ff,nelem,nmethod,vold,v,dtime,theta2,iexplicit,mi,
     &     ipvar,var,ipvarf,varf,dt)
!     
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 3 (correction **)
!     
!     ff: rhs 
!     
      implicit none
!     
      character*8 lakonl
!     
      integer konl(8),nk,nelem,i,j,i1,j1,nmethod,jj,jj1,kk,
     &     nope,mint3d,iflag,iexplicit,mi(*),ipvar(*),index,ii,
     &     ipvarf(*)
!     
      real*8 co(3,*),xl(3,8),shp(4,8),ff(78),xsjmod,vl(0:mi(2),8),
     &     vel(3),div,voldl(0:mi(2),8),v(0:mi(2),*),vold(0:mi(2),*),
     &     term,var(*),varf(*),dt(*),dtime,
     &     xi,et,ze,xsj,weight,shpv(8),theta2,dpress(3),ddpress(3)
!     
!     
!     
      include "gauss.f"
!     
      iflag=3
!     
      if(lakonl(4:4).eq.'4') then
        nope=4
        mint3d=1
      elseif(lakonl(4:4).eq.'6') then
        nope=6
        mint3d=2
      elseif(lakonl(4:5).eq.'8R') then
        nope=8
        mint3d=1
      elseif(lakonl(4:4).eq.'8') then
        nope=8
        mint3d=8
      endif
!     
!     initialisation for distributed forces
!     
      do i=1,3*nope
        ff(i)=0.d0
      enddo
!     
!     temperature, velocity and conservative variables
!     (rho*energy density, rho*velocity and rho)
!     
      if(iexplicit.eq.0) then
        do i1=1,nope
          vl(4,i1)=v(4,konl(i1))
        enddo
      endif
!     
!     computation of the matrix: loop over the Gauss points
!     
      index=ipvar(nelem)
      do kk=1,mint3d
        if(lakonl(4:5).eq.'8R') then
          weight=weight3d1(kk)
        elseif(lakonl(4:4).eq.'8') then
          weight=weight3d2(kk)
        elseif(lakonl(4:4).eq.'4') then
          weight=weight3d4(kk)
        elseif(lakonl(4:5).eq.'6 ') then
          weight=weight3d7(kk)
        endif
!     
!     copying the shape functions, their derivatives and the
!     Jacobian determinant from field var
!     
        do jj=1,nope
          do ii=1,4
            index=index+1
            shp(ii,jj)=var(index)
          enddo
        enddo
        index=index+1
        xsj=var(index)
!     
        xsjmod=dtime*xsj*weight
!     
!     calculating the temperature temp, the velocity vel, the
!     divergence of the velocity div and the divergence
!     of the shape function times the velocity shpv(*)
!     in the integration point
!     
        do i1=1,nope
          index=index+1
          shpv(i1)=var(index)
        enddo
        index=index+5
        do i1=1,3
          index=index+1
          dpress(i1)=var(index)
        enddo
        index=index+6
!     
!     only for the semi-implicit procedure: calculate ddpress
!     
        if(iexplicit.eq.0) then
          do i1=1,3
            ddpress(i1)=0.d0
          enddo
          do i1=1,nope
            do j1=1,3
              ddpress(j1)=ddpress(j1)+shp(j1,i1)*vl(4,i1)
            enddo
          enddo
        endif
!     
!     determination of rhs
!     
        if(iexplicit.eq.1) then
          jj1=1
          do jj=1,nope
            term=xsjmod*(shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)
            ff(jj1)=ff(jj1)-dpress(1)*term
            ff(jj1+1)=ff(jj1+1)-dpress(2)*term
            ff(jj1+2)=ff(jj1+2)-dpress(3)*term
            jj1=jj1+3
          enddo
        else
          jj1=1
          do jj=1,nope
!     
!     with stability term 
!
!     change on 20200719: the factor of 2 is wrong.
!            
            ff(jj1)=ff(jj1)-xsjmod*(dpress(1)*(shp(4,jj)+
     &           (1.d0-theta2)*dt(konl(jj))
     &           *shpv(jj))+theta2*shp(4,jj)*ddpress(1))
            ff(jj1+1)=ff(jj1+1)-xsjmod*(dpress(2)*(shp(4,jj)+
     &           (1.d0-theta2)*dt(konl(jj))
     &           *shpv(jj))+theta2*shp(4,jj)*ddpress(2))
            ff(jj1+2)=ff(jj1+2)-xsjmod*(dpress(3)*(shp(4,jj)+
     &           (1.d0-theta2)*dt(konl(jj))
     &           *shpv(jj))+theta2*shp(4,jj)*ddpress(3))
c            ff(jj1)=ff(jj1)-xsjmod*(dpress(1)*(shp(4,jj)+
c     &           (1.d0-theta2)*dt(konl(jj))
c     &           *shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(1))
c            ff(jj1+1)=ff(jj1+1)-xsjmod*(dpress(2)*(shp(4,jj)+
c     &           (1.d0-theta2)*dt(konl(jj))
c     &           *shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(2))
c            ff(jj1+2)=ff(jj1+2)-xsjmod*(dpress(3)*(shp(4,jj)+
c     &           (1.d0-theta2)*dt(konl(jj))
c     &           *shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(3))
            jj1=jj1+3
          enddo
        endif
!     
      enddo
!     
      return
      end
