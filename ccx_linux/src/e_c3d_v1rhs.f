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
      subroutine e_c3d_v1rhs(co,nk,konl,lakonl,p1,p2,omx,bodyfx,
     &     nbody,ff,nelem,nmethod,rhcon,nrhcon,ielmat,ntmat_,vold,vcon,
     &     idist,dtime,matname,mi,
     &     ttime,time,istep,iinc,shcon,nshcon,
     &     turbulent,vcontu,yy,nelemface,sideface,nface,compressible,
     &     ipvar,var,ipvarf,varf,sti,ithermal,dt)
!     
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 1 (correction *)
!     
!     ff: rhs 
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!     
      integer konl(8),ifaceq(8,6),nk,nbody,nelem,ithermal(*),mi(*),
     &     idist,i,j,k,i1,i2,j1,nmethod,ii,jj,jj1,id,ipointer,
     &     ig,kk,nrhcon(*),ielmat(mi(3),*),nshcon(*),ntmat_,nope,
     &     nopes,imat,turbulent,compressible,
     &     mint2d,mint3d,ifacet(6,4),ifacew(8,5),istep,iinc,
     &     iflag,k1,nelemface(*),nface,ipvar(*),index,ipvarf(*)
!     
      real*8 co(3,*),shp(4,8),dvi,p1(3),p2(3),dt(*),
     &     bodyfx(3),ff(78),bf(3),q(3),c1,c2,xsjmod,
     &     rhcon(0:1,ntmat_,*),vel(3),div,shcon(0:3,ntmat_,*),
     &     voldl(0:mi(2),8),xsj2(3),shp2(7,8),omcor,
     &     vold(0:mi(2),*),om,omx,const,xsj,temp,tt(3,3),
     &     vcon(0:4,*),vconl(0:4,8),rho,weight,shpv(8),t(3,3),
     &     vcontu(2,*),vcontul(2,8),cvel(3),vkl(3,3),corio(3),xkin,
     &     xtuf,vort,un,yy(*),yyl(8),y,f2,unt,umt,a1,arg2,
     &     var(*),varf(*),sti(6,mi(1),*),tu
!     
      real*8 dtime,ttime,time
!     
!     
!     
      include "gauss.f"
!     
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/),(/8,5/))
      iflag=3
      a1=0.31d0
!     
      imat=ielmat(1,nelem)
      amat=matname(imat)
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
!     temperature, velocity, conservative variables
!     (rho*velocity and rho) and if turbulent 
!     rho*turbulence variables
!     
      do i1=1,nope
        do i2=0,3
          voldl(i2,i1)=vold(i2,konl(i1))
        enddo
        do i2=1,4
          vconl(i2,i1)=vcon(i2,konl(i1))
        enddo
        if(turbulent.ne.0) then
          vcontul(1,i1)=vcontu(1,konl(i1))
          vcontul(2,i1)=vcontu(2,konl(i1))
          yyl(i1)=yy(konl(i1))
        endif
      enddo
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
!     calculating of
!     the temperature temp
!     the velocity vel
!     the auxiliary variable cvel
!     the velocity gradient vkl
!     the divergence of the velocity div 
!     the divergence of the shape function times the velocity shpv(*)
!     in the integration point
!     
        temp=0.d0
        do i1=1,3
          vel(i1)=0.d0
          cvel(i1)=0.d0
          do j1=1,3
            vkl(i1,j1)=0.d0
          enddo
        enddo
        do i1=1,nope
          temp=temp+shp(4,i1)*voldl(0,i1)
          do j1=1,3
            vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
            do k1=1,3
              vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
            enddo
          enddo
        enddo
        if(compressible.eq.1) then
          div=vkl(1,1)+vkl(2,2)+vkl(3,3)
        else
          div=0.d0
        endif
!     
        do i1=1,nope
          shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
     &         shp(3,i1)*vel(3)+shp(4,i1)*div
          do j1=1,3
            cvel(j1)=cvel(j1)+shpv(i1)*vconl(j1,i1)
          enddo
        enddo
!     
!     storing shpv, vel and temp
!     
        do i1=1,nope
          index=index+1
          var(index)=shpv(i1)
        enddo
        index=index+1
        var(index)=temp
!     
!     material data (density and dynamic viscosity)
!     
        call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     determining the dissipative stress 
!     
        do i1=1,3
          do j1=i1,3
            t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
          enddo
          if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
        enddo
!     
!     calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!     
        if(compressible.eq.1) then
!     
!     gas
!     
          if((nbody.ne.0).or.(turbulent.ne.0)) then
            rho=0.d0
            do i1=1,nope
              rho=rho+shp(4,i1)*vconl(4,i1)
            enddo
          endif
          if(turbulent.ne.0) then
            xkin=0.d0
            xtuf=0.d0
            y=0.d0
            do i1=1,nope
              xkin=xkin+shp(4,i1)*vcontul(1,i1)
              xtuf=xtuf+shp(4,i1)*vcontul(2,i1)
              y=y+shp(4,i1)*yyl(i1)
            enddo
            xkin=xkin/rho
            xtuf=xtuf/rho
            var(index+8)=rho
            var(index+9)=y
            var(index+10)=xkin
            var(index+11)=xtuf
          endif
        else
!     
!     liquid
!     
          if((nbody.ne.0).or.(turbulent.ne.0)) then
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
          endif
          if(turbulent.ne.0) then
            xkin=0.d0
            xtuf=0.d0
            y=0.d0
            do i1=1,nope
              xkin=xkin+shp(4,i1)*vcontul(1,i1)
              xtuf=xtuf+shp(4,i1)*vcontul(2,i1)
              y=y+shp(4,i1)*yyl(i1)
            enddo
            xkin=xkin/rho
            xtuf=xtuf/rho
            var(index+8)=rho
            var(index+9)=y
            var(index+10)=xkin
            var(index+11)=xtuf
          endif
        endif
!     
!     
!     adding the turbulent stress
!     
        if(turbulent.ne.0) then
!     
!     vorticity
!     
          vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &         (vkl(1,3)-vkl(3,1))**2+
     &         (vkl(2,1)-vkl(1,2))**2)
!     
!     kinematic viscosity
!     
          un=dvi/rho
!     
!     factor F2
!     
          c1=dsqrt(xkin)/(0.09d0*xtuf*y)
          c2=500.d0*un/(y*y*xtuf)
          arg2=max(2.d0*c1,c2)
          f2=dtanh(arg2*arg2)
!     
!     kinematic and dynamic turbulent viscosity
!     
          unt=a1*xkin/max(a1*xtuf,vort*f2)
c     unt=xkin/xtuf
          var(index+12)=unt
!     
          umt=unt*rho
!     
!     calculating the production (anisotropic part of
!     the turbulent stress is, apart from the dynamic
!     viscosity, identical to the viscous stress)
!     
c     tu=vort*vort
          tu=(t(1,1)*vkl(1,1)+t(1,2)*vkl(1,2)+t(1,3)*vkl(1,3)+
     &         t(1,2)*vkl(2,1)+t(2,2)*vkl(2,2)+t(2,3)*vkl(2,3)+
     &         t(1,3)*vkl(3,1)+t(2,3)*vkl(3,2)+t(3,3)*vkl(3,3))
!     
!     correction for compressible fluids
!     
          if(compressible.eq.1) then
            tu=tu-2.d0*xtuf*div/3.d0
          endif
          var(index+13)=tu
!     
!     calculating the turbulent stress
!     
          do i1=1,3
            do j1=i1,3
              tt(i1,j1)=t(i1,j1)
            enddo
            tt(i1,i1)=tt(i1,i1)-2.d0*xtuf/3.d0
          enddo
!     
!     adding the viscous stress
!     
          do i1=1,3
            do j1=i1,3
              t(i1,j1)=dvi*t(i1,j1)+umt*tt(i1,j1)
c     t(i1,j1)=dvi*t(i1,j1)
            enddo
          enddo
        else
!     
          do i1=1,3
            do j1=i1,3
              t(i1,j1)=dvi*t(i1,j1)
            enddo
          enddo
        endif
!     
!     storing the total dissipative stress x velocity
!     (viscous + turbulent)
!     
        index=index+1
        var(index)=t(1,1)*vel(1)+t(1,2)*vel(2)+t(1,3)*vel(3)
        index=index+1
        var(index)=t(1,2)*vel(1)+t(2,2)*vel(2)+t(2,3)*vel(3)
        index=index+1
        var(index)=t(1,3)*vel(1)+t(2,3)*vel(2)+t(3,3)*vel(3)
!     
!     determination of lhs and rhs
!     
        jj1=1
        do jj=1,nope
!     
!     convective + diffusive
!     
          ff(jj1)=ff(jj1)-xsjmod*
     &         (cvel(1)*(shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)
     &         +shp(1,jj)*t(1,1)+shp(2,jj)*t(1,2)+shp(3,jj)*t(1,3))
          ff(jj1+1)=ff(jj1+1)-xsjmod*
     &         (cvel(2)*(shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)
     &         +shp(1,jj)*t(1,2)+shp(2,jj)*t(2,2)+shp(3,jj)*t(2,3))
          ff(jj1+2)=ff(jj1+2)-xsjmod*
     &         (cvel(3)*(shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)
     &         +shp(1,jj)*t(1,3)+shp(2,jj)*t(2,3)+shp(3,jj)*t(3,3))
          jj1=jj1+3
        enddo
!     
!     computation of contribution due to body forces
!     
        if(nbody.ne.0) then
!     
!     initialisation for the body forces
!     
          om=omx*rho
          omcor=2.d0*rho*dsqrt(omx)
          do ii=1,3
            bf(ii)=bodyfx(ii)*rho
          enddo
!     
          if(om.gt.0.d0) then
            do i1=1,3
!     
!     computation of the global coordinates of the gauss
!     point
!     
              q(i1)=0.d0
              do j1=1,nope
                q(i1)=q(i1)+shp(4,j1)*co(i1,konl(j1))
              enddo
!     
              q(i1)=q(i1)-p1(i1)
            enddo
            const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!     
!     Coriolis forces
!     
            corio(1)=vel(2)*p2(3)-vel(3)*p2(2)
            corio(2)=vel(3)*p2(1)-vel(1)*p2(3)
            corio(3)=vel(1)*p2(2)-vel(2)*p2(1)
!     
!     inclusion of the centrifugal force into the body force
!     
            do i1=1,3
              bf(i1)=bf(i1)+(q(i1)-const*p2(i1))*om+
     &             corio(i1)*omcor
            enddo
          endif
!     
!     storing the body force
!     
          index=index+1
          var(index)=bf(1)*vel(1)+bf(2)*vel(2)+bf(3)*vel(3)
          index=index+9
!     
          jj1=1
          do jj=1,nope
            ff(jj1)=ff(jj1)+xsjmod*bf(1)*(shp(4,jj)+
     &           dt(konl(jj))*shpv(jj)/2.d0)
            ff(jj1+1)=ff(jj1+1)+xsjmod*bf(2)*(shp(4,jj)+
     &           dt(konl(jj))*shpv(jj)/2.d0)
            ff(jj1+2)=ff(jj1+2)+xsjmod*bf(3)*(shp(4,jj)+
     &           dt(konl(jj))*shpv(jj)/2.d0)
            jj1=jj1+3
          enddo
        else
          index=index+10
        endif
!     
      enddo
!     
      if(nface.ne.0) then
        index=ipvarf(nelem)
!     
!     free stream or solid surface boundaries
!     
        nopes=0
        call nident(nelemface,nelem,nface,id)
        do
          if((id.eq.0).or.(nelemface(id).ne.nelem)) exit
          ig=ichar(sideface(id)(1:1))-48
!     
          if(nopes.eq.0) then
            if(lakonl(4:4).eq.'4') then
              nopes=3
              mint2d=1
            elseif(lakonl(4:4).eq.'6') then
              mint2d=1
            elseif(lakonl(4:5).eq.'8R') then
              nopes=4
              mint2d=1
            elseif(lakonl(4:4).eq.'8') then
              nopes=4
              mint2d=4
            endif
          endif
!     
          if(lakonl(4:4).eq.'6') then
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
!     
          do i=1,mint2d
!     
!     local coordinates of the surface integration
!     point within the surface local coordinate system
!     
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
              weight=weight2d1(i)
            elseif(lakonl(4:4).eq.'8') then
              weight=weight2d2(i)
            elseif((lakonl(4:4).eq.'4').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
              weight=weight2d4(i)
            endif
!     
            do i1=1,nopes
              index=index+1
              shp2(4,i1)=varf(index)
            enddo
            do i1=1,3
              index=index+1
              xsj2(i1)=varf(index)
            enddo
!     
!     derivative of the volumetric shape functions
!     needed for the dissipative stress
!     
            do i1=1,nope
              do j1=1,4
                index=index+1
                shp(j1,i1)=varf(index)
              enddo
            enddo
!     
!     calculating of
!     the temperature temp
!     the velocity gradient vkl
!     in the integration point
!     
            temp=0.d0
            do i1=1,3
              vel(i1)=0.d0
              do j1=1,3
                vkl(i1,j1)=0.d0
              enddo
            enddo
            do i1=1,nope
              temp=temp+shp(4,i1)*voldl(0,i1)
              do j1=1,3
                vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
                do k1=1,3
                  vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                enddo
              enddo
            enddo
!     
!     storing the temperature
!     
            index=index+1
            varf(index)=temp
!     
            if(compressible.eq.1) then
              div=vkl(1,1)+vkl(2,2)+vkl(3,3)
            else
              div=0.d0
            endif
!     
!     material data (dynamic viscosity)
!     
            call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,
     &           dvi)
!     
!     determining the dissipative stress 
!     
            do i1=1,3
              do j1=i1,3
                t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
              enddo
              if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
            enddo
!     
!     calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!     
            if(turbulent.ne.0) then
              if(compressible.eq.1) then
!     
!     gas
!     
                rho=0.d0
                do i1=1,nope
                  rho=rho+shp(4,i1)*vconl(4,i1)
                enddo
              else
!     
!     liquid
!     
                call materialdata_rho(rhcon,nrhcon,imat,rho,
     &               temp,ntmat_,ithermal)
              endif
              y=0.d0
              xkin=0.d0
              xtuf=0.d0
              do i1=1,nope
                y=y+shp(4,i1)*yyl(i1)
                xkin=xkin+shp(4,i1)*vcontul(1,i1)
                xtuf=xtuf+shp(4,i1)*vcontul(2,i1)
              enddo
              xkin=xkin/rho
              xtuf=xtuf/rho
!     
              varf(index+4)=rho
              varf(index+5)=y
              varf(index+6)=xkin
              varf(index+7)=xtuf
            endif
!     
!     calculation of turbulent auxiliary variables
!     
            if(turbulent.ne.0) then
!     
!     vorticity
!     
              vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &             (vkl(1,3)-vkl(3,1))**2+
     &             (vkl(2,1)-vkl(1,2))**2)
!     
!     kinematic viscosity
!     
              un=dvi/rho
!     
!     factor F2
!     
              c1=dsqrt(xkin)/(0.09d0*xtuf*y)
              c2=500.d0*un/(y*y*xtuf)
              arg2=max(2.d0*c1,c2)
              f2=dtanh(arg2*arg2)
!     
!     kinematic and dynamic turbulent viscosity
!     
              unt=a1*xkin/max(a1*xtuf,vort*f2)
c     unt=xkin/xtuf
              varf(index+8)=unt
!     
              umt=unt*rho
!     
              do i1=1,3
                do j1=i1,3
                  t(i1,j1)=(dvi+umt)*t(i1,j1)
c     t(i1,j1)=(dvi)*t(i1,j1)
                enddo
                t(i1,i1)=t(i1,i1)-2.d0*rho*xkin/3.d0
              enddo
            else
!     
              do i1=1,3
                do j1=i1,3
                  t(i1,j1)=dvi*t(i1,j1)
                enddo
              enddo
            endif
!     
!     storing stress . velocity
!     
            index=index+1
            varf(index)=t(1,1)*vel(1)+t(1,2)*vel(2)+
     &           t(1,3)*vel(3)
            index=index+1
            varf(index)=t(1,2)*vel(1)+t(2,2)*vel(2)+
     &           t(2,3)*vel(3)
            index=index+1
            varf(index)=t(1,3)*vel(1)+t(2,3)*vel(2)+
     &           t(3,3)*vel(3)
            index=index+5
!     
            do k=1,nopes
              if(nope.eq.8) then
                ipointer=(ifaceq(k,ig)-1)*3
              elseif(nope.eq.4) then
                ipointer=(ifacet(k,ig)-1)*3
              else
                ipointer=(ifacew(k,ig)-1)*3
              endif
              ff(ipointer+1)=ff(ipointer+1)+shp2(4,k)*
     &             (t(1,1)*xsj2(1)+t(1,2)*xsj2(2)+t(1,3)*xsj2(3))*
     &             weight*dtime
              ff(ipointer+2)=ff(ipointer+2)+shp2(4,k)*
     &             (t(1,2)*xsj2(1)+t(2,2)*xsj2(2)+t(2,3)*xsj2(3))*
     &             weight*dtime
              ff(ipointer+3)=ff(ipointer+3)+shp2(4,k)*
     &             (t(1,3)*xsj2(1)+t(2,3)*xsj2(2)+t(3,3)*xsj2(3))*
     &             weight*dtime
            enddo
          enddo
          id=id-1
        enddo
      endif
!     
      return
      end
