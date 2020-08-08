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
      subroutine e_c3d_krhs(co,nk,konl,lakonl,ffk,fft,nelem,nmethod,
     &     rhcon,nrhcon,ielmat,ntmat_,vold,vcon,dtime,matname,mi,
     &     shcon,nshcon,vcontu,compressible,yy,nelemface,sideface,nface,
     &     iturbulent,ithermal,ipvar,var,ipvarf,varf,dt,ggk,ggt)
!     
!     computation of the turbulence element matrix and rhs for the
!     element with the topology in konl: step 4
!     
!     ffk and fft: rhs (x 2: kinetic term and turbulence frequency term):
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!     
      integer konl(20),ifaceq(8,6),nk,nelem,nload,i,j,k,i1,i2,j1,k1,
     &     nmethod,ii,jj,id,ipointer,ig,kk,nrhcon(*),mi(*),
     &     ielmat(mi(3),*),nshcon(*),ipvar(*),ipvarf(*),index,
     &     ntmat_,nope,nopes,imat,mint2d,mint3d,ifacet(6,4),nopev,
     &     ifacew(8,5),istep,iinc,layer,kspt,jltyp,iflag,iscale,
     &     ithermal(*),
     &     compressible,idf,igl,nelemface(*),nface,iturbulent
!     
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),dvi,prod,diss,
     &     ffk(78),xsjmod,vkl(3,3),rhcon(0:1,ntmat_,*),reltime,
     &     t(3,3),bfv,press,vel(3),div,shcon(0:3,ntmat_,*),pgauss(3),
     &     xkin,xtuf,voldl(0:mi(2),20),yyl(20),tvk(3),tvt(3),
     &     xl2(3,8),xsj2(3),shp2(7,8),vold(0:mi(2),*),tvnk,tvnt,
     &     om,omx,xi,et,ze,const,xsj,fft(78),dxkin(3),dt(*),
     &     temp,vcon(0:4,*),vconl(0:4,20),rho,dxtuf(3),ggk(78),
     &     weight,shpv(20),rhokin,rhotuf,y,vort,c1,c2,arg2,f2,
     &     a1,unt,umt,cdktuf,arg1,f1,skin,skin1,skin2,stuf,stuf1,
     &     stuf2,beta,beta1,beta2,betas,gamm,gamm1,xkappa,un,
     &     gamm2,umsk,umst,tu,tuk,tut,vcontu(2,*),vcontul(2,20),
     &     f1m,yy(*),xsjmodk,xsjmodt,xi3d,et3d,ze3d,xlocal20(3,9,6),
     &     xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),ggt(78),
     &     xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),var(*),varf(*)
!     
      real*8 dtime,ttime,time
!     
      include "gauss.f"
      include "xlocal.f"
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
!     
!     turbulence constants
!     
c     
c     following constant is for the SST Modell
      skin1=0.85d0
c     skin1=0.5d0
      skin2=1.d0
      stuf1=0.5d0
      stuf2=0.856d0
      beta1=0.075d0
      beta2=0.0828d0
      a1=0.31d0
      betas=0.09d0
      xkappa= 0.41d0
!     
      gamm1=beta1/betas-stuf1*xkappa*xkappa/dsqrt(betas)
      gamm2=beta2/betas-stuf2*xkappa*xkappa/dsqrt(betas)
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
      elseif(lakonl(4:5).eq.'10') then
        nope=10
        mint3d=4
      elseif(lakonl(4:5).eq.'15') then
        nope=15
        mint3d=9
      elseif(lakonl(4:6).eq.'20R') then
        nope=20
        mint3d=8
      elseif(lakonl(4:4).eq.'2') then
        nope=20
        mint3d=27
      else
        mint3d=0
      endif
!     
!     initialisation for distributed forces
!     
      do i=1,nope
        ffk(i)=0.d0
        fft(i)=0.d0
        ggk(i)=0.d0
        ggt(i)=0.d0
      enddo
!     
!     temperature, velocity and conservative variables
!     (rho*energy density, rho*velocity and rho)
!     
      do i1=1,nope
        vcontul(1,i1)=vcontu(1,konl(i1))
        vcontul(2,i1)=vcontu(2,konl(i1))
      enddo
!     
!     computation of the matrix: loop over the Gauss points
!     
      index=ipvar(nelem)
      do kk=1,mint3d
        if(lakonl(4:5).eq.'8R') then
          weight=weight3d1(kk)
        elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &         then
          weight=weight3d2(kk)
        elseif(lakonl(4:4).eq.'2') then
          weight=weight3d3(kk)
        elseif(lakonl(4:5).eq.'10') then
          weight=weight3d5(kk)
        elseif(lakonl(4:4).eq.'4') then
          weight=weight3d4(kk)
        elseif(lakonl(4:5).eq.'15') then
          weight=weight3d8(kk)
        elseif(lakonl(4:5).eq.'6 ') then
          weight=weight3d7(kk)
        elseif(lakonl(4:5).eq.'6R') then
          weight=weight3d11(kk)
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
!     retrieving shpv and temp
!     
        do i1=1,nope
          index=index+1
          shpv(i1)=var(index)
        enddo
        index=index+1
        temp=var(index)
!     
!     calculating of
!     rho times turbulent kinetic energy times shpv(*): rhokin
!     rho times turbulence frequency times shpv(*): rhotuf
!     
        rhokin=0.d0
        rhotuf=0.d0
        do i1=1,nope
        enddo
        do i1=1,nope
          rhokin=rhokin+shpv(i1)*vcontul(1,i1)
          rhotuf=rhotuf+shpv(i1)*vcontul(2,i1)
        enddo
!     
!     material data (density and dynamic viscosity)
!     
        call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     calculation of the spatial derivatives of the turbulent kinetic energy
!     and the turbulence frequency for gases and liquids
!     
        rho=var(index+8)
        y=var(index+9)
        xkin=var(index+10)
        xtuf=var(index+11)
!     
        if(compressible.eq.1) then
!     
!     gas
!     
          do j1=1,3
            dxkin(j1)=0.d0
            dxtuf(j1)=0.d0
          enddo
          do i1=1,nope
            do j1=1,3
              dxkin(j1)=dxkin(j1)+shp(j1,i1)*vcontul(1,i1)
              dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*vcontul(2,i1)
            enddo
          enddo
          do j1=1,3
            dxkin(j1)=dxkin(j1)/rho
            dxtuf(j1)=dxtuf(j1)/rho
          enddo
        else
!     
!     liquid
!     
          do j1=1,3
            dxkin(j1)=0.d0
            dxtuf(j1)=0.d0
          enddo
          do i1=1,nope
            do j1=1,3
              dxkin(j1)=dxkin(j1)+shp(j1,i1)*vcontul(1,i1)
              dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*vcontul(2,i1)
            enddo
          enddo
          do j1=1,3
            dxkin(j1)=dxkin(j1)/rho
            dxtuf(j1)=dxtuf(j1)/rho
          enddo
        endif
!     
!     calculation of turbulent auxiliary variables
!     
c     !        vorticity
c     !
c     vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
c     &              (vkl(1,3)-vkl(3,1))**2+
c     &              (vkl(2,1)-vkl(1,2))**2)
c     !
c     !        kinematic viscosity
c     !
        un=dvi/rho
c     !
c     !        factor F2
c     !
        c1=dsqrt(xkin)/(0.09d0*xtuf*y)
        c2=500.d0*un/(y*y*xtuf)
c     arg2=max(2.d0*c1,c2)
c     f2=dtanh(arg2*arg2)
!     
!     kinematic and dynamic turbulent viscosity
!     
c     unt=a1*xkin/max(a1*xtuf,vort*f2)
c     unt=xkin/xtuf
        unt=var(index+12)
c     write(*,*) 'e_c3d_krhs unt ',unt
c     
        umt=unt*rho
!     
!     factor F1
!     
        if(iturbulent.eq.1) then
!     
!     k-epsilon model
!     
          f1=0.d0
        elseif(iturbulent.eq.2) then
!     
!     k-omega model
!     
          f1=1.d0
        else
!     
!     SST model
!     
          cdktuf=max(2.d0*rho*stuf2*
     &         (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+dxkin(3)*dxtuf(3))/
     &         xtuf,1.d-20)
          arg1=min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
          f1=dtanh(arg1**4.d0)
        endif
        f1m=1.d0-f1
!     
!     interpolation of the constants
!     
        skin=f1*skin1+f1m*skin2
        stuf=f1*stuf1+f1m*stuf2
        beta=f1*beta1+f1m*beta2
        gamm=f1*gamm1+f1m*gamm2
!     
!     auxiliary quantities
!     
        umsk=dvi+skin*umt
        umst=dvi+stuf*umt
        tu=var(index+13)
        index=index+13
!     
c     c         prod=unt*tu
c     c         diss=betas*xtuf*xkin
c     c         write(*,*) 'e_c3d_krhs ',prod,20*diss
c     c         if(prod.gt.20.d0*diss) prod=20.d0*diss
        tuk=rho*(unt*tu-betas*xtuf*xkin)
        tut=rho*(gamm*tu-beta*xtuf*xtuf+2.d0*f1m*stuf2*
     &       (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+dxkin(3)*dxtuf(3))/
     &       xtuf)
!     
        do i1=1,3
          dxkin(i1)=dxkin(i1)*umsk
          dxtuf(i1)=dxtuf(i1)*umst
        enddo
!     
!     determination of lhs and rhs
!     
        do jj=1,nope
!     
          ffk(jj)=ffk(jj)-
     &         xsjmod*((shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)*
     &         (rhokin-tuk)+(shp(1,jj)*dxkin(1)+shp(2,jj)*dxkin(2)
     &         +shp(3,jj)*dxkin(3)))
          fft(jj)=fft(jj)-
     &         xsjmod*((shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)*
     &         (rhotuf-tut)+(shp(1,jj)*dxtuf(1)+shp(2,jj)*dxtuf(2)
     &         +shp(3,jj)*dxtuf(3)))
          ggk(jj)=ggk(jj)+betas*xsjmod*shp(4,jj)*xtuf
          ggt(jj)=ggt(jj)+2.d0*beta*xsjmod*shp(4,jj)*xtuf
c     write(*,*) 'e_c3d_krsh ',jj,ffk(jj),fft(jj)
        enddo
!     
      enddo
!     
      if(nface.ne.0) then
        index=ipvarf(nelem)
!     
!     free stream or solid surface boundaries
!     
        nopes=0
        call nident(nelemface,nelem,nface,idf)
        do
          if((idf.eq.0).or.(nelemface(idf).ne.nelem)) exit
          read(sideface(idf)(1:1),'(i1)') ig
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
            elseif(lakonl(4:5).eq.'10') then
              nopes=6
              mint2d=3
            elseif(lakonl(4:6).eq.'20R') then
              nopes=8
              mint2d=4
            elseif(lakonl(4:4).eq.'2') then
              nopes=8
              mint2d=9
            endif
          endif
!     
          if(lakonl(4:4).eq.'6') then
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          elseif(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
              nopes=6
              mint2d=3
            else
              nopes=8
              mint2d=4
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
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
              weight=weight2d2(i)
            elseif(lakonl(4:4).eq.'2') then
              weight=weight2d3(i)
            elseif((lakonl(4:5).eq.'10').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
              weight=weight2d5(i)
            elseif((lakonl(4:4).eq.'4').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
              weight=weight2d4(i)
            endif
!     
!     local surface normal
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
!     retrieving the temperature
!     
            index=index+1
            temp=varf(index)
!     
!     retrieving other variables
!     
            rho=varf(index+4)
            y=varf(index+5)
            xkin=varf(index+6)
            xtuf=varf(index+7)
!     
!     material data (dynamic viscosity)
!     
            call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     calculation of the spatial derivatives of the turbulent kinetic energy
!     and the turbulence frequency for gases and liquids
!     
            do j1=1,3
              dxkin(j1)=0.d0
              dxtuf(j1)=0.d0
            enddo
            do i1=1,nope
              do j1=1,3
                dxkin(j1)=dxkin(j1)+shp(j1,i1)*vcontul(1,i1)
                dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*vcontul(2,i1)
              enddo
            enddo
            do j1=1,3
              dxkin(j1)=dxkin(j1)/rho
              dxtuf(j1)=dxtuf(j1)/rho
            enddo
!     
!     calculation of turbulent auxiliary variables
c     !     
c     !     vorticity
c     !     
c     vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
c     &              (vkl(1,3)-vkl(3,1))**2+
c     &              (vkl(2,1)-vkl(1,2))**2)
!     
!     kinematic viscosity
!     
            un=dvi/rho
!     
!     factor F2
!     
            if(y.gt.0.d0) then
              c1=dsqrt(xkin)/(0.09d0*xtuf*y)
              c2=500.d0*un/(y*y*xtuf)
            endif
c     arg2=max(2.d0*c1,c2)
c     f2=dtanh(arg2*arg2)
!     
!     kinematic and dynamic turbulent viscosity
!     
c     unt=a1*xkin/max(a1*xtuf,vort*f2)
            index=index+8
            unt=varf(index)
            umt=unt*rho
!     
!     factor F1
!     
            if(iturbulent.eq.1) then
!     
!     k-epsilon model
!     
              f1=0.d0
            elseif(iturbulent.eq.2) then
!     
!     k-omega model
!     
              f1=1.d0
            else
!     
!     SST model
!     
              if(y.gt.0.d0) then
!     
!     finite distance from wall
!     
                cdktuf=max(2.d0*rho*stuf2*
     &               (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+
     &               dxkin(3)*dxtuf(3))/xtuf,1.d-20)
                arg1=
     &               min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
                f1=dtanh(arg1**4.d0)
              else
!     
!     wall
!     
                f1=1.d0
              endif
            endif
            f1m=1.d0-f1
!     
!     interpolation of the constants
!     
            skin=f1*skin1+f1m*skin2
            stuf=f1*stuf1+f1m*stuf2
!     
!     auxiliary quantities
!     
            umsk=dvi+skin*umt
            umst=dvi+stuf*umt
!     
!     determining the stress and and stress x velocity + conductivity x
!     temperature gradient
!     
            do i1=1,3
              tvk(i1)=umsk*dxkin(i1)
              tvt(i1)=umsk*dxtuf(i1)
            enddo
!     
            tvnk=tvk(1)*xsj2(1)+tvk(2)*xsj2(2)+tvk(3)*xsj2(3)
            tvnt=tvt(1)*xsj2(1)+tvt(2)*xsj2(2)+tvt(3)*xsj2(3)
!     
            xsjmodk=tvnk*weight*dtime
            xsjmodt=tvnt*weight*dtime
            do k=1,nopes
              if((nope.eq.20).or.(nope.eq.8)) then
                ipointer=ifaceq(k,ig)
              elseif((nope.eq.10).or.(nope.eq.4)) then
                ipointer=ifacet(k,ig)
              else
                ipointer=ifacew(k,ig)
              endif
              ffk(ipointer)=ffk(ipointer)+shp2(4,k)*xsjmodk
              fft(ipointer)=fft(ipointer)+shp2(4,k)*xsjmodt
            enddo
          enddo
          idf=idf-1
        enddo
      endif
!     
      return
      end

