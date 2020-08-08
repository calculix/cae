!     
!     CalculiX 3-dimensional finite element program
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
      subroutine e_c3d_trhs(co,nk,konl,lakonl,p1,p2,omx,bodyfx,
     &     nbody,ff,nelem,nmethod,rhcon,nrhcon,
     &     ielmat,ntmat_,vold,vcon,nelemload,
     &     sideload,xload,nload,idist,dtime,matname,mi,
     &     ttime,time,istep,iinc,xloadold,reltimef,shcon,nshcon,cocon,
     &     ncocon,physcon,nelemface,sideface,nface,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,compressible,
     &     vcontu,yy,turbulent,ipvar,var,ipvarf,varf,dt)
!     
!     computation of the energy element matrix and rhs for the element with
!     element with the topology in konl: step 4
!     
!     ff: rhs 
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
!     
      integer konl(8),ifaceq(8,6),nelemload(2,*),nk,nbody,nelem,
     &     nload,idist,i,j,k,i1,i2,j1,ncocon(2,*),k1,node,nfield,mi(*),
     &     nmethod,ii,jj,id,ipointer,ig,kk,nrhcon(*),ielmat(mi(3),*),
     &     nshcon(*),flux,compressible,
     &     ntmat_,nope,nopes,imat,mint2d,mint3d,ifacet(6,4),nopev,
     &     ifacew(8,5),istep,iinc,layer,kspt,jltyp,iflag,nelemface(*),
     &     nface,igl,idf,ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),
     &     iscale,turbulent,ipvar(*),index,ipvarf(*),iemchange
!     
      real*8 co(3,*),shp(4,8),xs2(3,7),dvi,dt(*),
     &     p1(3),p2(3),bodyf(3),bodyfx(3),ff(78),cond,enthalpy,
     &     bf(3),q(3),xsjmod,dtem(3),vkl(3,3),corio(3),sinktemp,
     &     rhcon(0:1,ntmat_,*),reltimef,t(3,3),tv(3),bfv,press,
     &     vel(3),div,shcon(0:3,ntmat_,*),pgauss(3),dxsj2,areaj,
     &     voldl(0:mi(2),8),xloadold(2,*),cocon(0:6,ntmat_,*),
     &     xl2(3,8),xsj2(3),shp2(7,8),vold(0:mi(2),*),xload(2,*),
     &     om,omx,xi,et,ze,const,xsj,field,physcon(*),tvn,
     &     temp,vcon(0:4,*),vconl(0:4,8),rho,xi3d,et3d,ze3d,
     &     weight,shpv(8),xlocal20(3,9,6),coefmpc(*),
     &     xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),vl(0:mi(2),8),
     &     xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),omcor,
     &     shpvnithi(8),vcontu(2,*),vcontul(2,8),yy(*),yyl(8),
     &     y,xtuf,xkin,vort,un,unt,umt,f2,arg2,c1,c2,a1,var(*),varf(*)
!     
      real*8 dtime,ttime,time,tvar(2),coords(3)
!     
!     
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
      a1=0.31d0
!     
      tvar(1)=time
      tvar(2)=ttime+dtime
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
      do i=1,nope
        ff(i)=0.d0
      enddo
!     
!     temperature, velocity and conservative variables
!     (rho*energy density, rho*velocity and rho)
!     
      do i1=1,nope
        voldl(0,i1)=vold(0,konl(i1))
        voldl(4,i1)=vold(4,konl(i1))
        vconl(0,i1)=vcon(0,konl(i1))
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
!     calculating of the temperature gradient dtem
!     in the integration point
!     
        enthalpy=0.d0
        do i1=1,3
          dtem(i1)=0.d0
        enddo
        do i1=1,nope
          do j1=1,3
            dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
          enddo
        enddo
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
!     calculating the total enthalpy     
!     in the integration point
!     
        do i1=1,nope
          enthalpy=enthalpy+shpv(i1)*(vconl(0,i1)+voldl(4,i1))
        enddo
!     
!     material data (conductivity)
!     
        call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,cond)
!     
!     retrieving stress x velocity
!     
        do i1=1,3
          index=index+1
          tv(i1)=var(index)
        enddo
!     
!     determining stress x velocity + conductivity x
!     temperature gradient
!     
        do i1=1,3
          tv(i1)=tv(i1)+cond*dtem(i1)
        enddo
!     
!     determination of the rhs of the energy equations
!     
        do jj=1,nope
          ff(jj)=ff(jj)-xsjmod*(
     &         (shp(4,jj)+dt(konl(jj))*shpv(jj)/2.d0)*enthalpy+
     &         shp(1,jj)*tv(1)+shp(2,jj)*tv(2)+shp(3,jj)*tv(3))
        enddo
!     
!     computation of contribution due to body forces
!     
        if(nbody.ne.0) then
!     
!     retrieving bfv (scalar product of the body force
!     with the velocity)
!     
          index=index+1
          bfv=var(index)
          index=index+9
!     
          do jj=1,nope
            ff(jj)=ff(jj)+xsjmod*(shp(4,jj)+
     &           dt(konl(jj))*shpv(jj)/2.d0)*bfv
          enddo
        else
          index=index+10
        endif
!     
!     distributed heat flux
!     
        if(nload.gt.0) then
          call nident2(nelemload,nelem,nload,id)
          areaj=xsj*weight
          do
            if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
            if(sideload(id)(1:2).ne.'BF') then
              id=id-1
              cycle
            endif
            if(sideload(id)(3:4).eq.'NU') then
              do j=1,3
                pgauss(j)=0.d0
                do i1=1,nope
                  pgauss(j)=pgauss(j)+
     &                 shp(4,i1)*co(j,konl(i1))
                enddo
              enddo
              jltyp=1
              iscale=1
              call dflux(xload(1,id),temp,istep,iinc,tvar,
     &             nelem,kk,pgauss,jltyp,temp,press,sideload(id),
     &             areaj,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,
     &             nmpc,ikmpc,ilmpc,iscale,mi)
              if((nmethod.eq.20).and.(iscale.ne.0))
     &             xload(1,id)=xloadold(1,id)+
     &             (xload(1,id)-xloadold(1,id))*reltimef
            endif
            do jj=1,nope
              ff(jj)=ff(jj)+xsjmod*(shp(4,jj)+
     &             dt(konl(jj))*shpv(jj)/2.d0)*xload(1,id)
            enddo
            exit
          enddo
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
        call nident(nelemface,nelem,nface,idf)
        do
          if((idf.eq.0).or.(nelemface(idf).ne.nelem)) exit
          ig=ichar(sideface(idf)(1:1))-48
!     
!     check for distributed flux
!     an adiabatic face must be declared as a face with
!     distributed flux zero!
!     
          flux=0
          call nident2(nelemload,nelem,nload,id)
          do
            if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
            if((sideload(id)(1:1).ne.'F').and.
     &           (sideload(id)(1:1).ne.'R').and.
     &           (sideload(id)(1:1).ne.'S')) then
              id=id-1
              cycle
            endif
            igl=ichar(sideload(id)(2:2))-48
            if(igl.ne.ig) then
              id=id-1
              cycle
            endif
            flux=1
            exit
          enddo
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
!     calculating of the temperature gradient dtem
!     in the integration point
!     
            do i1=1,3
              dtem(i1)=0.d0
            enddo
            do i1=1,nope
              do j1=1,3
                dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
              enddo
            enddo
!     
!     retrieving the temperature
!     
            index=index+1
            varf(index)=temp
!     
!     material data (conductivity)
!     
            call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &           cond)
!     
!     determining  stress x velocity + conductivity x
!     temperature gradient
!     
            do i1=1,3
              index=index+1
              tv(i1)=varf(index)
              if(flux.eq.0) then
                tv(i1)=tv(i1)+cond*dtem(i1)
              endif
            enddo
            index=index+5
!     
            tvn=tv(1)*xsj2(1)+tv(2)*xsj2(2)+tv(3)*xsj2(3)
!     
            if(flux.eq.1) then
              dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &             xsj2(3)*xsj2(3))
              areaj=dxsj2*weight
              sinktemp=xload(2,id)
!     
!     for nonuniform load: determine the coordinates of the
!     point (transferred into the user subroutine)
!     
              if((sideload(id)(3:4).eq.'NU').or.
     &             (sideload(id)(5:6).eq.'NU')) then
                if(nope.eq.8) then
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifaceq(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                elseif(nope.eq.4) then
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifacet(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                else
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifacew(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                endif
                jltyp=ichar(sideload(id)(2:2))-48
                jltyp=jltyp+10
                if(sideload(id)(1:1).eq.'S') then
                  iscale=1
                  call dflux(xload(1,id),temp,istep,iinc,tvar,
     &                 nelem,i,coords,jltyp,temp,press,
     &                 sideload(id),areaj,vold,co,lakonl,konl,
     &                 ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &                 iscale,mi)
                  if((nmethod.eq.1).and.(iscale.ne.0))
     &                 xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                elseif(sideload(id)(1:1).eq.'F') then
                  call film(xload(1,id),sinktemp,temp,istep,
     &                 iinc,tvar,nelem,i,coords,jltyp,field,
     &                 nfield,sideload(id),node,areaj,vold,mi)
                  if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                elseif(sideload(id)(1:1).eq.'R') then
                  call radiate(xload(1,id),xload(2,id),temp,istep,
     &                 iinc,tvar,nelem,i,coords,jltyp,field,
     &                 nfield,sideload(id),node,areaj,vold,mi,
     &                 iemchange)
                  if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                endif
              endif
!     
              if(sideload(id)(1:1).eq.'S') then
!     
!     flux INTO the face is positive (input deck convention)
!     this is different from the convention in the theory
!     
                tvn=tvn+xload(1,id)*dxsj2
              elseif(sideload(id)(1:1).eq.'F') then
                tvn=tvn-xload(1,id)*(temp-sinktemp)*dxsj2
              elseif(sideload(id)(1:1).eq.'R') then
                tvn=tvn-physcon(2)*
     &               xload(1,id)*((temp-physcon(1))**4-
     &               (xload(2,id)-physcon(1))**4)*dxsj2
              endif
            endif
!     
            xsjmod=tvn*weight*dtime
            do k=1,nopes
              if(nope.eq.8) then
                ipointer=ifaceq(k,ig)
              elseif(nope.eq.4) then
                ipointer=ifacet(k,ig)
              else
                ipointer=ifacew(k,ig)
              endif
              ff(ipointer)=ff(ipointer)+
     &             shp2(4,k)*xsjmod
            enddo
          enddo
          idf=idf-1
        enddo
      endif
!     
      return
      end
