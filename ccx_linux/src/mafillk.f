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
      subroutine mafillk(nef,ipnei,neifa,neiel,vfa,xxn,area,
     &  au,ad,jq,irow,nzs,b,vel,umfa,alet,ale,gradkfa,xxi,
     &  body,volume,ielfa,lakonf,ifabou,nbody,neq,
     &  dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,xrlfa,
     &  xxj,nactdohinv,a1,a2,a3,flux,nefa,nefb,iau6,xxni,xxnj,
     &  iturbulent,f1,of2,yy,umel,gradkel,gradoel,sc)
!
!     filling the matrix for the conservation of energy
!
      implicit none
!
      logical knownflux
!
      character*8 lakonf(*)
!
      integer i,nef,indexf,ipnei(*),j,ifa,iel,neifa(*),
     &  neiel(*),jq(*),irow(*),nzs,ielfa(4,*),nefa,nefb,
     &  ipointer,ifabou(*),nbody,neq,indexb,nactdohinv(*),
     &  iau6(6,*),iturbulent
!
      real*8 xflux,vfa(0:7,*),xxn(3,*),area(*),au(*),ad(*),b(neq),
     &  vel(nef,0:7),umfa(*),alet(*),ale(*),coef,gradkfa(3,*),
     &  xxi(3,*),body(0:3,*),volume(*),dtimef,velo(nef,0:7),
     &  veloo(nef,0:7),rhovol,cvel(*),gradvel(3,3,*),sk,
     &  cvfa(*),hcfa(*),div,xload(2,*),gamma(*),xrlfa(3,*),
     &  xxj(3,*),a1,a2,a3,flux(*),xxnj(3,*),xxni(3,*),difcoef,
     &  f1(*),of2(*),yy(*),umel(*),gradkel(3,*),gradoel(3,*),
     &  cd,arg1,sc(*),constant
!
!
!
      do i=nefa,nefb
!
         if(iturbulent.eq.1) then
!
!           k-epsilon model
!
            sk=1.d0
         elseif(iturbulent.eq.2) then
!
!           k-omega model
!
            sk=0.5d0
         else
!
!           BSL and SST model
!
            cd=max(2.d0*vel(i,5)*0.856d0*
     &             (gradkel(1,i)*gradoel(1,i)+
     &              gradkel(2,i)*gradoel(2,i)+
     &              gradkel(3,i)*gradoel(3,i))/vel(i,7),
     &           1.d-20)
            arg1=min(max(dsqrt(vel(i,6))/(0.09d0*vel(i,7)*yy(i)),
     &                   500.d0*umel(i)/(vel(i,5)*yy(i)**2*vel(i,7))),
     &               4.d0*vel(i,5)*0.856d0*vel(i,6)/(cd*yy(i)**2))
            f1(i)=dtanh(arg1**4)
            if(iturbulent.eq.3) then
!
!              BSL model
!
               sk=0.5d0*f1(i)+(1.d0-f1(i))
            else
!
!              SST model
!
               sk=0.85d0*f1(i)+(1.d0-f1(i))
            endif
         endif
!
         do indexf=ipnei(i)+1,ipnei(i+1)
!
!     convection
!
            ifa=neifa(indexf)
            iel=neiel(indexf)
            xflux=flux(indexf)
!
            if(xflux.ge.0.d0) then
!     
!     outflowing flux
!     
               ad(i)=ad(i)+xflux
!
               b(i)=b(i)-gamma(ifa)*(vfa(6,ifa)-vel(i,6))*xflux
!
            else
               if(iel.gt.0) then
!
!                    incoming flux from neighboring element
!
                  au(indexf)=au(indexf)+xflux
!
                  b(i)=b(i)-gamma(ifa)*(vfa(6,ifa)-vel(iel,6))*xflux
!
               else
!
!                    incoming flux through boundary
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if(((ifabou(indexb+1).ne.0).and.
     &                   (ifabou(indexb+2).ne.0).and.
     &                   (ifabou(indexb+3).ne.0)).or.
     &                    (dabs(xflux).lt.1.d-10)) then
                        b(i)=b(i)-vfa(6,ifa)*xflux
                     endif
                  endif
               endif
            endif
!
!           diffusion
!
            if(iturbulent.le.3) then
!
!              k-epsilon, k-omega or BSL model
!
               difcoef=umfa(ifa)+sk*vfa(5,ifa)*vfa(6,ifa)/vfa(7,ifa)
            else
!
!              SST model
!
               difcoef=umfa(ifa)+sk*vfa(5,ifa)*(0.31d0*vfa(6,ifa))/
     &                           max(0.31d0*vfa(7,ifa),of2(i))
            endif
!
            if(iel.ne.0) then
!     
!              neighboring element
!     
               coef=difcoef*alet(indexf)
               ad(i)=ad(i)+coef
               au(indexf)=au(indexf)-coef
!     
!              correction for non-orthogonal grid
!     
               b(i)=b(i)+difcoef*
     &              (gradkfa(1,ifa)*xxnj(1,indexf)+
     &               gradkfa(2,ifa)*xxnj(2,indexf)+
     &               gradkfa(3,ifa)*xxnj(3,indexf))
            else
!     
!              boundary; either temperature given or adiabatic
!              or outlet
!     
               knownflux=.false.
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  if((ifabou(ipointer+5).ne.0).or.
     &               (ifabou(ipointer+1).ne.0).or.
     &               (ifabou(ipointer+2).ne.0).or.
     &               (ifabou(ipointer+3).ne.0)) then
!     
!                    no outlet:
!                    (i.e. no wall || no sliding || at least one velocity given)
!                    turbulent variable is assumed fixed
!     
                     coef=difcoef*ale(indexf)
                     ad(i)=ad(i)+coef
                     b(i)=b(i)+coef*vfa(6,ifa)
                  else
!     
!                     outlet: no diffusion
!     
                  endif
               endif
!     
!              correction for non-orthogonal grid
!     
               if(.not.knownflux) then
                  b(i)=b(i)+difcoef*
     &                 (gradkfa(1,ifa)*xxnj(1,indexf)+
     &                  gradkfa(2,ifa)*xxnj(2,indexf)+
     &                  gradkfa(3,ifa)*xxnj(3,indexf))
               endif
            endif
         enddo
!
!           viscous dissipation
!     
         rhovol=vel(i,5)*volume(i)
!
!        sink terms are treated implicitly (lhs)
!
         ad(i)=ad(i)+rhovol*0.09d0*vel(i,7)
!
!        source terms are treated explicitly (rhs)
!
         if(iturbulent.le.3) then
!
!           k-epsilon, k-omega and BSL-model: the turbulent
!           viscosity=k/omega
!
            b(i)=b(i)+rhovol*vel(i,6)*((
     &           (2.d0*(gradvel(1,1,i)**2+gradvel(2,2,i)**2+
     &                  gradvel(3,3,i)**2)+
     &           (gradvel(1,2,i)+gradvel(2,1,i))**2+
     &           (gradvel(1,3,i)+gradvel(3,1,i))**2+
     &           (gradvel(2,3,i)+gradvel(3,2,i))**2))/vel(i,7))
         else
!
!           SST model: other definition of the turbulent
!           viscosity
!
            b(i)=b(i)+rhovol*0.31d0*vel(i,6)*((
     &           (2.d0*(gradvel(1,1,i)**2+gradvel(2,2,i)**2+
     &                  gradvel(3,3,i)**2)+
     &           (gradvel(1,2,i)+gradvel(2,1,i))**2+
     &           (gradvel(1,3,i)+gradvel(3,1,i))**2+
     &           (gradvel(2,3,i)+gradvel(3,2,i))**2))/
     &           max(0.31d0*vel(i,7),of2(i)))
         endif
!
!           transient term
!
         constant=rhovol*sc(i)
         b(i)=b(i)-(a2*velo(i,6)+a3*veloo(i,6))*constant
         constant=a1*constant
         ad(i)=ad(i)+constant
!     
      enddo
!     
      return
      end
