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
      subroutine mafillt(nef,ipnei,neifa,neiel,vfa,xxn,area,
     &  au,ad,jq,irow,nzs,b,vel,umel,alet,ale,gradtfa,xxi,
     &  body,volume,ielfa,lakonf,ifabou,nbody,neq,
     &  dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,xrlfa,
     &  xxj,nactdohinv,a1,a2,a3,flux,nefa,nefb,iau6,xxni,xxnj,
     &  iturbulent,of2,sc)
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
     &  vel(nef,0:7),umel(*),alet(*),ale(*),coef,gradtfa(3,*),
     &  xxi(3,*),body(0:3,*),volume(*),dtimef,velo(nef,0:7),
     &  veloo(nef,0:7),rhovol,constant,cvel(*),gradvel(3,3,*),
     &  cvfa(*),hcfa(*),div,xload(2,*),gamma(*),xrlfa(3,*),
     &  xxj(3,*),a1,a2,a3,flux(*),xxnj(3,*),xxni(3,*),visdis,
     &  four,of2(*),sc(*)
!
!
!
      do i=nefa,nefb
         do indexf=ipnei(i)+1,ipnei(i+1)
!
!     convection
!
            ifa=neifa(indexf)
            iel=neiel(indexf)
            xflux=flux(indexf)*cvfa(ifa)
!
            if(xflux.ge.0.d0) then
!     
!     outflowing flux
!     
               ad(i)=ad(i)+xflux
!
               b(i)=b(i)-gamma(ifa)*(vfa(0,ifa)-vel(i,0))*xflux
!
            else
               if(iel.gt.0) then
!
!                    incoming flux from neighboring element
!
                  au(indexf)=au(indexf)+xflux
!
                  b(i)=b(i)-gamma(ifa)*(vfa(0,ifa)-vel(iel,0))*xflux
!
               else
!
!                    incoming flux through boundary
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if((ifabou(indexb).ne.0).or.
     &                    (dabs(xflux).lt.1.d-10)) then
                        b(i)=b(i)-vfa(0,ifa)*xflux
                     endif
                  endif
               endif
            endif
!
!           diffusion
!
!           k_turbulent=cp*mu_turbulent/Prandl_turbulent
!           Prandl_turbulent=0.9
!     
            if(iturbulent.eq.0) then
               four=hcfa(ifa)
            elseif(iturbulent.le.3) then
!
!              k-epsilon, k-omega or BSL model
!
               four=hcfa(ifa)+cvfa(ifa)*vfa(5,ifa)*vfa(6,ifa)/
     &                        (vfa(7,ifa)*0.9d0)
            else
!
!              SST model
!
               four=hcfa(ifa)+cvfa(ifa)*vfa(5,ifa)*
     &              (0.31d0*vfa(6,ifa))/
     &              (max(0.31d0*vfa(7,ifa),of2(i))*0.9d0)
            endif
!
            if(iel.ne.0) then
!     
!              neighboring element
!     
               coef=four*alet(indexf)
               ad(i)=ad(i)+coef
               au(indexf)=au(indexf)-coef
!     
!              correction for non-orthogonal grid
!     
               b(i)=b(i)+four*
     &              (gradtfa(1,ifa)*xxnj(1,indexf)+
     &               gradtfa(2,ifa)*xxnj(2,indexf)+
     &               gradtfa(3,ifa)*xxnj(3,indexf))
            else
!     
!              boundary; either temperature given or adiabatic
!              or outlet
!     
               knownflux=.false.
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  if(ifabou(ipointer+6).gt.0) then
!     
!     heat flux is known
!     
                     b(i)=b(i)-xload(1,ifabou(ipointer+6))
                     knownflux=.true.
!     
c                  elseif((ifabou(ipointer).ne.0).or.
c     &                    (ifabou(ipointer+1).ne.0).or.
c     &                    (ifabou(ipointer+2).ne.0).or.
c     &                    (ifabou(ipointer+3).ne.0)) then
                  elseif(ifabou(ipointer).ne.0) then
!     
!                    temperature given or no outlet:
!                    temperature is assumed fixed
!     
                     coef=four*ale(indexf)
                     ad(i)=ad(i)+coef
                     b(i)=b(i)+coef*vfa(0,ifa)
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
                  b(i)=b(i)+four*
     &                 (gradtfa(1,ifa)*xxnj(1,indexf)+
     &                  gradtfa(2,ifa)*xxnj(2,indexf)+
     &                  gradtfa(3,ifa)*xxnj(3,indexf))
               endif
            endif
         enddo
!
!           viscous dissipation
!
         b(i)=b(i)+umel(i)*volume(i)*
     &        (2.d0*(gradvel(1,1,i)**2+gradvel(2,2,i)**2+
     &        gradvel(3,3,i)**2)+
     &        (gradvel(1,2,i)+gradvel(2,1,i))**2+
     &        (gradvel(1,3,i)+gradvel(3,1,i))**2+
     &        (gradvel(2,3,i)+gradvel(3,2,i))**2)
!     
!           body heat source and body sources
!     
         rhovol=vel(i,5)*volume(i)
!
         if(nbody.gt.0) then
            b(i)=b(i)+rhovol*body(0,i)
         endif
!
!           transient term
!
         constant=rhovol*cvel(i)*sc(i)
         b(i)=b(i)-(a2*velo(i,0)+a3*veloo(i,0))*constant
         constant=a1*constant
         ad(i)=ad(i)+constant
!     
      enddo
!     
      return
      end
