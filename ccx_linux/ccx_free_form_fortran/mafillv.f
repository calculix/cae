!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
      subroutine mafillv(nef,ipnei,neifa,neiel,vfa,xxn,area,&
        auv,adv,jq,irow,nzs,bv,vel,cosa,umfa,alet,ale,gradvfa,xxi,&
        body,volume,ielfa,lakonf,ifabou,nbody,&
        dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,a1,&
        a2,a3,flux,nefa,nefb,icyclic,c,ifatie,iau6,xxna,xxnj,&
        iturbulent,gradvel,of2,yy,umel,ncfd,inlet,sc)
      !
      implicit none
      !
      character*8 lakonf(*)
      !
      integer i,nef,indexf,ipnei(*),j,ifa,iel,neifa(*),icyclic,&
        neiel(*),jq(*),irow(*),nzs,iwall,ielfa(4,*),nefa,nefb,&
        ipointer,ifabou(*),nbody,k,indexb,nactdohinv(*),&
        ifatie(*),iau6(6,*),iturbulent,ncfd,inlet(*)
      !
      real*8 xflux,vfa(0:7,*),xxn(3,*),area(*),auv(*),adv(*),bv(nef,3),&
        vel(nef,0:7),cosa(*),umfa(*),alet(*),ale(*),coef,gradvfa(3,3,*),&
        xxi(3,*),body(0:3,*),volume(*),coef2,dtimef,velo(nef,0:7),&
        veloo(nef,0:7),rhovol,constant,sel(3,*),xrlfa(3,*),gamma(*),&
        xxj(3,*),a1,a2,a3,flux(*),c(3,3),xxna(3,*),xxnj(3,*),difcoef,&
        xl1,xl2,aa,bb,gradvel(3,3,*),of2(*),f2,arg2,yy(*),umel(*),sc(*)
      !
      intent(in) nef,ipnei,neifa,neiel,vfa,xxn,area,&
        jq,irow,nzs,vel,cosa,umfa,alet,ale,gradvfa,xxi,&
        body,volume,ielfa,lakonf,ifabou,nbody,&
        dtimef,velo,veloo,xrlfa,gamma,xxj,nactdohinv,a1,&
        a2,a3,flux,gradvel,yy,umel,iturbulent,ncfd,sc
      !
      intent(inout) adv,auv,bv,sel,of2
      !
      do i=nefa,nefb
         do indexf=ipnei(i)+1,ipnei(i+1)
            !
            !              convection
            !
            ifa=neifa(indexf)
            iel=neiel(indexf)
            xflux=flux(indexf)
            !
            if(xflux.ge.0.d0) then
               !
               !                 outflowing xflux
               !
               if(iel.eq.0) then
                  adv(i)=adv(i)+xflux
               else
                  adv(i)=adv(i)+xflux
                  do k=1,ncfd
                    !   retarded gamma
                    bv(i,k)=bv(i,k)-gamma(ifa)*(vfa(k,ifa)&
                             -vel(i,k))*xflux
                  enddo
               endif
            else
               if(iel.gt.0) then
                  !
                  !                    incoming flux from neighboring element
                  !
                  if((icyclic.eq.0).or.(ifatie(ifa).eq.0)) then
                     auv(indexf)=auv(indexf)+xflux
                     do k=1,ncfd
                        !    retarded gamma
                        bv(i,k)=bv(i,k)-gamma(ifa)*&
                                (vfa(k,ifa)-vel(iel,k))*xflux
                     enddo
                  else
                     do k=1,ncfd
                        bv(i,k)=bv(i,k)-vfa(k,ifa)*xflux
                     enddo
                  endif
               else
                  !
                  !                    incoming flux through boundary
                  !
                  if(ielfa(2,ifa).lt.0) then
                     if(inlet(ifa).ne.0) then
                        do k=1,ncfd
                           bv(i,k)=bv(i,k)-vfa(k,ifa)*xflux
                        enddo
                     endif
                  endif
               endif
            endif
            !
            !              diffusion (laminar + turbulent)
            !
            if(iturbulent.eq.0) then
               difcoef=umfa(ifa)
            elseif(iturbulent.le.3) then
               !
               !              k-epsilon, k-omega or BSL model
               !
               difcoef=umfa(ifa)+vfa(5,ifa)*vfa(6,ifa)/vfa(7,ifa)
            else
               !
               !              SST model
               !
               arg2=max(2.d0*dsqrt(vel(i,6))/(0.09d0*vel(i,7)*yy(i)),&
                        500.d0*umel(i)/(vel(i,5)*yy(i)**2*vel(i,7)))
               f2=dtanh(arg2**2)
               of2(i)=dsqrt((gradvel(3,2,i)-gradvel(2,3,i))**2+&
                            (gradvel(1,3,i)-gradvel(3,1,i))**2+&
                            (gradvel(2,1,i)-gradvel(1,2,i))**2)*f2
               difcoef=umfa(ifa)+vfa(5,ifa)*0.31d0*vfa(6,ifa)/&
                       max(0.31d0*vfa(7,ifa),of2(i))
            endif
            !
            if(iel.ne.0) then
               !
               !                 neighboring element
               !
               coef=difcoef*alet(indexf)
               adv(i)=adv(i)+coef
               if((icyclic.eq.0).or.(ifatie(ifa).eq.0)) then
                  auv(indexf)=auv(indexf)-coef
               elseif(ifatie(ifa).gt.0) then
                  !
                  !                 for cyclic symmetry the term is retarded, since
                  !                 otherwise the x-, y- and z- components are linked
                  !                 (i.e. the x-, y- and z- momentum equations cannot
                  !                  be solved separately any more)
                  !
                  do k=1,ncfd
                     bv(i,k)=bv(i,k)+&
                      (c(k,1)*vel(iel,1)+c(k,2)*vel(iel,2)+&
                       c(k,3)*vel(iel,3))*coef
                  enddo
               else
                  do k=1,ncfd
                     bv(i,k)=bv(i,k)+&
                      (c(1,k)*vel(iel,1)+c(2,k)*vel(iel,2)+&
                       c(3,k)*vel(iel,3))*coef
                  enddo
               endif
               !
               !                 correction for non-orthogonal grid
               !
               do k=1,ncfd
                  bv(i,k)=bv(i,k)+difcoef*&
                    (gradvfa(k,1,ifa)*xxnj(1,indexf)+&
                     gradvfa(k,2,ifa)*xxnj(2,indexf)+&
                     gradvfa(k,3,ifa)*xxnj(3,indexf))
               enddo
            else
               !
               !                 boundary; check whether wall (specified by user),
               !                           outlet (no velocity boundary conditions) or
               !                           none of those
               !
               iwall=0
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  iwall=ifabou(ipointer+5)
               endif
               if(iwall.eq.0) then
                  !
                  !                    external face, but no wall
                  !
                  if(ielfa(3,ifa).gt.0) then
                     !
                     !                       no outlet: face velocity fixed
                     !
                     coef=difcoef*ale(indexf)
                     adv(i)=adv(i)+coef
                     do k=1,ncfd
                        bv(i,k)=bv(i,k)+coef*vfa(k,ifa)
                     enddo
                  else
                  !
                  !                       outlet: no diffusion
                  !
                  endif
                  !
                  !                    correction for non-orthogonal grid
                  !
                  do k=1,ncfd
                     bv(i,k)=bv(i,k)+difcoef*&
                       (gradvfa(k,1,ifa)*xxnj(1,indexf)+&
                        gradvfa(k,2,ifa)*xxnj(2,indexf)+&
                        gradvfa(k,3,ifa)*xxnj(3,indexf))
                  enddo
               elseif(iwall.gt.0) then
                  !
                  !                    wall
                  !
                  coef=difcoef*cosa(indexf)
                  adv(i)=adv(i)+coef
                  !
                  !                    correction for non-orthogonal grid and nonzero
                  !                    wall velocity
                  !
                  coef2=(vel(i,1)*xxn(1,indexf)+&
                       vel(i,2)*xxn(2,indexf)+&
                       vel(i,3)*xxn(3,indexf))*coef
                  do k=1,ncfd
                     bv(i,k)=bv(i,k)+coef*vfa(k,ifa)+&
                       coef2*xxn(k,indexf)
                  enddo
               endif
            endif
            !
            !           pressure and turbulent kinetic energy
            !
            if(iturbulent.eq.0) then
               do k=1,ncfd
                  bv(i,k)=bv(i,k)&
                       -vfa(4,ifa)*xxna(k,indexf)
               enddo
            else
               do k=1,ncfd
                  bv(i,k)=bv(i,k)&
                       -(vfa(4,ifa)+2.d0*vfa(5,ifa)*vfa(6,ifa)/3.d0)&
                       *xxna(k,indexf)
               enddo
            endif
         enddo
         !
         !        body force
         !
         rhovol=vel(i,5)*volume(i)
         !
         if(nbody.gt.0) then
            do k=1,ncfd
               bv(i,k)=bv(i,k)+rhovol*body(k,i)
            enddo
         endif
         !
         !           transient term
         !
         constant=rhovol*sc(i)
         do k=1,ncfd
            bv(i,k)=bv(i,k)-(a2*velo(i,k)+a3*veloo(i,k))*constant
         enddo
         constant=a1*constant
         adv(i)=adv(i)+constant
      !
      !        pressure contribution to b
      !
      !          if(iturbulent.eq.0) then
      !             do indexf=ipnei(i)+1,ipnei(i+1)
      !                ifa=neifa(indexf)
      !                do k=1,ncfd
      !                   bv(i,k)=bv(i,k)
      !      &                 -vfa(4,ifa)*xxna(k,indexf)
      !                enddo
      !             enddo
      !          else
      !             do indexf=ipnei(i)+1,ipnei(i+1)
      !                ifa=neifa(indexf)
      !                do k=1,ncfd
      !                   bv(i,k)=bv(i,k)
      !      &               -(vfa(4,ifa)+2.d0*vfa(5,ifa)*vfa(6,ifa)/3.d0)
      !      &               *xxna(k,indexf)
      !                enddo
      !             enddo
      !          endif
      !
      enddo
      !
      return
      end
