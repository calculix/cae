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
      subroutine mafillpcomp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,
     &  xlen,cosb,ielmatf,mi,a1,a2,a3,velo,veloo,dtimef,shcon,
     &  ntmat_,vel,nactdohinv,xrlfa,flux,nefa,nefb,iau6,xxicn,
     &  gamma,inlet)
!
!     filling the lhs and rhs to calculate p
!
!     bc:
!     outflow, p known: diffusion (subsonic)
!              p unknown: convection (supersonic)
!     inflow, p known: none
!             p unknown: convection (subsonic)
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,nef,indexf,ipnei(*),neifa(*),inlet(*),
     &  neiel(*),iel,ifa,irow(*),ielfa(4,*),
     &  ifabou(*),neq,jq(*),indexb,
     &  neij(*),nzs,imat,nefa,nefb,iau6(6,*),
     &  mi(*),ielmatf(mi(3),*),ntmat_,nactdohinv(*)
!
      real*8 coef,vfa(0:7,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(neq),au(nzs),xle(*),xxn(3,*),ap(*),b(neq),cosb(*),
     &  hfa(3,*),gradpel(3,*),bp(*),xxi(3,*),xlen(*),r,
     &  xflux,velo(nef,0:7),veloo(nef,0:7),dtimef,
     &  shcon(0:3,ntmat_,*),vel(nef,0:7),a1,a2,a3,
     &  xrlfa(3,*),gamma(*),flux(*),xxicn(3,*)
!
!
!
      do i=nefa,nefb
         imat=ielmatf(1,i)
         r=shcon(3,1,imat)
         do indexf=ipnei(i)+1,ipnei(i+1)
!
!           loop over all element faces
!
            ifa=neifa(indexf)
            iel=neiel(indexf)
            xflux=flux(indexf)
!
!           convection (density correction)
!
            if(xflux.ge.0.d0) then
!
!              outflowing flux
!
               if(iel.gt.0) then
!
!                 internal face
!
                  ad(i)=ad(i)+xflux/(vfa(5,ifa)*r*vfa(0,ifa))
!
               elseif(ielfa(3,ifa).le.0) then
                  indexb=-ielfa(2,ifa)
                  if(indexb.gt.0) then
                     if((ifabou(indexb+4).eq.0).and.
     &                  (ifabou(indexb+5).eq.0)) then
!
!                       outflow, pressure not known
!
                        ad(i)=ad(i)+xflux/(vfa(5,ifa)*r*vfa(0,ifa))
                     endif
                  else
!
!                    outflow, pressure not known
!
                     ad(i)=ad(i)+xflux/(vfa(5,ifa)*r*vfa(0,ifa))
                  endif
               endif
            else
!
!              inflowing flux
!
               if(iel.gt.0) then
!
!                 internal face
!
                  au(indexf)=au(indexf)+xflux/(vfa(5,ifa)*r*vfa(0,ifa))
               else
!
!                 external face
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if((inlet(ifa).ne.0).and.
     &                  (ifabou(indexb+4).eq.0)) then
!
!                       inlet
!                       and no pressure given (typical subsonic inlet)
!
                        ad(i)=ad(i)+xflux/(vfa(5,ifa)*r*vfa(0,ifa))
                     endif
                  endif
               endif
            endif
!
!           diffusion (velocity correction)
!
            if(iel.gt.0) then
!
!              internal face
!               
               coef=vfa(5,ifa)*area(ifa)*advfa(ifa)/
     &              xlet(indexf)
               ad(i)=ad(i)+coef
               au(indexf)=au(indexf)-coef
            else
!
!              external face
!
               if(ielfa(2,ifa).lt.0) then
                  indexb=-ielfa(2,ifa)
                  if((xflux.ge.0.d0).and.
     &               (ifabou(indexb+4).ne.0)) then
!
!                    not all velocity components known
!                    pressure known (typical subsonic outlet)
!
                     coef=vfa(5,ifa)*area(ifa)*advfa(ifa)/
     &                    xle(indexf)
                     ad(i)=ad(i)+coef
                  endif
               endif
            endif
!
!           flux
!
            b(i)=b(i)-flux(indexf)
         enddo
!
!        transient term
!
         ad(i)=ad(i)+volume(i)/(r*vel(i,0)*dtimef)
         b(i)=b(i)+(velo(i,5)-vel(i,5))*volume(i)/dtimef
      enddo
!
      return
      end
