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
      subroutine rhspcomp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,
     &  xlen,cosb,ielmatf,mi,a1,a2,a3,velo,veloo,dtimef,shcon,
     &  ntmat_,vel,nactdohinv,xrlfa,flux,nefa,nefb,iau6,xxicn,
     &  gamma,xxnj,gradpcfa)
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
      integer i,nef,indexf,ipnei(*),neifa(*),
     &  neiel(*),iel,ifa,irow(*),ielfa(4,*),
     &  ifabou(*),neq,jq(*),indexb,
     &  neij(*),nzs,imat,nefa,nefb,iau6(6,*),
     &  mi(*),ielmatf(mi(3),*),ntmat_,nactdohinv(*)
!
      real*8 coef,vfa(0:7,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(neq),au(nzs),xle(*),xxn(3,*),ap(*),b(neq),cosb(*),
     &  hfa(3,*),gradpel(3,*),bp(*),xxi(3,*),xlen(*),r,xxnj(3,*),
     &  xflux,velo(nef,0:7),veloo(nef,0:7),dtimef,gradpcfa(3,*),
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
!           diffusion (velocity correction)
!
            if(iel.gt.0) then
!
!              internal face
!               
               b(i)=b(i)+vfa(5,ifa)*advfa(ifa)*
     &                   (gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                    gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                    gradpcfa(3,ifa)*xxnj(3,indexf))
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
!                    pressure known
!
                     b(i)=b(i)+vfa(5,ifa)*advfa(ifa)*
     &                    (gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                     gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                     gradpcfa(3,ifa)*xxnj(3,indexf))
                  endif
               endif
            endif
!
!           flux
!
            b(i)=b(i)-xflux
         enddo
!
!        transient term
!
         b(i)=b(i)+(velo(i,5)-vel(i,5))*volume(i)/dtimef
      enddo
!
      return
      end
