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
      subroutine mafillp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,
     &  xlen,cosb,nefa,nefb,iau6,xxicn,flux)
!
!     filling the lhs and rhs to calculate p
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,nef,indexf,ipnei(*),j,neifa(*),nefa,nefb,
     &  neiel(*),iel,ifa,irow(*),ielfa(4,*),compressible,
     &  ifabou(*),neq,jq(*),iel2,indexb,indexf2,
     &  j2,neij(*),nzs,k,iau6(6,*)
!
      real*8 coef,vfa(0:7,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(*),au(*),xle(*),xxn(3,*),ap(*),b(*),cosb(*),
     &  hfa(3,*),gradpel(3,*),bp(*),xxi(3,*),xlen(*),bp_ifa,
     &  xxicn(3,*),flux(*)
!
      do i=nefa,nefb
         indexf=ipnei(i)
         do j=1,ipnei(i+1)-indexf
!     
!     diffusion
! 
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) then
               coef=vfa(5,ifa)*area(ifa)*advfa(ifa)/
     &              (xlet(indexf)*cosb(indexf))
               ad(i)=ad(i)+coef
               if(i.gt.iel) au(iau6(j,i))=au(iau6(j,i))-coef
            else
               if(ielfa(2,ifa).lt.0) then
                  indexb=-ielfa(2,ifa)
                  if(((ifabou(indexb+1).eq.0).or.
     &                (ifabou(indexb+2).eq.0).or.
     &                ( ifabou(indexb+3).eq.0)).and.
     &               (ifabou(indexb+4).ne.0)) then
!     
!                    pressure given (only if not all velocity
!                    components are given)
!     
                     coef=vfa(5,ifa)*area(ifa)*advfa(ifa)/
     &                    (xle(indexf)*cosb(indexf))
                     ad(i)=ad(i)+coef
                  endif
               endif
            endif
!     
!           right hand side: sum of the flux
!     
            b(i)=b(i)-flux(indexf)
         enddo
      enddo
!     
      return
      end
