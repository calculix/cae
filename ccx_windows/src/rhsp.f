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
      subroutine rhsp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,xlen,nefa,nefb,
     &  xxicn,flux,xxnj,gradpcfa,cosb)
!
!     filling the lhs and rhs to calculate the first correction to the
!     pressure p'
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,nef,indexf,ipnei(*),j,neifa(*),
     &  neiel(*),iel,ifa,irow(*),ielfa(4,*),nefa,nefb,
     &  ifabou(*),neq,nzs,jq(*),iel2,indexb,knownflux,
     &  iatleastonepressurebc,j2,indexf2,neij(*)
!
      real*8 coef,vfa(0:7,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(*),au(*),xle(*),xxn(3,*),ap(*),b(*),bp(*),
     &  hfa(3,*),xxi(3,*),gradpel(3,*),xlen(*),bp_ifa,xxicn(3,*),
     &  flux(*),xxnj(3,*),gradpcfa(3,*),cosb(*)
!
      do i=nefa,nefb
         do indexf=ipnei(i)+1,ipnei(i+1)
            knownflux=0
!
!              diffusion
!
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) then
!
!              internal face
!
               b(i)=b(i)+vfa(5,ifa)*area(ifa)*advfa(ifa)*
     &                   (gradpcfa(1,ifa)*xxicn(1,indexf)+
     &                    gradpcfa(2,ifa)*xxicn(2,indexf)+
     &                    gradpcfa(3,ifa)*xxicn(3,indexf))
            else
!
!              external face
!
               iel2=ielfa(2,ifa)
               if(iel2.lt.0) then
                  if(ifabou(-iel2+4).gt.0) then
!     
!                    pressure given
!                        
                     b(i)=b(i)+vfa(5,ifa)*area(ifa)*advfa(ifa)*
     &                    (gradpcfa(1,ifa)*xxicn(1,indexf)+
     &                     gradpcfa(2,ifa)*xxicn(2,indexf)+
     &                     gradpcfa(3,ifa)*xxicn(3,indexf))
                  endif
               endif
            endif
!     
!           right hand side: sum of the flux
!     
            b(i)=b(i)-flux(indexf)
         enddo
      enddo
c         write(*,*) 'rhsp ',1,b(1)
!
      return
      end
