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
      subroutine applympc_dpel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
     &  xle,xxi,icyclic,xxn,ipnei,ifatie,
     &  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
     &  iflag,xxj,xlet)
!
!     applying the MPC's in the calculation of the pressure
!     correction
!
      implicit none
!
      character*20 labmpc(*)
!
      integer nface,ielfa(4,*),ifabou(*),i,iel1,iel2,nef,ibou,
     &  neifa(*),icyclic,ifa,indexf,l,m,ipnei(*),ifatie(*),
     &  is,ie,nmpc,ipompc(*),nodempc(3,*),ifaext(*),nfaext,nactdoh(*),
     &  iel,index,mpc,ipointer,k,ielorig,iface,iflag
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xbounact(*),xl1,xl2,
     &   vfap(0:7,nface),gradpcel(3,*),gradpcfa(3,*),rf(3,*),area(*),
     &   volume(*),xle(*),xxi(3,*),c(3,3),gradnor,xxn(3,*),coefmpc(*),
     &   coefnorm,sum,xxj(3,*),xlet(*),dd
!
!
!
!     Multiple point constraints
!
      if(nmpc.gt.0) then
         do k=1,nfaext
            i=ifaext(k)
            if(ielfa(2,i).ge.0) cycle
            ipointer=-ielfa(2,i)
            if(ifabou(ipointer+4).ge.0) cycle
            mpc=-ifabou(ipointer+4)
!     
            index=ipompc(mpc)
            sum=0.d0
            coefnorm=0.d0
            do
               if(index.eq.0) exit
               if(nodempc(1,index).lt.0) then
!
!                 a negative number refers to a boundary
!                 condition (fields nodeboun, ndirboun..)
!                 resulting from a SPC in local coordinates
!                  
                  sum=sum+0.d0
               else
!
!                 face term
!
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  sum=sum+coefmpc(index)
     &                 *vfa(nodempc(2,index),neifa(ipnei(iel)+iface))
                  coefnorm=coefnorm+coefmpc(index)**2
               endif
               index=nodempc(3,index)
            enddo
!
!           distribute the sum across all terms which are not
!           fixed by boundary conditions
!
            index=ipompc(mpc)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).gt.0) then
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  vfa(nodempc(2,index),neifa(ipnei(iel)+iface))=
     &                 vfa(nodempc(2,index),neifa(ipnei(iel)+iface))-
     &                 sum*coefmpc(index)/coefnorm
               endif
               index=nodempc(3,index)
            enddo
         enddo
      endif
!
      return
      end
