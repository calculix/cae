!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine applympc_hfa(nface,ielfa,is,ie,ifabou,ipompc,hfa,&
        coefmpc,nodempc,ipnei,neifa,labmpc,xbounact,nactdoh)
      !
      !     applies MPC's to the faces
      !
      implicit none
      !
      character*20 labmpc(*)
      !
      integer i,j,nface,ielfa(4,*),ipointer,is,ie,ifabou(*),mpc,&
        ipompc(*),index,iel,iface,nodempc(3,*),ipnei(*),neifa(*),&
        nactdoh(*),ielorig
      !
      real*8 coefmpc(*),denominator,hfa(3,*),sum,xbounact(*)
      !
      do i=1,nface
         if(ielfa(2,i).ge.0) cycle
         ipointer=-ielfa(2,i)
         do j=is,ie
            if(ifabou(ipointer+j).ge.0) cycle
            mpc=-ifabou(ipointer+j)
            index=ipompc(mpc)
            denominator=coefmpc(index)
            sum=0.d0
            index=nodempc(3,index)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).lt.0) then
                  !
                  !                 a negative number refers to a boundary
                  !                 condition (fields nodeboun, ndirboun..)
                  !                 resulting from a SPC in local coordinates
                  !
                  sum=sum+coefmpc(index)*xbounact(-nodempc(1,index))
               else
                  !
                  !                 face term
                  !
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  sum=sum+coefmpc(index)&
                       *hfa(nodempc(2,index),neifa(ipnei(iel)+iface))
               endif
               index=nodempc(3,index)
            enddo
            hfa(j,i)=-sum/denominator
         enddo
      enddo
      !
      return
      end
