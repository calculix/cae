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
      subroutine identifyface(konl1,konl2,iface1,iface2)
      !
      !     the topology of one and the same 8-node element is described
      !     in two different ways: konl1 and konl2. This routine identifies
      !     the internal face number iface2 of konl2 which corresponds to face
      !     number iface1 of konl1
      !
      implicit none
      !
      integer konl1(8),konl2(8),iface1,iface2,konl(3),i,j,kode,node,&
        ifaceq(8,6)
      !
      intent(in) konl1,konl2,iface1
      !
      intent(inout) iface2
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
                  5,6,7,8,13,14,15,16,&
                  1,2,6,5,9,18,13,17,&
                  2,3,7,6,10,19,14,18,&
                  3,4,8,7,11,20,15,19,&
                  4,1,5,8,12,17,16,20/
      !
      do j=1,3
         konl(j)=0
      enddo
      !
      !     the face is identified using the first three entries
      !     of topology konl2
      !
      do i=1,4
         node=konl1(ifaceq(i,iface1))
         do j=1,3
            if(konl2(j).eq.node) then
               konl(j)=1
               exit
            endif
         enddo
      enddo
      !
      kode=4*konl(1)+2*konl(2)+konl(3)
      if(kode.eq.7) then
         iface2=1
      elseif(kode.eq.0) then
         iface2=2
      elseif(kode.eq.6) then
         iface2=3
      elseif(kode.eq.3) then
         iface2=4
      elseif(kode.eq.1) then
         iface2=5
      else
         iface2=6
      endif
      !
      return
      end
