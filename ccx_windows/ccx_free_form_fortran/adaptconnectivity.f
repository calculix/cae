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
      subroutine adaptconnectivity(konl1,iface,konl2,konl2f,&
        neighbor)
      !
      !     the 8-node volumetric elements with topologies konl1 and
      !     konl2 have the face iface in common. iface is the local face
      !     number for topology konl1. A new topology for element 2
      !     is looked for, which is such that the local xi, eta and
      !     zeta axes point in the same direction for both elements
      !
      implicit none
      !
      integer konl1(8),konl2(8),iface,konl2f(8),ifaceq(8,6),&
        ifaceqtrans(4,6),konf1l(4),konf1r(4),konf2l(4),konf2r(4),&
        node,i,j,i1,i2,i3,neighbor(8,8,8)
      !
      intent(in) konl1,iface,konl2,neighbor
      !
      intent(inout) konl2f
      !
      !     facial nodes
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
                  5,6,7,8,13,14,15,16,&
                  1,2,6,5,9,18,13,17,&
                  2,3,7,6,10,19,14,18,&
                  3,4,8,7,11,20,15,19,&
                  4,1,5,8,12,17,16,20/
      !
      !     opposite (translated) nodes to ifaceq
      !
      data ifaceqtrans /8,7,6,5,&
                        1,2,3,4,&
                        4,3,7,8,&
                        1,4,8,5,&
                        2,1,5,6,&
                        3,2,6,7/
      !
      do i=1,4
         konf1l(i)=ifaceqtrans(i,iface)
         konf1r(i)=ifaceq(i,iface)
      enddo
      !
      do i=1,4
         node=konl1(konf1r(i))
         do j=1,8
            if(konl2(j).eq.node) then
               konf2l(i)=j
               exit
            endif
         enddo
      enddo
      !
      do i=1,4
         i2=konf2l(i)
         if(i.eq.1) then
            i1=konf2l(4)
         else
            i1=konf2l(i-1)
         endif
         if(i.eq.4) then
            i3=konf2l(1)
         else
            i3=konf2l(i+1)
         endif
         konf2r(i)=neighbor(i1,i2,i3)
      enddo
      !
      do i=1,4
         konl2f(konf1l(i))=konl2(konf2l(i))
         konl2f(konf1r(i))=konl2(konf2r(i))
      enddo
      !
      return
      end
