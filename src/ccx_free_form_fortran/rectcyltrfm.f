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
      subroutine rectcyltrfm(node,co,cs,icntrl,fin,fout)
      !
      !     transforms a vector fin(1..3) in node "node" with coordinates
      !     in co(1..3,node) from the global rectangular system into a local
      !     cyclic symmetric cylindrical system or vice versa. The result
      !     is stored in fout(1..3)
      !
      implicit none
      !
      integer node,icntrl,i
      !
      real*8 co(3,*),cs(17,*),csab(7),fin(3),fout(3),a(3,3)
      !
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
      !
      if(icntrl.eq.2) then
         !
         !        global rectangular -> cylindrical
         !
         call transformatrix(csab,co(1,node),a)
         !
         fout(1)=fin(1)*a(1,1)+fin(2)*a(2,1)+fin(3)*a(3,1)
         fout(2)=fin(1)*a(1,2)+fin(2)*a(2,2)+fin(3)*a(3,2)
         fout(3)=fin(1)*a(1,3)+fin(2)*a(2,3)+fin(3)*a(3,3)
      !
      elseif(icntrl.eq.-2) then
         !
         !        cylindrical -> global rectangular
         !
         call transformatrix(csab,co(1,node),a)
         !
         fout(1)=fin(1)*a(1,1)+fin(2)*a(1,2)+fin(3)*a(1,3)
         fout(2)=fin(1)*a(2,1)+fin(2)*a(2,2)+fin(3)*a(2,3)
         fout(3)=fin(1)*a(3,1)+fin(2)*a(3,2)+fin(3)*a(3,3)
      !
      endif
      !
      return
      end













