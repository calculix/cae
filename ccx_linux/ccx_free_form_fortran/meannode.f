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
      subroutine meannode(nk,inum,v)
      !
      !     taking the mean at the nodes
      !
      implicit none
      !
      integer nk,inum(*),i,j
      !
      real*8 v(0:4,*)
      !
      do i=1,nk
         if(inum(i).ne.0) then
            inum(i)=abs(inum(i))
            do j=0,4
               v(j,i)=v(j,i)/inum(i)
            enddo
         endif
         write(*,*) 'meannode ',i,(v(j,i),j=1,3)
      enddo
      !
      return
      end
