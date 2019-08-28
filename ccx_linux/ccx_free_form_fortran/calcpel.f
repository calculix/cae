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
      subroutine calcpel(ne,nactdoh,vel,b,nef)
      !
      !     stores the pressure into field vel
      !
      implicit none
      !
      integer ne,nactdoh(*),i,nef
      !
      real*8 vel(nef,0:7),b(*)
      !
      do i=1,ne
         vel(i,4)=vel(i,4)+b(i)
      !          write(*,*) i
      !          write(*,*) nactdoh(i)
      !          write(*,*) b(nactdoh(i))
      !          write(*,*) 'calcpel ',i,vel(i,4)
      enddo
      !       write(*,*) 'calcpel ',vel(1600,4)
      !       write(*,*)
      !
      return
      end
