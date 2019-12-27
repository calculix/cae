!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
      subroutine reorderlhs(au,am,iamorig,nz_numa,nz_numb)
      !
      !     reorders matrix elements into compressed row format
      !
      implicit none
      !
      integer iamorig(*),nz_numa,nz_numb,i
      !
      real*8 au(*),am(*)
      !
      intent(in) au,iamorig,nz_numa,nz_numb
      !
      intent(inout) am
      !
      do i=nz_numa,nz_numb
         am(i)=au(iamorig(i))
      enddo
      !
      return
      end
