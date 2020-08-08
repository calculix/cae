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
      subroutine reorderrhs(a,b,vold,neighblock,nneighblock)
!
!     reorders matrix elements into compressed row format
!
      implicit none
!
      integer neighblock(3,*),nneighblock,i,j,indexf,iel
!
      real*8 a(*),b(*),vold(*)
!
!
!
      do j=1,nneighblock
!
!        location in au/auv
!
         indexf=neighblock(1,j)
!
!        neighboring block element number
!
         iel=neighblock(2,j)
!
!        equation number
!
         i=neighblock(3,j)
!
         b(i)=b(i)-a(indexf)*vold(iel)
!
      enddo
!
      return
      end
