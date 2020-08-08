!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine extrapol_pel5(ielfa,gradpfa,nfacea,nfaceb,ncfd)
!     
!     taking the mean of the facial pressure gradient for faces in
!     between two elements
!     
      implicit none
!     
      integer ielfa(4,*),nfacea,nfaceb,i,l,iel2,ncfd
!     
      real*8 gradpfa(3,*)
!     
!     
!     
      do i=nfacea,nfaceb
        iel2=ielfa(2,i)
        if(iel2.gt.0) then
            do l=1,ncfd
              gradpfa(l,i)=gradpfa(l,i)/2.d0
            enddo
        endif
      enddo
!     
      return
      end
