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
!     y=A*x for real sparse symmetric matrices
!
!     storage of the matrix:
!        au: first lower triangle
!        ad: diagonal terms
!
      subroutine preconvert2slapcol(irow,ia,jq,ja,nzs,nef)
!
      implicit none
!
      integer irow(*),ia(*),jq(*),ja(*),nzs,nef,i,j,k
!
!     converting the CalculiX format into the SLAP column format
!
      k=nzs+nef
!
      do i=nef,1,-1
         do j=jq(i+1)-1,jq(i),-1
            ia(k)=irow(j)
            k=k-1
         enddo
         ia(k)=i
         k=k-1
      enddo
!
      ja(1)=1
      do i=2,nef+1
         ja(i)=jq(i)+i-1
      enddo
!
      return
      end
