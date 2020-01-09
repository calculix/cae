!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine writemeshinp(kontet,netet_,cotet,nktet,ij)
      !
      implicit none
      !
      character*1 ending
      character*132 fninp
      !
      integer kontet(4,*),netet_,i,j,nsom,nktet,ij
      !
      real*8 cotet(3,*)
      !
      !     storing the mesh in input format
      !
      write(ending,'(i1)') ij
      fninp='finemesh.inp'//ending
      open(2,file=fninp(1:13),status='unknown')
      !
      !     storing the nodes
      !
      write(2,*) '*NODE'
      do i=1,nktet
         write(2,100) i,(cotet(j,i),j=1,3)
      enddo
      !
      !     storing the elements
      !
      write(2,*) '*ELEMENT,TYPE=C3D4,ELSET=TET'
      nsom=0
      do i=1,netet_
         if(kontet(1,i).eq.0) cycle
         nsom=nsom+1
         write(2,101) i,(kontet(j,i),j=1,4)
      enddo
      write(*,*) 'number of tetrahedra= ',nsom
      close(2)
 !
 100  format(i10,',',e20.13,',',e20.13,',',e20.13)
 101  format(5(i10,','))
      !
      return
      end
