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
      subroutine writerandomfield(d,nev,abserr,relerr)
      !
      !     writes the eigenvalues and the error measures of
      !     the randomfield in the .dat file
      !
      implicit none
      !
      integer iobject,i,nev
      real*8 d(*),abserr,relerr
      !
      write(5,*)
      write(5,*) 'EIGENVALUES OF MODESHAPES OF RANDOMFIELD'  
      write(5,*)
      do i=nev,1,-1
         write(5,'(7x,e14.7)') d(i)
      enddo
      !
      write(5,*)
      write(5,*) 'ABSOLUTE ERROR W.R.T. THE VARIANCE OF THE RANDOMFIELD'
      write(5,*)
      write(5,'(7x,e14.7)') abserr
      write(5,*)
      write(5,*) 'RELATIVE ERROR W.R.T. THE VARIANCE OF THE RANDOMFIELD'
      write(5,*)
      write(5,'(7x,e14.7)') relerr

      return
      end

