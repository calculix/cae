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
      subroutine openfilefluid(jobname)
!
      implicit none
!
      character*132 jobname,fnfrd
      integer i
!
!     opening frd file
!
      do i=1,132
         if(jobname(i:i).eq.' ') exit
      enddo
      i=i-1
      if(i.gt.128) then
         write(*,*) 
     &     '*ERROR in openfilefluid: input file name is too long:'
         write(*,'(a132)') jobname(1:132)
         write(*,*) '       exceeds 128 characters'
         call exit(201)
      endif
!
      fnfrd=jobname(1:i)//'.frd'
      open(13,file=fnfrd(1:i+5),status='unknown',position='append')
!
      return
      end
