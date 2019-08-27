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
      subroutine posttransition(dgdxglob,nobject,nk,nodedesi,ndesi,&
         objectset)
      !
      !     Normalizing the sensitivities
      !
      implicit none
      !
      character*81 objectset(4,*)
      !
      integer nobject,nk,nodedesi(*),i,ndesi,m
      !
      real*8 dgdxglob(2,nk,nobject),dd
      !
      !     Scaling the greatest sensitivity value (absolute) to 1
      !
      do m=1,nobject
         if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
         dd=0.d0
         do i=1,ndesi
            dd=max(dd,abs(dgdxglob(2,nodedesi(i),m)))
         enddo
         do i=1,ndesi
            dgdxglob(2,nodedesi(i),m)=dgdxglob(2,nodedesi(i),m)/dd
         enddo
      enddo
      !
      return        
      end




