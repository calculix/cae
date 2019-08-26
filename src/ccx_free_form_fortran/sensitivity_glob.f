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
      subroutine sensitivity_glob(dgdx,dgdxglob,nobject,ndesi,&
        nodedesi,nk)
      !
      !    prepares the sensitivities for the output in the frd-file
      !
      implicit none
      !
      integer nobject,ndesi,nodedesi(*),nk,&
        iobject,node,idesvar
      !
      real*8 dgdx(ndesi,nobject),dgdxglob(2,nk,nobject)
      !
      !     copy the sensitivities in a global node vector
      !
      do idesvar=1,ndesi
         node=nodedesi(idesvar)
         do iobject=1,nobject
            dgdxglob(1,node,iobject)=dgdx(idesvar,iobject)
         enddo
      enddo
      !
      return        
      end
