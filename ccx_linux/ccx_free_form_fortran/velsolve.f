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
      subroutine velsolve(nef,ipnei,bv,auv,adv,vel,temp,neiel)
      !
      !
      implicit none
      !
      integer nef,ipnei(*),i,j,indexf,neiel(*),iel
      !
      real*8 bv(nef,3),auv(*),adv(*),vel(nef,0:7),temp(nef,0:7)
      !
      do i=1,nef
         do j=1,3
            temp(i,j)=bv(i,j)
         enddo
         do indexf=ipnei(i)+1,ipnei(i+1)
            iel=neiel(indexf)
            do j=1,3
               temp(i,j)=temp(i,j)-auv(indexf)*vel(iel,j)
            enddo
         enddo
         do j=1,3
            temp(i,j)=temp(i,j)/adv(i)
         enddo
      enddo
      !
      return
      end
