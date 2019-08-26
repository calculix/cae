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
      subroutine fill_neiel(nef,ipnei,neiel,neielcp)
      !
      !     copy neiel into neielcp, thereby substituting the zero's by
      !     neighboring values
      !
      implicit none
      !
      integer nef,ipnei(*),neiel(*),neielcp(*) ,i,j,indexf
      !
      intent(in) nef,ipnei,neiel
      !
      intent(inout) neielcp
      !
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            if(neiel(indexf).eq.0) then
               if(indexf.eq.ipnei(i)+1) then
                  do j=ipnei(i)+2,ipnei(i+1)
                     if(neiel(j).ne.0) then
                        neielcp(indexf)=neiel(j)
                        exit
                     endif
                  enddo
               else
                  neielcp(indexf)=neielcp(indexf-1)
               endif
            else
               neielcp(indexf)=neiel(indexf)
            endif
         enddo
      enddo
      !
      return
      end
