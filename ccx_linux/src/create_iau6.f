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
      subroutine create_iau6(nef,ipnei,neiel,jq,irow,nzs,iau6,lakonf)
!
!     sets up a field of pointers iau6(j,i) for neighbor j of
!     element (cell) i (for CFD-applications) into field auv,aup,aut
!
      implicit none
!
      character*8 lakonf(*)
!
      integer nef,ipnei(*),neiel(*),jq(*),irow(*),nzs,iau6(6,*),
     &  numfaces,id,i,j,iel,indexf
!
!
!
      do i=1,nef
         indexf=ipnei(i)
!
         do j=1,ipnei(i+1)-ipnei(i)
            indexf=indexf+1
            iel=neiel(indexf)
            if(iel.eq.0) cycle
            if(i.gt.iel) then
               call nident(irow(jq(iel)),i,jq(iel+1)-jq(iel),id)
               iau6(j,i)=jq(iel)+id-1
            elseif(i.lt.iel) then
               call nident(irow(jq(i)),iel,jq(i+1)-jq(i),id)
               iau6(j,i)=nzs+jq(i)+id-1
            endif
         enddo
      enddo
!     
      return
      end
