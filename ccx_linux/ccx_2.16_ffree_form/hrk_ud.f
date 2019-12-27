!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
      subroutine hrk_ud(ielfa,vel,ipnei,nef,flux,vfa,nfacea,nfaceb)
      !
      !     determine the facial temperature using upwind difference
      !
      implicit none
      !
      integer ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef,nfacea,&
        nfaceb
      !
      real*8 vel(nef,0:7),flux(*),vfa(0:7,*)
      !
      intent(in) ielfa,vel,ipnei,nef,flux,nfacea,nfaceb
      !
      intent(inout) vfa
      !
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
         !
         !        faces with only one neighbor need not be treated
         !        unless outlet
         !
         !          if((iel2.le.0).and.(ielfa(3,i).ge.0)) cycle
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
         !
         if(flux(indexf).ge.0.d0) then
            !
            !           outflow && (neighbor || outlet)
            !
            vfa(6,i)=vel(iel1,6)
         elseif(iel2.gt.0) then
            !
            !           inflow && neighbor
            !
            vfa(6,i)=vel(iel2,6)
         endif
      enddo
      !
      return
      end
