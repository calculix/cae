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
      subroutine postprojectgrad(ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,iconst,objectset,iconstacti,
     &   inameacti)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node,iconst,i,m,iconstacti(*),
     &   inameacti(*)
!
      real*8 dgdxglob(2,nk,nobject),scalar,dd,len
!
!     calculation of final projected gradient
!
      if(nactive.gt.0) then
         do irow=1,ndesi
            node=nodedesi(irow)
            dgdxglob(2,node,nobject)=dgdxglob(2,node,1)
     &          -dgdxglob(2,node,nobject)
         enddo
         objectset(1,nobject)(1:11)='PROJECTGRAD'        
!
         write(*,*)
         write(*,*) '*INFO in postprojectgrad:'
         write(*,*) '      at least 1 constraint active:'
         write(*,*) '      projected gradient has been '
         write(*,*) '      calculated based on the '
         write(*,*) '      constraints:'
         if(nnlconst.eq.nactive) then
            do i=1,nactive
               write(*,'(7x,a12)') objectset(1,ipoacti(i))
            enddo
            write(*,*)
         elseif((nnlconst.lt.nactive).and.(nnlconst.gt.0)) then
            do i=1,nnlconst
               write(*,'(7x,a12)') objectset(1,ipoacti(i))
            enddo
            write(*,'(7x,a12)') objectset(1,inameacti(i))
            write(*,*)
         elseif(nnlconst.eq.0) then
            write(*,'(7x,a12)') objectset(1,inameacti(i))
            write(*,*)
         endif
      else
         write(*,*)
         write(*,*) '*INFO in postprojectgrad:'
         write(*,*) '      no constraint active:'
         write(*,*) '      no projected gradient '
         write(*,*) '      calculated' 
         write(*,*)
      endif    
!
      return        
      end

