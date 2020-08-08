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
      subroutine scalesen(dgdxglob,nobject,nk,nodedesi,ndesi,
     &   objectset,iscaleflag,ipoacti,nnlconst,nactive)
!
!     Scaling the sensitivities      
!
!     iscaleflag=1: length of the vector is scaled to 1 
!     iscaleflag=2: one time filtered sensitivities (dgdxglob(2,i,j) are
!                   copied to dgdxglob(1,i,j). Preparation for the second 
!                   filter step   
!     iscaleflag=3: greatest vector value is scaled to 1
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer nobject,nk,nodedesi(*),i,ndesi,m,iscaleflag,node,
     &   ipoacti(*),nnlconst,nactive,istart,iend
!
      real*8 dgdxglob(2,nk,nobject),dd
!
!
      if(iscaleflag.eq.1) then
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            dd=0.d0
            do i=1,ndesi
               node=nodedesi(i)
               dd=dd+dgdxglob(2,node,m)**2
            enddo
            dd=dsqrt(dd)
            do i=1,ndesi
               node=nodedesi(i)
               dgdxglob(2,node,m)=dgdxglob(2,node,m)/dd
            enddo
         enddo
      elseif(iscaleflag.eq.2) then
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle  
            do i=1,ndesi
               node=nodedesi(i)
               dgdxglob(1,node,m)=dgdxglob(2,node,m)
            enddo
         enddo
      elseif(iscaleflag.eq.3) then
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            dd=0.d0
            do i=1,ndesi
               node=nodedesi(i)
               dd=max(dd,abs(dgdxglob(2,node,m)))
            enddo
            do i=1,ndesi
               node=nodedesi(i)
               dgdxglob(2,node,m)=dgdxglob(2,node,m)/dd
               if(objectset(1,m)(1:11).eq.'PROJECTGRAD') then
                  dgdxglob(1,node,m)=dgdxglob(1,node,m)/dd
               endif
            enddo
         enddo
!
!        set all sensitivities of design variables which are linearly 
!        constraint exactly to zero (avoids numeric noise in the 
!        design update) 
!   
         if(nactive.gt.0) then
            istart=nnlconst+1
            iend=nactive
            do m=1,nobject
               if(objectset(1,m)(1:11).eq.'PROJECTGRAD') then
                  do i=istart,iend
                     node=nodedesi(ipoacti(i))
                     dgdxglob(2,node,m)=0.d0
                  enddo
               endif
            enddo
         endif
      endif
!     
      return        
      end
