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
      subroutine scalesen(dgdxglob,nobject,nk,nodedesi,ndesi,&
         objectset,iscaleflag)
      !
      !     Scaling the sensitivities
      !
      !     iscaleflag=0: greatest vector value is scaled to 1
      !     iscaleflag=1: length of the vector is scaled to 1
      !
      implicit none
      !
      character*81 objectset(4,*)
      !
      integer nobject,nk,nodedesi(*),i,ndesi,m,iscaleflag,kflag,node
      !
      real*8 dgdxglob(2,nk,nobject),dd,len
      !
      !
      kflag=0
      if(iscaleflag.eq.0) then
         if(objectset(1,nobject)(1:11).eq.'PROJECTGRAD') then
            kflag=1
         endif
         !
         len=0.d0
         dd=0.d0
         if(kflag.eq.1) then
            do i=1,ndesi
               node=nodedesi(i)
               len=len+dgdxglob(2,node,nobject)**2
               dd=max(dd,abs(dgdxglob(1,node,nobject)))
            enddo
            if(dd.ne.0.d0) then
               do i=1,ndesi
                  node=nodedesi(i)
                  dgdxglob(1,node,nobject)=dgdxglob(1,node,nobject)/dd
               enddo
            endif
         else
            do i=1,ndesi
               node=nodedesi(i)
               len=len+dgdxglob(2,node,1)**2
            enddo
         endif
         len=dsqrt(len)
         !          write(5,*)
         !          write(5,*) 'LENGTH OF DESCENT GRADIENT VECTOR:'
         !          write(5,*)
         !          write(5,'(7x,e14.7)') len
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
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
      elseif(iscaleflag.eq.1) then
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
      endif
      !
      return        
      end
      



