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
      subroutine writematrix(au,ad,irow,jq,neq,number)
      ! >
      ! >     \brief writes an matrix to file (for debugging purposes)
      ! >
      implicit none
      !
      character*12 name
      character*14 name2
      !
      integer irow(*),jq(*),neq,i,j,k,number,ii,row,&
           column,idiff
      !
      real*8 au(*),ad(*),help,aij,aji,diff,&
           maxdiff,maxhelp
      !
      name='matrix_'//char(number+96)//'.out'
      name2='matrix_'//char(number+96)//'_t.out'
      open(10,file=name,status='unknown')
      write(10,*) 'matrix number ',number
      !
      do i=1,neq
         if(ad(i).gt.1.e-20 .or. ad(i).lt. -1.e-20)then
            write(10,*) 'row ',i,' value ',ad(i)
         endif
      enddo
      !
      diff=0.0
      idiff=0
      maxdiff=0.0
      !
      do i=1,neq
         k=jq(i+1)-jq(i)
         if(k.gt.0)write(10,*) 'comlumn ', i
         help=0.0
         do j=jq(i),jq(i+1)-1
            aij=0.0
            aji=0.0
            row=i
            column=irow(j)
            aij=au(j)
               !             do ii=jq(column),jq(column+1)-1
               !                if(irow(ii)==row)then
               !                   aji=au(ii)
               !                   exit
               !                endif
               !             enddo
               !             diff=abs(aij-aji)
               !             if(diff.lt.1.0d-17) then
               !                write(10,100) i,irow(j),au(j)
               !             else
               !                if(diff.gt.maxdiff) maxdiff=diff
               !                idiff=idiff+1
               write(10,100) i,irow(j),au(j)
         !             endif
         !             help=help+au(j)
         enddo
      enddo
 !       write(10,*)'maxdiffsymm',maxdiff,'idiff',idiff
 !       maxhelp=0.0
 !       do i=1,neq
 !          help=0.0
 !          if(ad(i).gt.0.0)then
 !             help=ad(i)+help
 !             write(10,*) help
 !             do ii=1,neq
 !                do j=jq(ii),jq(ii+1)-1
 !                   if(irow(j).eq.i) then
 !                      help=help+au(j)
 !                      write(10,*) ii,au(j),help
 !                   endif
 !                enddo
 !             enddo

      !             if(abs(help).gt.maxhelp) maxhelp=abs(help)
      !             write(10,*)'column',i,'diff',help
      !          endif
      !       enddo
      write(10,*)'maxdiff_impuls',maxhelp
      !
      close(10)
 100  format('column ',i10,1x,'row ', i10,1x,'value ',e15.8)      
 101  format('column ',i10,1x,'row ', i10,1x,'value ',e15.8,&
           1x,'diff',e15.8)
      return
      end
      
