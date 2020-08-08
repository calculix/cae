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
      subroutine writelm(iter,xlambd,nactive,nnlconst,objectset,
     &   nobject,ipoacti,iconstacti,inameacti)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer nactive,nobject,nnlconst,iter,ipos,i,ipoacti(*),
     &   iconstacti(*),inameacti(*)   
!
      real*8 xlambd(*)
!
      write(5,*)
      if(iter.eq.1) then      
         write(5,*) 'LAGRANGE MULTIPLIERS IN THE 1st ITERATION:'
      elseif(iter.eq.2) then
         write(5,*) 'LAGRANGE MULTIPLIERS IN THE 2nd ITERATION:'
      elseif(iter.eq.3) then
         write(5,*) 'LAGRANGE MULTIPLIERS IN THE 3rd ITERATION:'
      elseif((iter.gt.3).and.(iter.lt.10)) then
         write(5,'(a28,i1,a4,a13)') 'LAGRANGE MULTIPLIERS IN THE ',
     &            iter,'th ITERATION:' 
      else
         write(5,'(a28,i2,a4,a13)') 'LAGRANGE MULTIPLIERS IN THE ',
     &            iter,'th ITERATION:' 
      endif
      write(5,*)
      do i=1,nactive
         ipos=ipoacti(i)
         if(i.le.nnlconst) then
            if(iconstacti(i).eq.-1) then
               if(xlambd(i).gt.0.d0) then            
                  write(5,'(7x,a12,a4,e14.7,a12)') objectset(1,ipos),
     &                    'LE  ',xlambd(i),'active    '
               else
                  write(5,'(7x,a12,a4,e14.7,a12)') objectset(1,ipos),
     &                    'LE  ',xlambd(i),'not active'
               endif
            else
               if(xlambd(i).gt.0.d0) then
                  write(5,'(7x,a12,a4,e14.7,a12)') objectset(1,ipos),
     &                    'GE  ',xlambd(i),'not active' 
               else
                  write(5,'(7x,a12,a4,e14.7,a12)') objectset(1,ipos),
     &                    'GE  ',xlambd(i),'active    '
               endif
            endif   
         else
            if(iconstacti(i).eq.-1) then      
               if(xlambd(i).gt.0.d0) then   
                  write(5,'(7x,a12,a4,e14.7,a12)')
     &                    objectset(1,inameacti(i)),
     &                    'LE  ',xlambd(i), 'active    '
               else
                  write(5,'(7x,a12,a4,e14.7,a12)')
     &                    objectset(1,inameacti(i)),
     &                    'LE  ',xlambd(i), 'not active'
               endif
            else
               if(xlambd(i).gt.0.d0) then   
                  write(5,'(7x,a12,a4,e14.7,a12)')
     &                    objectset(1,inameacti(i)),
     &                    'GE  ',xlambd(i), 'not active'
               else
                  write(5,'(7x,a12,a4,e14.7,a12)')
     &                    objectset(1,inameacti(i)),
     &                    'GE  ',xlambd(i), 'active    '
               endif
            endif
         endif
      enddo  
!
      return        
      end




