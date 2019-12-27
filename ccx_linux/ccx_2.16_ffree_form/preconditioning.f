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
!     diagonal preconditioning of a matrix
!
      subroutine preconditioning(ad,au,b,neq,irow,jq,adaux)
      !
      implicit none
      !
      integer neq,irow(*),jq(*),i,ic,ir
      !
      real*8 ad(*),au(*),b(*),adaux(*),adc
      !
      intent(in) neq,irow,jq
      !
      intent(inout) ad,au,b,adaux
      !
      !     inverse of the square root of the diagonal
      !     the sign takes care that the diagonal term becomes 1
      !     (and not -1)
      !
      !     taking zero's on the diagonal into account (adaux(i)=1 in such case)
      !
      do i=1,neq
         !      write(*,*) 'preconditioning ',i,ad(i)
         if(dabs(ad(i)).lt.1.d-30) then
            adaux(i)=dsign(1.d0,ad(i))
         else
            adaux(i)=dsign(1.d0/dsqrt(dabs(ad(i))),ad(i))
         endif
      enddo
      !       do i=1,neq
      !          adaux(i)=1.d0/dsqrt(dabs(ad(i)))
      !       enddo
      !
      !     scaling the matrix and the right hand side
      !
      do ic=1,neq
         adc=dabs(adaux(ic))
         !
         !        scaling the diagonal
         !
         ad(ic)=ad(ic)*adc*adaux(ic)
         !
         !        scaling the off-diagonal terms
         !
         do i=jq(ic),jq(ic+1)-1
            ir=irow(i)
            !             write(*,*) 'au before',i,au(i),adc,ir,adaux(ir)
            au(i)=au(i)*adc*adaux(ir)
         !             write(*,*) 'au after',i,au(i)
         enddo
         !
         !        scaling the right hand side
         !
         b(ic)=b(ic)*adaux(ic)
      enddo
      !
      !     taking the absolute value
      !
      do i=1,neq
         adaux(i)=dabs(adaux(i))
      enddo
      !
      return
      end
