!
!     applying the steepest descent concept to the sensitivities
!     calculated in opt1.inp 
!
!     results in geometrical changes needed to decrease the value
!     of the objective function. These are written as boundary 
!     conditions for a subsequent linear elastic calculation
!
      implicit none
!
      character*80 text
      integer i,j
      real*8 xn(3,261),sen(261)
!
!     reading the sensitivities
!
      open(1,file='opt1.frd',status='old')
!
      do
         read(1,'(a80)',end=1) text
!
         if(text(6:10).eq.'NORMZ') then
            read(1,*)
            do i=1,261
               read(1,'(13x,3e12.5)') (xn(j,i),j=1,3)
            enddo
         endif
!
         if(text(6:12).eq.'PRJGRAD') then
            read(1,*)
            read(1,*)
            do i=1,261
               read(1,'(25x,e12.5)') sen(i)
            enddo
            exit
         endif
      enddo
!
!     the sensitivity is the change of the objective function while
!     thickening the body in local normal direction; in order to
!     decrease the objective function one has to thicken the body if
!     the sensitivity is negative and make the body more slender if
!     the sensitivity is positive.
!
      do i=1,261
         sen(i)=-sen(i)*1.d-1
      enddo
!
 1    close(1)
!
!     write thickness changes as boundary conditions in file
!
      open(2,file='opt1.bou',status='unknown')
      write(2,300)
 300  format('*BOUNDARY')
      do i=1,261
!
!        no change if the norm of the normal is not 1 or if the
!        sensitivity is zero
!
         if(xn(1,i)**2+xn(2,i)**2+xn(3,i)**2.lt.0.5d0) cycle
         if(dabs(sen(i)).le.0.d0) cycle
         write(2,301) i,xn(1,i)*sen(i)
         write(2,302) i,xn(2,i)*sen(i)
         write(2,303) i,xn(3,i)*sen(i)
      enddo
 301  format(i5,',1,1,',e20.13)
 302  format(i5,',2,2,',e20.13)
 303  format(i5,',3,3,',e20.13)
      close(2)
!
      return
      end

      
