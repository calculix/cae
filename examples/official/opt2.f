!
!     reads coordinates and displacement results from opt2.inp
!     stores new position of nodes to file
!
      implicit none
!
      character*80 text
      integer i,j
      real*8 co(3,261),disp(3,261)
!
!     read coordinates and displacements from opt2.frd
!
      open(1,file='opt2.frd',status='old')
!
      do
         read(1,'(a80)',end=1) text
!
         if(text(5:6).eq.'2C') then
            do i=1,261
               read(1,'(13x,3e12.5)') (co(j,i),j=1,3)
            enddo
         endif
!
         if(text(6:8).eq.'ALL') then
            do i=1,261
               read(1,'(13x,3e12.5)') (disp(j,i),j=1,3)
            enddo
            exit
         endif
      enddo
!
 1    close(1)
!
!     write new position of mesh nodes to file
!
      open(2,file='opt3.inc',status='unknown')
      write(2,300)
 300  format('*NODE')
      do i=1,261
         write(2,301) i,co(1,i)+disp(1,i),
     &                  co(2,i)+disp(2,i),
     &                  co(3,i)+disp(3,i)
      enddo
 301  format(i5,3(',',e20.13))
      close(2)
!
      return
      end
