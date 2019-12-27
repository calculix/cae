!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine writerefinemesh(kontet,netet_,cotet,nktet,jobnamec,&
                 ipkon,kon,lakon,iquad,iedtet,iedgmid,itreated,ne)
      !
      implicit none
      !
      character*8 lakon(*)
      character*132 fninp,jobnamec
      !
      integer kontet(4,*),netet_,i,j,nktet,ipkon(*),kon(*),&
        indexe,iquad,iedtet(6,*),iedgmid(*),itreated(*),ne
      !
      real*8 cotet(3,*)
      !
      intent(in) kontet,netet_,nktet,jobnamec,ne,&
                 ipkon,kon,lakon,iquad,iedtet,iedgmid
      !
      intent(inout) cotet
      !
      !     stores the refined mesh in input format
      !
      do i=1,132
         if(ichar(jobnamec(i:i)).eq.0) exit
      !          write(*,*) 'writerefinemesh ',i,ichar(jobnamec(i:i))
      enddo
      if(i.gt.129) then
         write(*,*) '*ERROR in writerefinemesh'
         write(*,*) '       jobname has more than 129 characters'
         call exit(201)
      endif
      fninp(1:i+3)=jobnamec(1:i-1)//'.fin'
      !
      !     storing the mesh in input format
      !
      open(2,file=fninp(1:i+3),status='unknown')
      !
      !     storing the nodes
      !
      write(2,102)
 102  format('*NODE')
      do i=1,nktet
         !
         !        setting too small numbers to zero (else the exponent in the
         !        output contains 3 digits and the letter "D" is omitted)
         !
         do j=1,3
            if(dabs(cotet(j,i)).lt.1.d-99) cotet(j,i)=0.d0
         enddo
         write(2,100) i,(cotet(j,i),j=1,3)
      enddo
      !
      !     storing the tetrahedral elements
      !
      if(iquad.eq.0) then
         write(2,103)
 103     format('*ELEMENT,TYPE=C3D4,ELSET=TET')
         do i=1,netet_
            if(kontet(1,i).ne.0) then
               write(2,101) i,(kontet(j,i),j=1,4)
            endif
         enddo
      else
         write(2,104)
 104     format('*ELEMENT,TYPE=C3D10,ELSET=TET')
         do i=1,netet_
            if(kontet(1,i).ne.0) then
               write(2,101) i,(kontet(j,i),j=1,4),&
                 (iedgmid(iedtet(j,i)),j=1,6)
            endif
         enddo
      endif
      !
      !     storing all other elements
      !
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if((lakon(i)(1:4).eq.'C3D4').or.&
            (lakon(i)(1:5).eq.'C3D10')) cycle
         if(lakon(i)(1:4).eq.'C3D6') then
            indexe=ipkon(i)
            write(2,105) lakon(i)
 105        format('*ELEMENT,TYPE=',a8)
            write(2,101) i,(kon(indexe+j),j=1,6)
         elseif(lakon(i)(1:4).eq.'C3D8') then
            indexe=ipkon(i)
            write(2,106) lakon(i)
 106        format('*ELEMENT,TYPE=',a8)
            write(2,101) i,(kon(indexe+j),j=1,8)
         elseif(lakon(i)(1:5).eq.'C3D15') then
            indexe=ipkon(i)
            write(2,107) lakon(i)
 107        format('*ELEMENT,TYPE=',a8)
            write(2,101) i,(kon(indexe+j),j=1,15)
         elseif(lakon(i)(1:5).eq.'C3D20') then
            indexe=ipkon(i)
            write(2,108) lakon(i)
 108        format('*ELEMENT,TYPE=',a8)
            write(2,101) i,(kon(indexe+j),j=1,15)
            write(2,101) (kon(indexe+j),j=16,20)
         endif
      enddo
      !
      !     write information whether projection took place on edge
      !     or on surface
      !
      !       write(2,109)
      !  109  format('*TEMPERATURE')
      !       do i=1,nktet
      !          write(2,100) i,1.d0*itreated(i)
      !       enddo
      !
      close(2)
 !
 100  format(i10,',',e20.13,',',e20.13,',',e20.13)
 101  format(11(i10,','))
      !
      return
      end
