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
      subroutine writerefinemesh(kontet,netet_,cotet,nktet,jobnamec,
     &     ipkon,kon,lakon,iquad,iedtet,iedgmid,ne,
     &     number,jfix,iparentel,nk)
!
      implicit none
!
      character*8 lakon(*)
      character*10 elestr
      character*132 fnrfn,jobnamec,el_header
!
      integer kontet(4,*),netet_,i,j,k,nktet,ipkon(*),kon(*),node,
     &     indexe,iquad,iedtet(6,*),iedgmid(*),ne,number(*),nk,
     &     jfix(*),iparentel(*)
!
      real*8 cotet(3,*)
!
!     give nodes of the unrefined mesh which were not fixed
!     a new node number in order to avoid collisions with the
!     refined mesh
!      
      do i=1,netet_
        if(kontet(1,i).ne.0) then
          do j=1,4
            node=kontet(j,i)
            if((jfix(node).ne.1).and.(node.le.nk)) then
              if(number(node).ne.0) then
                kontet(j,i)=number(node)
              else
                nktet=nktet+1
                number(node)=nktet
                kontet(j,i)=nktet
                do k=1,3
                  cotet(k,nktet)=cotet(k,node)
                enddo
              endif
            endif
          enddo
        endif
      enddo
!
!     stores the refined mesh in input format
!
      do i=1,132
         if(ichar(jobnamec(i:i)).eq.0) exit
      enddo
      if(i.gt.125) then
         write(*,*) '*ERROR in writerefinemesh'
         write(*,*) '       jobname has more than 124 characters'
         call exit(201)
      endif
      fnrfn(1:i+7)=jobnamec(1:i-1)//'.rfn.inp'
!
!     storing the mesh in input format
!
      open(2,file=fnrfn(1:i+7),status='unknown',position='append')
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
        do i=1,netet_
          if(kontet(1,i).ne.0) then
!     
!     keyword card
!     
            write(elestr,'(i10)') iparentel(i)
            do k=1,10
              if(elestr(k:k).ne.' ') exit
            enddo
            el_header='*ELEMENT,PARENT='//elestr(k:10)//
     &           ',TYPE=C3D4'               
            write(2,*) el_header(1:36-k+1)
!     
!     topology
!     
            write(2,101) i,(kontet(j,i),j=1,4)
          endif
        enddo
      else
        do i=1,netet_
          if(kontet(1,i).ne.0) then
!     
!     keyword card
!     
            write(elestr,'(i10)') iparentel(i)
            do k=1,11
              if(elestr(k:k).ne.' ') exit
            enddo
            el_header='*ELEMENT,PARENT='//elestr(k:10)//
     &           ',TYPE=C3D10'               
            write(2,*) el_header(1:37-k+1)
!     
!     topology
!     
            write(2,101) i,(kontet(j,i),j=1,4),
     &           (iedgmid(iedtet(j,i)),j=1,6)
          endif
        enddo
      endif
!     
!     storing all other elements
!     
c     do i=1,ne
c     if(ipkon(i).lt.0) cycle
c     if((lakon(i)(1:4).eq.'C3D4').or.
c     &      (lakon(i)(1:5).eq.'C3D10')) cycle
c     if(lakon(i)(1:4).eq.'C3D6') then
c     indexe=ipkon(i)
c     write(2,105) lakon(i)
c     105        format('*ELEMENT,TYPE=',a8)
c     write(2,101) i,(kon(indexe+j),j=1,6)
c     elseif(lakon(i)(1:4).eq.'C3D8') then
c     indexe=ipkon(i)
c     write(2,106) lakon(i)
c     106        format('*ELEMENT,TYPE=',a8)
c     write(2,101) i,(kon(indexe+j),j=1,8)
c     elseif(lakon(i)(1:5).eq.'C3D15') then
c     indexe=ipkon(i)
c     write(2,107) lakon(i)
c     107        format('*ELEMENT,TYPE=',a8)
c     write(2,101) i,(kon(indexe+j),j=1,15)
c     elseif(lakon(i)(1:5).eq.'C3D20') then
c     indexe=ipkon(i)
c     write(2,108) lakon(i)
c     108        format('*ELEMENT,TYPE=',a8)
c     write(2,101) i,(kon(indexe+j),j=1,15)
c     write(2,101) (kon(indexe+j),j=16,20)
c     elseif(lakon(i)(1:4).eq.'A3D4') then
c     indexe=ipkon(i)
c     write(2,103) 
c     write(2,101) i,(kon(indexe+j),j=1,4)
c     lakon(i)(1:1)='C'
c     elseif(lakon(i)(1:5).eq.'A3D10') then
c     indexe=ipkon(i)
c     write(2,104)
c     write(2,101) i,(kon(indexe+j),j=1,10)
c     lakon(i)(1:1)='C'
c     endif
c     enddo
!     
      close(2)
!     
 100  format(i10,',',e20.13,',',e20.13,',',e20.13)
 101  format(11(i10,','))
 103  format('*ELEMENT,TYPE=C3D4,ELSET=TET')
 104  format('*ELEMENT,TYPE=C3D10,ELSET=TET')
!     
      return
      end
