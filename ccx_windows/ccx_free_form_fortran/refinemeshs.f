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
      subroutine refinemeshs(inpc,textpart,filab,istep,istat,n,iline,&
        ipol,inl,ipoinp,inp,ipoinpc,ier)
      !
      !     reading *REFINE MESH in the input deck
      !
      implicit none
      !
      character*1 inpc(*)
      character*87 filab(*)
      character*132 textpart(16)
      !
      integer istep,istat,n,key,ii,ier,iline,ipol,inl,ipoinp(2,*),&
        inp(3,*),ipoinpc(0:*)
      !
      if(istep.lt.1) then
         write(*,*)&
      '*ERROR reading *REFINE MESH: *REFINE MESH'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         ier=1
         return
      endif
      !
      !     reset the refinement controlling field
      !
      filab(48)(1:6)='RM    '
      !
      do ii=2,n
        if(textpart(ii)(1:6).eq.'LIMIT=') then
           filab(48)(7:26)=textpart(ii)(7:26)
        elseif(textpart(ii)(1:4).eq.'USER') then
           filab(48)(1:6)='RMUSER'
        else
            write(*,*)&
                   '*WARNING reading *REFINE MESH:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',&
                       textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,&
      "*REFINE MESH %")
        endif
      enddo
      !
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,&
              ipoinp,inp,ipoinpc)
         if((key.eq.1).or.(istat.lt.0)) return
         !
         if(n.gt.1) then
            write(*,*) '*WARNING reading *REFINE MESH'
            write(*,*) '         only one refinement controlling'
            write(*,*) '         field allowed'
            call inputwarning(inpc,ipoinpc,iline,&
      "*REFINE MESH %")
            return
         endif
         !
         do ii=1,n
            if((textpart(ii)(1:4).eq.'ERR ').or.&
               (textpart(ii)(1:4).eq.'U   ').or.&
               (textpart(ii)(1:4).eq.'NT  ').or.&
               (textpart(ii)(1:4).eq.'S   ').or.&
               (textpart(ii)(1:4).eq.'E   ').or.&
               (textpart(ii)(1:4).eq.'ME  ').or.&
               (textpart(ii)(1:4).eq.'PEEQ').or.&
               (textpart(ii)(1:4).eq.'ENER').or.&
               (textpart(ii)(1:4).eq.'HFL')) then
               filab(48)(3:6)=textpart(ii)(1:4)
            else
               write(*,*)&
      '*WARNING reading *REFINE MESH: label not applicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline,&
      "*REFINE MESH %")
            endif
         enddo
      enddo
      !
      return
      end






