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
      subroutine equcons(inpc,textpart,istep,istat,n,iline,ipol,inl,&
              ipoinp,inp,ipoinpc,equcon,equbou)
      !
      !     reading the input deck: *EQUALITY CONSTRAINTS
      !
      !     criteria: MASS
      !               STRESS
      !               DISPLACEMENT
      !
      implicit none
      !
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 equcon(2)
      !
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),&
        inp(3,*),ipoinpc(0:*)
      !
      real*8 equbou
      !
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *EQUALITY CONSTRAINTS:&
      *EQUALITY CONSTRAINTS can only be used within a SENSITIVITY STEP'
         call exit(201)
      endif
      !
      do i=2,n
         if(textpart(i)(1:9).eq.'CRITERIA=') then
            read(textpart(i)(10:85),'(a80)',iostat=istat)&
                 equcon(1)(1:80)
         elseif(textpart(i)(1:7).eq.'ENTITY=') then
            read(textpart(i)(8:85),'(a80)',iostat=istat)&
                 equcon(2)(1:80)
         endif
      enddo
      !
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,&
              ipoinp,inp,ipoinpc)
      !
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) equbou
      !
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,&
              ipoinp,inp,ipoinpc)
      !
      return
      end

