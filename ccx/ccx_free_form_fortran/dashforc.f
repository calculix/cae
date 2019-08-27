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
      subroutine dashforc(xl,konl,vl,imat,elcon,nelcon,&
        elas,fnl,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,&
        plicon,nplicon,npmat_,vel,time,nmethod,mi)
      !
      !     calculates the force of the dashpot
      !
      implicit none
      !
      character*8 lakonl
      !
      integer konl(20),i,j,imat,ncmat_,ntmat_,nope,nmethod,&
        kode,nelcon(2,*),nplicon(0:ntmat_,*),npmat_,id,niso,mi(*)
      !
      real*8 xl(3,20),elas(21),t0l,t1l,vl(0:mi(2),20),plconloc(802),&
        pl(0:3,9),xn(3),al,dd,vel(1:3,20),time,&
        elcon(0:ncmat_,ntmat_,*),elconloc(21),xk,fk,fnl(3,9),&
        plicon(0:2*npmat_,ntmat_,*),xiso(200),yiso(200)
      !
      !     actual positions of the nodes belonging to the dashpot
      !
      do i=1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+vl(j,i)
         enddo
      enddo
      !
      dd=dsqrt((pl(1,2)-pl(1,1))**2&
           +(pl(2,2)-pl(2,1))**2&
           +(pl(3,2)-pl(3,1))**2)
      do i=1,3
         xn(i)=(pl(i,2)-pl(i,1))/dd
      enddo
      !
      al=(vel(1,2)-vel(1,1))*xn(1)+&
         (vel(2,2)-vel(2,1))*xn(2)+&
         (vel(3,2)-vel(3,1))*xn(3)
      !
      !     interpolating the material data
      !
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,&
           elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
      !
      !     calculating the dashpot force and the dashpot constant
      !
      if(kode.eq.2)then
         xk=elconloc(1)
         fk=xk*al
      else
         if(nmethod.ne.5) then
            write(*,*) '*ERROR in dashdamp: the damping coefficient'
            write(*,*) '       may depend on temperature and frequency'
            write(*,*) '       only; the latter is only allowed for'
            write(*,*) '       steady state dynamics calculations'
            call exit(201)
         endif
         niso=int(plconloc(801))
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,time,niso,id)
         if(id.eq.0) then
            xk=yiso(1)
         elseif(id.eq.niso) then
            xk=yiso(niso)
         else
            xk=yiso(id)+(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))&
               *(time-xiso(id))
         endif
         fk=xk*al
      endif
      !
      do i=1,3
         fnl(i,1)=-fk*xn(i)
         fnl(i,2)=fk*xn(i)
      enddo
      !
      return
      end
      
      
