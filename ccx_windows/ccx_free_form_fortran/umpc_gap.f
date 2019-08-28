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
      subroutine umpc_gap(x,u,f,a,jdof,n,force,iit,idiscon)
      !
      !     updates the coefficients in a gap mpc (name GAP)
      !
      !     a gap MPC is triggered by a *GAP definition applied to
      !     a GAPUNI element. The gap direction is stored in
      !     x(1..3,7), the clearance in x(1,8), which is also the
      !     constant term
      !
      !     INPUT:
      !
      !     x(3,n)             Carthesian coordinates of the nodes in the
      !                        user mpc.
      !     u(3,n)             Actual displacements of the nodes in the
      !                        user mpc.
      !     jdof               Actual degrees of freedom of the mpc terms
      !     n                  number of terms in the user mpc
      !     force              Actual value of the mpc force
      !     iit                iteration number
      !
      !     OUTPUT:
      !
      !     f                  Actual value of the mpc. If the mpc is
      !                        exactly satisfied, this value is zero
      !     a(n)               coefficients of the linearized mpc
      !     jdof               Corrected degrees of freedom of the mpc terms
      !     idiscon            0: no discontinuity
      !                        1: discontinuity
      !                        If a discontinuity arises the previous
      !                        results are not extrapolated at the start of
      !                        a new increment
      !
      implicit none
      !
      integer jdof(*),n,iit,ifix,idiscon
      !
      real*8 x(3,*),u(3,*),f,a(*),dist(3),xn(3),force
      !
      !       write(*,*) (jdof(i),i=1,7)
      if(jdof(7).eq.1) then
         ifix=1
      else
         ifix=0
         jdof(7)=2
      endif
      !
      dist(1)=u(1,4)-u(1,1)
      dist(2)=u(2,4)-u(2,1)
      dist(3)=u(3,4)-u(3,1)
      !
      !     gap direction
      !
      xn(1)=x(1,7)
      xn(2)=x(2,7)
      xn(3)=x(3,7)
      !
      f=dist(1)*xn(1)+dist(2)*xn(2)+dist(3)*xn(3)+x(1,8)
      !
      a(7)=-1.d0
      !
      !     only one change per increment is allowed
      !        (change= from free to linked or vice versa)
      !     ifix=0: free
      !     ifix=1: linked
      !
      if(ifix.eq.0) then
         !
         !        previous state: free
         !
         if(f.gt.0) then
            !
            !           new state: free
            !
            f=0.d0
         elseif(iit.le.1) then
            !
            !           new state: linked
            !
            write(*,*) 'switch to linked'
            write(*,*)
            jdof(7)=1
            idiscon=1
         else
            !
            !           new state: free
            !
            f=0.d0
         endif
      else
         !
         !        previous state: linked
         !
         if(force.ge.0.d0) then
         !
         !           new state: linked
         !
         elseif(iit.le.1) then
            !
            !           new state: free
            !
            write(*,*) 'switch to free'
            write(*,*)
            jdof(7)=2
            f=0.d0
            idiscon=1
         else
         !
         !           new state: linked
         !
         endif
      endif
      !
      if(dabs(xn(jdof(1))).gt.1.d-10) then
         a(1)=-xn(jdof(1))
         if(jdof(1).eq.1) then
            jdof(2)=2
            jdof(3)=3
         elseif(jdof(1).eq.2) then
            jdof(2)=3
            jdof(3)=1
         else
            jdof(2)=1
            jdof(3)=2
         endif
         a(2)=-xn(jdof(2))
         a(3)=-xn(jdof(3))
      else
         if(jdof(1).eq.3) then
            jdof(1)=1
         else
            jdof(1)=jdof(1)+1
         endif
         if(dabs(xn(jdof(1))).gt.1.d-10) then
            a(1)=-xn(jdof(1))
            if(jdof(1).eq.1) then
               jdof(2)=2
               jdof(3)=3
            elseif(jdof(1).eq.2) then
               jdof(2)=3
               jdof(3)=1
            else
               jdof(2)=1
               jdof(3)=2
            endif
            a(2)=-xn(jdof(2))
            a(3)=-xn(jdof(3))
         else
            if(jdof(1).eq.3) then
               jdof(1)=1
            else
               jdof(1)=jdof(1)+1
            endif
            if(dabs(xn(jdof(1))).gt.1.d-10) then
               a(1)=-xn(jdof(1))
               if(jdof(1).eq.1) then
                  jdof(2)=2
                  jdof(3)=3
               elseif(jdof(1).eq.2) then
                  jdof(2)=3
                  jdof(3)=1
               else
                  jdof(2)=1
                  jdof(3)=2
               endif
               a(2)=-xn(jdof(2))
               a(3)=-xn(jdof(3))
            endif
         endif
      endif
      !
      a(4)=xn(1)
      a(5)=xn(2)
      a(6)=xn(3)
      jdof(4)=1
      jdof(5)=2
      jdof(6)=3
      !
      !       write(*,*) 'jdof,a'
      !       do i=1,7
      !          write(*,*) jdof(i),a(i)
      !       enddo
      !       write(*,*) 'f ',f
      !
      return
      end













