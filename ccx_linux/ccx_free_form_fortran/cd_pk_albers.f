!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2019 Guido Dhondt
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
!     author: Yannick Muller
!
      subroutine cd_pk_albers(rad,d,xl,reynolds,p2,p1,beta,kappa,cd,u,&
           T1,R)
      !
      implicit none
      !
      real*8 rad,d,xl,reynolds,p2,p1,beta,kappa,&
           cd,R,u,T1
      !
      rad=rad
      d=d
      xl=xl
      reynolds=reynolds
      p2=p2
      p1=p1
      beta=beta
      kappa=kappa
      R=R
      u=u
      T1=T1

      
      cd=1.d0

      write(*,*) '*WARNING while using subroutine cd_pk_albers.f'
      write(*,*) 'cd implicitely taken equal to 1'
      !
      return
      end
      
