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
      subroutine materialdata_cfd1(nef,vel,shcon,nshcon,ielmatf,&
        ntmat_,mi,cvel,physcon,ithermal,umel,rhcon,nrhcon,nefa,nefb)
      !
      !     calculation of material properties at elements centers
      !     (incompressible fluids)
      !
      implicit none
      !
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal,&
        nshcon(2,*),nrhcon(*),nefa,nefb
      !
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),cvel(*),&
        physcon(*),umel(*),rhcon(0:1,ntmat_,*)
      !
      intent(in) nef,shcon,nshcon,ielmatf,ntmat_,mi,physcon,ithermal,&
        nefa,nefb
      !
      intent(inout) vel,cvel,umel
      !
      do i=nefa,nefb
         t1l=vel(i,0)
         imat=ielmatf(1,i)
         !
         !        density
         !
         call materialdata_rho(rhcon,nrhcon,imat,vel(i,5),t1l,ntmat_,&
                  ithermal)
         !
         !        heat capacity at constant volume
         !        (for liquids: =heat capacity at constant pressure)
         !
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cvel(i),&
             physcon)
         !
         !        dynamic viscosity
         !
         call materialdata_dvi(shcon,nshcon,imat,umel(i),t1l,ntmat_,&
                  ithermal)
      enddo
      !
      return
      end
