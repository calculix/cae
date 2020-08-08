!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine updatecon(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     updating the conservative variables
!
      implicit none
!
      integer convergence,compressible
!
      integer nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,mi(*),
     &  nshcon(*),ielmat(mi(3),*),nk,ithermal(*),i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),vcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  vcontu(2,*),vtu(2,*),co(3,*),factor
!     
!     volumetric energy density
!     
      if(ithermal(1).gt.1) then
         do i=1,nk
            vcon(0,i)=vcon(0,i)+v(0,i)
         enddo
      endif
!     
!     volumetric momentum density
!     pressure (liquid) or density (gas)
!     
      do i=1,nk
         if(inomat(i).eq.0) cycle
!
         do j=1,3
            vcon(j,i)=vcon(j,i)+v(j,i)
         enddo
!
         if(compressible.eq.1) then
            vcon(4,i)=vcon(4,i)+v(4,i)
         else
            vold(4,i)=vold(4,i)+v(4,i)
         endif
      enddo
!     
!     volumetric turbulent density
!     
      if(turbulent.ne.0) then
         do i=1,nk
            if(inomat(i).eq.0) cycle
            vcontu(1,i)=max(1.d-10,vcontu(1,i)+vtu(1,i))
            vcontu(2,i)=max(1.d-10,vcontu(2,i)+vtu(2,i))
c            vcontu(1,i)=vcontu(1,i)+vtu(1,i)
c            vcontu(2,i)=vcontu(2,i)+vtu(2,i)
         enddo
      endif
!     
      return
      end
      
