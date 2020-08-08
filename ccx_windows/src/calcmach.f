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
      subroutine calcmach(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     calculates 
!       vold (temperature,velocity and pressure)
!       vcon (volumetric energy density, volumetric momentum
!                density and density)
!       at the nodes    
!
!     prints if iout=1
!
      implicit none
!
      integer convergence,compressible,
     &  nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,mi(*),
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
!     calculate kappa (cp/cv) and store it in v(0,*)
!     calculate the Mach number and store it in v(1,*)
!     
      do i=1,nk
         imat=inomat(i)
         temp=vold(0,i)
         call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &        nshcon,cp,physcon)
         r=shcon(3,1,imat)
         vel2=vold(1,i)**2+vold(2,i)**2+vold(3,i)**2
         v(0,i)=cp/(cp-r)
         v(1,i)=dsqrt((vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
     &        /(v(0,i)*r*(temp-physcon(1))))
      enddo
!     
      return
      end
      
