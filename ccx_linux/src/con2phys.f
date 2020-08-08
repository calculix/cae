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
      subroutine con2phys(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     calculates the physical variable from the conservative
!     variables:
!       for gases: temperature,velocity and pressure 
!       for thermal liquids: temperature, velocity and density
!       for athermal liquids: velocity and density
!
      implicit none
!
      integer convergence,compressible,
     &  nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,mi(*),
     &  nshcon(*),ielmat(mi(3),*),nk,ithermal(*),i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun,nnn
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),vcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  vcontu(2,*),vtu(2,*),co(3,*),factor
!     
      if(ithermal(1).gt.1) then
!     
         do i=1,nk
            if(inomat(i).eq.0) cycle
            imat=inomat(i)
            temp=vold(0,i)
!     
            if(compressible.eq.1) then
!     
!              gas
!
!              calculation of the static temperature
!     
               rho=vcon(4,i)
               r=shcon(3,1,imat)
               c1=(vcon(0,i)-(vcon(1,i)**2+vcon(2,i)**2+
     &              vcon(3,i)**2)/(2.d0*rho))/rho
!
!              the energy density per volume unit cannot be
!              negative
!
               if(c1.lt.1.d-10) c1=1.d-10
!     
               temp0=temp
               j=0
               do
                  call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &                 nshcon,cp,physcon)
                  temp=c1/(cp-r)+physcon(1)
c                  write(*,*) 'con2phys ',j,temp
c     start shallow
c     temp=max(c1/(cp),1.d-2)+physcon(1)
c     end shallow
                  j=j+1
                  if((dabs(temp-temp0).lt.1.d-4*dabs(temp)).or.
     &                 (dabs(temp-temp0).lt.1.d-10)) then
                     vold(0,i)=temp
                     exit
                  endif
                  if(j.gt.100) then
                     write(*,*) 
     &                    '*ERROR in con2phys: too many iterations'
                     write(*,*) '       for node',i
                     write(*,*) '       actual temperature ',temp,' K'
                     stop
                  endif
                  temp0=temp
               enddo
!     
!              determining the pressure (gas equation)
!     
               vold(4,i)=rho*r*(temp-physcon(1))
c     start shallow
c     vold(4,i)=5.d0*(rho*rho-(0.005*co(1,i))**2)
c     end shallow
!     
!              calculating the velocity from the conservative variables
!     
               do k=1,3
                  if(nactdoh(k,i).ne.0) then
                     vold(k,i)=vcon(k,i)/rho
                  endif
               enddo
               cycle
!     
            else
!     
!              thermal liquid
!     
               c1=vcon(0,i)
               c2=(vcon(1,i)**2+vcon(2,i)**2+vcon(3,i)**2)/2.d0
               temp0=temp
               j=0
!     
!              iterating to find the temperature
!     
               do
                  call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &                 nshcon,cp,physcon)
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_,ithermal)
                  temp=(c1-c2/rho)/(rho*cp)+physcon(1)
                  j=j+1
                  if((dabs(temp-temp0).lt.1.d-4*dabs(temp)).or.
     &               (dabs(temp-temp0).lt.1.d-10)) then
                     vold(0,i)=temp
                     exit
                  endif
                  if(j.gt.100) then
                     write(*,*) 
     &                 '*ERROR in con2phys: too many iterations'
                     write(*,*) '       for node',i
                     write(*,*) '       actual temperature ',temp,' K'
                     stop
                  endif
                  temp0=temp
               enddo
!
!              storing the density
!
               vcon(4,i)=rho
!     
!              calculating the velocity
!     
               do k=1,3
                  if(nactdoh(k,i).ne.0) then
                     vold(k,i)=vcon(k,i)/rho
                  endif
               enddo
            endif
         enddo
      else
!     
!     athermal liquid calculation
!     
         do i=1,nk
            if(inomat(i).eq.0) cycle
            imat=inomat(i)
            temp=vold(0,i)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
!     
!           storing the density
!
            vcon(4,i)=rho
!
!           calculating the velocity
!     
            do k=1,3
               vold(k,i)=vcon(k,i)/rho
            enddo
         enddo
      endif
!     
      return
      end
      
      
