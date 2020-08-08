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
      subroutine hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
!     determine the critical depth
!     
      implicit none
!      
      real*8 xflow,rho,b,dg,sqrts0,hk,theta,tth,c1,xflow2,
     &  A,dBBdh,dAdh,BB,dhk
!
      hk=((xflow/(rho*b))**2/(dg*sqrts0))**(1.d0/3.d0)
!
      if(dabs(theta).lt.1.d-10) return
!
!     critical depth for trapezoid, non-rectangular cross section
!
      tth=dtan(theta)
      c1=rho*rho*dg*sqrts0
      xflow2=xflow*xflow
!
      do
         A=hk*(b+hk*tth)
         dBBdh=2.d0*tth
         dAdh=b+hk*dBBdh
         BB=dAdh
         dhk=(xflow2*BB-c1*A**3)/(xflow2*dBBdh-3.d0*c1*A*A*dAdh)
         if(dabs(dhk)/dhk.lt.1.d-3) exit
         hk=hk-dhk
      enddo
!
      return
      end
      

