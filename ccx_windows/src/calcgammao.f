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
      subroutine calcgammao(nface,ielfa,vel,gradoel,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef,flux)
!
!     determine gamma for the turbulent dissipation rate:
!        upwind difference: gamma=0
!        central difference: gamma=1
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:7),gradoel(3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*)
!
      do i=1,nface
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         vcd=vel(iel2,7)-vel(iel1,7)
         if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,7))) vcd=0.d0
!
         if(flux(indexf).ge.0.d0) then
            vud=2.d0*xlet(indexf)*
     &           (gradoel(1,iel1)*xxj(1,indexf)+
     &            gradoel(2,iel1)*xxj(2,indexf)+
     &            gradoel(3,iel1)*xxj(3,indexf))
         else
            vud=2.d0*xlet(indexf)*
     &           (gradoel(1,iel2)*xxj(1,indexf)+
     &            gradoel(2,iel2)*xxj(2,indexf)+
     &            gradoel(3,iel2)*xxj(3,indexf))
         endif
!
         if(dabs(vud).lt.1.d-20) then
            gamma(i)=0.d0
            cycle
         endif
!            
         phic=1.d0-vcd/vud
!
         if(phic.ge.1.d0) then
            gamma(i)=0.d0
         elseif(phic.le.0.d0) then
            gamma(i)=0.d0
         elseif(betam.le.phic) then
            gamma(i)=1.d0
         else
            gamma(i)=phic/betam
         endif
c         write(*,*) 'calcgammat',i,gamma(i)
      enddo
!            
      return
      end
