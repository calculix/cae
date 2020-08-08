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
      subroutine hrt_mod_smart(nface,ielfa,vel,gradtel,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef,flux,vfa)
!
!     use the modified smart scheme to determine the facial
!     temperature
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:7),gradtel(3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*),vfa(0:7,*)
!
      do i=1,nface
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!        unless outlet
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
!           outflow && (neighbor || outlet)
!
!           outlet
!
            if(iel2.le.0) then
               vfa(0,i)=vel(iel1,0)
               cycle
            endif
!
!           neighbor
!
            vcd=vel(iel1,0)-vel(iel2,0)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,0))) vcd=0.d0
!
            vud=2.d0*xlet(indexf)*
     &           (gradtel(1,iel1)*xxj(1,indexf)+
     &            gradtel(2,iel1)*xxj(2,indexf)+
     &            gradtel(3,iel1)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(0,i)=vel(iel1,0)
               cycle
            endif
!     
            phic=1.d0+vcd/vud
c            write(*,*) 'calcvfa1 ',i,phic
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(0,i)=vel(iel1,0)
            elseif(phic.le.1.d0/6.d0) then
               vfa(0,i)=3.d0*vel(iel1,0)-2.d0*vel(iel2,0)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(0,i)=3.d0*vel(iel1,0)/4.d0+vel(iel2,0)/4.d0+vud/8.d0
            else
               vfa(0,i)=vel(iel1,0)/3.d0+2.d0*vel(iel2,0)/3.d0
            endif
         elseif(iel2.gt.0) then
!
            vcd=vel(iel2,0)-vel(iel1,0)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel2,0))) vcd=0.d0
!
            vud=-2.d0*xlet(indexf)*
     &           (gradtel(1,iel2)*xxj(1,indexf)+
     &            gradtel(2,iel2)*xxj(2,indexf)+
     &            gradtel(3,iel2)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(0,i)=vel(iel2,0)
               cycle
            endif
!     
            phic=1.d0+vcd/vud
c            write(*,*) 'calcvfa2 ',i,phic
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(0,i)=vel(iel2,0)
            elseif(phic.le.1.d0/6.d0) then
               vfa(0,i)=3.d0*vel(iel2,0)-2.d0*vel(iel1,0)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(0,i)=3.d0*vel(iel2,0)/4.d0+vel(iel1,0)/4.d0+vud/8.d0
            else
               vfa(0,i)=vel(iel2,0)/3.d0+2.d0*vel(iel1,0)/3.d0
            endif
         endif
      enddo
!            
      return
      end
