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
      subroutine hrr_mod_smart(vfa,shcon,ielmat,ntmat_,
     &  mi,ielfa,ipnei,vel,nef,flux,gradpel,gradtel,xxj,
     &  xlet,nfacea,nfaceb)
!
!     calculation of the density at the face centers
!     (compressible fluids)
!
!     facial temperature and pressure is only used for external
!     faces
!
      implicit none
!
      integer i,j,imat,ntmat_,mi(*),ipnei(*),nef,iel1,iel2,
     &  ielmat(mi(3),*),ielfa(4,*),indexf,nfacea,nfaceb
!
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),vel(nef,0:7),flux(*),
     &  r,gradpel(3,*),gradtel(3,*),xxj(3,*),phic,vud,
     &  vcd,xlet(*)
!
!
!     
      do i=nfacea,nfaceb
!
!        take the material of the first adjacent element
!        (for the specific gas constant)        
!
         imat=ielmat(1,ielfa(1,i))
         r=shcon(3,1,imat)
!
!        calculate gamma
!
         iel2=ielfa(2,i)
!
!        no neighbor
!
         if(iel2.le.0) then
            vfa(5,i)=vfa(4,i)/(r*vfa(0,i))
            cycle
         endif
!         
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
!           outflow
!
            vcd=vel(iel1,4)-vel(iel2,4)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,4))) vcd=0.d0
!
            vud=2.d0*xlet(indexf)*
     &           (gradpel(1,iel1)*xxj(1,indexf)+
     &            gradpel(2,iel1)*xxj(2,indexf)+
     &            gradpel(3,iel1)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(5,i)=vel(iel1,4)/(r*vfa(0,i))
               cycle
            endif
!     
            phic=1.d0+vcd/vud
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(5,i)=vel(iel1,4)
            elseif(phic.le.1.d0/6.d0) then
               vfa(5,i)=3.d0*vel(iel1,4)-2.d0*vel(iel2,4)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(5,i)=3.d0*vel(iel1,4)/4.d0+vel(iel2,4)/4.d0+vud/8.d0
            else
               vfa(5,i)=vel(iel1,4)/3.d0+2.d0*vel(iel2,4)/3.d0
            endif
            vfa(5,i)=vfa(5,i)/(r*vfa(0,i))
         elseif(iel2.gt.0) then
!
            vcd=vel(iel2,4)-vel(iel1,4)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel2,4))) vcd=0.d0
!
            vud=-2.d0*xlet(indexf)*
     &           (gradpel(1,iel2)*xxj(1,indexf)+
     &            gradpel(2,iel2)*xxj(2,indexf)+
     &            gradpel(3,iel2)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(5,i)=vel(iel2,4)/(r*vfa(0,i))
               cycle
            endif
!     
            phic=1.d0+vcd/vud
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(5,i)=vel(iel2,4)
            elseif(phic.le.1.d0/6.d0) then
               vfa(5,i)=3.d0*vel(iel2,4)-2.d0*vel(iel1,4)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(5,i)=3.d0*vel(iel2,4)/4.d0+vel(iel1,4)/4.d0+vud/8.d0
            else
               vfa(5,i)=vel(iel2,4)/3.d0+2.d0*vel(iel1,4)/3.d0
            endif
            vfa(5,i)=vfa(5,i)/(r*vfa(0,i))
         endif
!
      enddo
!            
      return
      end
