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
! > regularization function for tangential mortar contact, iwan model
! > see phd-thesis Sitzmann Chapter 3.2.2., semi-smooth Newton for tangential contact
! >
! > @param [in] lambdan                 normal contact pressure
! > @param [in] utilt                     \f$ \tilde{u}_\tau=\tau \cdot (\hat{u}-\bar{u}) \f$
! > @param [in] bp                     friction bound
! > @param [in] atau2                     tangential stiffness
! > @param [out] resreg                     evaluated regularization function
! > @param [in] divmode                      indicates whether funtion (==0) or derivate (==1) should be called
! > @param [in] regmode                      selects used semi-Newton (==2 is active)
! > @param [in,out] lambdaiwan          Lagrange multiplier splitted to Iwan elements
! > @param [in,out] lambdaiwanini           Lagrange multiplier splitted to Iwan elements at start of increment
! > @param [in] inode                     slave node number
! > @param [in] n                     slave normal
! > @param [in] t                     slave tangent
! > @param [in] mu                     friction coefficient
! > @param [out] rslip                     matrix used in semi-smooth Newton
! > @param [out] ltslip                     matrix used in semi-smooth Newton
! > @param [out] ltu                     vector used in semi-smooth Newton
! > @param [out] yielded                     debugging parameter
! > @param [in] iit                     Newton-Raphson iteration
! > @param [in] debug                     debug output flag
! > @param [in] niwan                     number of iwan elements
! > @param [in] dut                     \f$ \Delta \tilde{u}_\tau \f$
! >
      subroutine regularization_slip_iwan(lambdan,&
           utilt,bp,atau2,resreg,divmode,regmode,lambdaiwan,&
           lambdaiwanini,inode,n,t,mu,rslip,ltslip,ltu,yielded,iit,&
           debug,niwan,dut)
      !
      !     regularization function of tangential contact
      !     implementation of iwan-model
      !     Author: Saskia Sitzmann
      implicit none
      !
      integer i,j,k,niwan,divmode,regmode,inode,yielded,iit,&
           debug,imodification
      !
      real*8 bp,atau2,sb,rb,kiwan,fstar(10),dut(*),&
           taulim,rbmax,nut,alpha(10),resreg(2),&
           rev,help(2),nhelp,scal,n(*),t(*),lambdan,&
           lpt(2),lptini(2),d(2),nd,lambdaiwan(3,niwan,*),&
           lambdaiwanini(3,niwan,*),lp(4),mp(4),fp(4),ep,mu,&
           ltslip(*),rslip(*),ltu(*),utilt(*),det,nlpt
      !
      !      number of iwan elements
      !      niwan=1
      !
      kiwan=atau2
      imodification=regmode
      do i=1,niwan
         alpha(i)=2.0*real(i)/(real(niwan)*real(niwan+1))
         fstar(i)=alpha(i)*bp*real(niwan)
      enddo
      yielded=0
      resreg(1)=0.0
      resreg(2)=0.0
      if(imodification==1)then
         if(divmode.eq.0)then
            !      update lambdaiwan
            do i=1,niwan
               lptini(1)=lambdaiwanini(1,i,inode)*t(1)&
                    +lambdaiwanini(2,i,inode)*t(2)&
                    +lambdaiwanini(3,i,inode)*t(3)
               lptini(2)=lambdaiwanini(1,i,inode)*t(4)&
                    +lambdaiwanini(2,i,inode)*t(5)&
                    +lambdaiwanini(3,i,inode)*t(6)
               d(1)=niwan*lptini(1)+kiwan*utilt(1)
               d(2)=niwan*lptini(2)+kiwan*utilt(2)
               nd=sqrt(d(1)*d(1)+d(2)*d(2))
               nhelp=sqrt(utilt(1)*utilt(1)+utilt(2)*utilt(2))
               if(debug.eq.1)then
                  write(*,*)'imodification',imodification
                  write(*,*)'lini',lptini(1),lptini(2)
                  write(*,*)'d',d(1),d(2)
               endif
               !      write(*,*)'l ',i,lptini(2),d(2)/real(niwan)
               if(nd.le.fstar(i))then
                  !      update ok
                  resreg(1)=resreg(1)+d(1)/real(niwan)
                  resreg(2)=resreg(2)+d(2)/real(niwan)
               else 
                  !      radial return mapping
                  yielded=yielded+1
                  d(1)=fstar(i)*d(1)/(nd)
                  d(2)=fstar(i)*d(2)/(nd)           
                  resreg(1)=resreg(1)+d(1)/real(niwan)
                  resreg(2)=resreg(2)+d(2)/real(niwan)
               endif
               if(debug.eq.1)then
                  write(*,*)'liwan',d(1)/real(niwan),d(2)/real(niwan)
               endif
            !            write(*,*)'lu',i,lptini(2),d(2)/real(niwan)
            enddo
         elseif(divmode.eq.1)then
            !      Newton iteration
            do i=1,6         
               rslip(i)=t(i) 
               ltslip(i)=0.0      
            enddo
            do i=1,2
               ltu(i)=0.0
            enddo
            do i=1,niwan
               lptini(1)=lambdaiwanini(1,i,inode)*t(1)&
                    +lambdaiwanini(2,i,inode)*t(2)&
                    +lambdaiwanini(3,i,inode)*t(3)
               lptini(2)=lambdaiwanini(1,i,inode)*t(4)&
                    +lambdaiwanini(2,i,inode)*t(5)&
                    +lambdaiwanini(3,i,inode)*t(6)
               d(1)=niwan*lptini(1)+kiwan*utilt(1)
               d(2)=niwan*lptini(2)+kiwan*utilt(2)
               nd=sqrt(d(1)*d(1)+d(2)*d(2))
               nhelp=sqrt(utilt(1)*utilt(1)+utilt(2)*utilt(2))
               resreg(1)=d(1)
               resreg(2)=d(2)
               if(nd.lt.fstar(i).or.(iit.eq.1))then
                  !      update ok
                  do j=1,3
                     do k=1,2
                        ltslip((k-1)*3+j)=ltslip((k-1)*3+j)&
                             +(kiwan/real(niwan))*t((k-1)*3+j)
                     enddo
                  enddo
                  ltu(1)=ltu(1)+d(1)/real(niwan)
                  ltu(2)=ltu(2)+d(2)/real(niwan)
               else 
                  !      radial return mapping
                  yielded=yielded+1
                  ep=fstar(i)/(nd*real(niwan))
                  fp(1)=d(1)*d(1)/(nd*nd)
                  fp(2)=d(1)*d(2)/(nd*nd)
                  fp(3)=d(2)*d(1)/(nd*nd)
                  fp(4)=d(2)*d(2)/(nd*nd)
                  mp(1)=-(ep*fp(1)-ep)
                  mp(2)=-ep*fp(2)
                  mp(3)=-ep*fp(3)   
                  mp(4)=-(ep*fp(4)-ep)
                  lp(1)=kiwan*mp(1)  
                  lp(2)=kiwan*mp(2)
                  lp(3)=kiwan*mp(3)
                  lp(4)=kiwan*mp(4)
                  if(debug.eq.1)then
                     write(*,*) 't1',t(1),t(2),t(3)
                     write(*,*) 't1',t(4),t(5),t(6)
                     write(*,*)'ep',ep
                     write(*,*) 'fp',fp(1),fp(2),fp(3),fp(4)
                     write(*,*) 'mp',mp(1),mp(2),mp(3),mp(4)
                     write(*,*) 'lp',lp(1),lp(2),lp(3),lp(4)
                     write(*,*)'rn',(d(1)/nd)*alpha(i)*mu*n(1),&
                          (d(1)/nd)*alpha(i)*mu*n(2),&
                          (d(1)/nd)*alpha(i)*mu*n(3)
                     write(*,*)'rn',(d(2)/nd)*alpha(i)*mu*n(1),&
                          (d(2)/nd)*alpha(i)*mu*n(2),&
                          (d(2)/nd)*alpha(i)*mu*n(3)
                     det=lp(1)*lp(4)-lp(2)*lp(3) 
                     write(*,*)'det',det
                  endif     
                  do j=1,3
                     do k=1,2
                        ltslip((k-1)*3+j)=ltslip((k-1)*3+j)&
                             +lp((k-1)*2+1)*t(j)&
                             +lp((k-1)*2+2)*t(j+3)
                        !
                        rslip((k-1)*3+j)=rslip((k-1)*3+j)&
                             -(d(k)/nd)*alpha(i)*mu*n(j)
                     enddo
                  enddo
               endif
            enddo      
         endif
      !
      else
         !      alternative Newton iteration
         if(debug.eq.1)then
            write(*,*)'imodification',imodification
         endif
         if(divmode.eq.0) then
            !      update lambdaiwan
            do i=1,niwan
               lptini(1)=lambdaiwanini(1,i,inode)*t(1)&
                    +lambdaiwanini(2,i,inode)*t(2)&
                    +lambdaiwanini(3,i,inode)*t(3)
               lptini(2)=lambdaiwanini(1,i,inode)*t(4)&
                    +lambdaiwanini(2,i,inode)*t(5)&
                    +lambdaiwanini(3,i,inode)*t(6)
               lpt(1)=lambdaiwan(1,i,inode)*t(1)&
                    +lambdaiwan(2,i,inode)*t(2)&
                 +lambdaiwan(3,i,inode)*t(3)
               lpt(2)=lambdaiwan(1,i,inode)*t(4)&
                    +lambdaiwan(2,i,inode)*t(5)&
                    +lambdaiwan(3,i,inode)*t(6)
               nlpt=sqrt(lpt(1)*lpt(1)+lpt(2)*lpt(2))
               d(1)=niwan*lptini(1)+kiwan*utilt(1)
               d(2)=niwan*lptini(2)+kiwan*utilt(2)
               nd=sqrt(d(1)*d(1)+d(2)*d(2))
               nhelp=sqrt(utilt(1)*utilt(1)+utilt(2)*utilt(2))
               if(debug.eq.1)then
                  write(*,*)'lini',lptini(1),lptini(2)
                  write(*,*)'d',d(1),d(2),nd
               endif
               !      write(*,*)'l ',i,lptini(2),d(2)/real(niwan)
               if(nd.le.fstar(i))then
                  !      update ok
                  d(1)=(d(1)+kiwan*dut(1))/real(niwan)
                  d(2)=(d(2)+kiwan*dut(2))/real(niwan)
                  lambdaiwan(1,i,inode)=d(1)*t(1)
                  lambdaiwan(2,i,inode)=d(1)*t(2)
                  lambdaiwan(3,i,inode)=d(1)*t(3)
                  lambdaiwan(1,i,inode)=lambdaiwan(1,i,inode)&
                       +d(2)*t(4)
                  lambdaiwan(2,i,inode)=lambdaiwan(2,i,inode)&
                       +d(2)*t(5)
                  lambdaiwan(3,i,inode)=lambdaiwan(3,i,inode)&
                       +d(2)*t(6)
                  resreg(1)=resreg(1)+d(1)
                  resreg(2)=resreg(2)+d(2)
               else 
                  !      radial return mapping
                  yielded=yielded+1
                  ep=fstar(i)/(nd*real(niwan))
                  fp(1)=lpt(1)*d(1)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(2)=lpt(1)*d(2)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(3)=lpt(2)*d(1)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(4)=lpt(2)*d(2)/(nd*max(fstar(i)/real(niwan),nlpt))
                  mp(1)=-(ep*fp(1)-ep)
                  mp(2)=-ep*fp(2)
                  mp(3)=-ep*fp(3)   
                  mp(4)=-(ep*fp(4)-ep)
                  lp(1)=kiwan*mp(1)  
                  lp(2)=kiwan*mp(2)
                  lp(3)=kiwan*mp(3)
                  lp(4)=kiwan*mp(4)
                  d(1)=mu*lambdan*d(1)*alpha(i)/nd&
                       +lp(1)*dut(1)+lp(2)*dut(2)
                  d(2)=mu*lambdan*d(2)*alpha(i)/nd&
                       +lp(3)*dut(1)+lp(4)*dut(2)
                  lambdaiwan(1,i,inode)=d(1)*t(1)
                  lambdaiwan(2,i,inode)=d(1)*t(2)
                  lambdaiwan(3,i,inode)=d(1)*t(3)
                  lambdaiwan(1,i,inode)=lambdaiwan(1,i,inode)&
                       +d(2)*t(4)
                  lambdaiwan(2,i,inode)=lambdaiwan(2,i,inode)&
                       +d(2)*t(5)
                  lambdaiwan(3,i,inode)=lambdaiwan(3,i,inode)&
                       +d(2)*t(6)
                  resreg(1)=resreg(1)+d(1)
                  resreg(2)=resreg(2)+d(2)
                  
               endif
               if(debug.eq.1)then
                  write(*,*)'liwan',d(1),d(2)
               endif
            enddo
         elseif(divmode.eq.1)then
            !      Newton iteration
            do i=1,6         
               rslip(i)=t(i) 
               ltslip(i)=0.0      
            enddo
            do i=1,2
               ltu(i)=0.0
            enddo
            do i=1,niwan
               lptini(1)=lambdaiwanini(1,i,inode)*t(1)&
                    +lambdaiwanini(2,i,inode)*t(2)&
                    +lambdaiwanini(3,i,inode)*t(3)
               lptini(2)=lambdaiwanini(1,i,inode)*t(4)&
                    +lambdaiwanini(2,i,inode)*t(5)&
                    +lambdaiwanini(3,i,inode)*t(6)
               lpt(1)=lambdaiwan(1,i,inode)*t(1)&
                    +lambdaiwan(2,i,inode)*t(2)&
                    +lambdaiwan(3,i,inode)*t(3)
               lpt(2)=lambdaiwan(1,i,inode)*t(4)&
                    +lambdaiwan(2,i,inode)*t(5)&
                    +lambdaiwan(3,i,inode)*t(6)
               d(1)=niwan*lptini(1)+kiwan*utilt(1)
               d(2)=niwan*lptini(2)+kiwan*utilt(2)
               nd=sqrt(d(1)*d(1)+d(2)*d(2))
               nlpt=sqrt(lpt(1)*lpt(1)+lpt(2)*lpt(2))
               nhelp=sqrt(utilt(1)*utilt(1)+utilt(2)*utilt(2))
               resreg(1)=d(1)
               resreg(2)=d(2)
               if(nd.lt.fstar(i).or.(iit.eq.1))then
                  !      update ok
                  do j=1,3
                     do k=1,2
                        ltslip((k-1)*3+j)=ltslip((k-1)*3+j)&
                             +(kiwan/real(niwan))*t((k-1)*3+j)
                        !      check for niwan>1
                        rslip((k-1)*3+j)=rslip((k-1)*3+j)&
                             +(mu*(lpt(k)-(d(k)/real(niwan)))/bp)*n(j)
                     enddo
                  enddo
                  ltu(1)=ltu(1)+lpt(1)
                  ltu(2)=ltu(2)+lpt(2)
               else 
                  !      radial return mapping
                  yielded=yielded+1
                  ep=fstar(i)/(nd*real(niwan))
                  fp(1)=lpt(1)*d(1)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(2)=lpt(1)*d(2)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(3)=lpt(2)*d(1)/(nd*max(fstar(i)/real(niwan),nlpt))
                  fp(4)=lpt(2)*d(2)/(nd*max(fstar(i)/real(niwan),nlpt))
                  mp(1)=-(ep*fp(1)-ep)
                  mp(2)=-ep*fp(2)
                  mp(3)=-ep*fp(3)   
                  mp(4)=-(ep*fp(4)-ep)
                  lp(1)=kiwan*mp(1)  
                  lp(2)=kiwan*mp(2)
                  lp(3)=kiwan*mp(3)
                  lp(4)=kiwan*mp(4)
                  if(debug.eq.1)then
                     write(*,*) 't1',t(1),t(2),t(3)
                     write(*,*) 't1',t(4),t(5),t(6)
                     write(*,*)'ep',ep
                     write(*,*) 'fp',fp(1),fp(2),fp(3),fp(4)
                     write(*,*) 'mp',mp(1),mp(2),mp(3),mp(4)
                     write(*,*) 'lp',lp(1),lp(2),lp(3),lp(4)
                     write(*,*)'rn',(d(1)/nd)*alpha(i)*mu*n(1),&
                          (d(1)/nd)*alpha(i)*mu*n(2),&
                          (d(1)/nd)*alpha(i)*mu*n(3)
                     write(*,*)'rn',(d(2)/nd)*alpha(i)*mu*n(1),&
                          (d(2)/nd)*alpha(i)*mu*n(2),&
                          (d(2)/nd)*alpha(i)*mu*n(3)
                     det=lp(1)*lp(4)-lp(2)*lp(3) 
                     write(*,*)'det',det
                  endif     
                  do j=1,3
                     do k=1,2
                        ltslip((k-1)*3+j)=ltslip((k-1)*3+j)&
                             +lp((k-1)*2+1)*t(j)&
                             +lp((k-1)*2+2)*t(j+3)
                        !
                        rslip((k-1)*3+j)=rslip((k-1)*3+j)&
                             -(d(k)/nd)*alpha(i)*mu*n(j)
                     enddo
                  enddo
               endif
            enddo
         endif
         
      endif
      !
      return
      end
      
