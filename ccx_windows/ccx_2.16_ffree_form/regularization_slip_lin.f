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
! > regularization function for tangential mortar contact, linear regularization
! > see phd-thesis Sitzmann Chapter 3.2.2., semi-smooth Newton for tangential contact
! >
! > @param [in] lambdan                 normal contact pressure
! > @param [in] utilt                     \f$ \tilde{u}_\tau=\tau \cdot (\hat{u}-\bar{u}) \f$
! > @param [in] bp                     friction bound
! > @param [in] atauinv                     inverse of tangential stiffness
! > @param [out] resreg                     evaluated regularization function
! > @param [in] divmode                      indicates whether funtion (==0) or derivate (==1) should be called
! > @param [in] regmode                      not used
! > @param [in] islavact                     (i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node)
! > @param [in] lambdat                     shear stress
! > @param [in] lambdainit              not used
! > @param [in] lambdatilt              shear stress-shear stress at start of the increment
! > @param [in] constantt              parameter in semi-smooth Newton
! > @param [in] debug                     debug output flag
! > @param [in] inode                     slave node number
! > @param [in] n                     slave normal
! > @param [in] n2                     transformed slave normal
! > @param [in] t                     slave tangent
! > @param [in] that                     transformed slave tangent
! > @param [in] mu                     friction coefficient
! > @param [out] rslip                     matrix used in semi-smooth Newton
! > @param [out] ltslip                     matrix used in semi-smooth Newton
! > @param [out] ltu                     vector used in semi-smooth Newton
! >
      subroutine regularization_slip_lin(lambdan,&
           utilt,bp,atauinv,resreg,divmode,regmode,&
           islavact,lambdat,lambdainit,lambdatilt,&
           constantt,debug,&
           inode,n,n2,t,that,mu,rslip,ltslip,ltu)
      !
      !     regularization function of tangential contact
      !     Author: Saskia Sitzmann
      implicit none
      !
      integer i,j,k,l,m,niwan,divmode,regmode,inode,&
           yielded,imodification,islavact(*),debug
      !
      real*8 bp,atauinv,sb,rb,kiwan,constantt,&
           taulim,rbmax,nut,resreg(2),&
           rev,help(2),nhelp,scal,n(*),t(*),lambdan,&
           lpt(2),lptini(2),d(2),nd,&
           lp(4),mp(4),fp(4),ep,mu,&
           ltslip(*),rslip(*),ltu(*),utilt(*),&
           lambdat(*),lambdainit(*),lambdatilt(*),nlambdat,&
           nlmhat,lmhat(2),alpha,beta,delta,det,mp2(4),&
           mp3(4),vp(2),hp(2),rp(2),that(*),n2(*)
      !
      imodification=1
      !
      resreg(1)=0.0
      resreg(2)=0.0
      !
      constantt=min(constantt,1.0/atauinv);       
      if(divmode.eq.0)then
         !
         lmhat(1)=lambdat(1)&
              +constantt*(utilt(1)-atauinv*lambdatilt(1))
         lmhat(2)=lambdat(2)&
              +constantt*(utilt(2)-atauinv*lambdatilt(2))
         nlmhat=sqrt(lmhat(1)*lmhat(1)+lmhat(2)*lmhat(2))
         if(islavact(inode).eq.1)then 
            resreg(1)=lmhat(1)
            resreg(2)=lmhat(2)
         else if(islavact(inode).eq.2)then
            resreg(1)=bp*lmhat(1)/nlmhat
            resreg(2)=bp*lmhat(2)/nlmhat           
         endif
      elseif(divmode.eq.1)then
         if(islavact(inode).eq.1)then
            !
            do l=1,3
               do m=1,2   
                  ltslip((m-1)*3+l)=t(3*(m-1)+l)
               enddo
            enddo       
            do l=1,3
               do m=1,2
                  rslip((m-1)*3+l)=atauinv*that((m-1)*3+l)&
                       -mu*(utilt(m)-atauinv*lambdatilt(m))*n2(l)/bp
               enddo      
            enddo   
            ltu(1)=atauinv*lambdat(1)
            ltu(2)=atauinv*lambdat(2)            
         else if(islavact(inode).eq.2)then
            lmhat(1)=lambdat(1)&
                 +constantt*(utilt(1)-atauinv*lambdatilt(1))
            lmhat(2)=lambdat(2)&
                 +constantt*(utilt(2)-atauinv*lambdatilt(2))
            nlmhat=sqrt(lmhat(1)*lmhat(1)+lmhat(2)*lmhat(2))
            nlambdat=sqrt(lambdat(1)*lambdat(1)+lambdat(2)*lambdat(2))
            ep=bp/nlmhat
            if(mu.gt.1.E-10)then            
               if(imodification.eq.1)then           
                  do l=1,2             
                     do m=1,2               
                        fp(2*(m-1)+l)=(1.0/(max(bp,nlambdat)*nlmhat))&
                             *lambdat(m)*lmhat(l)
                     enddo
                  enddo
               else if(imodification.eq.2)then
                  fp(1)= (1.0/(max(bp,nlambdat)*nlmhat))&
                       *lambdat(1)*lmhat(1)
                  fp(2)= (1.0/(max(bp,nlambdat)*nlmhat)*2)&
                       *(lambdat(2)*lmhat(1)+lambdat(1)*lmhat(2))
                  fp(3)= (1.0/(max(bp,nlambdat)*nlmhat)*2)&
                       *(lambdat(2)*lmhat(1)+lambdat(1)*lmhat(2))
                  fp(4)= (1.0/(max(bp,nlambdat)*nlmhat))&
                       *lambdat(2)*lmhat(2)
               else if(imodification.eq.3)then
                  fp(1)= (1.0/(nlmhat*nlmhat))*lmhat(1)*lmhat(1)
                  fp(2)= (1.0/(nlmhat*nlmhat))*lmhat(2)*lmhat(1)
                  fp(3)= (1.0/(nlmhat*nlmhat))*lmhat(1)*lmhat(2)
                  fp(4)= (1.0/(nlmhat*nlmhat))*lmhat(2)*lmhat(2)             
               endif 
               mp(1)=1.0
               mp(2)=0.0
               mp(3)=0.0
               mp(4)=1.0
               do l=1,2
                  do m=1,2
                     mp(2*(m-1)+l)=(mp(2*(m-1)+l)-fp(2*(m-1)+l))*ep
                  enddo
               enddo         
               alpha=(lambdat(1)*lmhat(1)+lambdat(2)*lmhat(2))&
                    /(nlambdat*nlmhat)
               delta=min(1.0,nlambdat/bp)
               if(imodification.eq.1)then
                  if(alpha.lt.0.0)then
                     beta=1.0/(1.0-alpha*delta)
                  else
                     beta=1.0
                  endif
               else if(imodification.eq.2)then
                  beta=2/(2-(alpha-1)*delta)
               else if(imodification.eq.3)then
                  beta=1.0
               endif  
               !              /* H_p==mp2 */
               mp2(1)=1.0
               mp2(2)=0.0
               mp2(3)=0.0
               mp2(4)=1.0
               do l=1,4
                  mp2(l)=mp2(l)-beta*(1.0-constantt*atauinv)*mp(l) 
               enddo 
               det=mp2(1)*mp2(4)-mp2(2)*mp2(3)
               mp3(1)=(1.0/det)*mp2(4)
               mp3(2)=-(1.0/det)*mp2(2)
               mp3(3)=-(1.0/det)*mp2(3)
               mp3(4)=(1.0/det)*mp2(1)           
               vp(1)=(1.0/(nlmhat))*(mp3(1)*lmhat(1)+mp3(2)*lmhat(2))
               vp(2)=(1.0/(nlmhat))*(mp3(3)*lmhat(1)+mp3(4)*lmhat(2)) 
               !     /* solving for du -> need more arrays*/
               lp(1)=constantt*(mp3(1)*mp(1)+mp3(2)*mp(3))
               lp(2)=constantt*(mp3(1)*mp(2)+mp3(2)*mp(4))
               lp(3)=constantt*(mp3(3)*mp(1)+mp3(4)*mp(3))
               lp(4)=constantt*(mp3(3)*mp(2)+mp3(4)*mp(4))              
               do l=1,3
                  do m=1,2   
                     ltslip((m-1)*3+l)=lp(2*(m-1)+1)*t(l)&
                          +lp((m-1)*2+2)*t(3+l)
                  enddo
               enddo       
               do l=1,3
                  do m=1,2
                     rslip((m-1)*3+l)=that((m-1)*3+l)&
                          -mu*vp(m)*n2(l)
                  enddo      
               enddo  
               ltu(1)=-((1.0/constantt)-atauinv)&
                    *(lp(1)*lambdat(1)+lp(2)*lambdat(2))
               ltu(2)=-((1.0/constantt)-atauinv)&
                    *(lp(3)*lambdat(1)+lp(4)*lambdat(2))
               if(debug.eq.1)then
                  write(*,*) 't1',t(1),t(2),t(3)
                  write(*,*) 't1',t(4),t(5),t(6)
                  write(*,*)'ep',ep,'vp',vp(1),vp(2)
                  write(*,*) 'fp',fp(1),fp(2),fp(3),fp(4)
                  write(*,*) 'mp',mp(1),mp(2),mp(3),mp(4)
                  write(*,*) 'mp2',mp2(1),mp2(2),mp2(3),mp2(4)
                  write(*,*) 'mp3',mp3(1),mp3(2),mp3(3),mp3(4)
                  write(*,*)'rsn',mu*vp(1)*n2(1),mu*vp(1)*n2(1),&
                       mu*vp(1)*n2(1)
                  write(*,*)'rsn',mu*vp(2)*n2(1),mu*vp(2)*n2(1),&
                       mu*vp(2)*n2(1)
                  write(*,*)'det',det,beta,alpha,delta
                  write(*,*) 'lp',lp(1),lp(2),lp(3),lp(4)
               endif 
            else
               !     ///contact tie without friction
               do l=1,3
                  do m=1,2
                     rslip((m-1)*3+l)=t((m-1)*3+l)       
                  enddo
               enddo
               do l=1,3
                  do m=1,2   
                     ltslip((m-1)*3+l)=0.0
                  enddo
               enddo 
               ltu(1)=0.0
               ltu(2)=0.0
            endif
         endif
      endif
      !
      return
      end
      
