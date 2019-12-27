!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2019 Guido Dhondt
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
!
! >
! >  \brief function extracting the contact constants
! >
! > @param [out] mu              friction coefficient
! > @param [out] regmode              regularization method in normal direction (=1 linear, =2 piece-wise liner,=3 exponential,=4 tied)
! > @param [out] regmodet       regularization method in tangnetial direction (=1 linear, =2 Iwan model)
! > @param [out] fkninv              inverse of normal stiffness \f$ \frac{1}{a_n} \f$
! > @param [out] fktauinv        inverse of tangential stiffness \f$ \frac{1}{a_\tau} \f$
! > @param [out] p0              parameter needed for exponential regularization
! > @param [out] beta              parameter needed for exponential regularization
! > @param [out] niwan              number of Iwan elements (1-10) for Iwan model
! >
      subroutine getcontactparams(mu,regmode,regmodet,fkninv,fktauinv,&
           p0,beta,tietol,elcon,itie,ncmat_,ntmat_,niwan)
       !     Author: Saskia Sitzmann
      
      implicit none
      !
      integer itie,imat,ncmat_,ntmat_,regmode,regmodet,niwan
      !
      real*8 mu,fkninv,fktauinv,p0,beta,tietol(3,*),&
           elcon(0:ncmat_,ntmat_,*)
      !
      itie=itie+1
      imat=int(tietol(2,itie))
      !
      if(ncmat_.lt.6)then
         mu=0.0
         fktauinv=0.0
         regmodet=1
         niwan=1
      else
         mu=elcon(6,1,imat)
         if(elcon(7,1,imat).le.0.0)then
            fktauinv=0.0
         else
            fktauinv=1.0/elcon(7,1,imat)
         endif
         regmodet=1
         if(elcon(8,1,imat).le.0.99)then
            regmodet=1
            niwan=1
         else
            regmodet=2
           niwan=int(elcon(8,1,imat))
            niwan=min(10,niwan)
         endif
         if(fktauinv.lt.1.e-8)then 
            regmodet=1
            niwan=1
         endif
      endif
      !
      !     exponential regularization
      !
      if(ncmat_.gt.2)then
         if(elcon(3,1,imat).gt.1.4 .and.&
              elcon(3,1,imat).lt.1.6 )then
            regmode=3
            p0=elcon(2,1,imat)
            beta=1.0/(elcon(1,1,imat))
            fkninv=0.0
            if(mu.gt.1.e-10)then
               write(*,*)'getcontactparams:'
               write(*,*)'exponential pressure overclosure',&
                    'with friction not yet supportes'
               call exit(201)
            endif
         !
         !     linear regularization
         else if(elcon(3,1,imat).gt.2.4 .and.&
                 elcon(3,1,imat).lt.2.6 )then
            regmode=1
            fkninv=1.0/elcon(2,1,imat)
            p0=0.0
            beta=0.0
         !
         !     piecewiese linear regularization
         else if(elcon(3,1,imat).gt.3.4 .and.&
                 elcon(3,1,imat).lt.3.6 )then
            regmode=2
            p0=0.0
            beta=0.0
            fkninv=0.0
         !
         !     tied contact
         else if(elcon(3,1,imat).gt.4.4 .and.&
                 elcon(3,1,imat).lt.4.6 )then
            regmode=4
            p0=0.0
            beta=0.0
            fkninv=0.0
            mu=0.0
            fktauinv=0.0
            regmodet=4
         else
            regmode=1
            fkninv=0.0
            p0=0.0
            beta=0.0
         endif
      else
         regmode=1
         fkninv=0.0
         p0=0.0
         beta=0.0
      endif 
      itie=itie-1
      !
      return
      end
      
