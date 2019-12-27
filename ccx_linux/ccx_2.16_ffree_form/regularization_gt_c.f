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
! >
! > \brief regularization function for tangential contact mortar (old, only for contactstress2)
! >
! > @param [in]       lambdatt       lambdatilde_tau=lambda_tau-bar{lambda}_tau
! > @param [in]       divmode        indicates whether funtion or derivate
! >                                    should be called
! >                                  =0 function called
! >                                  =1 derivative called
! > @param [in]       regmode               selects regularization funtion
! >                                  =1 perturbed Lagrange
! > @param [out] gtc               result regularization function
! > @param [in]  atauinvloc      stiffness constant for perturbed Lagrange
! > @param [in]        p0              parameter for exponential regularization
! > @param [in]        beta              parameter for exponential regularization
! > @param [in]        elcon              material parameters
! > @param [in]        nelcon              (1,i) number of elastic constants for material i (2,i) number of temperature points
! > @param [in]        itie              current tie
! > @param [in]        ntmat_               maximum number of temperature data points for any material
! > @param [in]        plicon              isotropic hardening curve or points for pressure-overclosure=tabular
! > @param [in]        nplicon              isotropic hardening curve.
! > @param [in]        npmat_              maximum number of data points for plicon
! > @param [in]        ncmat_              maximum number of elastic material constants
! > @param [in]        tietol              (1,i) tie tolerance (2,i) contant interaction material definition
! >
      subroutine regularization_gt_c(lambdatt,divmode,regmode,&
           gtc,atauinvloc,p0,beta,elcon,nelcon,itie,ntmat_,&
           plicon,nplicon,npmat_,ncmat_,tietol)
      !
      !     regularization function for tangential contact
      !     Author: Saskia Sitzmann
      !
      implicit none
      !
      integer divmode,regmode,i,ndata,kode,npmat_,ncmat_,&
           itie,ntmat_,nplicon(0:ntmat_,*),nelcon(2,*),&
           imat
      !
      !
      real*8 lambdatt(*),gtc(*),pn_d(40),gn_d(40),aninv1(40),t(40),&
           beta,p0,atauinvloc,elconloc(21),plconloc(82),t1l,&
           elcon(0:ncmat_,ntmat_,*),plicon(0:2*npmat_,ntmat_,*),&
           tietol(2,*)

      kode=-51
      t1l=0.0
       !      imat=int(tietol(2,itie+1))
      
      !
      !     perturbed Lagrange
      !
      if(regmode.eq.1)then
         if(divmode.eq.0)then
            gtc(1)=atauinvloc*lambdatt(1)
            gtc(2)=atauinvloc*lambdatt(2)
         elseif(divmode.eq.1)then
            gtc(1)=atauinvloc
            gtc(2)=atauinvloc
         else
            write(*,*)'error in regularzation_gt_c.f!'
            call exit(201)
         endif
      else
         gtc(1)=0.0
         gtc(2)=0.0
      endif
      
      return
      end
