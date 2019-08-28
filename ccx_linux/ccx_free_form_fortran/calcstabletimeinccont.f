!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine calcstabletimeinccont(ne,lakon,kon,ipkon,mi,&
        ielmat,elcon,mortar,adb,alpha,nactdof,springarea,&
        ne0,ntmat_,ncmat_,dtcont)
      !
      !     Calculates the critical time increment (CTI) for contact
      !     spring elements based on the Courant
      !     Criterion for Explicit Dynamics calculations.
      !
      implicit none
      !
      character*8 lakon(*),lakonl
      !
      integer j,ne,nope,kon(*),ipkon(*),indexe,nelem,&
        ncmat_,ntmat_,mi(*),ielmat(mi(3),*),imat,mortar,&
        indexn,nactdof(0:mi(2),*),ne0,nopem
      !
      real*8 elcon(0:ncmat_,ntmat_,*),safefac,xk,adb(*),alpha,bet,gam,&
        critom,damping,springms,springmm,springfac,dtcont,xmacont,&
        springarea(2,*),areaslav
      !
      dtcont=1.d30
      safefac=0.80d0
      !
      xmacont=0.0d0
      !
      damping=0.d0
      !
      bet=(1.d0-alpha)*(1.d0-alpha)/4.d0
      gam=0.5d0-alpha
      !
      !     Omega Critical
      !     Om_cr=dt*freq_max
      !
      critom=dsqrt(damping*damping*(1.d0+2.d0*alpha*(1.d0-gam))&
           *(1.d0+2.d0*alpha*(1.d0-gam))&
          +   2.d0*(gam+2.d0*alpha*(gam-bet)))
      critom=0.98d0*(-damping*(1.d0+2.d0*alpha*(1.d0-gam))+critom)&
           /(gam+2.d0*alpha*(gam-bet)) !eq 25 miranda
      !
      !     ** DO per element
      !
      do nelem=ne0+1,ne
         indexe=ipkon(nelem)
         lakonl=lakon(nelem)
         imat=ielmat(1,nelem)
               !
               xk=elcon(2,1,imat)
               !
               springmm=0.0d0
               springms=0.0d0
               !
               if(mortar.eq.0) then  
                  !
                  !                 node-to-face
                  !
                  !                 nope is the total number of nodes:
                  !                 master nodes+1 slave node
                  !                 (notice -47 instead of -48)
                  !
                  nope=ichar(lakonl(8:8))-47
                  springfac=1.d0  
                  !
                  do indexn=1,nope                   
                     xmacont=0.0d0
                     !
                     do j=1,3           
                       if(nactdof(j,kon(indexe+indexn)).gt.0)then
                           xmacont=max(xmacont,&
                               adb(nactdof(j,kon(indexe+indexn))))
                        endif
                     enddo
                     !
                     if(indexn.eq.nope)then
                        springms=springms+xmacont
                        !
                        !                       mean mass at master nodes
                        !
                        springmm=springmm/(nope-1.0d0)    
                     else
                        springmm=springmm+xmacont    
                     endif        
                  enddo
                  !
                  areaslav=springarea(1,kon(indexe+nope+1))
                  xk=xk*areaslav
                  !
                  !                 checking whether mass or slave side is constrained
                  !
                  if((springmm.le.0.d0).and.(springms.le.0.d0)) then
                     cycle
                  elseif(springmm.le.0.d0) then
                     springmm=springms
                  elseif(springms.le.0.d0) then
                     springms=springmm
                  endif
               !
               elseif(mortar.eq.1) then
                  !
                  !                 face-to-face
                  !
                  nopem=ichar(lakonl(8:8))-48
                  springfac=0.1d0
                  !
                  do indexn=1,kon(indexe) 
                     !
                     xmacont=0.0d0
                     do j=1,3           
                        if(nactdof(j,kon(indexe+indexn)).gt.0)then
                           xmacont=max(xmacont,&
                                adb(nactdof(j,kon(indexe+indexn))))
                        endif
                     enddo
                     if(indexn.gt.nopem)then !slave
                        springms=springms+xmacont    
                     else !master
                        springmm=springmm+xmacont    
                     endif        
                  enddo
                  !
                  !           mean mass at slave and master nodes
                  !
                  springms=springms/(kon(indexe)-nopem)
                  springmm=springmm/nopem
                  !
                  areaslav=springarea(1,kon(1+indexe+kon(indexe)))
                  xk=xk*areaslav
                  !
                  !                 checking whether mass or slave side is constrained
                  !
                  if((springmm.le.0.d0).and.(springms.le.0.d0)) then
                     cycle
                  elseif(springmm.le.0.d0) then
                     springmm=springms
                  elseif(springms.le.0.d0) then
                     springms=springmm
                  endif
               endif
               !
               !              linear pressure-overclosure
               !
               if(int(elcon(3,1,imat)).eq.2) then
                  !
                  springmm=springmm/2.0d0
                  springms=springms/2.0d0      
                  !
                  dtcont =&
                       min(dtcont,&
                       springfac*critom*dsqrt((springmm*springms)/&
                       ((springmm+springms)*xk)))
               !
               else
                  write(*,*) '*ERROR in calcstabletimeinccont:'
                  write(*,*) '       in explicit dynamic calculations'
                  write(*,*) '       only linear pressure-overclosure'
                  write(*,*) '       is allowed'
                  call exit(201)
               endif
      enddo
      !     ** ENDDO per element
      !
      dtcont=dtcont* safefac
      !
      return
      end
