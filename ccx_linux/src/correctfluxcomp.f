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
      subroutine correctfluxcomp(nef,ipnei,neifa,neiel,flux,vfa,advfa,
     &  area,vel,alet,ielfa,ale,ifabou,ielmatf,mi,shcon,
     &  ntmat_,nefa,nefb,xxnj,gradpcfa,inlet)
!
!     correction of v due to the balance of mass
!     the correction is in normal direction to the face
!
!     bc:
!     outflow, p known: diffusion (subsonic)
!              p unknown: convection (supersonic)
!     inflow, p known: none
!             p unknown: convection (subsonic)
!
      implicit none
!
      integer i,nef,indexf,ipnei(*),neifa(*),neiel(*),ielfa(4,*),
     &  iel,ifa,ifabou(*),mi(*),ielmatf(mi(3),*),ntmat_,imat,indexb,
     &  nefa,nefb,inlet(*)
!
      real*8 flux(*),vfa(0:7,*),advfa(*),area(*),vel(nef,0:7),alet(*),
     &  ale(*),r,xflux,shcon(0:3,ntmat_,*),xxnj(3,*),gradpcfa(3,*)
!
!
!
      do i=nefa,nefb
c         totflux=0.d0
         imat=ielmatf(1,i)
         r=shcon(3,1,imat)
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            iel=neiel(indexf)
            xflux=flux(indexf)
            if(xflux.ge.0.d0) then
!
!              outflowing flux
!
               if(iel.gt.0) then
!
!                 internal face (velocity and density contribution)
!
                  flux(indexf)=flux(indexf)+vfa(5,ifa)*advfa(ifa)*
     &                 ((vel(i,4)-vel(iel,4))*alet(indexf)
     &                 -(gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                   gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                   gradpcfa(3,ifa)*xxnj(3,indexf)))
     &                        +flux(indexf)*vel(i,4)/
     &                           (vfa(5,ifa)*r*vfa(0,ifa))
               else
!
!                 external face
!
                  if(ielfa(3,ifa).le.0) then
                     indexb=-ielfa(2,ifa)
                     if(indexb.gt.0) then
                        if(ifabou(indexb+4).ne.0) then
!
!                          outflow, pressure known: diffusion term
!                          (typical subsonic outlet)
!
                           flux(indexf)=flux(indexf)+vfa(5,ifa)*
     &                        advfa(ifa)*(vel(i,4)*ale(indexf)
     &                             -(gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                               gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                               gradpcfa(3,ifa)*xxnj(3,indexf)))
                        else
!
!                          outflow, pressure unknown: convection term
!
                           flux(indexf)=flux(indexf)*(1.d0+vel(i,4)/
     &                          (vfa(5,ifa)*r*vfa(0,ifa)))
                        endif
                     else
!
!                       outflow, pressure unknown: convection term
!
                        flux(indexf)=flux(indexf)*(1.d0+vel(i,4)/
     &                       (vfa(5,ifa)*r*vfa(0,ifa)))
                     endif
                  endif
               endif
            else
!
!              inflowing flux
!
               if(iel.gt.0) then
!
!                 internal face (velocity and density contribution)
!
                  flux(indexf)=flux(indexf)+vfa(5,ifa)*advfa(ifa)
     &                           *((vel(i,4)-vel(iel,4))*alet(indexf)
     &                            -(gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                              gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                              gradpcfa(3,ifa)*xxnj(3,indexf)))
     &                        +flux(indexf)*vel(iel,4)/
     &                           (vfa(5,ifa)*r*vfa(0,ifa))
               else
!
!                 external face
!
                  indexb=-ielfa(2,ifa)
                  if(indexb.gt.0) then
                     if((inlet(ifa).ne.0).and.
     &                  (ifabou(indexb+4).eq.0)) then
!
!                       inlet and pressure unknown
!                       typical subsonic inlet conditions
!                       
                        flux(indexf)=flux(indexf)*(1.d0+vel(iel,4)/
     &                        (vfa(5,ifa)*r*vfa(0,ifa)))
                     endif
                  endif
               endif
            endif
         enddo
      enddo
! 
      return
      end
