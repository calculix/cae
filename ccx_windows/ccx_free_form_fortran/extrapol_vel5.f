!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine extrapol_vel5(ielfa,ipnei,vel,xlet,gradvfa,xxj,&
        nef,nfacea,nfaceb)
      !
      !     correct the facial velocity gradients:
      !     Moukalled et al. p 289
      !
      implicit none
      !
      integer ielfa(4,*),ipnei(*),nef,nfacea,nfaceb,i,k,l,indexf,iel1,&
        iel2
      !
      real*8 vel(nef,0:7),xlet(*),gradvfa(3,3,*),xxj(3,*),dd
      !
      intent(in) ielfa,ipnei,vel,xlet,xxj,&
        nef,nfacea,nfaceb
      !
      intent(inout) gradvfa
      !
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            iel1=ielfa(1,i)
            indexf=ipnei(iel1)+ielfa(4,i)
            do k=1,3
               dd=(vel(iel2,k)-vel(iel1,k))/xlet(indexf)&
                    -gradvfa(k,1,i)*xxj(1,indexf)&
                    -gradvfa(k,2,i)*xxj(2,indexf)&
                    -gradvfa(k,3,i)*xxj(3,indexf)
               do l=1,3
                  gradvfa(k,l,i)=gradvfa(k,l,i)+dd*xxj(l,indexf)
               enddo
            enddo
         endif
      enddo
      !
      return
      end
