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
      subroutine extrapol_dpel5(ielfa,ipnei,vel,xlet,gradpcfa,xxj,
     &  nef,nfacea,nfaceb,ncfd)
!
      implicit none
!
      integer ielfa(4,*),ipnei(*),nef,nfacea,nfaceb,i,k,iel1,iel2,
     &  indexf,ncfd
!
      real*8 vel(nef,0:7),xlet(*),gradpcfa(3,*),xxj(3,*),dd
!
!
!
!     correct the facial pressure gradients:
!     Moukalled et al. p 289
!
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            iel1=ielfa(1,i)
            indexf=ipnei(iel1)+ielfa(4,i)
            dd=(vel(iel2,4)-vel(iel1,4))/xlet(indexf)
     &           -gradpcfa(1,i)*xxj(1,indexf)
     &           -gradpcfa(2,i)*xxj(2,indexf)
     &           -gradpcfa(3,i)*xxj(3,indexf)
            do k=1,ncfd
               gradpcfa(k,i)=gradpcfa(k,i)+dd*xxj(k,indexf)
            enddo
         endif
      enddo
!            
      return
      end
