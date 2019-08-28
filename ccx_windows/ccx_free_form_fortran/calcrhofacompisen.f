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
      subroutine calcrhofacompisen(nface,vfa,shcon,ielmat,ntmat_,&
        mi,ielfa,cvfa,velo,nef)
      !
      !     calculation of the density at the face centers
      !     (compressible fluids)
      !
      implicit none
      !
      integer nface,i,imat,ntmat_,mi(*),nef,velo(nef,0:7),&
        ielmat(mi(3),*),ielfa(4,*),j
      !
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),cvfa(*) 
      !
      do i=1,nface
         t1l=vfa(0,i)
         j=ielfa(1,i)
         !
         !        take the material of the first adjacent element
         !
         imat=ielmat(1,ielfa(1,i))
         vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*&
             !      &       (10.5d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
             (5.98696d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
      !      &       (1.41827d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
      !          vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*
      !      &       (velo(j,0)+(velo(j,1)**2+velo(j,2)**2+velo(j,3)**2)/2.d0-
      !      &(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
      enddo
      !
      return
      end
