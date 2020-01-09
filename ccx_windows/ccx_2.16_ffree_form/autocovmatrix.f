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
      subroutine autocovmatrix(co,ad,au,jqs,irows,ndesi,nodedesi,&
        physcon)
      !
      !     calculates the values of the autocovariance matrix
      !
      implicit none
      !
      integer jqs(*),irows(*),ndesi,nodedesi(*),idof,j,jdof,node1,&
        node2
      !
      real*8 co(3,*),ad(*),au(*),physcon(*),dist,corrlength,sigma
      !
      corrlength=physcon(13)
      sigma=physcon(12)
      !
      do idof=1,ndesi
         ad(idof)=sigma*sigma
         do j=jqs(idof),jqs(idof+1)-1
            jdof=irows(j)
            node1=nodedesi(idof)
            node2=nodedesi(jdof)
            dist=dsqrt((co(1,node1)-co(1,node2))**2+&
                       (co(2,node1)-co(2,node2))**2+&
                       (co(3,node1)-co(3,node2))**2)
            !
            !           assign the value to the autocovariance matrix
            !
            au(j)=ad(idof)*dexp(-(dist/corrlength)**2)
         enddo
      enddo
      !
      return        
      end




