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
!
! >    \brief Generate local transformation matrix \f$ T_e^{lin} \f$ for quad-lin mortar method
! >    see phd-thesis Sitzmann equation (4.4)
! >    Author:Saskia Sitzmann
! >
! > @param   [in]     ipkon       pointer into field kon
! > @param   [in]     kon         Field containing the connectivity of the elements in succesive order
! > @param   [in]     lakon       element label
! > @param   [in]     islavsurf   islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
! > @param   [out]    contr       field containing T_e contributions for current face
! > @param   [out]    icontr1     (i)  row  of contribution(i)
! > @param   [out]    icontr2     (i)  column of contribution(i)
! > @param   [out]    icounter    counter variable for contr
! > @param   [in]     lface        current slave face
! >
      subroutine createtele_lin(ipkon,kon,lakon,&
           islavsurf,&
           contr,icontr1,icontr2,icounter,lface)
      !
      !     Generate local transformation matrix T
      !
      !     Author: Sitzmann,Saskia ;
      !
      implicit none
      !
      logical debug
      !
      character*8 lakon(*)
      !
      integer ipkon(*),kon(*),konl(20),islavsurf(2,*),&
           lface,icounter,icontr1(*),icontr2(*),j,nope,&
           ifaces,nelems,jfaces,m,nopes,&
           ifac,getiface,lnode(2,8),modf,idummy
      !
      real*8 contr(*),alpha
      !
      debug=.false.
      alpha=1.0/2.0
      icounter=0
      ifaces = islavsurf(1,lface)
      nelems = int(ifaces/10)
      jfaces = ifaces - nelems*10
      call getnumberofnodes(nelems,jfaces,lakon,nope,&
           nopes,idummy)
      do j=1,nope
         konl(j)=kon(ipkon(nelems)+j)
      enddo
      do m=1,nopes
         ifac=getiface(m,jfaces,nope)
         lnode(1,m)=konl(ifac)
      enddo
      if(nopes.eq.8) then
         do j=1,4
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo
         do j=5,8
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo
         do j=1,4
            icounter=icounter+1
            contr(icounter)=alpha
            icontr1(icounter)=lnode(1,j+4)
            icontr2(icounter)=lnode(1,j)
            icounter=icounter+1
            contr(icounter)=alpha
            icontr1(icounter)=lnode(1,modf(4,j-1)+4)
            icontr2(icounter)=lnode(1,j)
         enddo           
      elseif(nopes.eq.4) then
         do j=1,4
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo                 
      elseif(nopes.eq.6) then
         do j=1,3
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo
         do j=4,6
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo
         do j=1,3
            icounter=icounter+1
            contr(icounter)=alpha
            icontr1(icounter)=lnode(1,j+3)
            icontr2(icounter)=lnode(1,j)
            icounter=icounter+1
            contr(icounter)=alpha
            icontr1(icounter)=lnode(1,modf(3,j-1)+3)
            icontr2(icounter)=lnode(1,j)
         enddo
      else
         do j=1,3
            icounter=icounter+1
            contr(icounter)=1.0
            icontr1(icounter)=lnode(1,j)
            icontr2(icounter)=lnode(1,j)
         enddo
      endif 
      !
      if(debug)then
         write(*,*) 'createtele: contri,iscontr,imcontr',lface
         do j=1, icounter
            write(*,*)contr(j),icontr1(j),icontr2(j)
         enddo             
      endif
      !
      return
      end
      
