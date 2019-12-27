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
      subroutine free_disc_pumping(node1,node2,nodem,nelem,lakon,kon,&
              ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,&
              nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,&
              nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,&
              iaxial,iplausi)
      !
      !     Free-Disc-Pump-Flow
      !
      implicit none
      !
      logical identity
      character*8 lakon(*)
      character*81 set(*)
      !
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,&
           ielprop(*),nodef(*),idirf(*),index,iflag,iaxial,&
           ipkon(*),kon(*),mi(*),nrhcon(*),ntmat_,nshcon(*),iplausi
      !
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),ttime,time,&
           physcon(*),co(3,*),vold(0:mi(2),*),r,cp,dvi,&
           shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*)
      !
      intent(in) node1,node2,nodem,nelem,lakon,kon,ipkon,&
              nactdog,ielprop,prop,iflag,v,&
              cp,r,physcon,dvi,set,mi,&
              ttime,time,iaxial
      !
      intent(inout) identity,xflow,idirf,nodef,numf,f,df,iplausi
      !
      if (iflag.eq.0) then
         identity=.true.
         !
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
      !
      elseif ((iflag.eq.1).or.(iflag.eq.2).or.(iflag.eq.3))then
      !
      !    User defined flow element
      !
      endif
      !
      return
      end
      

