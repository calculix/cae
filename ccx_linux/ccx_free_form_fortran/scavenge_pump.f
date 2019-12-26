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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine scavenge_pump(node1,node2,nodem,nelem,lakon,kon,ipkon,&
              nactdog,identity,ielprop,prop,iflag,v,xflow,f,&
              nodef,idirf,df,cp,r,physcon,dvi,numf,set,ntmat_,mi,&
              ttime,time,iaxial,iplausi)
      !
      !     scavenge pump element
      !
      !     author: Yannick Muller
      !
      implicit none
      !
      logical identity
      character*8 lakon(*)
      character*81 set(*)
      !
      integer nelem,nactdog(0:3,*),numf,node1,node2,nodem,&
           ielprop(*),nodef(5),idirf(5),index,iflag,&
           ipkon(*),kon(*),mi(*),ntmat_,iaxial,iplausi
      !
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),kappa,cp,physcon(*)&
           ,dvi,R,ttime,time
      !
      intent(in) node1,node2,nodem,nelem,lakon,kon,ipkon,&
              nactdog,ielprop,prop,iflag,v,cp,r,physcon,dvi,set,mi,&
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
      elseif (iflag.eq.1) then      
         if(v(1,nodem).ne.0.d0) return
      !
      elseif (iflag.eq.2) then
      !
      elseif (iflag.eq.3) then
      !
      endif
      return
      end
