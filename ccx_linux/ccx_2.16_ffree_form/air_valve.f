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
      subroutine air_valve(node1,node2,nodem,nelem,lakon,kon,ipkon,&
              nactdog,identity,ielprop,prop,iflag,v,xflow,f,&
              nodef,idirf,df,cp,r,physcon,dvi,numf,set,co,vold,mi,&
              ttime,time,iaxial,iplausi)
      !
      !     proprietary
      !
      implicit none
      !
      logical identity
      character*8 lakon(*)
      character*81 set(*)
      !
      integer node1,node2,nodem,nelem,kon(*),ipkon(*),nactdog(0:3,*),&
        ielprop(*),iflag,nodef(5),idirf(5),numf,mi(*),iaxial,&
        iplausi
      !
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),cp,r,physcon(*),dvi,&
        ttime,time,co(3,*),vold(0:mi(2),*)
      !
      intent(in) node1,node2,nodem,nelem,lakon,kon,ipkon,&
              nactdog,ielprop,prop,iflag,v,cp,r,physcon,dvi,set,mi,&
              ttime,time,iaxial,co,vold
      !
      intent(inout) identity,xflow,idirf,nodef,numf,f,df,iplausi
      !
      return
      end
      
      
      
     
      
      
      
    
