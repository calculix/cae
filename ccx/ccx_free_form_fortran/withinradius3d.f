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
      subroutine withinradius3d(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,n,&
        ir,r,nr,radius)
      !
      !     determines out of a cloud of n nodes with coordinates in
      !     (xo,yo,zo) those nodes which are within a distance "radius"
      !     from node (xp,yp,zp); The node numbers are stored in field
      !     ir, the distances in field r
      !
      !     x, y and z are xo, yo and zo sorted in ascending order,
      !     respectively. The order is kept in nx, ny and nz.
      !
      implicit none
      !
      integer n,nx(n),ny(n),nz(n),ir(n+6),nr,kflag,iflag,&
        i,j,id,idx,idy,idz,eight,idummy,node,l,nrprev
      !
      real*8 x(n),y(n),z(n),xo(n),yo(n),zo(n),xp,yp,zp,r(n+6),xr,yr,&
        zr,c(8),dd,xw,xe,ys,yn,zb,zt,sqrt_rmaxini,radius,radius2
      !
      intent(in) xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,n,radius
      !
      intent(inout) ir,r,nr
      !
      iflag=1
      kflag=2
      eight=8
      !
      radius2=radius*radius
      nrprev=0
      !
      !     identify position of xp, yp and zp
      !
      call ident(x,xp,n,idx)
      call ident(y,yp,n,idy)
      call ident(z,zp,n,idz)
      !
      sqrt_rmaxini=1.d30
      !
      !     initialization of the maximal distance in each direction
      !
      xw=0.d0
      xe=0.d0
      ys=0.d0
      yn=0.d0
      zb=0.d0
      zt=0.d0
      !
      i=1
      !
      do
         !
         nrprev=nr
         !
         !        west
         !
         id=idx+1-i
         if(id.gt.0) then
            node=nx(id)
            xw=xo(node)-xp
            yr=yo(node)-yp
            zr=zo(node)-zp
            dd=xw*xw+yr*yr+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.0) then
            xw=sqrt_rmaxini
         endif
         !
         !        east
         !
         id=idx+i
         if(id.le.n) then
            node=nx(id)
            xe=xo(node)-xp
            yr=yo(node)-yp
            zr=zo(node)-zp
            dd=xe*xe+yr*yr+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.n+1) then
            xe=sqrt_rmaxini
         endif
         !
         !        south
         !
         id=idy+1-i
         if(id.gt.0) then
            node=ny(id)
            xr=xo(node)-xp
            ys=yo(node)-yp
            zr=zo(node)-zp
            dd=xr*xr+ys*ys+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.0) then
            ys=sqrt_rmaxini
         endif
         !
         !        north
         !
         id=idy+i
         if(id.le.n) then
            node=ny(id)
            xr=xo(node)-xp
            yn=yo(node)-yp
            zr=zo(node)-zp
            dd=xr*xr+yn*yn+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.n+1) then
            yn=sqrt_rmaxini
         endif
         !
         !        bottom
         !
         id=idz+1-i
         if(id.gt.0) then
            node=nz(id)
            xr=xo(node)-xp
            yr=yo(node)-yp
            zb=zo(node)-zp
            dd=xr*xr+yr*yr+zb*zb
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.0) then
            zb=sqrt_rmaxini
         endif
         !
         !        top
         !
         id=idz+i
         if(id.le.n) then
            node=nz(id)
            xr=xo(node)-xp
            yr=yo(node)-yp
            zt=zo(node)-zp
            dd=xr*xr+yr*yr+zt*zt
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.n+1) then
            zt=sqrt_rmaxini
         endif
         !
         !        check the corners of the box
         !
         c(1)=xe*xe+yn*yn+zb*zb
         c(2)=xw*xw+yn*yn+zb*zb
         c(3)=xw*xw+ys*ys+zb*zb
         c(4)=xe*xe+ys*ys+zb*zb
         c(5)=xe*xe+yn*yn+zt*zt
         c(6)=xw*xw+yn*yn+zt*zt
         c(7)=xw*xw+ys*ys+zt*zt
         c(8)=xe*xe+ys*ys+zt*zt
         call dsort(c,idummy,eight,iflag)
         !
         !        check for new entries
         !
         if(nr.gt.nrprev) then
            call dsort(r,ir,nr,kflag)
            !
            !           reject equal entries
            !
            nrprev=1
            loop: do j=2,nr
               do l=nrprev,1,-1
                  if(ir(j).eq.ir(l)) cycle loop
               enddo
               nrprev=nrprev+1
               r(nrprev)=r(j)
               ir(nrprev)=ir(j)
            enddo loop
         endif
         !
         if(c(1).ge.radius2) then
            nr=nrprev
            exit
         endif
         !
         i=i+1
      !
      enddo
      !
      return
      end
