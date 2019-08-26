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
      subroutine shape26h(xi,et,ze,xl,xsj,shp,iflag,konl)
      !
      !     shape functions and derivatives for a 26-node quadratic
      !     isoparametric brick element. -1<=xi,et,ze<=1
      !
      !     iflag=1: calculate only the value of the shape functions
      !     iflag=2: calculate the value of the shape functions and
      !              the Jacobian determinant
      !     iflag=3: calculate the value of the shape functions, the
      !              value of their derivatives w.r.t. the global
      !              coordinates and the Jacobian determinant
      !
      implicit none
      !
      integer i,j,k,iflag,konl(26),jf(3,6),ifaceq(8,6)
      !
      real*8 shp(4,26),xs(3,3),xsi(3,3),xl(3,26),shpe(4,26),dd,&
        dd1,dd2,dd3,xi,et,ze,xsj,omg,omh,omr,opg,oph,opr,&
        tpgphpr,tmgphpr,tmgmhpr,tpgmhpr,tpgphmr,tmgphmr,tmgmhmr,tpgmhmr,&
        omgopg,omhoph,omropr,omgmopg,omhmoph,omrmopr,fxi(3),fet(3),&
        fze(3),dfxi(3),dfet(3),dfze(3)
      !
      intent(in) xi,et,ze,xl,iflag,konl
      !
      intent(out) shp,xsj
      !
      jf=reshape((/2,2,1,2,2,3,2,1,2,3,2,2,2,3,2,1,2,2/),(/3,6/))
      !
      ifaceq=reshape((/4,3,2,1,11,10,9,12,&
                        5,6,7,8,13,14,15,16,&
                        1,2,6,5,9,18,13,17,&
                        2,3,7,6,10,19,14,18,&
                        3,4,8,7,11,20,15,19,&
                        4,1,5,8,12,17,16,20/),(/8,6/))
      !
      !     shape functions in one dimension
      !
      fxi(1)=xi*(xi-1.d0)/2.d0
      fxi(2)=(1.d0-xi)*(1.d0+xi)
      fxi(3)=xi*(xi+1.d0)/2.d0
      !
      fet(1)=et*(et-1.d0)/2.d0
      fet(2)=(1.d0-et)*(1.d0+et)
      fet(3)=et*(et+1.d0)/2.d0
      !
      fze(1)=ze*(ze-1.d0)/2.d0
      fze(2)=(1.d0-ze)*(1.d0+ze)
      fze(3)=ze*(ze+1.d0)/2.d0
      !
      !     shape functions and their glocal derivatives
      !
      omg=1.d0-xi
      omh=1.d0-et
      omr=1.d0-ze
      opg=1.d0+xi
      oph=1.d0+et
      opr=1.d0+ze
      tpgphpr=opg+oph+ze
      tmgphpr=omg+oph+ze
      tmgmhpr=omg+omh+ze
      tpgmhpr=opg+omh+ze
      tpgphmr=opg+oph-ze
      tmgphmr=omg+oph-ze
      tmgmhmr=omg+omh-ze
      tpgmhmr=opg+omh-ze
      omgopg=omg*opg/4.d0
      omhoph=omh*oph/4.d0
      omropr=omr*opr/4.d0
      omgmopg=(omg-opg)/4.d0
      omhmoph=(omh-oph)/4.d0
      omrmopr=(omr-opr)/4.d0
      !
      !     shape functions
      !
      shp(4, 1)=-omg*omh*omr*tpgphpr/8.d0
      shp(4, 2)=-opg*omh*omr*tmgphpr/8.d0
      shp(4, 3)=-opg*oph*omr*tmgmhpr/8.d0
      shp(4, 4)=-omg*oph*omr*tpgmhpr/8.d0
      shp(4, 5)=-omg*omh*opr*tpgphmr/8.d0
      shp(4, 6)=-opg*omh*opr*tmgphmr/8.d0
      shp(4, 7)=-opg*oph*opr*tmgmhmr/8.d0
      shp(4, 8)=-omg*oph*opr*tpgmhmr/8.d0
      shp(4, 9)=omgopg*omh*omr
      shp(4,10)=omhoph*opg*omr
      shp(4,11)=omgopg*oph*omr
      shp(4,12)=omhoph*omg*omr
      shp(4,13)=omgopg*omh*opr
      shp(4,14)=omhoph*opg*opr
      shp(4,15)=omgopg*oph*opr
      shp(4,16)=omhoph*omg*opr
      shp(4,17)=omropr*omg*omh
      shp(4,18)=omropr*opg*omh
      shp(4,19)=omropr*opg*oph
      shp(4,20)=omropr*omg*oph
      !
      !     correction for the extra nodes in the middle of the faces
      !
      do i=1,6
         if(konl(20+i).eq.konl(20)) then
            !
            !           no extra node in this face
            !
            shp(4,20+i)=0.d0
         else
            shp(4,20+i)=fxi(jf(1,i))*fet(jf(2,i))*fze(jf(3,i))
            do j=1,4
               shp(4,ifaceq(j,i))=shp(4,ifaceq(j,i))+shp(4,20+i)/4.d0
            enddo
            do j=5,8
               shp(4,ifaceq(j,i))=shp(4,ifaceq(j,i))-shp(4,20+i)/2.d0
            enddo
         endif
      enddo
      !
      if(iflag.eq.1) return
      !
      !     derivative of the shape functions in one dimension
      !
      dfxi(1)=(2.d0*xi-1.d0)/2.d0
      dfxi(2)=-2.d0*xi
      dfxi(3)=(2.d0*xi+1.d0)/2.d0
      !
      dfet(1)=(2.d0*et-1.d0)/2.d0
      dfet(2)=-2.d0*et
      dfet(3)=(2.d0*et+1.d0)/2.d0
      !
      dfze(1)=(2.d0*ze-1.d0)/2.d0
      dfze(2)=-2.d0*ze
      dfze(3)=(2.d0*ze+1.d0)/2.d0
      !
      !     local derivatives of the shape functions: xi-derivative
      !
      shpe(1, 1)=omh*omr*(tpgphpr-omg)/8.d0
      shpe(1, 2)=(opg-tmgphpr)*omh*omr/8.d0
      shpe(1, 3)=(opg-tmgmhpr)*oph*omr/8.d0
      shpe(1, 4)=oph*omr*(tpgmhpr-omg)/8.d0
      shpe(1, 5)=omh*opr*(tpgphmr-omg)/8.d0
      shpe(1, 6)=(opg-tmgphmr)*omh*opr/8.d0
      shpe(1, 7)=(opg-tmgmhmr)*oph*opr/8.d0
      shpe(1, 8)=oph*opr*(tpgmhmr-omg)/8.d0
      shpe(1, 9)=omgmopg*omh*omr
      shpe(1,10)=omhoph*omr
      shpe(1,11)=omgmopg*oph*omr
      shpe(1,12)=-omhoph*omr
      shpe(1,13)=omgmopg*omh*opr
      shpe(1,14)=omhoph*opr
      shpe(1,15)=omgmopg*oph*opr
      shpe(1,16)=-omhoph*opr
      shpe(1,17)=-omropr*omh
      shpe(1,18)=omropr*omh
      shpe(1,19)=omropr*oph
      shpe(1,20)=-omropr*oph
      !
      !     correction for the extra nodes in the middle of the faces
      !
      do i=1,6
         if(konl(20+i).eq.konl(20)) then
            !
            !           no extra node in this face
            !
            shpe(1,20+i)=0.d0
         else
            shpe(1,20+i)=dfxi(jf(1,i))*fet(jf(2,i))*fze(jf(3,i))
            do j=1,4
               shpe(1,ifaceq(j,i))=shpe(1,ifaceq(j,i))+shpe(1,20+i)/4.d0
            enddo
            do j=5,8
               shpe(1,ifaceq(j,i))=shpe(1,ifaceq(j,i))-shpe(1,20+i)/2.d0
            enddo
         endif
      enddo
      !
      !     local derivatives of the shape functions: eta-derivative
      !
      shpe(2, 1)=omg*omr*(tpgphpr-omh)/8.d0
      shpe(2, 2)=opg*omr*(tmgphpr-omh)/8.d0
      shpe(2, 3)=opg*(oph-tmgmhpr)*omr/8.d0
      shpe(2, 4)=omg*(oph-tpgmhpr)*omr/8.d0
      shpe(2, 5)=omg*opr*(tpgphmr-omh)/8.d0
      shpe(2, 6)=opg*opr*(tmgphmr-omh)/8.d0
      shpe(2, 7)=opg*(oph-tmgmhmr)*opr/8.d0
      shpe(2, 8)=omg*(oph-tpgmhmr)*opr/8.d0
      shpe(2, 9)=-omgopg*omr
      shpe(2,10)=omhmoph*opg*omr
      shpe(2,11)=omgopg*omr
      shpe(2,12)=omhmoph*omg*omr
      shpe(2,13)=-omgopg*opr
      shpe(2,14)=omhmoph*opg*opr
      shpe(2,15)=omgopg*opr
      shpe(2,16)=omhmoph*omg*opr
      shpe(2,17)=-omropr*omg
      shpe(2,18)=-omropr*opg
      shpe(2,19)=omropr*opg
      shpe(2,20)=omropr*omg
      !
      !     correction for the extra nodes in the middle of the faces
      !
      do i=1,6
         if(konl(20+i).eq.konl(20)) then
            !
            !           no extra node in this face
            !
            shpe(2,20+i)=0.d0
         else
            shpe(2,20+i)=fxi(jf(1,i))*dfet(jf(2,i))*fze(jf(3,i))
            do j=1,4
               shpe(2,ifaceq(j,i))=shpe(2,ifaceq(j,i))+shpe(2,20+i)/4.d0
            enddo
            do j=5,8
               shpe(2,ifaceq(j,i))=shpe(2,ifaceq(j,i))-shpe(2,20+i)/2.d0
            enddo
         endif
      enddo
      !
      !     local derivatives of the shape functions: zeta-derivative
      !
      shpe(3, 1)=omg*omh*(tpgphpr-omr)/8.d0
      shpe(3, 2)=opg*omh*(tmgphpr-omr)/8.d0
      shpe(3, 3)=opg*oph*(tmgmhpr-omr)/8.d0
      shpe(3, 4)=omg*oph*(tpgmhpr-omr)/8.d0
      shpe(3, 5)=omg*omh*(opr-tpgphmr)/8.d0
      shpe(3, 6)=opg*omh*(opr-tmgphmr)/8.d0
      shpe(3, 7)=opg*oph*(opr-tmgmhmr)/8.d0
      shpe(3, 8)=omg*oph*(opr-tpgmhmr)/8.d0
      shpe(3, 9)=-omgopg*omh
      shpe(3,10)=-omhoph*opg
      shpe(3,11)=-omgopg*oph
      shpe(3,12)=-omhoph*omg
      shpe(3,13)=omgopg*omh
      shpe(3,14)=omhoph*opg
      shpe(3,15)=omgopg*oph
      shpe(3,16)=omhoph*omg
      shpe(3,17)=omrmopr*omg*omh
      shpe(3,18)=omrmopr*opg*omh
      shpe(3,19)=omrmopr*opg*oph
      shpe(3,20)=omrmopr*omg*oph
      !
      !     correction for the extra nodes in the middle of the faces
      !
      do i=1,6
         if(konl(20+i).eq.konl(20)) then
            !
            !           no extra node in this face
            !
            shpe(3,20+i)=0.d0
         else
            shpe(3,20+i)=fxi(jf(1,i))*fet(jf(2,i))*dfze(jf(3,i))
            do j=1,4
               shpe(3,ifaceq(j,i))=shpe(3,ifaceq(j,i))+shpe(3,20+i)/4.d0
            enddo
            do j=5,8
               shpe(3,ifaceq(j,i))=shpe(3,ifaceq(j,i))-shpe(3,20+i)/2.d0
            enddo
         endif
      enddo
      !
      !     computation of the local derivative of the global coordinates
      !     (xs)
      !
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,26
            xs(i,j)=xs(i,j)+xl(i,k)*shpe(j,k)
          enddo
        enddo
      enddo
      !
      !     computation of the jacobian determinant
      !
      dd1=xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2)
      dd2=xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3)
      dd3=xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1)
      xsj=xs(1,1)*dd1+xs(1,2)*dd2+xs(1,3)*dd3
      !
      if(iflag.eq.2) return
      !
      dd=1.d0/xsj
      !
      !     computation of the global derivative of the local coordinates
      !     (xsi) (inversion of xs)
      !
      xsi(1,1)=dd1*dd
      xsi(1,2)=(xs(1,3)*xs(3,2)-xs(1,2)*xs(3,3))*dd
      xsi(1,3)=(xs(1,2)*xs(2,3)-xs(2,2)*xs(1,3))*dd
      xsi(2,1)=dd2*dd
      xsi(2,2)=(xs(1,1)*xs(3,3)-xs(3,1)*xs(1,3))*dd
      xsi(2,3)=(xs(1,3)*xs(2,1)-xs(1,1)*xs(2,3))*dd
      xsi(3,1)=dd3*dd
      xsi(3,2)=(xs(1,2)*xs(3,1)-xs(1,1)*xs(3,2))*dd
      xsi(3,3)=(xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2))*dd
      !
      !     computation of the global derivatives of the shape functions
      !
      do k=1,26
        do j=1,3
          shp(j,k)=shpe(1,k)*xsi(1,j)+shpe(2,k)*xsi(2,j)&
                +shpe(3,k)*xsi(3,j)
        enddo
      enddo
      !
      return
      end
