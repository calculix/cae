!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine calcstabletimeincvol(ne0,lakon,co,kon,ipkon,mi,&
        ielmat,dtvol,alpha,wavespeed)
      !
      !     Calculates the critical time increment (CTI) based on the Courant
      !     Criterion for Explicit Dynamics calculations. Temperature is
      !     assumed averaged from the centroid of the element and material
      !     wave propagation speeds must be calculated before.
      !
      implicit none
      !
      character*8 lakon(*),lakonl
      !
      integer i,j,ne0,nope,kon(*),ipkon(*),indexe,konl(26),nelem,&
        iflag,nopes,nfaces,ig,ifaceq(8,6),ifacet(6,4),ifacew(8,5),&
        mi(*),ielmat(mi(3),*),imat,elemmin
      !
      real*8 xi,et,ze,weight,co(3,*),xl(3,26),xsj,shp(4,26),xl2(3,9),&
        xsj2(3),xs2(3,7),shp2(7,9),hmin,area,volume,&
        wavspd,dtvol,safefac,alpha,bet,gam,critom,damping,&
        wavespeed(*),geomfac,quadfac
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
           5,6,7,8,13,14,15,16,&
           1,2,6,5,9,18,13,17,&
           2,3,7,6,10,19,14,18,&
           3,4,8,7,11,20,15,19,&
           4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,&
           1,2,4,5,9,8,&
           2,3,4,6,10,9,&
           1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,&
           4,5,6,10,11,12,0,0,&
           1,2,5,4,7,14,10,13,&
           2,3,6,5,8,15,11,14,&
           4,6,3,1,12,15,9,13/
      !
      include "gauss.f"
      !
      iflag=2
      dtvol=1.d30
      safefac=0.80d0
      quadfac=0.3d0
      !
      damping=0 
      !
      bet=(1.d0-alpha)*(1.d0-alpha)/4.d0
      gam=0.5d0-alpha
      !
      !     Omega Critical
      !     Om_cr=dt*freq_max
      !
      critom=dsqrt(damping*damping*(1.d0+2.d0*alpha*(1.d0-gam))&
           *(1.d0+2.d0*alpha*(1.d0-gam))&
           +2.d0*(gam+2.d0*alpha*(gam-bet)) )
      critom=0.98d0*(-damping*(1.d0+2.d0*alpha*(1.d0-gam))+critom)&
           /(gam+2.d0*alpha*(gam-bet)) !eq 25 miranda
      !
      !     ** DO per element
      do nelem=1,ne0
         if(ipkon(nelem).lt.0) cycle
         indexe=ipkon(nelem)
         !
         lakonl=lakon(nelem)
         imat=ielmat(1,nelem)
         !
         geomfac=1.d0
         !
         if(lakon(nelem)(4:5).eq.'20')then 
            nope=20     
            nopes=8          
            nfaces=6
            geomfac=quadfac
         elseif(lakon(nelem)(1:5).eq.'C3D8I')then
            nope=8 
            nopes=4
            nfaces=6
            geomfac=quadfac
            geomfac=0.5d0
         elseif(lakon(nelem)(4:4).eq.'8') then
            nope=8
            nopes=4
            nfaces=6
         elseif(lakon(nelem)(4:5).eq.'10')then
            nope=10
            nopes=6
            nfaces=4
            geomfac=quadfac
         elseif(lakon(nelem)(4:4).eq.'4') then
            nope=4
            nopes=3
            nfaces=4
         elseif(lakon(nelem)(4:5).eq.'15')then
            nope=15
            nfaces=5 
            geomfac=quadfac 
         elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
            nfaces=5 
         else           
            cycle   
         endif
         !
         !     Find center of the element for avg temp value on the element  to
         !     get properties later
         !     if HEX
         if((lakon(nelem)(4:5).eq.'20').or.&
              (lakon(nelem)(4:4).eq.'8')) then
            xi=0.d0
            et=0.d0
            ze=0.d0
            weight=8.d0
         !     if TET
         elseif((lakon(nelem)(4:5).eq.'10').or.&
                 (lakon(nelem)(4:4).eq.'4')) then
            xi=gauss3d4(1,1)
            et=gauss3d4(2,1)
            ze=gauss3d4(3,1)
            weight=weight3d4(1)
         !     elseif WEDGES
         elseif((lakonl(4:5).eq.'15').or.&
                (lakonl(4:4).eq.'6'))then
            xi=1.d0/3.d0
            et=1.d0/3.d0
            ze=0.d0
            weight=1.d0
         endif   
         !
         do i=1,nope
            konl(i)=kon(indexe+i)
            do j=1,3
               xl(j,i)=co(j,konl(i))
            enddo
         enddo
         !
         if   (nope.eq.20)then
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.10)then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.4) then
            call shape4tet (xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.15)then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
         else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
         endif
         !
         wavspd=wavespeed(imat)
         !
         !     Divides volume accordingly per geometry of element
         !     Carlo MT proposal
         !     if HEX
         if((lakon(nelem)(4:5).eq.'20').or.&
              (lakon(nelem)(4:4).eq.'8')) then
            volume=weight*xsj
         !     if TET
         elseif((lakon(nelem)(4:5).eq.'10').or.&
                 (lakon(nelem)(4:4).eq.'4')) then
            volume=weight*xsj/3.d0
         !     if WEDGES
         elseif ( (lakonl(4:5).eq.'15').or.&
                 (lakonl(4:4).eq.'6'))then
            volume=weight*xsj/2.d0
         endif
         !
         hmin=1.d30
         !
         !     DO over sides
         do ig=1,nfaces
            if(lakon(nelem)(4:4).eq.'6')then
               if(ig.le.2)then
                  nopes=3
               else
                  nopes=4
               endif
            elseif(lakon(nelem)(4:5).eq.'15')then
               if(ig.le.2)then
                  nopes=6
               else
                  nopes=8
               endif
            endif
            !
            if((nope.eq.20).or.(nope.eq.8))then
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4))then
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifacet(i,ig)))
                  enddo
               enddo
            else
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifacew(i,ig)))
                  enddo
               enddo
            endif
            !
            if((nopes.eq.4).or.(nopes.eq.8))then
               xi=0.d0
               et=0.d0
               weight=4.d0
            else
               xi=1.d0/3.d0 
               et=1.d0/3.d0 
               weight=0.5d0
            endif
            !
            if    (nopes.eq.8) then
               call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            !             elseif(nopes.eq.7) then
            !                call shape7tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
               call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
            !
            area=weight*dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+&
                 xsj2(3)*xsj2(3))
            hmin=min(hmin,(volume/area))
         !
         enddo
         !     ENDDO over sides
         !
         if(critom/2*hmin/wavspd*geomfac.lt.dtvol)then
            elemmin=nelem
         endif

         dtvol=min(dtvol,critom/2* hmin/wavspd*geomfac) 
      enddo
      !     ** ENDDO per element
      !
      dtvol=dtvol*safefac
      !
      return
      end
