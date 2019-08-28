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
      subroutine calcmatwavspeed(ne0,elcon,nelcon,&
                 rhcon,nrhcon,alcon,nalcon,orab,ntmat_,ithermal,alzero,&
                 plicon,nplicon,plkcon,nplkcon,npmat_,mi,dtime,&
                 xstiff,ncmat_,vold,ielmat,t0,t1,&
                 matname,lakon,wavespeed,nmat,ipkon)
 !
 !     **************
 !     Calculates the propagation wave speed in a material, selecting
 !     appropiate procedure between isotropic, single crystals, and
 !     anisotropic materials. All other cases of orthotropy are treated
 !     as anisotropic.

      !     Based on the procedure in:
      !     C. Lane. The Development of a 2D Ultrasonic Array Inspection
      !     for Single Crystal Turbine Blades.
      !     Switzerland: Springer International Publishing, 2014.
      !
      implicit none
      !
      character*80 matname(*),amat
      character*8 lakon(*),lakonl
      !
      integer i,i1,nelem,ne0,nelcon(2,*),nrhcon(*),nalcon(2,*),imat,&
           ntmat_,ithermal,mattyp,ihyper,istiff,kode,mi(*),kk,&
           nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,iorth,&
           ielmat(mi(3),*),nope,iorien,ipkon(*),null,&
           konl(26),nopered,npmat_,nmat
      !
      real*8 elas(21),wavespeed(*),rhcon(0:1,ntmat_,*)  ,&
           alcon(0:6,ntmat_,*),coords(3),orab(7,*),rho,alzero(*),&
           t0l,t1l,elconloc(21),eth(6),plicon(0:2*npmat_,ntmat_,*),&
           plkcon(0:2*npmat_,ntmat_,*),plconloc(802),dtime,&
           xstiff(27,mi(1),*),elcon(0:ncmat_,ntmat_,*),&
           t0(*),t1(*),shp(4,26),vold(0:mi(2),*),tt,&
           e,un,um,al,wavspd
      !
      write(*,*)'++CMT: Calculating Material Wave Speeds...'
      write(*,*)
      !
      null=0
      !
      !     initialization of wavespeed
      !
      do i=1,nmat
         wavespeed(i)=-1.d0
      enddo
      !
      !     ** DO per element
      !
      do nelem=1,ne0 
         if(ipkon(nelem).lt.0) cycle
         !
         lakonl=lakon(nelem)
         imat=ielmat(1,nelem)
         amat=matname(imat)
         wavspd=wavespeed(imat)
         !
         if(lakonl(1:2).ne.'ES') then
            !
            !     orientation is not important for the wave speed
            !
            iorien=0
            !
            if(nelcon(1,imat).lt.0) then
               ihyper=1
            else
               ihyper=0
            endif
            !
            !     calculating the temperature in the integration
            !     point
            !
            kk=1                ! Only 1 integration point is considered, CENTER
            t0l=0.d0
            t1l=0.d0
            if(ithermal.eq.1) then
               if(lakonl(4:5).eq.'8 ') then
                  do i1=1,nope
                     t0l=t0l+t0(konl(i1))/8.d0
                     t1l=t1l+t1(konl(i1))/8.d0
                  enddo
               elseif(lakonl(4:6).eq.'20 ')then
                  nopered=20
                  call lintemp(t0,t1,konl,nopered,kk,t0l,t1l)
               elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(t0,konl,t0l,null,shp)
                  call linscal10(t1,konl,t1l,null,shp)
               else
                  do i1=1,nope
                     t0l=t0l+shp(4,i1)*t0(konl(i1))
                     t1l=t1l+shp(4,i1)*t1(konl(i1))
                  enddo
               endif
            elseif(ithermal.ge.2) then
               if(lakonl(4:5).eq.'8 ') then
                  do i1=1,nope
                     t0l=t0l+t0(konl(i1))/8.d0
                     t1l=t1l+vold(0,konl(i1))/8.d0
                  enddo
               elseif(lakonl(4:6).eq.'20 ')then
                  nopered=20
                  call lintemp_th(t0,vold,konl,nopered,kk,t0l,t1l,mi)
               elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(t0,konl,t0l,null,shp)
                  call linscal10(vold,konl,t1l,mi(2),shp)
               else
                  do i1=1,nope
                     t0l=t0l+shp(4,i1)*t0(konl(i1))
                     t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                  enddo
               endif
            endif
            tt=t1l-t0l
            !
            kode=nelcon(1,imat)         
            
            !     material data and local stiffness matrix
            !
            istiff=1
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,&
                 imat,amat,iorien,coords,orab,ntmat_,elas,rho,&
                 nelem,ithermal,alzero,mattyp,t0l,t1l,&
                 ihyper,istiff,elconloc,eth,kode,plicon,&
                 nplicon,plkcon,nplkcon,npmat_,&
                 plconloc,mi(1),dtime,kk,&
                 xstiff,ncmat_)
            !
            if(mattyp.eq.1) then
               e=elas(1)
               un=elas(2)
               um=e/(1.d0+un)
               al=un*um/(1.d0-2.d0*un)
               um=um/2.d0
               wavspd=max(wavspd,&
                   dsqrt((e*(1-un))/((1+un)*(1-2*un)*rho)))
            elseif(mattyp.eq.2) then
               !
               !     single crystal
               !
               if(((elas(1).eq.elas(3)).and.(elas(1).eq.elas(6)).and.&
                    (elas(3).eq.elas(6))).and.&
                    ((elas(2).eq.elas(4)).and.(elas(2).eq.elas(5)).and.&
                    (elas(4).eq.elas(5))).and.&
                    ((elas(7).eq.elas(8)).and.(elas(7).eq.elas(9)).and.&
                    (elas(8).eq.elas(9)))) then
                  wavspd=max(wavspd,dsqrt((1/3.d0)*(elas(1)+2.0*elas(2)+&
                       4.0d0*elas(7))/rho))
                  wavspd=max(wavspd,dsqrt((1/2.d0)*&
                       (elas(1)+elas(2)+2.0*elas(7))/rho))
                  wavspd=max(wavspd,dsqrt(elas(1)/rho))
               else
                  iorth=1
                  call anisomaxwavspd(elas,rho,iorth,wavspd )
               endif
            elseif(mattyp.eq.3) then
               iorth=0
               call anisomaxwavspd(elas,rho,iorth,wavspd)
            endif
            wavespeed(imat)=wavspd
         endif
      enddo     
      !
      do i=1,nmat
         write(*,*) 'Wave Speed for mat. ',i,wavespeed(i)   
      enddo
      write(*,*)
      !
      return
      end
