!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine calcstressheatfluxfem(kon,lakon,ipkon,
     &  ielmat,ntmat_,vold,matname,mi,shcon,nshcon,
     &  turbulent,compressible,ipvar,var,sti,qfx,cocon,
     &  ncocon,ne,isti,iqfx)
!
!     calculating the stresses and heat flow in the integration points
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*),amat
!
      integer kon(*),nelem,mi(*),ne,konl(20),ipkon(*),j,i1,i2,j1,ii,
     &  jj,
     &  indexe,isti,iqfx,kk,ielmat(mi(3),*),nshcon(*),ntmat_,nope,
     &  imat,ncocon(2,*),mint3d,k1,ipvar(*),index,turbulent,compressible
!
      real*8 shp(4,20),dvi,cond,cocon(0:6,ntmat_,*),div,
     &  shcon(0:3,ntmat_,*),voldl(0:mi(2),20),vold(0:mi(2),*),temp,
     &  tt(3,3),rho,t(3,3),vkl(3,3),xkin,unt,umt,
     &  var(*),sti(6,mi(1),*),dtem(3),qfx(3,mi(1),*)
!
      do nelem=1,ne
!
         if(ipkon(nelem).lt.0) cycle
         if(lakon(nelem)(1:1).ne.'F') cycle
         indexe=ipkon(nelem)
!     
         imat=ielmat(1,nelem)
         amat=matname(imat)
!     
         if(lakon(nelem)(4:4).eq.'4') then
            nope=4
            mint3d=1
         elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
            mint3d=2
         elseif(lakon(nelem)(4:5).eq.'8R') then
            nope=8
            mint3d=1
         elseif(lakon(nelem)(4:4).eq.'8') then
            nope=8
            mint3d=8
         elseif(lakon(nelem)(4:5).eq.'10') then
            nope=10
            mint3d=4
         elseif(lakon(nelem)(4:5).eq.'15') then
            nope=15
            mint3d=9
         elseif(lakon(nelem)(4:6).eq.'20R') then
            nope=20
            mint3d=8
         elseif(lakon(nelem)(4:4).eq.'2') then
            nope=20
            mint3d=27
         else
            nope=0
            mint3d=0
         endif
!     
         do j=1,nope
            konl(j)=kon(indexe+j) 
         enddo
!     
!     storing the local temperature and velocity
!     
         do i1=1,nope
            do i2=0,3
               voldl(i2,i1)=vold(i2,konl(i1))
            enddo
         enddo
!     
!     computation of the matrix: loop over the Gauss points
!     
         index=ipvar(nelem)
c                  write(*,*) 'index ',index
         do kk=1,mint3d
!     
!     copying the shape functions and their derivatives from field var
!     
            do jj=1,nope
               do ii=1,4
                  index=index+1
                  shp(ii,jj)=var(index)
c                  write(*,*) 'shp ',jj,ii,shp(ii,jj)
               enddo
            enddo
!     
!     calculating of
!     the velocity gradient vkl
!     the temperature gradient dtem
!   
            if(isti.gt.0) then
               do i1=1,3
                  do j1=1,3
                     vkl(i1,j1)=0.d0
                  enddo
               enddo
               do i1=1,nope
                  do j1=1,3
                     do k1=1,3
                        vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                     enddo
                  enddo
               enddo
               if(compressible.eq.1) div=vkl(1,1)+vkl(2,2)+vkl(3,3)
            endif
c               write(*,*) 'calcstressheatflux ',kk,nelem,
c     &             ((vkl(i1,j1),j1=1,3),i1=1,3)
!     
            if(iqfx.gt.0) then
               do i1=1,3
                  dtem(i1)=0.d0
               enddo
               do i1=1,nope
                  do j1=1,3
                     dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
                  enddo
               enddo
            endif
!     
!     storing shpv, vel and temp
!     
            index=index+nope+2
            temp=var(index)
c               write(*,*) 'calcstressheatflux ',kk,nelem,
c     &             temp
!     
!     determining the dissipative stress 
!     
            if(isti.gt.0) then
              call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,
     &             dvi)
               do i1=1,3
                  do j1=i1,3
                     t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
                  enddo
                  if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
               enddo
!     
!     calculating the stress
!     
               if(turbulent.ne.0) then
!     
                  rho=var(index+8)
                  xkin=var(index+10)
                  unt=var(index+12)
!     
                  umt=unt*rho
!     
!     calculating the turbulent stress
!     
                  do i1=1,3
                     do j1=i1,3
                        tt(i1,j1)=umt*t(i1,j1)
                     enddo
                     tt(i1,i1)=tt(i1,i1)-2.d0*rho*xkin/3.d0
                  enddo
!     
!     adding the viscous stress
!     
                  do i1=1,3
                     do j1=i1,3
                        t(i1,j1)=dvi*t(i1,j1)+tt(i1,j1)
                     enddo
                  enddo
               else
!     
                  do i1=1,3
                     do j1=i1,3
                        t(i1,j1)=dvi*t(i1,j1)
                     enddo
                  enddo
               endif
!     
!     storing the total stress
!     
               sti(1,kk,nelem)=t(1,1)
               sti(2,kk,nelem)=t(2,2)
               sti(3,kk,nelem)=t(3,3)
               sti(4,kk,nelem)=t(1,2)
               sti(5,kk,nelem)=t(1,3)
               sti(6,kk,nelem)=t(2,3)
c               write(*,*) 'calcstressheatflux ',kk,nelem,
c     &             (sti(i1,kk,nelem),i1=1,6)
            endif
!     
!     storing the heat flow
!     
            if(iqfx.gt.0) then
               call materialdata_cond(imat,ntmat_,temp,cocon,
     &               ncocon,cond)
               qfx(1,kk,nelem)=-cond*dtem(1)
               qfx(2,kk,nelem)=-cond*dtem(2)
               qfx(3,kk,nelem)=-cond*dtem(3)
            endif
!
            index=index+13
!     
         enddo
      enddo
!     
      return
      end

