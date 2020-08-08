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
      subroutine newtincf(ithermal,dtimef,compressible,vel,
     &     hcel,umel,cvel,h,sc,iturbulent,ipkonf,nmethod,nef,lakonf,
     &     xxn,ipnei)
!
!     updates tincf 
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,j,k,ithermal(*),compressible,iturbulent,ipkonf(*),
     &  nmethod,nef,indexf,ipnei(*)
!
      real*8 dtimef,vel(nef,0:7),hcel(*),h(*),umel(*),cvel(*),sc(*),
     &  u,scmax,xn(3),hconv,hdiff,xxn(3,*),dd
!
c      if((nmethod.eq.1).and.(compressible.eq.1)) then
      if(nmethod.eq.1) then
!
!        compressible stationary flow: mass scaling
!
         dtimef=1.d30
!
!        determining the time increment needed by convection
!
         do i=1,nef
            if(ipkonf(i).lt.0) cycle
            if(lakonf(i)(1:1).ne.'F') cycle
!
!           norm of the velocity
!
            u=dsqrt(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)
!
!           calculating the convective height hconv and the diffusive
!           height hdiff
!
!           the height h(*) is calculated in initialcfd as the ratio
!           of the volume divided by the area of one of the element
!           faces. It is supposed to be orthogonal to the face. The
!           convective height is the length of a vector in the direction
!           of the velocity such that its projection on the surface    
!           normal has the value h(*)
!
            hconv=1.d30
!            
            indexf=ipnei(i)
            do j=1,ipnei(i+1)-ipnei(i)
               indexf=indexf+1
               do k=1,3
                  xn(k)=xxn(k,indexf)
               enddo
               dd=dabs(xn(1)*vel(i,1)+xn(2)*vel(i,2)+xn(3)*vel(i,3))
               if(dd.gt.1.d-30) then
                  hconv=min(hconv,h(indexf)*u/dd)
               endif
            enddo
!
!           convection
!            
            if(u.gt.1.d-30) then
               dtimef=min(dtimef,hconv/u)
            endif
         enddo
c         write(*,*) 'newtincf convection',dtimef
!
!        determining the mass scaling based on the time
!        increment needed by diffusion
!
         scmax=1.d0
         do i=1,nef
            if(ipkonf(i).lt.0) cycle
            if(lakonf(i)(1:1).ne.'F') cycle
!
!            determining the diffusive height =
!           the minimum height in the element
!            
            hdiff=1.d30
!            
            indexf=ipnei(i)
            do j=1,ipnei(i+1)-ipnei(i)
               indexf=indexf+1
               hdiff=min(hdiff,h(indexf))
            enddo
!
            hdiff=hdiff*hdiff
!            
            sc(i)=1.d0
!
!           viscous diffision
!
            if(iturbulent.eq.0) then
               if(vel(i,5).gt.0.d0) then
                  sc(i)=max(sc(i),dtimef*2.d0*umel(i)/
     &                 (vel(i,5)*hdiff))
               endif
            elseif(vel(i,5)*vel(i,7).gt.0.d0) then
               sc(i)=max(sc(i),dtimef*2.d0*
     &              (umel(i)/vel(i,5)+vel(i,6)/vel(i,7))/hdiff)
            endif
!
!           thermal diffusion
!
            if(ithermal(1).gt.0) then
               if(vel(i,5)*cvel(i)*hcel(i).gt.0.d0) then
                  sc(i)=max(sc(i),dtimef*2.d0*hcel(i)/
     &                 (vel(i,5)*cvel(i)*hdiff))
               endif
            endif
c            write(*,*) 'newtincf diffusion scale ',i,sc(i)
            if(sc(i).gt.scmax) then
               scmax=sc(i)
c               write(*,*) 'newtincf diffusion scale ',i,scmax
c               write(*,*) i,dtimef,umel(i),vel(i,5),hdiff
            endif
         enddo
         write(*,*) 'largest diffusion scaling factor ',scmax
!
c         dtimef=0.9d0*dtimef
         write(*,*) 'new value of tinc: ',dtimef
      else
!
!        incompressible flow or transient compressible flow:
!        no mass scaling         
!
         dtimef=1.d30
         do i=1,nef
            if(ipkonf(i).lt.0) cycle
            if(lakonf(i)(1:1).ne.'F') cycle
!
!           norm of the velocity
!
            u=dsqrt(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)
!
!           calculating the convective height hconv and the diffusive
!           height hdiff
!
!           the height h(*) is calculated in initialcfd as the ratio
!           of the volume divided by the area of one of the element
!           faces. It is supposed to be orthogonal to the face. The
!           convective height is the length of a vector in the direction
!           of the velocity such that its projection on the surface    
!           normal has the value h(*)
!
            hconv=1.d30
            hdiff=1.d30
!            
            indexf=ipnei(i)
            do j=1,ipnei(i+1)-ipnei(i)
               indexf=indexf+1
               do k=1,3
                  xn(k)=xxn(k,indexf)
               enddo
               dd=dabs(xn(1)*vel(i,1)+xn(2)*vel(i,2)+xn(3)*vel(i,3))
               if(dd.gt.1.d-30) then
                  hconv=min(hconv,h(indexf)*u/dd)
               endif
               hdiff=min(hdiff,h(indexf))
            enddo
!
            hdiff=hdiff*hdiff
!
!           convection
!            
            if(u.gt.1.d-30) then
               dtimef=min(dtimef,hconv/u)
            endif
!
!           viscous diffusion
!
            if(iturbulent.eq.0) then
               if(vel(i,5).gt.0.d0) then
                  dtimef=min(dtimef,vel(i,5)*hdiff/
     &                 (2.d0*umel(i)))
               endif
            elseif(vel(i,5)*vel(i,7).gt.0.d0) then
               dtimef=min(dtimef,hdiff/
     &              (2.d0*(umel(i)/vel(i,5)+vel(i,6)/vel(i,7))))
            endif
!
!           thermal diffision
!
            if(ithermal(1).gt.0) then
               if(vel(i,5)*cvel(i)*hcel(i).gt.0.d0) then
                  dtimef=min(dtimef,
     &                 vel(i,5)*cvel(i)*hdiff/(2.d0*hcel(i)))
               endif
            endif
         enddo
!
c         dtimef=0.9d0*dtimef
         write(*,*) 'new value of tinc: ',dtimef
!        
      endif
!
      return
      end
