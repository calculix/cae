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
      subroutine checkprojectgrad(nactiveold,nactive,ipoacti,ipoactiold,&
         objectset,xlambd,nnlconst,iconstacti,iconstactiold,inameacti,&
         inameactiold,g0,nobject)
      !
      !     checks the lagrange multipliers and reduces the number of active
      !     constraints if possible and update the values of the response
      !     functions for the linear constraints
      !
      implicit none
      !
      character*81 objectset(4,*)
      !
      integer nactiveold,nactive,ipoacti(*),ipoactiold(*),i,j,&
         nnlconst,nnlconstold,iconstacti(*),iconstactiold(*),&
         inameacti(*),inameactiold(*),iact,nobject
      !
      real*8 xlambd(*),g0(*)
      !
      !
      do i=1,nactive
         ipoactiold(i)=ipoacti(i)
         iconstactiold(i)=iconstacti(i)
         inameactiold(i)=inameacti(i)
      enddo
      !
      nactiveold=nactive
      nnlconstold=nnlconst
      nactive=0
      nnlconst=0
      !
      do i=1,nactiveold
         if(i.le.nnlconstold) then           
            if(iconstactiold(i).eq.-1) then
               if(xlambd(i).gt.0.d0) then
                  nactive=nactive+1
                  nnlconst=nnlconst+1
                  ipoacti(nactive)=ipoactiold(i)
                  iconstacti(nactive)=iconstactiold(i)
                  inameacti(nactive)=inameactiold(i)
               endif
            elseif(iconstactiold(i).eq.1) then
               if(xlambd(i).lt.0.d0) then
                  nactive=nactive+1
                  nnlconst=nnlconst+1
                  ipoacti(nactive)=ipoactiold(i)
                  iconstacti(nactive)=iconstactiold(i)
                  inameacti(nactive)=inameactiold(i)
               endif
            endif
         else
            if(iconstactiold(i).eq.-1) then      
               if(xlambd(i).gt.0.d0) then
                  nactive=nactive+1
                  ipoacti(nactive)=ipoactiold(i)
                  iconstacti(nactive)=iconstactiold(i)
                  inameacti(nactive)=inameactiold(i)
               endif
            elseif(iconstactiold(i).eq.1) then
               if(xlambd(i).lt.0.d0) then
                  nactive=nactive+1
                  ipoacti(nactive)=ipoactiold(i)
                  iconstacti(nactive)=iconstactiold(i)
                  inameacti(nactive)=inameactiold(i)
               endif      
            endif
         endif
      enddo
      !
      !     update the values of the response functions for the linear
      !     constraint
      !
      do i=1,nobject
         if((objectset(1,i)(1:9).eq.'THICKNESS').or.&
            (objectset(1,i)(1:12).eq.'FIXSHRINKAGE')) then
            iact=0
            do j=1,nactive
               if(inameacti(j).eq.i) then
                  iact=iact+1
               endif
            enddo
            g0(i)=iact
         endif
      enddo
      !
      return        
      end
