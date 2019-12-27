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
      subroutine actideacti(set,nset,istartset,iendset,ialset,&
                 objectset,ipkon,iobject,ne)
      !
      !
      implicit none
      !
      character*81 objectset(4,*),set(*)
      !
      integer i,j,k,nset,istartset(*),iendset(*),ialset(*),ipkon(*),&
        iobject,ne
      !
      intent(in) set,nset,istartset,iendset,ialset,&
                 objectset,iobject,ne
      !
      intent(inout) ipkon
      !
      !     determining the set
      !
      do i=1,nset
         if(objectset(3,iobject).eq.set(i)) exit
      enddo
      !
      if(i.le.nset) then
         !
         !        deactivate all elements
         !
         do j=1,ne
            if(ipkon(j).lt.0) cycle
            ipkon(j)=-2-ipkon(j)
         enddo
         !
         !        reactivate the elements belonging to the set
         !
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               ipkon(ialset(j))=-ipkon(ialset(j))-2
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  ipkon(k)=-ipkon(k)-2
               enddo
            endif
         enddo
      endif
      !
      return
      end
