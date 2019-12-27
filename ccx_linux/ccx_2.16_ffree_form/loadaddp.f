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
      subroutine loadaddp(nelement,label,nelemload,sideload,&
        xload,nload,nload_,iamload,iamplitude,nam,node)
      !
      !     adds a facial dload condition to the data base
      !
      implicit none
      !
      character*20 label,sideload(*)
      !
      integer nelemload(2,*),iamload(2,*),nelement,nload,nload_,j,&
          iamplitude,nam,node,id
      !
      real*8 xload(2,*)
      !
      call nident2(nelemload,nelement,nload,id)
      if(id.gt.0) then
         !
         !        it is possible that several *DLOAD, *FILM or
         !        *RADIATE boundary conditions are applied to one
         !        and the same element
         !
         if(nelemload(1,id).eq.nelement) then
            do
               if(sideload(id).eq.label) then
                  !
                  !     loading on same element face and sector
                  !     detected: values are replaced
                  !
                  xload(1,id)=0.d0
                  xload(2,id)=0.d0
                  if(nam.gt.0) then
                     iamload(1,id)=iamplitude
                     iamload(2,id)=0.d0
                  endif
                  return
               endif
               id=id-1
               if((id.eq.0).or.(nelemload(1,id).ne.nelement)) then
                  exit
               endif
            enddo
         endif
      endif
      !
      !     loading a element face on which no previous loading
      !     was applied
      !
      !     loading conditions on one and the same element are
      !     alphabetized based on field sideload
      !
      nload=nload+1
      if(nload.gt.nload_) then
         write(*,*) '*ERROR in loadadd: increase nload_'
         call exit(201)
      endif
      !
      !     shifting existing loading
      !
      do j=nload,id+2,-1
         nelemload(1,j)=nelemload(1,j-1)
         nelemload(2,j)=nelemload(2,j-1)
         sideload(j)=sideload(j-1)
         xload(1,j)=xload(1,j-1)
         xload(2,j)=xload(2,j-1)
         if(nam.gt.0) then
            iamload(1,j)=iamload(1,j-1)
            iamload(2,j)=iamload(2,j-1)
         endif
      enddo
      !
      !     inserting new loading
      !
      nelemload(1,id+1)=nelement
      nelemload(2,id+1)=node
      sideload(id+1)=label
      xload(1,id+1)=0.d0
      xload(2,id+1)=0.d0
      if(nam.gt.0) then
         iamload(1,id+1)=iamplitude
         iamload(2,id+1)=0
      endif
      !
      return
      end

