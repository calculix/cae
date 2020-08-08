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
      subroutine cataloguenodes(ipofano,ifano,ifreefa,ielfa,ifabou,
     &  ipkon,konf,lakon,nface,nk)
!
!     catalogues the nodes for interpolation
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipofano(*),ifano(2,*),ifreefa,ielfa(4,*),indexb,nk,
     &  ifabou(*),iel,indexe,ipkon(*),k,node,ifaceq(8,6),nface,i,j,
     &  konf(*),ifacet(7,4),ifacew(8,5),nf(5)
!
!     nodes belonging to the cell faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data nf /3,3,4,4,4/
!
      ifreefa=0
!
!     taking into account the walls (first priority)
!
      do i=1,nface
         indexb=ielfa(2,i)
         if(indexb.ge.0) cycle
c         if(ifabou(-indexb+5).eq.1) then
         if(ifabou(-indexb+5).gt.0) then
            iel=ielfa(1,i)
            indexe=ipkon(iel)
            j=ielfa(4,i)
!
            if(lakon(iel)(4:4).eq.'8') then
               do k=1,4
                  node=konf(indexe+ifaceq(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            elseif(lakon(iel)(4:4).eq.'6') then
               do k=1,nf(j)
                  node=konf(indexe+ifacew(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            else
               do k=1,3
                  node=konf(indexe+ifacet(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            endif
         endif
      enddo
!
!     marking the nodes used so far by putting a minus sign in front
!
      do i=1,nk
         if(ipofano(i).gt.0) ipofano(i)=-ipofano(i)
      enddo
!
!     taking into account the nodes at which all velocity
!     components are prescribed (second priority)
!
      do i=1,nface
         indexb=ielfa(2,i)
         if(indexb.ge.0) cycle
         if((ifabou(-indexb+1).ne.0).and.
     &      (ifabou(-indexb+2).ne.0).and.
     &      (ifabou(-indexb+3).ne.0)) then
            iel=ielfa(1,i)
            indexe=ipkon(iel)
            j=ielfa(4,i)
!
            if(lakon(iel)(4:4).eq.'8') then
               do k=1,4
                  node=konf(indexe+ifaceq(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            elseif(lakon(iel)(4:4).eq.'6') then
               do k=1,nf(j)
                  node=konf(indexe+ifacew(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            else
               do k=1,3
                  node=konf(indexe+ifacet(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            endif
         endif
      enddo
!
!     marking the nodes used so far by putting a minus sign in front
!
      do i=1,nk
         if(ipofano(i).gt.0) ipofano(i)=-ipofano(i)
      enddo
!
!     taking into account the nodes in which the pressure is
!     prescribed (third priority)
!
      do i=1,nface
         indexb=ielfa(2,i)
         if(indexb.ge.0) cycle
         if(ifabou(-indexb+4).ne.0) then
            iel=ielfa(1,i)
            indexe=ipkon(iel)
            j=ielfa(4,i)
!
            if(lakon(iel)(4:4).eq.'8') then
               do k=1,4
                  node=konf(indexe+ifaceq(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            elseif(lakon(iel)(4:4).eq.'6') then
               do k=1,nf(j)
                  node=konf(indexe+ifacew(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            else
               do k=1,3
                  node=konf(indexe+ifacet(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            endif
         endif
      enddo
!
!     marking the nodes used so far by putting a minus sign in front
!
      do i=1,nk
         if(ipofano(i).gt.0) ipofano(i)=-ipofano(i)
      enddo
!
!     taking into account faces with given mass flow (fourth priority)
!
      do i=1,nface
         indexb=ielfa(2,i)
         if(indexb.ge.0) cycle
c         if(ifabou(-indexb+5).eq.2) then
         if(ifabou(-indexb+5).lt.0) then
            iel=ielfa(1,i)
            indexe=ipkon(iel)
            j=ielfa(4,i)
!
            if(lakon(iel)(4:4).eq.'8') then
               do k=1,4
                  node=konf(indexe+ifaceq(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            elseif(lakon(iel)(4:4).eq.'6') then
               do k=1,nf(j)
                  node=konf(indexe+ifacew(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            else
               do k=1,3
                  node=konf(indexe+ifacet(k,j))
                  if(ipofano(node).eq.0) then
                     ifreefa=ifreefa+1
                     ipofano(node)=ifreefa
                     ifano(1,ifreefa)=i
                  elseif(ipofano(node).gt.0) then
                     ifreefa=ifreefa+1
                     ifano(2,ifreefa)=ipofano(node)
                     ifano(1,ifreefa)=i
                     ipofano(node)=ifreefa
                  endif
               enddo
            endif
         endif
      enddo
!
!     marking the nodes used so far by putting a minus sign in front
!
      do i=1,nk
         if(ipofano(i).gt.0) ipofano(i)=-ipofano(i)
      enddo
!
!     taking into account the rest (last priority)
!
      do i=1,nface
         indexb=ielfa(2,i)
         iel=ielfa(1,i)
         indexe=ipkon(iel)
         j=ielfa(4,i)
!     
         if(lakon(iel)(4:4).eq.'8') then
            do k=1,4
               node=konf(indexe+ifaceq(k,j))
               if(ipofano(node).eq.0) then
                  ifreefa=ifreefa+1
                  ipofano(node)=ifreefa
                  ifano(1,ifreefa)=i
               elseif(ipofano(node).gt.0) then
                  ifreefa=ifreefa+1
                  ifano(2,ifreefa)=ipofano(node)
                  ifano(1,ifreefa)=i
                  ipofano(node)=ifreefa
               endif
            enddo
         elseif(lakon(iel)(4:4).eq.'6') then
            do k=1,nf(j)
               node=konf(indexe+ifacew(k,j))
               if(ipofano(node).eq.0) then
                  ifreefa=ifreefa+1
                  ipofano(node)=ifreefa
                  ifano(1,ifreefa)=i
               elseif(ipofano(node).gt.0) then
                  ifreefa=ifreefa+1
                  ifano(2,ifreefa)=ipofano(node)
                  ifano(1,ifreefa)=i
                  ipofano(node)=ifreefa
               endif
            enddo
         else
            do k=1,3
               node=konf(indexe+ifacet(k,j))
               if(ipofano(node).eq.0) then
                  ifreefa=ifreefa+1
                  ipofano(node)=ifreefa
                  ifano(1,ifreefa)=i
               elseif(ipofano(node).gt.0) then
                  ifreefa=ifreefa+1
                  ifano(2,ifreefa)=ipofano(node)
                  ifano(1,ifreefa)=i
                  ipofano(node)=ifreefa
               endif
            enddo
         endif
      enddo
!
!     removing all markings
!
      do i=1,nk
         if(ipofano(i).lt.0) ipofano(i)=-ipofano(i)
      enddo
!     
      return
      end
