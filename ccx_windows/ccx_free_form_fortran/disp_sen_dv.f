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
      subroutine disp_sen_dv(nodeset,istartset,iendset,ialset,iobject,&
        mi,nactdof,dgdu,vold,objectset,nactdofinv,neq)
      !
      !     calculates the sum of the square of the displacements of a node
      !     set and its derivative w.r.t. the coordinates of the mesh
      !
      implicit none
      !
      character*81 objectset(4,*)
      !
      integer istartset(*),iendset(*),ialset(*),nodeset,idir,&
        idof,iobject,mi(*),nactdof(0:mi(2),*),j,k,nactdofinv(*),&
        inode,node,neq,mt
      !
      real*8 dgdu(*),vold(0:mi(2),*)
      !
      intent(in) nodeset,istartset,iendset,ialset,iobject,mi,&
        nactdof,vold
      !
      intent(inout) dgdu
      !
      mt=mi(2)+1
      !
      !     check for the existence of a set, else take the complete mesh
      !
      if(nodeset.eq.0) then
         do idof=1,neq
            inode=nactdofinv(idof)               
            idir=inode-mt*(inode/mt);
            node=inode/mt+1;
            if(objectset(1,iobject)(1:12).eq.'DISPLACEMENT') then
               dgdu(idof)=2.d0*vold(idir,node)
               elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                  if(idir.eq.1) then
                     dgdu(idof)=2.d0*vold(idir,node)
                  endif   
               elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                  if(idir.eq.2) then
                     dgdu(idof)=2.d0*vold(idir,node)
                  endif   
               elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                  if(idir.eq.3) then
                     dgdu(idof)=2.d0*vold(idir,node)
                  endif   
               endif                     
            enddo
      else
         do j=istartset(nodeset),iendset(nodeset)      
            if(ialset(j).gt.0) then
               if(objectset(1,iobject)(1:12).eq.'DISPLACEMENT') then
                  do idir=1,3
                     idof=nactdof(idir,ialset(j))
                     if(idof.gt.0) then
                        dgdu(idof)=2.d0*vold(idir,ialset(j))
                     endif
                  enddo
               elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                  idof=nactdof(1,ialset(j))
                  if(idof.gt.0) then
                     dgdu(idof)=2.d0*vold(1,ialset(j))
                  endif   
               elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                  idof=nactdof(2,ialset(j))
                  if(idof.gt.0) then
                     dgdu(idof)=2.d0*vold(2,ialset(j))
                  endif   
               elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                  idof=nactdof(3,ialset(j))
                  if(idof.gt.0) then
                     dgdu(idof)=2.d0*vold(3,ialset(j))
                  endif   
               endif    
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(objectset(1,iobject)(1:12).eq.'DISPLACEMENT') then
                     do idir=1,3
                        idof=nactdof(idir,ialset(j))
                        if(idof.gt.0) then
                           dgdu(idof)=2.d0*vold(idir,ialset(j))
                        endif
                     enddo
                  elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                     idof=nactdof(1,ialset(j))
                     if(idof.gt.0) then
                        dgdu(idof)=2.d0*vold(1,ialset(j))
                     endif     
                  elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                     idof=nactdof(2,ialset(j))
                     if(idof.gt.0) then
                        dgdu(idof)=2.d0*vold(2,ialset(j))
                     endif     
                  elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                     idof=nactdof(3,ialset(j))
                     if(idof.gt.0) then
                        dgdu(idof)=2.d0*vold(3,ialset(j))
                     endif     
                  endif      
               enddo
            endif
         enddo  
      endif
      !
      return
      end
      
