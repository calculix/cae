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
      subroutine objective_disp_dx(nodeset,istartset,iendset,ialset,&
        nk,idesvarc,iobject,mi,nactdof,dgdx,ndesi,nobject,vold,b,&
        objectset)
      !
      !     calculates the sum of the square of the displacements of a node
      !     set and its derivative w.r.t. the coordinates of the mesh
      !
      implicit none
      !
      character*81 objectset(4,*)
      !
      integer nk,istartset(*),iendset(*),ialset(*),nodeset,idir,&
        idof,idesvarc,iobject,mi(*),nactdof(0:mi(2),*),j,k,ndesi,&
        nobject,idesvar
      !
      real*8 dgdx(ndesi,nobject),vold(0:mi(2),*),b(*)
      !
      intent(in) nodeset,istartset,iendset,ialset,&
        nk,idesvarc,iobject,mi,nactdof,ndesi,nobject,vold,b
      !
      intent(inout) dgdx
      !
      idesvar=idesvarc+1
      !
      !     check for the existence of a set, else take the complete mesh
      !
      if(nodeset.eq.0) then
         do j=1,nk
            if(objectset(1,iobject)(1:12).eq.'DISPLACEMENT') then
               do idir=1,3
                  idof=nactdof(idir,j)
                  if(idof.gt.0) then
                     dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                          +2.d0*vold(idir,j)*b(idof)
                  endif
               enddo
            elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
               idof=nactdof(1,j)
               if(idof.gt.0) then
                  dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                       +2.d0*vold(1,j)*b(idof)
               endif
            elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
               idof=nactdof(2,j)
               if(idof.gt.0) then
                  dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                       +2.d0*vold(2,j)*b(idof)
               endif
            elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
               idof=nactdof(3,j)
               if(idof.gt.0) then
                  dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                       +2.d0*vold(3,j)*b(idof)
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
                        dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                             +2.d0*vold(idir,ialset(j))*b(idof)
                     endif
                  enddo
               elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                  idof=nactdof(1,ialset(j))
                  if(idof.gt.0) then
                     dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                          +2.d0*vold(1,ialset(j))*b(idof)
                  endif
               elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                  idof=nactdof(2,ialset(j))
                  if(idof.gt.0) then
                     dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                          +2.d0*vold(2,ialset(j))*b(idof)
                  endif
               elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                  idof=nactdof(3,ialset(j))
                  if(idof.gt.0) then
                     dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                          +2.d0*vold(3,ialset(j))*b(idof)
                  endif
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(objectset(1,iobject)(1:12).eq.'DISPLACEMENT') then
                     do idir=1,3
                        idof=nactdof(idir,k)
                        if(idof.gt.0) then
                           dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                                +2.d0*vold(idir,k)*b(idof)
                        endif
                     enddo
                  elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                     idof=nactdof(1,k)
                     if(idof.gt.0) then
                        dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                             +2.d0*vold(1,k)*b(idof)
                     endif
                  elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                     idof=nactdof(2,k)
                     if(idof.gt.0) then
                        dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                             +2.d0*vold(2,k)*b(idof)
                     endif
                  elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                     idof=nactdof(3,k)
                     if(idof.gt.0) then
                        dgdx(idesvar,iobject)=dgdx(idesvar,iobject)&
                             +2.d0*vold(3,k)*b(idof)
                     endif
                  endif
               enddo
            endif
         enddo
      endif
      !
      return
      end
      
