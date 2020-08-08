!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine applybounv(nodeboun,ndirboun,nboun,v,nmpc,nodempc,
     &     ipompc,coefmpc,inomat,mi)
!
!     boundary conditions for Delta V* (CFD, CBS method)
!
!     1) applies velocity SPC's
!     2) applies MPC's for the conservative variables v(1..3,*)    
!     
      implicit none
!     
      integer mi(*),nodeboun(*),ndirboun(*),i,nboun,node,index,
     &     nmpc,nodempc(3,*),ipompc(*),ndir,inomat(*)
!     
      real*8 v(0:mi(2),*),coefmpc(*),residu,correction
!     
      do i=1,nboun
!     
!     check whether fluid node
!     
        node=nodeboun(i)
        if(inomat(node).eq.0) cycle
        if((ndirboun(i).ge.1).and.(ndirboun(i).le.3)) then
          v(ndirboun(i),node)=0.d0
        endif
      enddo
!     
!     MPC's are treated by distributing the residual proportional to
!     the coefficients
!     
!     Right now it is assumed that the MPC's are independent of each 
!     other, i.e. degrees of freedom used in one MPC are not used in 
!     any other MPC
!     
      do i=1,nmpc
        index=ipompc(i)
!     
!     check whether fluid node
!     
        node=nodempc(1,index)
        if(inomat(node).eq.0) cycle
!     
!     pressure multiple point constraints for incompressible materials
!     have already been treated
!     
        ndir=nodempc(2,index)
        if((ndir.ge.1).and.(ndir.le.3)) then
          residu=coefmpc(index)*v(ndir,node)
          if(index.eq.0) cycle
          do
            index=nodempc(3,index)
            if(index.eq.0) exit
            node=nodempc(1,index)
            ndir=nodempc(2,index)
            residu=residu+coefmpc(index)*v(ndir,node)
          enddo
!     
!     correcting all terms of the MPC
!     
          index=ipompc(i)
          do
            node=nodempc(1,index)
            ndir=nodempc(2,index)
!     
            correction=-residu*coefmpc(index)
            v(ndir,node)=v(ndir,node)+correction
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        endif
      enddo
!     
      return
      end
      
