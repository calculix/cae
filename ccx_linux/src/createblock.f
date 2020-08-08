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
      subroutine createblock(nef,ipnei,neiel,iam,jam,iamorig,
     &  nflnei,nz_num,num_cpus,nestart,ineighblock,neighblock,
     &  icyclic,neifa,ifatie,nneighblock)
!
!     1) removes zero neighbors from field ia
!     2) removes cyclic symmetric neighbors from field ia
!     3) generates num_cpu blocks
!
      integer i,j,nef,ipnei(*),neiel(*),iam(*),jam(*),iamorig(*),
     &  indexf,nflnei,nz_num,kflag,n,num_cpus,nestart(*),numd,
     &  ifreeneigh,ineighblock(*),neighblock(3,*),neifa(*),ifatie(*),
     &  icyclic,ifa,iel,isubtract,nneighblock
! 
!     nestart(i) points to the element before the block for which
!     cpu i is responsible
!
      numd=int(nef/num_cpus)+1
      nestart(1)=0
      do i=2,num_cpus
         nestart(i)=nestart(i-1)+numd
      enddo
      nestart(num_cpus+1)=nef
!
      ifreeneigh=0
      j=0
!
      do k=1,num_cpus
         ineighblock(k)=ifreeneigh
         do i=nestart(k)+1,nestart(k+1)
            jam(i)=j
            do indexf=ipnei(i)+1,ipnei(i+1)
               iel=neiel(indexf)
!
!              no neighbor
!
               if(iel.eq.0) cycle
!
!              cyclic symmetric neighbor
!
               if(icyclic.eq.1) then
                  ifa=neifa(indexf)
                  if(ifatie(ifa).ne.0) cycle
               endif
!
!              neighbor belongs to adjacent block
!
               if((iel.le.nestart(k)).or.
     &            (iel.gt.nestart(k+1))) then
                  ifreeneigh=ifreeneigh+1
!
!                 location in au/auv
!
                  neighblock(1,ifreeneigh)=indexf
!
!                 neighboring block element number 
!
                  neighblock(2,ifreeneigh)=iel
!
!                 equation number
!
                  neighblock(3,ifreeneigh)=i
                  cycle
               endif
!
!              genuine element of current block
!               
               j=j+1
               iam(j)=iel
               iamorig(j)=indexf
            enddo
!     
!     adding the diagonal element
!     
            j=j+1
            iam(j)=i
            iamorig(j)=nflnei+i
         enddo
      enddo
!
      jam(nef+1)=j
      ineighblock(num_cpus+1)=ifreeneigh
!
      nz_num=j
      nneighblock=ifreeneigh
!
!     sorting the column numbers within each row
!     localizing the diagonal elements and storing them in uam
!
      kflag=2
      do i=1,nef
         n=jam(i+1)-jam(i)
         call isortii(iam(jam(i)+1),iamorig(jam(i)+1),n,kflag)
      enddo
!
!     subtracting from iam the number of elements belonging
!     to the preceding blocks
!
      do k=2,num_cpus
         isubtract=nestart(k)
         do i=nestart(k)+1,nestart(k+1)
            do j=jam(i)+1,jam(i+1)
               iam(j)=iam(j)-isubtract
            enddo
         enddo
      enddo
!
      return
      end
