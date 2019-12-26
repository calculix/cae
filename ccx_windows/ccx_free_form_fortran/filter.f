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
      subroutine filter(dgdxglob,nobject,nk,nodedesi,ndesi,&
                 objectset,xo,yo,zo,x,y,z,nx,ny,nz,neighbor,r,&
                 ndesia,ndesib,xdesi,distmin)
      !
      !     Filtering of sensitivities
      !
      implicit none
      !
      character*81 objectset(4,*)

      integer nobject,nk,nodedesi(*),nnodesinside,i,actdir,&
              ndesi,j,m,neighbor(ndesi+6),nx(ndesi),&
              ny(ndesi),nz(ndesi),istat,ndesia,ndesib
      !
      real*8 dgdxglob(2,nk,nobject),xo(ndesi),yo(ndesi),zo(ndesi),&
             x(ndesi),y(ndesi),z(ndesi),filterrad,r(ndesi+6),&
             filterval(ndesi),nominator,denominator,distmin,&
             xdesi(3,ndesi),scalar,pi,sigma
      !
      !     Calculate filtered sensitivities
      !
      !     Check if direction weighting is turned on
      !
      if(objectset(2,1)(14:16).eq.'DIR') then
         actdir=1
      else
         actdir=0
      endif
      !
      !     Assign filter radius (taken from first defined object function)
      !
      read(objectset(2,1)(21:40),'(f20.0)',iostat=istat) filterrad     
      !
      !     For the GAUSS filter search in the 3sigma distance
      !
      if(objectset(2,1)(1:5).eq.'GAUSS') then
         sigma=filterrad
         filterrad=3*filterrad
      endif
      !
      do j=ndesia,ndesib
         !
         call near3d_se(xo,yo,zo,x,y,z,nx,ny,nz,xo(j),yo(j),zo(j),&
              ndesi,neighbor,r,nnodesinside,filterrad)
         !
         !        Calculate function value of the filterfunction CONST
         !
         if(objectset(2,1)(1:5).eq.'CONST') then
            do i=1,nnodesinside
               filterval(i)=1.0d0
            enddo
         !
         !        Calculate function value of the filterfunction LIN
         !
         elseif(objectset(2,1)(1:3).eq.'LIN') then
            do i=1,nnodesinside
               filterval(i)=(1-dsqrt(r(i))/filterrad)*filterrad
            enddo
         !
         !        Calculate function value of the filterfunction QUAD
         !
         elseif(objectset(2,1)(1:4).eq.'QUAD') then
            do i=1,nnodesinside
               filterval(i)=(-1*(1+dsqrt(r(i))/filterrad)&
                            *(-1+dsqrt(r(i))/filterrad))*filterrad
            enddo
         !
         !        Calculate function value of the filterfunction CUB
         !
         elseif(objectset(2,1)(1:5).eq.'CUBIC') then
            do i=1,nnodesinside
               filterval(i)=(2*(dsqrt(r(i))/filterrad)**3&
                            -3*(dsqrt(r(i))/filterrad)**2+1)*filterrad
            enddo 
         !
         !        Calculate function value of the filterfunction GAUSS
         !
         elseif(objectset(2,1)(1:5).eq.'GAUSS') then
            pi=4.d0*datan(1.d0)
            do i=1,nnodesinside
               filterval(i)=1/(dsqrt(2*pi)*sigma)*exp(-(dsqrt(r(i))**2)&
                            /(2*sigma**2))
            enddo 
         endif
         !
         !        Calculate filtered sensitivity
         !
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle             
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle             
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            nominator=0.d0
            denominator=0.d0
            do i=1,nnodesinside
               if(actdir.eq.1) then          
                  scalar=(xdesi(1,j)*xdesi(1,neighbor(i))&
                      +xdesi(2,j)*xdesi(2,neighbor(i))&
                      +xdesi(3,j)*xdesi(3,neighbor(i)))/(distmin**2)
                  !                   if(objectset(1,m)(1:4).eq.'MASS') then
                  !                      scalar=1.d0
                  !                   endif
                  if(scalar.lt.0.d0) then
                     scalar=0.d0
                  endif
                  nominator=nominator+filterval(i)*scalar*&
                      dgdxglob(1,nodedesi(neighbor(i)),m)
                  denominator=denominator+filterval(i)      
               else
                  nominator=nominator+filterval(i)*&
                      dgdxglob(1,nodedesi(neighbor(i)),m)
                  denominator=denominator+filterval(i)
               endif
            enddo 
            dgdxglob(2,nodedesi(j),m)=nominator/denominator
         enddo
      enddo
      !
      return        
      end




