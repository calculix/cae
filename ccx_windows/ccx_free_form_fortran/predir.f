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
!     y=A*x for real sparse matrices (symmetric and non-symmetric)
!
!     storage of the matrix in a:
!        - first the lower triangular terms
!        - then, if the matrix is non-symmetric, the upper triangular terms
!        - finally the diagonal terms
!
      subroutine predir(n,b,x,nelt,ia,ja,a,isym,itol,tol,itmax,iter,&
        err,ierr,iunit,r,z,dz,rwork,iwork)
      !
      implicit none
      !
      integer n,nelt,ia(*),ja(*),isym,itol,itmax,iter,ierr,&
        iunit,iwork(*)
      !
      real*8 b(*),x(*),a(*),tol,err,r(*),z(*),dz(*),&
        rwork(*)
      !
      external matvec,msolve
      !
      call dir(n,b,x,nelt,ia,ja,a,isym,matvec,msolve,itol,tol,itmax,&
        iter,err,ierr,iunit,r,z,dz,rwork,iwork)
      !
      return
      end
