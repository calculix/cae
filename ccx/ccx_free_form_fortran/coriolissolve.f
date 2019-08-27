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
      subroutine coriolissolve(cc,nev,a,b,x,eiga,eigb,eigcorio,iter,d,&
        temp)
      !
      !     solves for the complex eigenfrequencies due to Coriolis
      !     forces
      !
      implicit none
      !
      logical wantx
      !
      integer nev,iter(*),n,i,j,k,l,m,iii,ksmall,nincrement,&
        nevcomplex,&
        ninc,id,ipos
      !
      real*8 cc(nev,*),d(*),size
      !
      complex*16 a(nev,*),b(nev,*),x(nev,*),eiga(*),eigb(*),eigcorio(*),&
        temp(nev,*),deig,eig,eignew
      !
      n=nev
      wantx=.false.
      nincrement=10
      nevcomplex=0
      ninc=2
      ipos=0
      !
      !     loop over all modes
      !
      loop2: do i=1,nev
      if(d(i).lt.0.d0) cycle
         !
         !        initial values are taken at regular intervals between
         !        the eigenvalues without Coriolis
         !
         loop1: do m=1,nincrement
            !
            !           initial frequency
            !
            if(ipos.eq.0) then
               eig=dsqrt(d(i))*(1.d0,0.d0)*m/nincrement
            else
               eig=(dsqrt(d(i-1))+&
                (dsqrt(d(i))-dsqrt(d(i-1)))*m/nincrement)*(1.d0,0.d0)
            endif
            !
            iii=0
            do
               do k=1,nev
                  do l=1,nev
                     a(k,l)=eig*cc(k,l)*(0.d0,1.d0)
                     b(k,l)=-cc(k,l)*(0.d0,1.d0)
                  enddo
                  a(k,k)=a(k,k)-eig*eig+d(k)
                  b(k,k)=b(k,k)+2.d0*eig
               enddo
               !
               !     solving for the complex eigenvalues
               !
               call dlzhes(n,a,n,b,n,temp,n,wantx)
               call dlzit(n,a,n,b,n,temp,n,wantx,iter,eiga,eigb)
               !
               if(iter(1).eq.-1) then
                  write(*,*) '*ERROR in coriolissolve: fatal error'
                  write(*,*) '       in dlzit'
                  call exit(201)
               elseif(cdabs(eigb(1)).lt.1.d-10) then
                  write(*,*) '*ERROR in coriolissolve: eigenvalue'
                  write(*,*) '       out of bounds'
                  deig=0.d0
               else
                  size=1.d30
                  do k=1,nev
                     deig=eiga(k)/eigb(k)
                     if(cdabs(deig).lt.size) then
                        ksmall=k
                        size=cdabs(deig)
                     endif
                  enddo
                  deig=eiga(ksmall)/eigb(ksmall)
               endif
               iii=iii+1
               !                write(*,*) 'coriolissolve ',iii,eig
               !
               if((cdabs(deig).lt.1.d-3*cdabs(eig)).or.&
                    (cdabs(deig).lt.1.d-3).or.(iii.eq.100)) then
                  !                   write(*,*) 'csolve ',i,m,(eig+deig)
                  !
                  id=0
                  eignew=eig+deig
                  !
                  if((i.ne.1).or.(m.ne.1)) then
                     call ident2(eigcorio,eignew,nevcomplex,ninc,id)
                     if(id.ne.0) then
                        if((cdabs(eignew-eigcorio(id)).lt.&
                             2.d-3*cdabs(eignew)).or.&
                           (cdabs(eignew).lt.1.d-3)) cycle loop1
                     endif
                     if(id.ne.nevcomplex) then
                        if((cdabs(eignew-eigcorio(id+1)).lt.&
                             2.d-3*cdabs(eignew)).or.&
                           (cdabs(eignew).lt.1.d-3)) cycle loop1
                     endif
                  endif
                  !
                  !     solve for the eigenvalues AND eigenvectors
                  !
                  do k=1,nev
                     do l=1,nev
                        a(k,l)=eig*cc(k,l)*(0.d0,1.d0)
                        b(k,l)=-cc(k,l)*(0.d0,1.d0)
                     enddo
                     a(k,k)=a(k,k)-eig*eig+d(k)
                     b(k,k)=b(k,k)+2.d0*eig
                  enddo
                  wantx=.true.
                  call dlzhes(n,a,n,b,n,temp,n,wantx)
                  call dlzit(n,a,n,b,n,temp,n,wantx,iter,eiga,eigb)
                  !
                  !     copying the eigenvector into the right place
                  !
                  nevcomplex=nevcomplex+1
                  !
                  do k=nevcomplex,id+2,-1
                     eigcorio(k)=eigcorio(k-1)
                     do j=1,nev
                        x(j,k)=x(j,k-1)
                     enddo
                  enddo
                  !
                  do j=1,nev
                     x(j,id+1)=temp(j,ksmall)
                  enddo
                  eigcorio(id+1)=eignew
                  !
                  if(nevcomplex.eq.nev) exit loop2
                  exit
               else
                  eig=eig+deig
               endif
            !
            enddo
         enddo loop1
         ipos=1
      enddo loop2
      !
      return
      end
      
