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
      subroutine orthonl(w,vo,elas,s,ii1,jj1,weight)
      !
      !     This routine replaces the following lines in e_c3d.f for
      !     an orthotropic material
      !
      !                      do i1=1,3
      !                        iii1=ii1+i1-1
      !                        do j1=1,3
      !                          jjj1=jj1+j1-1
      !                          do k1=1,3
      !                            do l1=1,3
      !                              s(iii1,jjj1)=s(iii1,jjj1)
      !     &                              +anisox(i1,k1,j1,l1)*w(k1,l1)
      !                              do m1=1,3
      !                                s(iii1,jjj1)=s(iii1,jjj1)
      !     &                              +anisox(i1,k1,m1,l1)*w(k1,l1)
      !     &                                 *vo(j1,m1)
      !     &                              +anisox(m1,k1,j1,l1)*w(k1,l1)
      !     &                                 *vo(i1,m1)
      !                                do n1=1,3
      !                                  s(iii1,jjj1)=s(iii1,jjj1)
      !     &                                  +anisox(m1,k1,n1,l1)
      !     &                                  *w(k1,l1)*vo(i1,m1)*vo(j1,n1)
      !                                enddo
      !                              enddo
      !                            enddo
      !                          enddo
      !                        enddo
      !                      enddo
      !
      implicit none
      !
      integer ii1,jj1
      !
      real*8 w(3,3),vo(3,3),elas(21),s(60,60),weight
      !
      intent(in) w,vo,elas,ii1,jj1,weight
      !
      intent(inout) s
      !
      s(ii1,jj1)=s(ii1,jj1)+((elas( 1)+elas( 1)*vo(1,1)&
      +(elas( 1)+elas( 1)*vo(1,1)&
      )*vo(1,1)+(elas( 7)*vo(1,2))*vo(1,2)&
      +(elas( 8)*vo(1,3))&
      *vo(1,3))*w(1,1)&
      +(elas( 2)*vo(1,2)+(elas( 2)*vo(1,2))*vo(1,1)+(elas( 7)&
      +elas( 7)*vo(1,1))*vo(1,2)&
      )*w(1,2)&
      +(elas( 4)*vo(1,3)+(elas( 4)*vo(1,3))*vo(1,1)&
      +(elas( 8)+elas( 8)*vo(1,1))&
      *vo(1,3))*w(1,3)&
      +(elas( 7)*vo(1,2)+(elas( 7)*vo(1,2))*vo(1,1)+(elas( 2)&
      +elas( 2)*vo(1,1))*vo(1,2)&
      )*w(2,1)&
      +(elas( 7)+elas( 7)*vo(1,1)&
      +(elas( 7)+elas( 7)*vo(1,1)&
      )*vo(1,1)+(elas( 3)*vo(1,2))*vo(1,2)&
      +(elas( 9)*vo(1,3))&
      *vo(1,3))*w(2,2)&
      +((elas( 5)*vo(1,3))*vo(1,2)&
      +(elas( 9)*vo(1,2))&
      *vo(1,3))*w(2,3)&
      +(elas( 8)*vo(1,3)+(elas( 8)*vo(1,3))*vo(1,1)&
      +(elas( 4)+elas( 4)*vo(1,1))&
      *vo(1,3))*w(3,1)&
      +((elas( 9)*vo(1,3))*vo(1,2)&
      +(elas( 5)*vo(1,2))&
      *vo(1,3))*w(3,2)&
      +(elas( 8)+elas( 8)*vo(1,1)&
      +(elas( 8)+elas( 8)*vo(1,1)&
      )*vo(1,1)+(elas( 9)*vo(1,2))*vo(1,2)&
      +(elas( 6)*vo(1,3))&
      *vo(1,3))*w(3,3))*weight
      s(ii1,jj1+1)=s(ii1,jj1+1)+((elas( 1)*vo(2,1)&
      +(elas( 1)*vo(2,1)&
      )*vo(1,1)+(elas( 7)&
      +elas( 7)*vo(2,2))*vo(1,2)&
      +(elas( 8)*vo(2,3))&
      *vo(1,3))*w(1,1)&
      +(elas( 2)&
      +elas( 2)*vo(2,2)+(elas( 2)&
      +elas( 2)*vo(2,2))*vo(1,1)+(elas( 7)*vo(2,1))*vo(1,2)&
      )*w(1,2)&
      +(elas( 4)*vo(2,3)+(elas( 4)*vo(2,3))*vo(1,1)&
      +(elas( 8)*vo(2,1))&
      *vo(1,3))*w(1,3)&
      +(elas( 7)&
      +elas( 7)*vo(2,2)+(elas( 7)&
      +elas( 7)*vo(2,2))*vo(1,1)+(elas( 2)*vo(2,1))*vo(1,2)&
      )*w(2,1)&
      +(elas( 7)*vo(2,1)&
      +(elas( 7)*vo(2,1)&
      )*vo(1,1)+(elas( 3)&
      +elas( 3)*vo(2,2))*vo(1,2)&
      +(elas( 9)*vo(2,3))&
      *vo(1,3))*w(2,2)&
      +((elas( 5)*vo(2,3))*vo(1,2)&
      +(elas( 9)+elas( 9)*vo(2,2))&
      *vo(1,3))*w(2,3)&
      +(elas( 8)*vo(2,3)+(elas( 8)*vo(2,3))*vo(1,1)&
      +(elas( 4)*vo(2,1))&
      *vo(1,3))*w(3,1)&
      +((elas( 9)*vo(2,3))*vo(1,2)&
      +(elas( 5)+elas( 5)*vo(2,2))&
      *vo(1,3))*w(3,2)&
      +(elas( 8)*vo(2,1)&
      +(elas( 8)*vo(2,1)&
      )*vo(1,1)+(elas( 9)&
      +elas( 9)*vo(2,2))*vo(1,2)&
      +(elas( 6)*vo(2,3))&
      *vo(1,3))*w(3,3))*weight
      s(ii1,jj1+2)=s(ii1,jj1+2)+((elas( 1)*vo(3,1)&
      +(elas( 1)*vo(3,1)&
      )*vo(1,1)+(elas( 7)*vo(3,2))*vo(1,2)&
      +(elas( 8)+elas( 8)*vo(3,3))&
      *vo(1,3))*w(1,1)&
      +(elas( 2)*vo(3,2)&
      +(elas( 2)*vo(3,2))*vo(1,1)+(elas( 7)*vo(3,1))*vo(1,2)&
      )*w(1,2)&
      +(elas( 4)&
      +elas( 4)*vo(3,3)+(elas( 4)&
      +elas( 4)*vo(3,3))*vo(1,1)&
      +(elas( 8)*vo(3,1))&
      *vo(1,3))*w(1,3)&
      +(elas( 7)*vo(3,2)+(elas( 7)*vo(3,2))*vo(1,1)&
      +(elas( 2)*vo(3,1))*vo(1,2)&
      )*w(2,1)&
      +(elas( 7)*vo(3,1)&
      +(elas( 7)*vo(3,1)&
      )*vo(1,1)+(elas( 3)*vo(3,2))*vo(1,2)&
      +(elas( 9)+elas( 9)*vo(3,3))&
      *vo(1,3))*w(2,2)&
      +((elas( 5)&
      +elas( 5)*vo(3,3))*vo(1,2)&
      +(elas( 9)*vo(3,2))&
      *vo(1,3))*w(2,3)&
      +(elas( 8)&
      +elas( 8)*vo(3,3)+(elas( 8)&
      +elas( 8)*vo(3,3))*vo(1,1)&
      +(elas( 4)*vo(3,1))&
      *vo(1,3))*w(3,1)&
      +((elas( 9)&
      +elas( 9)*vo(3,3))*vo(1,2)&
      +(elas( 5)*vo(3,2))&
      *vo(1,3))*w(3,2)&
      +(elas( 8)*vo(3,1)&
      +(elas( 8)*vo(3,1)&
      )*vo(1,1)+(elas( 9)*vo(3,2))*vo(1,2)&
      +(elas( 6)+elas( 6)*vo(3,3))&
      *vo(1,3))*w(3,3))*weight
      s(ii1+1,jj1)=s(ii1+1,jj1)+((elas( 7)*vo(1,2)&
      +(elas( 1)+elas( 1)*vo(1,1)&
      )*vo(2,1)+(elas( 7)*vo(1,2))*vo(2,2)&
      +(elas( 8)*vo(1,3))&
      *vo(2,3))*w(1,1)&
      +(elas( 7)+elas( 7)*vo(1,1)&
      +(elas( 2)*vo(1,2))*vo(2,1)+(elas( 7)&
      +elas( 7)*vo(1,1))*vo(2,2)&
      )*w(1,2)&
      +((elas( 4)*vo(1,3))*vo(2,1)&
      +(elas( 8)+elas( 8)*vo(1,1))&
      *vo(2,3))*w(1,3)&
      +(elas( 2)+elas( 2)*vo(1,1)&
      +(elas( 7)*vo(1,2))*vo(2,1)+(elas( 2)&
      +elas( 2)*vo(1,1))*vo(2,2)&
      )*w(2,1)&
      +(elas( 3)*vo(1,2)+(elas( 7)+elas( 7)*vo(1,1)&
      )*vo(2,1)+(elas( 3)*vo(1,2))*vo(2,2)&
      +(elas( 9)*vo(1,3))&
      *vo(2,3))*w(2,2)&
      +(elas( 5)*vo(1,3)+(elas( 5)*vo(1,3))*vo(2,2)&
      +(elas( 9)*vo(1,2))&
      *vo(2,3))*w(2,3)&
      +((elas( 8)*vo(1,3))*vo(2,1)&
      +(elas( 4)+elas( 4)*vo(1,1))&
      *vo(2,3))*w(3,1)&
      +(elas( 9)*vo(1,3)+(elas( 9)*vo(1,3))*vo(2,2)&
      +(elas( 5)*vo(1,2))&
      *vo(2,3))*w(3,2)&
      +(elas( 9)*vo(1,2)+(elas( 8)+elas( 8)*vo(1,1)&
      )*vo(2,1)+(elas( 9)*vo(1,2))*vo(2,2)&
      +(elas( 6)*vo(1,3))&
      *vo(2,3))*w(3,3))*weight
      s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+((elas( 7)&
      +elas( 7)*vo(2,2)+(elas( 1)*vo(2,1)&
      )*vo(2,1)+(elas( 7)&
      +elas( 7)*vo(2,2))*vo(2,2)&
      +(elas( 8)*vo(2,3))&
      *vo(2,3))*w(1,1)&
      +(elas( 7)*vo(2,1)&
      +(elas( 2)&
      +elas( 2)*vo(2,2))*vo(2,1)+(elas( 7)*vo(2,1))*vo(2,2)&
      )*w(1,2)&
      +((elas( 4)*vo(2,3))*vo(2,1)&
      +(elas( 8)*vo(2,1))&
      *vo(2,3))*w(1,3)&
      +(elas( 2)*vo(2,1)&
      +(elas( 7)&
      +elas( 7)*vo(2,2))*vo(2,1)+(elas( 2)*vo(2,1))*vo(2,2)&
      )*w(2,1)&
      +(elas( 3)&
      +elas( 3)*vo(2,2)+(elas( 7)*vo(2,1)&
      )*vo(2,1)+(elas( 3)&
      +elas( 3)*vo(2,2))*vo(2,2)&
      +(elas( 9)*vo(2,3))&
      *vo(2,3))*w(2,2)&
      +(elas( 5)*vo(2,3)+(elas( 5)*vo(2,3))*vo(2,2)&
      +(elas( 9)+elas( 9)*vo(2,2))&
      *vo(2,3))*w(2,3)&
      +((elas( 8)*vo(2,3))*vo(2,1)&
      +(elas( 4)*vo(2,1))&
      *vo(2,3))*w(3,1)&
      +(elas( 9)*vo(2,3)+(elas( 9)*vo(2,3))*vo(2,2)&
      +(elas( 5)+elas( 5)*vo(2,2))&
      *vo(2,3))*w(3,2)&
      +(elas( 9)&
      +elas( 9)*vo(2,2)+(elas( 8)*vo(2,1)&
      )*vo(2,1)+(elas( 9)&
      +elas( 9)*vo(2,2))*vo(2,2)&
      +(elas( 6)*vo(2,3))&
      *vo(2,3))*w(3,3))*weight
      s(ii1+1,jj1+2)=s(ii1+1,jj1+2)+((elas( 7)*vo(3,2)+(elas( 1)*vo(3,1)&
      )*vo(2,1)+(elas( 7)*vo(3,2))*vo(2,2)&
      +(elas( 8)+elas( 8)*vo(3,3))&
      *vo(2,3))*w(1,1)&
      +(elas( 7)*vo(3,1)&
      +(elas( 2)*vo(3,2))*vo(2,1)+(elas( 7)*vo(3,1))*vo(2,2)&
      )*w(1,2)&
      +((elas( 4)&
      +elas( 4)*vo(3,3))*vo(2,1)&
      +(elas( 8)*vo(3,1))&
      *vo(2,3))*w(1,3)&
      +(elas( 2)*vo(3,1)&
      +(elas( 7)*vo(3,2))*vo(2,1)+(elas( 2)*vo(3,1))*vo(2,2)&
      )*w(2,1)&
      +(elas( 3)*vo(3,2)+(elas( 7)*vo(3,1)&
      )*vo(2,1)+(elas( 3)*vo(3,2))*vo(2,2)&
      +(elas( 9)+elas( 9)*vo(3,3))&
      *vo(2,3))*w(2,2)&
      +(elas( 5)&
      +elas( 5)*vo(3,3)+(elas( 5)&
      +elas( 5)*vo(3,3))*vo(2,2)&
      +(elas( 9)*vo(3,2))&
      *vo(2,3))*w(2,3)&
      +((elas( 8)&
      +elas( 8)*vo(3,3))*vo(2,1)&
      +(elas( 4)*vo(3,1))&
      *vo(2,3))*w(3,1)&
      +(elas( 9)&
      +elas( 9)*vo(3,3)+(elas( 9)&
      +elas( 9)*vo(3,3))*vo(2,2)&
      +(elas( 5)*vo(3,2))&
      *vo(2,3))*w(3,2)&
      +(elas( 9)*vo(3,2)+(elas( 8)*vo(3,1)&
      )*vo(2,1)+(elas( 9)*vo(3,2))*vo(2,2)&
      +(elas( 6)+elas( 6)*vo(3,3))&
      *vo(2,3))*w(3,3))*weight
      s(ii1+2,jj1)=s(ii1+2,jj1)+((elas( 8)*vo(1,3)&
      +(elas( 1)+elas( 1)*vo(1,1)&
      )*vo(3,1)+(elas( 7)*vo(1,2))*vo(3,2)&
      +(elas( 8)*vo(1,3))&
      *vo(3,3))*w(1,1)&
      +((elas( 2)*vo(1,2))*vo(3,1)+(elas( 7)&
      +elas( 7)*vo(1,1))*vo(3,2)&
      )*w(1,2)&
      +(elas( 8)+elas( 8)*vo(1,1)&
      +(elas( 4)*vo(1,3))*vo(3,1)&
      +(elas( 8)+elas( 8)*vo(1,1))&
      *vo(3,3))*w(1,3)&
      +((elas( 7)*vo(1,2))*vo(3,1)+(elas( 2)&
      +elas( 2)*vo(1,1))*vo(3,2)&
      )*w(2,1)&
      +(elas( 9)*vo(1,3)+(elas( 7)+elas( 7)*vo(1,1)&
      )*vo(3,1)+(elas( 3)*vo(1,2))*vo(3,2)&
      +(elas( 9)*vo(1,3))&
      *vo(3,3))*w(2,2)&
      +(elas( 9)*vo(1,2)+(elas( 5)*vo(1,3))*vo(3,2)&
      +(elas( 9)*vo(1,2))&
      *vo(3,3))*w(2,3)&
      +(elas( 4)+elas( 4)*vo(1,1)&
      +(elas( 8)*vo(1,3))*vo(3,1)&
      +(elas( 4)+elas( 4)*vo(1,1))&
      *vo(3,3))*w(3,1)&
      +(elas( 5)*vo(1,2)+(elas( 9)*vo(1,3))*vo(3,2)&
      +(elas( 5)*vo(1,2))&
      *vo(3,3))*w(3,2)&
      +(elas( 6)*vo(1,3)+(elas( 8)+elas( 8)*vo(1,1)&
      )*vo(3,1)+(elas( 9)*vo(1,2))*vo(3,2)&
      +(elas( 6)*vo(1,3))&
      *vo(3,3))*w(3,3))*weight
      s(ii1+2,jj1+1)=s(ii1+2,jj1+1)+((elas( 8)*vo(2,3)&
      +(elas( 1)*vo(2,1)&
      )*vo(3,1)+(elas( 7)&
      +elas( 7)*vo(2,2))*vo(3,2)&
      +(elas( 8)*vo(2,3))&
      *vo(3,3))*w(1,1)&
      +((elas( 2)&
      +elas( 2)*vo(2,2))*vo(3,1)+(elas( 7)*vo(2,1))*vo(3,2)&
      )*w(1,2)&
      +(elas( 8)*vo(2,1)&
      +(elas( 4)*vo(2,3))*vo(3,1)&
      +(elas( 8)*vo(2,1))&
      *vo(3,3))*w(1,3)&
      +((elas( 7)&
      +elas( 7)*vo(2,2))*vo(3,1)+(elas( 2)*vo(2,1))*vo(3,2)&
      )*w(2,1)&
      +(elas( 9)*vo(2,3)+(elas( 7)*vo(2,1)&
      )*vo(3,1)+(elas( 3)&
      +elas( 3)*vo(2,2))*vo(3,2)&
      +(elas( 9)*vo(2,3))&
      *vo(3,3))*w(2,2)&
      +(elas( 9)&
      +elas( 9)*vo(2,2)+(elas( 5)*vo(2,3))*vo(3,2)&
      +(elas( 9)+elas( 9)*vo(2,2))&
      *vo(3,3))*w(2,3)&
      +(elas( 4)*vo(2,1)&
      +(elas( 8)*vo(2,3))*vo(3,1)&
      +(elas( 4)*vo(2,1))&
      *vo(3,3))*w(3,1)&
      +(elas( 5)&
      +elas( 5)*vo(2,2)+(elas( 9)*vo(2,3))*vo(3,2)&
      +(elas( 5)+elas( 5)*vo(2,2))&
      *vo(3,3))*w(3,2)&
      +(elas( 6)*vo(2,3)+(elas( 8)*vo(2,1)&
      )*vo(3,1)+(elas( 9)&
      +elas( 9)*vo(2,2))*vo(3,2)&
      +(elas( 6)*vo(2,3))&
      *vo(3,3))*w(3,3))*weight
      s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+((elas( 8)&
      +elas( 8)*vo(3,3)+(elas( 1)*vo(3,1)&
      )*vo(3,1)+(elas( 7)*vo(3,2))*vo(3,2)&
      +(elas( 8)+elas( 8)*vo(3,3))&
      *vo(3,3))*w(1,1)&
      +((elas( 2)*vo(3,2))*vo(3,1)+(elas( 7)*vo(3,1))*vo(3,2)&
      )*w(1,2)&
      +(elas( 8)*vo(3,1)&
      +(elas( 4)&
      +elas( 4)*vo(3,3))*vo(3,1)&
      +(elas( 8)*vo(3,1))&
      *vo(3,3))*w(1,3)&
      +((elas( 7)*vo(3,2))*vo(3,1)+(elas( 2)*vo(3,1))*vo(3,2)&
      )*w(2,1)&
      +(elas( 9)&
      +elas( 9)*vo(3,3)+(elas( 7)*vo(3,1)&
      )*vo(3,1)+(elas( 3)*vo(3,2))*vo(3,2)&
      +(elas( 9)+elas( 9)*vo(3,3))&
      *vo(3,3))*w(2,2)&
      +(elas( 9)*vo(3,2)+(elas( 5)&
      +elas( 5)*vo(3,3))*vo(3,2)&
      +(elas( 9)*vo(3,2))&
      *vo(3,3))*w(2,3)&
      +(elas( 4)*vo(3,1)&
      +(elas( 8)&
      +elas( 8)*vo(3,3))*vo(3,1)&
      +(elas( 4)*vo(3,1))&
      *vo(3,3))*w(3,1)&
      +(elas( 5)*vo(3,2)+(elas( 9)&
      +elas( 9)*vo(3,3))*vo(3,2)&
      +(elas( 5)*vo(3,2))&
      *vo(3,3))*w(3,2)&
      +(elas( 6)&
      +elas( 6)*vo(3,3)+(elas( 8)*vo(3,1)&
      )*vo(3,1)+(elas( 9)*vo(3,2))*vo(3,2)&
      +(elas( 6)+elas( 6)*vo(3,3))&
      *vo(3,3))*w(3,3))*weight
      !
      return
      end
