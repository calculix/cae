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
      subroutine defplas(elconloc,elas,emec,ithermal,icmd,
     &  beta,stre,ckl,vj)
!
!     calculates stiffness and stresses for the deformation plasticity
!     material law
!
!     icmd=3: calculates stress at mechanical strain
!     else: calculates stress and stiffness matrix at mechanical strain
!
      implicit none
!
      integer cauchy
!
      integer ithermal(*),icmd,i,j,k,l,m,n,ii,istart,iend,nt,kk(84)
!
      real*8 elconloc(*),elas(*),emec(*),beta(*),s(6),al,
     &  ee,un,s0,xn,stre(*),eq,c0,c1,c2,c3,dkl(3,3),ekl(3,3),
     &  q,dq,pp,el(6),ckl(3,3),vj
!
      kk=(/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &  1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,3,3,1,3,
     &  1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3,2,3,
     &  2,3,2,3/)
!
      cauchy=1
!
      istart=1
      iend=1
!
!     determining linear elastic material constants
!
      ee=elconloc(1)
      un=elconloc(2)
      s0=elconloc(3)
      xn=elconloc(4)
      al=elconloc(5)
!
      do i=1,6
         el(i)=emec(i)
      enddo
!
!     major loop
!
      do ii=istart,iend
!
         c0=(el(1)+el(2)+el(3))/3.d0
!
         el(1)=el(1)-c0
         el(2)=el(2)-c0
         el(3)=el(3)-c0
!
!        equivalent deviatoric strain
!
         eq=dsqrt(2.d0/3.d0*(el(1)*el(1)+el(2)*el(2)+
     &                       el(3)*el(3)+2.d0*(el(4)*el(4)+
     &                       el(5)*el(5)+el(6)*el(6))))
!
!        initial value of the Mises equivalent stress (q)
!
         c1=3.d0*ee*eq/(2.d0*(1.d0+un))
!
         if(c1.le.s0) then
            q=c1
         else
            q=(s0**(xn-1)*ee*eq/al)**(1.d0/xn)
         endif
!
!        determining the Mises equivalent stress q
!
         c1=2.d0*(1.d0+un)/3.d0
         do
            c2=al*(q/s0)**(xn-1.d0)
            dq=(ee*eq-(c1+c2)*q)/(c1+xn*c2)
            if((dabs(dq).lt.q*1.d-4).or.(dabs(dq).lt.1.d-10)) exit
            q=q+dq
         enddo
!
         if(icmd.ne.3) then
!
!           calculating the tangent stiffness matrix
!
!           initialization of the Delta Dirac function
!
            do i=1,3
               do j=1,3
                  dkl(i,j)=0.d0
               enddo
            enddo
            do i=1,3
               dkl(i,i)=1.d0
            enddo
!
            ekl(1,1)=el(1)
            ekl(2,2)=el(2)
            ekl(3,3)=el(3)
            ekl(1,2)=el(4)
            ekl(1,3)=el(5)
            ekl(2,3)=el(6)
            ekl(2,1)=ekl(1,2)
            ekl(3,1)=ekl(1,3)
            ekl(3,2)=ekl(2,3)
!
            if(eq.lt.1.d-10) then
               c1=ee/(1.d0+un)
               c2=0.d0
            else
               c1=2.d0/(3.d0*eq)
               c2=c1*(1.d0/eq-1.d0/(eq+(xn-1.d0)*c2*q/ee))
               c1=c1*q
            endif
            c3=(ee/(1.d0-2.d0*un)-c1)/3.d0
!
            nt=0
            do i=1,21
               k=kk(nt+1)
               l=kk(nt+2)
               m=kk(nt+3)
               n=kk(nt+4)
               nt=nt+4
               elas(i)=c1*((dkl(k,m)*dkl(l,n)+dkl(k,n)*dkl(l,m))/2.d0
     &                     -c2*ekl(k,l)*ekl(m,n))
     &                 +c3*dkl(k,l)*dkl(m,n)
            enddo
!
!           conversion of the stiffness matrix from spatial coordinates
!           coordinates into material coordinates
!
            call stiff2mat(elas,ckl,vj,cauchy)
!
         endif
!
!        calculating the stress
!
         pp=-ee*c0/(1.d0-2.d0*un)
!
         if(eq.lt.1.d-10) then
            c1=0.d0
         else
            c1=2.d0*q/(3.d0*eq)
         endif
!
         do i=1,6
            s(i)=el(i)*c1
         enddo
         do i=1,3
            s(i)=s(i)-pp
         enddo
!
         do i=1,6
            stre(i)=s(i)
         enddo
!
!        converting the stress into the material frame of
!        reference
!
         call str2mat(stre,ckl,vj,cauchy)
!
      enddo
!
      return
      end
