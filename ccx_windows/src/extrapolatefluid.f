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
      subroutine extrapolatefluid(nk,ipofano,ifano,inum,vfa,v,ielfa,
     &  ithermal,imach,ikappa,xmach,xkappa,shcon,nshcon,ntmat_,ielmat,
     &  physcon,mi,iturb,xturb,gradtfa,gradvfa,gradpfa,gradkfa,gradofa,
     &  co,cofa,ifabou)
!
!     extrapolates the field values at the center of the faces to
!     the nodes
!
      implicit none
!
      logical gradient
!
      integer nk,ipofano(*),ifano(2,*),inum(*),ielfa(4,*),i,l,indexf,
     &  iface,ithermal(*),imach,ikappa,imat,nshcon(*),ntmat_,mi(*),
     &  ielmat(mi(3),*),iturb,ifabou(*)
!
      real*8 vfa(0:7,*),v(0:4,*),cp,r,xk,xmach(*),xkappa(*),t1l,
     &  shcon(0:3,ntmat_,*),physcon(*),xturb(2,*),gradtfa(3,*),
     &  gradvfa(3,3,*),gradpfa(3,*),gradkfa(3,*),gradofa(3,*),
     &  vfal(0:7),co(3,*),cofa(3,*)
!
c      sum=0.d0
!
      do i=1,nk
!
!        athermal calculations
!
         if(ithermal(1).eq.0) then
            do l=1,4
               v(l,i)=0.d0
            enddo
            inum(i)=0
            indexf=ipofano(i)
            do
               if(indexf.eq.0) exit
               iface=ifano(1,indexf)
               if(ielfa(2,iface).gt.0) then
                  do l=1,3
                     vfal(l)=vfa(l,iface)+
     &                    gradvfa(l,1,iface)*(co(1,i)-cofa(1,iface))+
     &                    gradvfa(l,2,iface)*(co(2,i)-cofa(2,iface))+
     &                    gradvfa(l,3,iface)*(co(3,i)-cofa(3,iface))
                  enddo
                  vfal(4)=vfa(4,iface)+
     &                 gradpfa(1,iface)*(co(1,i)-cofa(1,iface))+
     &                 gradpfa(2,iface)*(co(2,i)-cofa(2,iface))+
     &                 gradpfa(3,iface)*(co(3,i)-cofa(3,iface))
               else
                  do l=1,3
                     vfal(l)=vfa(l,iface)
                  enddo
                  vfal(4)=vfa(4,iface)
               endif
               do l=1,4
                  v(l,i)=v(l,i)+vfal(l)
               enddo
               inum(i)=inum(i)+1
               indexf=ifano(2,indexf)
            enddo
            if(inum(i).gt.0) then
               do l=1,4
                  v(l,i)=v(l,i)/inum(i)
               enddo
            endif
         else
!
!           1) thermal incompressible calculations 
!           2) compressible calculations
!
            do l=0,4
               v(l,i)=0.d0
            enddo
            inum(i)=0
            indexf=ipofano(i)
            do
               if(indexf.eq.0) exit
               iface=ifano(1,indexf)
               if(ielfa(2,iface).gt.0) then
!
!                 internal node
!
                  gradient=.true.
               elseif(ielfa(2,iface).eq.0) then
!
!                 external node without boundary condition
!
                  gradient=.false.
               elseif(ielfa(2,iface).lt.0) then
!
!                 external node with boundary condition
!
                  if(ifabou(-ielfa(2,iface)+5).lt.0) then
!
!                    sliding conditions
!
                     gradient=.true.
                  else
!
!                    other conditions
!
                     gradient=.false.
                  endif
               endif
               if(gradient) then
                  vfal(0)=vfa(0,iface)+
     &                 gradtfa(1,iface)*(co(1,i)-cofa(1,iface))+
     &                 gradtfa(2,iface)*(co(2,i)-cofa(2,iface))+
     &                 gradtfa(3,iface)*(co(3,i)-cofa(3,iface))
                  do l=1,3
                     vfal(l)=vfa(l,iface)+
     &                    gradvfa(l,1,iface)*(co(1,i)-cofa(1,iface))+
     &                    gradvfa(l,2,iface)*(co(2,i)-cofa(2,iface))+
     &                    gradvfa(l,3,iface)*(co(3,i)-cofa(3,iface))
                  enddo
                  vfal(4)=vfa(4,iface)+
     &                 gradpfa(1,iface)*(co(1,i)-cofa(1,iface))+
     &                 gradpfa(2,iface)*(co(2,i)-cofa(2,iface))+
     &                 gradpfa(3,iface)*(co(3,i)-cofa(3,iface))
               else
                  vfal(0)=vfa(0,iface)
                  do l=1,3
                     vfal(l)=vfa(l,iface)
                  enddo
                  vfal(4)=vfa(4,iface)
               endif
               do l=0,4
                  v(l,i)=v(l,i)+vfal(l)
               enddo
               if(imach.eq.1) then
                  t1l=vfal(0)
                  imat=ielmat(1,ielfa(1,iface))
                  r=shcon(3,1,imat)
                  call materialdata_cp_sec(imat,ntmat_,t1l,
     &                 shcon,nshcon,cp,physcon)
                  xk=cp/(cp-r)
                  xmach(i)=xmach(i)+dsqrt((vfal(1)**2+
     &               vfal(2)**2+vfal(3)**2)/(xk*r*t1l))
               endif
               if(ikappa.eq.1) then
                  xkappa(i)=xkappa(i)+xk
               endif
!
               inum(i)=inum(i)+1
               indexf=ifano(2,indexf)
            enddo
            if(inum(i).gt.0) then
               do l=0,4
                  v(l,i)=v(l,i)/inum(i)
               enddo
               if(imach.eq.1) xmach(i)=xmach(i)/inum(i)
               if(ikappa.eq.1) xkappa(i)=xkappa(i)/inum(i)
            endif
         endif
!
!        turbulence output
!
         if(iturb.ne.0) then
            do l=1,2
               xturb(l,i)=0.d0
            enddo
            indexf=ipofano(i)
            do
               if(indexf.eq.0) exit
               iface=ifano(1,indexf)
               if(ielfa(2,iface).gt.0) then
                  vfal(6)=vfa(6,iface)+
     &                 gradkfa(1,iface)*(co(1,i)-cofa(1,iface))+
     &                 gradkfa(2,iface)*(co(2,i)-cofa(2,iface))+
     &                 gradkfa(3,iface)*(co(3,i)-cofa(3,iface))
c                  write(*,*) 'extrapolatefluid i1',i,vfal(6)
c                  if(sum.lt.vfal(6)) sum=vfal(6)
                  vfal(7)=vfa(7,iface)+
     &                 gradofa(1,iface)*(co(1,i)-cofa(1,iface))+
     &                 gradofa(2,iface)*(co(2,i)-cofa(2,iface))+
     &                 gradofa(3,iface)*(co(3,i)-cofa(3,iface))
               else
                  vfal(6)=vfa(6,iface)
c                  write(*,*) 'extrapolatefluid i2',i,vfal(6)
c                  if(sum.lt.vfal(6)) sum=vfal(6)
                  vfal(7)=vfa(7,iface)
               endif
               do l=1,2
                  xturb(l,i)=xturb(l,i)+vfal(l+5)
               enddo
c               inum(i)=inum(i)+1
               indexf=ifano(2,indexf)
            enddo
            if(inum(i).gt.0) then
               do l=1,2
                  xturb(l,i)=xturb(l,i)/inum(i)
c                  write(*,*) 'extrapolatesum ',i,xturb(1,i)
c                  if(sum.lt.abs(xturb(1,i))) sum=abs(xturb(1,i))
               enddo
!
!              check for values with an exponent consisting of 3 digits
!              (are not correctly displayed in cgx)
!
               if(xturb(1,i).lt.1.e-30) xturb(1,i)=0.d0
            endif
         endif
      enddo
c      write(*,*) 'extrapolatefluid sum',sum
!  
      return
      end
