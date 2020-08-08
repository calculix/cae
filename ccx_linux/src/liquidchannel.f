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
      subroutine liquidchannel(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,rho,g,co,dvi,numf,mi,ipkon,kon,iplausi)
!
!     open channel for incompressible media
!
!     SG: sluice gate
!     SO: sluice opening
!     WE: weir
!     WO: weir opening
!     DS: discontinuous slope
!     DO: discontinuous slope opening
!       : default channel mit linearly varying trapezoid cross
!         section
!     
      implicit none
!     
      logical identity,bresse,jump
      character*8 lakon(*)
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,indexup,i,
     &     ielprop(*),nodef(*),idirf(*),index,iflag,mi(*),nsol,
     &     inv,numf,nodesg,nelemdown,nelemup,node0,kon(*),ipkon(*),
     &     iplausi
!      
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),b,d,c,p,
     &     h1,h2,rho,dvi,friction,reynolds,dg,b1,b2,
     &     g(3),dl,xks,z1,z2,co(3,*),xflow2,dyg3dbj,dyg4dbj,
     &     s0,sqrts0,hk,form_fact,h1ns,h2ns,h0,dyg3deta,dyg4deta,
     &     dh3dh1,dh4dh2,dh3dm,dh4dm,eta,dA3deta,dA4deta,bj,
     &     theta,cth,tth,um1,um2,A1,A2,P1,P2,D1,D2,dA1dh1,dA2dh2,
     &     dP1dh1,dP2dh2,dD1dh1,dD2dh2,h3,h4,dh3deta,xn1,xn2,xt1,xt2,
     &     dh4deta,yg3,yg4,dyg3dh3,dyg4dh4,A3,A4,dA3dh3,dA4dh4,
     &     dum1dh1,dum2dh2,c1,c2,dbds,dbjdeta,e0,e1,e2,e3,
     &     dyg3dm,dyg4dm,dA3dm,dA4dm,dyg3dh1,dyg4dh2,
     &     dA3dh1,dA4dh2,solreal(3),solimag(3),dist
!
!
!
!     iflag=0: check whether all parameters in the element equation
!              are known => equation is not needed
!     iflag=1: calculation of the initial flux
!     iflag=2: evaluate the element equation and all derivatives
!     iflag=3: correct the channel depth in order to move a jump
!
      if (iflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif((iflag.eq.1).or.(iflag.eq.2))then
!     
         index=ielprop(nelem)
!     
         h1=v(2,node1)
         h2=v(2,node2)
!     
         z1=-g(1)*co(1,node1)-g(2)*co(2,node1)-g(3)*co(3,node1)
         z2=-g(1)*co(1,node2)-g(2)*co(2,node2)-g(3)*co(3,node2)
!
         dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!     
         if(iflag.eq.1) then
!
!           in a first call of liquidchannel the flow is determined,
!           in a second call the channel depth is calculated
!
            if(lakon(nelem)(6:7).eq.'SG') then
!
!              sluice gate
!
               b=prop(index+1)
               s0=prop(index+2)
               if(s0.lt.-1.d0) then
                  s0=dasin((z1-z2)/dl)
               endif
               sqrts0=dsqrt(1.d0-s0*s0)
               theta=0.d0
               h2=prop(index+3)
!
               if(dabs(xflow).lt.1.d-30) then
!
!                 determine initial mass flow
!                  
                  if(nactdog(2,node1).ne.0) then
!
!                    upstream level not known
!
                     xflow=0.d0
                  else
                     xflow=2.d0*dg*(rho*b*h2)**2*(h1-h2*sqrts0)
                     if(xflow.lt.0.d0) then
                        write(*,*)'*ERROR in liquidchannel: water level'
                        write(*,*) '       upstream of sluice gate is '
                        write(*,*) '       smaller than downstream heigh
     &t'
                        call exit(201)
                     else
                        xflow=dsqrt(xflow)
                     endif
                  endif
               else
!
!                 determine the downstream depth
!                 and the upstream depth if not defined as BC
!
                  call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
                  v(3,node2)=hk
                  if(h2.gt.hk) then
!
!                    for initial conditions
!
                     if(nactdog(2,node1).ne.0) v(2,node1)=3.d0*hk/2.d0
                     v(2,node2)=hk
                  else
!
!                    for initial conditions
!
                     if(nactdog(2,node1).ne.0) v(2,node1)=
     &                  xflow**2/(2.d0*dg*(rho*b*h2)**2)+h2*sqrts0
                     v(2,node2)=h2
                  endif
               endif
            elseif(lakon(nelem)(6:7).eq.'WE') then
!
!              weir
!
               b=prop(index+1)
               p=prop(index+2)
               c=prop(index+3)
               sqrts0=1.d0
               theta=0.d0
!
               if(dabs(xflow).lt.1.d-30) then
!
!                 determine initial mass flow
!                  
                  if(nactdog(2,node1).ne.0) then
!
!                    upstream level unknown
!
                     xflow=0.d0
                  else
                     if(h1.le.p) then
                        write(*,*) '*ERROR in liquidchannel'
                        write(*,*) '       weir height exceeds'
                        write(*,*) '       upstream level'
                        call exit(201)
                      endif
                      xflow=rho*c*b*(h1-p)**(1.5d0)
                   endif
               else
!
!                 determine the downstream depth
!                 and the upstream depth if not defined as BC
!
                  call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
                  v(3,node2)=hk
!
!                 for initial conditions
!
                  if(nactdog(2,node1).ne.0) v(2,node1)=p+3.d0*hk/2.d0
!
!                 next value is used for downstream initial values
!
                  v(2,node2)=hk
               endif
!
            elseif(lakon(nelem)(6:7).eq.'DS') then
               if(dabs(xflow).lt.1.d-30) then
!
!                 initial mass flow cannot be determined for this
!                 type of element
!                  
                  xflow=0.d0
               else
!     
!              determine the downstream depth
!
                  b=prop(index+1)
                  s0=prop(index+2)
                  if(s0.lt.-1.d0) then
                     s0=dasin((z1-z2)/dl)
                  endif
                  sqrts0=dsqrt(1.d0-s0*s0)
                  theta=prop(index+4)
!     
                  call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
                  v(3,node2)=hk
!
!                 initial condition for fluid depth
!                 supercritical value
!
                  v(2,node2)=hk/2.d0
               endif
!               
            endif
         else
!
!           calculating f and its derivatives
!
            bresse=.false.
            jump=.false.
!
            xflow2=xflow*xflow
!
!           element properties
!
            if((lakon(nelem)(6:7).eq.'SG').or.
     &         (lakon(nelem)(6:7).eq.'SO').or.
     &         (lakon(nelem)(6:7).eq.'WO').or.
     &         (lakon(nelem)(6:7).eq.'RE').or.
     &         (lakon(nelem)(6:7).eq.'  ').or.
     &         (lakon(nelem)(6:7).eq.'DS').or.
     &         (lakon(nelem)(6:7).eq.'DO')) then
               b=prop(index+1)
               s0=prop(index+2)
               if(s0.lt.-1.d0) then
                  s0=dasin((z1-z2)/dl)
               endif
               sqrts0=dsqrt(1.d0-s0*s0)
               if(lakon(nelem)(6:7).ne.'SG') then
                  dl=prop(index+3)
                  theta=prop(index+4)
                  xks=prop(index+5)
                  if(dl.le.0.d0) then
                     dl=dsqrt((co(1,node2)-co(1,node1))**2+
     &                    (co(2,node2)-co(2,node1))**2+
     &                    (co(3,node2)-co(3,node1))**2)
                  endif
               else
                  theta=0.d0
               endif
            elseif(lakon(nelem)(6:7).eq.'WE') then
               b=prop(index+1)
               p=prop(index+2)
               c=prop(index+3)
               sqrts0=1.d0
               theta=0.d0
            elseif((lakon(nelem)(6:7).eq.'CO').or.
     &             (lakon(nelem)(6:7).eq.'EL')) then
               b1=prop(index+1)
!
               s0=prop(index+2)
               if(s0.lt.-1.d0) then
                  s0=0.d0
               endif
               sqrts0=dsqrt(1.d0-s0*s0)
!
               dl=prop(index+3)
               if(dl.le.0.d0) then
                  dl=dsqrt((co(1,node2)-co(1,node1))**2+
     &                 (co(2,node2)-co(2,node1))**2+
     &                 (co(3,node2)-co(3,node1))**2)
               endif
!
               b2=prop(index+4)
               b=(b1+b2)/2.d0
               theta=0.d0
               xks=0.d0
            elseif((lakon(nelem)(6:7).eq.'ST').or.
     &             (lakon(nelem)(6:7).eq.'DR')) then
               b=prop(index+1)
!
               s0=prop(index+2)
               if(s0.lt.-1.d0) then
                  s0=0.d0
               endif
               sqrts0=dsqrt(1.d0-s0*s0)
!
               dl=prop(index+3)
               if(dl.le.0.d0) then
                  dl=dsqrt((co(1,node2)-co(1,node1))**2+
     &                 (co(2,node2)-co(2,node1))**2+
     &                 (co(3,node2)-co(3,node1))**2)
               endif
!
               d=prop(index+4)
               b1=b
               b2=b
               theta=0.d0
               xks=0.d0
            endif
!
            if(xflow.ge.0.d0) then
               inv=1
            else
               inv=-1
            endif
!
!           standard element equation: unknowns are the mass flow
!           and the depth upstream and downstream
!
            numf=3
            nodef(1)=node1
            nodef(2)=nodem
            nodef(3)=node2
            idirf(1)=2
            idirf(2)=1
            idirf(3)=2
!
            if(lakon(nelem)(6:7).eq.'SG') then
!
!              sluice gate
!              1-SG-2-SO-3
!
!              h2 cannot exceed HKmax
!
               h2=prop(index+3)
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
               if(h2.gt.hk) h2=hk
!
               nelemdown=nint(prop(index+5))
               h3=v(2,kon(ipkon(nelemdown)+3))
               call hns(b,theta,rho,dg,sqrts0,xflow,h2,h2ns)
               if(h3.lt.h2ns) then
!
!                 Q=f_SG(h1,h2): sluice gate equation between
!                 1 and 2
!
!                 next line for output only
!
                  v(2,node2)=h2
c                  write(30,*) 'SG: sluice gate equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3,'h2ns= ',h2ns
                  df(1)=2.d0*dg*(rho*b*h2)**2
                  df(2)=-2.d0*xflow
                  f=df(1)*(h1-h2*sqrts0)
                  df(3)=2.d0*f/h2-df(1)*sqrts0
                  f=f-xflow2
               else
!
!                 fake equation
!                  
c                  write(30,*) 'SG: fake equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3,'h2ns= ',h2ns
                  numf=1
                  nodef(1)=nodem
                  idirf(1)=3
                  f=prop(index+4)-0.5d0
                  df(1)=1.d0
               endif
            elseif(lakon(nelem)(6:7).eq.'SO') then
!
!              sluice opening (element streamdown of sluice gate)
!              0-SG-1-SO-2
!
               nelemup=nint(prop(index+6))
               node0=kon(ipkon(nelemup)+1)
               h0=v(2,node0)
               h1=prop(ielprop(nelemup)+3)
!
!              h1 cannot exceed HKmax
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
               if(h1.gt.hk) h1=hk
!
               call hns(b,theta,rho,dg,sqrts0,xflow,h1,h1ns)
               if(h2.lt.h1ns) then
!
!                 bresse (frontwater)
!
c                  write(30,*) 'SO: Bresse equation '
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2,'h1ns= ',h1ns
                  bresse=.true.
               else
!
!                 Q=f_SG(h0,h2): sluice gate equation between 0 and 2
!                 (backwater)
!
!                 reset gate height
!                  
                  h1=prop(ielprop(nelemup)+3)
!
c                  write(30,*) 'SO: Sluice gate eqn. between 0 and 2 '
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2,'h1ns= ',h1ns
                  numf=4
                  nodef(4)=node0
                  idirf(4)=2
!
                  if(h2.gt.h1) then
!
!                    gate flow (water touches gate)
!                    section = b * h1
!
!                    next line for output only
!
                     v(2,node1)=h1
                     df(4)=2.d0*dg*(rho*b*h1)**2
                     df(3)=-df(4)*sqrts0
                     df(2)=-2.d0*xflow
                     f=df(4)*(h0-h2*sqrts0)
                     df(1)=2.d0*f/h1
                  else
!
!                    incomplete inflexion (water does not touch gate)
!                    section = b * h2
!
!                    next line for output only
!
                     v(2,node1)=h2
                     df(4)=2.d0*dg*(rho*b*h2)**2
                     df(3)=-df(4)*sqrts0
                     df(2)=-2.d0*xflow
                     f=df(4)*(h0-h2*sqrts0)
                     df(3)=df(3)+2.d0*f/h2
                     df(1)=0.d0
                  endif
                  f=f-xflow2
               endif
            elseif(lakon(nelem)(6:7).eq.'WE') then
!
!              weir
!              1-WE-2-WO-3
!
               nelemdown=nint(prop(index+5))
               h3=v(2,kon(ipkon(nelemdown)+3))
!
!              default depth for weir is hk
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(h3.lt.p+hk) then
!
!                 output only
!
                  v(2,node2)=p+hk
!
!                 Q=f_WE(h1): weir equation
!
c                  write(30,*) 'WE: weir equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3,'hk= ',hk
                  f=rho*c*b*(h1-p)**(1.5d0)
                  df(1)=3.d0*f/(2.d0*(h1-p))
                  f=f-xflow
                  df(2)=-1.d0
                  df(3)=0.d0
               else
!
!                 fake equation
!                  
c                  write(30,*) 'WE: weir equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3,'hk= ',hk
                  numf=1
                  nodef(1)=nodem
                  idirf(1)=3
                  f=prop(index+4)-0.5d0
                  df(1)=1.d0
               endif
            elseif(lakon(nelem)(6:7).eq.'WO') then
!
!              weir opening (element streamdown of weir)
!              0-WE-1-WO-2
!
               nelemup=nint(prop(index+6))
               node0=kon(ipkon(nelemup)+1)
               h0=v(2,node0)
!
               p=prop(ielprop(nelemup)+2)
!
!              default depth for weir is hk
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(h2.lt.p+hk) then
!
!                 bresse between 1 and 2
!
                  h1=hk
c                  write(30,*) 'WO: Bresse equation '
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2,'hk= ',hk
                  p=prop(ielprop(nelemup)+2)
                  s0=dasin(p/dsqrt(dl**2+p**2))
c                  write(*,*) 's0=',p,dl,s0
                  sqrts0=dsqrt(1.d0-s0*s0)
                  bresse=.true.
               else
!
!                 output only
!
                  v(2,node1)=h2
!
!                 bresse between 0 and 2
!
c                  write(30,*) 'WO: Bresse eqn. between 0 and 2 '
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2,'hk= ',hk
                  nodef(1)=node0
                  h1=h0
                  bresse=.true.
               endif
            elseif(lakon(nelem)(6:7).eq.'DS') then
!     
!     discontinuous slope
!     1-DS-2-DO-3
!     
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!     
               if(h1.gt.hk) then
                  nelemdown=nint(prop(index+8))
                  h3=v(2,kon(ipkon(nelemdown)+3))
                  if(h3.le.hk) then
!     
!                    upstream: backwater curve
!                    downstream: frontwater curve
!     
                     h2=hk
                     bresse=.true.
c                  write(30,*) 'DS:  back/front bresse'
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3
!
!                    for output purposes
!
                     v(2,node2)=h2
                  else
!     
!                    both curves are backwater curves
!                    fake equation
!     
c                  write(30,*) 'DS:  back/back fake equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2,'h3= ',h3
                     numf=1
                     nodef(1)=nodem
                     idirf(1)=3
                     f=prop(index+7)-0.5d0
                     df(1)=1.d0
                  endif
               else
!     
!                 both curves are frontwater curves
!                 fake equation
!     
c                  write(30,*) 'DS:  front/front fake equation '
c                  write(30,*)'h1= ',h1,'h2= ',h2
                  nelemup=nint(prop(index+6))
                  numf=1
                  nodef(1)=kon(ipkon(nelemup)+2)
                  idirf(1)=3
                  f=prop(index+7)-0.5d0
                  df(1)=1.d0
               endif
            elseif(lakon(nelem)(6:7).eq.'DO') then
!     
!              discontinuous slope opening 
!              (element streamdown of discontinuous slope)
!              0-DS-1-DO-2
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               nelemup=nint(prop(index+6))
               node0=kon(ipkon(nelemup)+1)
               h0=v(2,node0)
!
               if(h0.gt.hk) then
                  if(h2.le.hk) then
!
!                    upstream: backwater curve
!                    downstream: frontwater curve
!                    bresse between 1 and 2
!     
                     h1=hk
c                  write(30,*) 'DO: back/front bresse 1-2'
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2
                     bresse=.true.
                  else
!     
!                    both curves are backwater curves
!                    bresse between 0 and 2
!     
c                  write(30,*) 'DO: back/back bresse 0-2'
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2
                     nodef(1)=node0
                     h1=h0
                     bresse=.true.
!
!                    output purposes
!
                     v(2,node1)=(h0+h2)/2.d0
                  endif
               else
!     
!                 both curves are frontwater curves
!                 bresse between 0 and 2
!     
c                  write(30,*) 'DO: front/front bresse 0-2'
c                  write(30,*)'h0= ',h0,'h1= ',h1,'h2= ',h2
                  nodef(1)=node0
                  h1=h0
                  bresse=.true.
!
!                 output purposes
!
                  v(2,node1)=(h0+h2)/2.d0
               endif
            elseif(lakon(nelem)(6:7).eq.'RE') then
!
!              element upstream of a reservoir
!              calculating the critical depth
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
               if(h1.ge.hk) then
!
!                 backwater curve
!
                  if(h2.lt.hk) h2=hk
c                  write(30,*) 'RE: Bresse downstream equation '
c                  write(30,*) 'h1= ',h1,'h2= ',h2,'hk= ',hk
                  bresse=.true.
               else
!
!                 frontwater curve
!
                  call hns(b,theta,rho,dg,sqrts0,xflow,h1,h1ns)
                  if(h2.le.h1ns) then
c                  write(30,*) 'RE: fake equation '
c                  write(30,*) 'h1= ',h1,'h2= ',h2,'h1ns= ',h1ns
!
!                    fake equation
!                     
                     nelemup=nint(prop(index+6))
                     nodesg=kon(ipkon(nelemup)+2)
                     numf=1
                     nodef(1)=nodesg
                     idirf(1)=3
!
!                    retrieving previous value of eta
!
                     index=ielprop(nelemup)
                     if(lakon(nelemup)(6:7).eq.'SG') then
                        f=prop(index+4)-0.5d0
                     elseif(lakon(nelemup)(6:7).eq.'WE') then
                        f=prop(index+4)-0.5d0
                     elseif(lakon(nelemup)(6:7).eq.'DS') then
                        f=prop(index+7)-0.5d0
                     endif
                     df(1)=1.d0
                  else
c                  write(30,*) 'RE: Bresse downstream equation '
c                  write(30,*) 'h1= ',h1,'h2= ',h2,'h1ns= ',h1ns
                     bresse=.true.
                  endif
               endif
            elseif(lakon(nelem)(6:7).eq.'CO') then
c               write(30,*) 'CO: contraction '
c               write(30,*)'h1= ',h1,'h2= ',h2
!
               call hcrit(xflow,rho,b2,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(inv.eq.-1) then
                  if((h1.gt.hk).and.(h2.lt.hk)) then
                     jump=.true.
                  endif
               else
                  if((h1.lt.hk).and.(h2.gt.hk)) then
                     jump=.true.
                  endif
               endif
!
c               write(*,*) 'CO ',jump
!
               if(.not.jump) then
                  c1=rho*rho*dg
                  c2=b1*b2*h1*h2
                  df(1)=b1*(2.d0*xflow2+c1*b1*b2*h2**3)
                  df(3)=b2*(2.d0*xflow2+c1*b1*b1*h1**3)
                  f=h1*df(1)-h2*df(3)
                  df(1)=df(1)-3.d0*c1*c2*b1*h1
                  df(3)=3.d0*c1*c2*b1*h2-df(3)
                  df(2)=4.d0*(b1*h1-b2*h2)*xflow
               endif
            elseif(lakon(nelem)(6:7).eq.'EL') then
c               write(30,*) 'EL: enlargement '
c               write(30,*)'h1= ',h1,'h2= ',h2
!
               call hcrit(xflow,rho,b2,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(inv.eq.-1) then
                  if((h1.gt.hk).and.(h2.lt.hk)) then
                     jump=.true.
                  endif
               else
                  if((h1.lt.hk).and.(h2.gt.hk)) then
                     jump=.true.
                  endif
               endif
!
c               write(*,*) 'EL ',jump
!
               if(.not.jump) then
                  c1=rho*rho*dg
                  c2=b1*b2*h1*h2
                  df(1)=b1*(2.d0*xflow2+c1*b2*b2*h2**3)
                  df(3)=b2*(2.d0*xflow2+c1*b1*b2*h1**3)
                  f=h1*df(1)-h2*df(3)
                  df(1)=df(1)-3.d0*c1*c2*b2*h1
                  df(3)=3.d0*c1*c2*b2*h2-df(3)
                  df(2)=4.d0*(b1*h1-b2*h2)*xflow
               endif
            elseif(lakon(nelem)(6:7).eq.'DR') then
c               write(30,*) 'DR: drop '
c               write(30,*)'h1= ',h1,'h2= ',h2
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(inv.eq.-1) then
                  if((h1.gt.hk).and.(h2.lt.hk)) then
                     jump=.true.
                  endif
               else
                  if((h1.lt.hk).and.(h2.gt.hk)) then
                     jump=.true.
                  endif
               endif
!
               if(.not.jump) then
                  c1=rho*rho*dg
                  df(1)=2.d0*xflow2+c1*b*b*h2**3
                  df(3)=2.d0*xflow2+c1*b*b*h1*(h1+d)**2
                  f=h1*df(1)-h2*df(3)
                  df(1)=df(1)-c1*b*b*h2*(3.d0*h1+d)*(h1+d)
                  df(3)=3.d0*c1*b*b*h1*h2*h2-df(3)
                  df(2)=4.d0*(h1-h2)*xflow
               endif
            elseif(lakon(nelem)(6:7).eq.'ST') then
c               write(30,*) 'ST: step '
c               write(30,*)'h1= ',h1,'h2= ',h2
!
               call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
               v(3,node2)=hk
!
               if(inv.eq.-1) then
                  if((h1.gt.hk).and.(h2.lt.hk)) then
                     jump=.true.
                  endif
               else
                  if((h1.lt.hk).and.(h2.gt.hk)) then
                     jump=.true.
                  endif
               endif
!
               if(.not.jump) then
                  c1=rho*rho*dg
                  df(1)=2.d0*xflow2+c1*b*b*h2*(h2+d)**2
                  df(3)=2.d0*xflow2+c1*b*b*h1**3
                  f=h1*df(1)-h2*df(3)
                  df(1)=df(1)-3.d0*c1*b*b*h1*h1*h2
                  df(3)=c1*b*b*h1*(3.d0*h2+d)*(h2+d)-df(3)
                  df(2)=4.d0*(h1-h2)*xflow
               endif
            elseif(lakon(nelem)(6:7).eq.'  ') then
               bresse=.true.
c                  write(30,*)  'straight: Bresse equation '
c                  write(30,*) 'h1= ',h1,'h2= ',h2
            endif
!
!           bresse equation
!
            if((bresse).or.(jump)) then
!
               if(xks.gt.0.d0) then
!
!                 White-Coolebrook
!
!                 hydraulic diameter
!
                  d=2.d0*(h1+h2)
                  reynolds=4.d0*xflow/(b*dvi)
                  form_fact=1.d0
                  call friction_coefficient(dl,d,xks,reynolds,form_fact,
     &                 friction)
               endif
!
               if(bresse) then
                  call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
                  v(3,node2)=hk
                  if(inv.eq.-1) then
                     if((h1.gt.hk).and.(h2.lt.hk)) then
                        jump=.true.
                     endif
                  else
                     if((h1.lt.hk).and.(h2.gt.hk)) then
                        jump=.true.
                     endif
                  endif
                  b1=b
                  b2=b
               endif
!
!              geometric data
!
               cth=dcos(theta)
               tth=dtan(theta)
!
!              nonprismatic cross section
!
               if(lakon(nelem)(6:7).eq.'  ') then
                  dbds=prop(index+7)
               else
                  dbds=0.d0
               endif
!
!              width at water surface
!
               dD1dh1=2.d0*tth
               dD2dh2=dD1dh1
               D1=b1+h1*dD1dh1
               D2=b2+dl*dbds+h2*dD2dh2
!
!              cross section
!
               A1=h1*(b1+h1*tth)
               A2=h2*(b2+dl*dbds+h2*tth)
               dA1dh1=D1
               dA2dh2=D2
!
!              perimeter
!
               P1=b1+2.d0*h1/cth
               P2=b2+dl*dbds+2.d0*h2/cth
               dP1dh1=2.d0/cth
               dP2dh2=dP1dh1
!
!              factor for friction
!
               if(xks.gt.0.d0) then
!                 White-Coolebrook
                  um1=friction/8.d0
                  um2=um1
                  dum1dh1=0.d0
                  dum2dh2=0.d0
               else
!                 Manning
                  um1=xks*xks*dg*(P1/A1)**(1.d0/3.d0)
                  um2=xks*xks*dg*(P2/A2)**(1.d0/3.d0)
                  dum1dh1=xks*xks*dg*
     &                    (P1**(-2.d0/3.d0)*dP1dh1*A1**(1.d0/3.d0)-
     &                     A1**(-2.d0/3.d0)*dA1dh1*P1**(1.d0/3.d0))/
     &                     (3.d0*A1**(2.d0/3d0))
                  dum2dh2=xks*xks*dg*
     &                    (P2**(-2.d0/3.d0)*dP2dh2*A2**(1.d0/3.d0)-
     &                     A2**(-2.d0/3.d0)*dA2dh2*P2**(1.d0/3.d0))/
     &                     (3.d0*A2**(2.d0/3d0))
               endif
!
!              constants
!
               c1=rho*rho*dg
               c2=c1*sqrts0
               c1=c1*s0
!
!              hydraulic jump
!
               if(jump) then
c                  write(30,*) 
c     &              'liquidchannel: jump in element,hk ',nelem,hk
                  nelemup=prop(index+6)
                  indexup=ielprop(nelemup)
                  if(lakon(nelemup)(6:7).eq.'SG') then
                     eta=prop(indexup+4)
                     prop(indexup+7)=nelem
                  elseif(lakon(nelemup)(6:7).eq.'WE') then
                     eta=prop(indexup+4)
                     prop(indexup+7)=nelem
                  elseif(lakon(nelemup)(6:7).eq.'DS') then
                     eta=prop(indexup+7)
                     prop(indexup+9)=nelem
                  endif
!
!                 determining h3, h4 and derivatives
!
!                 numerator
!
                  xt1=c1*A1**3+(h1*dbds-um1*P1)*xflow2
                  xt2=c1*A2**3+(h2*dbds-um2*P2)*xflow2
!
!                 denominator
!
                  xn1=c2*A1**3-D1*xflow2
                  xn2=c2*A2**3-D2*xflow2
!
!                 h3 and h4
!
                  h3=h1+dl*xt1/xn1*eta
                  h4=h2-dl*xt2/xn2*(1.d0-eta)
c                  write(30,*) 
c     &              'liquidchannel: h3,h4,eta ',h3,h4,eta
!
                  if(bresse) then
!     
!                    width at jump
!     
                     bj=b+dbds*eta*dl
!     
!                    cross sections and derivatives
!
                     A3=h3*(bj+h3*tth)
                     A4=h4*(bj+h4*tth)
                     dA3dh3=bj+2.d0*h3*tth
                     dA4dh4=bj+2.d0*h4*tth
!
!                    center of gravity and derivatives
!
                     yg3=h3*(3.d0*bj+2.d0*h3*tth)/(6.d0*(bj+h3*tth))
                     yg4=h4*(3.d0*bj+2.d0*h4*tth)/(6.d0*(bj+h4*tth))
                     dyg3dh3=((3.d0*bj+4.d0*h3*tth)*(bj+tth)
     &                    -tth*h3*(3.d0*bj+2.d0*h3*tth))/
     &                    (6.d0*(bj+h3*tth)**2)
                     dyg4dh4=((3.d0*bj+4.d0*h4*tth)*(bj+tth)
     &                    -tth*h4*(3.d0*bj+2.d0*h4*tth))/
     &                    (6.d0*(bj+h4*tth)**2)
                     dyg3dbj=h3*h3*tth/(6.d0*(bj+h3*tth)**2)
                     dyg4dbj=h4*h4*tth/(6.d0*(bj+h4*tth)**2)
                  endif
!     
!                 derivative of h3 w.r.t. h1 and of h4 w.r.t. h2
!     
                  dh3dh1=1.d0+((3.d0*c1*A1*A1*dA1dh1
     &                   +(dbds-dum1dh1*P1-um1*dP1dh1)*xflow2)*xn1
     &                   -(3.d0*c2*A1*A1*dA1dh1-dD1dh1*xflow2)*xt1)/
     &                    (xn1*xn1)*eta*dl
                  dh4dh2=1.d0-((3.d0*c1*A2*A2*dA2dh2
     &                   +(dbds-dum2dh2*P2-um2*dP2dh2)*xflow2)*xn2
     &                   -(3.d0*c2*A2*A2*dA2dh2-dD2dh2*xflow2)*xt2)/
     &                    (xn2*xn2)*(1.d0-eta)*dl
!
                  if(bresse) then
                     dA3dh1=dA3dh3*dh3dh1
                     dA4dh2=dA4dh4*dh4dh2
                     dyg3dh1=dyg3dh3*dh3dh1
                     dyg4dh2=dyg4dh4*dh4dh2
                  endif
!
!                 derivative of h3 and h4 w.r.t. the mass flow
!
                  dh3dm=((dbds*h1-um1*P1)*xn1+D1*xt1)*2.d0*xflow/
     &                  (xn1*xn1)*eta*dl
                  dh4dm=-((dbds*h2-um2*P2)*xn2+D2*xt2)*2.d0*xflow/
     &                  (xn2*xn2)*(1.d0-eta)*dl
!
                  if(bresse) then
                     dA3dm=dA3dh3*dh3dm
                     dA4dm=dA4dh4*dh4dm
                     dyg3dm=dyg3dh3*dh3dm
                     dyg4dm=dyg4dh4*dh4dm
                  endif
!
!                 derivative of h3 and h4 w.r.t. eta
!
                  dh3deta=dl*xt1/xn1
                  dh4deta=dl*xt2/xn2
!
                  if(bresse) then
                     dbjdeta=dbds*dl
!
!                    derivative of A3, A4, yg3 and yg4 w.r.t. eta
!
                     dA3deta=dA3dh3*dh3deta+h3*dbjdeta
                     dA4deta=dA4dh4*dh4deta+h4*dbjdeta
                     dyg3deta=dyg3dh3*dh3deta+dyg3dbj*dbjdeta
                     dyg4deta=dyg4dh4*dh4deta+dyg4dbj*dbjdeta
                  endif
!
                  numf=4
                  nodef(4)=kon(ipkon(nelemup)+2)
                  idirf(4)=3
!
                  if(bresse) then
                     f=A4*xflow2+c2*(A3*A3*A4*yg3-A3*A4*A4*yg4)
     &                    -A3*xflow2
                     df(1)=c2*(2.d0*A3*dA3dh1*A4*yg3+A3*A3*A4*dyg3dh1
     &                    -dA3dh1*A4*A4*yg4)-dA3dh1*xflow2
                     df(2)=2.d0*xflow*(A4-A3)+
     &                    (c2*(2.d0*A3*A4*yg3-A4*A4*yg4)-xflow2)*dA3dm+
     &                    (c2*(A3*A3*yg3-2.d0*A3*A4*yg4)+xflow2)*dA4dm+
     &                    c2*A3*A3*A4*dyg3dm-c2*A3*A4*A4*dyg4dm
                     df(3)=c2*(A3*A3*dA4dh2*yg3-2.d0*A3*A4*dA4dh2*yg4
     &                    -A3*A4*A4*dyg4dh2)+dA4dh2*xflow2
                     df(4)=dA4deta*xflow2+
     &                    c2*(2.d0*A3*dA3deta*A4*yg3+A3*A3*dA4deta*yg3
     &                    +A3*A3*A4*dyg3deta-dA3deta*A4*A4*yg4
     &                    -A3*2.d0*A4*dA4deta*yg4-A3*A4*A4*dyg4deta)
     &                    -dA3deta*xflow2
                  elseif(lakon(nelem)(6:7).eq.'CO') then
                     f=b2*h4*(2.d0*xflow2+c2*b1*b1*h3**3)-
     &                 b1*h3*(2.d0*xflow2+c2*b1*b2*h4**3)
!                    dfdh3
                     df(1)=3.d0*b2*h4*c2*b1*b1*h3*h3-
     &                 b1*(2.d0*xflow2+c2*b1*b2*h4**3)
!                    dfdh4
                     df(3)=b2*(2.d0*xflow2+c2*b1*b1*h3**3)-
     &                     3.d0*b1*h3*c2*b1*b2*h4*h4
!                    dfdm
                     df(2)=4.d0*xflow*(b2*h4-b1*h3)+
     &                     df(1)*dh3dm+df(3)*dh4dm
!                    dfdeta
                     df(4)=df(1)*dh3deta+df(3)*dh4deta
!                    dfdh1
                     df(1)=df(1)*dh3dh1
!                    dfdh2
                     df(3)=df(3)*dh4dh2
                  elseif(lakon(nelem)(6:7).eq.'EL') then
                     f=b2*h4*(2.d0*xflow2+c2*b1*b2*h3**3)-
     &                 b1*h3*(2.d0*xflow2+c2*b2*b2*h4**3)
!                    dfdh3
                     df(1)=3.d0*b2*h4*c2*b1*b2*h3*h3-
     &                 b1*(2.d0*xflow2+c2*b2*b2*h4**3)
!                    dfdh4
                     df(3)=b2*(2.d0*xflow2+c2*b1*b2*h3**3)-
     &                     3.d0*b1*h3*c2*b2*b2*h4*h4
!                    dfdm
                     df(2)=4.d0*xflow*(b2*h4-b1*h3)+
     &                     df(1)*dh3dm+df(3)*dh4dm
!                    dfdeta
                     df(4)=df(1)*dh3deta+df(3)*dh4deta
!                    dfdh1
                     df(1)=df(1)*dh3dh1
!                    dfdh2
                     df(3)=df(3)*dh4dh2
                  elseif(lakon(nelem)(6:7).eq.'DR') then
                     f=h4*(2.d0*xflow2+c2*b*b*h3*(h3+d)**2)-
     &                 h3*(2.d0*xflow2+c2*b*b*h4**3)
!                    dfdh3
                     df(1)=h4*c2*b*b*(3.d0*h3+d)*(h3+d)-
     &                     (2.d0*xflow2+c2*b*b*h4**3)
!                    dfdh4
                     df(3)=(2.d0*xflow2+c2*b*b*h3*(h3+d)**2)-
     &                     3.d0*h3*c2*b*b*h4*h4
!                    dfdm
                     df(2)=4.d0*xflow*(h4-h3)+
     &                     df(1)*dh3dm+df(3)*dh4dm
!                    dfdeta
                     df(4)=df(1)*dh3deta+df(3)*dh4deta
!                    dfdh1
                     df(1)=df(1)*dh3dh1
!                    dfdh2
                     df(3)=df(3)*dh4dh2
                  elseif(lakon(nelem)(6:7).eq.'ST') then
                     f=h4*(2.d0*xflow2+c2*b*b*h3**3)-
     &                 h3*(2.d0*xflow2+c2*b*b*h4*(h4+d)**2)
!                    dfdh3
                     df(1)=3.d0*h4*c2*b*b*h3*h3-
     &                     (2.d0*xflow2+c2*b*b*h4*(h4+d)**2)
!                    dfdh4
                     df(3)=(2.d0*xflow2+c2*b*b*h3**3)-
     &                     h3*c2*b*b*(3.d0*h4+d)*(h4+d)
!                    dfdm
                     df(2)=4.d0*xflow*(h4-h3)+
     &                     df(1)*dh3dm+df(3)*dh4dm
!                    dfdeta
                     df(4)=df(1)*dh3deta+df(3)*dh4deta
!                    dfdh1
                     df(1)=df(1)*dh3dh1
!                    dfdh2
                     df(3)=df(3)*dh4dh2
                  endif
               else
!
!                 regular Bresse equation
!
                  f=c2*(A1**3+A2**3)-xflow2*(D1+D2)
                  df(1)=-f+(h2-h1)*(c2*dA1dh1*3.d0*A1*A1-xflow2*dD1dh1)
     &                  -dl*(c1*3.d0*A1*A1*dA1dh1
     &                       -(dum1dh1*P1+um1*dP1dh1-dbds)*xflow2)
                  df(2)=(-(h2-h1)*(D1+D2)
     &                   +dl*(um1*P1+um2*P2-(h1+h2)*dbds))*2.d0*xflow
                  df(3)=f+(h2-h1)*(c2*dA2dh2*3.d0*A2*A2-xflow2*dD2dh2)
     &                  -dl*(c1*3.d0*A2*A2*dA2dh2
     &                       -(dum2dh2*P2+um2*dP2dh2-dbds)*xflow2)
                  f=(h2-h1)*f-dl*(c1*(A1**3+A2**3)
     &                            -(um1*P1+um2*P2-(h1+h2)*dbds)*xflow2)
               endif
            endif
         endif
      elseif(iflag.eq.3) then
!
!        only if called from resultnet in case the element contains
!        a hydraulic jump and eta<0 or eta>1. This means that the 
!        jump does not take place in the element itself. By adjusting
!        h1 or h2 the jump is forced into a neighboring element
!
         index=ielprop(nelem)
c         write(30,*) 'iflag=3, nelem',nelem,lakon(nelem)
!     
         h1=v(2,node1)
         h2=v(2,node2)
!     
         z1=-g(1)*co(1,node1)-g(2)*co(2,node1)-g(3)*co(3,node1)
         z2=-g(1)*co(1,node2)-g(2)*co(2,node2)-g(3)*co(3,node2)
!
         dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!
         xflow2=xflow*xflow
!
!        determine eta (present location of jump)
!
         nelemup=prop(index+6)
         indexup=ielprop(nelemup)
         if(lakon(nelemup)(6:7).eq.'SG') then
            eta=prop(indexup+4)
            prop(indexup+4)=0.5d0
            prop(indexup+7)=0.d0
         elseif(lakon(nelemup)(6:7).eq.'WE') then
            eta=prop(indexup+4)
            prop(indexup+4)=0.5d0
            prop(indexup+7)=0.d0
         elseif(lakon(nelemup)(6:7).eq.'DS') then
            eta=prop(indexup+7)
            prop(indexup+7)=0.5d0
            prop(indexup+9)=0.d0
         endif
!     
!     element properties
!     
         if((lakon(nelem)(6:7).eq.'SG').or.
     &        (lakon(nelem)(6:7).eq.'SO').or.
     &        (lakon(nelem)(6:7).eq.'RE').or.
     &        (lakon(nelem)(6:7).eq.'  ').or.
     &        (lakon(nelem)(6:7).eq.'DS').or.
     &        (lakon(nelem)(6:7).eq.'DO')) then
            b=prop(index+1)
            s0=prop(index+2)
            if(s0.lt.-1.d0) then
               s0=dasin((z1-z2)/dl)
            endif
            sqrts0=dsqrt(1.d0-s0*s0)
            if(lakon(nelem)(6:7).ne.'SG') then
               dl=prop(index+3)
               theta=prop(index+4)
               xks=prop(index+5)
               if(dl.le.0.d0) then
                  dl=dsqrt((co(1,node2)-co(1,node1))**2+
     &                 (co(2,node2)-co(2,node1))**2+
     &                 (co(3,node2)-co(3,node1))**2)
               endif
            else
               theta=0.d0
            endif
         elseif(lakon(nelem)(6:7).eq.'WE') then
            b=prop(index+1)
            p=prop(index+2)
            c=prop(index+3)
         elseif((lakon(nelem)(6:7).eq.'CO').or.
     &          (lakon(nelem)(6:7).eq.'EL')) then
            b1=prop(index+1)
            s0=prop(index+2)
            if(s0.lt.-1.d0) then
               s0=dasin((z1-z2)/dl)
            endif
            sqrts0=dsqrt(1.d0-s0*s0)
            b2=prop(index+4)
         elseif((lakon(nelem)(6:7).eq.'DR').or.
     &          (lakon(nelem)(6:7).eq.'ST'))then
            b=prop(index+1)
            s0=prop(index+2)
            if(s0.lt.-1.d0) then
               s0=dasin((z1-z2)/dl)
            endif
            sqrts0=dsqrt(1.d0-s0*s0)
            d=prop(index+4)
         endif
!
!        contraction, enlargement, drop and step:
!        adjust h1 or h2 by solving the appropriate
!        momentum equation
!
         if((lakon(nelem)(6:7).eq.'CO').or.
     &          (lakon(nelem)(6:7).eq.'EL').or.
     &          (lakon(nelem)(6:7).eq.'DR').or.
     &          (lakon(nelem)(6:7).eq.'ST'))then
            c2=rho*rho*dg*sqrts0
!
            if(eta.gt.1.d0) then
!
!              h1 is given, h2 is unknown
!
               if(lakon(nelem)(6:7).eq.'CO') then
                  e3=b1*h1*c2*b1*b2
                  e0=2.d0*b1*h1*xflow2/e3
                  e1=-(2.d0*xflow2+c2*b1*b1*h1**3)*b2/e3
                  e2=0.d0
               elseif(lakon(nelem)(6:7).eq.'EL') then
                  e3=b1*h1*c2*b2*b2
                  e0=2.d0*b1*h1*xflow2/e3
                  e1=-(2.d0*xflow2+c2*b1*b2*h1**3)*b2/e3
                  e2=0.d0
               elseif(lakon(nelem)(6:7).eq.'DR') then
                  e3=h1*c2*b*b
                  e0=h1*2.d0*xflow2/e3
                  e1=-(2.d0*xflow2+c2*b*b*h1*(h1+d)**2)/e3
                  e2=0.d0
               elseif(lakon(nelem)(6:7).eq.'ST') then
                  e3=h1*c2*b*b
                  e0=h1*2.d0*xflow2/e3
                  e1=(h1*c2*b*b*d*d-(2.d0*xflow2+c2*b*b*h1**3))/e3
                  e2=h1*c2*b*b*2.d0*d/e3
               endif
!
!              solve the cubic equation
!
               call cubic(e0,e1,e2,solreal,solimag,nsol)
!
!              determine the real solution closest to h1
!               
               dist=1.d30
               do i=1,nsol
                  if(dabs(solreal(i)-h1).lt.dist) then
                     dist=dabs(solreal(i)-h1)
                     h2=solreal(i)
                  endif
               enddo
               if(nactdog(2,node2).ne.0) v(2,node2)=h2
            elseif(eta.lt.0.d0) then
!
!              h2 is given, h1 is unknown
!
               if(lakon(nelem)(6:7).eq.'CO') then
                  e3=c2*b1*b1*b2*h2
                  e0=2.d0*xflow2*b2*h2/e3
                  e1=-b1*(2.d0*xflow2+c2*b1*b2*h2**3)/e3
                  e2=0.d0
               elseif(lakon(nelem)(6:7).eq.'EL') then
                  e3=c2*b1*b2*b2*h2
                  e0=2.d0*xflow2*b2*h2/e3
                  e1=-b1*(2.d0*xflow2+c2*b2*b2*h2**3)/e3
                  e2=0.d0
               elseif(lakon(nelem)(6:7).eq.'DR') then
                  e3=c2*b*b*h2
                  e0=2.d0*xflow2*h2/e3
                  e1=(c2*b*b*d*d*h2-(2.d0*xflow2+c2*b*b*h2**3))/e3
                  e2=c2*b*b*2.d0*d*h2/e3
               elseif(lakon(nelem)(6:7).eq.'ST') then
                  e3=c2*b*b*h2
                  e0=2.d0*xflow2*h2/e3
                  e1=-(2.d0*xflow2+c2*b*b*h2*(h2+d)**2)/e3
                  e2=0.d0
               endif   
!
!              solve the cubic equation
!
               call cubic(e0,e1,e2,solreal,solimag,nsol)
c               write(30,*) 'check ',solreal(1)**3+e1*solreal(1)+e0
!
c               write(30,*) 'nsol',nsol
c               write(30,*) 'solreal',(solreal(i),i=1,3)
c               write(30,*) 'solimag',(solimag(i),i=1,3)
!
!              determine the real solution closest to h2
!               
               dist=1.d30
               do i=1,nsol
                  if(dabs(solreal(i)-h2).lt.dist) then
                     dist=dabs(solreal(i)-h2)
                     h1=solreal(i)
                  endif
               enddo
               if(nactdog(2,node1).ne.0) v(2,node1)=h1
            endif
            return
         endif
!
         if(xks.gt.0.d0) then
!     
!     White-Coolebrook
!     
!     hydraulic diameter
!     
            d=2.d0*(h1+h2)
            reynolds=4.d0*xflow/(b*dvi)
            form_fact=1.d0
            call friction_coefficient(dl,d,xks,reynolds,form_fact,
     &           friction)
         endif
!     
!     geometric data
!     
         cth=dcos(theta)
         tth=dtan(theta)
!     
!     nonprismatic cross section
!     
         if(lakon(nelem)(6:7).eq.'  ') then
            dbds=prop(index+7)
         else
            dbds=0.d0
         endif
!     
!     width at water surface
!     
         dD1dh1=2.d0*tth
         dD2dh2=dD1dh1
         D1=b+h1*dD1dh1
         D2=b+dl*dbds+h2*dD2dh2
!     
!     cross section
!     
         A1=h1*(b+h1*tth)
         A2=h2*(b+dl*dbds+h2*tth)
!     
!     perimeter
!     
         P1=b+2.d0*h1/cth
         P2=b+dl*dbds+2.d0*h2/cth
!     
!     factor for friction
!     
         if(xks.gt.0.d0) then
!     White-Coolebrook
            um1=friction/8.d0
            um2=um1
         else
!     Manning
            um1=xks*xks*dg*(P1/A1)**(1.d0/3.d0)
            um2=xks*xks*dg*(P2/A2)**(1.d0/3.d0)
         endif
!     
!     constants
!     
         c1=rho*rho*dg
         c2=c1*sqrts0
         c1=c1*s0
!
         if(eta.gt.1.d0) then
            xt1=c1*A1**3+(h1*dbds-um1*P1)*xflow2
            xn1=c2*A2**3-D2*xflow2
            if(nactdog(2,node2).ne.0) v(2,node2)=h1+dl*xt1/xn1
c            write(30,*) 'move jump:  h1 h2,h2new ',h1,h2,v(2,node2)
         elseif(eta.lt.0.d0) then
            xt2=c1*A2**3+(h2*dbds-um2*P2)*xflow2
            xn2=c2*A2**3-D2*xflow2
            if(nactdog(2,node1).ne.0) 
     &           v(2,node1)=h2-dl*xt2/xn2
c            write(30,*) 'move jump: h1 h1new h2 ',h1,v(2,node1),h2
         endif
      endif
!     
      return
      end
      

