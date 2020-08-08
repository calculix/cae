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
!     calculate the initial conditions for the gas 
!         - the initial pressure
!         - identifying the chambers and gas pipe nodes
!           for gas networks
!         - the initial flow
!         - calculating the static temperature for gas networks
!     
      subroutine initialchannel(itg,ieg,ntg,ac,bc,lakon,v,
     &     ipkon,kon,nflow,ikboun,nboun,prop,ielprop,
     &     nactdog,ndirboun,nodeboun,xbounact,
     &     ielmat,ntmat_,shcon,nshcon,physcon,ipiv,nteq,
     &     rhcon,nrhcon,ipobody,ibody,xbodyact,co,nbody,network,
     &     iin_abs,vold,set,istep,iit,mi,ineighe,ilboun,ttime,
     &     time,iaxial)
!     
      implicit none
!     
      logical identity,calcinitialpressure,gravity,gaspipe
!
      character*8 lakon(*)
      character*81 set(*)
!           
      integer mi(*),ieg(*),nflow,i,j,ntg,ielmat(mi(3),*),ntmat_,id,
     &     node1,node2,
     &     nelem,index,nshcon(*),ipkon(*),kon(*),ikboun(*),nboun,idof,
     &     nodem,idirf(8),nactdog(0:3,*),imat,ielprop(*),id1,id2,
     &     nodef(8),ndirboun(*),nodeboun(*),itg(*),node,kflag,ipiv(*),
     &     nrhs,info,idof1,idof2,nteq,nrhcon(*),ipobody(2,*),ibody(3,*),
     &     nbody,numf,network,iin_abs,icase,index2,index1,nelem1,nelem2,
     &     node11,node21,node12,node22,istep,iit,ineighe(*),
     &     ilboun(*),nelemup,k,node2up,ider,iaxial,iplausi
!     
      real*8 ac(nteq,nteq), bc(nteq),prop(*),shcon(0:3,ntmat_,*),
     &     f,df(8),xflow,xbounact(*),v(0:mi(2),*),cp,r,tg1,
     &     tg2,gastemp,physcon(*),pressmin,dvi,rho,g(3),z1,z2,
     &     rhcon(0:1,ntmat_,*),co(3,*),xbodyact(7,*),kappa,
     &     a,Tt,Ts,pressmax,constant,vold(0:mi(2),*),href,
     &     ttime,time
!
      kflag=1
!
!     applying the boundary conditions
!
      do j=1,nboun
         v(ndirboun(j),nodeboun(j))=xbounact(j)
      enddo
!     
!     determining the initial pressure (not for purely thermal networks)
!     and identifying the chamber and gas pipe nodes (only for gas
!     networks)
!
      if(network.gt.2) then
!   
!        determining whether pressure initial conditions 
!        are provided for all nodes
!     
         pressmin=-1.d0
         pressmax=0.d0
         constant=1.55d0

         do i=1,ntg
            node=itg(i)
            if(v(2,node).lt.1.d-10) then
               v(2,node)=0.d0
            else
               if(pressmin.lt.0.d0) then
                  pressmin=v(2,node)
               elseif(v(2,node).lt.pressmin) then
                  pressmin=v(2,node)
               endif
!
               if(v(2,node).gt.pressmax)then
                  pressmax=v(2,node)
               endif
!
            endif
         enddo
!         
         if(pressmin.lt.0.d0) then
            write(*,*) 
     &        '*ERROR in initialchannel: minimum initial pressure'
            write(*,*) '       is smaller than zero'
            call exit(201)
         endif
!
!        in nodes in which no initial pressure is given v(2,*)
!        is replaced by -n, where n is the number of elements the
!        node belongs to: allows to find boundary nodes of the 
!        network
!
         gaspipe=.false.
         calcinitialpressure=.false.
!
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            call nident(itg,node1,ntg,id1)
            call nident(itg,node2,ntg,id2)
!     
            if (((lakon(nelem)(1:5).eq.'DGAPF').and.(iin_abs.eq.0))
     &           .or.((lakon(nelem)(1:3).eq.'DRE')
     &           .and.(lakon(nelem)(1:7).ne.'DREWAOR')
     &           .and.(iin_abs.eq.0))) then 
!     
!     In the case of a element of type GASPIPE or RESTRICTOR 
!     (except TYPE= RESTRICTOR WALL ORIFICE)
!     the number of pipes connected to node 1 and 2
!     are computed and stored in ineighe(id1)
!     respectively ineighe(id2)
!     
               gaspipe=.true.
               if(node1.ne.0) then
                  if (ineighe(id1).ge.0) then
!
                     if(node2.ne.0)then
                        ineighe(id1)=ineighe(id1)+1
                     endif
                  endif
               endif
               if(node2.ne.0) then
                  if (ineighe(id2).ge.0) then
                     if(node1.ne.0) then
                        ineighe(id2)=ineighe(id2)+1
                     endif
                  endif
               endif
            else
               if(iin_abs.eq.0) then
!     
!     for all other elements (different from GASPIPE or 
!     RESTRICTOR), including RESTRICTOR WALL ORIFICE 
!     ineighe(idi)=-1
!     which means that they are connected to chambers
!     i.e. static and total values are equal
!     
                  if (node1.ne.0) then
                     ineighe(id1)=-1
                  endif
                  if(node2.ne.0) then
                     ineighe(id2)=-1
                  endif
               endif
            endif
!     
            if((node1.eq.0).or.(node2.eq.0)) cycle
            if(v(2,node1).lt.1.d-10) then
               v(2,node1)=v(2,node1)-1.d0
               calcinitialpressure=.true.
            endif
            if(v(2,node2).lt.1.d-10) then
               v(2,node2)=v(2,node2)-1.d0
               calcinitialpressure=.true.
            endif
         enddo
!
!        for each end node i: if ineighe(i)<0: chamber
!                         else: ineighe(i)=number of pipe connections
!
      else
!
!       identifying the chamber nodes for purely thermal
!       gas networks (needed to determine the static
!       temperature which is used for the material properties)
!         
         do i=1,nflow
            nelem=ieg(i)
            if((lakon(nelem)(2:3).eq.'LP').or.
     &         (lakon(nelem)(2:3).eq.'LI')) cycle
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            if(node1.ne.0) then
               call nident(itg,node1,ntg,id1)
               ineighe(id1)=-1
            endif
            if(node2.ne.0) then
               call nident(itg,node2,ntg,id2)
               ineighe(id2)=-1
            endif
         enddo
      endif
!
!     temperature initial conditions
!
      do i=1,ntg
         node=itg(i)
         if (nactdog(0,node).eq.0) cycle
         if (v(0,node)-physcon(1).lt.1.d-10) then
            write(*,*)
     &           '*WARNING in initialchannel : the initial temperature f
     &or node',node
            write(*,*) 
     &           'is O Kelvin or less; the default is taken (293 K)'
            write(*,*)
            v(0,node)=293.d0-physcon(1)
         endif
      enddo
!    
!     initialisation of bc
!     
      do i=1,nteq
         bc(i)=0.d0
      enddo
!  
!     determining the initial mass flow in those nodes for which no
!     flux boundary conditions are defined
!     liquid channels are treated separately
!   
      if(network.gt.2) then
!     
!     calculate the initial mass flow
!     
!     check whether the mass flow is given as a boundary condition
!     
         do j=1,nflow
            nelem=ieg(j)
            index=ipkon(nelem)
            nodem=kon(index+2)
            if(nactdog(1,nodem).eq.0) then
               idof=8*(nodem-1)+1
               call nident(ikboun,idof,nboun,id)
               if(id.gt.0) then
                  if(ikboun(id).eq.idof) then
                     xflow=xbounact(ilboun(id))
                     if(dabs(xflow).gt.1.d-30) exit
                  endif
               endif
            endif
         enddo
!     
         if(dabs(xflow).gt.1.d-30) then
!     
!     if nonzero: set all mass flow to this value
!     
            do j=1,nflow
               nelem=ieg(j)
               index=ipkon(nelem)
               nodem=kon(index+2)
               if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow
            enddo
         else
!     
!     calculate the mass flow: look for a sluice gate or weir        
!     
            do j=1,nflow
               nelem=ieg(j)
               if((lakon(nelem)(6:7).ne.'SG').and.
     &              (lakon(nelem)(6:7).ne.'WE')) cycle
               index=ipkon(nelem)
               node1=kon(index+1)
               node2=kon(index+3)
               if((node1.eq.0).or.(node2.eq.0)) cycle
               nodem=kon(index+2)
!     
!     determine the gravity vector
!     
               gravity=.false.
               do k=1,3
                  g(k)=0.d0
               enddo
               if(nbody.gt.0) then
                  index=nelem
                  do
                     k=ipobody(1,index)
                     if(k.eq.0) exit
                     if(ibody(1,k).eq.2) then
                        g(1)=g(1)+xbodyact(1,k)*xbodyact(2,k)
                        g(2)=g(2)+xbodyact(1,k)*xbodyact(3,k)
                        g(3)=g(3)+xbodyact(1,k)*xbodyact(4,k)
                        gravity=.true.
                     endif
                     index=ipobody(2,index)
                     if(index.eq.0) exit
                  enddo
               endif
               if(.not.gravity) then
                  write(*,*)
     &              '*ERROR in initialchannel: no gravity vector'
                  write(*,*) '       was defined for liquid element',
     &                 nelem
                  call exit(201)
               endif
!     
               tg1=v(0,node1)
               tg2=v(0,node2)
               gastemp=(tg1+tg2)/2.d0
               imat=ielmat(1,nelem)
               call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,
     &              cp,r,dvi,rhcon,nrhcon,rho)
!     
               call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &              nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &              nodef,idirf,df,cp,r,rho,physcon,g,co,dvi,numf,
     &              vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &              ttime,time,iaxial,iplausi)
!     
               if(dabs(xflow).gt.1.d-30) exit
            enddo
!     
            if(dabs(xflow).gt.1.d-30) then
!     
!     if nonzero: set all mass flow to this value
!     
               do j=1,nflow
                  nelem=ieg(j)
                  index=ipkon(nelem)
                  nodem=kon(index+2)
                  if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow
               enddo
            else
               write(*,*) '*ERROR in initialchannel: initial mass flow'
               write(*,*) '       cannot be determined'
               call exit(201)
            endif
         endif
!     
!     calculate the depth
!     
         if(calcinitialpressure) then
!     
!     determine the streamdown depth for sluice gates,
!     weirs and discontinuous slopes
!     
            do j=1,nflow
               nelem=ieg(j)
               if((lakon(nelem)(6:7).ne.'SG').and.
     &              (lakon(nelem)(6:7).ne.'WE').and.
     &              (lakon(nelem)(6:7).ne.'DS')) cycle
               index=ipkon(nelem)
               node1=kon(index+1)
               node2=kon(index+3)
               if((node1.eq.0).or.(node2.eq.0)) cycle
               nodem=kon(index+2)
!     
!     determine the gravity vector
!     
               gravity=.false.
               do k=1,3
                  g(k)=0.d0
               enddo
               if(nbody.gt.0) then
                  index=nelem
                  do
                     k=ipobody(1,index)
                     if(k.eq.0) exit
                     if(ibody(1,k).eq.2) then
                        g(1)=g(1)+xbodyact(1,k)*xbodyact(2,k)
                        g(2)=g(2)+xbodyact(1,k)*xbodyact(3,k)
                        g(3)=g(3)+xbodyact(1,k)*xbodyact(4,k)
                        gravity=.true.
                     endif
                     index=ipobody(2,index)
                     if(index.eq.0) exit
                  enddo
               endif
               if(.not.gravity) then
                  write(*,*)
     &              '*ERROR in initialchannel: no gravity vector'
                  write(*,*) '       was defined for liquid element',
     &                 nelem
                  call exit(201)
               endif
!     
               tg1=v(0,node1)
               tg2=v(0,node2)
               gastemp=(tg1+tg2)/2.d0
               imat=ielmat(1,nelem)
               call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,
     &              cp,r,dvi,rhcon,nrhcon,rho)
!     
               call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &              nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &              nodef,idirf,df,cp,r,rho,physcon,g,co,dvi,numf,
     &              vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &              ttime,time,iaxial,iplausi)
!     
            enddo
!     
!     for all other elements the depth is taken to be
!     0.9 of the depth in the downstream node of the
!     streamup reference element
!     
            do j=1,nflow
               nelem=ieg(j)
               if((lakon(nelem)(6:7).eq.'SG').or.
     &              (lakon(nelem)(6:7).eq.'WE').or.
     &              (lakon(nelem)(6:7).eq.'DS')) cycle
!     
               index=ipkon(nelem)
               node1=kon(index+1)
               node2=kon(index+3)
               if((node1.eq.0).or.(node2.eq.0)) cycle
!     
               index=ielprop(nelem)
               nelemup=nint(prop(index+6))
               node2up=kon(ipkon(nelemup)+3)
               href=0.9d0*v(2,node2up)
               if(nactdog(2,node1).ne.0) 
     &              v(2,node1)=href
               if(nactdog(2,node2).ne.0) 
     &              v(2,node2)=href
            enddo
!     
!     reapplying the boundary conditions (the depth of the
!     sluice gate may have changed if it exceeded the critical
!     value
!     
            do j=1,nboun
               v(ndirboun(j),nodeboun(j))=xbounact(j)
            enddo
         endif
!     
!     calculating the static temperature for nodes belonging to gas pipes
!     and restrictors (except RESTRICTOR WALL ORIFICE)
!     
      endif
!     
!     for chambers the static temperature equals the total
!     temperature
!     
      do i=1,ntg
         if(ineighe(i).eq.-1) v(3,itg(i))=v(0,itg(i))
      enddo
!     
      return
      end
      
