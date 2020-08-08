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
!     construction of the B matrix      
!
!     author: Yannick Muller
!     
      subroutine flowoutput(itg,ieg,ntg,nteq,bc,lakon,ntmat_,
     &     v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
     &     ielmat,prop,ielprop,nactdog,nacteq,iin,physcon,
     &     camt,camf,camp,rhcon,nrhcon,vold,jobnamef,set,istartset,
     &     iendset,ialset,nset,mi,iaxial,istep,iit)
!     
      implicit none
!     
      logical identity
!
      character*8 lakon(*)
      character*81 set(*)
      character*132 jobnamef(*),fnnet
!     
      integer mi(*),itg(*),ieg(*),ntg,nflow,ielmat(mi(3),*),i,
     &     nrhcon(*),node,iaxial,ider,idirf(8),ieq,imat,kflag,
     &     ntmat_,nteq,nshcon(*),nelem,index,ipkon(*),kon(*),iin,
     &     nactdog(0:3,*),nacteq(0:3,*),ielprop(*),node1,nodem,node2,
     &     istartset(*),iendset(*),ialset(*),nset,nodef(8),numf,
     &     istep,iit,iplausi
!     
      real*8 physcon(*),v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),
     &     prop(*),dtime,ttime,time,xflow,camp(*),camt(*),camf(*),
     &     rhcon(0:1,ntmat_,*),vold(0:mi(2),*),eta,
     &     bc(*),cp,dvi,df(8),gastemp,f,g(3),r,rho,ts1,ts2
!
!     
      do i=1,132
         if(jobnamef(1)(i:i).eq.' ') exit
      enddo
      i=i-1
      fnnet=jobnamef(1)(1:i)//'.net'
      open(1,file=fnnet,status='unknown')
!
      kflag=3
!
      do i=1,nflow
         nelem=ieg(i)
!     
!        output for gas networks
!
         if((lakon(nelem)(2:5).ne.'LIPI').and.
     &      (lakon(nelem)(2:5).ne.'LICH')) then
!     
            index=ipkon(nelem)
            node1=kon(index+1)
            nodem=kon(index+2)
            node2=kon(index+3)
!
            xflow=v(1,nodem)
!
            if(lakon(nelem)(2:3).ne.'LP') then
!
!              compressible
!
               if(node1.eq.0) then
                  ts1=v(3,node2)
                  ts2=ts1
               elseif(node2.eq.0) then
                  ts1=v(3,node1)
                  ts2=ts1
               else
                  ts1=v(3,node1)
                  ts2=v(3,node2)
               endif
               gastemp=(ts1+ts2)/2.d0
            else
!
!              incompressible
!
               if(xflow.gt.0) then
                  gastemp=v(3,node1)
               else
                  gastemp=v(3,node2)
               endif
            endif
!
            imat=ielmat(1,nelem)
!
            call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,
     &         dvi,rhcon,nrhcon,rho)
!
            if(nacteq(2,nodem).ne.0) then
               ieq=nacteq(2,nodem)
!
!              dummy set number
!
               numf=1
!
               call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &              nactdog,identity,
     &              ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &              cp,r,rho,physcon,g,co,dvi,numf,vold,set,shcon,
     &              nshcon,rhcon,nrhcon,ntmat_,mi,ider,ttime,time,
     &              iaxial,iplausi)
            endif
         endif
!            
         if(lakon(ieg(i))(2:5).eq.'LICH') then
            if((lakon(ieg(i))(6:7).eq.'SG').or.
     &           (lakon(ieg(i))(6:7).eq.'WE').or.
     &           (lakon(ieg(i))(6:7).eq.'DS')) then
               index=ipkon(ieg(i))
               node=kon(index+2)
               if(nactdog(3,node).eq.0) cycle
               index=ielprop(ieg(i))
               if(lakon(ieg(i))(6:7).eq.'SG') then
                  eta=prop(index+4)
                  nelem=int(prop(index+7))      
               elseif(lakon(ieg(i))(6:7).eq.'WE') then
                  eta=prop(index+4)
                  nelem=int(prop(index+7))      
               elseif(lakon(ieg(i))(6:7).eq.'DS') then
                  eta=prop(index+7)
                  nelem=int(prop(index+9))      
               endif
               if(nelem.ne.0) then
                  write(*,*) '     *INFO in flowoutput: hydraulic jump'
                  write(*,*) '           in element ',nelem,'.'
                  write(*,*) '           relative location:',eta
                  write(*,*)
               endif
            endif
         endif
      enddo
!
      close(1)
!
      return
      end
