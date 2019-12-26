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
      subroutine tempload_em(xforcold,xforc,xforcact,iamforc,nforc,&
        xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,&
        xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,&
        amta,namta,nam,ampli,time,reltime,ttime,dtime,ithermal,nmethod,&
        xbounold,xboun,xbounact,iamboun,nboun,&
        nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,&
        co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,&
        ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,&
        iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,&
        h0scale,inomat,ipobody,iponoel,inoel,ipkon,kon,lakon,&
        ielprop,prop,ielmat,shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,&
        ncocon)
      !
      !     calculates the loading at a given time
      !
      implicit none
      !
      logical gasnode
      !
      character*1 entity
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 amname(*)
      character*81 tieset(3,*)
      !
      integer iamforc(*),iamload(2,*),iamt1(*),nelemload(2,*),&
        nam,i,istart,iend,id,nforc,nload,nk,namta(3,*),ithermal,&
        nmethod,iamt1i,iamboun(*),nboun,iamforci,iambouni,&
        iamloadi1,iamloadi2,ibody(3,*),itg(*),ntg,idof,one,&
        nbody,iambodyi,nodeboun(*),ndirboun(*),nodeforc(2,*),&
        ndirforc(*),istep,iinc,msecpt,node,j,ikboun(*),ilboun(*),&
        ipresboun,mi(*),ntrans,inotr(2,*),idummy,integerglob(*),&
        istartset(*),iendset(*),ialset(*),ntie,iselect(1),&
        nmpc,ikmpc(*),ilmpc(*),nodempc(3,*),k,ist,index,ipompc(*),&
        inomat(*),ipobody(2,*),iponoel(*),inoel(2,*),ipkon(*),&
        kon(*),ielprop(*),ielmat(mi(3),*),nshcon(*),nrhcon(*),&
        ncocon(2,*),ntmat_
      !
      real*8 xforc(*),xforcact(*),xload(2,*),xloadact(2,*),&
        t1(*),t1act(*),amta(2,*),ampli(*),time,fixed_temp,&
        xforcold(*),xloadold(2,*),t1old(*),reltime,coefmpc(*),&
        xbounold(*),xboun(*),xbounact(*),ttime,dtime,reftime,&
        xbody(7,*),xbodyold(7,*),xbodyact(7,*),co(3,*),&
        vold(0:mi(2),*),abqtime(2),coords(3),trab(7,*),&
        veold(0:mi(2),*),ddummy,doubleglob(*),h0scale,prop(*),&
        shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),cocon(0:6,ntmat_,*)
      !
      intent(in) xforcold,xforc,iamforc,nforc,&
        xloadold,xload,iamload,nload,ibody,xbody,nbody,&
        xbodyold,t1old,t1,iamt1,nk,ipkon,kon,lakon,&
        amta,namta,nam,time,reltime,ttime,dtime,ithermal,nmethod,&
        xbounold,xboun,iamboun,nboun,ielprop,prop,&
        nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,ielmat,&
        co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,&
        ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,&
        iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,&
        inomat,ipobody,iponoel,inoel,shcon,nshcon,rhcon,&
        nrhcon,ntmat_,cocon,ncocon
      !
      intent(inout) xloadact,xforcact,xbounact,h0scale,ampli,t1act,&
        xbodyact
      !
      data msecpt /1/
      data one /1/
      !
      h0scale=1.d0
      !
      !     if an amplitude is active, the loading is scaled according to
      !     the actual time. If no amplitude is active, then the load is
      !     - scaled according to the relative time for a static step
      !     - applied as a step loading for a dynamic step
      !
      !     calculating all amplitude values for the current time
      !
      do i=1,nam
         if(namta(3,i).lt.0) then
            reftime=ttime+time
         else
            reftime=time
         endif
         if(abs(namta(3,i)).ne.i) then
            reftime=reftime-amta(1,namta(1,i))
            istart=namta(1,abs(namta(3,i)))
            iend=namta(2,abs(namta(3,i)))
            if(istart.eq.0) then
               call uamplitude(reftime,amname(namta(3,i)),ampli(i))
               cycle
            endif
         else
            istart=namta(1,i)
            iend=namta(2,i)
            if(istart.eq.0) then
               call uamplitude(reftime,amname(i),ampli(i))
               cycle
            endif
         endif
         call identamta(amta,reftime,istart,iend,id)
         if(id.lt.istart) then
            ampli(i)=amta(2,istart)
         elseif(id.eq.iend) then
            ampli(i)=amta(2,iend)
         else
            ampli(i)=amta(2,id)+(amta(2,id+1)-amta(2,id))&
                 *(reftime-amta(1,id))/(amta(1,id+1)-amta(1,id))
         endif
      enddo
      !
      !     scaling the boundary conditions
      !
      do i=1,nboun
         if((xboun(i).lt.1.2357111318d0).and.&
              (xboun(i).gt.1.2357111316d0)) then
            !
            !           user subroutine for boundary conditions
            !
            node=nodeboun(i)
            !
            !           check whether node is a gasnode
            !
            gasnode=.false.
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
               if(itg(id).eq.node) then
                  gasnode=.true.
               endif
            endif
            !
            abqtime(1)=time
            abqtime(2)=ttime+time
            !
            !           a gasnode cannot move (displacement DOFs are used
            !           for other purposes, e.g. mass flow and pressure)
            !
            if(gasnode) then
               do j=1,3
                  coords(j)=co(j,node)
               enddo
            else
               do j=1,3
                  coords(j)=co(j,node)+vold(j,node)
               enddo
            endif
            !
            if(ndirboun(i).eq.0) then
               call utemp(xbounact(i),msecpt,istep,iinc,abqtime,node,&
                    coords,vold,mi,iponoel,inoel,&
                    ipobody,xbodyact,ibody,ipkon,kon,&
                    lakon,ielprop,prop,ielmat,&
                    shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
            else
               call uboun(xbounact(i),istep,iinc,abqtime,node,&
                    ndirboun(i),coords,vold,mi,iponoel,inoel,&
                    ipobody,xbodyact,ibody,ipkon,kon,&
                    lakon,ielprop,prop,ielmat,&
                    shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
            endif
            cycle
         endif
         if((xboun(i).lt.1.9232931375d0).and.&
              (xboun(i).gt.1.9232931373d0)) then
            !
            !           boundary conditions for submodel
            !
            node=nodeboun(i)
            !
            !           check whether node is a gasnode
            !
            gasnode=.false.
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
               if(itg(id).eq.node) then
                  gasnode=.true.
               endif
            endif
            !
            !           for the interpolation of submodels the undeformed
            !           geometry is taken
            !
            do j=1,3
               coords(j)=co(j,node)
            enddo
            !
            entity='N'
            one=1
            iselect(1)=ndirboun(i)+1
            call interpolsubmodel(integerglob,doubleglob,xbounact(i),&
                 coords,iselect,one,node,tieset,istartset,iendset,&
                 ialset,ntie,entity)
            !             write(*,*) 'tempload ',node,ndirboun(i),xbounact(i)
            !
            if(nmethod.eq.1) then
               xbounact(i)=xbounold(i)+&
                    (xbounact(i)-xbounold(i))*reltime
            endif
            cycle
         endif
         !
         if(nam.gt.0) then
            iambouni=iamboun(i)
         else
            iambouni=0
         endif
         if(iambouni.gt.0) then
            xbounact(i)=xboun(i)*ampli(iambouni)
            !
            !           scaling of the voltage across a coil in an
            !           electromagnetic calculation
            !
            node=nodeboun(i)
            if(inomat(node).eq.0) h0scale=ampli(iambouni)
         !
         elseif(nmethod.eq.1) then
            xbounact(i)=xbounold(i)+&
               (xboun(i)-xbounold(i))*reltime
         else
            xbounact(i)=xboun(i)
         endif
      enddo
      !
      !     scaling the loading
      !
      do i=1,nforc
         if(ndirforc(i).eq.0) then
            if((xforc(i).lt.1.2357111318d0).and.&
               (xforc(i).gt.1.2357111316d0)) then
               !
               !              user subroutine for the concentrated heat flux
               !
               node=nodeforc(1,i)
               !
               !              check whether node is a gasnode
               !
               gasnode=.false.
               call nident(itg,node,ntg,id)
               if(id.gt.0) then
                  if(itg(id).eq.node) then
                     gasnode=.true.
                  endif
               endif
               !
               abqtime(1)=time
               abqtime(2)=ttime+time
               !
               !              a gasnode cannot move (displacement DOFs are used
               !              for other purposes, e.g. mass flow and pressure)
               !
               if(gasnode) then
                  do j=1,3
                     coords(j)=co(j,node)
                  enddo
               else
                  do j=1,3
                     coords(j)=co(j,node)+vold(j,node)
                  enddo
               endif
               !
               call cflux(xforcact(i),msecpt,istep,iinc,abqtime,node,&
                    coords,vold,mi,ipkon,kon,lakon,iponoel,inoel,&
                    ielprop,prop,ielmat,shcon,nshcon,rhcon,nrhcon,&
                    ntmat_,cocon,ncocon)
               cycle
            endif
         else
            if((xforc(i).lt.1.2357111318d0).and.&
               (xforc(i).gt.1.2357111316d0)) then
               !
               !              user subroutine for the concentrated load
               !
               node=nodeforc(1,i)
               !
               abqtime(1)=time
               abqtime(2)=ttime+time
               !
               do j=1,3
                  coords(j)=co(j,node)+vold(j,node)
               enddo
               !
               call cload(xforcact(i),istep,iinc,abqtime,node,&
                    ndirforc(i),coords,vold,mi,ntrans,trab,inotr,veold)
               cycle
            elseif((xforc(i).lt.1.9232931375d0).and.&
                   (xforc(i).gt.1.9232931373d0)) then
               !
               !     boundary conditions for submodel
               !
               node=nodeforc(1,i)
               !
               !     for the interpolation of submodels the undeformed
               !     geometry is taken
               !
               do j=1,3
                  coords(j)=co(j,node)
               enddo
               !
               entity='N'
               one=1
               iselect(1)=ndirforc(i)+10
               call interpolsubmodel(integerglob,doubleglob,xforcact(i),&
                    coords,iselect,one,node,tieset,istartset,iendset,&
                    ialset,ntie,entity)
               !
               if(nmethod.eq.1) then
                  xforcact(i)=xforcold(i)+&
                       (xforcact(i)-xforcold(i))*reltime
               endif
               cycle
            endif
         endif
         !
         if(nam.gt.0) then
            iamforci=iamforc(i)
         else
            iamforci=0
         endif
         if(iamforci.gt.0) then
            if(amname(iamforci).eq.'RAMP12357111317') then
               xforcact(i)=xforcold(i)+&
                    (xforc(i)-xforcold(i))*reltime
            else
               xforcact(i)=xforc(i)*ampli(iamforci)
            endif
         elseif(nmethod.eq.1) then
            xforcact(i)=xforcold(i)+&
               (xforc(i)-xforcold(i))*reltime
         else
            xforcact(i)=xforc(i)
         endif
      enddo
      !
      do i=1,nload
         ipresboun=0
         !
         !        check for pressure boundary conditions
         !
         if(sideload(i)(3:4).eq.'NP') then
            node=nelemload(2,i)
            idof=8*(node-1)+2
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.idof) then
                  ipresboun=1
                  xloadact(1,i)=xbounact(ilboun(id))
               endif
            endif
         endif
         !
         if(ipresboun.eq.0) then
            if(nam.gt.0) then
               iamloadi1=iamload(1,i)
               iamloadi2=iamload(2,i)
            else
               iamloadi1=0
               iamloadi2=0
            endif
            if(iamloadi1.gt.0) then
               if(amname(iamloadi1).eq.'RAMP12357111317') then
                  xloadact(1,i)=xloadold(1,i)+&
                       (xload(1,i)-xloadold(1,i))*reltime
               else
                  xloadact(1,i)=xload(1,i)*ampli(iamloadi1)
               endif
            elseif(nmethod.eq.1) then
               xloadact(1,i)=xloadold(1,i)+&
                    (xload(1,i)-xloadold(1,i))*reltime
            else
               xloadact(1,i)=xload(1,i)
            endif
            if(iamloadi2.gt.0) then
               xloadact(2,i)=xload(2,i)*ampli(iamloadi2)
            !             elseif(nmethod.eq.1) then
            !                xloadact(2,i)=xload(2,i)
            else
               xloadact(2,i)=xload(2,i)
            endif
         endif
      enddo
      !
      do i=1,nbody
         if(nam.gt.0) then
            iambodyi=ibody(2,i)
         else
            iambodyi=0
         endif
         if(iambodyi.gt.0) then
            if(amname(iambodyi).eq.'RAMP12357111317') then
               xbodyact(1,i)=xbodyold(1,i)+&
                    (xbody(1,i)-xbodyold(1,i))*reltime
            else
               xbodyact(1,i)=xbody(1,i)*ampli(iambodyi)
            endif
         elseif(nmethod.eq.1) then
            xbodyact(1,i)=xbodyold(1,i)+&
                 (xbody(1,i)-xbodyold(1,i))*reltime
         else
            xbodyact(1,i)=xbody(1,i)
         endif
      enddo
      !
      !     scaling the temperatures
      !
      if(ithermal.eq.1) then
         do i=1,nk
            if((t1(i).lt.1.2357111318d0).and.&
                 (t1(i).gt.1.2357111316d0)) then
               !
               abqtime(1)=time
               abqtime(2)=ttime+time
               !
               do j=1,3
                  coords(j)=co(j,i)+vold(j,i)
               enddo
               call utemp(t1act(i),msecpt,istep,iinc,abqtime,i,&
                    coords,vold,mi,iponoel,inoel,&
                    ipobody,xbodyact,ibody,ipkon,kon,&
                    lakon,ielprop,prop,ielmat,&
                    shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
               cycle
            endif
            !
            if((t1(i).lt.1.9232931375d0).and.&
                 (t1(i).gt.1.9232931373d0)) then
               !
               !              for the interpolation of submodels the undeformed
               !              geometry is taken
               !
               do j=1,3
                  coords(j)=co(j,i)
               enddo
               !
               entity='N'
               one=1
               iselect(1)=1
               call interpolsubmodel(integerglob,doubleglob,t1act(i),&
                    coords,iselect,one,i,tieset,istartset,iendset,&
                    ialset,ntie,entity)
               !
               if(nmethod.eq.1) then
                  t1act(i)=t1old(i)+(t1act(i)-t1old(i))*reltime
               endif
               cycle
            endif
            !
            if(nam.gt.0) then
               iamt1i=iamt1(i)
            else
               iamt1i=0
            endif
            if(iamt1i.gt.0) then
               if(amname(iamt1i).eq.'RAMP12357111317') then
                  t1act(i)=t1old(i)+(t1(i)-t1old(i))*reltime
               else
                  t1act(i)=t1(i)*ampli(iamt1i)
               endif
            elseif(nmethod.eq.1) then
               t1act(i)=t1old(i)+(t1(i)-t1old(i))*reltime
            else
               t1act(i)=t1(i)
            endif
         enddo
         !
         !        taking temperature MPC's into account
         !
         do j=1,nmpc
            k=mod(ikmpc(j),8)
            if(k.ne.0) cycle
            i=ilmpc(j)
            ist=ipompc(i)
            node=nodempc(1,ist)
            index=nodempc(3,ist)
            fixed_temp=0.d0
            if(index.ne.0) then
               do
                  fixed_temp=fixed_temp-&
                    coefmpc(index)*t1act(nodempc(1,index))
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
            t1act(node)=fixed_temp/coefmpc(ist)
         enddo
      endif
      !       write(*,*) 'nboun'
      !       do i=1,nboun
      !          write(*,'(i7,1x,e11.4,1x,e11.4)') i,xbounact(i),xboun(i)
      !       enddo
      !       write(*,*) 'nforc'
      !       do i=1,nforc
      !          write(*,'(i7,1x,e11.4,1x,e11.4)') i,xforcact(i),xforc(i)
      !       enddo
      !       write(*,*) 'nload'
      !       do i=1,nload
      !          write(*,'(i7,1x,e11.4,1x,e11.4)') i,xloadact(1,i),xload(1,i)
      !       enddo
      !
      return
      end
