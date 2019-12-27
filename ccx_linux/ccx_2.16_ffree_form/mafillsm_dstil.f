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
! >
! > \brief computation of the stiffness/mass matrix and rhs for the element with transformed basis functions
! > see mafillsm
! > Author: Saskia Sitzmann
! >
! > @param [in] co       coordinates of nodes
! > @param [in] nk       number of nodes
! > @param [in] kon      storing the connectivity list of elem. in succ. order
! > @param [in] ipkon        pointer into field kon...
! > @param [in] lakon        (i) label for element i
! > @param [in] ne       number of elements
! > @param [in] nodeboun         (i) node of SPC i
! > @param [in] ndirboun     (i) direction of SPC i
! > @param [in] xboun            (i) value of SPC i
! > @param [in] nboun            number of SPCs
! > @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
! > @param [in] nodempc          nodes and directions of MPCs
! > @param [in] coefmpc          coefficients of MPCs
! > @param [in] nmpc     number of mpcs
! > @param [in] nodeforc     point force, node
! > @param [in] ndirforc     point force, dir
! > @param [in] xforc        point force, value
! > @param [in] nforc        number of point forces
! > @param [in] nelemload    element faces of which are loaded
! > @param [in] sideload     load label
! > @param [in] xload        concentrated load
! > @param [in] nload        number of facial distributed loads
! > @param [in] xbody        body forces
! > @param [in] ipobody      pointer to xbody...
! > @param [in] nbody        number of mechanical body loads
! > @param [in] cgr      gravity forces
! > @param [out] ad      diagonal terms of stiffness matrix K
! > @param [out] au      off-diagonal terms of stiffness matrix K
! > @param [out] fext        external forces
! > @param [in] nactdof      (i,j) actual degree of freedom for direction i of node j
! > @param [in] icol     column numbers for stiffness/mass matrix
! > @param [in] jq       column pointer to field irow
! > @param [in] irow     row numbers for stiffness/mass matrix
! > @param [in] neq      number of active degrees of freedom
! > @param [out] nzl     highes column numbeer with non-zero entry
! > @param [in] nmethod      analysis method
! > @param [in] ikmpc        sorted dofs idof=8*(node-1)+dir for MPCs
! > @param [in] ilmpc        SPC numbers for sorted dofs
! > @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
! > @param [in] ilboun           SPC numbers for sorted dofs
! > @param [in] elcon        material parameters
! > @param [in] nelcon           (1,i) number of elastic constants for material i (2,i) number of temperature points
! > @param [in] rhcon        (0,j,i) temperature (1,j,i) density at the density temperature point j of material i
! > @param [in] nrhcon       (i) number of temperature data points for the density material i
! > @param [in] alcon        (0,j,i) temperature, (k,j,i) expansion coefficient k at expansion temperature point j of material i
! > @param [in] nalcon       (1,i) number of expansion constants  (2,i) number of temperature data points for expansion coefficients for material i
! > @param [in] alzero       used in material data
! > @param [in] ielmat           (j,i) material number of layer j
! > @param [in] ielorien     (j,i) orientation number of layer j
! > @param [in] norien       number of orientations
! > @param [in] orab     (7,*) description of local coordinate system
! > @param [in] ntmat_       maximum number of temperature data points for any material
! > @param [in] t0       needed for spring elements
! > @param [in] t1       needed for spring elements
! > @param [in] ithermal     >0 thermal effects are taken into account
! > @param [in] prestr       NOT USED
! > @param [in] iprestr      NOT USED
! > @param [in] vold     displacement of nodes
! > @param [in] iperturb     flag indicating what geomatric method to use
! > @param [in] sti      (1:6,k,i) 2nd order stress of element i in integration point k
! > @param [in] nzs      number of non-zero entries of mass/stiffness matrix
! > @param [in] stx      (1:6,k,i) buckling stress of element i in integration point k
! > @param [out] adb     diagonal terms of mass matrix M
! > @param [out] aub     off-diagonal terms of mass matrix M
! > @param [in] iexpl        flag indicating whether explicit (=1) or implicit (=0) method is used
! > @param [in] plicon       isotropic hardening curve
! > @param [in] nplicon          pointer to isotropic hardening curve
! > @param [in] plkcon       kinematic hardening curve
! > @param [in] nplkcon      pointer to kinematik hardening curve
! > @param [in] xstiff       (1:27,k,i) stiffness matrix of element i in integration point k for last increment
! > @param [in] npmat_       maximum number of data points for plicon
! > @param [in] dtime        real time increment size
! > @param [in] matname      (i) name of material i
! > @param [in] mi       (1) max # of integration points per element (2) max degree of freedom per element
! > @param [in] ncmat_       maximum number of elastic material constants
! > @param [in] mass     flag indicating whether to calculate the mass matrix
! > @param [in] stiffness    flag indicating whether to calculate the stiffness matrix
! > @param [in] buckling     flag indicating whether to calculate the buckling stiffness
! > @param [in] rhsi     flag indicating whether to calculate the right hand side
! > @param [in] intscheme    flag indicating what integration scheme to use
! > @param [in] physcon      used for thermal calculation
! > @param [in] shcon        used for thermal calculation, @see materialdata_th
! > @param [in] nshcon       used for thermal calculation, @see materialdata_th
! > @param [in] cocon        used for thermal calculation, @see materialdata_th
! > @param [in] ncocon       used for thermal calculation, @see materialdata_th
! > @param [in] ttime        real time size of all previous steps
! > @param [in] time     real time size of all previous increments including the present increment
! > @param [in] istep        step number
! > @param [in] iinc     increment number
! > @param [in] coriolis     flag indicating whether to calculate the coriolis matrix
! > @param [in] ibody        used for assinging body forces
! > @param [in] xloadold     magnitude of load at start of step
! > @param [in] reltime      theta+dtheta
! > @param [in] veold        current velocity
! > @param [in] springarea   used in spring elements
! > @param [in] nstate_      maximum number of state variables
! > @param [in] xstateini    state variables at start if the increment
! > @param [in] xstate       current state variables
! > @param [in] thicke       thickness matrix
! > @param [in] integerglob  used for submodel
! > @param [in] doubleglob   used for submodel
! > @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface
! > @param [in] istartset        (i) pointer to ialset containing the first set member
! > @param [in] iendset          (i) pointer to ialset containing the last set member
! > @param [in] ialset           set members
! > @param [in] ntie     number of ties
! > @param [in] nasym        flag indicating whether matrix is symmetric or not
! > @param [in] pslavsurf    field storing  position xil, etal and weight for integration point on slave side
! > @param [in] pmastsurf    field storing position and etal for integration points on master side
! > @param [in] mortar       param indicating what contact formulation is used (=0 NTS penalty, =1 GPTS penalty , >1 STS mortar)
! > @param [in] clearini     used for spring elements
! > @param [in] ielprop      used for beam elements
! > @param [in] prop     used for beam elements
! > @param [in] nslavnode    (i)pointer into field isalvnode for contact tie i
! > @param [in] islavnode    field storing the nodes of the slave surface
! > @param [in] islavsurf    islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i
! > @param [in] islavelinv       (i)==0 if there is no slave node in the element, >0 otherwise
! > @param [in] autloc       transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$
! > @param [in] irowtloc     field containing row numbers of autloc
! > @param [in] jqtloc           pointer into field irowtloc
! >
      subroutine mafillsm_dstil(co,nk,kon,ipkon,lakon,ne,nodeboun,&
        ndirboun,&
        xboun,nboun,&
        ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,&
        nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,&
        ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,&
        ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,&
        nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,&
        t0,t1,ithermal,prestr,&
        iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,&
        nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,&
        matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,&
        physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,&
        coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,&
        xstateini,xstate,thicke,integerglob,doubleglob,&
        tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,&
        mortar,clearini,ielprop,prop,ne0,fnext,nea,neb,kscale,&
        iponoel,inoel,network,&
        nslavnode,islavnode,islavsurf,islavelinv,&
        autloc, irowtloc,jqtloc)
      !
      !     filling the stiffness matrix in spare matrix format (sm)
      !
      implicit none
      !
      integer mass(2),stiffness,buckling,rhsi,stiffonly(2),coriolis
      !
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*)
      !
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),&
        nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),ikmpc(*),&
        ilmpc(*),ikboun(*),ilboun(*),mi(*),nstate_,ne0,nasym,&
        nactdof(0:mi(2),*),irow(*),icolumn,ialset(*),ielprop(*),&
        nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),ntie,&
        ielorien(mi(3),*),integerglob(*),istartset(*),iendset(*),&
        ipkon(*),intscheme,ncocon(2,*),nshcon(*),ipobody(2,*),nbody,&
        ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq(2),nzl,nmethod,&
        ithermal(2),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,&
        ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,&
        mpc1,mpc2,index1,index2,jdof,node1,node2,kflag,icalccg,&
        ntmat_,indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,&
        nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,mortar,&
        nea,neb,kscale,iponoel(*),inoel(2,*),network,ndof,&
        islavelinv(*),nslavnode(*),islavnode(*),&
        islavsurf(2,*),jqtloc(*),irowtloc(*),ii,jqtloc1(21),&
        irowtloc1(96),i1,j1,j2,konl(26)
      !
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),&
        p2(3),ad(*),au(*),bodyf(3),fext(*),xloadold(2,*),reltime,&
        t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),&
        ff(60),fnext(0:mi(2),*),&
        sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),&
        elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),springarea(2,*),&
        alcon(0:6,ntmat_,*),physcon(*),cocon(0:6,ntmat_,*),prop(*),&
        xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),&
        shcon(0:3,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),&
        plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),&
        xstiff(27,mi(1),*),veold(0:mi(2),*),om,valu2,value,dtime,ttime,&
        time,thicke(mi(3),*),doubleglob(*),clearini(3,9,*),&
        pslavsurf(3,*),pmastsurf(6,*),&
        autloc(*),autloc1(96)
      !
      intent(in) co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,&
        xboun,nboun,&
        ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,&
        nforc,nelemload,sideload,nload,xbody,ipobody,nbody,&
        nactdof,icol,jq,irow,neq,nzl,&
        ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,&
        nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,&
        t0,t1,ithermal,prestr,&
        iprestr,vold,iperturb,sti,nzs,stx,iexpl,plicon,&
        nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,&
        matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,&
        physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,&
        coriolis,ibody,xloadold,reltime,veold,nstate_,&
        xstateini,thicke,integerglob,doubleglob,&
        tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,&
        mortar,clearini,ielprop,prop,ne0,nea,neb
      !
      intent(inout) fext,ad,au,adb,aub,xload,nmethod,cgr,springarea,&
        xstate
      !
      kflag=2
      i0=0
      icalccg=0
      !       write(*,*) loc(kflag)
      !       write(*,*) loc(s)
      !       write(*,*) loc(sm)
      !       write(*,*) loc(ff)
      !       write(*,*) loc(index1)
      !
      if((stiffness.eq.1).and.(mass(1).eq.0).and.(buckling.eq.0)) then
         stiffonly(1)=1
      else
         stiffonly(1)=0
      endif
      if((stiffness.eq.1).and.(mass(2).eq.0).and.(buckling.eq.0)) then
         stiffonly(2)=1
      else
         stiffonly(2)=0
      endif
      !
      if(rhsi.eq.1) then
         !
         !        distributed forces (body forces or thermal loads or
         !        residual stresses or distributed face loads)
         !
         if((nbody.ne.0).or.(ithermal(1).ne.0).or.&
            (iprestr.ne.0).or.(nload.ne.0)) then
            idist=1
         else
            idist=0
         endif
      !
      endif
      !
      if((ithermal(1).le.1).or.(ithermal(1).eq.3)) then
      !
      !     mechanical analysis: loop over all elements
      !
      do i=nea,neb
        !
        if((ipkon(i).lt.0).or.(lakon(i)(1:1).eq.'F')) cycle
        indexe=ipkon(i)
        !      Bernhardi start
        if(lakon(i)(1:5).eq.'C3D8I') then
           nope=11
           ndof=3
        elseif(lakon(i)(4:5).eq.'20') then
           !      Bernhardi end
           nope=20
           ndof=3
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
           ndof=3
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
           ndof=3
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
           ndof=3
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
           ndof=3
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
           ndof=3
        elseif((lakon(i)(1:2).eq.'ES').and.(lakon(i)(7:7).ne.'F')) then
           !
           !          spring and contact spring elements (NO dashpot elements
           !          = ED... elements)
           !
           nope=ichar(lakon(i)(8:8))-47
           ndof=3
           !
           !          local contact spring number
           !          if friction is involved, the contact spring element
           !          matrices are determined in mafillsmas.f
           !
           if(lakon(i)(7:7).eq.'C') then
              if(nasym.eq.1) cycle
              if(mortar.eq.1) nope=kon(indexe)
           endif
        elseif(lakon(i)(1:4).eq.'MASS') then
           nope=1
           ndof=3
        elseif(lakon(i)(1:1).eq.'U') then
           !            if(lakon(i)(7:7).eq.' ') then
           !               nope=ichar(lakon(i)(8:8))-48
           !            else
           !               nope=10*(ichar(lakon(i)(7:7))-48)
           !      &             +ichar(lakon(i)(8:8))-48
           !            endif
           !            ndof=ichar(lakon(i)(6:6))-48
           ndof=ichar(lakon(i)(7:7))
           nope=ichar(lakon(i)(8:8))
        else
           cycle
        endif
        !
        !
        !     mortar start
        !
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
        !
        !     mortar end
        !
        om=0.d0
        !
        if((nbody.gt.0).and.(lakon(i)(1:1).ne.'E')) then
           !
           !          assigning centrifugal forces
           !
           bodyf(1)=0.d0
           bodyf(2)=0.d0
           bodyf(3)=0.d0
           !
           index=i
           do
              j=ipobody(1,index)
              if(j.eq.0) exit
              if(ibody(1,j).eq.1) then
                 om=xbody(1,j)
                 p1(1)=xbody(2,j)
                 p1(2)=xbody(3,j)
                 p1(3)=xbody(4,j)
                 p2(1)=xbody(5,j)
                 p2(2)=xbody(6,j)
                 p2(3)=xbody(7,j)
              !
              !          assigning gravity forces
              !
              elseif(ibody(1,j).eq.2) then
                 bodyf(1)=bodyf(1)+xbody(1,j)*xbody(2,j)
                 bodyf(2)=bodyf(2)+xbody(1,j)*xbody(3,j)
                 bodyf(3)=bodyf(3)+xbody(1,j)*xbody(4,j)
              !
              !          assigning newton gravity forces
              !
              elseif(ibody(1,j).eq.3) then
                 call newton(icalccg,ne,ipkon,lakon,kon,t0,co,rhcon,&
                      nrhcon,ntmat_,physcon,i,cgr,bodyf,ielmat,ithermal,&
                      vold,mi)
              endif
              index=ipobody(2,index)
              if(index.eq.0) exit
           enddo
        endif
        !
        !
        !     mortar start
        !
        !     generate local transformation matrix for current element
        !
        if(islavelinv(i).gt.0)then
          if(nope.eq.20 .or. nope.eq.10 .or. nope.eq.15)then
             jqtloc1(1)=1
          ii=1
          do i1=1,nope
                 node1=konl(i1)
                 do j1=jqtloc(node1),jqtloc(node1+1)-1
                    node2=irowtloc(j1)
                    do j2=1,nope
                       if(konl(j2).eq.node2)then
                          autloc1(ii)=autloc(j1)
                          irowtloc1(ii)=j2
                          ii=ii+1
                       endif
                    enddo
                 enddo
                 jqtloc1(i1+1)=ii
          enddo
           else
             jqtloc1(1)=1
             ii=1
             do i1=1,nope
               jqtloc1(i1+1)=ii
             enddo       
           endif
        endif
        !
        !      mortar end
        !
        !         WRITE(*,*) ' mafillsm_dstil: elem',i
        if(lakon(i)(1:1).ne.'U') then
           call e_c3d_dstil(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,&
                ff,i,&
                nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,&
                alzero,ielmat,ielorien,norien,orab,ntmat_,&
                t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,&
                nload,idist,sti,stx,iexpl,plicon,&
                nplicon,plkcon,nplkcon,xstiff,npmat_,&
                dtime,matname,mi(1),ncmat_,mass(1),stiffness,buckling,&
                rhsi,intscheme,ttime,time,istep,iinc,coriolis,xloadold,&
                reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,&
                springarea,nstate_,xstateini,xstate,ne0,ipkon,thicke,&
                integerglob,doubleglob,tieset,istartset,&
                iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,mortar,&
                clearini,ielprop,prop,kscale,&
                islavelinv,islavsurf,&
                autloc1,irowtloc1,jqtloc1)
        else
           call e_c3d_u(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,i,&
                nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,&
                alzero,ielmat,ielorien,norien,orab,ntmat_,&
                t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,&
                nload,idist,sti,stx,iexpl,plicon,&
                nplicon,plkcon,nplkcon,xstiff,npmat_,&
                dtime,matname,mi(1),ncmat_,mass(1),stiffness,buckling,&
                rhsi,intscheme,ttime,time,istep,iinc,coriolis,xloadold,&
                reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,&
                ne0,ipkon,thicke,&
                integerglob,doubleglob,tieset,istartset,&
                iendset,ialset,ntie,nasym,&
                ielprop,prop)
        endif
        !
        !         WRITE(*,*) ' mafillsm_dstil: elem',i,'ready'
        do jj=1,ndof*nope
          !
          j=(jj-1)/ndof+1
          k=jj-ndof*(j-1)
          !
          node1=kon(indexe+j)
          jdof1=nactdof(k,node1)
          !
          do ll=jj,ndof*nope
            !
            l=(ll-1)/ndof+1
            m=ll-ndof*(l-1)
            !
            node2=kon(indexe+l)
            jdof2=nactdof(m,node2)
            !
            !           check whether one of the DOF belongs to a SPC or MPC
            !
            if((jdof1.gt.0).and.(jdof2.gt.0)) then
               if(stiffonly(1).eq.1) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,&
                       s(jj,ll),jj,ll)
               else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,&
                       s(jj,ll),sm(jj,ll),jj,ll)
               endif
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
               !
               !              idof1: genuine DOF
               !              idof2: nominal DOF of the SPC/MPC
               !
               if(jdof1.le.0) then
                  idof1=jdof2
                  idof2=jdof1
               else
                  idof1=jdof1
                  idof2=jdof2
               endif
               if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
                     !
                     !                    regular DOF / MPC
                     !
                     id=(-idof2+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdof(nodempc(2,index),nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof1.eq.idof2) value=2.d0*value
                        if(idof2.gt.0) then
                           if(stiffonly(1).eq.1) then
                              call add_sm_st(au,ad,jq,irow,idof1,&
                                   idof2,value,i0,i0)
                           else
                              valu2=-coefmpc(index)*sm(jj,ll)/&
                                     coefmpc(ist)
                              !
                              if(idof1.eq.idof2) valu2=2.d0*valu2
                              !
                              call add_sm_ei(au,ad,aub,adb,jq,irow,&
                                   idof1,idof2,value,valu2,i0,i0)
                           endif
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
               !
               !              regular DOF / SPC
               !
               if(rhsi.eq.1) then
               elseif(nmethod.eq.2) then
                  value=s(jj,ll)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
               endif
            else
               idof1=jdof1
               idof2=jdof2
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
                     !
                     !                    MPC id1 / MPC id1
                     !
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),&
                                      nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdof(nodempc(2,index2),&
                                         nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*&
                                s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(1).eq.1) then
                                 call add_sm_st(au,ad,jq,irow,&
                                   idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*&
                                   sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_ei(au,ad,aub,adb,jq,&
                                   irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
                           !
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
                     !
                     !                    MPC id1 / MPC id2
                     !
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),&
                                      nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdof(nodempc(2,index2),&
                                         nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*&
                                s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                           if(idof1.eq.idof2) value=2.d0*value
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(1).eq.1) then
                                 call add_sm_st(au,ad,jq,irow,&
                                   idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*&
                                   sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                                 !
                                 if(idof1.eq.idof2) valu2=2.d0*valu2
                                 !
                                 call add_sm_ei(au,ad,aub,adb,jq,&
                                   irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
                           !
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
          !
          !           WRITE(*,*)'                  add sm ready'
          if(rhsi.eq.1) then
             !
             !            distributed forces
             !
             if(idist.ne.0) then
                !
                !               updating the external force vector for dynamic
                !               calculations
                !
                if(nmethod.eq.4) fnext(k,node1)=fnext(k,node1)+ff(jj)
                !
                if(jdof1.le.0) then
                   if(nmpc.ne.0) then
                      idof1=jdof1
                      if(idof1.ne.2*(idof1/2)) then
                         id=(-idof1+1)/2
                         ist=ipompc(id)
                         index=nodempc(3,ist)
                         if(index.eq.0) cycle
                         do
                            jdof1=nactdof(nodempc(2,index),&
                                 nodempc(1,index))
                            if(jdof1.gt.0) then
                               fext(jdof1)=fext(jdof1)&
                                    -coefmpc(index)*ff(jj)&
                                    /coefmpc(ist)
                            endif
                            index=nodempc(3,index)
                            if(index.eq.0) exit
                         enddo
                      endif
                   endif
                   cycle
                endif
                fext(jdof1)=fext(jdof1)+ff(jj)
             endif
          endif
        !
        enddo
      enddo
      !
      endif
      !       WRITE(*,*) ' mafillsm_dstil: mech ready'
      if(ithermal(1).gt.1) then
      !
      !     thermal analysis: loop over all elements
      !
      do i=nea,neb
        !
        if((ipkon(i).lt.0).or.(lakon(i)(1:1).eq.'F')) cycle
        indexe=ipkon(i)
        if(lakon(i)(4:5).eq.'20') then
           nope=20
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
         elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).ne.'A')) then
            !
            !          contact spring and advection elements
            !
            nope=ichar(lakon(i)(8:8))-47
           !
           !          local contact spring number
           !
           if(lakon(i)(7:7).eq.'C') then
              if(mortar.eq.1) nope=kon(indexe)
           endif
        elseif((lakon(i)(1:2).eq.'D ').or.&
               ((lakon(i)(1:1).eq.'D').and.(network.eq.1))) then
           !
           !          asymmetrical contribution -> mafillsmas.f
           !
           cycle
        else
           cycle
        endif
        !
        !        TODO update e_c3d_th_dstil.f
        call e_c3d_th_dstil(co,nk,kon,lakon(i),s,sm,&
        ff,i,nmethod,rhcon,nrhcon,ielmat,ielorien,norien,orab,&
        ntmat_,t0,t1,ithermal,vold,iperturb,nelemload,&
        sideload,xload,nload,idist,iexpl,dtime,&
        matname,mi(1),mass(2),stiffness,buckling,rhsi,intscheme,&
        physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,&
        xstiff,xloadold,reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,&
        ilmpc,springarea,plkcon,nplkcon,npmat_,ncmat_,elcon,nelcon,&
        lakon,pslavsurf,pmastsurf,mortar,clearini,plicon,nplicon,&
        ipkon,ielprop,prop,iponoel,inoel,sti,xstateini,xstate,&
        nstate_,network,ipobody,xbody,ibody,&
        islavelinv,islavsurf,autloc1,irowtloc1,jqtloc1)
        !
        do jj=1,nope
          !
          j=jj
          !
          node1=kon(indexe+j)
          jdof1=nactdof(0,node1)
          !
          do ll=jj,nope
            !
            l=ll
            !
            node2=kon(indexe+l)
            jdof2=nactdof(0,node2)
            !
            !           check whether one of the DOF belongs to a SPC or MPC
            !
            if((jdof1.gt.0).and.(jdof2.gt.0)) then
               if(stiffonly(2).eq.1) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,&
                       s(jj,ll),jj,ll)
               else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,&
                       s(jj,ll),sm(jj,ll),jj,ll)
               endif
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
               !
               !              idof1: genuine DOF
               !              idof2: nominal DOF of the SPC/MPC
               !
               if(jdof1.le.0) then
                  idof1=jdof2
                  idof2=jdof1
               else
                  idof1=jdof1
                  idof2=jdof2
               endif
               if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
                     !
                     !                    regular DOF / MPC
                     !
                     id=(-idof2+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdof(nodempc(2,index),nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof1.eq.idof2) value=2.d0*value
                        if(idof2.gt.0) then
                           if(stiffonly(2).eq.1) then
                              call add_sm_st(au,ad,jq,irow,idof1,&
                                   idof2,value,i0,i0)
                           else
                              valu2=-coefmpc(index)*sm(jj,ll)/&
                                     coefmpc(ist)
                              !
                              if(idof1.eq.idof2) valu2=2.d0*valu2
                              !
                              call add_sm_ei(au,ad,aub,adb,jq,irow,&
                                   idof1,idof2,value,valu2,i0,i0)
                           endif
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
               !
               !              regular DOF / SPC
               !
               if(rhsi.eq.1) then
               elseif(nmethod.eq.2) then
                  value=s(jj,ll)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
               endif
            else
               idof1=jdof1
               idof2=jdof2
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
                     !
                     !                    MPC id1 / MPC id1
                     !
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),&
                                      nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdof(nodempc(2,index2),&
                                         nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*&
                                s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(2).eq.1) then
                                 call add_sm_st(au,ad,jq,irow,&
                                   idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*&
                                   sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_ei(au,ad,aub,adb,jq,&
                                   irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
                           !
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
                     !
                     !                    MPC id1 / MPC id2
                     !
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),&
                                      nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdof(nodempc(2,index2),&
                                         nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*&
                                s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                           if(idof1.eq.idof2) value=2.d0*value
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(2).eq.1) then
                                 call add_sm_st(au,ad,jq,irow,&
                                   idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*&
                                   sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                                 !
                                 if(idof1.eq.idof2) valu2=2.d0*valu2
                                 !
                                 call add_sm_ei(au,ad,aub,adb,jq,&
                                   irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
                           !
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
          !
          !           WRITE(*,*) ' mafillsm_dstil: df'
          if(rhsi.eq.1) then
             !
             !            distributed forces
             !
             if(idist.ne.0) then
                if(jdof1.le.0) then
                   if(nmpc.ne.0) then
                      idof1=jdof1
                      if(idof1.ne.2*(idof1/2)) then
                         id=(-idof1+1)/2
                         ist=ipompc(id)
                         index=nodempc(3,ist)
                         if(index.eq.0) cycle
                         do
                            jdof1=nactdof(nodempc(2,index),&
                                 nodempc(1,index))
                            if(jdof1.gt.0) then
                               fext(jdof1)=fext(jdof1)&
                                    -coefmpc(index)*ff(jj)&
                                    /coefmpc(ist)
                            endif
                            index=nodempc(3,index)
                            if(index.eq.0) exit
                         enddo
                      endif
                   endif
                   cycle
                endif
                fext(jdof1)=fext(jdof1)+ff(jj)
             endif
          endif
        !
        enddo
      enddo
      !
      endif
      !
      !       WRITE(*,*) ' mafillsm_dstil ready: elems',nea,neb
      return
      end
