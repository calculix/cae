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
      subroutine resultsmech_u(co,kon,ipkon,lakon,ne,v,&
        stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,&
        ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,&
        iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,&
        veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,&
        xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,&
        ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,&
        reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,&
        ikin,nal,ne0,thicke,emeini,nelem,ielprop,prop)
      !
      !     calculates nal,qa,fn,xstiff,ener,eme,eei,stx for user elements
      !
      implicit none
      !
      character*8 lakon(*)
      character*80 matname(*)
      !
      integer kon(*),mi(*),&
        nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),&
        ielorien(mi(3),*),ntmat_,ipkon(*),ne0,&
        istep,iinc,ne,ithermal(2),iprestr,&
        nener,norien,iperturb(*),iout,&
        nal,icmd,nmethod,ielas,&
        ncmat_,nstate_,ikin,ielprop(*),&
        nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,&
        calcul_cauchy,calcul_qa,nelem
      !
      real*8 co(3,*),v(0:mi(2),*),stiini(6,mi(1),*),&
        stx(6,mi(1),*),prop(*),&
        elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),&
        alcon(0:6,ntmat_,*),vini(0:mi(2),*),&
        alzero(*),orab(7,*),fn(0:mi(2),*),&
        t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),&
        vold(0:mi(2),*),ener(mi(1),*),eei(6,mi(1),*),enerini(mi(1),*),&
        veold(0:mi(2),*),qa(*),dtime,time,ttime,&
        plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),&
        xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),&
        xstateini(nstate_,mi(1),*),reltime,&
        thicke(mi(3),*),emeini(6,mi(1),*)
      !
      intent(in) co,kon,ipkon,lakon,ne,v,&
        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,&
        ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,&
        iprestr,iperturb,iout,vold,nmethod,&
        veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,&
        xstateini,xstate,npmat_,matname,mi,ielas,icmd,&
        ncmat_,nstate_,stiini,vini,enerini,istep,iinc,&
        reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,&
        ikin,ne0,thicke,emeini,nelem,ielprop,prop
      !
      intent(inout) nal,qa,fn,xstiff,ener,eme,eei,stx
      !
      if(lakon(nelem)(2:2).eq.'1') then
         call resultsmech_u1(co,kon,ipkon,lakon,ne,v,&
              stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,&
              ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,&
              iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,&
              veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,&
              xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,&
              ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,&
              reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,&
              ikin,nal,ne0,thicke,emeini,nelem,ielprop,prop)
      else
         write(*,*) '*ERROR in resultsmech_u.f: user element'
         write(*,*) '       ',lakon(nelem)(1:5),' is not defined'
         call exit(201)
      endif
      !
      return
      end
