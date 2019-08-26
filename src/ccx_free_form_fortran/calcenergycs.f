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
      subroutine calcenergycs(co,kon,ipkon,lakon,ne,v,stx,ielmat,&
        ielorien,norien,orab,t1,ithermal,eme,fn,vold,veold,time,ttime,&
        xstate,mi,nstate_,ener,eei,ikin,ne0,thicke,mortar,ielprop,prop,&
        prlab,nprint,nload,nelemload,xload,sideload,trab,inotr,ntrans,&
        set,nset,istartset,iendset,ialset,prset,qfx,orname,nk,islavsurf,&
        filab,inum,enern)
      !
      !     calculates the real and imaginary part of the strain energy
      !     for cyclic symmetric frequency calculations and stores the
      !     result (if requested) in the .dat file
      !     (relevant requests: ENER or ELSE on the *EL PRINT card,
      !      CELS on the *CONTACT PRINT card)
      !
      !     extrapolates the values to the nodes for frd-output
      !     (relevant requests: ENER on the *EL FILE/ELEMENT OUTPUT card,
      !      CELS on the *CONTACT FILE/CONTACT OUTPUT card)
      !
      implicit none
      !
      logical force
      !
      character*1 cflag
      character*6 prlab(*),prlabcp(nprint)
      character*8 lakon(*),lakonl
      character*20 sideload(*)
      character*80 orname(*)
      character*81 set(*),prset(*)
      character*87 filab(*)
      !
      integer kon(*),mi(*),ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),&
        ne0,null,ne,ithermal(2),i,k,jj,norien,mint3d,nstate_,ikin,&
        nlayer,ielprop(*),mortar,ntrans,nset,nload,nelemload(2,*),nk,&
        istartset(*),inotr(2,*),islavsurf(2,*),iendset(*),ialset(*),&
        inum(*),nprint,ndim,nfield,iorienloc
      !
      real*8 co(3,*),v(0:mi(2),*),stx(6,mi(1),*),prop(*),orab(7,*),&
        fn(0:mi(2),*),t1(*),eme(6,mi(1),*),vold(0:mi(2),*),&
        veold(0:mi(2),*),ener(mi(1),*),eei(6,mi(1),*),xi,et,ze,&
        weight,time,ttime,xstate(nstate_,mi(1),*),thicke(mi(3),*),&
        xload(2,*),trab(7,*),qfx(3,mi(1),*),enern(*)
      !
      intent(in) co,kon,ipkon,lakon,ne,v,ielmat,ielorien,norien,orab,t1,&
        ithermal,vold,veold,time,ttime,xstate,mi,nstate_,ikin,ne0,&
        thicke,mortar,ielprop,prop,xload,trab,set,sideload,qfx,prset,&
        orname,ntrans,nset,nload,nelemload,nk,istartset,inotr,islavsurf,&
        iendset,ialset,filab,inum,prlab,nprint,enern
      !
      intent(inout) fn,ener,eme,eei,stx
      !
      null=0
      !
      !     determining the real and imaginary part of the strain energy
      !
      do i=1,ne
         !
         lakonl=lakon(i)
         !
         !        determining the number of integration points
         !
         if(lakonl(1:1).eq.'U') then
            mint3d=ichar(lakonl(6:6))
         elseif(lakonl(4:5).eq.'8R') then
            mint3d=1
         elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
               mint3d=50
            else
               call beamintscheme(lakonl,mint3d,ielprop(i),prop,&
                    null,xi,et,ze,weight)
            endif
         elseif((lakonl(4:4).eq.'8').or.&
                (lakonl(4:6).eq.'20R')) then
            if(lakonl(7:8).eq.'LC') then
               mint3d=8*nlayer
            else
               mint3d=8
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).eq.'LC') then
               mint3d=6*nlayer
            else
               mint3d=9
            endif
         elseif(lakonl(4:4).eq.'6') then
            mint3d=2
         !
         !           contact elements
         !
         elseif(i.gt.ne0) then
            mint3d=1
         elseif(lakonl(1:1).eq.'E') then
            mint3d=0
         endif
         !
         k=ne+i
         !
         !        energy_r=eme_r*stx_r-eme_i*stx_i
         !        energy_i=eme_r*stx_i+eme_i*stx_r
         !
         do jj=1,mint3d
            ener(jj,i)=((eme(1,jj,i)*stx(1,jj,i)+&
                         eme(2,jj,i)*stx(2,jj,i)+&
                         eme(3,jj,i)*stx(3,jj,i))/2.d0+&
                         eme(4,jj,i)*stx(4,jj,i)+&
                         eme(5,jj,i)*stx(5,jj,i)+&
                         eme(6,jj,i)*stx(6,jj,i))-&
                       ((eme(1,jj,k)*stx(1,jj,k)+&
                         eme(2,jj,k)*stx(2,jj,k)+&
                         eme(3,jj,k)*stx(3,jj,k))/2.d0+&
                         eme(4,jj,k)*stx(4,jj,k)+&
                         eme(5,jj,k)*stx(5,jj,k)+&
                         eme(6,jj,k)*stx(6,jj,k))
            ener(jj,k)=((eme(1,jj,i)*stx(1,jj,k)+&
                         eme(2,jj,i)*stx(2,jj,k)+&
                         eme(3,jj,i)*stx(3,jj,k))/2.d0+&
                         eme(4,jj,i)*stx(4,jj,k)+&
                         eme(5,jj,i)*stx(5,jj,k)+&
                         eme(6,jj,i)*stx(6,jj,k))+&
                       ((eme(1,jj,k)*stx(1,jj,i)+&
                         eme(2,jj,k)*stx(2,jj,i)+&
                         eme(3,jj,k)*stx(3,jj,i))/2.d0+&
                         eme(4,jj,k)*stx(4,jj,i)+&
                         eme(5,jj,k)*stx(5,jj,i)+&
                         eme(6,jj,k)*stx(6,jj,i))
         enddo
      enddo
      !
      !     storage in the .dat file
      !
      !     real values
      !
      call writere()
      !       write(5,*)
      !       write(5,*) ' REAL VALUES '
      !
      do i=1,nprint
         if((prlab(i)(1:4).eq.'ENER').or.&
            (prlab(i)(1:4).eq.'ELSE').or.&
            (prlab(i)(1:4).eq.'CELS')) then
            prlabcp(i)=prlab(i)
         else
            prlabcp(i)='      '
         endif
      enddo
      !
      call printout(set,nset,istartset,iendset,ialset,nprint,&
        prlabcp,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,&
        mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,&
        orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,thicke,&
        eme,islavsurf,mortar,time,ielprop,prop,veold,orname,&
        nelemload,nload,sideload,xload)
      !
      !     imaginary values
      !
      call writeim()
      !       write(5,*)
      !       write(5,*) ' IMAGINARY VALUES '
      !
      call printout(set,nset,istartset,iendset,ialset,nprint,&
        prlabcp,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener(1,ne+1),&
        mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,&
        orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,thicke,&
        eme,islavsurf,mortar,time,ielprop,prop,veold,orname,&
        nelemload,nload,sideload,xload)
      !
      !     determining the internal energy in the nodes
      !     for output in frd format
      !
      if(filab(7)(1:4).eq.'ENER') then
         nfield=1
         ndim=1
         iorienloc=0
         cflag=filab(7)(5:5)
         force=.false.
         call extrapolate(ener,enern,ipkon,inum,kon,lakon,nfield,nk,&
              ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,&
              vold,force,ielmat,thicke,ielprop,prop)
         call extrapolate(ener(1,ne+1),enern(nk+1),ipkon,inum,kon,&
              lakon,nfield,nk,&
              ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,&
              vold,force,ielmat,thicke,ielprop,prop)
      endif
      !
      return
      end
