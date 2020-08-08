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
      subroutine submodels(inpc,textpart,set,istartset,iendset,ialset,
     &  nset,nset_,nalset,nalset_,nk,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,nsubmodel,tieset,tietol,ntie,
     &  ntie_,jobnamec,amta,namta,nam,nam_,namtot_,ier)
!
!     reading the input deck: *SUBMODEL
!
      implicit none
!
      character*1 type,inpc(*)
      character*81 set(*),noset,submset,globset,facialset,
     &  tieset(3,*)
      character*132 textpart(16),jobnamec(*)
!
      integer nset,nset_,nalset,nalset_,istep,istat,n,key,i,nk,
     &  j,istartset(*),iendset(*),ialset(*),ipos,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),l,k,kstart,kend,ipoinpc(0:*),
     &  iglobset,kincrement,nsubmodel,kflag,idummy,nalsetold,ntie,ntie_,
     &  nlength,namta(3,*),nam,nam_,namtot_,ier
!
      real*8 tietol(3,*),amta(2,*)
!
      data kflag /1/
!
      if(istep.gt.0) then
         write(*,*) 
     &          '*ERROR reading *SUBMODEL: *SURFACE should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR reading *SUBMODEL: increase ntie_'
         ier=1
         return
      endif
!
      tietol(1,ntie)=-1.d0
      iglobset=0
!
!     if no *AMPLITUDE in the input deck, add a fake amplitude
!     iamboun and iamload are used to store the global step used
!     for interpolation; these fields are not allocated if the
!     number of amplitudes is zero
!
      if(nam.eq.0) then
         nam=1
         if(nam_.lt.1) then
            write(*,*) '*ERROR reading *SUBMODEL: increase nam_'
            ier=1
            return
         endif
         if(namtot_.lt.2) then
            write(*,*) 
     &           '*ERROR reading *SUBMODEL: increase namtot_'
            ier=1
            return
         endif
         amta(1,1)=0.d0
         amta(2,1)=1.d0
         amta(1,2)=1.d0
         amta(2,2)=1.d0
         namta(1,1)=1
         namta(2,1)=2
         namta(3,1)=1
      endif
!
      type='N'
c      globset(1:1)=' '
      nsubmodel=nsubmodel+1
      submset(1:8)='SUBMODEL'
      if(nsubmodel.lt.10) then
         submset(9:10)='00'
         write(submset(11:11),'(i1)') nsubmodel
      elseif(nsubmodel.lt.100) then
         submset(9:9)='0'
         write(submset(10:11),'(i2)') nsubmodel
      elseif(nsubmodel.lt.1000) then
         write(submset(9:11),'(i3)') nsubmodel
      else
         write(*,*) '*ERROR reading *SUBMODEL: no more than 999'
         write(*,*) '       submodels allowed'
         ier=1
         return
      endif
      do i=12,81
         submset(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:12).eq.'SURFACE') then
               type='T'
            elseif(textpart(i)(6:9).eq.'NODE') then
               type='N'
            else
               write(*,*) 
     &          '*ERROR reading *SUBMODEL: unknown type'
               ier=1
               return
            endif
         elseif(textpart(i)(1:12).eq.'GLOBALELSET=') then
            globset(1:80)=textpart(i)(13:92)
            globset(81:81)=' '
            ipos=index(globset,' ')
            globset(ipos:ipos)='E'
            do iglobset=1,nset
               if(set(iglobset).eq.globset) exit
            enddo
            if(iglobset.gt.nset) then
               write(*,*) 
     &            '*ERROR reading *SUBMODEL: global element set ',
     &            globset(1:ipos-1)
               write(*,*) '       does not exist'
               ier=1
               return
            endif
            do j=istartset(iglobset),iendset(iglobset)
               if(ialset(j).lt.0) then
                  write(*,*) 
     &               '*ERROR reading *SUBMODEL: global element set ',
     &               globset(1:ipos-1)
                  write(*,*) '       was defined using GENERATE;'
                  write(*,*) '       this is not allowed'
                  ier=1
                  return
               endif
            enddo
         elseif(textpart(i)(1:6).eq.'INPUT=') then
            jobnamec(4)(1:126)=textpart(i)(7:132)
            jobnamec(4)(127:132)='      '
            loop1: do j=1,126
               if(jobnamec(4)(j:j).eq.'"') then
                  do k=j+1,126
                     if(jobnamec(4)(k:k).eq.'"') then
                        do l=k-1,126
                           jobnamec(4)(l:l)=' '
                           exit loop1
                        enddo
                     endif
                     jobnamec(4)(k-1:k-1)=jobnamec(4)(k:k)
                  enddo
                  jobnamec(4)(126:126)=' '
               endif
            enddo loop1
         else
            write(*,*) 
     &        '*WARNING reading *SUBMODEL: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SUBMODEL%")
         endif
      enddo
!
!     check whether a global elset was defined
!
c      if(globset(1:1).eq.' ') then
c         write(*,*) '*ERROR reading *SUBMODEL: no global elset'
c         write(*,*) '       defined'
c         ier=1
c         return
c      endif
!
c      ipos=index(submset,' ')
c      if(ipos.eq.1) then
c         write(*,*) '*ERROR reading *SUBMODEL: no name specified'
c         ier=1
c         return
c      endif
      submset(12:12)=type
!
!     submodel set (nodal set or element face set)
!
      nset=nset+1
      if(nset.gt.nset_) then
         write(*,*) '*ERROR reading *SUBMODEL: increase nset_'
         ier=1
         return
      endif
      set(nset)=submset
      istartset(nset)=nalset+1
      iendset(nset)=0
!
      if(type.eq.'N') then
!
!        node surface
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               if(iendset(nset).eq.0) then
                  write(*,*) '*ERROR reading *SUBMODEL: nodal'
                  write(*,*) '       submodel contains no nodes'
                  ier=1
                  return
               endif
               exit
            endif
!
            do l=1,n
               if(nalset+1.gt.nalset_) then
                  write(*,*) 
     &               '*ERROR reading *SUBMODEL: increase nalset_'
                  ier=1
                  return
               endif
!     
               read(textpart(l)(1:10),'(i10)',iostat=istat)
     &             ialset(nalset+1)
               if(istat.gt.0) then
                  noset=textpart(l)(1:80)
                  noset(81:81)=' '
                  ipos=index(noset,' ')
                  noset(ipos:ipos)='N'
                  do i=1,nset
                     if(set(i).eq.noset) then
                        do j=istartset(i),iendset(i)
                           if(ialset(j).gt.0) then
                              nalset=nalset+1
                              if(nalset.gt.nalset_) then
                                 write(*,*) 
     &                   '*ERROR reading *SUBMODEL: increase nalset_'
                                 ier=1
                                 return
                              endif
                              ialset(nalset)=ialset(j)
                           else
                              kstart=ialset(nalset-1)
                              kend=ialset(nalset)
                              nalset=nalset-1
                              kincrement=-ialset(j)
                              do k=kstart+kincrement,kend,kincrement
                                 nalset=nalset+1
                                 if(nalset.gt.nalset_) then
                                    write(*,*) 
     &                     '*ERROR reading *SUBMODEL: increase nalset_'
                                    ier=1
                                    return
                                 endif
                                 ialset(nalset)=k
                              enddo
                           endif
                        enddo
                        iendset(nset)=nalset
                        exit
                     endif
                  enddo
                  if(i.gt.nset-1) then
                     noset(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *SUBMODEL: node set ',
     &                     noset
                     write(*,*) '       does not exist'
                     ier=1
                     return
                  endif
               else
                  if(ialset(nalset+1).gt.nk) then
                     write(*,*) '*WARNING reading *SUBMODEL: value ',
     &                    ialset(nalset+1)
                     write(*,*) '         in set ',set(nset),' > nk'
                  else
                     nalset=nalset+1
                     iendset(nset)=nalset
                  endif
               endif
            enddo
         enddo
!     
      else
!
!        facial submodel interface surface
!
         loop:do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
!
!           check whether any facial surface was defined
!
            if((istat.lt.0).or.(key.eq.1)) then
               if(iendset(nset).eq.0) then
                  write(*,*) '*ERROR reading *SUBMODEL: facial'
                  write(*,*) '       submodel contains no faces'
                  ier=1
                  return
               endif
               exit
            endif
!
            if(n.gt.1) then
               write(*,*) '*ERROR reading *SUBMODEL: no more than'
               write(*,*) '       one entry per line allowed for'
               write(*,*) '       a facial submodel interface'
               ier=1
               return
            endif
!
            facialset(1:80)=textpart(1)(1:80)
            facialset(81:81)=' '
            ipos=index(facialset,' ')
            facialset(ipos:ipos)='T'
            do i=1,nset
               if(set(i).eq.facialset) then
                  do j=istartset(i),iendset(i)
                     nalset=nalset+1
                     if(nalset.gt.nalset_) then
                        write(*,*)
     &                   '*ERROR in reading *SUBMODEL: increase nalset_'
                        ier=1
                        return
                     endif
                     ialset(nalset)=ialset(j)
                  enddo
                  iendset(nset)=nalset
                  cycle loop
               endif
            enddo
         enddo loop
      endif
!
      tieset(1,ntie)(81:81)='S'
!
!     sorting the local (nodal or element face) set
!
      write(tieset(2,ntie)(1:10),'(i10)') nset
      tieset(2,ntie)(11:11)=type
      nlength=iendset(nset)-istartset(nset)+1
      call isortii(ialset(istartset(nset)),idummy,nlength,kflag)
!
!     expanding and sorting the global element set
!
      write(tieset(3,ntie)(1:10),'(i10)') iglobset
!
!     if iglobset equals zero (i.e. no global element set was
!     defined), the complete global model is taken as the master model
!
      if(iglobset.eq.0) return
!
!     expanding the global element set (no "generate" allowed)
!
      nalsetold=nalset+1
      do j=istartset(iglobset),iendset(iglobset)
         if(ialset(j).gt.0) then
            nalset=nalset+1
            if(nalset.gt.nalset_) then
               write(*,*) 
     &              '*ERROR reading *SUBMODEL: increase nalset_'
               ier=1
               return
            endif
            ialset(nalset)=ialset(j)
         else
            kstart=ialset(nalset-1)
            kend=ialset(nalset)
            nalset=nalset-1
            kincrement=-ialset(j)
            do k=kstart+kincrement,kend,kincrement
               nalset=nalset+1
               if(nalset.gt.nalset_) then
                  write(*,*) 
     &                 '*ERROR reading *SUBMODEL: increase nalset_'
                  ier=1
                  return
               endif
               ialset(nalset)=k
            enddo
         endif
      enddo
      istartset(iglobset)=nalsetold
      iendset(iglobset)=nalset
!
!     sorting the global element set
!  
      nlength=iendset(iglobset)-istartset(iglobset)+1
      call isortii(ialset(istartset(iglobset)),idummy,nlength,kflag)
!     
      return
      end

