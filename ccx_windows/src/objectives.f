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
      subroutine objectives(inpc,textpart,istep,istat,n,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc,nener,nobject,objectset,objective_flag,
     &        set,nset,ntie,tieset,ier)        
!
!     reading the input deck: *OBJECTIVE
!
!     criteria: DISPLACEMENT
!               X-DISP
!               Y-DISP
!               Z-DISP
!               EIGENFREQUENCY
!               GREEN
!               MASS
!               STRAIN ENERGY
!               STRESS
!            
      implicit none
!
      logical objective_flag
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objectset(4,*),set(*),tieset(3,*)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),nset,
     &  inp(3,*),ipoinpc(0:*),nener,nobject,k,ipos,nconst,icoordinate,
     &  ntie,ier
!
      real*8 rho,stress
!
!     initialization
!
      rho=0.d0
      stress=0.d0
      icoordinate=0
!
!     check whether the design variables are the coordinates
!
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'D') then
            if(tieset(1,i)(1:10).eq.'COORDINATE') then
               icoordinate=1
               exit
            elseif(tieset(1,i)(1:11).eq.'ORIENTATION') then
               exit
            endif
         endif
      enddo
!
      if(objective_flag) then
         write(*,*) '*ERROR reading *OBJECTIVE'
         write(*,*) '       no more than one *OBJECTIVE keyword'
         write(*,*) '       is allowed per *SENSITIVITY step'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OBJECTIVE: *OBJECTIVE can
     &only be used within a SENSITIVITY STEP'     
         ier=1
         return
      endif
!
!     check if constraints are already defined
!
      nconst=0
      if(nobject.gt.0) then
         nconst=nobject
         nobject=0
      endif
!
!     check if it is a minimization or maximization problem
!
      if(textpart(2)(1:7).eq.'TARGET=') then
         if(textpart(1)(8:10).eq.'MIN') then
            objectset(2,1)(17:19)='MIN'
         elseif(textpart(2)(8:10).eq.'MAX') then 
            objectset(2,1)(17:19)='MAX'  
         else
            write(*,*) 
     &        '*WARNING optimization TARGET not known.'
            write(*,*)
     &        '         Minimization problem assumed as default.'
            objectset(2,1)(17:19)='MIN' 
         endif              
      else
         write(*,*) 
     &     '*WARNING optimization TARGET not specified.'
         write(*,*)
     &        '      Minimization problem assumed as default.'
         objectset(2,1)(17:19)='MIN' 
      endif        
!
!     reading the objectives
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      do
         if(textpart(1)(1:12).eq.'DISPLACEMENT') then
            nobject=nobject+1
            objectset(1,nobject)(1:12)='DISPLACEMENT'
            do k=13,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='N'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *OBJECTIVE: node set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
         elseif(textpart(1)(1:6).eq.'X-DISP') then
            nobject=nobject+1
            objectset(1,nobject)(1:6)='X-DISP'
            do k=7,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='N'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *OBJECTIVE: node set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
         elseif(textpart(1)(1:6).eq.'Y-DISP') then
            nobject=nobject+1
            objectset(1,nobject)(1:6)='Y-DISP'
            do k=7,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='N'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *OBJECTIVE: node set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
         elseif(textpart(1)(1:6).eq.'Z-DISP') then
            nobject=nobject+1
            objectset(1,nobject)(1:6)='Z-DISP'
            do k=7,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='N'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *OBJECTIVE: node set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
         elseif(textpart(1)(1:14).eq.'EIGENFREQUENCY') then
            nobject=nobject+1
            objectset(1,nobject)(1:14)='EIGENFREQUENCY'
            do k=15,20
               objectset(1,nobject)(k:k)=' '
            enddo
         elseif(textpart(1)(1:5).eq.'GREEN') then
            nobject=nobject+1
            objectset(1,nobject)(1:5)='GREEN'
            do k=6,20
               objectset(1,nobject)(k:k)=' '
            enddo
         elseif(textpart(1)(1:4).eq.'MASS') then
            nobject=nobject+1
            objectset(1,nobject)(1:4)='MASS'
            do k=5,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='E'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) 
     &                  '*ERROR reading *OBJECTIVE: element set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
         elseif(textpart(1)(1:12).eq.'STRAINENERGY') then
            nobject=nobject+1
            objectset(1,nobject)(1:12)='STRAINENERGY'
            do k=13,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='E'
!
!              check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) 
     &                 '*ERROR reading *OBJECTIVE: element set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
            nener=1
         elseif(textpart(1)(1:6).eq.'STRESS') then
            nobject=nobject+1
            objectset(1,nobject)(1:6)='STRESS'
            do k=7,20
               objectset(1,nobject)(k:k)=' '
            enddo
            if(n.ge.2) then
               read(textpart(2)(1:80),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               if(ipos.ne.1) then
                  objectset(3,nobject)(ipos:ipos)='N'
!
!                 check the existence of the set
!
                  do i=1,nset
                     if(set(i).eq.objectset(3,nobject)) exit
                  enddo
                  if(i.gt.nset) then
                     objectset(3,nobject)(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *OBJECTIVE: node set ',
     &                    objectset(3,nobject)
                     write(*,*) '       has not yet been defined. '
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*OBJECTIVE%",ier)
                     return
                  endif
               endif
            endif
!
!           rho for the Kreisselmeier-Steinhauser function
!
            if(n.ge.3) then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat) rho
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*OBJECTIVE%",ier)
                  return
               endif
               objectset(2,nobject)(41:60)=textpart(3)(1:20)
            endif
!
            if(icoordinate.eq.1) then
               if(rho.lt.1.d0) then
                  write(*,*) '*ERROR reading *OBJECTIVE'
                  write(*,*) '       first Kreisselmeier-Steinhauser'
                  write(*,*) '       parameter rho cannot be less'
                  write(*,*) '       than 1'
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*OBJECTIVE%",ier)
                  return
               endif
            endif
!
!           the target stress for the Kreisselmeier-Steinhauser function
!
            if(n.ge.4) then
               read(textpart(4)(1:20),'(f20.0)',iostat=istat) stress
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*OBJECTIVE%",ier)
                  return
               endif
               objectset(2,nobject)(61:80)=textpart(4)(1:20)
            endif
!
            if(stress.le.0.d0) then
               if(icoordinate.eq.1) then
                  write(*,*) '*ERROR reading *OBJECTIVE'
                  write(*,*) '       the target stress in the'
                  write(*,*) '       Kreisselmeier-Steinhauser function'
                  write(*,*) '       must be strictly positive'
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*OBJECTIVE%",ier)
                  return
               endif
            endif
         else
            write(*,*) '*ERROR reading *OBJECTIVE'
            write(*,*) '       objective function not known'
            call inputerror(inpc,ipoinpc,iline,
     &           "*OBJECTIVE%",ier)
            return
         endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!
!        if constraints are defined only one objective function is allowed
!
         if((nconst.gt.0).and.(nobject.gt.1)) then
            write(*,*) '*ERROR reading *OBJECTIVE'
            write(*,*) '       in the case constraints are defined,'
            write(*,*) '       the definition of only 1 objective' 
            write(*,*) '       is allowed'
            call inputerror(inpc,ipoinpc,iline,
     &           "*OBJECTIVE%",ier)
            return
         endif
             
      enddo
!  
!     In the case constraints are already defined, the actual number of
!     nobjects is restored
!   
      if(nconst.gt.0) then
         nobject=nconst
      endif
!
      return
      end
      
      
