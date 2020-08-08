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
!     check wether the property array contains values or IDs
!     if first gas iteration (iin =0) then
!     for each property of each fluid element
!     if IDs: the IDs contained in prop array are stored in prop_store array;
!     the IDs are interpreted and values are set/stored in prop array
!     for each property of each fluid element
!     else (after convergence)
!     if IDs: the IDs stored in prop_store array are copied back in
!     prop array
!     
!     2016.05.19 Y. Muller
!     
      subroutine checkinputvaluesnet(ieg,nflow,prop,ielprop,lakon)
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer ieg(*),nflow,i,ielprop(*),index,nelem
!     
      real*8 prop(*)
!     
      do i=1,nflow
        nelem=ieg(i)
        index=ielprop(nelem)
        if(index.lt.0) cycle
!     
!     modifying the prop array (formerly done in fluidsections.f)
!     
        if(lakon(nelem)(2:8).eq.'REBRJI2') then
          if(1.d0-(prop(index+5)+prop(index+6))/
     &         prop(index+4).gt.0.01d0)then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       in element type RESTRICTOR
     &BRANCH JOINT IDELCHIK2'
            write(*,*) '       A0 ist not equal to A1+A2'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
!     
        elseif((lakon(nelem)(2:3).eq.'OR').and.
     &         (lakon(nelem)(2:5).ne.'ORC1').and.
     &         (lakon(nelem)(2:5).ne.'ORB1').and.
     &         (lakon(nelem)(2:5).ne.'ORB2').and.
     &         (lakon(nelem)(2:5).ne.'ORBT').and.
     &         (lakon(nelem)(2:5).ne.'ORPN').and.
     &         (lakon(nelem)(2:5).ne.'ORFL'))then
          if(prop(index+2).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: diameter '
            write(*,*) '       of the orifice is not positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+3).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '        length of the orifice is ',
     &           'not positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((lakon(nelem)(2:5).ne.'ORC1').and.
     &         (lakon(nelem)(2:5).ne.'ORBF').and.
     &         (lakon(nelem)(2:5).ne.'ORRG').and.
     &         (lakon(nelem)(2:5).ne.'ORSG').and.
     &         (lakon(nelem)(2:5).ne.'ORGA').and.
     &         (lakon(nelem)(2:5).ne.'ORBO'))then
            if((prop(index+4).gt.1.d-20).and.
     &           (prop(index+5).gt.1.d-20)) then
              write(*,*)
     &             '*ERROR in checkinputvaluesnet:'
              write(*,*) 'either the radius of ',
     &             'the orifice must be zero or the chamfer angle'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
          endif
          if(prop(index+4)/prop(index+2).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) 
     &           '        r/d of an orifice must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+4)/prop(index+2).gt.0.82d0) then
            write(*,*)
     &           '*WARNING in checkinputvaluesnet: '
            write(*,*) '        r/d of an orifice ',
     &           'must not exceed 0.82'
            write(*,*) '         element number: ',nelem
          endif
          if(prop(index+5).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the chamfer angle of an ',
     &           'orifice must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+5).gt.90.d0) then
            write(*,*)
     &           '*ERROR in checkinputvaluesnet: '
            write(*,*) '       the chamfer angle of an ',
     &           'orifice must not exceed 90°'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+6).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '      d/D (orifice)must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+6).gt.0.5d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       d/D (orifice)must not exceed 0.5'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+3)/prop(index+2).lt.0.d0) then
            write(*,*)
     &           '*ERROR in checkinputvaluesnet:'
            write(*,*) '        L/d of an orifice ',
     &           'must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(lakon(nelem)(4:4).eq.'P') then
            if(prop(index+3)/prop(index+2).gt.2.d0) then
              write(*,*)
     &             '*WARNING in checkinputvaluesnet: '
              write(*,*) '       L/d of an orifice ',
     &             'with Parker must not exceed 2'
              write(*,*) '       element number: ',nelem
!     call exit(201)
            endif
          endif
        elseif((lakon(nelem)(2:5).eq.'ORBT'))then
          if(prop(index+2).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '        ps1/pt1 (bleedtapping) ',
     &           'must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+2).gt.1.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       ps1/pt1 (bleed tapping) ',
     &           'must not exceed 1'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
        elseif((lakon(nelem)(2:5).eq.'ORPN'))then
          if(prop(index+2).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       theta (preswirlnozzle) ',
     &           'must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+2).gt.90.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       theta (preswirl nozzle) ',
     &           'must not exceed 90°'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+3).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       k_phi (preswirl nozzle) ',
     &           'must not be negative'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+3).gt.1.05d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       k_phi (preswirlnozzle) ',
     &           'must not exceed 1.05'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
        elseif((lakon(nelem)(2:5).eq.'ORBG'))then
          if(prop(index+1).lt.0.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       section area is not positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if(prop(index+2).lt.0.d0 .or.
     &         prop(index+2).ge.1.d0) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       using Bragg Method'
            write(*,*) '       Cd by crtitical pressure ratio '
            write(*,*) '       *FLUID SECTIONS position 2'
            write(*,*) '       0 < Cd _crit < 1'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
        elseif(lakon(nelem)(2:5).eq.'ORB1')then
          if(dabs(1.d0-prop(index+4)/50.d0).lt.0.00001d0.and.
     &         dabs(1.d0-prop(index+5)/0.00015d0).gt.0.00001d0)then    
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 50 the bristle diameter'
            write(*,*) '       has to be 0.00015m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+4)/100.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+5)/0.00007d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 100 the bristle diameter'
            write(*,*) '       has to be 0.00007m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+4)/140.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+5)/0.0001d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 140 the bristle diameter'
            write(*,*) '       has to be 0.00010m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+4)/200.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+5)/0.00007d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 200 the bristle diameter'
            write(*,*) '       has to be 0.00007m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)       
          endif
          if(dabs(1.d0-prop(index+4)/50.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+4)/100.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+4)/140.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+4)/200.d0).gt.0.00001d0)then    
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       Only a brsitle density of'
            write(*,*) '       50, 100, 140 and 200 is supported'
            write(*,*) '       with corresponding bristle diameter'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+4)
            call exit(201)
          endif
!     
        elseif(lakon(nelem)(2:5).eq.'ORB2')then
          if(dabs(1.d0-prop(index+3)/50.d0).lt.0.00001d0.and.
     &         dabs(1.d0-prop(index+4)/0.00015d0).gt.0.00001d0)then    
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 50 the bristle diameter'
            write(*,*) '       has to be 0.00015m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+3)/100.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+4)/0.00007d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 100 the bristle diameter'
            write(*,*) '       has to be 0.00007m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+3)/140.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+4)/0.0001d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 140 the bristle diameter'
            write(*,*) '       has to be 0.00010m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)
          elseif(dabs(1.d0-prop(index+3)/200.d0).lt.0.00001d0.and.
     &           dabs(1.d0-prop(index+4)/0.00007d0).gt.0.00001d0)then
            write(*,*) '*ERROR in checkinputvaluesnet: '                
            write(*,*) '       For a brush seal with a bristle'
            write(*,*) '       density of 200 the bristle diameter'
            write(*,*) '       has to be 0.00007m!'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+5)
            call exit(201)       
          endif
          if(dabs(1.d0-prop(index+3)/50.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+3)/100.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+3)/140.d0).gt.0.00001d0.and.
     &         dabs(1.d0-prop(index+3)/200.d0).gt.0.00001d0)then    
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       Only a brsitle density of'
            write(*,*) '       50, 100, 140 and 200 is supported'
            write(*,*) '       with corresponding bristle diameter'
            write(*,*) '       For element number:', nelem
            write(*,*) '       the value is:',prop(index+4)
            call exit(201)
          endif
!     
        elseif((lakon(nelem)(2:4).eq.'LAB').and.
     &         (lakon(nelem)(2:5).ne.'LABF').and.
     &         (lakon(nelem)(2:5).ne.'LABD')) then
          if((prop(index+1).gt.1000.d0)
     &         .or.(prop(index+1).lt.0.d0)) then
            write(*,*)
     &           '*ERROR in checkinputvaluesnet: ',
     &           'the selected pitch t'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=t<=1000 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+2).gt.100.d0)
     &         .or.(prop(index+2).lt.0.d0)) then
            write(*,*)
     &           '*ERROR in checkinputvaluesnet: ',
     &           'the selected gap s'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=s<=100 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+3).gt.5000.d0)
     &         .or.(prop(index+3).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected diameter d'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=d<=5000'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+5).gt.9.d0)
     &         .or.(prop(index+5).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*)'       the selected spike number n'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=n<=9'
            write(*,*) prop(index+4)
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+5).gt.100.d0)
     &         .or.(prop(index+5).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected spike breadth'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=b<=100 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+6).gt.9.d0)
     &         .or.(prop(index+6).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected spike height'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=b<=20 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+7).gt.4.d0)
     &         .or.(prop(index+7).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected Honeycomb cell width'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=L<=4 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+8).gt.5.d0)
     &         .or.(prop(index+8).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected edge radius'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=r<=5 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+9).gt.100.d0)
     &         .or.(prop(index+9).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected position of the spike'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=X<=0.1 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
          if((prop(index+10).gt.100.d0)
     &         .or.(prop(index+10).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       the selected height of the step'
            write(*,*) '       for the labyrinth is not correct'
            write(*,*) '       0<=Hst<=0.1 mm'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
!     
        elseif((lakon(nelem)(2:6).eq.'CARBS'))then
          if(lakon(nelem)(2:8).eq.'CARBSGE') then
            if(prop(index+1).le.0.d0) then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       the selected diameter of ',
     &             'the carbon seal'
              write(*,*) 
     &             '       has been defined as less or equal to 0'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
          else
            if((prop(index+1).le.0.d0)
     &           .or.(prop(index+2).le.0.d0)
     &           .or.(prop(index+3).le.0.d0)) then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       the selected diameter'
              write(*,*) '       or the selected length'
              write(*,*) 
     &             '       or the selected gap of the carbon seal'
              write(*,*) 
     &             '       has been defined as less or equal to 0'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
          endif
!     
        elseif((lakon(nelem)(2:4).eq.'RCVL'))then
          if(prop(index+2).lt.(prop(index+3))) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       element TYPE=ROTATING CAVITY ',
     &           '(Radial inflow)'
            write(*,*) '       the specified upstream radius is ',
     &           'smaller than'
            write(*,*) '       the specified downstream radius!'
            write(*,*) '       Please check the element definition.'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
        elseif((lakon(nelem)(2:4).eq.'RCVN'))then
          if(prop(index+1).lt.(prop(index+2))) then
            write(*,*) '*ERROR in checkinputvaluesnet: '
            write(*,*) '       element TYPE=ROTATING CAVITY ',
     &           '(Radial inflow)'
            write(*,*) '       the specified upstream radius is ',
     &           'smaller than'
            write(*,*) '       the specified downstream radius!'
            write(*,*) '       Please check the element definition.'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
!     
        elseif(lakon(nelem)(2:3).eq.'RE') then
          if(((prop(index+1).le.0.d0)
     &         .or.(prop(index+2).le.0.d0)
     &         .or.(prop(index+3).le.0.d0))
     &         .and.(lakon(nelem)(2:5).ne.'REBR')
     &         .and.(lakon(nelem)(2:5).ne.'REEX')
     &         .and.(lakon(nelem)(2:7).ne.'REWAOR')
     &         .and.(lakon(nelem)(2:5).ne.'REEN'))then
            write(*,*) '*ERROR in checkinputvaluesnet:'
            write(*,*) '       A1,A2 or Dh less or equal 0'
            write(*,*) '       element number: ',nelem
            call exit(201)
!     
          elseif((lakon(nelem)(2:5).eq.'REEL'))then
            if(prop(index+1).ge.(prop(index+2))) then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) 
     &             '       Section A1 is greater than section A2'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
!     
          elseif((lakon(nelem)(2:5).eq.'RECO'))then
            if(prop(index+1).lt.(prop(index+2))) then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) 
     &             '       Section A2 is greater than section A1'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
          endif
!     
          if((lakon(nelem)(2:5).eq.'REBR'))then
            if((prop(index+1).le.0.d0)
     &           .or.(prop(index+2).le.0.d0)
     &           .or.(prop(index+3).le.0.d0))then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       trying to define a branch '
              write(*,*) '       all three elements must be ',
     &             'different from 0'
              write(*,*) '       element number: ',nelem
              call exit(201)
!     
            elseif((prop(index+4).le.0.d0)
     &             .or.(prop(index+5).le.0.d0)
     &             .or.(prop(index+6).le.0.d0))then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       trying to define a branch '
              write(*,*) '       all sections must be positive'
              write(*,*) '       element number: ',nelem
              call exit(201)
!     
            elseif((prop(index+7).lt.0)
     &             .or.(prop(index+8).lt.0))then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       trying to define a branch '
              write(*,*) '       alpha1 & 2 cannot be negative'
              write(*,*) '       element number: ',nelem
              call exit(201)
!     
            elseif((prop(index+7).gt.90)
     &             .or.(prop(index+8).gt.90)) then
              write(*,*) '*ERROR in checkinputvaluesnet:'
              write(*,*) '       trying to define a branch '
              write(*,*) '       alpha1 & 2 cannot greater than ',
     &             '90 gegrees'
              write(*,*) '       element number: ',nelem
              call exit(201)
!     
            elseif((lakon(nelem)(6:8).eq.'SI1')
     &             .or.(lakon(nelem)(6:8).eq.'JI1')
     &             .or.(lakon(nelem)(6:8).eq.'JI2')) then
              if(prop(index+7).gt.0) then
                write(*,*) '*ERROR in checkinputvaluesnet:'
                write(*,*) '       trying to define a branch '
                write(*,*) 
     &               '       Type IDELCHIK JOINT1 or SPLIT1&2'
                write(*,*) '       alpha1 must be 0 degrees'
                write(*,*) '       element number: ',nelem
                call exit(201)
              endif
!     
            elseif((lakon(nelem)(6:8).eq.'SI1').or.
     &             (lakon(nelem)(6:8).eq.'JI1'))then
              if(prop(index+4).ne.(prop(index+5)))then
                write(*,*) '*ERROR in checkinputvaluesnet:'
                write(*,*) '       trying to define a branch '
                write(*,*) '       Type IDELCHIK SPLIT1 or JOINT1 '
                write(*,*) '       A1=A0'
                write(*,*) '       element number: ',nelem
                call exit(201)
              endif
!     
            elseif(lakon(nelem)(6:8).eq.'JI2') then
              if((prop(index+5)+(prop(index+6)))
     &             .ne.prop(index+4))then
                write(*,*) '*ERROR in checkinputvaluesnet:'
                write(*,*) '       trying to define a branch '
                write(*,*) '       Type IDELCHIK JOINT2 '
                write(*,*) '       A1+A2 must be equal to A0'
                write(*,*) '       element number: ',nelem
                call exit(201)
              endif
            endif
          endif
!     
!     General Vortex
!     
        elseif((lakon(nelem)(2:3).eq.'VO'))then
!     
!     inner and outer radius less or equal to 0
!     
          if((prop(index+1).le.0.d0) .or.
     &         (prop(index+2).le.0.d0)) then
            write(*,*)'*ERROR in checkinputvaluesnet:'
            write(*,*)'       trying to define a VORTEX'
            write(*,*)'       R1 and R2 must be positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
!     
          elseif(prop(index+3).le.0.d0) then
            write(*,*)'*ERROR in checkinputvaluesnet:'
            write(*,*)'       trying to define a VORTEX'
            write(*,*)'       eta must be different positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif
!     
!     FREE VORTEX
!     
          if((lakon(nelem)(2:5).eq.'VOFR'))then
!     
!     the swirl comes from another upstream  element
!     
            if(prop(index+6).ne.0d0) then
!     
!     the rotation speed must be 0
!     
              if(prop(index+7).ne.0.d0) then
                write(*,*)'*ERROR in checkinputvaluesnet:'
                write(*,*)'       trying to define a FREE VORTEX'
                write(*,*)
     &               '       rotation speed and upstream element'
                write(*,*)'       cannot be simultaneously used '
                write(*,*) '       element number: ',nelem
                call exit(201)
              endif
            endif
!     
!     FORCED VORTEX
!     
          elseif((lakon(nelem)(2:5).eq.'VOFO'))then
!     
!     Core swirl ratio must be defined and positive
!     
            if(prop(index+4).le.0.d0) then
              write(*,*)'*ERROR in checkinputvaluesnet:'
              write(*,*)'       trying to define a FORCED VORTEX'
              write(*,*)'       Core swirl ratio Kr is strictly ',
     &             'positive'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
            if(prop(index+4).gt.1.d0) then
              write(*,*)'*WARNING in checkinputvaluesnet:'
              write(*,*)'       trying to define a FORCED VORTEX'
              write(*,*)'       Core swirl ratio Kr is ',
     &             'greater than 1'
              write(*,*) '       element number: ',nelem
!     call exit(201)
            endif
!     
!     Rotation speed must be defined and positive
!     
            if(prop(index+5).le.0.d0) then
              write(*,*)'*ERROR in checkinputvaluesnet:'
              write(*,*)'       trying to define a FORCED VORTEX'
              write(*,*)
     &             '       Rotation speed n is strictly positive'
              write(*,*) '       element number: ',nelem
              call exit(201)
            endif
          endif
!     
!     Absolute/relative system
!     
        elseif((lakon(nelem)(2:4).eq.'ATR').or.
     &         (lakon(nelem)(2:4).eq.'RTA')) then
          if(prop(index+1).le.0.d0) then
            write(*,*)'*ERROR in checkinputvaluesnet:'
            write(*,*)'       trying to define an element'
            write(*,*)'       TYPE= ABSOLUTE TO RELATIVE or'
            write(*,*)'       TYPE= RELATIVE TO ABSOLUTE'
            write(*,*)'       Rotation speed is strictly positive'
            write(*,*)'       element number: ',nelem
            call exit(201)
          elseif(prop(index+3).ne.0.d0) then
            if(prop(index+2).ne.0.d0) then
              write(*,*)'*ERROR in checkinputvaluesnet:'
              write(*,*)'       trying to define an element'
              write(*,*)'       TYPE= ABSOLUTE TO RELATIVE or'
              write(*,*)'       TYPE= RELATIVE TO ABSOLUTE'
              write(*,*)'       reference element has been provided'
              write(*,*)'       but tangential velocity ',
     &             'is already defined'
              write(*,*)'       element number: ',nelem
              call exit(201)
            endif
          endif
!     
!     Air Valve 
!     
        elseif(lakon(nelem)(2:5).eq.'AVLV')then
!     
          if((dabs(prop(index+1)-1.d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-1.25d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-1.5d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-2.d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-2.5d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-3.d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-4.d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-6.d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-6.5d0).gt.1.d-8).and.
     &         (dabs(prop(index+1)-8.d0).gt.1.d-8).and.
     &         nint(prop(index+3)).eq.1)then
!     
            write(*,*) '*ERROR in air_valve: '
            write(*,*) '  Specified diameter is: ',
     &           prop(index+1)
            write(*,*) '  inch(s) diameter should be a value'
            write(*,*) '  in inch in1, 1.25, 1.5, 2, 2.5, 3,' 
            write(*,*) '  4, 6, 6.5, 8,'    
            call exit(201)
!     
          endif
!     
          if((prop(index+1).gt.90d0).or.
     &         (prop(index+1).le.0)) then 
!     
            write(*,*) '*ERROR in air_valve: '
            write(*,*) '  valve opening should be a value'
            write(*,*) '  comprised between 0 grad and 90 grad'
            call exit(201)
!     
          endif
!     
        elseif((lakon(nelem)(2:5).ne.'LIPU').and.
     &         (lakon(nelem).ne.'       ').and.
     &         (lakon(nelem)(2:4).ne.'LAB').and.
     &         (lakon(nelem)(2:5).ne.'CHAR').and.
     &         (lakon(nelem)(2:6).ne.'CARBS'))then
          if((prop(index+1).lt.0.d0)) then
            write(*,*) '*ERROR in checkinputvaluesnet: section area'
            write(*,*) '       is not positive'
            write(*,*) '       element number: ',nelem
            call exit(201)
          endif               
!     
        endif
!     
      enddo
!     
      return
      end
