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
      subroutine extrapolatefem(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,nelemload,
     &  nload,nodeboun,nboun,ndirboun,vold,ithermal,force,
     &  icfdout,ielmat,thicke,filab)
!
!     extrapolates field values at the integration points to the 
!     nodes
!
!     the C3D20RB element has 50 integration points, however, the
!     first 8 integration points coincide with those of a C3D20R
!     element. In this routine and in errorestimator.f the C3D20RB
!     element is treated as an ordinary C3D20R element
!
!     the number of internal state variables is limited to 999
!     (cfr. array field)
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*),lakonl
      character*87 filab(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,indexe,nope,
     &  nonei20(3,12),nfield,nonei10(3,6),nk,i,j,k,l,ndim,
     &  nonei15(3,9),iorienloc,iorien,ielorien(mi(3),*),konl,
     &  mint3d,m,iflag,nelemload(2,*),nload,nboun,jj,ll,
     &  nodeboun(*),ndirboun(*),ithermal(*),icfdout,ielmat(mi(3),*),
     &  nlayer,nopeexp,ilayer,kk,mint2d,nopes,kl,ki,iel,ielstart,
     &  ielemremesh,itet(4),iwedge(2,9),mesh_in_original_form
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),field(999,20*mi(3)),a8(8,8),
     &  a4(4,4),a27(20,27),a9(6,9),a2(6,2),orab(7,*),co(3,*),
     &  coords(3,27),xi,et,ze,xl(3,20),xsj,shp(4,20),weight,
     &  yiloc(6,27),a(3,3),b(3,3),c(3,3),vold(0:mi(2),*),tlayer(4),
     &  dlayer(4),xlayer(mi(3),4),thickness,xs2(3,7),xl2(3,8),
     &  xsj2(3),shp2(7,8),thicke(mi(3),*),coloc(3,8),scale,
     &  xwedge(2,2,9),a14(8,14)
!
      include "gauss.f"
!
      data coloc /-1.d0,-1.d0,-1.d0,
     &             1.d0,-1.d0,-1.d0,
     &            -1.d0,1.d0,-1.d0,
     &             1.d0,1.d0,-1.d0,
     &            -1.d0,-1.d0,1.d0,
     &             1.d0,-1.d0,1.d0,
     &            -1.d0,1.d0,1.d0,
     &             1.d0,1.d0,1.d0/
!
!     a 10-node tet is remeshed into 10 4-node tets at contact
!     interfaces; itet contains the linear tet elements to which
!     the integration points of the parent 10-node tet belong
!
      data itet /1,2,3,10/
!
      data iwedge /1,0,2,0,3,0,1,5,2,6,3,7,5,0,6,0,7,0/
!
      data xwedge /0.975615382715435242d0,0.243846172845647580d-1,
     &             0.d0,0.d0,
     &             0.975615382715435242d0,0.243846172845647580d-1,
     &             0.d0,0.d0,
     &             0.975615382715435242d0,0.243846172845647580d-1,
     &             0.d0,0.d0,
     &             0.d0,.5d0,.5d0,0.d0,
     &             0.d0,.5d0,.5d0,0.d0,
     &             0.d0,.5d0,.5d0,0.d0,
     &             0.243846172845647580d-1,0.975615382715435242d0,
     &             0.d0,0.d0,
     &             0.243846172845647580d-1,0.975615382715435242d0,
     &             0.d0,0.d0,
     &             0.243846172845647580d-1,0.975615382715435242d0,
     &             0.d0,0.d0/
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
      data a2 /  1.1455,-0.1455,1.1455,-0.1455,1.1455,-0.1455,
     &           -0.1455,1.1455,-0.1455,1.1455,-0.1455,1.1455/
      data a4 /  1.92705, -0.30902, -0.30902, -0.30902,
     &          -0.30902,  1.92705, -0.30902, -0.30902,
     &          -0.30902, -0.30902,  1.92705, -0.30902,
     &          -0.30902, -0.30902, -0.30902,  1.92705/
      data a9 / 1.63138,-0.32628,-0.32628,-0.52027, 0.10405, 0.10405,
     &         -0.32628, 1.63138,-0.32628, 0.10405,-0.52027, 0.10405,
     &         -0.32628,-0.32628, 1.63138, 0.10405, 0.10405,-0.52027,
     &          0.55556,-0.11111,-0.11111, 0.55556,-0.11111,-0.11111,
     &         -0.11111, 0.55556,-0.11111,-0.11111, 0.55556,-0.11111,
     &         -0.11111,-0.11111, 0.55556,-0.11111,-0.11111, 0.55556,
     &         -0.52027, 0.10405, 0.10405, 1.63138,-0.32628,-0.32628,
     &          0.10405,-0.52027, 0.10405,-0.32628, 1.63138,-0.32628,
     &          0.10405, 0.10405,-0.52027,-0.32628,-0.32628, 1.63138/
!
!     extrapolation from a 2x2x2=8 integration point scheme in a hex to
!     the vertex nodes
!    
      data a8 /2.549,-.683,.183,-.683,-.683,.183,
     &        -.04904,.183,-.683,2.549,-.683,.183,
     &        .183,-.683,.183,-.04904,-.683,.183,
     &        -.683,2.549,.183,-.04904,.183,-.683,
     &        .183,-.683,2.549,-.683,-.04904,.183,
     &        -.683,.183,-.683,.183,-.04904,.183,
     &        2.549,-.683,.183,-.683,.183,-.683,
     &        .183,-.04904,-.683,2.549,-.683,.183,
     &        .183,-.04904,.183,-.683,-.683,.183,
     &        -.683,2.549,-.04904,.183,-.683,.183,
     &        .183,-.683,2.549,-.683/  
!
!     extrapolation from a 14 integration point scheme in a hex to
!     the vertex nodes
!    
      data a14 /
     &  0.1396E+01,-0.3026E+00,0.1124E-01,-0.3026E+00,
     &  -0.3026E+00,0.1124E-01,0.4901E-01,0.1124E-01,
     &  -0.3026E+00,0.1396E+01,-0.3026E+00,0.1124E-01,
     &  0.1124E-01,-0.3026E+00,0.1124E-01,0.4901E-01,
     &  0.1124E-01,-0.3026E+00,0.1396E+01,-0.3026E+00,
     &  0.4901E-01,0.1124E-01,-0.3026E+00,0.1124E-01,
     &  -0.3026E+00,0.1124E-01,-0.3026E+00,0.1396E+01,
     &  0.1124E-01,0.4901E-01,0.1124E-01,-0.3026E+00,
     &  -0.3026E+00,0.1124E-01,0.4901E-01,0.1124E-01,
     &  0.1396E+01,-0.3026E+00,0.1124E-01,-0.3026E+00,
     &  0.1124E-01,-0.3026E+00,0.1124E-01,0.4901E-01,
     &  -0.3026E+00,0.1396E+01,-0.3026E+00,0.1124E-01,
     &  0.4901E-01,0.1124E-01,-0.3026E+00,0.1124E-01,
     &  0.1124E-01,-0.3026E+00,0.1396E+01,-0.3026E+00,
     &  0.1124E-01,0.4901E-01,0.1124E-01,-0.3026E+00,
     &  -0.3026E+00,0.1124E-01,-0.3026E+00,0.1396E+01,
     &  0.2069E+00,0.2069E+00,-0.6408E-01,-0.6408E-01,
     &  0.2069E+00,0.2069E+00,-0.6408E-01,-0.6408E-01,
     &  -0.6408E-01,0.2069E+00,0.2069E+00,-0.6408E-01,
     &  -0.6408E-01,0.2069E+00,0.2069E+00,-0.6408E-01,
     &  -0.6408E-01,-0.6408E-01,0.2069E+00,0.2069E+00,
     &  -0.6408E-01,-0.6408E-01,0.2069E+00,0.2069E+00,
     &  0.2069E+00,-0.6408E-01,-0.6408E-01,0.2069E+00,
     &  0.2069E+00,-0.6408E-01,-0.6408E-01,0.2069E+00,
     &  0.2069E+00,0.2069E+00,0.2069E+00,0.2069E+00,
     &  -0.6408E-01,-0.6408E-01,-0.6408E-01,-0.6408E-01,
     &  -0.6408E-01,-0.6408E-01,-0.6408E-01,-0.6408E-01,
     &  0.2069E+00,0.2069E+00,0.2069E+00,0.2069E+00/
!
!     extrapolation from a 3x3x3=27 integration point scheme in a hex to
!     the all nodes in a 20-node element
!    
      data a27 /
     &  2.37499,-0.12559,-0.16145,-0.12559,-0.12559,-0.16145, 0.11575,
     & -0.16145, 0.32628, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     & -0.10405, 0.11111, 0.32628, 0.11111,-0.10405, 0.11111,-0.31246,
     & -0.31246, 0.31481, 0.31481, 0.31481, 0.31481,-0.16902,-0.16902,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.19444, 0.15961,-0.00661,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.12559, 2.37499,
     & -0.12559,-0.16145,-0.16145,-0.12559,-0.16145, 0.11575, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111,-0.10405,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.31481,-0.16902,-0.16902, 0.31481,-0.27072,-0.19444,
     & -0.27072, 1.28439, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.22898, 0.22898, 0.22898, 0.22898, 0.05556, 0.05556, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556,-0.22222,-0.22222,
     & -0.22222,-0.22222, 0.31481,-0.31246,-0.31246, 0.31481,-0.16902,
     &  0.31481, 0.31481,-0.16902,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.19444, 0.15961,-0.00661, 0.15961,-0.27072,-0.27072,
     &  0.15961,-0.12559,-0.16145,-0.12559, 2.37499,-0.16145, 0.11575,
     & -0.16145,-0.12559, 0.11111, 0.11111, 0.32628, 0.32628,-0.10405,
     & -0.10405, 0.11111, 0.11111, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.16902,-0.16902, 0.31481,
     &  0.31481,-0.19444,-0.27072, 1.28439,-0.27072,-0.00661, 0.15961,
     & -0.19444, 0.15961, 0.15961, 0.15961,-0.27072,-0.27072,-0.16145,
     & -0.12559, 2.37499,-0.12559, 0.11575,-0.16145,-0.12559,-0.16145,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.11111,
     & -0.10405,-0.10405, 0.11111, 0.32628, 0.11111,-0.31246, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     &  1.28439,-0.19444,-0.00661,-0.19444,-0.48824,-0.48824, 0.22898,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898, 0.05556,-0.22222,
     &  0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.31481,-0.31246, 0.31481,-0.16902,
     &  0.31481,-0.31246, 0.31481,-0.16902,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.19444, 1.28439,
     & -0.19444,-0.00661,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898, 0.22898,-0.48824,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,
     & -0.29630,-0.29630,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     & -0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.31246,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     & -0.27072,-0.19444,-0.00661,-0.19444, 1.28439, 0.22898, 0.22898,
     & -0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222,
     &  0.05556, 0.05556, 0.05556, 0.05556,-0.16902, 0.31481,-0.31246,
     &  0.31481,-0.16902, 0.31481,-0.31246, 0.31481, 0.15961,-0.27072,
     & -0.27072, 0.15961, 0.15961,-0.27072,-0.27072, 0.15961,-0.00661,
     & -0.19444, 1.28439,-0.19444,-0.12559,-0.16145, 0.11575,-0.16145,
     &  2.37499,-0.12559,-0.16145,-0.12559, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.31481, 0.31481,-0.16902,-0.16902,-0.31246,
     & -0.31246, 0.31481, 0.31481,-0.19444, 0.15961,-0.00661, 0.15961,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.16145,-0.12559,-0.16145, 0.11575,-0.12559, 2.37499,
     & -0.12559,-0.16145, 0.11111, 0.11111,-0.10405,-0.10405, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     &  0.31481,-0.16902,-0.16902, 0.31481,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,-0.19444,
     & -0.27072, 1.28439,-0.27072, 0.15961, 0.15961,-0.27072, 0.22898,
     &  0.22898, 0.22898, 0.22898,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.22222,-0.22222,-0.22222,-0.22222,-0.16902, 0.31481,
     &  0.31481,-0.16902, 0.31481,-0.31246,-0.31246, 0.31481, 0.15961,
     & -0.19444, 0.15961,-0.00661,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.27072,-0.27072, 0.15961,-0.16145, 0.11575,-0.16145,
     & -0.12559,-0.12559,-0.16145,-0.12559, 2.37499,-0.10405,-0.10405,
     &  0.11111, 0.11111, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.32628,-0.16902,-0.16902, 0.31481, 0.31481,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.00661, 0.15961,-0.19444,
     &  0.15961,-0.19444,-0.27072, 1.28439,-0.27072, 0.15961, 0.15961,
     & -0.27072,-0.27072, 0.11575,-0.16145,-0.12559,-0.16145,-0.16145,
     & -0.12559, 2.37499,-0.12559,-0.10405, 0.11111, 0.11111,-0.10405,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.11111/
!
      data iflag /1/
!
      if(filab(1)(3:3).eq.'C') then
         mesh_in_original_form=0
      else
         mesh_in_original_form=1
      endif
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,nk
         do j=1,nfield
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).eq.-1) then
            cycle
         elseif((lakon(i)(5:6).eq.'IC').and.
     &          (mesh_in_original_form.eq.1)) then
            cycle
         elseif((ipkon(i).lt.-1).and.
     &          (mesh_in_original_form.eq.1)) then
!
!           for C3D20(R) elements remeshed due to contact:
!           the values at the integration points of the 
!           corresponding C3D20R element are interpolated
!           from the values at the integration points of the
!           substituting C3D8I elements
!
!           This means that a C3D20 element 
!           element is converted into a C3D20R element (having
!           8 integration points instead of 27)
!
            ielstart=kon(-2-ipkon(i)+1)-1
            if(lakon(i)(4:5).eq.'20') then
               lakon(i)(6:6)='R'
               scale=2.d0-dsqrt(3.d0)
               do j=1,8
                  iel=ielstart+j
                  xi=coloc(1,j)*scale
                  et=coloc(2,j)*scale
                  ze=coloc(3,j)*scale
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
                  do k=1,ndim
                     yi(k,j,i)=0.d0
                     do l=1,8
                        yi(k,j,i)=yi(k,j,i)+shp(4,l)*yi(k,l,iel)
                     enddo
                  enddo
               enddo
c            elseif(lakon(i)(4:5).eq.'10') then
c               do j=1,4
c                  iel=ielstart+itet(j)
c                  do k=1,ndim
c                     yi(k,j,i)=yi(k,1,iel)
c                  enddo
c               enddo
c            elseif(lakon(i)(4:5).eq.'15') then
c               do j=1,9
c                  do k=1,ndim
c                     yi(k,j,i)=0.d0
c                     do l=1,2
c                        iel=ielstart+iwedge(l,j)
c                        if(iwedge(l,j).eq.0) cycle
c                        do m=1,2
c                           yi(k,j,i)=yi(k,j,i)+xwedge(m,l,j)*
c     &                                     yi(k,m,iel)
c                        enddo
c                     enddo
c                  enddo
c               enddo
            else
               cycle
            endif
            indexe=-ipkon(i)-2
            ielemremesh=kon(indexe+1)
            kon(indexe+1)=kon(ipkon(ielemremesh)+1)
         elseif(ipkon(i).lt.-1) then
            cycle
         else
            indexe=ipkon(i)
         endif
!
c         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') then
            nlayer=0
            do j=1,mi(3)
               if(ielmat(j,i).gt.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               nopeexp=28
            elseif(lakonl(4:5).eq.'15') then
               nopeexp=21
            endif
         endif
!
!        if icfdout=1 cfd-output is requested (i.e. this routine
!        is being called from compfluid.c)
!
         if((lakonl(1:1).eq.'F').and.(icfdout.ne.1)) then
            cycle
         elseif((lakonl(1:1).ne.'F').and.(icfdout.eq.1)) then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).eq.'A'))then
            inum(kon(indexe+1))=inum(kon(indexe+1))+1
            inum(kon(indexe+2))=inum(kon(indexe+2))+1
            cycle
         else
            cycle
         endif
!
!     storage in local coordinates
!
!     calculation of the integration point coordinates for
!     output in the local system
!
         if((iorienloc.ne.0).and.
     &        ((lakonl(7:8).eq.'LC').or.(ielorien(1,i).ne.0))) then
!
            if(lakonl(7:8).ne.'LC') then
               iorien=ielorien(1,i)
            else
!     
!     composite materials
!     
               mint2d=4
               nopes=8
!     determining the number of layers
!     
               nlayer=0
               do k=1,mi(3)
                  if(ielmat(k,i).ne.0) then
                     nlayer=nlayer+1
                  endif
               enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
               indexe=ipkon(i)
               do kk=1,mint2d
                  xi=gauss3d2(1,kk)
                  et=gauss3d2(2,kk)
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  tlayer(kk)=0.d0
                  do k=1,nlayer
                     thickness=0.d0
                     do j=1,nopes
                        thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
                     enddo
                     tlayer(kk)=tlayer(kk)+thickness
                     xlayer(k,kk)=thickness
                  enddo
               enddo
!     
               ilayer=0
               do k=1,4
                  dlayer(k)=0.d0
               enddo
            endif
!     
            if(lakon(i)(4:5).eq.'8R') then
               mint3d=1
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R').or.
     &             (lakon(i)(4:6).eq.'26R')) then
               if(lakonl(7:8).eq.'LC') then
                  mint3d=8*nlayer
               else
                  mint3d=8
               endif
            elseif(lakon(i)(4:4).eq.'2') then
               mint3d=27
            elseif((lakon(i)(4:5).eq.'10').or.
     &             (lakon(i)(4:5).eq.'14')) then
               mint3d=4
            elseif(lakon(i)(4:4).eq.'4') then
               mint3d=1
            elseif(lakon(i)(4:5).eq.'15') then
               mint3d=9
            elseif(lakon(i)(4:4).eq.'6') then
               mint3d=2
            endif
!
            do j=1,nope
               konl=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl)
               enddo
            enddo
!
            do j=1,mint3d
               if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,j)
                  et=gauss3d1(2,j)
                  ze=gauss3d1(3,j)
                  weight=weight3d1(j)
               elseif((lakon(i)(4:4).eq.'8').or.
     &                (lakon(i)(4:6).eq.'20R').or.
     &                (lakon(i)(4:6).eq.'26R'))
     &                 then
                  if(lakonl(7:8).ne.'LC') then
                     xi=gauss3d2(1,j)
                     et=gauss3d2(2,j)
                     ze=gauss3d2(3,j)
                     weight=weight3d2(j)
                  else
                     kl=mod(jj,8)
                     if(kl.eq.0) kl=8
!     
                     xi=gauss3d2(1,kl)
                     et=gauss3d2(2,kl)
                     ze=gauss3d2(3,kl)
                     weight=weight3d2(kl)
!     
                     ki=mod(jj,4)
                     if(ki.eq.0) ki=4
!     
                     if(kl.eq.1) then
                        ilayer=ilayer+1
                        if(ilayer.gt.1) then
                           do k=1,4
                              dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                           enddo
                        endif
                     endif
                     ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*
     &                    xlayer(ilayer,ki))/tlayer(ki)-1.d0
                     weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!                    material and orientation
!     
                     iorien=ielorien(ilayer,i)
                  endif
               elseif(lakon(i)(4:4).eq.'2') then
                  xi=gauss3d3(1,j)
                  et=gauss3d3(2,j)
                  ze=gauss3d3(3,j)
                  weight=weight3d3(j)
               elseif((lakon(i)(4:5).eq.'10').or.
     &                (lakon(i)(4:5).eq.'14')) then
                  xi=gauss3d5(1,j)
                  et=gauss3d5(2,j)
                  ze=gauss3d5(3,j)
                  weight=weight3d5(j)
               elseif(lakon(i)(4:4).eq.'4') then
                  xi=gauss3d4(1,j)
                  et=gauss3d4(2,j)
                  ze=gauss3d4(3,j)
                  weight=weight3d4(j)
               elseif(lakon(i)(4:5).eq.'15') then
                  xi=gauss3d8(1,j)
                  et=gauss3d8(2,j)
                  ze=gauss3d8(3,j)
                  weight=weight3d8(j)
               elseif(lakon(i)(4:4).eq.'6') then
                  xi=gauss3d7(1,j)
                  et=gauss3d7(2,j)
                  ze=gauss3d7(3,j)
                  weight=weight3d7(j)
               endif
!
               if(nope.eq.20) then
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.8) then
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.10) then
                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.4) then
                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.15) then
                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
               endif
!
!              layer without orientation in a composite
!
               if(iorien.eq.0) then
                  if(nfield.eq.3) then
                     do k=1,3
                        yiloc(k,j)=yi(k,j,i)
                     enddo
                  elseif(nfield.eq.6) then
                     do k=1,6
                        yiloc(k,j)=yi(k,j,i)
                     enddo
                  endif
                  cycle
               endif
!
               do k=1,3
                  coords(k,j)=0.d0
                  do l=1,nope
                     coords(k,j)=coords(k,j)+xl(k,l)*shp(4,l)
                  enddo
               enddo
!
               if(nfield.eq.3) then
                  call transformatrix(orab(1,iorien),coords(1,j),a)
                  yiloc(1,j)=yi(1,j,i)*a(1,1)+yi(2,j,i)*a(2,1)+
     &                     yi(3,j,i)*a(3,1)
                  yiloc(2,j)=yi(1,j,i)*a(1,2)+yi(2,j,i)*a(2,2)+
     &                     yi(3,j,i)*a(3,2)
                  yiloc(3,j)=yi(1,j,i)*a(1,3)+yi(2,j,i)*a(2,3)+
     &                     yi(3,j,i)*a(3,3)
               elseif(nfield.eq.6) then
                  call transformatrix(orab(1,iorien),coords(1,j),a)
                  b(1,1)=yi(1,j,i)
                  b(2,2)=yi(2,j,i)
                  b(3,3)=yi(3,j,i)
                  b(1,2)=yi(4,j,i)
                  b(1,3)=yi(5,j,i)
                  b(2,3)=yi(6,j,i)
                  b(2,1)=b(1,2)
                  b(3,1)=b(1,3)
                  b(3,2)=b(2,3)
                  do k=1,3
                     do l=1,3
                        c(k,l)=0.d0
                        do m=1,3
                           c(k,l)=c(k,l)+b(k,m)*a(m,l)
                        enddo
                     enddo
                  enddo
                  do k=1,3
                     do l=k,3
                        b(k,l)=0.d0
                        do m=1,3
                           b(k,l)=b(k,l)+a(m,k)*c(m,l)
                        enddo
                     enddo
                  enddo
                  yiloc(1,j)=b(1,1)
                  yiloc(2,j)=b(2,2)
                  yiloc(3,j)=b(3,3)
                  yiloc(4,j)=b(1,2)
                  yiloc(5,j)=b(1,3)
                  yiloc(6,j)=b(2,3)
               endif
            enddo
!
c     Bernhardi start
            if((lakonl(4:6).eq.'20R').or.(lakonl(4:6).eq.'26R').or.
     &         (lakonl(4:5).eq.'8 ').or.(lakonl(4:5).eq.'8I')) then
c     Bernhardi end
               if(lakonl(7:8).ne.'LC') then
                  do j=1,8
                     do k=1,nfield
                        field(k,j)=0.d0
                        do l=1,8
                           field(k,j)=field(k,j)+a8(j,l)*yiloc(k,l)
                        enddo
                     enddo
                  enddo
               else
                  do m=1,nlayer
                     jj=20*(m-1)
                     ll=8*(m-1)
                     do j=1,8
                        do k=1,nfield
                           field(k,jj+j)=0.d0
                           do l=1,8
                              field(k,jj+j)=
     &                           field(k,jj+j)+a8(j,l)*yiloc(k,ll+l)
                           enddo
                        enddo
                     enddo
                  enddo
               endif
            elseif(lakonl(4:4).eq.'8') then
               do j=1,8
                  do k=1,nfield
                     field(k,j)=yiloc(k,1)
                  enddo
               enddo
            elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,4
                        field(k,j)=field(k,j)+a4(j,l)*yiloc(k,l)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'2') then
               do j=1,20
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,27
                        field(k,j)=field(k,j)+a27(j,l)*yiloc(k,l)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'4') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=yiloc(k,1)
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'1') then
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,9
                        field(k,j)=field(k,j)+a9(j,l)*yiloc(k,l)
                     enddo
                  enddo
               enddo
            else
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,2
                        field(k,j)=field(k,j)+a2(j,l)*yiloc(k,l)
                     enddo
                  enddo
               enddo
            endif
         else
!
!        storage in global coordinates
!
!        determining the field values in the vertex nodes
!        for C3D20R and C3D8: trilinear extrapolation (= use of the
!                             C3D8 linear interpolation functions)
!        for C3D8R: constant field value in each element
!        for C3D10: use of the C3D4 linear interpolation functions
!        for C3D4: constant field value in each element
!        for C3D15: use of the C3D6 linear interpolation functions
!        for C3D6: use of a linear interpolation function
!
c     Bernhardi start
            if((lakonl(4:6).eq.'20R').or.(lakonl(4:6).eq.'26R').or.
     &         (lakonl(4:5).eq.'8 ').or.(lakonl(4:5).eq.'8I')) then
c     Bernhardi end
               if(lakonl(7:8).ne.'LC') then
                  do j=1,8
                     do k=1,nfield
                        field(k,j)=0.d0
                        do l=1,8
                           field(k,j)=field(k,j)+a8(j,l)*yi(k,l,i)
                        enddo
                     enddo
                  enddo
               else
                  do m=1,nlayer
                     jj=20*(m-1)
                     ll=8*(m-1)
                     do j=1,8
                        do k=1,nfield
                           field(k,jj+j)=0.d0
                           do l=1,8
                              field(k,jj+j)=
     &                           field(k,jj+j)+a8(j,l)*yi(k,ll+l,i)
                           enddo
                        enddo
                     enddo
                  enddo
               endif
            elseif(lakonl(4:4).eq.'8') then
               do j=1,8
                  do k=1,nfield
                     field(k,j)=yi(k,1,i)
                  enddo
               enddo
            elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,4
                        field(k,j)=field(k,j)+a4(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'2') then
               do j=1,20
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,27
                        field(k,j)=field(k,j)+a27(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'4') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=yi(k,1,i)
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'1') then
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,9
                        field(k,j)=field(k,j)+a9(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            else
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,2
                        field(k,j)=field(k,j)+a2(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            endif
         endif
!
!        determining the field values in the midside nodes
!
         if((lakonl(4:6).eq.'20R').or.(lakonl(4:6).eq.'26R')) then
            if(lakonl(7:8).ne.'LC') then
               do j=9,20
                  do k=1,nfield
                     field(k,j)=(field(k,nonei20(2,j-8))+
     &                    field(k,nonei20(3,j-8)))/2.d0
                  enddo
               enddo
            else
               do m=1,nlayer
                  jj=20*(m-1)
                  do j=9,20
                     do k=1,nfield
                        field(k,jj+j)=(field(k,jj+nonei20(2,j-8))
     &                     +field(k,jj+nonei20(3,j-8)))/2.d0
                     enddo
                  enddo
               enddo
            endif
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
            do j=5,10
               do k=1,nfield
                  field(k,j)=(field(k,nonei10(2,j-4))+
     &                 field(k,nonei10(3,j-4)))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            do j=7,15
               do k=1,nfield
                  field(k,j)=(field(k,nonei15(2,j-6))+
     &                 field(k,nonei15(3,j-6)))/2.d0
               enddo
            enddo
         endif
!
!        transferring the field values into yn
!
         if(lakonl(7:8).ne.'LC') then
            do j=1,nope
               do k=1,nfield
                  yn(k,kon(indexe+j))=yn(k,kon(indexe+j))+
     &                 field(k,j)
               enddo
               inum(kon(indexe+j))=inum(kon(indexe+j))+1
            enddo
         else
            do j=1,nope*nlayer
               do k=1,nfield
                  yn(k,kon(indexe+nopeexp+j))=
     &            yn(k,kon(indexe+nopeexp+j))+field(k,j)
               enddo
               inum(kon(indexe+nopeexp+j))=inum(kon(indexe+nopeexp+j))+1
            enddo
         endif
!
c     Bernhardi start
c        incompatible modes elements
         if(lakonl(1:5).eq.'C3D8I') then
            do j=1,3
               do k=1,nfield
                  yn(k,kon(indexe+nope+j))=0.0d0
               enddo
cccc               inum(kon(indexe+nope+j))=inum(kon(indexe+nope+j))+1
            enddo
         endif
c     Bernhardi end
!
         if((ipkon(i).lt.-1).and.
     &      (mesh_in_original_form.eq.1)) kon(indexe+1)=ielemremesh
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!
!     for 1d and 2d elements only:
!     finding the solution in the original nodes
!
      if((cflag.ne.' ').and.(cflag.ne.'E')) then
         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,ne,cflag,co,
     &         vold,force,mi)
      endif
!
      return
      end
