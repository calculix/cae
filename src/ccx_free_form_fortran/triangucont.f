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
      subroutine triangucont(ncont,ntie,tieset,nset,set,istartset,&
        iendset,ialset,itietri,lakon,ipkon,kon,koncont,kind1,kind2,&
        co,nk)
      !
      !     generate a triangulation of the contact master surfaces
      !
      implicit none
      !
      character*1 kind1,kind2,c
      character*3 m1,m2,m3
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*81 tieset(3,*),mastset,set(*)
      character*88 fntria
      !
      integer ncont,ntie,i,j,k,l,nset,istartset(*),iendset(*),&
        ialset(*),itrifac9(3,8),itrifac7(3,6),&
        iright,itietri(2,ntie),nelem,jface,indexe,ipkon(*),nope,m,one,&
        ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,ilen,&
        ntrifac,itrifac3(3,1),itrifac4(3,2),itrifac6(3,4),itrifac8(3,6),&
        itrifac(3,8),nnodelem,nface,nodef(9),kon(*),koncont(4,*),nk,&
        ncontini
      !
      real*8 co(3,*)
      !
      !     nodes per face for hex elements
      !
      data ifaceq /4,3,2,1,11,10,9,12,&
                  5,6,7,8,13,14,15,16,&
                  1,2,6,5,9,18,13,17,&
                  2,3,7,6,10,19,14,18,&
                  3,4,8,7,11,20,15,19,&
                  4,1,5,8,12,17,16,20/
      !
      !     nodes per face for tet elements
      !
      data ifacet /1,3,2,7,6,5,&
                   1,2,4,5,9,8,&
                   2,3,4,6,10,9,&
                   1,4,3,8,10,7/
      !
      !     nodes per face for linear wedge elements
      !
      data ifacew1 /1,3,2,0,&
                   4,5,6,0,&
                   1,2,5,4,&
                   2,3,6,5,&
                   3,1,4,6/
      !
      !     nodes per face for quadratic wedge elements
      !
      data ifacew2 /1,3,2,9,8,7,0,0,&
                   4,5,6,10,11,12,0,0,&
                   1,2,5,4,7,14,10,13,&
                   2,3,6,5,8,15,11,14,&
                   3,1,4,6,9,13,12,15/
      !
      !     triangulation for three-node face
      !
      data itrifac3 /1,2,3/
      !
      !     triangulation for four-node face
      !
      data itrifac4 /1,2,4,2,3,4/
      !
      !     triangulation for six-node face
      !
      data itrifac6 /1,4,6,4,2,5,6,5,3,4,5,6/
      !
      !     triangulation for seven-node face
      !
      data itrifac7 /1,4,7,4,2,7,2,5,7,5,3,7,3,6,7,6,1,7/
      !
      !     triangulation for eight-node face
      !
      data itrifac8 /1,5,8,5,2,6,7,6,3,8,7,4,8,5,7,5,6,7/
      !
      !     triangulation for nine-node face
      !
      data itrifac9 /1,5,9,5,2,9,2,6,9,6,3,9,3,7,9,7,4,9,4,8,9,8,1,9/
      !
      ncont=0
      !
      do i=1,ntie
         !
         !        check for contact conditions
         !
         ncontini=ncont
         !
         if((tieset(1,i)(81:81).eq.kind1).or.&
            (tieset(1,i)(81:81).eq.kind2)) then
            mastset=tieset(3,i)
            !
            !           determining the master surface
            !
            do j=1,nset
               if(set(j).eq.mastset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in triangucont: master surface',&
                     mastset
               write(*,*) '       does not exist'
               call exit(201)
            endif
            iright=j
            !
            itietri(1,i)=ncont+1
            !
            do j=istartset(iright),iendset(iright)
               !
               nelem=int(ialset(j)/10.d0)
               jface=ialset(j)-10*nelem
               !
               indexe=ipkon(nelem)
               !
               if(lakon(nelem)(4:5).eq.'20') then
                  nnodelem=8
                  nface=6
               elseif(lakon(nelem)(4:4).eq.'8') then
                  nnodelem=4
                  nface=6
               elseif(lakon(nelem)(4:5).eq.'10') then
                  nnodelem=6
                  nface=4
               elseif(lakon(nelem)(4:4).eq.'4') then
                  nnodelem=3
                  nface=4
               elseif(lakon(nelem)(4:5).eq.'15') then
                  if(jface.le.2) then
                     nnodelem=6
                  else
                     nnodelem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelem)(4:4).eq.'6') then
                  if(jface.le.2) then
                     nnodelem=3
                  else
                     nnodelem=4
                  endif
                  nface=5
                  nope=6
               else
                  cycle
               endif
               !
               !     determining the nodes of the face
               !
               if(nface.eq.4) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacet(k,jface))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew1(k,jface))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew2(k,jface))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifaceq(k,jface))
                  enddo
               endif
               !
               !     number of triangles
               !
               if(nnodelem.eq.3) then
                  ntrifac=1
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac3(k,l)
                     enddo
                  enddo
               elseif(nnodelem.eq.4) then
                  ntrifac=2
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac4(k,l)
                     enddo
                  enddo
               elseif(nnodelem.eq.6) then
                  ntrifac=4
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac6(k,l)
                     enddo
                  enddo
               elseif(nnodelem.eq.7) then
                  ntrifac=6
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac7(k,l)
                     enddo
                  enddo
               elseif(nnodelem.eq.8) then
                  ntrifac=6
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac8(k,l)
                     enddo
                  enddo
               elseif(nnodelem.eq.9) then
                  ntrifac=8
                  do l=1,ntrifac
                     do k=1,3
                        itrifac(k,l)=itrifac9(k,l)
                     enddo
                  enddo
               endif
               !
               !     storing the topology of the triangles
               !
               do l=1,ntrifac
                  !
                  ncont=ncont+1
                  do k=1,3
                     node=nodef(itrifac(k,l))
                     koncont(k,ncont)=node
                  enddo
                  !
                  koncont(4,ncont)=ialset(j)
               !
               enddo
            !
            enddo
            !
            itietri(2,i)=ncont
         !
         endif
      !
      !     storing the triangulation in .frd format
      !
      !          ilen=index(tieset(3,i),' ')
      !          fntria(1:3)='Tri'
      !          do j=4,ilen+2
      !             fntria(j:j)=tieset(3,i)(j-3:j-3)
      !          enddo
      !          fntria(ilen+3:ilen+6)='.frd'
      !          do j=ilen+7,88
      !             fntria(j:j)=' '
      !          enddo
      ! !
      !          open(70,file=fntria,status='unknown')
      !          c='C'
      !          m1=' -1'
      !          m2=' -2'
      !          m3=' -3'
      !          p0='    0'
      !          p1='    1'
      !          p2='    2'
      !          p3='    3'
      !          p7='    7'
      !          p9999=' 9999'
      !          one=1
      !          write(70,'(a5,a1)') p1,c
      !          write(70,'(a5,a1,67x,i1)') p2,c,one
      !          do j=1,nk
      !             write(70,'(a3,i10,1p,3e12.5)') m1,j,(co(k,j),k=1,3)
      !          enddo
      !          write(70,'(a3)') m3
      !          write(70,'(a5,a1,67x,i1)') p3,c,one
      !          do j=ncontini+1,ncont
      !             write(70,'(a3,i10,2a5)')m1,j,p7,p0
      !             write(70,'(a3,3i10)') m2,(koncont(k,j),k=1,3)
      !          enddo
      !          write(70,'(a3)') m3
      !          write(70,'(a5)') p9999
      !          close(70)
      !
      enddo
      !
      return
      end

