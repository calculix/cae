!
!     SLATEC: public domain
!
! deck isort
      subroutine isortiiddc(ix1,ix2,dy1,dy2,cy,n,kflag)
      !
      !     modified to make the same interchanges in an integer (ix2), two
      !     double (dy1 and dy2) and a char*20 aray (cy)
      !
      ! ***BEGIN PROLOGUE  ISORT
      ! ***PURPOSE  Sort an array and optionally make the same interchanges in
      !             an auxiliary array.  The array may be sorted in increasing
      !             or decreasing order.  A slightly modified QUICKSORT
      !             algorithm is used.
      ! ***LIBRARY   SLATEC
      ! ***CATEGORY  N6A2A
      ! ***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
      ! ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
      ! ***AUTHOR  Jones, R. E., (SNLA)
      !            Kahaner, D. K., (NBS)
      !            Wisniewski, J. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    ISORT sorts array IX1 and optionally makes the same interchanges in
      !    array IY.  The array IX1 may be sorted in increasing order or
      !    decreasing order.  A slightly modified quicksort algorithm is used.
      !
      !    Description of Parameters
      !       IX1 - integer array of values to be sorted
      !       IY - integer array to be (optionally) carried along
      !       N  - number of values in integer array IX1 to be sorted
      !       KFLAG - control parameter
      !             =  2  means sort IX1 in increasing order and carry IY along.
      !             =  1  means sort IX1 in increasing order (ignoring IY)
      !             = -1  means sort IX1 in decreasing order (ignoring IY)
      !             = -2  means sort IX1 in decreasing order and carry IY along.
      !
      ! ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
      !                  for sorting with minimal storage, Communications of
      !                  the ACM, 12, 3 (1969), pp. 185-187.
      ! ***ROUTINES CALLED  XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    761118  DATE WRITTEN
      !    810801  Modified by David K. Kahaner.
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891009  Removed unreferenced statement labels.  (WRB)
      !    891009  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
      !    901012  Declared all variables; changed X,Y to IX1,IY. (M. McClain)
      !    920501  Reformatted the REFERENCES section.  (DWL, WRB)
      !    920519  Clarified error messages.  (DWL)
      !    920801  Declarations section rebuilt and code restructured to use
      !            IF-THEN-ELSE-ENDIF.  (RWC, WRB)
      !   100411  changed the dimension of IL and IU from 21 to 31.
      !
      !     field IL and IU have the dimension 31. This is log2 of the largest
      !     array size to be sorted. If arrays larger than 2**31 in length have
      !     to be sorted, this dimension has to be modified accordingly
      !
      ! ***END PROLOGUE  ISORT
      !      .. Scalar Arguments ..
      implicit none
      !
      integer kflag, n,iside,istat
      !      .. Array Arguments ..
      integer ix1(2,*),ix2(2,*)
      real*8 dy1(2,*),dy2(2,*)
      character*20 cy(*)
      !      .. Local Scalars ..
      real r
      integer i, ij, j, k, kk, l, m, nn, t, tt,tx21,tx12,tx22,&
        ttx21,ttx12,ttx22
      real*8 tty11,tty12,ty11,ty12,tty21,tty22,ty21,ty22
      character*20 uuy,uy
      !      .. Local Arrays ..
      integer il(31), iu(31)
      !      .. External Subroutines ..
      !      EXTERNAL XERMSG
      !      .. Intrinsic Functions ..
      intrinsic abs, int
      ! ***FIRST EXECUTABLE STATEMENT  ISORT
      !
      do i=1,n
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix1(1,i)=10*ix1(1,i)+iside
      enddo
      !
      nn = n
      if (nn .lt. 1) then
         !         CALL XERMSG ('SLATEC', 'ISORT',
         !     +      'The number of values to be sorted is not positive.', 1, 1)
         return
      endif
      !
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
         !         CALL XERMSG ('SLATEC', 'ISORT',
         !     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
         !     +      1)
         return
      endif
      !
      !      Alter array IX1 to get decreasing order if needed
      !
      if (kflag .le. -1) then
         do 10 i=1,nn
            ix1(1,i) = -ix1(1,i)
   10    continue
      endif
      !
      if (kk .eq. 2) go to 100
      !
      !      Sort IX1 only
      !
      m = 1
      i = 1
      j = nn
      r = 0.375e0
   !
   20 if (i .eq. j) go to 60
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
   !
   30 k = i
      !
      !      Select a central element of the array and save it in location T
      !
      ij = i + int((j-i)*r)
      t = ix1(1,ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (ix1(1,i) .gt. t) then
         ix1(1,ij) = ix1(1,i)
         ix1(1,i) = t
         t = ix1(1,ij)
      endif
      l = j
      !
      !      If last element of array is less than than T, interchange with T
      !
      if (ix1(1,j) .lt. t) then
         ix1(1,ij) = ix1(1,j)
         ix1(1,j) = t
         t = ix1(1,ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (ix1(1,i) .gt. t) then
            ix1(1,ij) = ix1(1,i)
            ix1(1,i) = t
            t = ix1(1,ij)
         endif
      endif
   !
   !      Find an element in the second half of the array which is smaller
   !      than T
   !
   40 l = l-1
      if (ix1(1,l) .gt. t) go to 40
   !
   !      Find an element in the first half of the array which is greater
   !      than T
   !
   50 k = k+1
      if (ix1(1,k) .lt. t) go to 50
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = ix1(1,l)
         ix1(1,l) = ix1(1,k)
         ix1(1,k) = tt
         go to 40
      endif
      !
      !      Save upper and lower subscripts of the array yet to be sorted
      !
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 70
   !
   !      Begin again on another portion of the unsorted array
   !
   60 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
   !
   70 if (j-i .ge. 1) go to 30
      if (i .eq. 1) go to 20
      i = i-1
   !
   80 i = i+1
      if (i .eq. j) go to 60
      t = ix1(1,i+1)
      if (ix1(1,i) .le. t) go to 80
      k = i
   !
   90 ix1(1,k+1) = ix1(1,k)
      k = k-1
      if (t .lt. ix1(1,k)) go to 90
      ix1(1,k+1) = t
      go to 80
  !
  !      Sort IX1 and carry IY along
  !
  100 m = 1
      i = 1
      j = nn
      r = 0.375e0
  !
  110 if (i .eq. j) go to 150
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
  !
  120 k = i
      !
      !      Select a central element of the array and save it in location T
      !
      ij = i + int((j-i)*r)
      t = ix1(1,ij)
      ty11 = dy1(1,ij)
      ty21 = dy1(2,ij)
      ty12 = dy2(1,ij)
      ty22 = dy2(2,ij)
      tx21 = ix1(2,ij)
      tx12=ix2(1,ij)
      tx22=ix2(2,ij)
      uy = cy(ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (ix1(1,i) .gt. t) then
         ix1(1,ij) = ix1(1,i)
         ix1(1,i) = t
         t = ix1(1,ij)
         dy1(1,ij) = dy1(1,i)
         dy1(2,ij) = dy1(2,i)
         dy2(1,ij) = dy2(1,i)
         dy2(2,ij) = dy2(2,i)
         ix1(2,ij) = ix1(2,i)
         ix2(1,ij)=ix2(1,i)
         ix2(2,ij)=ix2(2,i)
         cy(ij) = cy(i)
         dy1(1,i) = ty11
         dy1(2,i) = ty21
         dy2(1,i) = ty12
         dy2(2,i) = ty22
         ix1(2,i) = tx21
         ix2(1,i)=tx12
         ix2(2,i)=tx22
         cy(i) = uy
         ty11 = dy1(1,ij)
         ty21 = dy1(2,ij)
         ty12 = dy2(1,ij)
         ty22 = dy2(2,ij)
         tx21 = ix1(2,ij)
         tx12=ix2(1,ij)
         tx22=ix2(2,ij)
         uy = cy(ij)
      endif
      l = j
      !
      !      If last element of array is less than T, interchange with T
      !
      if (ix1(1,j) .lt. t) then
         ix1(1,ij) = ix1(1,j)
         ix1(1,j) = t
         t = ix1(1,ij)
         dy1(1,ij) = dy1(1,j)
         dy1(2,ij) = dy1(2,j)
         dy2(1,ij) = dy2(1,j)
         dy2(2,ij) = dy2(2,j)
         ix1(2,ij) = ix1(2,j)
         ix2(1,ij)=ix2(1,j)
         ix2(2,ij)=ix2(2,j)
         cy(ij) = cy(j)
         dy1(1,j) = ty11
         dy1(2,j) = ty21
         dy2(1,j) = ty12
         dy2(2,j) = ty22
         ix1(2,j) = tx21
         ix2(1,j)=tx12
         ix2(2,j)=tx22
         cy(j) = uy
         ty11 = dy1(1,ij)
         ty21 = dy1(2,ij)
         ty12 = dy2(1,ij)
         ty22 = dy2(2,ij)
         tx21 = ix1(2,ij)
         tx12=ix2(1,ij)
         tx22=ix2(2,ij)
         uy = cy(ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (ix1(1,i) .gt. t) then
            ix1(1,ij) = ix1(1,i)
            ix1(1,i) = t
            t = ix1(1,ij)
            dy1(1,ij) = dy1(1,i)
            dy1(2,ij) = dy1(2,i)
            dy2(1,ij) = dy2(1,i)
            dy2(2,ij) = dy2(2,i)
            ix1(2,ij) = ix1(2,i)
            ix2(1,ij)=ix2(1,i)
            ix2(2,ij)=ix2(2,i)
            cy(ij) = cy(i)
            dy1(1,i) = ty11
            dy1(2,i) = ty21
            dy2(1,i) = ty12
            dy2(2,i) = ty22
            ix1(2,i) = tx21
            ix2(1,i)=tx12
            ix2(2,i)=tx22
            cy(i) = uy
            ty11 = dy1(1,ij)
            ty21 = dy1(2,ij)
            ty12 = dy2(1,ij)
            ty22 = dy2(2,ij)
            tx21 = ix1(2,ij)
            tx12=ix2(1,ij)
            tx22=ix2(2,ij)
            uy = cy(ij)
         endif
      endif
  !
  !      Find an element in the second half of the array which is smaller
  !      than T
  !
  130 l = l-1
      if (ix1(1,l) .gt. t) go to 130
  !
  !      Find an element in the first half of the array which is greater
  !      than T
  !
  140 k = k+1
      if (ix1(1,k) .lt. t) go to 140
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = ix1(1,l)
         ix1(1,l) = ix1(1,k)
         ix1(1,k) = tt
         tty11 = dy1(1,l)
         tty21 = dy1(2,l)
         tty12 = dy2(1,l)
         tty22 = dy2(2,l)
         ttx21 = ix1(2,l)
         ttx12=ix2(1,l)
         ttx22=ix2(2,l)
         uuy = cy(l)
         dy1(1,l) = dy1(1,k)
         dy1(2,l) = dy1(2,k)
         dy2(1,l) = dy2(1,k)
         dy2(2,l) = dy2(2,k)
         ix1(2,l) = ix1(2,k)
         ix2(1,l)=ix2(1,k)
         ix2(2,l)=ix2(2,k)
         cy(l) = cy(k)
         dy1(1,k) = tty11
         dy1(2,k) = tty21
         dy2(1,k) = tty12
         dy2(2,k) = tty22
         ix1(2,k) = ttx21
         ix2(1,k)=ttx12
         ix2(2,k)=ttx22
         cy(k) = uuy
         go to 130
      endif
      !
      !      Save upper and lower subscripts of the array yet to be sorted
      !
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 160
  !
  !      Begin again on another portion of the unsorted array
  !
  150 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
  !
  160 if (j-i .ge. 1) go to 120
      if (i .eq. 1) go to 110
      i = i-1
  !
  170 i = i+1
      if (i .eq. j) go to 150
      t = ix1(1,i+1)
      ty11 = dy1(1,i+1)
      ty21 = dy1(2,i+1)
      ty12 = dy2(1,i+1)
      ty22 = dy2(2,i+1)
      tx21 = ix1(2,i+1)
      tx12=ix2(1,i+1)
      tx22=ix2(2,i+1)
      uy = cy(i+1)
      if (ix1(1,i) .le. t) go to 170
      k = i
  !
  180 ix1(1,k+1) = ix1(1,k)
      dy1(1,k+1) = dy1(1,k)
      dy1(2,k+1) = dy1(2,k)
      dy2(1,k+1) = dy2(1,k)
      dy2(2,k+1) = dy2(2,k)
      ix1(2,k+1) = ix1(2,k)
      ix2(1,k+1)=ix2(1,k)
      ix2(2,k+1)=ix2(2,k)
      cy(k+1) = cy(k)
      k = k-1
      if (t .lt. ix1(1,k)) go to 180
      ix1(1,k+1) = t
      dy1(1,k+1) = ty11
      dy1(2,k+1) = ty21
      dy2(1,k+1) = ty12
      dy2(2,k+1) = ty22
      ix1(2,k+1) = tx21
      ix2(1,k+1)=tx12
      ix2(2,k+1)=tx22
      cy(k+1) = uy
      go to 170
  !
  !      Clean up
  !
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            ix1(1,i) = -ix1(1,i)
  200    continue
      endif
      !
      do i=1,nn
         read(cy(i)(2:2),'(i1)',iostat=istat) iside 
         if(istat.gt.0) iside=0
         ix1(1,i)=(ix1(1,i)-iside)/10
      enddo
      !
      return
      end
