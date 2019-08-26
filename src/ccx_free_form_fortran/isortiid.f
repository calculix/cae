!
!     SLATEC: public domain
!
! deck isort
      subroutine isortiid (ix,iy,dy,n,kflag)
      !
      !     modified to make the same interchanges in an integer (iy)
      !     and double (dy) array!
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
      !    ISORT sorts array IX and optionally makes the same interchanges in
      !    array IY.  The array IX may be sorted in increasing order or
      !    decreasing order.  A slightly modified quicksort algorithm is used.
      !
      !    Description of Parameters
      !       IX - integer array of values to be sorted
      !       IY - integer array to be (optionally) carried along
      !       N  - number of values in integer array IX to be sorted
      !       KFLAG - control parameter
      !             =  2  means sort IX in increasing order and carry IY along.
      !             =  1  means sort IX in increasing order (ignoring IY)
      !             = -1  means sort IX in decreasing order (ignoring IY)
      !             = -2  means sort IX in decreasing order and carry IY along.
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
      !    901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
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
      !
      implicit none
      !      .. Scalar Arguments ..
      integer kflag, n
      !      .. Array Arguments ..
      integer ix(*)
      real*8 dy(*)
      integer iy(*)
      !      .. Local Scalars ..
      real r
      integer i, ij, j, k, kk, l, m, nn, t, tt
      real*8 tty,ty
      integer uuy,uy
      !      .. Local Arrays ..
      integer il(31), iu(31)
      !      .. External Subroutines ..
      !      EXTERNAL XERMSG
      !      .. Intrinsic Functions ..
      intrinsic abs, int
      ! ***FIRST EXECUTABLE STATEMENT  ISORT
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
      !      Alter array IX to get decreasing order if needed
      !
      if (kflag .le. -1) then
         do 10 i=1,nn
            ix(i) = -ix(i)
   10    continue
      endif
      !
      if (kk .eq. 2) go to 100
      !
      !      Sort IX only
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
      t = ix(ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (ix(i) .gt. t) then
         ix(ij) = ix(i)
         ix(i) = t
         t = ix(ij)
      endif
      l = j
      !
      !      If last element of array is less than than T, interchange with T
      !
      if (ix(j) .lt. t) then
         ix(ij) = ix(j)
         ix(j) = t
         t = ix(ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (ix(i) .gt. t) then
            ix(ij) = ix(i)
            ix(i) = t
            t = ix(ij)
         endif
      endif
   !
   !      Find an element in the second half of the array which is smaller
   !      than T
   !
   40 l = l-1
      if (ix(l) .gt. t) go to 40
   !
   !      Find an element in the first half of the array which is greater
   !      than T
   !
   50 k = k+1
      if (ix(k) .lt. t) go to 50
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = ix(l)
         ix(l) = ix(k)
         ix(k) = tt
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
      t = ix(i+1)
      if (ix(i) .le. t) go to 80
      k = i
   !
   90 ix(k+1) = ix(k)
      k = k-1
      if (t .lt. ix(k)) go to 90
      ix(k+1) = t
      go to 80
  !
  !      Sort IX and carry IY along
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
      t = ix(ij)
      ty = dy(ij)
      uy = iy(ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (ix(i) .gt. t) then
         ix(ij) = ix(i)
         ix(i) = t
         t = ix(ij)
         dy(ij) = dy(i)
         iy(ij) = iy(i)
         dy(i) = ty
         iy(i) = uy
         ty = dy(ij)
         uy = iy(ij)
      endif
      l = j
      !
      !      If last element of array is less than T, interchange with T
      !
      if (ix(j) .lt. t) then
         ix(ij) = ix(j)
         ix(j) = t
         t = ix(ij)
         dy(ij) = dy(j)
         iy(ij) = iy(j)
         dy(j) = ty
         iy(j) = uy
         ty = dy(ij)
         uy = iy(ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (ix(i) .gt. t) then
            ix(ij) = ix(i)
            ix(i) = t
            t = ix(ij)
            dy(ij) = dy(i)
            iy(ij) = iy(i)
            dy(i) = ty
            iy(i) = uy
            ty = dy(ij)
            uy = iy(ij)
         endif
      endif
  !
  !      Find an element in the second half of the array which is smaller
  !      than T
  !
  130 l = l-1
      if (ix(l) .gt. t) go to 130
  !
  !      Find an element in the first half of the array which is greater
  !      than T
  !
  140 k = k+1
      if (ix(k) .lt. t) go to 140
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = ix(l)
         ix(l) = ix(k)
         ix(k) = tt
         tty = dy(l)
         uuy = iy(l)
         dy(l) = dy(k)
         iy(l) = iy(k)
         dy(k) = tty
         iy(k) = uuy
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
      t = ix(i+1)
      ty = dy(i+1)
      uy = iy(i+1)
      if (ix(i) .le. t) go to 170
      k = i
  !
  180 ix(k+1) = ix(k)
      dy(k+1) = dy(k)
      iy(k+1) = iy(k)
      k = k-1
      if (t .lt. ix(k)) go to 180
      ix(k+1) = t
      dy(k+1) = ty
      iy(k+1) = uy
      go to 170
  !
  !      Clean up
  !
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            ix(i) = -ix(i)
  200    continue
      endif
      return
      end
