!
!     SLATEC: public domain
!
! deck dsort
      subroutine dsort (dx, iy, n, kflag)
      !
      !     slight change: XERMSG was removed; error messages are
      !                    led to the screen;
      !
      ! ***BEGIN PROLOGUE  DSORT
      ! ***PURPOSE  Sort an array and optionally make the same interchanges in
      !             an auxiliary array.  The array may be sorted in increasing
      !             or decreasing order.  A slightly modified QUICKSORT
      !             algorithm is used.
      ! ***LIBRARY   SLATEC
      ! ***CATEGORY  N6A2B
      ! ***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
      ! ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
      ! ***AUTHOR  Jones, R. E., (SNLA)
      !            Wisniewski, J. A., (SNLA)
      ! ***ROUTINES CALLED  XERMSG
      ! ***DESCRIPTION
      !
      !    DSORT sorts array DX and optionally makes the same interchanges in
      !    array IY.  The array DX may be sorted in increasing order or
      !    decreasing order.  A slightly modified quicksort algorithm is used.
      !
      !    Description of Parameters
      !       DX - array of values to be sorted   (usually abscissas)
      !       IY - array to be (optionally) carried along
      !       N  - number of values in array DX to be sorted
      !       KFLAG - control parameter
      !             =  2  means sort DX in increasing order and carry IY along.
      !             =  1  means sort DX in increasing order (ignoring IY)
      !             = -1  means sort DX in decreasing order (ignoring IY)
      !             = -2  means sort DX in decreasing order and carry IY along.
      !
      ! ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
      !                  for sorting with minimal storage, Communications of
      !                  the ACM, 12, 3 (1969), pp. 185-187.
      ! ***REVISION HISTORY  (YYMMDD)
      !    761101  DATE WRITTEN
      !    761118  Modified to use the Singleton quicksort algorithm.  (JAW)
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891009  Removed unreferenced statement labels.  (WRB)
      !    891024  Changed category.  (WRB)
      !    891024  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
      !    901012  Declared all variables; changed X,Y to DX,IY; changed
      !            code to parallel SSORT. (M. McClain)
      !    920501  Reformatted the REFERENCES section.  (DWL, WRB)
      !    920519  Clarified error messages.  (DWL)
      !    920801  Declarations section rebuilt and code restructured to use
      !            IF-THEN-ELSE-ENDIF.  (RWC, WRB)
      !   100411  changed the dimension of IL and IU from 21 to 31.
      !   150514  inserted intent statements
      !
      !     field IL and IU have the dimension 31. This is log2 of the largest
      !     array size to be sorted. If arrays larger than 2**31 in length have
      !     to be sorted, this dimension has to be modified accordingly
      !
      ! ***END PROLOGUE  DSORT
      implicit none
      !      .. Scalar Arguments ..
      integer kflag, n,iy(*),ty,tty
      !      .. Array Arguments ..
      double precision dx(*)
      !      .. Local Scalars ..
      double precision r, t, tt
      integer i, ij, j, k, kk, l, m, nn
      !      .. Local Arrays ..
      integer il(31), iu(31)
      !      .. External Subroutines ..
      !       EXTERNAL XERMSG
      !      .. Intrinsic Functions ..
      intrinsic abs, int
      !
      intent(in) n,kflag 
      !
      intent(inout) dx,iy
      ! ***FIRST EXECUTABLE STATEMENT  DSORT
      nn = n
      if (nn .lt. 1) then
         write(*,*) '*error in dsort: the number of values to be'
         write(*,*) '       sorted is not positive'
         call exit(201)
      endif
      !
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
         write(*,*) '*error in dsort: the sort control parameter is'
         write(*,*) '       not 2, 1, -1 or -2'
         call exit(201)
      endif
      !
      !      Alter array DX to get decreasing order if needed
      !
      if (kflag .le. -1) then
         do 10 i=1,nn
            dx(i) = -dx(i)
   10    continue
      endif
      !
      if (kk .eq. 2) go to 100
      !
      !      Sort DX only
      !
      m = 1
      i = 1
      j = nn
      r = 0.375d0
   !
   20 if (i .eq. j) go to 60
      if (r .le. 0.5898437d0) then
         r = r+3.90625d-2
      else
         r = r-0.21875d0
      endif
   !
   30 k = i
      !
      !      Select a central element of the array and save it in location T
      !
      ij = i + int((j-i)*r)
      t = dx(ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (dx(i) .gt. t) then
         dx(ij) = dx(i)
         dx(i) = t
         t = dx(ij)
      endif
      l = j
      !
      !      If last element of array is less than than T, interchange with T
      !
      if (dx(j) .lt. t) then
         dx(ij) = dx(j)
         dx(j) = t
         t = dx(ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (dx(i) .gt. t) then
            dx(ij) = dx(i)
            dx(i) = t
            t = dx(ij)
         endif
      endif
   !
   !      Find an element in the second half of the array which is smaller
   !      than T
   !
   40 l = l-1
      if (dx(l) .gt. t) go to 40
   !
   !      Find an element in the first half of the array which is greater
   !      than T
   !
   50 k = k+1
      if (dx(k) .lt. t) go to 50
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = dx(l)
         dx(l) = dx(k)
         dx(k) = tt
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
      t = dx(i+1)
      if (dx(i) .le. t) go to 80
      k = i
   !
   90 dx(k+1) = dx(k)
      k = k-1
      if (t .lt. dx(k)) go to 90
      dx(k+1) = t
      go to 80
  !
  !      Sort DX and carry IY along
  !
  100 m = 1
      i = 1
      j = nn
      r = 0.375d0
  !
  110 if (i .eq. j) go to 150
      if (r .le. 0.5898437d0) then
         r = r+3.90625d-2
      else
         r = r-0.21875d0
      endif
  !
  120 k = i
      !
      !      Select a central element of the array and save it in location T
      !
      ij = i + int((j-i)*r)
      t = dx(ij)
      ty = iy(ij)
      !
      !      If first element of array is greater than T, interchange with T
      !
      if (dx(i) .gt. t) then
         dx(ij) = dx(i)
         dx(i) = t
         t = dx(ij)
         iy(ij) = iy(i)
         iy(i) = ty
         ty = iy(ij)
      endif
      l = j
      !
      !      If last element of array is less than T, interchange with T
      !
      if (dx(j) .lt. t) then
         dx(ij) = dx(j)
         dx(j) = t
         t = dx(ij)
         iy(ij) = iy(j)
         iy(j) = ty
         ty = iy(ij)
         !
         !         If first element of array is greater than T, interchange with T
         !
         if (dx(i) .gt. t) then
            dx(ij) = dx(i)
            dx(i) = t
            t = dx(ij)
            iy(ij) = iy(i)
            iy(i) = ty
            ty = iy(ij)
         endif
      endif
  !
  !      Find an element in the second half of the array which is smaller
  !      than T
  !
  130 l = l-1
      if (dx(l) .gt. t) go to 130
  !
  !      Find an element in the first half of the array which is greater
  !      than T
  !
  140 k = k+1
      if (dx(k) .lt. t) go to 140
      !
      !      Interchange these elements
      !
      if (k .le. l) then
         tt = dx(l)
         dx(l) = dx(k)
         dx(k) = tt
         tty = iy(l)
         iy(l) = iy(k)
         iy(k) = tty
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
      t = dx(i+1)
      ty = iy(i+1)
      if (dx(i) .le. t) go to 170
      k = i
  !
  180 dx(k+1) = dx(k)
      iy(k+1) = iy(k)
      k = k-1
      if (t .lt. dx(k)) go to 180
      dx(k+1) = t
      iy(k+1) = ty
      go to 170
  !
  !      Clean up
  !
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            dx(i) = -dx(i)
  200    continue
      endif
      return
      end
