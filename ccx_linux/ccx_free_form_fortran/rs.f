      !
      !      EISPACK -> SLATEC: public domain (cf. gams.nist.gov)
      !
      real*8 function pythag(a,b)
      implicit none
      real*8 a,b
      !
      !      finds dsqrt(a**2+b**2) without overflow or destructive underflow
      !
      real*8 p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
!
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
      !
      implicit none
      integer n,nm,ierr,matz
      real*8 a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!      this subroutine calls the recommended sequence of
!      subroutines from the eigensystem subroutine package (eispack)
      !      to find the eigenvalues and eigenvectors (if desired)
      !      of a real symmetric matrix.
      !
      !      on input
      !
      !         nm  must be set to the row dimension of the two-dimensional
      !         array parameters as declared in the calling program
      !         dimension statement.
      !
      !         n  is the order of the matrix  a.
      !
      !         a  contains the real symmetric matrix.
      !
      !         matz  is an integer variable set equal to zero if
      !         only eigenvalues are desired.  otherwise it is set to
      !         any non-zero integer for both eigenvalues and eigenvectors.
      !
      !      on output
      !
      !         w  contains the eigenvalues in ascending order.
      !
      !         z  contains the eigenvectors if matz is not zero.
      !
      !         ierr  is an integer output variable set equal to an error
      !            completion code described in the documentation for tqlrat
      !            and tql2.  the normal completion code is zero.
      !
      !         fv1  and  fv2  are temporary storage arrays.
      !
      !      questions and comments should be directed to burton s. garbow,
      !      mathematics and computer science div, argonne national laboratory
      !
      !      this version dated august 1983.
      !
      !      ------------------------------------------------------------------
      !
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
   !
   10 if (matz .ne. 0) go to 20
      !      .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      !   tqlrat encounters catastrophic underflow on the Vax
      !      call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
   !      .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
!
      subroutine tql1(n,d,e,ierr)
      !
      implicit none
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      real*8 d(n),e(n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!      this subroutine is a translation of the algol procedure tql1,
!      num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!      wilkinson.
!      handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!      this subroutine finds the eigenvalues of a symmetric
      !      tridiagonal matrix by the ql method.
      !
      !      on input
      !
      !         n is the order of the matrix.
      !
      !         d contains the diagonal elements of the input matrix.
      !
      !         e contains the subdiagonal elements of the input matrix
      !           in its last n-1 positions.  e(1) is arbitrary.
      !
      !       on output
      !
      !         d contains the eigenvalues in ascending order.  if an
      !           error exit is made, the eigenvalues are correct and
      !           ordered for indices 1,2,...ierr-1, but may not be
      !           the smallest eigenvalues.
      !
      !         e has been destroyed.
      !
      !         ierr is set to
      !           zero       for normal return,
      !           j          if the j-th eigenvalue has not been
      !                      determined after 30 iterations.
      !
      !      calls pythag for  dsqrt(a*a + b*b) .
      !
      !      questions and comments should be directed to burton s. garbow,
      !      mathematics and computer science div, argonne national laboratory
      !
      !      this version dated august 1983.
      !
      !      ------------------------------------------------------------------
      !
      ierr = 0
      if (n .eq. 1) go to 1001
      !
      do 100 i = 2, n
  100 e(i-1) = e(i)
      !
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      !
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
         !      .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  !      .......... e(n) is always zero, so there is no exit
  !                 through the bottom of the loop ..........
  110    continue
  !
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
         !      .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         !
         do 140 i = l2, n
  140    d(i) = d(i) - h
  !
  145    f = f + h
         !      .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         !      .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
         !
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
         !      .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
         !      .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
  !
  250    i = 1
  270    d(i) = p
  290 continue
      !
      go to 1001
 !      .......... set error -- no convergence to an
 !                 eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
!
      subroutine tql2(nm,n,d,e,z,ierr)
      !
      implicit none
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!      this subroutine is a translation of the algol procedure tql2,
!      num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!      wilkinson.
!      handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!      this subroutine finds the eigenvalues and eigenvectors
      !      of a symmetric tridiagonal matrix by the ql method.
      !      the eigenvectors of a full symmetric matrix can also
      !      be found if  tred2  has been used to reduce this
      !      full matrix to tridiagonal form.
      !
      !      on input
      !
      !         nm must be set to the row dimension of two-dimensional
      !           array parameters as declared in the calling program
      !           dimension statement.
      !
      !         n is the order of the matrix.
      !
      !         d contains the diagonal elements of the input matrix.
      !
      !         e contains the subdiagonal elements of the input matrix
      !           in its last n-1 positions.  e(1) is arbitrary.
      !
      !         z contains the transformation matrix produced in the
      !           reduction by  tred2, if performed.  if the eigenvectors
      !           of the tridiagonal matrix are desired, z must contain
      !           the identity matrix.
      !
      !       on output
      !
      !         d contains the eigenvalues in ascending order.  if an
      !           error exit is made, the eigenvalues are correct but
      !           unordered for indices 1,2,...,ierr-1.
      !
      !         e has been destroyed.
      !
      !         z contains orthonormal eigenvectors of the symmetric
      !           tridiagonal (or full) matrix.  if an error exit is made,
      !           z contains the eigenvectors associated with the stored
      !           eigenvalues.
      !
      !         ierr is set to
      !           zero       for normal return,
      !           j          if the j-th eigenvalue has not been
      !                      determined after 30 iterations.
      !
      !      calls pythag for  dsqrt(a*a + b*b) .
      !
      !      questions and comments should be directed to burton s. garbow,
      !      mathematics and computer science div, argonne national laboratory
      !
      !      this version dated august 1983.
      !
      !      ------------------------------------------------------------------
      !
      ierr = 0
      if (n .eq. 1) go to 1001
      !
      do 100 i = 2, n
  100 e(i-1) = e(i)
      !
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      !
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
         !      .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  !      .......... e(n) is always zero, so there is no exit
  !                 through the bottom of the loop ..........
  110    continue
  !
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
         !      .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         !
         do 140 i = l2, n
  140    d(i) = d(i) - h
  !
  145    f = f + h
         !      .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         !      .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            !      .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
  !
  200    continue
         !
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
      !      .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         !
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
         !
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         !
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
  !
  300 continue
      !
      go to 1001
 !      .......... set error -- no convergence to an
 !                 eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
!
      subroutine tred1(nm,n,a,d,e,e2)
      !
      implicit none
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),e2(n)
      real*8 f,g,h,scale
!
!      this subroutine is a translation of the algol procedure tred1,
!      num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!      handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!      this subroutine reduces a real symmetric matrix
      !      to a symmetric tridiagonal matrix using
      !      orthogonal similarity transformations.
      !
      !      on input
      !
      !         nm must be set to the row dimension of two-dimensional
      !           array parameters as declared in the calling program
      !           dimension statement.
      !
      !         n is the order of the matrix.
      !
      !         a contains the real symmetric input matrix.  only the
      !           lower triangle of the matrix need be supplied.
      !
      !      on output
      !
      !         a contains information about the orthogonal trans-
      !           formations used in the reduction in its strict lower
      !           triangle.  the full upper triangle of a is unaltered.
      !
      !         d contains the diagonal elements of the tridiagonal matrix.
      !
      !         e contains the subdiagonal elements of the tridiagonal
      !           matrix in its last n-1 positions.  e(1) is set to zero.
      !
      !         e2 contains the squares of the corresponding elements of e.
      !           e2 may coincide with e if the squares are not needed.
      !
      !      questions and comments should be directed to burton s. garbow,
      !      mathematics and computer science div, argonne national laboratory
      !
      !      this version dated august 1983.
      !
      !      ------------------------------------------------------------------
      !
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
      !      .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
         !      .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
         !
         if (scale .ne. 0.0d0) go to 140
         !
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
  !
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
  !
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
         !
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
         !      .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
         !
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
            !
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
  !
  220       e(j) = g
  240    continue
         !      .......... form p ..........
         f = 0.0d0
         !
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
         !
         h = f / (h + h)
         !      .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
         !      .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
            !
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
  !
  280    continue
  !
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
  !
  300 continue
      !
      return
      end
!
      subroutine tred2(nm,n,a,d,e,z)
      !
      implicit none
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),z(nm,n)
      real*8 f,g,h,hh,scale
!
!      this subroutine is a translation of the algol procedure tred2,
!      num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!      handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!      this subroutine reduces a real symmetric matrix to a
      !      symmetric tridiagonal matrix using and accumulating
      !      orthogonal similarity transformations.
      !
      !      on input
      !
      !         nm must be set to the row dimension of two-dimensional
      !           array parameters as declared in the calling program
      !           dimension statement.
      !
      !         n is the order of the matrix.
      !
      !         a contains the real symmetric input matrix.  only the
      !           lower triangle of the matrix need be supplied.
      !
      !      on output
      !
      !         d contains the diagonal elements of the tridiagonal matrix.
      !
      !         e contains the subdiagonal elements of the tridiagonal
      !           matrix in its last n-1 positions.  e(1) is set to zero.
      !
      !         z contains the orthogonal transformation matrix
      !           produced in the reduction.
      !
      !         a and z may coincide.  if distinct, a is unaltered.
      !
      !      questions and comments should be directed to burton s. garbow,
      !      mathematics and computer science div, argonne national laboratory
      !
      !      this version dated august 1983.
      !
      !      ------------------------------------------------------------------
      !
      do 100 i = 1, n
         !
         do 80 j = i, n
   80    z(j,i) = a(j,i)
         !
         d(i) = a(n,i)
  100 continue
      !
      if (n .eq. 1) go to 510
      !      .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
         !      .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
         !
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
         !
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
         !
         go to 290
  !
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
         !
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         !      .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
         !
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
            !
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
  !
  220       e(j) = g
  240    continue
         !      .......... form p ..........
         f = 0.0d0
         !
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
         !
         hh = f / (h + h)
         !      .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
         !      .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
            !
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
            !
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
  !
  290    d(i) = h
  300 continue
      !      .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
         !
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
         !
         do 360 j = 1, l
            g = 0.0d0
            !
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
            !
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
  !
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
  !
  500 continue
  !
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
      !
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
