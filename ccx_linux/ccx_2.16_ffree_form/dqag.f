!
!      CMLIB - public domain
!
      subroutine dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,&
          limit,lenw,last,iwork,work,phi,lambda1,zk0,Pup,Tup,rurd,xflow,&
           kup)
      ! ***begin prologue  dqag
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a1
      ! ***keywords  automatic integrator, general-purpose,
      !              integrand examinator, globally adaptive,
      !              gauss-kronrod
      ! ***author  piessens,robert,appl. math. & progr. div - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  the routine calculates an approximation result to a given
      !             definite integral i = integral of f over (a,b),
      !             hopefully satisfying following claim for accuracy
      !             abs(i-result)le.max(epsabs,epsrel*abs(i)).
      ! ***description
      !
      !         computation of a definite integral
      !         standard fortran subroutine
      !         double precision version
      !
      !             f      - double precision
      !                      function subprogam defining the integrand
      !                      function f(x). the actual name for f needs to be
      !                      declared e x t e r n a l in the driver program.
      !
      !             a      - double precision
      !                      lower limit of integration
      !
      !             b      - double precision
      !                      upper limit of integration
      !
      !             epsabs - double precision
      !                      absolute accoracy requested
      !             epsrel - double precision
      !                      relative accuracy requested
      !                      if  epsabs.le.0
      !                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
      !                      the routine will end with ier = 6.
      !
      !             key    - integer
      !                      key for choice of local integration rule
      !                      a gauss-kronrod pair is used with
      !                        7 - 15 points if key.lt.2,
      !                       10 - 21 points if key = 2,
      !                       15 - 31 points if key = 3,
      !                       20 - 41 points if key = 4,
      !                       25 - 51 points if key = 5,
      !                       30 - 61 points if key.gt.5.
      !
      !          on return
      !             result - double precision
      !                      approximation to the integral
      !
      !             abserr - double precision
      !                      estimate of the modulus of the absolute error,
      !                      which should equal or exceed abs(i-result)
      !
      !             neval  - integer
      !                      number of integrand evaluations
      !
      !             ier    - integer
      !                      ier = 0 normal and reliable termination of the
      !                              routine. it is assumed that the requested
      !                              accuracy has been achieved.
      !                      ier.gt.0 abnormal termination of the routine
      !                              the estimates for result and error are
      !                              less reliable. it is assumed that the
      !                              requested accuracy has not been achieved.
      !                       error messages
      !                      ier = 1 maximum number of subdivisions allowed
      !                              has been achieved. one can allow more
      !                              subdivisions by increasing the value of
      !                              limit (and taking the according dimension
      !                              adjustments into account). however, if
      !                              this yield no improvement it is advised
      !                              to analyze the integrand in order to
      !                              determine the integration difficulaties.
      !                              if the position of a local difficulty can
      !                              be determined (i.e.singularity,
      !                              discontinuity within the interval) one
      !                              will probably gain from splitting up the
      !                              interval at this point and calling the
      !                              integrator on the subranges. if possible,
      !                              an appropriate special-purpose integrator
      !                              should be used which is designed for
      !                              handling the type of difficulty involved.
      !                          = 2 the occurrence of roundoff error is
      !                              detected, which prevents the requested
      !                              tolerance from being achieved.
      !                          = 3 extremely bad integrand behaviour occurs
      !                              at some points of the integration
      !                              interval.
      !                          = 6 the input is invalid, because
      !                              (epsabs.le.0 and
      !                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
      !                              or limit.lt.1 or lenw.lt.limit*4.
      !                              result, abserr, neval, last are set
      !                              to zero.
      !                              except when lenw is invalid, iwork(1),
      !                              work(limit*2+1) and work(limit*3+1) are
      !                              set to zero, work(1) is set to a and
      !                              work(limit+1) to b.
      !
      !          dimensioning parameters
      !             limit - integer
      !                     dimensioning parameter for iwork
      !                     limit determines the maximum number of subintervals
      !                     in the partition of the given integration interval
      !                     (a,b), limit.ge.1.
      !                     if limit.lt.1, the routine will end with ier = 6.
      !
      !             lenw  - integer
      !                     dimensioning parameter for work
      !                     lenw must be at least limit*4.
      !                     if lenw.lt.limit*4, the routine will end with
      !                     ier = 6.
      !
      !             last  - integer
      !                     on return, last equals the number of subintervals
      !                     produced in the subdiviosion process, which
      !                     determines the number of significant elements
      !                     actually in the work arrays.
      !
      !          work arrays
      !             iwork - integer
      !                     vector of dimension at least limit, the first k
      !                     elements of which contain pointers to the error
      !                     estimates over the subintervals, such that
      !                     work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
      !                     form a decreasing sequence with k = last if
      !                     last.le.(limit/2+2), and k = limit+1-last otherwise
      !
      !             work  - double precision
      !                     vector of dimension at least lenw
      !                     on return
      !                     work(1), ..., work(last) contain the left end
      !                     points of the subintervals in the partition of
      !                      (a,b),
      !                     work(limit+1), ..., work(limit+last) contain the
      !                      right end points,
      !                     work(limit*2+1), ..., work(limit*2+last) contain
      !                      the integral approximations over the subintervals,
      !                     work(limit*3+1), ..., work(limit*3+last) contain
      !                      the error estimates.
      !
      ! ***references  (none)
      ! ***routines called  dqage,xerror
      ! ***end prologue  dqag
      real*8 a,abserr,b,epsabs,epsrel,f,result,work,d1mach(4),&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup
      integer ier,iwork,key,last,lenw,limit,lvl,l1,l2,l3,neval
      !
      dimension iwork(limit),work(lenw)
      !
      external f
      !
      !          check validity of lenw.
      !
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      ! ***first executable statement  dqag
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
      !
      !          prepare call for dqage.
      !
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      !
      call dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,&
        ier,work(1),work(l1),work(l2),work(l3),iwork,last,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      !
      !          call error handler if necessary.
      !
      lvl = 0
10    if(ier.eq.6) lvl = 1
      !      if(ier.ne.0) call xerror(26habnormal return from dqag ,26,ier,lvl)
      return
      end
      subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,&
         neval,ier,alist,blist,rlist,elist,iord,last,phi,lambda1,zk0,&
           Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqage
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a1
      ! ***keywords  automatic integrator, general-purpose,
      !              integrand examinator, globally adaptive,
      !              gauss-kronrod
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  the routine calculates an approximation result to a given
      !             definite integral   i = integral of f over (a,b),
      !             hopefully satisfying following claim for accuracy
      !             abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
      ! ***description
      !
      !         computation of a definite integral
      !         standard fortran subroutine
      !         double precision version
      !
      !         parameters
      !          on entry
      !             f      - double precision
      !                      function subprogram defining the integrand
      !                      function f(x). the actual name for f needs to be
      !                      declared e x t e r n a l in the driver program.
      !
      !             a      - double precision
      !                      lower limit of integration
      !
      !             b      - double precision
      !                      upper limit of integration
      !
      !             epsabs - double precision
      !                      absolute accuracy requested
      !             epsrel - double precision
      !                      relative accuracy requested
      !                      if  epsabs.le.0
      !                      and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
      !                      the routine will end with ier = 6.
      !
      !             key    - integer
      !                      key for choice of local integration rule
      !                      a gauss-kronrod pair is used with
      !                           7 - 15 points if key.lt.2,
      !                          10 - 21 points if key = 2,
      !                          15 - 31 points if key = 3,
      !                          20 - 41 points if key = 4,
      !                          25 - 51 points if key = 5,
      !                          30 - 61 points if key.gt.5.
      !
      !             limit  - integer
      !                      gives an upperbound on the number of subintervals
      !                      in the partition of (a,b), limit.ge.1.
      !
      !          on return
      !             result - double precision
      !                      approximation to the integral
      !
      !             abserr - double precision
      !                      estimate of the modulus of the absolute error,
      !                      which should equal or exceed abs(i-result)
      !
      !             neval  - integer
      !                      number of integrand evaluations
      !
      !             ier    - integer
      !                      ier = 0 normal and reliable termination of the
      !                              routine. it is assumed that the requested
      !                              accuracy has been achieved.
      !                      ier.gt.0 abnormal termination of the routine
      !                              the estimates for result and error are
      !                              less reliable. it is assumed that the
      !                              requested accuracy has not been achieved.
      !             error messages
      !                      ier = 1 maximum number of subdivisions allowed
      !                              has been achieved. one can allow more
      !                              subdivisions by increasing the value
      !                              of limit.
      !                              however, if this yields no improvement it
      !                              is rather advised to analyze the integrand
      !                              in order to determine the integration
      !                              difficulties. if the position of a local
      !                              difficulty can be determined(e.g.
      !                              singularity, discontinuity within the
      !                              interval) one will probably gain from
      !                              splitting up the interval at this point
      !                              and calling the integrator on the
      !                              subranges. if possible, an appropriate
      !                              special-purpose integrator should be used
      !                              which is designed for handling the type of
      !                              difficulty involved.
      !                          = 2 the occurrence of roundoff error is
      !                              detected, which prevents the requested
      !                              tolerance from being achieved.
      !                          = 3 extremely bad integrand behaviour occurs
      !                              at some points of the integration
      !                              interval.
      !                          = 6 the input is invalid, because
      !                              (epsabs.le.0 and
      !                               epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
      !                              result, abserr, neval, last, rlist(1) ,
      !                              elist(1) and iord(1) are set to zero.
      !                              alist(1) and blist(1) are set to a and b
      !                              respectively.
      !
      !             alist   - double precision
      !                       vector of dimension at least limit, the first
      !                        last  elements of which are the left
      !                       end points of the subintervals in the partition
      !                       of the given integration range (a,b)
      !
      !             blist   - double precision
      !                       vector of dimension at least limit, the first
      !                        last  elements of which are the right
      !                       end points of the subintervals in the partition
      !                       of the given integration range (a,b)
      !
      !             rlist   - double precision
      !                       vector of dimension at least limit, the first
      !                        last  elements of which are the
      !                       integral approximations on the subintervals
      !
      !             elist   - double precision
      !                       vector of dimension at least limit, the first
      !                        last  elements of which are the moduli of the
      !                       absolute error estimates on the subintervals
      !
      !             iord    - integer
      !                       vector of dimension at least limit, the first k
      !                       elements of which are pointers to the
      !                       error estimates over the subintervals,
      !                       such that elist(iord(1)), ...,
      !                       elist(iord(k)) form a decreasing sequence,
      !                       with k = last if last.le.(limit/2+2), and
      !                       k = limit+1-last otherwise
      !
      !             last    - integer
      !                       number of subintervals actually produced in the
      !                       subdivision process
      !
      ! ***references  (none)
      ! ***routines called  d1mach,dqk15,dqk21,dqk31,
      !                     dqk41,dqk51,dqk61,dqpsrt
      ! ***end prologue  dqage
      !
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,&
        blist,b1,b2,dabs,defabs,defab1,defab2,d1mach(4),elist,&
       epmach,epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,&
        resabs,result,rlist,uflow,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,&
        nrmax
      !
      dimension alist(limit),blist(limit),elist(limit),iord(limit),&
        rlist(limit)
      !
      external f


      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      !             list of major variables
      !             -----------------------
      !
      !            alist     - list of left end points of all subintervals
      !                        considered up to now
      !            blist     - list of right end points of all subintervals
      !                        considered up to now
      !            rlist(i)  - approximation to the integral over
      !                       (alist(i),blist(i))
      !            elist(i)  - error estimate applying to rlist(i)
      !            maxerr    - pointer to the interval with largest
      !                        error estimate
      !            errmax    - elist(maxerr)
      !            area      - sum of the integrals over the subintervals
      !            errsum    - sum of the errors over the subintervals
      !            errbnd    - requested accuracy max(epsabs,epsrel*
      !                        abs(result))
      !            *****1    - variable for the left subinterval
      !            *****2    - variable for the right subinterval
      !            last      - index for subdivision
      !
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach  is the largest relative spacing.
      !            uflow  is the smallest positive magnitude.
      !
      ! ***first executable statement  dqage
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      !            test on validity of parameters
      !            ------------------------------
      !
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.&
        epsrel.lt.max(0.5d+02*epmach,0.5d-28)) ier = 6
      if(ier.eq.6) go to 999
      !
      !            first approximation to the integral
      !            -----------------------------------
      !
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs,&
          phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      !
      !            test on accuracy.
      !
      errbnd = max(epsabs,epsrel*dabs(result))
      if(abserr.le.0.5d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)&
        .or.abserr.eq.0.0d+00) go to 60
      !
      !            initialization
      !            --------------
      !
      !
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
      !
      !            main do-loop
      !            ------------
      !
      do 30 last = 2,limit
        !
        !            bisect the subinterval with the largest error estimate.
        !
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2,&
            phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup )
        !
        !            improve previous approximations to integral
        !            and error and test for accuracy.
        !
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(dabs(rlist(maxerr)-area12).le.0.1d-04*dabs(area12)&
        .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 8
        !
        !            test for roundoff error and eventually set error flag.
        !
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
        !
        !            set error flag in the case that the number of subintervals
        !            equals limit.
        !
        if(last.eq.limit) ier = 1
        !
        !            set error flag in the case of bad integrand behaviour
        !            at a point of the integration range.
        !
        if(max(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*&
        epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
    !
    !            append the newly-created intervals to the list.
    !
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
   !
   !            call subroutine dqpsrt to maintain the descending ordering
   !            in the list of error estimates and select the subinterval
   !            with the largest error estimate (to be bisected next).
   !
   20   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax,phi,&
             lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        !  ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
   !
   !            compute final result.
   !            ---------------------
   !
   40 result = 0.0d+00
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
      subroutine dqk15(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk15
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  15-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b), with error
      !                            estimate
      !                        j = integral of abs(f) over (a,b)
      ! ***description
      !
      !            integration rules
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters
      !             on entry
      !               f      - double precision
      !                        function subprogram defining the integrand
      !                        function f(x). the actual name for f needs to be
      !                        declared e x t e r n a l in the calling program.
      !
      !               a      - double precision
      !                        lower limit of integration
      !
      !               b      - double precision
      !                        upper limit of integration
      !
      !             on return
      !               result - double precision
      !                        approximation to the integral i
      !                        result is computed by applying the 15-point
      !                        kronrod rule (resk) obtained by optimal addition
      !                        of abscissae to the7-point gauss rule(resg).
      !
      !               abserr - double precision
      !                        estimate of the modulus of the absolute error,
      !                        which should not exceed abs(i-result)
      !
      !               resabs - double precision
      !                        approximation to the integral j
      !
      !               resasc - double precision
      !                        approximation to the integral of abs(f-i/(b-a))
      !                        over (a,b)
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk15
      !
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1 
      external f
      !
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21


      !
      !            the abscissae and weights are given for the interval (-1,1).
      !            because of symmetry only the positive abscissae and their
      !            corresponding weights are given.
      !
      !            xgk    - abscissae of the 15-point kronrod rule
      !                     xgk(2), xgk(4), ...  abscissae of the 7-point
      !                     gauss rule
      !                     xgk(1), xgk(3), ...  abscissae which are optimally
      !                     added to the 7-point gauss rule
      !
      !            wgk    - weights of the 15-point kronrod rule
      !
      !            wg     - weights of the 7-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.129484966168869693270611432679082d0/
      data wg  (  2) /0.279705391489276667901467771423780d0/
      data wg  (  3) /0.381830050505118944950369775488975d0/
      data wg  (  4) /0.417959183673469387755102040816327d0/
      !
      data xgk (  1) /0.991455371120812639206854697526329d0/
      data xgk (  2) /0.949107912342758524526189684047851d0/
      data xgk (  3) /0.864864423359769072789712788640926d0/
      data xgk (  4) /0.741531185599394439863864773280788d0/
      data xgk (  5) /0.586087235467691130294144838258730d0/
      data xgk (  6) /0.405845151377397166906606412076961d0/
      data xgk (  7) /0.207784955007898467600689403773245d0/
      data xgk (  8) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.022935322010529224963732008058970d0/
      data wgk (  2) /0.063092092629978553290700663189204d0/
      data wgk (  3) /0.104790010322250183839876322541518d0/
      data wgk (  4) /0.140653259715525918745189590510238d0/
      data wgk (  5) /0.169004726639267902826583426598550d0/
      data wgk (  6) /0.190350578064785409913256402421014d0/
      data wgk (  7) /0.204432940075298892414161999234649d0/
      data wgk (  8) /0.209482141084727828012999174891714d0/
      !
      !
      !            list of major variables
      !            -----------------------
      !
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            absc   - abscissa
      !            fval*  - function value
      !            resg   - result of the 7-point gauss formula
      !            resk   - result of the 15-point kronrod formula
      !            reskh  - approximation to the mean value of f over (a,b),
      !                     i.e. to i/(b-a)
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      !
      ! ***first executable statement  dqk15
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 15-point kronrod approximation to
      !            the integral, and estimate the absolute error.
      !
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk21(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk21
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  21-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b), with error
      !                            estimate
      !                        j = integral of abs(f) over (a,b)
      ! ***description
      !
      !            integration rules
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters
      !             on entry
      !               f      - double precision
      !                        function subprogram defining the integrand
      !                        function f(x). the actual name for f needs to be
      !                        declared e x t e r n a l in the driver program.
      !
      !               a      - double precision
      !                        lower limit of integration
      !
      !               b      - double precision
      !                        upper limit of integration
      !
      !             on return
      !               result - double precision
      !                        approximation to the integral i
      !                        result is computed by applying the 21-point
      !                        kronrod rule (resk) obtained by optimal addition
      !                        of abscissae to the 10-point gauss rule (resg).
      !
      !               abserr - double precision
      !                        estimate of the modulus of the absolute error,
      !                        which should not exceed abs(i-result)
      !
      !               resabs - double precision
      !                        approximation to the integral j
      !
      !               resasc - double precision
      !                        approximation to the integral of abs(f-i/(b-a))
      !                        over (a,b)
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk21
      !
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
        zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
      !
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
      
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      !            the abscissae and weights are given for the interval (-1,1).
      !            because of symmetry only the positive abscissae and their
      !            corresponding weights are given.
      !
      !            xgk    - abscissae of the 21-point kronrod rule
      !                     xgk(2), xgk(4), ...  abscissae of the 10-point
      !                     gauss rule
      !                     xgk(1), xgk(3), ...  abscissae which are optimally
      !                     added to the 10-point gauss rule
      !
      !            wgk    - weights of the 21-point kronrod rule
      !
      !            wg     - weights of the 10-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.066671344308688137593568809893332d0/
      data wg  (  2) /0.149451349150580593145776339657697d0/
      data wg  (  3) /0.219086362515982043995534934228163d0/
      data wg  (  4) /0.269266719309996355091226921569469d0/
      data wg  (  5) /0.295524224714752870173892994651338d0/
      !
      data xgk (  1) /0.995657163025808080735527280689003d0/
      data xgk (  2) /0.973906528517171720077964012084452d0/
      data xgk (  3) /0.930157491355708226001207180059508d0/
      data xgk (  4) /0.865063366688984510732096688423493d0/
      data xgk (  5) /0.780817726586416897063717578345042d0/
      data xgk (  6) /0.679409568299024406234327365114874d0/
      data xgk (  7) /0.562757134668604683339000099272694d0/
      data xgk (  8) /0.433395394129247190799265943165784d0/
      data xgk (  9) /0.294392862701460198131126603103866d0/
      data xgk ( 10) /0.148874338981631210884826001129720d0/
      data xgk ( 11) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.011694638867371874278064396062192d0/
      data wgk (  2) /0.032558162307964727478818972459390d0/
      data wgk (  3) /0.054755896574351996031381300244580d0/
      data wgk (  4) /0.075039674810919952767043140916190d0/
      data wgk (  5) /0.093125454583697605535065465083366d0/
      data wgk (  6) /0.109387158802297641899210590325805d0/
      data wgk (  7) /0.123491976262065851077958109831074d0/
      data wgk (  8) /0.134709217311473325928054001771707d0/
      data wgk (  9) /0.142775938577060080797094273138717d0/
      data wgk ( 10) /0.147739104901338491374841515972068d0/
      data wgk ( 11) /0.149445554002916905664936468389821d0/
      !
      !
      !            list of major variables
      !            -----------------------
      !
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            absc   - abscissa
      !            fval*  - function value
      !            resg   - result of the 10-point gauss formula
      !            resk   - result of the 21-point kronrod formula
      !            reskh  - approximation to the mean value of f over (a,b),
      !                     i.e. to i/(b-a)
      !
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      !
      ! ***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 21-point kronrod approximation to
      !            the integral, and estimate the absolute error.
      !
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk31(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk31
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  31-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b) with error
      !                            estimate
      !                        j = integral of abs(f) over (a,b)
      ! ***description
      !
      !            integration rules
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters
      !             on entry
      !               f      - double precision
      !                        function subprogram defining the integrand
      !                        function f(x). the actual name for f needs to be
      !                        declared e x t e r n a l in the calling program.
      !
      !               a      - double precision
      !                        lower limit of integration
      !
      !               b      - double precision
      !                        upper limit of integration
      !
      !             on return
      !               result - double precision
      !                        approximation to the integral i
      !                        result is computed by applying the 31-point
      !                        gauss-kronrod rule (resk), obtained by optimal
      !                        addition of abscissae to the 15-point gauss
      !                        rule (resg).
      !
      !               abserr - double precison
      !                        estimate of the modulus of the modulus,
      !                        which should not exceed abs(i-result)
      !
      !               resabs - double precision
      !                        approximation to the integral j
      !
      !               resasc - double precision
      !                        approximation to the integral of abs(f-i/(b-a))
      !                        over (a,b)
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk31
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
        zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
      !
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      !            the abscissae and weights are given for the interval (-1,1).
      !            because of symmetry only the positive abscissae and their
      !            corresponding weights are given.
      !
      !            xgk    - abscissae of the 31-point kronrod rule
      !                     xgk(2), xgk(4), ...  abscissae of the 15-point
      !                     gauss rule
      !                     xgk(1), xgk(3), ...  abscissae which are optimally
      !                     added to the 15-point gauss rule
      !
      !            wgk    - weights of the 31-point kronrod rule
      !
      !            wg     - weights of the 15-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.030753241996117268354628393577204d0/
      data wg  (  2) /0.070366047488108124709267416450667d0/
      data wg  (  3) /0.107159220467171935011869546685869d0/
      data wg  (  4) /0.139570677926154314447804794511028d0/
      data wg  (  5) /0.166269205816993933553200860481209d0/
      data wg  (  6) /0.186161000015562211026800561866423d0/
      data wg  (  7) /0.198431485327111576456118326443839d0/
      data wg  (  8) /0.202578241925561272880620199967519d0/
      !
      data xgk (  1) /0.998002298693397060285172840152271d0/
      data xgk (  2) /0.987992518020485428489565718586613d0/
      data xgk (  3) /0.967739075679139134257347978784337d0/
      data xgk (  4) /0.937273392400705904307758947710209d0/
      data xgk (  5) /0.897264532344081900882509656454496d0/
      data xgk (  6) /0.848206583410427216200648320774217d0/
      data xgk (  7) /0.790418501442465932967649294817947d0/
      data xgk (  8) /0.724417731360170047416186054613938d0/
      data xgk (  9) /0.650996741297416970533735895313275d0/
      data xgk ( 10) /0.570972172608538847537226737253911d0/
      data xgk ( 11) /0.485081863640239680693655740232351d0/
      data xgk ( 12) /0.394151347077563369897207370981045d0/
      data xgk ( 13) /0.299180007153168812166780024266389d0/
      data xgk ( 14) /0.201194093997434522300628303394596d0/
      data xgk ( 15) /0.101142066918717499027074231447392d0/
      data xgk ( 16) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.005377479872923348987792051430128d0/
      data wgk (  2) /0.015007947329316122538374763075807d0/
      data wgk (  3) /0.025460847326715320186874001019653d0/
      data wgk (  4) /0.035346360791375846222037948478360d0/
      data wgk (  5) /0.044589751324764876608227299373280d0/
      data wgk (  6) /0.053481524690928087265343147239430d0/
      data wgk (  7) /0.062009567800670640285139230960803d0/
      data wgk (  8) /0.069854121318728258709520077099147d0/
      data wgk (  9) /0.076849680757720378894432777482659d0/
      data wgk ( 10) /0.083080502823133021038289247286104d0/
      data wgk ( 11) /0.088564443056211770647275443693774d0/
      data wgk ( 12) /0.093126598170825321225486872747346d0/
      data wgk ( 13) /0.096642726983623678505179907627589d0/
      data wgk ( 14) /0.099173598721791959332393173484603d0/
      data wgk ( 15) /0.100769845523875595044946662617570d0/
      data wgk ( 16) /0.101330007014791549017374792767493d0/
      !
      !
      !            list of major variables
      !            -----------------------
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            absc   - abscissa
      !            fval*  - function value
      !            resg   - result of the 15-point gauss formula
      !            resk   - result of the 31-point kronrod formula
      !            reskh  - approximation to the mean value of f over (a,b),
      !                     i.e. to i/(b-a)
      !
      !            machine dependent constants
      !            ---------------------------
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      ! ***first executable statement  dqk31
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 31-point kronrod approximation to
      !            the integral, and estimate the absolute error.
      !
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(16)*dabs(fc-reskh)
      do 20 j=1,15
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk41(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk41
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  41-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b), with error
      !                            estimate
      !                        j = integral of abs(f) over (a,b)
      ! ***description
      !
      !            integration rules
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters
      !             on entry
      !               f      - double precision
      !                        function subprogram defining the integrand
      !                        function f(x). the actual name for f needs to be
      !                        declared e x t e r n a l in the calling program.
      !
      !               a      - double precision
      !                        lower limit of integration
      !
      !               b      - double precision
      !                        upper limit of integration
      !
      !             on return
      !               result - double precision
      !                        approximation to the integral i
      !                        result is computed by applying the 41-point
      !                        gauss-kronrod rule (resk) obtained by optimal
      !                        addition of abscissae to the 20-point gauss
      !                        rule (resg).
      !
      !               abserr - double precision
      !                        estimate of the modulus of the absolute error,
      !                        which should not exceed abs(i-result)
      !
      !               resabs - double precision
      !                        approximation to the integral j
      !
      !               resasc - double precision
      !                        approximation to the integal of abs(f-i/(b-a))
      !                        over (a,b)
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk41
      !
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
        zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
      !
      dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      !            the abscissae and weights are given for the interval (-1,1).
      !            because of symmetry only the positive abscissae and their
      !            corresponding weights are given.
      !
      !            xgk    - abscissae of the 41-point gauss-kronrod rule
      !                     xgk(2), xgk(4), ...  abscissae of the 20-point
      !                     gauss rule
      !                     xgk(1), xgk(3), ...  abscissae which are optimally
      !                     added to the 20-point gauss rule
      !
      !            wgk    - weights of the 41-point gauss-kronrod rule
      !
      !            wg     - weights of the 20-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.017614007139152118311861962351853d0/
      data wg  (  2) /0.040601429800386941331039952274932d0/
      data wg  (  3) /0.062672048334109063569506535187042d0/
      data wg  (  4) /0.083276741576704748724758143222046d0/
      data wg  (  5) /0.101930119817240435036750135480350d0/
      data wg  (  6) /0.118194531961518417312377377711382d0/
      data wg  (  7) /0.131688638449176626898494499748163d0/
      data wg  (  8) /0.142096109318382051329298325067165d0/
      data wg  (  9) /0.149172986472603746787828737001969d0/
      data wg  ( 10) /0.152753387130725850698084331955098d0/
      !
      data xgk (  1) /0.998859031588277663838315576545863d0/
      data xgk (  2) /0.993128599185094924786122388471320d0/
      data xgk (  3) /0.981507877450250259193342994720217d0/
      data xgk (  4) /0.963971927277913791267666131197277d0/
      data xgk (  5) /0.940822633831754753519982722212443d0/
      data xgk (  6) /0.912234428251325905867752441203298d0/
      data xgk (  7) /0.878276811252281976077442995113078d0/
      data xgk (  8) /0.839116971822218823394529061701521d0/
      data xgk (  9) /0.795041428837551198350638833272788d0/
      data xgk ( 10) /0.746331906460150792614305070355642d0/
      data xgk ( 11) /0.693237656334751384805490711845932d0/
      data xgk ( 12) /0.636053680726515025452836696226286d0/
      data xgk ( 13) /0.575140446819710315342946036586425d0/
      data xgk ( 14) /0.510867001950827098004364050955251d0/
      data xgk ( 15) /0.443593175238725103199992213492640d0/
      data xgk ( 16) /0.373706088715419560672548177024927d0/
      data xgk ( 17) /0.301627868114913004320555356858592d0/
      data xgk ( 18) /0.227785851141645078080496195368575d0/
      data xgk ( 19) /0.152605465240922675505220241022678d0/
      data xgk ( 20) /0.076526521133497333754640409398838d0/
      data xgk ( 21) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.003073583718520531501218293246031d0/
      data wgk (  2) /0.008600269855642942198661787950102d0/
      data wgk (  3) /0.014626169256971252983787960308868d0/
      data wgk (  4) /0.020388373461266523598010231432755d0/
      data wgk (  5) /0.025882133604951158834505067096153d0/
      data wgk (  6) /0.031287306777032798958543119323801d0/
      data wgk (  7) /0.036600169758200798030557240707211d0/
      data wgk (  8) /0.041668873327973686263788305936895d0/
      data wgk (  9) /0.046434821867497674720231880926108d0/
      data wgk ( 10) /0.050944573923728691932707670050345d0/
      data wgk ( 11) /0.055195105348285994744832372419777d0/
      data wgk ( 12) /0.059111400880639572374967220648594d0/
      data wgk ( 13) /0.062653237554781168025870122174255d0/
      data wgk ( 14) /0.065834597133618422111563556969398d0/
      data wgk ( 15) /0.068648672928521619345623411885368d0/
      data wgk ( 16) /0.071054423553444068305790361723210d0/
      data wgk ( 17) /0.073030690332786667495189417658913d0/
      data wgk ( 18) /0.074582875400499188986581418362488d0/
      data wgk ( 19) /0.075704497684556674659542775376617d0/
      data wgk ( 20) /0.076377867672080736705502835038061d0/
      data wgk ( 21) /0.076600711917999656445049901530102d0/
      !
      !
      !            list of major variables
      !            -----------------------
      !
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            absc   - abscissa
      !            fval*  - function value
      !            resg   - result of the 20-point gauss formula
      !            resk   - result of the 41-point kronrod formula
      !            reskh  - approximation to mean value of f over (a,b), i.e.
      !                     to i/(b-a)
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      !
      ! ***first executable statement  dqk41
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 41-point gauss-kronrod approximation to
      !            the integral, and estimate the absolute error.
      !
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(21)*fc
      resabs = dabs(resk)
      do 10 j=1,10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(21)*dabs(fc-reskh)
      do 20 j=1,20
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk51(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk51
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  51-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math & progr. div. - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b) with error
      !                            estimate
      !                        j = integral of abs(f) over (a,b)
      ! ***description
      !
      !            integration rules
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters
      !             on entry
      !               f      - double precision
      !                        function subroutine defining the integrand
      !                        function f(x). the actual name for f needs to be
      !                        declared e x t e r n a l in the calling program.
      !
      !               a      - double precision
      !                        lower limit of integration
      !
      !               b      - double precision
      !                        upper limit of integration
      !
      !             on return
      !               result - double precision
      !                        approximation to the integral i
      !                        result is computed by applying the 51-point
      !                        kronrod rule (resk) obtained by optimal addition
      !                        of abscissae to the 25-point gauss rule (resg).
      !
      !               abserr - double precision
      !                        estimate of the modulus of the absolute error,
      !                        which should not exceed abs(i-result)
      !
      !               resabs - double precision
      !                        approximation to the integral j
      !
      !               resasc - double precision
      !                        approximation to the integral of abs(f-i/(b-a))
      !                        over (a,b)
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk51
      !
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
           resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
      !
      dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)
      d1mach(1)=1E21
      d1mach(2)=0d0
      d1mach(3)=0d0
      d1mach(4)=1E-21
      !
      !            the abscissae and weights are given for the interval (-1,1).
      !            because of symmetry only the positive abscissae and their
      !            corresponding weights are given.
      !
      !            xgk    - abscissae of the 51-point kronrod rule
      !                     xgk(2), xgk(4), ...  abscissae of the 25-point
      !                     gauss rule
      !                     xgk(1), xgk(3), ...  abscissae which are optimally
      !                     added to the 25-point gauss rule
      !
      !            wgk    - weights of the 51-point kronrod rule
      !
      !            wg     - weights of the 25-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.011393798501026287947902964113235d0/
      data wg  (  2) /0.026354986615032137261901815295299d0/
      data wg  (  3) /0.040939156701306312655623487711646d0/
      data wg  (  4) /0.054904695975835191925936891540473d0/
      data wg  (  5) /0.068038333812356917207187185656708d0/
      data wg  (  6) /0.080140700335001018013234959669111d0/
      data wg  (  7) /0.091028261982963649811497220702892d0/
      data wg  (  8) /0.100535949067050644202206890392686d0/
      data wg  (  9) /0.108519624474263653116093957050117d0/
      data wg  ( 10) /0.114858259145711648339325545869556d0/
      data wg  ( 11) /0.119455763535784772228178126512901d0/
      data wg  ( 12) /0.122242442990310041688959518945852d0/
      data wg  ( 13) /0.123176053726715451203902873079050d0/
      !
      data xgk (  1) /0.999262104992609834193457486540341d0/
      data xgk (  2) /0.995556969790498097908784946893902d0/
      data xgk (  3) /0.988035794534077247637331014577406d0/
      data xgk (  4) /0.976663921459517511498315386479594d0/
      data xgk (  5) /0.961614986425842512418130033660167d0/
      data xgk (  6) /0.942974571228974339414011169658471d0/
      data xgk (  7) /0.920747115281701561746346084546331d0/
      data xgk (  8) /0.894991997878275368851042006782805d0/
      data xgk (  9) /0.865847065293275595448996969588340d0/
      data xgk ( 10) /0.833442628760834001421021108693570d0/
      data xgk ( 11) /0.797873797998500059410410904994307d0/
      data xgk ( 12) /0.759259263037357630577282865204361d0/
      data xgk ( 13) /0.717766406813084388186654079773298d0/
      data xgk ( 14) /0.673566368473468364485120633247622d0/
      data xgk ( 15) /0.626810099010317412788122681624518d0/
      data xgk ( 16) /0.577662930241222967723689841612654d0/
      data xgk ( 17) /0.526325284334719182599623778158010d0/
      data xgk ( 18) /0.473002731445714960522182115009192d0/
      data xgk ( 19) /0.417885382193037748851814394594572d0/
      data xgk ( 20) /0.361172305809387837735821730127641d0/
      data xgk ( 21) /0.303089538931107830167478909980339d0/
      data xgk ( 22) /0.243866883720988432045190362797452d0/
      data xgk ( 23) /0.183718939421048892015969888759528d0/
      data xgk ( 24) /0.122864692610710396387359818808037d0/
      data xgk ( 25) /0.061544483005685078886546392366797d0/
      data xgk ( 26) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.001987383892330315926507851882843d0/
      data wgk (  2) /0.005561932135356713758040236901066d0/
      data wgk (  3) /0.009473973386174151607207710523655d0/
      data wgk (  4) /0.013236229195571674813656405846976d0/
      data wgk (  5) /0.016847817709128298231516667536336d0/
      data wgk (  6) /0.020435371145882835456568292235939d0/
      data wgk (  7) /0.024009945606953216220092489164881d0/
      data wgk (  8) /0.027475317587851737802948455517811d0/
      data wgk (  9) /0.030792300167387488891109020215229d0/
      data wgk ( 10) /0.034002130274329337836748795229551d0/
      data wgk ( 11) /0.037116271483415543560330625367620d0/
      data wgk ( 12) /0.040083825504032382074839284467076d0/
      data wgk ( 13) /0.042872845020170049476895792439495d0/
      data wgk ( 14) /0.045502913049921788909870584752660d0/
      data wgk ( 15) /0.047982537138836713906392255756915d0/
      data wgk ( 16) /0.050277679080715671963325259433440d0/
      data wgk ( 17) /0.052362885806407475864366712137873d0/
      data wgk ( 18) /0.054251129888545490144543370459876d0/
      data wgk ( 19) /0.055950811220412317308240686382747d0/
      data wgk ( 20) /0.057437116361567832853582693939506d0/
      data wgk ( 21) /0.058689680022394207961974175856788d0/
      data wgk ( 22) /0.059720340324174059979099291932562d0/
      data wgk ( 23) /0.060539455376045862945360267517565d0/
      data wgk ( 24) /0.061128509717053048305859030416293d0/
      data wgk ( 25) /0.061471189871425316661544131965264d0/
      !        note: wgk (26) was calculated from the values of wgk(1..25)
      data wgk ( 26) /0.061580818067832935078759824240066d0/
      !
      !
      !            list of major variables
      !            -----------------------
      !
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            absc   - abscissa
      !            fval*  - function value
      !            resg   - result of the 25-point gauss formula
      !            resk   - result of the 51-point kronrod formula
      !            reskh  - approximation to the mean value of f over (a,b),
      !                     i.e. to i/(b-a)
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      !
      ! ***first executable statement  dqk51
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 51-point kronrod approximation to
      !            the integral, and estimate the absolute error.
      !
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = dabs(resk)
      do 10 j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+absc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(26)*dabs(fc-reskh)
      do 20 j=1,25
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk61(f,a,b,result,abserr,resabs,resasc,phi,lambda1,&
           zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqk61
      ! ***date written   800101   (yymmdd)
      ! ***revision date  830518   (yymmdd)
      ! ***category no.  h2a1a2
      ! ***keywords  61-point gauss-kronrod rules
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  to compute i = integral of f over (a,b) with error
      !                            estimate
      !                        j = integral of dabs(f) over (a,b)
      ! ***description
      !
      !         integration rule
      !         standard fortran subroutine
      !         double precision version
      !
      !
      !         parameters
      !          on entry
      !            f      - double precision
      !                     function subprogram defining the integrand
      !                     function f(x). the actual name for f needs to be
      !                     declared e x t e r n a l in the calling program.
      !
      !            a      - double precision
      !                     lower limit of integration
      !
      !            b      - double precision
      !                     upper limit of integration
      !
      !          on return
      !            result - double precision
      !                     approximation to the integral i
      !                     result is computed by applying the 61-point
      !                     kronrod rule (resk) obtained by optimal addition of
      !                     abscissae to the 30-point gauss rule (resg).
      !
      !            abserr - double precision
      !                     estimate of the modulus of the absolute error,
      !                     which should equal or exceed dabs(i-result)
      !
      !            resabs - double precision
      !                     approximation to the integral j
      !
      !            resasc - double precision
      !                     approximation to the integral of dabs(f-i/(b-a))
      !
      !
      ! ***references  (none)
      ! ***routines called  d1mach
      ! ***end prologue  dqk61
      !
      double precision a,dabsc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,&
        d1mach(4),epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk,phi,lambda1,&
        zk0,Pup,Tup,rurd,xflow,kup
      integer j,jtw,jtwm1
      external f
      !
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
      d1mach(1)=1E21
      d1mach(2)=0
      d1mach(3)=0
      d1mach(4)=1E-21
      !
      !            the abscissae and weights are given for the
      !            interval (-1,1). because of symmetry only the positive
      !            abscissae and their corresponding weights are given.
      !
      !            xgk   - abscissae of the 61-point kronrod rule
      !                    xgk(2), xgk(4)  ... abscissae of the 30-point
      !                    gauss rule
      !                    xgk(1), xgk(3)  ... optimally added abscissae
      !                    to the 30-point gauss rule
      !
      !            wgk   - weights of the 61-point kronrod rule
      !
      !            wg    - weigths of the 30-point gauss rule
      !
      !
      !  gauss quadrature weights and kronron quadrature abscissae and weights
      !  as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
      !  bell labs, nov. 1981.
      !
      data wg  (  1) /0.007968192496166605615465883474674d0/
      data wg  (  2) /0.018466468311090959142302131912047d0/
      data wg  (  3) /0.028784707883323369349719179611292d0/
      data wg  (  4) /0.038799192569627049596801936446348d0/
      data wg  (  5) /0.048402672830594052902938140422808d0/
      data wg  (  6) /0.057493156217619066481721689402056d0/
      data wg  (  7) /0.065974229882180495128128515115962d0/
      data wg  (  8) /0.073755974737705206268243850022191d0/
      data wg  (  9) /0.080755895229420215354694938460530d0/
      data wg  ( 10) /0.086899787201082979802387530715126d0/
      data wg  ( 11) /0.092122522237786128717632707087619d0/
      data wg  ( 12) /0.096368737174644259639468626351810d0/
      data wg  ( 13) /0.099593420586795267062780282103569d0/
      data wg  ( 14) /0.101762389748405504596428952168554d0/
      data wg  ( 15) /0.102852652893558840341285636705415d0/
      !
      data xgk (  1) /0.999484410050490637571325895705811d0/
      data xgk (  2) /0.996893484074649540271630050918695d0/
      data xgk (  3) /0.991630996870404594858628366109486d0/
      data xgk (  4) /0.983668123279747209970032581605663d0/
      data xgk (  5) /0.973116322501126268374693868423707d0/
      data xgk (  6) /0.960021864968307512216871025581798d0/
      data xgk (  7) /0.944374444748559979415831324037439d0/
      data xgk (  8) /0.926200047429274325879324277080474d0/
      data xgk (  9) /0.905573307699907798546522558925958d0/
      data xgk ( 10) /0.882560535792052681543116462530226d0/
      data xgk ( 11) /0.857205233546061098958658510658944d0/
      data xgk ( 12) /0.829565762382768397442898119732502d0/
      data xgk ( 13) /0.799727835821839083013668942322683d0/
      data xgk ( 14) /0.767777432104826194917977340974503d0/
      data xgk ( 15) /0.733790062453226804726171131369528d0/
      data xgk ( 16) /0.697850494793315796932292388026640d0/
      data xgk ( 17) /0.660061064126626961370053668149271d0/
      data xgk ( 18) /0.620526182989242861140477556431189d0/
      data xgk ( 19) /0.579345235826361691756024932172540d0/
      data xgk ( 20) /0.536624148142019899264169793311073d0/
      data xgk ( 21) /0.492480467861778574993693061207709d0/
      data xgk ( 22) /0.447033769538089176780609900322854d0/
      data xgk ( 23) /0.400401254830394392535476211542661d0/
      data xgk ( 24) /0.352704725530878113471037207089374d0/
      data xgk ( 25) /0.304073202273625077372677107199257d0/
      data xgk ( 26) /0.254636926167889846439805129817805d0/
      data xgk ( 27) /0.204525116682309891438957671002025d0/
      data xgk ( 28) /0.153869913608583546963794672743256d0/
      data xgk ( 29) /0.102806937966737030147096751318001d0/
      data xgk ( 30) /0.051471842555317695833025213166723d0/
      data xgk ( 31) /0.000000000000000000000000000000000d0/
      !
      data wgk (  1) /0.001389013698677007624551591226760d0/
      data wgk (  2) /0.003890461127099884051267201844516d0/
      data wgk (  3) /0.006630703915931292173319826369750d0/
      data wgk (  4) /0.009273279659517763428441146892024d0/
      data wgk (  5) /0.011823015253496341742232898853251d0/
      data wgk (  6) /0.014369729507045804812451432443580d0/
      data wgk (  7) /0.016920889189053272627572289420322d0/
      data wgk (  8) /0.019414141193942381173408951050128d0/
      data wgk (  9) /0.021828035821609192297167485738339d0/
      data wgk ( 10) /0.024191162078080601365686370725232d0/
      data wgk ( 11) /0.026509954882333101610601709335075d0/
      data wgk ( 12) /0.028754048765041292843978785354334d0/
      data wgk ( 13) /0.030907257562387762472884252943092d0/
      data wgk ( 14) /0.032981447057483726031814191016854d0/
      data wgk ( 15) /0.034979338028060024137499670731468d0/
      data wgk ( 16) /0.036882364651821229223911065617136d0/
      data wgk ( 17) /0.038678945624727592950348651532281d0/
      data wgk ( 18) /0.040374538951535959111995279752468d0/
      data wgk ( 19) /0.041969810215164246147147541285970d0/
      data wgk ( 20) /0.043452539701356069316831728117073d0/
      data wgk ( 21) /0.044814800133162663192355551616723d0/
      data wgk ( 22) /0.046059238271006988116271735559374d0/
      data wgk ( 23) /0.047185546569299153945261478181099d0/
      data wgk ( 24) /0.048185861757087129140779492298305d0/
      data wgk ( 25) /0.049055434555029778887528165367238d0/
      data wgk ( 26) /0.049795683427074206357811569379942d0/
      data wgk ( 27) /0.050405921402782346840893085653585d0/
      data wgk ( 28) /0.050881795898749606492297473049805d0/
      data wgk ( 29) /0.051221547849258772170656282604944d0/
      data wgk ( 30) /0.051426128537459025933862879215781d0/
      data wgk ( 31) /0.051494729429451567558340433647099d0/
      !
      !            list of major variables
      !            -----------------------
      !
      !            centr  - mid point of the interval
      !            hlgth  - half-length of the interval
      !            dabsc  - abscissa
      !            fval*  - function value
      !            resg   - result of the 30-point gauss rule
      !            resk   - result of the 61-point kronrod rule
      !            reskh  - approximation to the mean value of f
      !                     over (a,b), i.e. to i/(b-a)
      !
      !            machine dependent constants
      !            ---------------------------
      !
      !            epmach is the largest relative spacing.
      !            uflow is the smallest positive magnitude.
      !
      epmach = d1mach(4)
      uflow = d1mach(1)
      !
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      !
      !            compute the 61-point kronrod approximation to the
      !            integral, and estimate the absolute error.
      !
      ! ***first executable statement  dqk61
      resg = 0.0d+00
      fc = f(centr,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      resk = wgk(31)*fc
      resabs = dabs(resk)
      do 10 j=1,15
        jtw = j*2
        dabsc = hlgth*xgk(jtw)
        fval1 = f(centr-dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = f(centr-dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fval2 = f(centr+dabsc,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
  15    continue
      reskh = resk*0.5d+00
      resasc = wgk(31)*dabs(fc-reskh)
      do 20 j=1,30
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)&
        abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
        ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax,&
           phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
      ! ***begin prologue  dqpsrt
      ! ***refer to  dqage,dqagie,dqagpe,dqawse
      ! ***routines called  (none)
      ! ***revision date  810101   (yymmdd)
      ! ***keywords  sequential sorting
      ! ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
      !            de doncker,elise,appl. math. & progr. div. - k.u.leuven
      ! ***purpose  this routine maintains the descending ordering in the
      !             list of the local error estimated resulting from the
      !             interval subdivision process. at each call two error
      !             estimates are inserted using the sequential search
      !             method, top-down for the largest error estimate and
      !             bottom-up for the smallest error estimate.
      ! ***description
      !
      !            ordering routine
      !            standard fortran subroutine
      !            double precision version
      !
      !            parameters (meaning at output)
      !               limit  - integer
      !                        maximum number of error estimates the list
      !                        can contain
      !
      !               last   - integer
      !                        number of error estimates currently in the list
      !
      !               maxerr - integer
      !                        maxerr points to the nrmax-th largest error
      !                        estimate currently in the list
      !
      !               ermax  - double precision
      !                        nrmax-th largest error estimate
      !                        ermax = elist(maxerr)
      !
      !               elist  - double precision
      !                        vector of dimension last containing
      !                        the error estimates
      !
      !               iord   - integer
      !                        vector of dimension last, the first k elements
      !                        of which contain pointers to the error
      !                        estimates, such that
      !                        elist(iord(1)),...,  elist(iord(k))
      !                        form a decreasing sequence, with
      !                        k = last if last.le.(limit/2+2), and
      !                        k = limit+1-last otherwise
      !
      !               nrmax  - integer
      !                        maxerr = iord(nrmax)
      !
      ! ***end prologue  dqpsrt
      !
      double precision elist,ermax,errmax,errmin,phi,lambda1,zk0,&
           Pup,Tup,rurd,xflow,kup
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,&
        nrmax
      dimension elist(last),iord(last)
      !
      !            check whether the list contains more than
      !            two error estimates.
      !
      ! ***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
   !
   !            this part of the routine is only executed if, due to a
   !            difficult integrand, subdivision increased the error
   !            estimate. in the normal case the insert procedure should
   !            start after the nrmax-th largest error estimate.
   !
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
        !  ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
   !
   !            compute the number of elements in the list to be maintained
   !            in descending order. this number depends on the number of
   !            subdivisions still allowed.
   !
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
      !
      !            insert errmax by traversing the list top-down,
      !            starting comparison from the element elist(iord(nrmax+1)).
      !
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
        !  ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
   !
   !            insert errmin by traversing the list bottom-up.
   !
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
        !  ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
   !
   !            set maxerr and ermax.
   !
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
