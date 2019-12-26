      !
      !    SLATEC: public domain
      !
      ! DECK DDEABM
      SUBROUTINE DDEABM (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,&
         RWORK, LRW, IWORK, LIW, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DDEABM
      ! ***PURPOSE  Solve an initial value problem in ordinary differential
      !             equations using an Adams-Bashforth method.
      ! ***LIBRARY   SLATEC (DEPAC)
      ! ***CATEGORY  I1A1B
      ! ***TYPE      DOUBLE PRECISION (DEABM-S, DDEABM-D)
      ! ***KEYWORDS  ADAMS-BASHFORTH METHOD, DEPAC, INITIAL VALUE PROBLEMS,
      !              ODE, ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
      ! ***AUTHOR  Shampine, L. F., (SNLA)
      !            Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    This is the Adams code in the package of differential equation
      !    solvers DEPAC, consisting of the codes DDERKF, DDEABM, and DDEBDF.
      !    Design of the package was by L. F. Shampine and H. A. Watts.
      !    It is documented in
      !         SAND79-2374 , DEPAC - Design of a User Oriented Package of ODE
      !                               Solvers.
      !    DDEABM is a driver for a modification of the code ODE written by
      !              L. F. Shampine and M. K. Gordon
      !              Sandia Laboratories
      !              Albuquerque, New Mexico 87185
      !
      !  **********************************************************************
      !  * ABSTRACT *
      !  ************
      !
      !    Subroutine DDEABM uses the Adams-Bashforth-Moulton
      !    Predictor-Corrector formulas of orders one through twelve to
      !    integrate a system of NEQ first order ordinary differential
      !    equations of the form
      !                          DU/DX = DF(X,U)
      !    when the vector Y(*) of initial values for U(*) at X=T is given.
      !    The subroutine integrates from T to TOUT. It is easy to continue the
      !    integration to get results at additional TOUT.  This is the interval
      !    mode of operation.  It is also easy for the routine to return with
      !    the solution at each intermediate step on the way to TOUT.  This is
      !    the intermediate-output mode of operation.
      !
      !    DDEABM uses subprograms DDES, DSTEPS, DINTP, DHSTRT, DHVNRM,
      !    D1MACH, and the error handling routine XERMSG.  The only machine
      !    dependent parameters to be assigned appear in D1MACH.
      !
      !  **********************************************************************
      !  * Description of The Arguments To DDEABM (An Overview) *
      !  **********************************************************************
      !
      !    The Parameters are
      !
      !       DF -- This is the name of a subroutine which you provide to
      !              define the differential equations.
      !
      !       NEQ -- This is the number of (first order) differential
      !              equations to be integrated.
      !
      !       T -- This is a DOUBLE PRECISION value of the independent
      !            variable.
      !
      !       Y(*) -- This DOUBLE PRECISION array contains the solution
      !               components at T.
      !
      !       TOUT -- This is a DOUBLE PRECISION point at which a solution is
      !               desired.
      !
      !       INFO(*) -- The basic task of the code is to integrate the
      !              differential equations from T to TOUT and return an
      !              answer at TOUT.  INFO(*) is an INTEGER array which is used
      !              to communicate exactly how you want this task to be
      !              carried out.
      !
      !       RTOL, ATOL -- These DOUBLE PRECISION quantities represent
      !                     relative and absolute error tolerances which you
      !                     provide to indicate how accurately you wish the
      !                     solution to be computed.  You may choose them to be
      !                     both scalars or else both vectors.
      !
      !       IDID -- This scalar quantity is an indicator reporting what
      !              the code did.  You must monitor this INTEGER variable to
      !              decide what action to take next.
      !
      !       RWORK(*), LRW -- RWORK(*) is a DOUBLE PRECISION work array of
      !              length LRW which provides the code with needed storage
      !              space.
      !
      !       IWORK(*), LIW -- IWORK(*) is an INTEGER work array of length LIW
      !              which provides the code with needed storage space and an
      !              across call flag.
      !
      !       RPAR, IPAR -- These are DOUBLE PRECISION and INTEGER parameter
      !              arrays which you can use for communication between your
      !              calling program and the DF subroutine.
      !
      !   Quantities which are used as input items are
      !              NEQ, T, Y(*), TOUT, INFO(*),
      !              RTOL, ATOL, RWORK(1), LRW and LIW.
      !
      !   Quantities which may be altered by the code are
      !              T, Y(*), INFO(1), RTOL, ATOL,
      !              IDID, RWORK(*) and IWORK(*).
      !
      !  **********************************************************************
      !  * INPUT -- What To Do On The First Call To DDEABM *
      !  **********************************************************************
      !
      !    The first call of the code is defined to be the start of each new
      !    problem.  Read through the descriptions of all the following items,
      !    provide sufficient storage space for designated arrays, set
      !    appropriate variables for the initialization of the problem, and
      !    give information about how you want the problem to be solved.
      !
      !
      !       DF -- Provide a subroutine of the form
      !                                DF(X,U,UPRIME,RPAR,IPAR)
      !              to define the system of first order differential equations
      !              which is to be solved.  For the given values of X and the
      !              vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
      !              evaluate the NEQ components of the system of differential
      !              equations  DU/DX=DF(X,U)  and store the derivatives in the
      !              array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
      !              equations I=1,...,NEQ.
      !
      !              Subroutine DF must NOT alter X or U(*).  You must declare
      !              the name df in an external statement in your program that
      !              calls DDEABM.  You must dimension U and UPRIME in DF.
      !
      !              RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
      !              arrays which you can use for communication between your
      !              calling program and subroutine DF. They are not used or
      !              altered by DDEABM.  If you do not need RPAR or IPAR,
      !              ignore these parameters by treating them as dummy
      !              arguments. If you do choose to use them, dimension them in
      !              your calling program and in DF as arrays of appropriate
      !              length.
      !
      !       NEQ -- Set it to the number of differential equations.
      !              (NEQ .GE. 1)
      !
      !       T -- Set it to the initial point of the integration.
      !              You must use a program variable for T because the code
      !              changes its value.
      !
      !       Y(*) -- Set this vector to the initial values of the NEQ solution
      !              components at the initial point.  You must dimension Y at
      !              least NEQ in your calling program.
      !
      !       TOUT -- Set it to the first point at which a solution
      !              is desired.  You can take TOUT = T, in which case the code
      !              will evaluate the derivative of the solution at T and
      !              return. Integration either forward in T  (TOUT .GT. T)  or
      !              backward in T  (TOUT .LT. T)  is permitted.
      !
      !              The code advances the solution from T to TOUT using
      !              step sizes which are automatically selected so as to
      !              achieve the desired accuracy.  If you wish, the code will
      !              return with the solution and its derivative following
      !              each intermediate step (intermediate-output mode) so that
      !              you can monitor them, but you still must provide TOUT in
      !              accord with the basic aim of the code.
      !
      !              The first step taken by the code is a critical one
      !              because it must reflect how fast the solution changes near
      !              the initial point.  The code automatically selects an
      !              initial step size which is practically always suitable for
      !              the problem. By using the fact that the code will not step
      !              past TOUT in the first step, you could, if necessary,
      !              restrict the length of the initial step size.
      !
      !              For some problems it may not be permissible to integrate
      !              past a point TSTOP because a discontinuity occurs there
      !              or the solution or its derivative is not defined beyond
      !              TSTOP.  When you have declared a TSTOP point (see INFO(4)
      !              and RWORK(1)), you have told the code not to integrate
      !              past TSTOP.  In this case any TOUT beyond TSTOP is invalid
      !              input.
      !
      !       INFO(*) -- Use the INFO array to give the code more details about
      !              how you want your problem solved.  This array should be
      !              dimensioned of length 15 to accommodate other members of
      !              DEPAC or possible future extensions, though DDEABM uses
      !              only the first four entries.  You must respond to all of
      !              the following items which are arranged as questions.  The
      !              simplest use of the code corresponds to answering all
      !              questions as YES ,i.e. setting ALL entries of INFO to 0.
      !
      !         INFO(1) -- This parameter enables the code to initialize
      !                itself.  You must set it to indicate the start of every
      !                new problem.
      !
      !             **** Is this the first call for this problem ...
      !                   YES -- set INFO(1) = 0
      !                    NO -- not applicable here.
      !                          See below for continuation calls.  ****
      !
      !         INFO(2) -- How much accuracy you want of your solution
      !                is specified by the error tolerances RTOL and ATOL.
      !                The simplest use is to take them both to be scalars.
      !                To obtain more flexibility, they can both be vectors.
      !                The code must be told your choice.
      !
      !             **** Are both error tolerances RTOL, ATOL scalars ...
      !                   YES -- set INFO(2) = 0
      !                          and input scalars for both RTOL and ATOL
      !                    NO -- set INFO(2) = 1
      !                          and input arrays for both RTOL and ATOL ****
      !
      !         INFO(3) -- The code integrates from T in the direction
      !                of TOUT by steps.  If you wish, it will return the
      !                computed solution and derivative at the next
      !                intermediate step (the intermediate-output mode) or
      !                TOUT, whichever comes first.  This is a good way to
      !                proceed if you want to see the behavior of the solution.
      !                If you must have solutions at a great many specific
      !                TOUT points, this code will compute them efficiently.
      !
      !             **** Do you want the solution only at
      !                  TOUT (and not at the next intermediate step) ...
      !                   YES -- set INFO(3) = 0
      !                    NO -- set INFO(3) = 1 ****
      !
      !         INFO(4) -- To handle solutions at a great many specific
      !                values TOUT efficiently, this code may integrate past
      !                TOUT and interpolate to obtain the result at TOUT.
      !                Sometimes it is not possible to integrate beyond some
      !                point TSTOP because the equation changes there or it is
      !                not defined past TSTOP.  Then you must tell the code
      !                not to go past.
      !
      !             **** Can the integration be carried out without any
      !                  Restrictions on the independent variable T ...
      !                   YES -- set INFO(4)=0
      !                    NO -- set INFO(4)=1
      !                          and define the stopping point TSTOP by
      !                          setting RWORK(1)=TSTOP ****
      !
      !       RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
      !              error tolerances to tell the code how accurately you want
      !              the solution to be computed.  They must be defined as
      !              program variables because the code may change them.  You
      !              have two choices --
      !                   Both RTOL and ATOL are scalars. (INFO(2)=0)
      !                   Both RTOL and ATOL are vectors. (INFO(2)=1)
      !              In either case all components must be non-negative.
      !
      !              The tolerances are used by the code in a local error test
      !              at each step which requires roughly that
      !                      ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
      !              for each vector component.
      !              (More specifically, a Euclidean norm is used to measure
      !              the size of vectors, and the error test uses the magnitude
      !              of the solution at the beginning of the step.)
      !
      !              The true (global) error is the difference between the true
      !              solution of the initial value problem and the computed
      !              approximation.  Practically all present day codes,
      !              including this one, control the local error at each step
      !              and do not even attempt to control the global error
      !              directly.  Roughly speaking, they produce a solution Y(T)
      !              which satisfies the differential equations with a
      !              residual R(T),    DY(T)/DT = DF(T,Y(T)) + R(T)   ,
      !              and, almost always, R(T) is bounded by the error
      !              tolerances.  Usually, but not always, the true accuracy of
      !              the computed Y is comparable to the error tolerances. This
      !              code will usually, but not always, deliver a more accurate
      !              solution if you reduce the tolerances and integrate again.
      !              By comparing two such solutions you can get a fairly
      !              reliable idea of the true error in the solution at the
      !              bigger tolerances.
      !
      !              Setting ATOL=0.D0 results in a pure relative error test on
      !              that component. Setting RTOL=0. results in a pure absolute
      !              error test on that component.  A mixed test with non-zero
      !              RTOL and ATOL corresponds roughly to a relative error
      !              test when the solution component is much bigger than ATOL
      !              and to an absolute error test when the solution component
      !              is smaller than the threshold ATOL.
      !
      !              Proper selection of the absolute error control parameters
      !              ATOL  requires you to have some idea of the scale of the
      !              solution components.  To acquire this information may mean
      !              that you will have to solve the problem more than once. In
      !              the absence of scale information, you should ask for some
      !              relative accuracy in all the components (by setting  RTOL
      !              values non-zero) and perhaps impose extremely small
      !              absolute error tolerances to protect against the danger of
      !              a solution component becoming zero.
      !
      !              The code will not attempt to compute a solution at an
      !              accuracy unreasonable for the machine being used.  It will
      !              advise you if you ask for too much accuracy and inform
      !              you as to the maximum accuracy it believes possible.
      !
      !       RWORK(*) -- Dimension this DOUBLE PRECISION work array of length
      !              LRW in your calling program.
      !
      !       RWORK(1) -- If you have set INFO(4)=0, you can ignore this
      !              optional input parameter.  Otherwise you must define a
      !              stopping point TSTOP by setting   RWORK(1) = TSTOP.
      !              (for some problems it may not be permissible to integrate
      !              past a point TSTOP because a discontinuity occurs there
      !              or the solution or its derivative is not defined beyond
      !              TSTOP.)
      !
      !       LRW -- Set it to the declared length of the RWORK array.
      !              You must have  LRW .GE. 130+21*NEQ
      !
      !       IWORK(*) -- Dimension this INTEGER work array of length LIW in
      !              your calling program.
      !
      !       LIW -- Set it to the declared length of the IWORK array.
      !              You must have  LIW .GE. 51
      !
      !       RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
      !              INTEGER type, respectively.  You can use them for
      !              communication between your program that calls DDEABM and
      !              the  DF subroutine.  They are not used or altered by
      !              DDEABM.  If you do not need RPAR or IPAR, ignore these
      !              parameters by treating them as dummy arguments.  If you do
      !              choose to use them, dimension them in your calling program
      !              and in DF as arrays of appropriate length.
      !
      !  **********************************************************************
      !  * OUTPUT -- After Any Return From DDEABM *
      !  **********************************************************************
      !
      !    The principal aim of the code is to return a computed solution at
      !    TOUT, although it is also possible to obtain intermediate results
      !    along the way.  To find out whether the code achieved its goal
      !    or if the integration process was interrupted before the task was
      !    completed, you must check the IDID parameter.
      !
      !
      !       T -- The solution was successfully advanced to the
      !              output value of T.
      !
      !       Y(*) -- Contains the computed solution approximation at T.
      !              You may also be interested in the approximate derivative
      !              of the solution at T.  It is contained in
      !              RWORK(21),...,RWORK(20+NEQ).
      !
      !       IDID -- Reports what the code did
      !
      !                          *** Task Completed ***
      !                    Reported by positive values of IDID
      !
      !              IDID = 1 -- A step was successfully taken in the
      !                        intermediate-output mode.  The code has not
      !                        yet reached TOUT.
      !
      !              IDID = 2 -- The integration to TOUT was successfully
      !                        completed (T=TOUT) by stepping exactly to TOUT.
      !
      !              IDID = 3 -- The integration to TOUT was successfully
      !                        completed (T=TOUT) by stepping past TOUT.
      !                        Y(*) is obtained by interpolation.
      !
      !                          *** Task Interrupted ***
      !                    Reported by negative values of IDID
      !
      !              IDID = -1 -- A large amount of work has been expended.
      !                        (500 steps attempted)
      !
      !              IDID = -2 -- The error tolerances are too stringent.
      !
      !              IDID = -3 -- The local error test cannot be satisfied
      !                        because you specified a zero component in ATOL
      !                        and the corresponding computed solution
      !                        component is zero.  Thus, a pure relative error
      !                        test is impossible for this component.
      !
      !              IDID = -4 -- The problem appears to be stiff.
      !
      !              IDID = -5,-6,-7,..,-32  -- Not applicable for this code
      !                        but used by other members of DEPAC or possible
      !                        future extensions.
      !
      !                          *** Task Terminated ***
      !                    Reported by the value of IDID=-33
      !
      !              IDID = -33 -- The code has encountered trouble from which
      !                        it cannot recover.  A message is printed
      !                        explaining the trouble and control is returned
      !                        to the calling program. For example, this occurs
      !                        when invalid input is detected.
      !
      !       RTOL, ATOL -- These quantities remain unchanged except when
      !              IDID = -2.  In this case, the error tolerances have been
      !              increased by the code to values which are estimated to be
      !              appropriate for continuing the integration.  However, the
      !              reported solution at T was obtained using the input values
      !              of RTOL and ATOL.
      !
      !       RWORK, IWORK -- Contain information which is usually of no
      !              interest to the user but necessary for subsequent calls.
      !              However, you may find use for
      !
      !              RWORK(11)--which contains the step size H to be
      !                         attempted on the next step.
      !
      !              RWORK(12)--if the tolerances have been increased by the
      !                         code (IDID = -2) , they were multiplied by the
      !                         value in RWORK(12).
      !
      !              RWORK(13)--Which contains the current value of the
      !                         independent variable, i.e. the farthest point
      !                         integration has reached. This will be different
      !                         from T only when interpolation has been
      !                         performed (IDID=3).
      !
      !              RWORK(20+I)--Which contains the approximate derivative
      !                         of the solution component Y(I).  In DDEABM, it
      !                         is obtained by calling subroutine DF to
      !                         evaluate the differential equation using T and
      !                         Y(*) when IDID=1 or 2, and by interpolation
      !                         when IDID=3.
      !
      !  **********************************************************************
      !  * INPUT -- What To Do To Continue The Integration *
      !  *             (calls after the first)             *
      !  **********************************************************************
      !
      !         This code is organized so that subsequent calls to continue the
      !         integration involve little (if any) additional effort on your
      !         part. You must monitor the IDID parameter in order to determine
      !         what to do next.
      !
      !         Recalling that the principal task of the code is to integrate
      !         from T to TOUT (the interval mode), usually all you will need
      !         to do is specify a new TOUT upon reaching the current TOUT.
      !
      !         Do not alter any quantity not specifically permitted below,
      !         in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or
      !         the differential equation in subroutine DF. Any such alteration
      !         constitutes a new problem and must be treated as such, i.e.
      !         you must start afresh.
      !
      !         You cannot change from vector to scalar error control or vice
      !         versa (INFO(2)) but you can change the size of the entries of
      !         RTOL, ATOL.  Increasing a tolerance makes the equation easier
      !         to integrate.  Decreasing a tolerance will make the equation
      !         harder to integrate and should generally be avoided.
      !
      !         You can switch from the intermediate-output mode to the
      !         interval mode (INFO(3)) or vice versa at any time.
      !
      !         If it has been necessary to prevent the integration from going
      !         past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
      !         code will not integrate to any TOUT beyond the currently
      !         specified TSTOP.  Once TSTOP has been reached you must change
      !         the value of TSTOP or set INFO(4)=0.  You may change INFO(4)
      !         or TSTOP at any time but you must supply the value of TSTOP in
      !         RWORK(1) whenever you set INFO(4)=1.
      !
      !         The parameter INFO(1) is used by the code to indicate the
      !         beginning of a new problem and to indicate whether integration
      !         is to be continued.  You must input the value  INFO(1) = 0
      !         when starting a new problem.  You must input the value
      !         INFO(1) = 1  if you wish to continue after an interrupted task.
      !         Do not set  INFO(1) = 0  on a continuation call unless you
      !         want the code to restart at the current T.
      !
      !                          *** Following A Completed Task ***
      !          If
      !              IDID = 1, call the code again to continue the integration
      !                      another step in the direction of TOUT.
      !
      !              IDID = 2 or 3, define a new TOUT and call the code again.
      !                      TOUT must be different from T. You cannot change
      !                      the direction of integration without restarting.
      !
      !                          *** Following An Interrupted Task ***
      !                      To show the code that you realize the task was
      !                      interrupted and that you want to continue, you
      !                      must take appropriate action and reset INFO(1) = 1
      !          If
      !              IDID = -1, the code has attempted 500 steps.
      !                      If you want to continue, set INFO(1) = 1 and
      !                      call the code again. An additional 500 steps
      !                      will be allowed.
      !
      !              IDID = -2, the error tolerances RTOL, ATOL have been
      !                      increased to values the code estimates appropriate
      !                      for continuing.  You may want to change them
      !                      yourself.  If you are sure you want to continue
      !                      with relaxed error tolerances, set INFO(1)=1 and
      !                      call the code again.
      !
      !              IDID = -3, a solution component is zero and you set the
      !                      corresponding component of ATOL to zero.  If you
      !                      are sure you want to continue, you must first
      !                      alter the error criterion to use positive values
      !                      for those components of ATOL corresponding to zero
      !                      solution components, then set INFO(1)=1 and call
      !                      the code again.
      !
      !              IDID = -4, the problem appears to be stiff.  It is very
      !                      inefficient to solve such problems with DDEABM.
      !                      The code DDEBDF in DEPAC handles this task
      !                      efficiently.  If you are absolutely sure you want
      !                      to continue with DDEABM, set INFO(1)=1 and call
      !                      the code again.
      !
      !              IDID = -5,-6,-7,..,-32  --- cannot occur with this code
      !                      but used by other members of DEPAC or possible
      !                      future extensions.
      !
      !                          *** Following A Terminated Task ***
      !          If
      !              IDID = -33, you cannot continue the solution of this
      !                      problem.  An attempt to do so will result in your
      !                      run being terminated.
      !
      !  **********************************************************************
      !  *Long Description:
      !
      !  **********************************************************************
      !  *             DEPAC Package Overview           *
      !  **********************************************************************
      !
      !  ....   You have a choice of three differential equation solvers from
      !  ....   DEPAC. The following brief descriptions are meant to aid you in
      !  ....   choosing the most appropriate code for your problem.
      !
      !  ....   DDERKF is a fifth order Runge-Kutta code. It is the simplest of
      !  ....   the three choices, both algorithmically and in the use of the
      !  ....   code. DDERKF is primarily designed to solve non-stiff and
      !  ....   mildly stiff differential equations when derivative evaluations
      !  ....   are not expensive. It should generally not be used to get high
      !  ....   accuracy results nor answers at a great many specific points.
      !  ....   Because DDERKF has very low overhead costs, it will usually
      !  ....   result in the least expensive integration when solving
      !  ....   problems requiring a modest amount of accuracy and having
      !  ....   equations that are not costly to evaluate. DDERKF attempts to
      !  ....   discover when it is not suitable for the task posed.
      !
      !  ....   DDEABM is a variable order (one through twelve) Adams code.
      !  ....   Its complexity lies somewhere between that of DDERKF and
      !  ....   DDEBDF.  DDEABM is primarily designed to solve non-stiff and
      !  ....   mildly stiff differential equations when derivative evaluations
      !  ....   are expensive, high accuracy results are needed or answers at
      !  ....   many specific points are required. DDEABM attempts to discover
      !  ....   when it is not suitable for the task posed.
      !
      !  ....   DDEBDF is a variable order (one through five) backward
      !  ....   differentiation formula code. it is the most complicated of
      !  ....   the three choices. DDEBDF is primarily designed to solve stiff
      !  ....   differential equations at crude to moderate tolerances.
      !  ....   If the problem is very stiff at all, DDERKF and DDEABM will be
      !  ....   quite inefficient compared to DDEBDF. However, DDEBDF will be
      !  ....   inefficient compared to DDERKF and DDEABM on non-stiff problems
      !  ....   because it uses much more storage, has a much larger overhead,
      !  ....   and the low order formulas will not give high accuracies
      !  ....   efficiently.
      !
      !  ....   The concept of stiffness cannot be described in a few words.
      !  ....   If you do not know the problem to be stiff, try either DDERKF
      !  ....   or DDEABM. Both of these codes will inform you of stiffness
      !  ....   when the cost of solving such problems becomes important.
      !
      !  *********************************************************************
      !
      ! ***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
      !                  oriented package of ODE solvers, Report SAND79-2374,
      !                  Sandia Laboratories, 1979.
      ! ***ROUTINES CALLED  DDES, XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891006  Cosmetic changes to prologue.  (WRB)
      !    891024  Changed references from DVNORM to DHVNRM.  (WRB)
      !    891024  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900510  Convert XERRWV calls to XERMSG calls.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DDEABM
      !
      INTEGER IALPHA, IBETA, IDELSN, IDID, IFOURU, IG, IHOLD,&
            INFO, IP, IPAR, IPHI, IPSI, ISIG, ITOLD, ITSTAR, ITWOU,&
            IV, IW, IWORK, IWT, IYP, IYPOUT, IYY, LIW, LRW, NEQ
      DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
      LOGICAL START,PHASE1,NORND,STIFF,INTOUT
      !
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),&
                RPAR(*),IPAR(*)
      !
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
      !
      EXTERNAL DF
      !
      !      CHECK FOR AN APPARENT INFINITE LOOP
      !
      ! ***FIRST EXECUTABLE STATEMENT  DDEABM
      IF ( INFO(1) .EQ. 0 ) IWORK(LIW) = 0
      IF (IWORK(LIW) .GE. 5) THEN
         IF (T .EQ. RWORK(21 + NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDEABM',&
               'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //&
               'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //&
               ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //&
               'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //&
               'CODE, PARTICULARLY INFO(1).', 13, 2)
            RETURN
         ENDIF
      ENDIF
      !
      !      CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
      !
      IDID=0
      IF (LRW .LT. 130+21*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE RWORK ' //&
            'ARRAY MUST BE AT LEAST 130 + 21*NEQ.$$' //&
            'YOU HAVE CALLED THE CODE WITH LRW = ' // XERN1, 1, 1)
         IDID=-33
      ENDIF
      !
      IF (LIW .LT. 51) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE IWORK ' //&
            'ARRAY MUST BE AT LEAST 51.$$YOU HAVE CALLED THE CODE ' //&
            'WITH LIW = ' // XERN1, 2, 1)
         IDID=-33
      ENDIF
      !
      !      COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
      !
      IYPOUT = 21
      ITSTAR = NEQ + 21
      IYP = 1 + ITSTAR
      IYY = NEQ + IYP
      IWT = NEQ + IYY
      IP = NEQ + IWT
      IPHI = NEQ + IP
      IALPHA = (NEQ*16) + IPHI
      IBETA = 12 + IALPHA
      IPSI = 12 + IBETA
      IV = 12 + IPSI
      IW = 12 + IV
      ISIG = 12 + IW
      IG = 13 + ISIG
      IGI = 13 + IG
      IXOLD = 11 + IGI
      IHOLD = 1 + IXOLD
      ITOLD = 1 + IHOLD
      IDELSN = 1 + ITOLD
      ITWOU = 1 + IDELSN
      IFOURU = 1 + ITWOU
      !
      RWORK(ITSTAR) = T
      IF (INFO(1) .EQ. 0) GO TO 50
      START = IWORK(21) .NE. (-1)
      PHASE1 = IWORK(22) .NE. (-1)
      NORND = IWORK(23) .NE. (-1)
      STIFF = IWORK(24) .NE. (-1)
      INTOUT = IWORK(25) .NE. (-1)
 !
 50   CALL DDES(DF,NEQ,T,Y,TOUT,INFO,RTOL,ATOL,IDID,RWORK(IYPOUT),&
               RWORK(IYP),RWORK(IYY),RWORK(IWT),RWORK(IP),RWORK(IPHI),&
               RWORK(IALPHA),RWORK(IBETA),RWORK(IPSI),RWORK(IV),&
               RWORK(IW),RWORK(ISIG),RWORK(IG),RWORK(IGI),RWORK(11),&
               RWORK(12),RWORK(13),RWORK(IXOLD),RWORK(IHOLD),&
               RWORK(ITOLD),RWORK(IDELSN),RWORK(1),RWORK(ITWOU),&
               RWORK(IFOURU),START,PHASE1,NORND,STIFF,INTOUT,IWORK(26),&
               IWORK(27),IWORK(28),IWORK(29),IWORK(30),IWORK(31),&
               IWORK(32),IWORK(33),IWORK(34),IWORK(35),IWORK(45),&
               RPAR,IPAR)
      !
      IWORK(21) = -1
      IF (START) IWORK(21) = 1
      IWORK(22) = -1
      IF (PHASE1) IWORK(22) = 1
      IWORK(23) = -1
      IF (NORND) IWORK(23) = 1
      IWORK(24) = -1
      IF (STIFF) IWORK(24) = 1
      IWORK(25) = -1
      IF (INTOUT) IWORK(25) = 1
      !
      IF (IDID .NE. (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T .NE. RWORK(ITSTAR)) IWORK(LIW) = 0
      !
      RETURN
      END
      ! DECK DDES
      SUBROUTINE DDES (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,&
         YPOUT, YP, YY, WT, P, PHI, ALPHA, BETA, PSI, V, W, SIG, G, GI,&
         H, EPS, X, XOLD, HOLD, TOLD, DELSGN, TSTOP, TWOU, FOURU, START,&
         PHASE1, NORND, STIFF, INTOUT, NS, KORD, KOLD, INIT, KSTEPS,&
         KLE4, IQUIT, KPREV, IVC, IV, KGI, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DDES
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEABM
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (DES-S, DDES-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DDEABM merely allocates storage for DDES to relieve the user of the
      !    inconvenience of a long call list.  Consequently  DDES  is used as
      !    described in the comments for  DDEABM .
      !
      ! ***SEE ALSO  DDEABM
      ! ***ROUTINES CALLED  D1MACH, DINTP, DSTEPS, XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    900510  Convert XERRWV calls to XERMSG calls, cvt GOTOs to
      !            IF-THEN-ELSE.  (RWC)
      !    910722  Updated AUTHOR section.  (ALS)
      ! ***END PROLOGUE  DDES
      !
      INTEGER IDID, INFO, INIT, IPAR, IQUIT, IV, IVC, K, KGI, KLE4,&
            KOLD, KORD, KPREV, KSTEPS, L, LTOL, MAXNUM, NATOLP, NEQ,&
            NRTOLP, NS
      DOUBLE PRECISION A, ABSDEL, ALPHA, ATOL, BETA, D1MACH,&
            DEL, DELSGN, DT, EPS, FOURU, G, GI, H,&
            HA, HOLD, P, PHI, PSI, RPAR, RTOL, SIG, T, TOLD, TOUT,&
            TSTOP, TWOU, U, V, W, WT, X, XOLD, Y, YP, YPOUT, YY
      LOGICAL STIFF,CRASH,START,PHASE1,NORND,INTOUT
      !
      DIMENSION Y(*),YY(*),WT(*),PHI(NEQ,16),P(*),YP(*),&
        YPOUT(*),PSI(12),ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),&
        GI(11),IV(10),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
      !
      EXTERNAL DF
      !
      ! .......................................................................
      !
      !   THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
      !   NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
      !   IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
      !   WORK.
      !
      SAVE MAXNUM
      DATA MAXNUM/500/
      !
      ! .......................................................................
      !
      ! ***FIRST EXECUTABLE STATEMENT  DDES
      IF (INFO(1) .EQ. 0) THEN
         !
         !  ON THE FIRST CALL , PERFORM INITIALIZATION --
         !         DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
         !         FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE
         !         VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
         !
         U=D1MACH(4)
         !                        -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
         TWOU=2.D0*U
         FOURU=4.D0*U
         !                        -- SET TERMINATION FLAG
         IQUIT=0
         !                        -- SET INITIALIZATION INDICATOR
         INIT=0
         !                        -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS=0
         !                        -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT= .FALSE.
         !                        -- SET INDICATOR FOR STIFFNESS DETECTION
         STIFF= .FALSE.
         !                        -- SET STEP COUNTER FOR STIFFNESS DETECTION
         KLE4=0
         !                        -- SET INDICATORS FOR STEPS CODE
         START= .TRUE.
         PHASE1= .TRUE.
         NORND= .TRUE.
         !                        -- RESET INFO(1) FOR SUBSEQUENT CALLS
         INFO(1)=1
      ENDIF
      !
      ! .......................................................................
      !
      !       CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
      !
      IF (INFO(1) .NE. 0  .AND.  INFO(1) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(1) MUST BE ' //&
            'SET TO 0 FOR THE START OF A NEW PROBLEM, AND MUST BE ' //&
            'SET TO 1 FOLLOWING AN INTERRUPTED TASK.  YOU ARE ' //&
            'ATTEMPTING TO CONTINUE THE INTEGRATION ILLEGALLY BY ' //&
            'CALLING THE CODE WITH INFO(1) = ' // XERN1, 3, 1)
         IDID=-33
      ENDIF
      !
      IF (INFO(2) .NE. 0  .AND.  INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(2) MUST BE ' //&
            '0 OR 1 INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //&
            'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //&
            XERN1, 4, 1)
         IDID=-33
      ENDIF
      !
      IF (INFO(3) .NE. 0  .AND.  INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(3) MUST BE ' //&
            '0 OR 1 INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT ' //&
            'MODE OF INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED ' //&
            'THE CODE WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID=-33
      ENDIF
      !
      IF (INFO(4) .NE. 0  .AND.  INFO(4) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(4)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(4) MUST BE ' //&
            '0 OR 1 INDICATING WHETHER OR NOT THE INTEGRATION ' //&
            'INTERVAL IS TO BE RESTRICTED BY A POINT TSTOP.  YOU ' //&
            'HAVE CALLED THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
         IDID=-33
      ENDIF
      !
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM,  THE NUMBER OF ' //&
            'EQUATIONS NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE ' //&
            'CALLED THE CODE WITH  NEQ = ' // XERN1, 6, 1)
         IDID=-33
      ENDIF
      !
      NRTOLP = 0
      NATOLP = 0
      DO 90 K=1,NEQ
         IF (NRTOLP .EQ. 0 .AND. RTOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') RTOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE RELATIVE ' //&
               'ERROR TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU ' //&
               'HAVE CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' //&
               XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //&
               'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
            IDID = -33
            NRTOLP = 1
         ENDIF
         !
         IF (NATOLP .EQ. 0 .AND. ATOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE ABSOLUTE ' //&
               'ERROR TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU ' //&
               'HAVE CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' //&
               XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //&
               'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID = -33
            NATOLP = 1
         ENDIF
         !
         IF (INFO(2) .EQ. 0) GO TO 100
         IF (NATOLP.GT.0 .AND. NRTOLP.GT.0) GO TO 100
   90 CONTINUE
  !
  100 IF (INFO(4) .EQ. 1) THEN
         IF (SIGN(1.D0,TOUT-T) .NE. SIGN(1.D0,TSTOP-T)&
            .OR. ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CALLED THE CODE WITH  TOUT = ' // XERN3 // ' BUT ' //&
               'YOU HAVE ALSO TOLD THE CODE (INFO(4) = 1) NOT TO ' //&
               'INTEGRATE PAST THE POINT TSTOP = ' // XERN4 //&
               ' THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
      !
      !      CHECK SOME CONTINUATION POSSIBILITIES
      !
      IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CALLED THE CODE WITH  T = TOUT = ' // XERN3 //&
               '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 9, 1)
            IDID=-33
         ENDIF
         !
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CHANGED THE VALUE OF T FROM ' // XERN3 // ' TO ' //&
               XERN4 //'  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',&
               10, 1)
            IDID=-33
         ENDIF
         !
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, BY ' //&
                  'CALLING THE CODE WITH TOUT = ' // XERN3 //&
                  ' YOU ARE ATTEMPTING TO CHANGE THE DIRECTION OF ' //&
                  'INTEGRATION.$$THIS IS NOT ALLOWED WITHOUT ' //&
                  'RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
      !
      !      INVALID INPUT DETECTED
      !
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
            IQUIT = -33
            INFO(1) = -1
         ELSE
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INVALID ' //&
               'INPUT WAS DETECTED ON SUCCESSIVE ENTRIES.  IT IS ' //&
               'IMPOSSIBLE TO PROCEED BECAUSE YOU HAVE NOT ' //&
               'CORRECTED THE PROBLEM, SO EXECUTION IS BEING ' //&
               'TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
      !
      ! .......................................................................
      !
      !      RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
      !      ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
      !      THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
      !      FOURU WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
      !
      DO 180 K=1,NEQ
        IF (RTOL(K)+ATOL(K) .GT. 0.D0) GO TO 170
        RTOL(K)=FOURU
        IDID=-2
  170   IF (INFO(2) .EQ. 0) GO TO 190
  180   CONTINUE
  !
  190 IF (IDID .NE. (-2)) GO TO 200
      !                        RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
      !                                                 SMALL POSITIVE VALUE
      INFO(1)=-1
      RETURN
  !
  !      BRANCH ON STATUS OF INITIALIZATION INDICATOR
  !             INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
  !                    AND DIRECTION NOT YET SET
  !             INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
  !             INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
  !
  200 IF (INIT .EQ. 0) GO TO 210
      IF (INIT .EQ. 1) GO TO 220
      GO TO 240
  !
  ! .......................................................................
  !
  !      MORE INITIALIZATION --
  !                          -- EVALUATE INITIAL DERIVATIVES
  !
  210 INIT=1
      A=T
      CALL DF(A,Y,YP,RPAR,IPAR)
      IF (T .NE. TOUT) GO TO 220
      IDID=2
      DO 215 L = 1,NEQ
  215    YPOUT(L) = YP(L)
      TOLD=T
      RETURN
  !
  !                          -- SET INDEPENDENT AND DEPENDENT VARIABLES
  !                                               X AND YY(*) FOR STEPS
  !                          -- SET SIGN OF INTEGRATION DIRECTION
  !                          -- INITIALIZE THE STEP SIZE
  !
  220 INIT = 2
      X = T
      DO 230 L = 1,NEQ
  230   YY(L) = Y(L)
      DELSGN = SIGN(1.0D0,TOUT-T)
      H = SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X)
  !
  ! .......................................................................
  !
  !    ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
  !    OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
  !
  240 DEL = TOUT - T
      ABSDEL = ABS(DEL)
  !
  ! .......................................................................
  !
  !    IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
  !
  250 IF(ABS(X-T) .LT. ABSDEL) GO TO 260
      CALL DINTP(X,YY,TOUT,Y,YPOUT,NEQ,KOLD,PHI,IVC,IV,KGI,GI,&
                                              ALPHA,G,W,XOLD,P)
      IDID = 3
      IF (X .NE. TOUT) GO TO 255
      IDID = 2
      INTOUT = .FALSE.
  255 T = TOUT
      TOLD = T
      RETURN
  !
  !    IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
  !    EXTRAPOLATE AND RETURN
  !
  260 IF (INFO(4) .NE. 1) GO TO 280
      IF (ABS(TSTOP-X) .GE. FOURU*ABS(X)) GO TO 280
      DT = TOUT - X
      DO 270 L = 1,NEQ
  270   Y(L) = YY(L) + DT*YP(L)
      CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
      IDID = 3
      T = TOUT
      TOLD = T
      RETURN
  !
  280 IF (INFO(3) .EQ. 0  .OR.  .NOT.INTOUT) GO TO 300
      !
      !    INTERMEDIATE-OUTPUT MODE
      !
      IDID = 1
      DO 290 L = 1,NEQ
        Y(L)=YY(L)
  290   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INTOUT = .FALSE.
      RETURN
  !
  ! .......................................................................
  !
  !      MONITOR NUMBER OF STEPS ATTEMPTED
  !
  300 IF (KSTEPS .LE. MAXNUM) GO TO 330
      !
      !                        A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      IDID=-1
      KSTEPS=0
      IF (.NOT. STIFF) GO TO 310
      !
      !                        PROBLEM APPEARS TO BE STIFF
      IDID=-4
      STIFF= .FALSE.
      KLE4=0
  !
  310 DO 320 L = 1,NEQ
        Y(L) = YY(L)
  320   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
  !
  ! .......................................................................
  !
  !    LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP
  !
  330 HA = ABS(H)
      IF (INFO(4) .NE. 1) GO TO 340
      HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
      EPS = 1.0D0
      LTOL = 1
      DO 350 L = 1,NEQ
        IF (INFO(2) .EQ. 1) LTOL = L
        WT(L) = RTOL(LTOL)*ABS(YY(L)) + ATOL(LTOL)
        IF (WT(L) .LE. 0.0D0) GO TO 360
  350   CONTINUE
      GO TO 380
  !
  !                        RELATIVE ERROR CRITERION INAPPROPRIATE
  360 IDID = -3
      DO 370 L = 1,NEQ
        Y(L) = YY(L)
  370   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
  !
  380 CALL DSTEPS(DF,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,CRASH,PHI,P,&
                 YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,NORND,KSTEPS,&
                 TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,RPAR,IPAR)
      !
      ! .......................................................................
      !
      IF(.NOT.CRASH) GO TO 420
      !
      !                        TOLERANCES TOO SMALL
      IDID = -2
      RTOL(1) = EPS*RTOL(1)
      ATOL(1) = EPS*ATOL(1)
      IF (INFO(2) .EQ. 0) GO TO 400
      DO 390 L = 2,NEQ
        RTOL(L) = EPS*RTOL(L)
  390   ATOL(L) = EPS*ATOL(L)
  400 DO 410 L = 1,NEQ
        Y(L) = YY(L)
  410   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
  !
  !    (STIFFNESS TEST) COUNT NUMBER OF CONSECUTIVE STEPS TAKEN WITH THE
  !    ORDER OF THE METHOD BEING LESS OR EQUAL TO FOUR
  !
  420 KLE4 = KLE4 + 1
      IF(KOLD .GT. 4) KLE4 = 0
      IF(KLE4 .GE. 50) STIFF = .TRUE.
      INTOUT = .TRUE.
      GO TO 250
      END
      ! DECK DINTP
      SUBROUTINE DINTP (X, Y, XOUT, YOUT, YPOUT, NEQN, KOLD, PHI, IVC,&
         IV, KGI, GI, ALPHA, OG, OW, OX, OY)
      ! ***BEGIN PROLOGUE  DINTP
      ! ***PURPOSE  Approximate the solution at XOUT by evaluating the
      !             polynomial computed in DSTEPS at XOUT.  Must be used in
      !             conjunction with DSTEPS.
      ! ***LIBRARY   SLATEC (DEPAC)
      ! ***CATEGORY  I1A1B
      ! ***TYPE      DOUBLE PRECISION (SINTRP-S, DINTP-D)
      ! ***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
      !              ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
      !              SMOOTH INTERPOLANT
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    The methods in subroutine  DSTEPS  approximate the solution near  X
      !    by a polynomial.  Subroutine  DINTP  approximates the solution at
      !    XOUT  by evaluating the polynomial there.  Information defining this
      !    polynomial is passed from  DSTEPS  so  DINTP  cannot be used alone.
      !
      !    Subroutine DSTEPS is completely explained and documented in the text
      !    "Computer Solution of Ordinary Differential Equations, the Initial
      !    Value Problem"  by L. F. Shampine and M. K. Gordon.
      !
      !    Input to DINTP --
      !
      !    The user provides storage in the calling program for the arrays in
      !    the call list
      !       DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
      !                 AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
      !    and defines
      !       XOUT -- point at which solution is desired.
      !    The remaining parameters are defined in  DSTEPS  and passed to
      !    DINTP  from that subroutine
      !
      !    Output from  DINTP --
      !
      !       YOUT(*) -- solution at  XOUT
      !       YPOUT(*) -- derivative of solution at  XOUT
      !    The remaining parameters are returned unaltered from their input
      !    values.  Integration with  DSTEPS  may be continued.
      !
      ! ***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
      !                  II, Report SAND84-0293, Sandia Laboratories, 1984.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    840201  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DINTP
      !
      INTEGER I, IQ, IV, IVC, IW, J, JQ, KGI, KOLD, KP1, KP2,&
              L, M, NEQN
      DOUBLE PRECISION ALP, ALPHA, C, G, GDI, GDIF, GI, GAMMA, H, HI,&
             HMU, OG, OW, OX, OY, PHI, RMU, SIGMA, TEMP1, TEMP2, TEMP3,&
             W, X, XI, XIM1, XIQ, XOUT, Y, YOUT, YPOUT
      !
      DIMENSION Y(*),YOUT(*),YPOUT(*),PHI(NEQN,16),OY(*)
      DIMENSION G(13),C(13),W(13),OG(13),OW(12),ALPHA(12),GI(11),IV(10)
      !
      ! ***FIRST EXECUTABLE STATEMENT  DINTP
      KP1 = KOLD + 1
      KP2 = KOLD + 2
      !
      HI = XOUT - OX
      H = X - OX
      XI = HI/H
      XIM1 = XI - 1.D0
      !
      !    INITIALIZE W(*) FOR COMPUTING G(*)
      !
      XIQ = XI
      DO 10 IQ = 1,KP1
        XIQ = XI*XIQ
        TEMP1 = IQ*(IQ+1)
 10     W(IQ) = XIQ/TEMP1
      !
      !    COMPUTE THE DOUBLE INTEGRAL TERM GDI
      !
      IF (KOLD .LE. KGI) GO TO 50
      IF (IVC .GT. 0) GO TO 20
      GDI = 1.0D0/TEMP1
      M = 2
      GO TO 30
 20   IW = IV(IVC)
      GDI = OW(IW)
      M = KOLD - IW + 3
 30   IF (M .GT. KOLD) GO TO 60
      DO 40 I = M,KOLD
 40     GDI = OW(KP2-I) - ALPHA(I)*GDI
      GO TO 60
 50   GDI = GI(KOLD)
 !
 !    COMPUTE G(*) AND C(*)
 !
 60   G(1) = XI
      G(2) = 0.5D0*XI*XI
      C(1) = 1.0D0
      C(2) = XI
      IF (KOLD .LT. 2) GO TO 90
      DO 80 I = 2,KOLD
        ALP = ALPHA(I)
        GAMMA = 1.0D0 + XIM1*ALP
        L = KP2 - I
        DO 70 JQ = 1,L
 70       W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
        G(I+1) = W(1)
 80     C(I+1) = GAMMA*C(I)
 !
 !    DEFINE INTERPOLATION PARAMETERS
 !
 90   SIGMA = (W(2) - XIM1*W(1))/GDI
      RMU = XIM1*C(KP1)/GDI
      HMU = RMU/H
      !
      !    INTERPOLATE FOR THE SOLUTION -- YOUT
      !    AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
      !
      DO 100 L = 1,NEQN
        YOUT(L) = 0.0D0
 100    YPOUT(L) = 0.0D0
      DO 120 J = 1,KOLD
        I = KP2 - J
        GDIF = OG(I) - OG(I-1)
        TEMP2 = (G(I) - G(I-1)) - SIGMA*GDIF
        TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
        DO 110 L = 1,NEQN
          YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 110      YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 120    CONTINUE
      DO 130 L = 1,NEQN
        YOUT(L) = ((1.0D0 - SIGMA)*OY(L) + SIGMA*Y(L)) +&
                   H*(YOUT(L) + (G(1) - SIGMA*OG(1))*PHI(L,1))
 130    YPOUT(L) = HMU*(OY(L) - Y(L)) +&
                      (YPOUT(L) + (C(1) + RMU*OG(1))*PHI(L,1))
      !
      RETURN
      END
      ! DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
      ! ***BEGIN PROLOGUE  XERMSG
      ! ***PURPOSE  Process error messages for SLATEC and other libraries.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3C
      ! ***TYPE      ALL (XERMSG-A)
      ! ***KEYWORDS  ERROR MESSAGE, XERROR
      ! ***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
      ! ***DESCRIPTION
      !
      !    XERMSG processes a diagnostic message in a manner determined by the
      !    value of LEVEL and the current value of the library error control
      !    flag, KONTRL.  See subroutine XSETF for details.
      !
      !     LIBRAR   A character constant (or character variable) with the name
      !              of the library.  This will be 'SLATEC' for the SLATEC
      !              Common Math Library.  The error handling package is
      !              general enough to be used by many libraries
      !              simultaneously, so it is desirable for the routine that
      !              detects and reports an error to identify the library name
      !              as well as the routine name.
      !
      !     SUBROU   A character constant (or character variable) with the name
      !              of the routine that detected the error.  Usually it is the
      !              name of the routine that is calling XERMSG.  There are
      !              some instances where a user callable library routine calls
      !              lower level subsidiary routines where the error is
      !              detected.  In such cases it may be more informative to
      !              supply the name of the routine the user called rather than
      !              the name of the subsidiary routine that detected the
      !              error.
      !
      !     MESSG    A character constant (or character variable) with the text
      !              of the error or warning message.  In the example below,
      !              the message is a character constant that contains a
      !              generic message.
      !
      !                    CALL XERMSG ('SLATEC', 'MMPY',
      !                   *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
      !                   *3, 1)
      !
      !              It is possible (and is sometimes desirable) to generate a
      !              specific message--e.g., one that contains actual numeric
      !              values.  Specific numeric values can be converted into
      !              character strings using formatted WRITE statements into
      !              character variables.  This is called standard Fortran
      !              internal file I/O and is exemplified in the first three
      !              lines of the following example.  You can also catenate
      !              substrings of characters to construct the error message.
      !              Here is an example showing the use of both writing to
      !              an internal file and catenating character strings.
      !
      !                    CHARACTER*5 CHARN, CHARL
      !                    WRITE (CHARN,10) N
      !                    WRITE (CHARL,10) LDA
      !                 10 FORMAT(I5)
      !                    CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
      !                   *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
      !                   *   CHARL, 3, 1)
      !
      !              There are two subtleties worth mentioning.  One is that
      !              the // for character catenation is used to construct the
      !              error message so that no single character constant is
      !              continued to the next line.  This avoids confusion as to
      !              whether there are trailing blanks at the end of the line.
      !              The second is that by catenating the parts of the message
      !              as an actual argument rather than encoding the entire
      !              message into one large character variable, we avoid
      !              having to know how long the message will be in order to
      !              declare an adequate length for that large character
      !              variable.  XERMSG calls XERPRN to print the message using
      !              multiple lines if necessary.  If the message is very long,
      !              XERPRN will break it into pieces of 72 characters (as
      !              requested by XERMSG) for printing on multiple lines.
      !              Also, XERMSG asks XERPRN to prefix each line with ' *  '
      !              so that the total line length could be 76 characters.
      !              Note also that XERPRN scans the error message backwards
      !              to ignore trailing blanks.  Another feature is that
      !              the substring '$$' is treated as a new line sentinel
      !              by XERPRN.  If you want to construct a multiline
      !              message without having to count out multiples of 72
      !              characters, just use '$$' as a separator.  '$$'
      !              obviously must occur within 72 characters of the
      !              start of each line to have its intended effect since
      !              XERPRN is asked to wrap around at 72 characters in
      !              addition to looking for '$$'.
      !
      !     NERR     An integer value that is chosen by the library routine's
      !              author.  It must be in the range -99 to 999 (three
      !              printable digits).  Each distinct error should have its
      !              own error number.  These error numbers should be described
      !              in the machine readable documentation for the routine.
      !              The error numbers need be unique only within each routine,
      !              so it is reasonable for each routine to start enumerating
      !              errors from 1 and proceeding to the next integer.
      !
      !     LEVEL    An integer value in the range 0 to 2 that indicates the
      !              level (severity) of the error.  Their meanings are
      !
      !             -1  A warning message.  This is used if it is not clear
      !                 that there really is an error, but the user's attention
      !                 may be needed.  An attempt is made to only print this
      !                 message once.
      !
      !              0  A warning message.  This is used if it is not clear
      !                 that there really is an error, but the user's attention
      !                 may be needed.
      !
      !              1  A recoverable error.  This is used even if the error is
      !                 so serious that the routine cannot return any useful
      !                 answer.  If the user has told the error package to
      !                 return after recoverable errors, then XERMSG will
      !                 return to the Library routine which can then return to
      !                 the user's routine.  The user may also permit the error
      !                 package to terminate the program upon encountering a
      !                 recoverable error.
      !
      !              2  A fatal error.  XERMSG will not return to its caller
      !                 after it receives a fatal error.  This level should
      !                 hardly ever be used; it is much better to allow the
      !                 user a chance to recover.  An example of one of the few
      !                 cases in which it is permissible to declare a level 2
      !                 error is a reverse communication Library routine that
      !                 is likely to be called repeatedly until it integrates
      !                 across some interval.  If there is a serious error in
      !                 the input such that another step cannot be taken and
      !                 the Library routine is called again without the input
      !                 error having been corrected by the caller, the Library
      !                 routine will probably be called forever with improper
      !                 input.  In this case, it is reasonable to declare the
      !                 error to be fatal.
      !
      !     Each of the arguments to XERMSG is input; none will be modified by
      !     XERMSG.  A routine may make multiple calls to XERMSG with warning
      !     level messages; however, after a call to XERMSG with a recoverable
      !     error, the routine should return to the user.  Do not try to call
      !     XERMSG with a second recoverable error after the first recoverable
      !     error because the error package saves the error number.  The user
      !     can retrieve this error number by calling another entry point in
      !     the error handling package and then clear the error number when
      !     recovering from the error.  Calling XERMSG in succession causes the
      !     old error number to be overwritten by the latest error number.
      !     This is considered harmless for error numbers associated with
      !     warning messages but must not be done for error numbers of serious
      !     errors.  After a call to XERMSG with a recoverable error, the user
      !     must be given a chance to call NUMXER or XERCLR to retrieve or
      !     clear the error number.
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
      ! ***REVISION HISTORY  (YYMMDD)
      !    880101  DATE WRITTEN
      !    880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
      !            THERE ARE TWO BASIC CHANGES.
      !            1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
      !                PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
      !                INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
      !                ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
      !                ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
      !                ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
      !                72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
      !                LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
      !            2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
      !                FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
      !                OF LOWER CASE.
      !    880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
      !            THE PRINCIPAL CHANGES ARE
      !            1.  CLARIFY COMMENTS IN THE PROLOGUES
      !            2.  RENAME XRPRNT TO XERPRN
      !            3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
      !                SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
      !                CHARACTER FOR NEW RECORDS.
      !    890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
      !            CLEAN UP THE CODING.
      !    890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
      !            PREFIX.
      !    891013  REVISED TO CORRECT COMMENTS.
      !    891214  Prologue converted to Version 4.0 format.  (WRB)
      !    900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
      !            NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
      !            LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
      !            XERCTL to XERCNT.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
      ! ***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
      !
      !        LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
      !        MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
      !           SHOULD BE PRINTED.
      !
      !        WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
      !           CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
      !           AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
      !
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.&
         LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //&
            'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//&
            'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
      !
      !        RECORD THE MESSAGE.
      !
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
      !
      !        HANDLE PRINT-ONCE WARNING MESSAGES.
      !
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
      !
      !        ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
      !
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
      !
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
      !
      !        SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
      !        ZERO AND THE ERROR IS NOT FATAL.
      !
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
      !
      !        ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
      !        MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
      !        AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
      !        IS NOT ZERO.
      !
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
      !
      !        IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
      !        PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
      !        FROM EACH OF THE FOLLOWING THREE OPTIONS.
      !        1.  LEVEL OF THE MESSAGE
      !               'INFORMATIVE MESSAGE'
      !               'POTENTIALLY RECOVERABLE ERROR'
      !               'FATAL ERROR'
      !        2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
      !               'PROG CONTINUES'
      !               'PROG ABORTED'
      !        3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
      !            MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
      !            WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
      !               'TRACEBACK REQUESTED'
      !               'TRACEBACK NOT REQUESTED'
      !        NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
      !        EXCEED 74 CHARACTERS.
      !        WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
      !
      IF (LKNTRL .GT. 0) THEN
         !
         !        THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
         !
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
         !
         !        THEN WHETHER THE PROGRAM WILL CONTINUE.
         !
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.&
             (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
         !
         !        FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
         !
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
      !
      !        NOW SEND OUT THE MESSAGE.
      !
      CALL XERPRN (' *  ', -1, MESSG, 72)
      !
      !        IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
      !           TRACEBACK.
      !
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
   !
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
      !
      !        IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
      !
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
   !
   !        IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
   !        CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
   !
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
      !
      !        THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
      !        FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
      !        SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
      !
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN&
               (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
      ! DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
      ! ***BEGIN PROLOGUE  XERPRN
      ! ***SUBSIDIARY
      ! ***PURPOSE  Print error messages processed by XERMSG.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3C
      ! ***TYPE      ALL (XERPRN-A)
      ! ***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
      ! ***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
      ! ***DESCRIPTION
      !
      !  This routine sends one or more lines to each of the (up to five)
      !  logical units to which error messages are to be sent.  This routine
      !  is called several times by XERMSG, sometimes with a single line to
      !  print and sometimes with a (potentially very long) message that may
      !  wrap around into multiple lines.
      !
      !  PREFIX  Input argument of type CHARACTER.  This argument contains
      !          characters to be put at the beginning of each line before
      !          the body of the message.  No more than 16 characters of
      !          PREFIX will be used.
      !
      !  NPREF   Input argument of type INTEGER.  This argument is the number
      !          of characters to use from PREFIX.  If it is negative, the
      !          intrinsic function LEN is used to determine its length.  If
      !          it is zero, PREFIX is not used.  If it exceeds 16 or if
      !          LEN(PREFIX) exceeds 16, only the first 16 characters will be
      !          used.  If NPREF is positive and the length of PREFIX is less
      !          than NPREF, a copy of PREFIX extended with blanks to length
      !          NPREF will be used.
      !
      !  MESSG   Input argument of type CHARACTER.  This is the text of a
      !          message to be printed.  If it is a long message, it will be
      !          broken into pieces for printing on multiple lines.  Each line
      !          will start with the appropriate prefix and be followed by a
      !          piece of the message.  NWRAP is the number of characters per
      !          piece; that is, after each NWRAP characters, we break and
      !          start a new line.  In addition the characters '$$' embedded
      !          in MESSG are a sentinel for a new line.  The counting of
      !          characters up to NWRAP starts over for each new line.  The
      !          value of NWRAP typically used by XERMSG is 72 since many
      !          older error messages in the SLATEC Library are laid out to
      !          rely on wrap-around every 72 characters.
      !
      !  NWRAP   Input argument of type INTEGER.  This gives the maximum size
      !          piece into which to break MESSG for printing on multiple
      !          lines.  An embedded '$$' ends a line, and the count restarts
      !          at the following character.  If a line break does not occur
      !          on a blank (it would split a word) that word is moved to the
      !          next line.  Values of NWRAP less than 16 will be treated as
      !          16.  Values of NWRAP greater than 132 will be treated as 132.
      !          The actual line length will be NPREF + NWRAP after NPREF has
      !          been adjusted to fall between 0 and 16 and NWRAP has been
      !          adjusted to fall between 16 and 132.
      !
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  I1MACH, XGETUA
      ! ***REVISION HISTORY  (YYMMDD)
      !    880621  DATE WRITTEN
      !    880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
      !            JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
      !            THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
      !            SLASH CHARACTER IN FORMAT STATEMENTS.
      !    890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
      !            STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
      !            LINES TO BE PRINTED.
      !    890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
      !            CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
      !    891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
      !    891214  Prologue converted to Version 4.0 format.  (WRB)
      !    900510  Added code to break messages between words.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
      ! ***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
      !
      !        A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
      !        ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
      !        ERROR MESSAGE UNIT.
      !
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
      !
      !        LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
      !        BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
      !        THE REST OF THIS ROUTINE.
      !
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
      !
      !        LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
      !        TIME FROM MESSG TO PRINT ON ONE LINE.
      !
      LWRAP = MAX(16, MIN(132, NWRAP))
      !
      !        SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
      !
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
      !
      !        IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
      !
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
      !
      !        SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
      !        STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
      !        WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
      !        WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
      !
      !        WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
      !        INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
      !        OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
      !        OF THE SECOND ARGUMENT.
      !
      !        THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
      !        FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
      !        OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
      !        POSITION NEXTC.
      !
      !        LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
      !                        REMAINDER OF THE CHARACTER STRING.  LPIECE
      !                        SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
      !                        WHICHEVER IS LESS.
      !
      !        LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
      !                        NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
      !                        PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
      !                        BLANK LINES.  THIS TAKES CARE OF THE SITUATION
      !                        WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
      !                        EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
      !                        SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
      !                        SHOULD BE INCREMENTED BY 2.
      !
      !        LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
      !
      !        ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
      !                        RESET LPIECE = LPIECE-1.  NOTE THAT THIS
      !                        PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
      !                        LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
      !                        AT THE END OF A LINE.
      !
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
         !
         !        THERE WAS NO NEW LINE SENTINEL FOUND.
         !
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
         !
         !        WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
         !        DON'T PRINT A BLANK LINE.
         !
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
         !
         !        LPIECE SHOULD BE SET DOWN TO LWRAP.
         !
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
         !
         !        IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
         !        WE SHOULD DECREMENT LPIECE BY ONE.
         !
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
      !
      !        PRINT
      !
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
      !
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
      ! DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,&
         ICOUNT)
      ! ***BEGIN PROLOGUE  XERSVE
      ! ***SUBSIDIARY
      ! ***PURPOSE  Record that an error has occurred.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3
      ! ***TYPE      ALL (XERSVE-A)
      ! ***KEYWORDS  ERROR, XERROR
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !  *Usage:
      !
      !         INTEGER  KFLAG, NERR, LEVEL, ICOUNT
      !         CHARACTER * (len) LIBRAR, SUBROU, MESSG
      !
      !         CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
      !
      !  *Arguments:
      !
      !         LIBRAR :IN    is the library that the message is from.
      !         SUBROU :IN    is the subroutine that the message is from.
      !         MESSG  :IN    is the message to be saved.
      !         KFLAG  :IN    indicates the action to be performed.
      !                       when KFLAG > 0, the message in MESSG is saved.
      !                       when KFLAG=0 the tables will be dumped and
      !                       cleared.
      !                       when KFLAG < 0, the tables will be dumped and
      !                       not cleared.
      !         NERR   :IN    is the error number.
      !         LEVEL  :IN    is the error severity.
      !         ICOUNT :OUT   the number of times this message has been seen,
      !                       or zero if the table has overflowed and does not
      !                       contain this message specifically.  When KFLAG=0,
      !                       ICOUNT will not be altered.
      !
      !  *Description:
      !
      !    Record that this error occurred and possibly dump and clear the
      !    tables.
      !
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  I1MACH, XGETUA
      ! ***REVISION HISTORY  (YYMMDD)
      !    800319  DATE WRITTEN
      !    861211  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900413  Routine modified to remove reference to KFLAG.  (WRB)
      !    900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
      !            sequence, use IF-THEN-ELSE, make number of saved entries
      !            easily changeable, changed routine name from XERSAV to
      !            XERSVE.  (RWC)
      !    910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
      ! ***FIRST EXECUTABLE STATEMENT  XERSVE
      !
      IF (KFLAG.LE.0) THEN
         !
         !         Dump the table.
         !
         IF (NMSG.EQ.0) RETURN
         !
         !         Print to each unit.
         !
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            !
            !            Print the table header.
            !
            WRITE (IUNIT,9000)
            !
            !            Print body of table.
            !
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),&
                  NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
            !
            !            Print number of other errors.
            !
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
         !
         !         Clear the error tables.
         !
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
         !
         !         PROCESS A MESSAGE...
         !         SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
         !         OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
         !
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.&
               MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.&
               LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
         !
         IF (NMSG.LT.LENTAB) THEN
            !
            !            Empty slot found for new message.
            !
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
            !
            !            Table is full.
            !
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
 !
 !      Formats.
 !
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /&
         ' LIBRARY    SUBROUTINE MESSAGE START             NERR',&
         '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
      ! DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
      ! ***BEGIN PROLOGUE  D1MACH
      ! ***PURPOSE  Return floating point machine dependent constants.
      ! ***LIBRARY   SLATEC
      ! ***CATEGORY  R1
      ! ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
      ! ***KEYWORDS  MACHINE CONSTANTS
      ! ***AUTHOR  Fox, P. A., (Bell Labs)
      !            Hall, A. D., (Bell Labs)
      !            Schryer, N. L., (Bell Labs)
      ! ***DESCRIPTION
      !
      !    D1MACH can be used to obtain machine-dependent parameters for the
      !    local machine environment.  It is a function subprogram with one
      !    (input) argument, and can be referenced as follows:
      !
      !         D = D1MACH(I)
      !
      !    where I=1,...,5.  The (output) value of D above is determined by
      !    the (input) value of I.  The results for various values of I are
      !    discussed below.
      !
      !    D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
      !    D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
      !    D1MACH( 3) = B**(-T), the smallest relative spacing.
      !    D1MACH( 4) = B**(1-T), the largest relative spacing.
      !    D1MACH( 5) = LOG10(B)
      !
      !    Assume double precision numbers are represented in the T-digit,
      !    base-B form
      !
      !               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
      !
      !    where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
      !    EMIN .LE. E .LE. EMAX.
      !
      !    The values of B, T, EMIN and EMAX are provided in I1MACH as
      !    follows:
      !    I1MACH(10) = B, the base.
      !    I1MACH(14) = T, the number of base-B digits.
      !    I1MACH(15) = EMIN, the smallest exponent E.
      !    I1MACH(16) = EMAX, the largest exponent E.
      !
      !    To alter this function for a particular environment, the desired
      !    set of DATA statements should be activated by removing the C from
      !    column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
      !    checked for consistency with the local operating system.
      !
      ! ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
      !                  a portable library, ACM Transactions on Mathematical
      !                  Software 4, 2 (June 1978), pp. 177-188.
      ! ***ROUTINES CALLED  XERMSG
      ! ***REVISION HISTORY  (YYMMDD)
      !    750101  DATE WRITTEN
      !    890213  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
      !    900618  Added DEC RISC constants.  (WRB)
      !    900723  Added IBM RS 6000 constants.  (WRB)
      !    900911  Added SUN 386i constants.  (WRB)
      !    910710  Added HP 730 constants.  (SMR)
      !    911114  Added Convex IEEE constants.  (WRB)
      !    920121  Added SUN -r8 compiler option constants.  (WRB)
      !    920229  Added Touchstone Delta i860 constants.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      !    920625  Added CONVEX -p8 and -pd8 compiler option constants.
      !            (BKS, WRB)
      !    930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
      !    010817  Elevated IEEE to highest importance; see next set of
      !            comments below.  (DWL)
      ! ***END PROLOGUE  D1MACH
      !
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
      !
      !  Initial data here correspond to the IEEE standard.  The values for
      !  DMACH(1), DMACH(3) and DMACH(4) are slight upper bounds.  The value
      !  for DMACH(2) is a slight lower bound.  The value for DMACH(5) is
      !  a 20-digit approximation.  If one of the sets of initial data below
      !  is preferred, do the necessary commenting and uncommenting. (DWL)
      DOUBLE PRECISION DMACH(5)
      DATA DMACH / 2.23D-308, 1.79D+308, 1.111D-16, 2.222D-16,&
                   0.30102999566398119521D0 /
      SAVE DMACH
      !
      ! c      EQUIVALENCE (DMACH(1),SMALL(1))
      ! c      EQUIVALENCE (DMACH(2),LARGE(1))
      ! c      EQUIVALENCE (DMACH(3),RIGHT(1))
      ! c      EQUIVALENCE (DMACH(4),DIVER(1))
      ! c      EQUIVALENCE (DMACH(5),LOG10(1))
      !
      !      MACHINE CONSTANTS FOR THE AMIGA
      !      ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
      !
      !      DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
      !      DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
      !      DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
      !      DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
      !      DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE AMIGA
      !      ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
      !
      !      DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
      !      DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
      !      DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
      !      DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
      !      DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE APOLLO
      !
      !      DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
      !      DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
      !      DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
      !      DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
      !      DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
      !
      !      DATA SMALL(1) / ZC00800000 /
      !      DATA SMALL(2) / Z000000000 /
      !      DATA LARGE(1) / ZDFFFFFFFF /
      !      DATA LARGE(2) / ZFFFFFFFFF /
      !      DATA RIGHT(1) / ZCC5800000 /
      !      DATA RIGHT(2) / Z000000000 /
      !      DATA DIVER(1) / ZCC6800000 /
      !      DATA DIVER(2) / Z000000000 /
      !      DATA LOG10(1) / ZD00E730E7 /
      !      DATA LOG10(2) / ZC77800DC0 /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
      !
      !      DATA SMALL(1) / O1771000000000000 /
      !      DATA SMALL(2) / O0000000000000000 /
      !      DATA LARGE(1) / O0777777777777777 /
      !      DATA LARGE(2) / O0007777777777777 /
      !      DATA RIGHT(1) / O1461000000000000 /
      !      DATA RIGHT(2) / O0000000000000000 /
      !      DATA DIVER(1) / O1451000000000000 /
      !      DATA DIVER(2) / O0000000000000000 /
      !      DATA LOG10(1) / O1157163034761674 /
      !      DATA LOG10(2) / O0006677466732724 /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
      !
      !      DATA SMALL(1) / O1771000000000000 /
      !      DATA SMALL(2) / O7770000000000000 /
      !      DATA LARGE(1) / O0777777777777777 /
      !      DATA LARGE(2) / O7777777777777777 /
      !      DATA RIGHT(1) / O1461000000000000 /
      !      DATA RIGHT(2) / O0000000000000000 /
      !      DATA DIVER(1) / O1451000000000000 /
      !      DATA DIVER(2) / O0000000000000000 /
      !      DATA LOG10(1) / O1157163034761674 /
      !      DATA LOG10(2) / O0006677466732724 /
      !
      !      MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
      !
      !      DATA SMALL(1) / Z"3001800000000000" /
      !      DATA SMALL(2) / Z"3001000000000000" /
      !      DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
      !      DATA LARGE(2) / Z"4FFE000000000000" /
      !      DATA RIGHT(1) / Z"3FD2800000000000" /
      !      DATA RIGHT(2) / Z"3FD2000000000000" /
      !      DATA DIVER(1) / Z"3FD3800000000000" /
      !      DATA DIVER(2) / Z"3FD3000000000000" /
      !      DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
      !      DATA LOG10(2) / Z"3FFFF7988F8959AC" /
      !
      !      MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
      !
      !      DATA SMALL(1) / 00564000000000000000B /
      !      DATA SMALL(2) / 00000000000000000000B /
      !      DATA LARGE(1) / 37757777777777777777B /
      !      DATA LARGE(2) / 37157777777777777777B /
      !      DATA RIGHT(1) / 15624000000000000000B /
      !      DATA RIGHT(2) / 00000000000000000000B /
      !      DATA DIVER(1) / 15634000000000000000B /
      !      DATA DIVER(2) / 00000000000000000000B /
      !      DATA LOG10(1) / 17164642023241175717B /
      !      DATA LOG10(2) / 16367571421742254654B /
      !
      !      MACHINE CONSTANTS FOR THE CELERITY C1260
      !
      !      DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
      !      DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
      !      DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
      !      DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
      !      DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -fn OR -pd8 COMPILER OPTION
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CC0000000000000' /
      !      DATA DMACH(4) / Z'3CD0000000000000' /
      !      DATA DMACH(5) / Z'3FF34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -fi COMPILER OPTION
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CA0000000000000' /
      !      DATA DMACH(4) / Z'3CB0000000000000' /
      !      DATA DMACH(5) / Z'3FD34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -p8 COMPILER OPTION
      !
      !      DATA DMACH(1) / Z'00010000000000000000000000000000' /
      !      DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3F900000000000000000000000000000' /
      !      DATA DMACH(4) / Z'3F910000000000000000000000000000' /
      !      DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
      !
      !      MACHINE CONSTANTS FOR THE CRAY
      !
      !      DATA SMALL(1) / 201354000000000000000B /
      !      DATA SMALL(2) / 000000000000000000000B /
      !      DATA LARGE(1) / 577767777777777777777B /
      !      DATA LARGE(2) / 000007777777777777774B /
      !      DATA RIGHT(1) / 376434000000000000000B /
      !      DATA RIGHT(2) / 000000000000000000000B /
      !      DATA DIVER(1) / 376444000000000000000B /
      !      DATA DIVER(2) / 000000000000000000000B /
      !      DATA LOG10(1) / 377774642023241175717B /
      !      DATA LOG10(2) / 000007571421742254654B /
      !
      !      MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
      !      NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
      !      STATIC DMACH(5)
      !
      !      DATA SMALL /    20K, 3*0 /
      !      DATA LARGE / 77777K, 3*177777K /
      !      DATA RIGHT / 31420K, 3*0 /
      !      DATA DIVER / 32020K, 3*0 /
      !      DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
      !
      !      MACHINE CONSTANTS FOR THE DEC ALPHA
      !      USING G_FLOAT
      !
      !      DATA DMACH(1) / '0000000000000010'X /
      !      DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
      !      DATA DMACH(3) / '0000000000003CC0'X /
      !      DATA DMACH(4) / '0000000000003CD0'X /
      !      DATA DMACH(5) / '79FF509F44133FF3'X /
      !
      !      MACHINE CONSTANTS FOR THE DEC ALPHA
      !      USING IEEE_FORMAT
      !
      !      DATA DMACH(1) / '0010000000000000'X /
      !      DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
      !      DATA DMACH(3) / '3CA0000000000000'X /
      !      DATA DMACH(4) / '3CB0000000000000'X /
      !      DATA DMACH(5) / '3FD34413509F79FF'X /
      !
      !      MACHINE CONSTANTS FOR THE DEC RISC
      !
      !      DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
      !      DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
      !      DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
      !      DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
      !      DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
      !
      !      MACHINE CONSTANTS FOR THE DEC VAX
      !      USING D_FLOATING
      !      (EXPRESSED IN INTEGER AND HEXADECIMAL)
      !      THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
      !      THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
      !
      !      DATA SMALL(1), SMALL(2) /        128,           0 /
      !      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
      !      DATA RIGHT(1), RIGHT(2) /       9344,           0 /
      !      DATA DIVER(1), DIVER(2) /       9472,           0 /
      !      DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
      !
      !      DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
      !      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
      !      DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
      !      DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
      !      DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
      !
      !      MACHINE CONSTANTS FOR THE DEC VAX
      !      USING G_FLOATING
      !      (EXPRESSED IN INTEGER AND HEXADECIMAL)
      !      THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
      !      THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
      !
      !      DATA SMALL(1), SMALL(2) /         16,           0 /
      !      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
      !      DATA RIGHT(1), RIGHT(2) /      15552,           0 /
      !      DATA DIVER(1), DIVER(2) /      15568,           0 /
      !      DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
      !
      !      DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
      !      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
      !      DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
      !      DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
      !      DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
      !
      !      MACHINE CONSTANTS FOR THE ELXSI 6400
      !      (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
      !
      !      DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
      !      DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
      !      DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
      !      DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
      !      DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
      !
      !      MACHINE CONSTANTS FOR THE HARRIS 220
      !
      !      DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
      !      DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
      !      DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
      !      DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
      !      DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
      !
      !      MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
      !
      !      DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
      !      DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
      !      DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
      !      DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
      !      DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
      !
      !      MACHINE CONSTANTS FOR THE HP 730
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CA0000000000000' /
      !      DATA DMACH(4) / Z'3CB0000000000000' /
      !      DATA DMACH(5) / Z'3FD34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE HP 2100
      !      THREE WORD DOUBLE PRECISION OPTION WITH FTN4
      !
      !      DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
      !      DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
      !      DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
      !      DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
      !      DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
      !
      !      MACHINE CONSTANTS FOR THE HP 2100
      !      FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
      !
      !      DATA SMALL(1), SMALL(2) /  40000B,       0 /
      !      DATA SMALL(3), SMALL(4) /       0,       1 /
      !      DATA LARGE(1), LARGE(2) /  77777B, 177777B /
      !      DATA LARGE(3), LARGE(4) / 177777B, 177776B /
      !      DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
      !      DATA RIGHT(3), RIGHT(4) /       0,    225B /
      !      DATA DIVER(1), DIVER(2) /  40000B,       0 /
      !      DATA DIVER(3), DIVER(4) /       0,    227B /
      !      DATA LOG10(1), LOG10(2) /  46420B,  46502B /
      !      DATA LOG10(3), LOG10(4) /  76747B, 176377B /
      !
      !      MACHINE CONSTANTS FOR THE HP 9000
      !
      !      DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
      !      DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
      !      DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
      !      DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
      !      DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
      !
      !      MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
      !      THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
      !      THE PERKIN ELMER (INTERDATA) 7/32.
      !
      !      DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
      !      DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
      !      DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
      !      DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
      !      DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
      !
      !      MACHINE CONSTANTS FOR THE IBM PC
      !      ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
      !      ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
      !
      !      DATA SMALL(1) / 2.23D-308  /
      !      DATA LARGE(1) / 1.79D+308  /
      !      DATA RIGHT(1) / 1.11D-16   /
      !      DATA DIVER(1) / 2.22D-16   /
      !      DATA LOG10(1) / 0.301029995663981195D0 /
      !
      !      MACHINE CONSTANTS FOR THE IBM RS 6000
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CA0000000000000' /
      !      DATA DMACH(4) / Z'3CB0000000000000' /
      !      DATA DMACH(5) / Z'3FD34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE INTEL i860
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CA0000000000000' /
      !      DATA DMACH(4) / Z'3CB0000000000000' /
      !      DATA DMACH(5) / Z'3FD34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
      !
      !      DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
      !      DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
      !      DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
      !      DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
      !      DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
      !
      !      MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
      !
      !      DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
      !      DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
      !      DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
      !      DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
      !      DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
      !
      !      MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
      !      32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
      !
      !      DATA SMALL(1), SMALL(2) /    8388608,           0 /
      !      DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
      !      DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
      !      DATA DIVER(1), DIVER(2) /  620756992,           0 /
      !      DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
      !
      !      DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
      !      DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
      !      DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
      !      DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
      !      DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
      !
      !      MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
      !      16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
      !
      !      DATA SMALL(1), SMALL(2) /    128,      0 /
      !      DATA SMALL(3), SMALL(4) /      0,      0 /
      !      DATA LARGE(1), LARGE(2) /  32767,     -1 /
      !      DATA LARGE(3), LARGE(4) /     -1,     -1 /
      !      DATA RIGHT(1), RIGHT(2) /   9344,      0 /
      !      DATA RIGHT(3), RIGHT(4) /      0,      0 /
      !      DATA DIVER(1), DIVER(2) /   9472,      0 /
      !      DATA DIVER(3), DIVER(4) /      0,      0 /
      !      DATA LOG10(1), LOG10(2) /  16282,   8346 /
      !      DATA LOG10(3), LOG10(4) / -31493, -12296 /
      !
      !      DATA SMALL(1), SMALL(2) / O000200, O000000 /
      !      DATA SMALL(3), SMALL(4) / O000000, O000000 /
      !      DATA LARGE(1), LARGE(2) / O077777, O177777 /
      !      DATA LARGE(3), LARGE(4) / O177777, O177777 /
      !      DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
      !      DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
      !      DATA DIVER(1), DIVER(2) / O022400, O000000 /
      !      DATA DIVER(3), DIVER(4) / O000000, O000000 /
      !      DATA LOG10(1), LOG10(2) / O037632, O020232 /
      !      DATA LOG10(3), LOG10(4) / O102373, O147770 /
      !
      !      MACHINE CONSTANTS FOR THE SILICON GRAPHICS
      !
      !      DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
      !      DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
      !      DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
      !      DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
      !      DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE SUN
      !
      !      DATA DMACH(1) / Z'0010000000000000' /
      !      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3CA0000000000000' /
      !      DATA DMACH(4) / Z'3CB0000000000000' /
      !      DATA DMACH(5) / Z'3FD34413509F79FF' /
      !
      !      MACHINE CONSTANTS FOR THE SUN
      !      USING THE -r8 COMPILER OPTION
      !
      !      DATA DMACH(1) / Z'00010000000000000000000000000000' /
      !      DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
      !      DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
      !      DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
      !      DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
      !
      !      MACHINE CONSTANTS FOR THE SUN 386i
      !
      !      DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
      !      DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
      !      DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
      !      DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
      !      DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
      !
      !      MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
      !
      !      DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
      !      DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
      !      DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
      !      DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
      !      DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
      !
      ! ***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',&
         'I OUT OF BOUNDS', 1, 2)
      !
      D1MACH = DMACH(I)
      RETURN
      !
      END
      ! DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
      ! ***BEGIN PROLOGUE  XGETUA
      ! ***PURPOSE  Return unit number(s) to which error messages are being
      !             sent.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3C
      ! ***TYPE      ALL (XGETUA-A)
      ! ***KEYWORDS  ERROR, XERROR
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !      Abstract
      !         XGETUA may be called to determine the unit number or numbers
      !         to which error messages are being sent.
      !         These unit numbers may have been set by a call to XSETUN,
      !         or a call to XSETUA, or may be a default value.
      !
      !      Description of Parameters
      !       --Output--
      !         IUNIT - an array of one to five unit numbers, depending
      !                 on the value of N.  A value of zero refers to the
      !                 default unit, as defined by the I1MACH machine
      !                 constant routine.  Only IUNIT(1),...,IUNIT(N) are
      !                 defined by XGETUA.  The values of IUNIT(N+1),...,
      !                 IUNIT(5) are not defined (for N .LT. 5) or altered
      !                 in any way by XGETUA.
      !         N     - the number of units to which copies of the
      !                 error messages are being sent.  N will be in the
      !                 range from 1 to 5.
      !
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  J4SAVE
      ! ***REVISION HISTORY  (YYMMDD)
      !    790801  DATE WRITTEN
      !    861211  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
      ! ***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      ! DECK DSTEPS
      SUBROUTINE DSTEPS (DF, NEQN, Y, X, H, EPS, WT, START, HOLD, K,&
         KOLD, CRASH, PHI, P, YP, PSI, ALPHA, BETA, SIG, V, W, G,&
         PHASE1, NS, NORND, KSTEPS, TWOU, FOURU, XOLD, KPREV, IVC, IV,&
         KGI, GI, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DSTEPS
      ! ***PURPOSE  Integrate a system of first order ordinary differential
      !             equations one step.
      ! ***LIBRARY   SLATEC (DEPAC)
      ! ***CATEGORY  I1A1B
      ! ***TYPE      DOUBLE PRECISION (STEPS-S, DSTEPS-D)
      ! ***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
      !              ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
      ! ***AUTHOR  Shampine, L. F., (SNLA)
      !            Gordon, M. K., (SNLA)
      !              MODIFIED BY H.A. WATTS
      ! ***DESCRIPTION
      !
      !    Written by L. F. Shampine and M. K. Gordon
      !
      !    Abstract
      !
      !    Subroutine  DSTEPS  is normally used indirectly through subroutine
      !    DDEABM .  Because  DDEABM  suffices for most problems and is much
      !    easier to use, using it should be considered before using  DSTEPS
      !    alone.
      !
      !    Subroutine DSTEPS integrates a system of  NEQN  first order ordinary
      !    differential equations one step, normally from X to X+H, using a
      !    modified divided difference form of the Adams Pece formulas.  Local
      !    extrapolation is used to improve absolute stability and accuracy.
      !    The code adjusts its order and step size to control the local error
      !    per unit step in a generalized sense.  Special devices are included
      !    to control roundoff error and to detect when the user is requesting
      !    too much accuracy.
      !
      !    This code is completely explained and documented in the text,
      !    Computer Solution of Ordinary Differential Equations, The Initial
      !    Value Problem  by L. F. Shampine and M. K. Gordon.
      !    Further details on use of this code are available in "Solving
      !    Ordinary Differential Equations with ODE, STEP, and INTRP",
      !    by L. F. Shampine and M. K. Gordon, SLA-73-1060.
      !
      !
      !    The parameters represent --
      !       DF -- subroutine to evaluate derivatives
      !       NEQN -- number of equations to be integrated
      !       Y(*) -- solution vector at X
      !       X -- independent variable
      !       H -- appropriate step size for next step.  Normally determined by
      !            code
      !       EPS -- local error tolerance
      !       WT(*) -- vector of weights for error criterion
      !       START -- logical variable set .TRUE. for first step,  .FALSE.
      !            otherwise
      !       HOLD -- step size used for last successful step
      !       K -- appropriate order for next step (determined by code)
      !       KOLD -- order used for last successful step
      !       CRASH -- logical variable set .TRUE. when no step can be taken,
      !            .FALSE. otherwise.
      !       YP(*) -- derivative of solution vector at  X  after successful
      !            step
      !       KSTEPS -- counter on attempted steps
      !       TWOU -- 2.*U where U is machine unit roundoff quantity
      !       FOURU -- 4.*U where U is machine unit roundoff quantity
      !       RPAR,IPAR -- parameter arrays which you may choose to use
      !             for communication between your program and subroutine F.
      !             They are not altered or used by DSTEPS.
      !    The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G,
      !    W,P,IV and GI are required for the interpolation subroutine SINTRP.
      !    The remaining variables and arrays are included in the call list
      !    only to eliminate local retention of variables between calls.
      !
      !    Input to DSTEPS
      !
      !       First call --
      !
      !    The user must provide storage in his calling program for all arrays
      !    in the call list, namely
      !
      !      DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
      !     1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
      !     2  RPAR(*),IPAR(*)
      !
      !     **Note**
      !
      !    The user must also declare  START ,  CRASH ,  PHASE1  and  NORND
      !    logical variables and  DF  an EXTERNAL subroutine, supply the
      !    subroutine  DF(X,Y,YP)  to evaluate
      !       DY(I)/DX = YP(I) = DF(X,Y(1),Y(2),...,Y(NEQN))
      !    and initialize only the following parameters.
      !       NEQN -- number of equations to be integrated
      !       Y(*) -- vector of initial values of dependent variables
      !       X -- initial value of the independent variable
      !       H -- nominal step size indicating direction of integration
      !            and maximum size of step.  Must be variable
      !       EPS -- local error tolerance per step.  Must be variable
      !       WT(*) -- vector of non-zero weights for error criterion
      !       START -- .TRUE.
      !       YP(*) -- vector of initial derivative values
      !       KSTEPS -- set KSTEPS to zero
      !       TWOU -- 2.*U where U is machine unit roundoff quantity
      !       FOURU -- 4.*U where U is machine unit roundoff quantity
      !    Define U to be the machine unit roundoff quantity by calling
      !    the function routine  D1MACH,  U = D1MACH(4), or by
      !    computing U so that U is the smallest positive number such
      !    that 1.0+U .GT. 1.0.
      !
      !    DSTEPS  requires that the L2 norm of the vector with components
      !    LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The
      !    array  WT  allows the user to specify an error test appropriate
      !    for his problem.  For example,
      !       WT(L) = 1.0  specifies absolute error,
      !             = ABS(Y(L))  error relative to the most recent value of the
      !                  L-th component of the solution,
      !             = ABS(YP(L))  error relative to the most recent value of
      !                  the L-th component of the derivative,
      !             = MAX(WT(L),ABS(Y(L)))  error relative to the largest
      !                  magnitude of L-th component obtained so far,
      !             = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed
      !                  relative-absolute test where  RELERR  is relative
      !                  error,  ABSERR  is absolute error and  EPS =
      !                  MAX(RELERR,ABSERR) .
      !
      !       Subsequent calls --
      !
      !    Subroutine  DSTEPS  is designed so that all information needed to
      !    continue the integration, including the step size  H  and the order
      !    K , is returned with each step.  With the exception of the step
      !    size, the error tolerance, and the weights, none of the parameters
      !    should be altered.  The array  WT  must be updated after each step
      !    to maintain relative error tests like those above.  Normally the
      !    integration is continued just beyond the desired endpoint and the
      !    solution interpolated there with subroutine  SINTRP .  If it is
      !    impossible to integrate beyond the endpoint, the step size may be
      !    reduced to hit the endpoint since the code will not take a step
      !    larger than the  H  input.  Changing the direction of integration,
      !    i.e., the sign of  H , requires the user set  START = .TRUE. before
      !    calling  DSTEPS  again.  This is the only situation in which  START
      !    should be altered.
      !
      !    Output from DSTEPS
      !
      !       Successful Step --
      !
      !    The subroutine returns after each successful step with  START  and
      !    CRASH  set .FALSE. .  X  represents the independent variable
      !    advanced one step of length  HOLD  from its value on input and  Y
      !    the solution vector at the new value of  X .  All other parameters
      !    represent information corresponding to the new  X  needed to
      !    continue the integration.
      !
      !       Unsuccessful Step --
      !
      !    When the error tolerance is too small for the machine precision,
      !    the subroutine returns without taking a step and  CRASH = .TRUE. .
      !    An appropriate step size and error tolerance for continuing are
      !    estimated and all other information is restored as upon input
      !    before returning.  To continue with the larger tolerance, the user
      !    just calls the code again.  A restart is neither required nor
      !    desirable.
      !
      ! ***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary
      !                  differential equations with ODE, STEP, and INTRP,
      !                  Report SLA-73-1060, Sandia Laboratories, 1973.
      ! ***ROUTINES CALLED  D1MACH, DHSTRT
      ! ***REVISION HISTORY  (YYMMDD)
      !    740101  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DSTEPS
      !
      INTEGER I, IFAIL, IM1, IP1, IPAR, IQ, J, K, KM1, KM2, KNEW,&
            KOLD, KP1, KP2, KSTEPS, L, LIMIT1, LIMIT2, NEQN, NS, NSM2,&
            NSP1, NSP2
      DOUBLE PRECISION ABSH, ALPHA, BETA, BIG, D1MACH,&
            EPS, ERK, ERKM1, ERKM2, ERKP1, ERR,&
            FOURU, G, GI, GSTR, H, HNEW, HOLD, P, P5EPS, PHI, PSI, R,&
            REALI, REALNS, RHO, ROUND, RPAR, SIG, TAU, TEMP1,&
            TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TWO, TWOU, U, V, W, WT,&
            X, XOLD, Y, YP
      LOGICAL START,CRASH,PHASE1,NORND
      DIMENSION Y(*),WT(*),PHI(NEQN,16),P(*),YP(*),PSI(12),&
        ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),&
        RPAR(*),IPAR(*)
      DIMENSION TWO(13),GSTR(13)
      EXTERNAL DF
      SAVE TWO, GSTR
      !
      DATA TWO(1),TWO(2),TWO(3),TWO(4),TWO(5),TWO(6),TWO(7),TWO(8),&
           TWO(9),TWO(10),TWO(11),TWO(12),TWO(13)&
           /2.0D0,4.0D0,8.0D0,16.0D0,32.0D0,64.0D0,128.0D0,256.0D0,&
            512.0D0,1024.0D0,2048.0D0,4096.0D0,8192.0D0/
      DATA GSTR(1),GSTR(2),GSTR(3),GSTR(4),GSTR(5),GSTR(6),GSTR(7),&
           GSTR(8),GSTR(9),GSTR(10),GSTR(11),GSTR(12),GSTR(13)&
           /0.5D0,0.0833D0,0.0417D0,0.0264D0,0.0188D0,0.0143D0,0.0114D0,&
            0.00936D0,0.00789D0,0.00679D0,0.00592D0,0.00524D0,0.00468D0/
      !
      !        ***     BEGIN BLOCK 0     ***
      !    CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
      !    PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
      !    STARTING STEP SIZE.
      !                    ***
      !
      !    IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
      !
      ! ***FIRST EXECUTABLE STATEMENT  DSTEPS
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5D0*EPS
      !
      !    IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
      !
      ROUND = 0.0D0
      DO 10 L = 1,NEQN
 10     ROUND = ROUND + (Y(L)/WT(L))**2
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15
      EPS = 2.0D0*ROUND*(1.0D0 + FOURU)
      RETURN
 15   CRASH = .FALSE.
      G(1) = 1.0D0
      G(2) = 0.5D0
      SIG(1) = 1.0D0
      IF(.NOT.START) GO TO 99
      !
      !    INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
      !
      !      CALL DF(X,Y,YP,RPAR,IPAR)
      !      SUM = 0.0
      DO 20 L = 1,NEQN
        PHI(L,1) = YP(L)
   20   PHI(L,2) = 0.0D0
      ! 20     SUM = SUM + (YP(L)/WT(L))**2
      !      SUM = SQRT(SUM)
      !      ABSH = ABS(H)
      !      IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM)
      !      H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
      !
      U = D1MACH(4)
      BIG = SQRT(D1MACH(2))
      CALL DHSTRT(DF,NEQN,X,X+H,Y,YP,WT,1,U,BIG,&
                   PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H)
      !
      HOLD = 0.0D0
      K = 1
      KOLD = 0
      KPREV = 0
      START = .FALSE.
      PHASE1 = .TRUE.
      NORND = .TRUE.
      IF(P5EPS .GT. 100.0D0*ROUND) GO TO 99
      NORND = .FALSE.
      DO 25 L = 1,NEQN
 25     PHI(L,15) = 0.0D0
 99   IFAIL = 0
 !        ***     END BLOCK 0     ***
 !
 !        ***     BEGIN BLOCK 1     ***
 !    COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
 !    THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
 !                    ***
 !
 100  KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      KM2 = K-2
      !
      !    NS IS THE NUMBER OF DSTEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
      !    ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
      !
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1
      NSP1 = NS+1
      IF (K .LT. NS) GO TO 199
      !
      !    COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
      !    ARE CHANGED
      !
      BETA(NS) = 1.0D0
      REALNS = NS
      ALPHA(NS) = 1.0D0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0D0
      IF(K .LT. NSP1) GO TO 110
      DO 105 I = NSP1,K
        IM1 = I-1
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H
        ALPHA(I) = H/TEMP1
        REALI = I
 105    SIG(I+1) = REALI*ALPHA(I)*SIG(I)
 110  PSI(K) = TEMP1
      !
      !    COMPUTE COEFFICIENTS G(*)
      !
      !    INITIALIZE V(*) AND SET W(*).
      !
      IF(NS .GT. 1) GO TO 120
      DO 115 IQ = 1,K
        TEMP3 = IQ*(IQ+1)
        V(IQ) = 1.0D0/TEMP3
 115    W(IQ) = V(IQ)
      IVC = 0
      KGI = 0
      IF (K .EQ. 1) GO TO 140
      KGI = 1
      GI(1) = W(2)
      GO TO 140
 !
 !    IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
 !
 120  IF(K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122
      JV = KP1 - IV(IVC)
      IVC = IVC - 1
      GO TO 123
 122  JV = 1
      TEMP4 = K*KP1
      V(K) = 1.0D0/TEMP4
      W(K) = V(K)
      IF (K .NE. 2) GO TO 123
      KGI = 1
      GI(1) = W(2)
 123  NSM2 = NS-2
      IF(NSM2 .LT. JV) GO TO 130
      DO 125 J = JV,NSM2
        I = K-J
        V(I) = V(I) - ALPHA(J+1)*V(I+1)
 125    W(I) = V(I)
      IF (I .NE. 2) GO TO 130
      KGI = NS - 1
      GI(KGI) = W(2)
 !
 !    UPDATE V(*) AND SET W(*)
 !
 130  LIMIT1 = KP1 - NS
      TEMP5 = ALPHA(NS)
      DO 135 IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1)
 135    W(IQ) = V(IQ)
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1)
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1
      IV(IVC) = LIMIT1 + 2
 !
 !    COMPUTE THE G(*) IN THE WORK VECTOR W(*)
 !
 140  NSP2 = NS + 2
      KPREV = K
      IF(KP1 .LT. NSP2) GO TO 199
      DO 150 I = NSP2,KP1
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO 145 IQ = 1,LIMIT2
 145      W(IQ) = W(IQ) - TEMP6*W(IQ+1)
 150    G(I) = W(1)
 199    CONTINUE
      !        ***     END BLOCK 1     ***
      !
      !        ***     BEGIN BLOCK 2     ***
      !    PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
      !    SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
      !    K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
      !                    ***
      !
      !    INCREMENT COUNTER ON ATTEMPTED DSTEPS
      !
      KSTEPS = KSTEPS + 1
      !
      !    CHANGE PHI TO PHI STAR
      !
      IF(K .LT. NSP1) GO TO 215
      DO 210 I = NSP1,K
        TEMP1 = BETA(I)
        DO 205 L = 1,NEQN
 205      PHI(L,I) = TEMP1*PHI(L,I)
 210    CONTINUE
 !
 !    PREDICT SOLUTION AND DIFFERENCES
 !
 215  DO 220 L = 1,NEQN
        PHI(L,KP2) = PHI(L,KP1)
        PHI(L,KP1) = 0.0D0
 220    P(L) = 0.0D0
      DO 230 J = 1,K
        I = KP1 - J
        IP1 = I+1
        TEMP2 = G(I)
        DO 225 L = 1,NEQN
          P(L) = P(L) + TEMP2*PHI(L,I)
 225      PHI(L,I) = PHI(L,I) + PHI(L,IP1)
 230    CONTINUE
      IF(NORND) GO TO 240
      DO 235 L = 1,NEQN
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU
 235    PHI(L,16) = (P(L) - Y(L)) - TAU
      GO TO 250
 240  DO 245 L = 1,NEQN
 245    P(L) = Y(L) + H*P(L)
 250  XOLD = X
      X = X + H
      ABSH = ABS(H)
      CALL DF(X,P,YP,RPAR,IPAR)
      !
      !    ESTIMATE ERRORS AT ORDERS K,K-1,K-2
      !
      ERKM2 = 0.0D0
      ERKM1 = 0.0D0
      ERK = 0.0D0
      DO 265 L = 1,NEQN
        TEMP3 = 1.0D0/WT(L)
        TEMP4 = YP(L) - PHI(L,1)
        IF(KM2)265,260,255
 255    ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
 260    ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
 265    ERK = ERK + (TEMP4*TEMP3)**2
      IF(KM2)280,275,270
 270  ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2)
 275  ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1)
 280  TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1))
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
      !
      !    TEST IF ORDER SHOULD BE LOWERED
      !
      IF(KM2)299,290,285
 285  IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      GO TO 299
 290  IF(ERKM1 .LE. 0.5D0*ERK) KNEW = KM1
 !
 !    TEST IF STEP SUCCESSFUL
 !
 299  IF(ERR .LE. EPS) GO TO 400
      !        ***     END BLOCK 2     ***
      !
      !        ***     BEGIN BLOCK 3     ***
      !    THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
      !    IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE
      !    THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
      !    TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
      !    PRECISION.
      !                    ***
      !
      !    RESTORE X, PHI(*,*) AND PSI(*)
      !
      PHASE1 = .FALSE.
      X = XOLD
      DO 310 I = 1,K
        TEMP1 = 1.0D0/BETA(I)
        IP1 = I+1
        DO 305 L = 1,NEQN
 305      PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
 310    CONTINUE
      IF(K .LT. 2) GO TO 320
      DO 315 I = 2,K
 315    PSI(I-1) = PSI(I) - H
 !
 !    ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
 !    SIZE
 !
 320  IFAIL = IFAIL + 1
      TEMP2 = 0.5D0
      IF(IFAIL - 3) 335,330,325
 325  IF(P5EPS .LT. 0.25D0*ERK) TEMP2 = SQRT(P5EPS/ERK)
 330  KNEW = 1
 335  H = TEMP2*H
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS
      RETURN
 340  GO TO 100
 !        ***     END BLOCK 3     ***
 !
 !        ***     BEGIN BLOCK 4     ***
 !    THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
 !    THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
 !    DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
 !                    ***
 400  KOLD = K
      HOLD = H
      !
      !    CORRECT AND EVALUATE
      !
      TEMP1 = H*G(KP1)
      IF(NORND) GO TO 410
      DO 405 L = 1,NEQN
        TEMP3 = Y(L)
        RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
        Y(L) = P(L) + RHO
        PHI(L,15) = (Y(L) - P(L)) - RHO
 405    P(L) = TEMP3
      GO TO 420
 410  DO 415 L = 1,NEQN
        TEMP3 = Y(L)
        Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
 415    P(L) = TEMP3
 420  CALL DF(X,Y,YP,RPAR,IPAR)
      !
      !    UPDATE DIFFERENCES FOR NEXT STEP
      !
      DO 425 L = 1,NEQN
        PHI(L,KP1) = YP(L) - PHI(L,1)
 425    PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      DO 435 I = 1,K
        DO 430 L = 1,NEQN
 430      PHI(L,I) = PHI(L,I) + PHI(L,KP1)
 435    CONTINUE
      !
      !    ESTIMATE ERROR AT ORDER K+1 UNLESS:
      !      IN FIRST PHASE WHEN ALWAYS RAISE ORDER,
      !      ALREADY DECIDED TO LOWER ORDER,
      !      STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
      !
      ERKP1 = 0.0D0
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE.
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455
      IF(KP1 .GT. NS) GO TO 460
      DO 440 L = 1,NEQN
 440    ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
      !
      !    USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER
      !    FOR NEXT STEP
      !
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. 0.5D0*ERK) GO TO 460
      GO TO 450
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460
 !
 !    HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE
 !    BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
 !
 !    RAISE ORDER
 !
 450  K = KP1
      ERK = ERKP1
      GO TO 460
 !
 !    LOWER ORDER
 !
 455  K = KM1
      ERK = ERKM1
 !
 !    WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
 !
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = K+1
      R = (P5EPS/ERK)**(1.0D0/TEMP2)
      HNEW = ABSH*MAX(0.5D0,MIN(0.9D0,R))
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
 465  H = HNEW
      RETURN
      !        ***     END BLOCK 4     ***
      END
      ! DECK FDUMP
      SUBROUTINE FDUMP
      ! ***BEGIN PROLOGUE  FDUMP
      ! ***PURPOSE  Symbolic dump (should be locally written).
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3
      ! ***TYPE      ALL (FDUMP-A)
      ! ***KEYWORDS  ERROR, XERMSG
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !         ***Note*** Machine Dependent Routine
      !         FDUMP is intended to be replaced by a locally written
      !         version which produces a symbolic dump.  Failing this,
      !         it should be replaced by a version which prints the
      !         subprogram nesting list.  Note that this dump must be
      !         printed on each of up to five files, as indicated by the
      !         XGETUA routine.  See XSETUA and XGETUA for details.
      !
      !      Written by Ron Jones, with SLATEC Common Math Library Subcommittee
      !
      ! ***REFERENCES  (NONE)
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    790801  DATE WRITTEN
      !    861211  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      ! ***END PROLOGUE  FDUMP
      ! ***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      ! DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
      ! ***BEGIN PROLOGUE  I1MACH
      ! ***PURPOSE  Return integer machine dependent constants.
      ! ***LIBRARY   SLATEC
      ! ***CATEGORY  R1
      ! ***TYPE      INTEGER (I1MACH-I)
      ! ***KEYWORDS  MACHINE CONSTANTS
      ! ***AUTHOR  Fox, P. A., (Bell Labs)
      !            Hall, A. D., (Bell Labs)
      !            Schryer, N. L., (Bell Labs)
      ! ***DESCRIPTION
      !
      !    I1MACH can be used to obtain machine-dependent parameters for the
      !    local machine environment.  It is a function subprogram with one
      !    (input) argument and can be referenced as follows:
      !
      !         K = I1MACH(I)
      !
      !    where I=1,...,16.  The (output) value of K above is determined by
      !    the (input) value of I.  The results for various values of I are
      !    discussed below.
      !
      !    I/O unit numbers:
      !      I1MACH( 1) = the standard input unit.
      !      I1MACH( 2) = the standard output unit.
      !      I1MACH( 3) = the standard punch unit.
      !      I1MACH( 4) = the standard error message unit.
      !
      !    Words:
      !      I1MACH( 5) = the number of bits per integer storage unit.
      !      I1MACH( 6) = the number of characters per integer storage unit.
      !
      !    Integers:
      !      assume integers are represented in the S-digit, base-A form
      !
      !                 sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
      !
      !                 where 0 .LE. X(I) .LT. A for I=0,...,S-1.
      !      I1MACH( 7) = A, the base.
      !      I1MACH( 8) = S, the number of base-A digits.
      !      I1MACH( 9) = A**S - 1, the largest magnitude.
      !
      !    Floating-Point Numbers:
      !      Assume floating-point numbers are represented in the T-digit,
      !      base-B form
      !                 sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
      !
      !                 where 0 .LE. X(I) .LT. B for I=1,...,T,
      !                 0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
      !      I1MACH(10) = B, the base.
      !
      !    Single-Precision:
      !      I1MACH(11) = T, the number of base-B digits.
      !      I1MACH(12) = EMIN, the smallest exponent E.
      !      I1MACH(13) = EMAX, the largest exponent E.
      !
      !    Double-Precision:
      !      I1MACH(14) = T, the number of base-B digits.
      !      I1MACH(15) = EMIN, the smallest exponent E.
      !      I1MACH(16) = EMAX, the largest exponent E.
      !
      !    To alter this function for a particular environment, the desired
      !    set of DATA statements should be activated by removing the C from
      !    column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
      !    checked for consistency with the local operating system.
      !
      ! ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
      !                  a portable library, ACM Transactions on Mathematical
      !                  Software 4, 2 (June 1978), pp. 177-188.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    750101  DATE WRITTEN
      !    891012  Added VAX G-floating constants.  (WRB)
      !    891012  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900618  Added DEC RISC constants.  (WRB)
      !    900723  Added IBM RS 6000 constants.  (WRB)
      !    901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
      !            (RWC)
      !    910710  Added HP 730 constants.  (SMR)
      !    911114  Added Convex IEEE constants.  (WRB)
      !    920121  Added SUN -r8 compiler option constants.  (WRB)
      !    920229  Added Touchstone Delta i860 constants.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      !    920625  Added Convex -p8 and -pd8 compiler option constants.
      !            (BKS, WRB)
      !    930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
      !    930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
      !            options.  (DWL, RWC and WRB).
      !    010817  Elevated IEEE to highest importance; see next set of
      !            comments below.  (DWL)
      ! ***END PROLOGUE  I1MACH
      !
      !  Initial data here correspond to the IEEE standard.  If one of the
      !  sets of initial data below is preferred, do the necessary commenting
      !  and uncommenting. (DWL)
      INTEGER IMACH(16),OUTPUT
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -126 /
      DATA IMACH(13) /        127 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1022 /
      DATA IMACH(16) /       1023 /
      SAVE IMACH
      ! c      EQUIVALENCE (IMACH(4),OUTPUT)
      !
      !      MACHINE CONSTANTS FOR THE AMIGA
      !      ABSOFT COMPILER
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -126 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1022 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE APOLLO
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        129 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1025 /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
      !
      !      DATA IMACH( 1) /          7 /
      !      DATA IMACH( 2) /          2 /
      !      DATA IMACH( 3) /          2 /
      !      DATA IMACH( 4) /          2 /
      !      DATA IMACH( 5) /         36 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         33 /
      !      DATA IMACH( 9) / Z1FFFFFFFF /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -256 /
      !      DATA IMACH(13) /        255 /
      !      DATA IMACH(14) /         60 /
      !      DATA IMACH(15) /       -256 /
      !      DATA IMACH(16) /        255 /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         48 /
      !      DATA IMACH( 6) /          6 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         39 /
      !      DATA IMACH( 9) / O0007777777777777 /
      !      DATA IMACH(10) /          8 /
      !      DATA IMACH(11) /         13 /
      !      DATA IMACH(12) /        -50 /
      !      DATA IMACH(13) /         76 /
      !      DATA IMACH(14) /         26 /
      !      DATA IMACH(15) /        -50 /
      !      DATA IMACH(16) /         76 /
      !
      !      MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         48 /
      !      DATA IMACH( 6) /          6 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         39 /
      !      DATA IMACH( 9) / O0007777777777777 /
      !      DATA IMACH(10) /          8 /
      !      DATA IMACH(11) /         13 /
      !      DATA IMACH(12) /        -50 /
      !      DATA IMACH(13) /         76 /
      !      DATA IMACH(14) /         26 /
      !      DATA IMACH(15) /     -32754 /
      !      DATA IMACH(16) /      32780 /
      !
      !      MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         64 /
      !      DATA IMACH( 6) /          8 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         63 /
      !      DATA IMACH( 9) / 9223372036854775807 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         47 /
      !      DATA IMACH(12) /      -4095 /
      !      DATA IMACH(13) /       4094 /
      !      DATA IMACH(14) /         94 /
      !      DATA IMACH(15) /      -4095 /
      !      DATA IMACH(16) /       4094 /
      !
      !      MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /    6LOUTPUT/
      !      DATA IMACH( 5) /         60 /
      !      DATA IMACH( 6) /         10 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         48 /
      !      DATA IMACH( 9) / 00007777777777777777B /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         47 /
      !      DATA IMACH(12) /       -929 /
      !      DATA IMACH(13) /       1070 /
      !      DATA IMACH(14) /         94 /
      !      DATA IMACH(15) /       -929 /
      !      DATA IMACH(16) /       1069 /
      !
      !      MACHINE CONSTANTS FOR THE CELERITY C1260
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          0 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / Z'7FFFFFFF' /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -126 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1022 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -fn COMPILER OPTION
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1023 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -fi COMPILER OPTION
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -p8 COMPILER OPTION
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         64 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         63 /
      !      DATA IMACH( 9) / 9223372036854775807 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         53 /
      !      DATA IMACH(12) /      -1023 /
      !      DATA IMACH(13) /       1023 /
      !      DATA IMACH(14) /        113 /
      !      DATA IMACH(15) /     -16383 /
      !      DATA IMACH(16) /      16383 /
      !
      !      MACHINE CONSTANTS FOR THE CONVEX
      !      USING THE -pd8 COMPILER OPTION
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         64 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         63 /
      !      DATA IMACH( 9) / 9223372036854775807 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         53 /
      !      DATA IMACH(12) /      -1023 /
      !      DATA IMACH(13) /       1023 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1023 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE CRAY
      !      USING THE 46 BIT INTEGER COMPILER OPTION
      !
      !      DATA IMACH( 1) /        100 /
      !      DATA IMACH( 2) /        101 /
      !      DATA IMACH( 3) /        102 /
      !      DATA IMACH( 4) /        101 /
      !      DATA IMACH( 5) /         64 /
      !      DATA IMACH( 6) /          8 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         46 /
      !      DATA IMACH( 9) / 1777777777777777B /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         47 /
      !      DATA IMACH(12) /      -8189 /
      !      DATA IMACH(13) /       8190 /
      !      DATA IMACH(14) /         94 /
      !      DATA IMACH(15) /      -8099 /
      !      DATA IMACH(16) /       8190 /
      !
      !      MACHINE CONSTANTS FOR THE CRAY
      !      USING THE 64 BIT INTEGER COMPILER OPTION
      !
      !      DATA IMACH( 1) /        100 /
      !      DATA IMACH( 2) /        101 /
      !      DATA IMACH( 3) /        102 /
      !      DATA IMACH( 4) /        101 /
      !      DATA IMACH( 5) /         64 /
      !      DATA IMACH( 6) /          8 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         63 /
      !      DATA IMACH( 9) / 777777777777777777777B /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         47 /
      !      DATA IMACH(12) /      -8189 /
      !      DATA IMACH(13) /       8190 /
      !      DATA IMACH(14) /         94 /
      !      DATA IMACH(15) /      -8099 /
      !      DATA IMACH(16) /       8190 /
      !
      !      MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
      !
      !      DATA IMACH( 1) /         11 /
      !      DATA IMACH( 2) /         12 /
      !      DATA IMACH( 3) /          8 /
      !      DATA IMACH( 4) /         10 /
      !      DATA IMACH( 5) /         16 /
      !      DATA IMACH( 6) /          2 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         15 /
      !      DATA IMACH( 9) /      32767 /
      !      DATA IMACH(10) /         16 /
      !      DATA IMACH(11) /          6 /
      !      DATA IMACH(12) /        -64 /
      !      DATA IMACH(13) /         63 /
      !      DATA IMACH(14) /         14 /
      !      DATA IMACH(15) /        -64 /
      !      DATA IMACH(16) /         63 /
      !
      !      MACHINE CONSTANTS FOR THE DEC ALPHA
      !      USING G_FLOAT
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1023 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE DEC ALPHA
      !      USING IEEE_FLOAT
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE DEC RISC
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE DEC VAX
      !      USING D_FLOATING
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         56 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE DEC VAX
      !      USING G_FLOATING
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1023 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE ELXSI 6400
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         32 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -126 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1022 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE HARRIS 220
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          0 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         24 /
      !      DATA IMACH( 6) /          3 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         23 /
      !      DATA IMACH( 9) /    8388607 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         23 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         38 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /         43 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         36 /
      !      DATA IMACH( 6) /          6 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         35 /
      !      DATA IMACH( 9) / O377777777777 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         27 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         63 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE HP 730
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE HP 2100
      !      3 WORD DOUBLE PRECISION OPTION WITH FTN4
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          4 /
      !      DATA IMACH( 4) /          1 /
      !      DATA IMACH( 5) /         16 /
      !      DATA IMACH( 6) /          2 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         15 /
      !      DATA IMACH( 9) /      32767 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         23 /
      !      DATA IMACH(12) /       -128 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         39 /
      !      DATA IMACH(15) /       -128 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE HP 2100
      !      4 WORD DOUBLE PRECISION OPTION WITH FTN4
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          4 /
      !      DATA IMACH( 4) /          1 /
      !      DATA IMACH( 5) /         16 /
      !      DATA IMACH( 6) /          2 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         15 /
      !      DATA IMACH( 9) /      32767 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         23 /
      !      DATA IMACH(12) /       -128 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         55 /
      !      DATA IMACH(15) /       -128 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE HP 9000
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          7 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         32 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -126 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1015 /
      !      DATA IMACH(16) /       1017 /
      !
      !      MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
      !      THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
      !      THE PERKIN ELMER (INTERDATA) 7/32.
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          7 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) /  Z7FFFFFFF /
      !      DATA IMACH(10) /         16 /
      !      DATA IMACH(11) /          6 /
      !      DATA IMACH(12) /        -64 /
      !      DATA IMACH(13) /         63 /
      !      DATA IMACH(14) /         14 /
      !      DATA IMACH(15) /        -64 /
      !      DATA IMACH(16) /         63 /
      !
      !      MACHINE CONSTANTS FOR THE IBM PC
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          0 /
      !      DATA IMACH( 4) /          0 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE IBM RS 6000
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          0 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE INTEL i860
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         36 /
      !      DATA IMACH( 6) /          5 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         35 /
      !      DATA IMACH( 9) / "377777777777 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         27 /
      !      DATA IMACH(12) /       -128 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         54 /
      !      DATA IMACH(15) /       -101 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         36 /
      !      DATA IMACH( 6) /          5 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         35 /
      !      DATA IMACH( 9) / "377777777777 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         27 /
      !      DATA IMACH(12) /       -128 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         62 /
      !      DATA IMACH(15) /       -128 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
      !      32-BIT INTEGER ARITHMETIC.
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         56 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
      !      16-BIT INTEGER ARITHMETIC.
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          5 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         16 /
      !      DATA IMACH( 6) /          2 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         15 /
      !      DATA IMACH( 9) /      32767 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         56 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      !      MACHINE CONSTANTS FOR THE SILICON GRAPHICS
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE SUN
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -125 /
      !      DATA IMACH(13) /        128 /
      !      DATA IMACH(14) /         53 /
      !      DATA IMACH(15) /      -1021 /
      !      DATA IMACH(16) /       1024 /
      !
      !      MACHINE CONSTANTS FOR THE SUN
      !      USING THE -r8 COMPILER OPTION
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          6 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         32 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         31 /
      !      DATA IMACH( 9) / 2147483647 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         53 /
      !      DATA IMACH(12) /      -1021 /
      !      DATA IMACH(13) /       1024 /
      !      DATA IMACH(14) /        113 /
      !      DATA IMACH(15) /     -16381 /
      !      DATA IMACH(16) /      16384 /
      !
      !      MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
      !
      !      DATA IMACH( 1) /          5 /
      !      DATA IMACH( 2) /          6 /
      !      DATA IMACH( 3) /          1 /
      !      DATA IMACH( 4) /          6 /
      !      DATA IMACH( 5) /         36 /
      !      DATA IMACH( 6) /          4 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         35 /
      !      DATA IMACH( 9) / O377777777777 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         27 /
      !      DATA IMACH(12) /       -128 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         60 /
      !      DATA IMACH(15) /      -1024 /
      !      DATA IMACH(16) /       1023 /
      !
      !      MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
      !
      !      DATA IMACH( 1) /          1 /
      !      DATA IMACH( 2) /          1 /
      !      DATA IMACH( 3) /          0 /
      !      DATA IMACH( 4) /          1 /
      !      DATA IMACH( 5) /         16 /
      !      DATA IMACH( 6) /          2 /
      !      DATA IMACH( 7) /          2 /
      !      DATA IMACH( 8) /         15 /
      !      DATA IMACH( 9) /      32767 /
      !      DATA IMACH(10) /          2 /
      !      DATA IMACH(11) /         24 /
      !      DATA IMACH(12) /       -127 /
      !      DATA IMACH(13) /        127 /
      !      DATA IMACH(14) /         56 /
      !      DATA IMACH(15) /       -127 /
      !      DATA IMACH(16) /        127 /
      !
      ! ***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
      !
      I1MACH = IMACH(I)
      RETURN
   !
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
      !
      !      CALL FDUMP
      !
      STOP
      END
      ! DECK DHSTRT
      SUBROUTINE DHSTRT (DF, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL,&
         BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
      ! ***BEGIN PROLOGUE  DHSTRT
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (HSTART-S, DHSTRT-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DHSTRT computes a starting step size to be used in solving initial
      !    value problems in ordinary differential equations.
      !
      !  **********************************************************************
      !   ABSTRACT
      !
      !      Subroutine DHSTRT computes a starting step size to be used by an
      !      initial value method in solving ordinary differential equations.
      !      It is based on an estimate of the local Lipschitz constant for the
      !      differential equation   (lower bound on a norm of the Jacobian) ,
      !      a bound on the differential equation  (first derivative) , and
      !      a bound on the partial derivative of the equation with respect to
      !      the independent variable.
      !      (all approximated near the initial point A)
      !
      !      Subroutine DHSTRT uses a function subprogram DHVNRM for computing
      !      a vector norm. The maximum norm is presently utilized though it
      !      can easily be replaced by any other vector norm. It is presumed
      !      that any replacement norm routine would be carefully coded to
      !      prevent unnecessary underflows or overflows from occurring, and
      !      also, would not alter the vector or number of components.
      !
      !  **********************************************************************
      !   On input you must provide the following
      !
      !       DF -- This is a subroutine of the form
      !                                DF(X,U,UPRIME,RPAR,IPAR)
      !              which defines the system of first order differential
      !              equations to be solved. For the given values of X and the
      !              vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
      !              evaluate the NEQ components of the system of differential
      !              equations  DU/DX=DF(X,U)  and store the derivatives in the
      !              array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
      !              equations I=1,...,NEQ.
      !
      !              Subroutine DF must not alter X or U(*). You must declare
      !              the name DF in an external statement in your program that
      !              calls DHSTRT. You must dimension U and UPRIME in DF.
      !
      !              RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
      !              arrays which you can use for communication between your
      !              program and subroutine DF. They are not used or altered by
      !              DHSTRT. If you do not need RPAR or IPAR, ignore these
      !              parameters by treating them as dummy arguments. If you do
      !              choose to use them, dimension them in your program and in
      !              DF as arrays of appropriate length.
      !
      !       NEQ -- This is the number of (first order) differential equations
      !              to be integrated.
      !
      !       A -- This is the initial point of integration.
      !
      !       B -- This is a value of the independent variable used to define
      !              the direction of integration. A reasonable choice is to
      !              set  B  to the first point at which a solution is desired.
      !              You can also use  B, if necessary, to restrict the length
      !              of the first integration step because the algorithm will
      !              not compute a starting step length which is bigger than
      !              ABS(B-A), unless  B  has been chosen too close to  A.
      !              (it is presumed that DHSTRT has been called with  B
      !              different from  A  on the machine being used. Also see the
      !              discussion about the parameter  SMALL.)
      !
      !       Y(*) -- This is the vector of initial values of the NEQ solution
      !              components at the initial point  A.
      !
      !       YPRIME(*) -- This is the vector of derivatives of the NEQ
      !              solution components at the initial point  A.
      !              (defined by the differential equations in subroutine DF)
      !
      !       ETOL -- This is the vector of error tolerances corresponding to
      !              the NEQ solution components. It is assumed that all
      !              elements are positive. Following the first integration
      !              step, the tolerances are expected to be used by the
      !              integrator in an error test which roughly requires that
      !                         ABS(LOCAL ERROR)  .LE.  ETOL
      !              for each vector component.
      !
      !       MORDER -- This is the order of the formula which will be used by
      !              the initial value method for taking the first integration
      !              step.
      !
      !       SMALL -- This is a small positive machine dependent constant
      !              which is used for protecting against computations with
      !              numbers which are too small relative to the precision of
      !              floating point arithmetic.  SMALL  should be set to
      !              (approximately) the smallest positive DOUBLE PRECISION
      !              number such that  (1.+SMALL) .GT. 1.  on the machine being
      !              used. The quantity  SMALL**(3/8)  is used in computing
      !              increments of variables for approximating derivatives by
      !              differences.  Also the algorithm will not compute a
      !              starting step length which is smaller than
      !              100*SMALL*ABS(A).
      !
      !       BIG -- This is a large positive machine dependent constant which
      !              is used for preventing machine overflows. A reasonable
      !              choice is to set big to (approximately) the square root of
      !              the largest DOUBLE PRECISION number which can be held in
      !              the machine.
      !
      !       SPY(*),PV(*),YP(*),SF(*) -- These are DOUBLE PRECISION work
      !              arrays of length NEQ which provide the routine with needed
      !              storage space.
      !
      !       RPAR,IPAR -- These are parameter arrays, of DOUBLE PRECISION and
      !              INTEGER type, respectively, which can be used for
      !              communication between your program and the DF subroutine.
      !              They are not used or altered by DHSTRT.
      !
      !  **********************************************************************
      !   On Output  (after the return from DHSTRT),
      !
      !       H -- is an appropriate starting step size to be attempted by the
      !              differential equation method.
      !
      !            All parameters in the call list remain unchanged except for
      !            the working arrays SPY(*),PV(*),YP(*), and SF(*).
      !
      !  **********************************************************************
      !
      ! ***SEE ALSO  DDEABM, DDEBDF, DDERKF
      ! ***ROUTINES CALLED  DHVNRM
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891024  Changed references from DVNORM to DHVNRM.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      ! ***END PROLOGUE  DHSTRT
      !
      INTEGER IPAR, J, K, LK, MORDER, NEQ
      DOUBLE PRECISION A, ABSDX, B, BIG, DA, DELF, DELY,&
            DFDUB, DFDXB, DHVNRM,&
            DX, DY, ETOL, FBND, H, PV, RELPER, RPAR, SF, SMALL, SPY,&
            SRYDPB, TOLEXP, TOLMIN, TOLP, TOLSUM, Y, YDPB, YP, YPRIME
      DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),&
                SF(*),RPAR(*),IPAR(*)
      EXTERNAL DF
         !
         !      ..................................................................
         !
         !      BEGIN BLOCK PERMITTING ...EXITS TO 160
         ! ***FIRST EXECUTABLE STATEMENT  DHSTRT
         DX = B - A
         ABSDX = ABS(DX)
         RELPER = SMALL**0.375D0
         !
         !         ...............................................................
         !
         !              COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
         !              DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
         !              INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW.
         !              ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE
         !              LOCALLY.
         !
         DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),&
                          100.0D0*SMALL*ABS(A)),DX)
         IF (DA .EQ. 0.0D0) DA = RELPER*DX
         CALL DF(A+DA,Y,SF,RPAR,IPAR)
         DO 10 J = 1, NEQ
            YP(J) = SF(J) - YPRIME(J)
   10    CONTINUE
         DELF = DHVNRM(YP,NEQ)
         DFDXB = BIG
         IF (DELF .LT. BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
         FBND = DHVNRM(SF,NEQ)
         !
         !         ...............................................................
         !
         !              COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ
         !              CONSTANT FOR THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS
         !              ALSO REPRESENTS AN ESTIMATE OF THE NORM OF THE JACOBIAN
         !              LOCALLY.  THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO
         !              ESTIMATE THE LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES.
         !              THE FIRST PERTURBATION VECTOR IS BASED ON THE INITIAL
         !              DERIVATIVES AND DIRECTION OF INTEGRATION. THE SECOND
         !              PERTURBATION VECTOR IS FORMED USING ANOTHER EVALUATION OF
         !              THE DIFFERENTIAL EQUATION.  THE THIRD PERTURBATION VECTOR
         !              IS FORMED USING PERTURBATIONS BASED ONLY ON THE INITIAL
         !              VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS CHANGED TO
         !              NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
         !              INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT
         !              COMPONENTS OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE
         !              CONSISTENT WITH THE SLOPES OF LOCAL SOLUTION CURVES.
         !              ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST
         !              DERIVATIVE.
         !
         !                                PERTURBATION VECTOR SIZE IS HELD
         !                                CONSTANT FOR ALL ITERATIONS. COMPUTE
         !                                THIS CHANGE FROM THE
         !                                        SIZE OF THE VECTOR OF INITIAL
         !                                        VALUES.
         DELY = RELPER*DHVNRM(Y,NEQ)
         IF (DELY .EQ. 0.0D0) DELY = RELPER
         DELY = SIGN(DELY,DX)
         DELF = DHVNRM(YPRIME,NEQ)
         FBND = MAX(FBND,DELF)
         IF (DELF .EQ. 0.0D0) GO TO 30
            !            USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
            DO 20 J = 1, NEQ
               SPY(J) = YPRIME(J)
               YP(J) = YPRIME(J)
   20       CONTINUE
         GO TO 50
   30    CONTINUE
            !            CANNOT HAVE A NULL PERTURBATION VECTOR
            DO 40 J = 1, NEQ
               SPY(J) = 0.0D0
               YP(J) = 1.0D0
   40       CONTINUE
            DELF = DHVNRM(YP,NEQ)
   50    CONTINUE
         !
         DFDUB = 0.0D0
         LK = MIN(NEQ+1,3)
         DO 140 K = 1, LK
            !            DEFINE PERTURBED VECTOR OF INITIAL VALUES
            DO 60 J = 1, NEQ
               PV(J) = Y(J) + DELY*(YP(J)/DELF)
   60       CONTINUE
            IF (K .EQ. 2) GO TO 80
               !               EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
               !               VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
               CALL DF(A,PV,YP,RPAR,IPAR)
               DO 70 J = 1, NEQ
                  PV(J) = YP(J) - YPRIME(J)
   70          CONTINUE
            GO TO 100
   80       CONTINUE
               !               USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
               !                                     IN COMPUTING ONE ESTIMATE
               CALL DF(A+DA,PV,YP,RPAR,IPAR)
               DO 90 J = 1, NEQ
                  PV(J) = YP(J) - SF(J)
   90          CONTINUE
  100       CONTINUE
            !            CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE
            !                           AND A LOCAL LIPSCHITZ CONSTANT
            FBND = MAX(FBND,DHVNRM(YP,NEQ))
            DELF = DHVNRM(PV,NEQ)
            !         ...EXIT
            IF (DELF .GE. BIG*ABS(DELY)) GO TO 150
            DFDUB = MAX(DFDUB,DELF/ABS(DELY))
            !      ......EXIT
            IF (K .EQ. LK) GO TO 160
            !            CHOOSE NEXT PERTURBATION VECTOR
            IF (DELF .EQ. 0.0D0) DELF = 1.0D0
            DO 130 J = 1, NEQ
               IF (K .EQ. 2) GO TO 110
                  DY = ABS(PV(J))
                  IF (DY .EQ. 0.0D0) DY = DELF
               GO TO 120
  110          CONTINUE
                  DY = Y(J)
                  IF (DY .EQ. 0.0D0) DY = DELY/RELPER
  120          CONTINUE
               IF (SPY(J) .EQ. 0.0D0) SPY(J) = YP(J)
               IF (SPY(J) .NE. 0.0D0) DY = SIGN(DY,SPY(J))
               YP(J) = DY
  130       CONTINUE
            DELF = DHVNRM(YP,NEQ)
  140    CONTINUE
  150    CONTINUE
         !
         !         PROTECT AGAINST AN OVERFLOW
         DFDUB = BIG
  160 CONTINUE
      !
      !      ..................................................................
      !
      !           COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
      !
      YDPB = DFDXB + DFDUB*FBND
      !
      !      ..................................................................
      !
      !           DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP
      !           SIZE IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR
      !           TOLERANCE RANGE IS SELECTED.
      !
      TOLMIN = BIG
      TOLSUM = 0.0D0
      DO 170 K = 1, NEQ
         TOLEXP = LOG10(ETOL(K))
         TOLMIN = MIN(TOLMIN,TOLEXP)
         TOLSUM = TOLSUM + TOLEXP
  170 CONTINUE
      TOLP = 10.0D0**(0.5D0*(TOLSUM/NEQ + TOLMIN)/(MORDER+1))
      !
      !      ..................................................................
      !
      !           COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND
      !           SECOND DERIVATIVE INFORMATION
      !
      !                             RESTRICT THE STEP LENGTH TO BE NOT BIGGER
      !                             THAN ABS(B-A).   (UNLESS  B  IS TOO CLOSE
      !                             TO  A)
      H = ABSDX
      !
      IF (YDPB .NE. 0.0D0 .OR. FBND .NE. 0.0D0) GO TO 180
         !
         !         BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
         !                      DERIVATIVE TERM (YDPB) ARE ZERO
         IF (TOLP .LT. 1.0D0) H = ABSDX*TOLP
      GO TO 200
  180 CONTINUE
      !
      IF (YDPB .NE. 0.0D0) GO TO 190
         !
         !         ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
         IF (TOLP .LT. FBND*ABSDX) H = TOLP/FBND
      GO TO 200
  190 CONTINUE
         !
         !         SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
         SRYDPB = SQRT(0.5D0*YDPB)
         IF (TOLP .LT. SRYDPB*ABSDX) H = TOLP/SRYDPB
  200 CONTINUE
      !
      !      FURTHER RESTRICT THE STEP LENGTH TO BE NOT
      !                                BIGGER THAN  1/DFDUB
      IF (H*DFDUB .GT. 1.0D0) H = 1.0D0/DFDUB
      !
      !      FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
      !      SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
      !      A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
      !      THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
      !                                      STEP LENGTH.
      H = MAX(H,100.0D0*SMALL*ABS(A))
      IF (H .EQ. 0.0D0) H = SMALL*ABS(B)
      !
      !      NOW SET DIRECTION OF INTEGRATION
      H = SIGN(H,DX)
      !
      RETURN
      END
      ! DECK DHVNRM
      DOUBLE PRECISION FUNCTION DHVNRM (V, NCOMP)
      ! ***BEGIN PROLOGUE  DHVNRM
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !      Compute the maximum norm of the vector V(*) of length NCOMP and
      !      return the result as DHVNRM
      !
      ! ***SEE ALSO  DDEABM, DDEBDF, DDERKF
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891024  Changed references from DVNORM to DHVNRM.  (WRB)
      !    891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      ! ***END PROLOGUE  DHVNRM
      !
      INTEGER K, NCOMP
      DOUBLE PRECISION V
      DIMENSION V(*)
      ! ***FIRST EXECUTABLE STATEMENT  DHVNRM
      DHVNRM = 0.0D0
      DO 10 K = 1, NCOMP
         DHVNRM = MAX(DHVNRM,ABS(V(K)))
   10 CONTINUE
      RETURN
      END
      ! DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
      ! ***BEGIN PROLOGUE  J4SAVE
      ! ***SUBSIDIARY
      ! ***PURPOSE  Save or recall global variables needed by error
      !             handling routines.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***TYPE      INTEGER (J4SAVE-I)
      ! ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !      Abstract
      !         J4SAVE saves and recalls several global variables needed
      !         by the library error handling routines.
      !
      !      Description of Parameters
      !       --Input--
      !         IWHICH - Index of item desired.
      !                 = 1 Refers to current error number.
      !                 = 2 Refers to current error control flag.
      !                 = 3 Refers to current unit number to which error
      !                     messages are to be sent.  (0 means use standard.)
      !                 = 4 Refers to the maximum number of times any
      !                      message is to be printed (as set by XERMAX).
      !                 = 5 Refers to the total number of units to which
      !                      each error message is to be written.
      !                 = 6 Refers to the 2nd unit for error messages
      !                 = 7 Refers to the 3rd unit for error messages
      !                 = 8 Refers to the 4th unit for error messages
      !                 = 9 Refers to the 5th unit for error messages
      !         IVALUE - The value to be set for the IWHICH-th parameter,
      !                  if ISET is .TRUE. .
      !         ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
      !                  given the value, IVALUE.  If ISET=.FALSE., the
      !                  IWHICH-th parameter will be unchanged, and IVALUE
      !                  is a dummy parameter.
      !       --Output--
      !         The (old) value of the IWHICH-th parameter will be returned
      !         in the function value, J4SAVE.
      !
      ! ***SEE ALSO  XERMSG
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    790801  DATE WRITTEN
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900205  Minor modifications to prologue.  (WRB)
      !    900402  Added TYPE section.  (WRB)
      !    910411  Added KEYWORDS section.  (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
      ! ***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      ! DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
      ! ***BEGIN PROLOGUE  XERCNT
      ! ***SUBSIDIARY
      ! ***PURPOSE  Allow user control over handling of errors.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3C
      ! ***TYPE      ALL (XERCNT-A)
      ! ***KEYWORDS  ERROR, XERROR
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !      Abstract
      !         Allows user control over handling of individual errors.
      !         Just after each message is recorded, but before it is
      !         processed any further (i.e., before it is printed or
      !         a decision to abort is made), a call is made to XERCNT.
      !         If the user has provided his own version of XERCNT, he
      !         can then override the value of KONTROL used in processing
      !         this message by redefining its value.
      !         KONTRL may be set to any value from -2 to 2.
      !         The meanings for KONTRL are the same as in XSETF, except
      !         that the value of KONTRL changes only for this message.
      !         If KONTRL is set to a value outside the range from -2 to 2,
      !         it will be moved back into that range.
      !
      !      Description of Parameters
      !
      !       --Input--
      !         LIBRAR - the library that the routine is in.
      !         SUBROU - the subroutine that XERMSG is being called from
      !         MESSG  - the first 20 characters of the error message.
      !         NERR   - same as in the call to XERMSG.
      !         LEVEL  - same as in the call to XERMSG.
      !         KONTRL - the current value of the control flag as set
      !                  by a call to XSETF.
      !
      !       --Output--
      !         KONTRL - the new value of KONTRL.  If KONTRL is not
      !                  defined, it will remain at its original value.
      !                  This changed value of control affects only
      !                  the current occurrence of the current message.
      !
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    790801  DATE WRITTEN
      !    861211  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900206  Routine changed from user-callable to subsidiary.  (WRB)
      !    900510  Changed calling sequence to include LIBRARY and SUBROUTINE
      !            names, changed routine name from XERCTL to XERCNT.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      ! ***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
      ! DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
      ! ***BEGIN PROLOGUE  XERHLT
      ! ***SUBSIDIARY
      ! ***PURPOSE  Abort program execution and print error message.
      ! ***LIBRARY   SLATEC (XERROR)
      ! ***CATEGORY  R3C
      ! ***TYPE      ALL (XERHLT-A)
      ! ***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
      ! ***AUTHOR  Jones, R. E., (SNLA)
      ! ***DESCRIPTION
      !
      !      Abstract
      !         ***Note*** machine dependent routine
      !         XERHLT aborts the execution of the program.
      !         The error message causing the abort is given in the calling
      !         sequence, in case one needs it for printing on a dayfile,
      !         for example.
      !
      !      Description of Parameters
      !         MESSG is as in XERMSG.
      !
      ! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
      !                  Error-handling Package, SAND82-0800, Sandia
      !                  Laboratories, 1982.
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    790801  DATE WRITTEN
      !    861211  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900206  Routine changed from user-callable to subsidiary.  (WRB)
      !    900510  Changed calling sequence to delete length of character
      !            and changed routine name from XERABT to XERHLT.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
      ! ***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
