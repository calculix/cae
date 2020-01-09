      !
      !    SLATEC: public domain
      !
      ! DECK DDEBDF
      SUBROUTINE DDEBDF (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,&
         RWORK, LRW, IWORK, LIW, RPAR, IPAR, DJAC)
      ! ***BEGIN PROLOGUE  DDEBDF
      ! ***PURPOSE  Solve an initial value problem in ordinary differential
      !             equations using backward differentiation formulas.  It is
      !             intended primarily for stiff problems.
      ! ***LIBRARY   SLATEC (DEPAC)
      ! ***CATEGORY  I1A2
      ! ***TYPE      DOUBLE PRECISION (DEBDF-S, DDEBDF-D)
      ! ***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DEPAC,
      !              INITIAL VALUE PROBLEMS, ODE,
      !              ORDINARY DIFFERENTIAL EQUATIONS, STIFF
      ! ***AUTHOR  Shampine, L. F., (SNLA)
      !            Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    This is the backward differentiation code in the package of
      !    differential equation solvers DEPAC, consisting of the codes
      !    DDERKF, DDEABM, and DDEBDF.  Design of the package was by
      !    L. F. Shampine and H. A. Watts.  It is documented in
      !         SAND-79-2374 , DEPAC - Design of a User Oriented Package of ODE
      !                               Solvers.
      !    DDEBDF is a driver for a modification of the code LSODE written by
      !              A. C. Hindmarsh
      !              Lawrence Livermore Laboratory
      !              Livermore, California 94550
      !
      !  **********************************************************************
      !  **             DEPAC PACKAGE OVERVIEW           **
      !  **********************************************************************
      !
      !         You have a choice of three differential equation solvers from
      !         DEPAC.  The following brief descriptions are meant to aid you
      !         in choosing the most appropriate code for your problem.
      !
      !         DDERKF is a fifth order Runge-Kutta code. It is the simplest of
      !         the three choices, both algorithmically and in the use of the
      !         code. DDERKF is primarily designed to solve non-stiff and mild-
      !         ly stiff differential equations when derivative evaluations are
      !         not expensive.  It should generally not be used to get high
      !         accuracy results nor answers at a great many specific points.
      !         Because DDERKF has very low overhead costs, it will usually
      !         result in the least expensive integration when solving
      !         problems requiring a modest amount of accuracy and having
      !         equations that are not costly to evaluate.  DDERKF attempts to
      !         discover when it is not suitable for the task posed.
      !
      !         DDEABM is a variable order (one through twelve) Adams code. Its
      !         complexity lies somewhere between that of DDERKF and DDEBDF.
      !         DDEABM is primarily designed to solve non-stiff and mildly
      !         stiff differential equations when derivative evaluations are
      !         expensive, high accuracy results are needed or answers at
      !         many specific points are required.  DDEABM attempts to discover
      !         when it is not suitable for the task posed.
      !
      !         DDEBDF is a variable order (one through five) backward
      !         differentiation formula code.  It is the most complicated of
      !         the three choices.  DDEBDF is primarily designed to solve stiff
      !         differential equations at crude to moderate tolerances.
      !         If the problem is very stiff at all, DDERKF and DDEABM will be
      !         quite inefficient compared to DDEBDF.  However, DDEBDF will be
      !         inefficient compared to DDERKF and DDEABM on non-stiff problems
      !         because it uses much more storage, has a much larger overhead,
      !         and the low order formulas will not give high accuracies
      !         efficiently.
      !
      !         The concept of stiffness cannot be described in a few words.
      !         If you do not know the problem to be stiff, try either DDERKF
      !         or DDEABM.  Both of these codes will inform you of stiffness
      !         when the cost of solving such problems becomes important.
      !
      !  **********************************************************************
      !  ** ABSTRACT **
      !  **********************************************************************
      !
      !    Subroutine DDEBDF uses the backward differentiation formulas of
      !    orders one through five to integrate a system of NEQ first order
      !    ordinary differential equations of the form
      !                          DU/DX = DF(X,U)
      !    when the vector Y(*) of initial values for U(*) at X=T is given.
      !    The subroutine integrates from T to TOUT. It is easy to continue the
      !    integration to get results at additional TOUT. This is the interval
      !    mode of operation. It is also easy for the routine to return with
      !    the solution at each intermediate step on the way to TOUT. This is
      !    the intermediate-output mode of operation.
      !
      !  **********************************************************************
      !  * Description of The Arguments To DDEBDF (An Overview) *
      !  **********************************************************************
      !
      !    The Parameters are:
      !
      !       DF -- This is the name of a subroutine which you provide to
      !             define the differential equations.
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
      !       RTOL, ATOL -- These DOUBLE PRECISION quantities
      !              represent relative and absolute error tolerances which you
      !              provide to indicate how accurately you wish the solution
      !              to be computed.  You may choose them to be both scalars
      !              or else both vectors.
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
      !              calling program and the DF subroutine (and the DJAC
      !              subroutine).
      !
      !       DJAC -- This is the name of a subroutine which you may choose to
      !              provide for defining the Jacobian matrix of partial
      !              derivatives DF/DU.
      !
      !   Quantities which are used as input items are
      !              NEQ, T, Y(*), TOUT, INFO(*),
      !              RTOL, ATOL, RWORK(1), LRW,
      !              IWORK(1), IWORK(2), and LIW.
      !
      !   Quantities which may be altered by the code are
      !              T, Y(*), INFO(1), RTOL, ATOL,
      !              IDID, RWORK(*) and IWORK(*).
      !
      !  **********************************************************************
      !  * INPUT -- What To Do On The First Call To DDEBDF *
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
      !              which is to be solved. For the given values of X and the
      !              vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
      !              evaluate the NEQ components of the system of differential
      !              equations  DU/DX=DF(X,U)  and store the derivatives in the
      !              array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
      !              equations I=1,...,NEQ.
      !
      !              Subroutine DF must not alter X or U(*). You must declare
      !              the name DF in an external statement in your program that
      !              calls DDEBDF. You must dimension U and UPRIME in DF.
      !
      !              RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
      !              arrays which you can use for communication between your
      !              calling program and subroutine DF. They are not used or
      !              altered by DDEBDF.  If you do not need RPAR or IPAR,
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
      !       TOUT -- Set it to the first point at which a solution is desired.
      !              You can take TOUT = T, in which case the code
      !              will evaluate the derivative of the solution at T and
      !              return.  Integration either forward in T  (TOUT .GT. T)
      !              or backward in T  (TOUT .LT. T)  is permitted.
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
      !              the problem.  By using the fact that the code will not
      !              step past TOUT in the first step, you could, if necessary,
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
      !              DEPAC or possible future extensions, though DDEBDF uses
      !              only the first six entries.  You must respond to all of
      !              the following items which are arranged as questions.  The
      !              simplest use of the code corresponds to answering all
      !              questions as YES ,i.e. setting all entries of INFO to 0.
      !
      !         INFO(1) -- This parameter enables the code to initialize
      !                itself.  You must set it to indicate the start of every
      !                new problem.
      !
      !             **** Is this the first call for this problem ...
      !                   YES -- Set INFO(1) = 0
      !                    NO -- Not applicable here.
      !                          See below for continuation calls.  ****
      !
      !         INFO(2) -- How much accuracy you want of your solution
      !                is specified by the error tolerances RTOL and ATOL.
      !                The simplest use is to take them both to be scalars.
      !                To obtain more flexibility, they can both be vectors.
      !                The code must be told your choice.
      !
      !             **** Are both error tolerances RTOL, ATOL scalars ...
      !                   YES -- Set INFO(2) = 0
      !                          and input scalars for both RTOL and ATOL
      !                    NO -- Set INFO(2) = 1
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
      !                  TOUT (and NOT at the next intermediate step) ...
      !                   YES -- Set INFO(3) = 0
      !                    NO -- Set INFO(3) = 1 ****
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
      !                  restrictions on the independent variable T ...
      !                   YES -- Set INFO(4)=0
      !                    NO -- Set INFO(4)=1
      !                          and define the stopping point TSTOP by
      !                          setting RWORK(1)=TSTOP ****
      !
      !         INFO(5) -- To solve stiff problems it is necessary to use the
      !                Jacobian matrix of partial derivatives of the system
      !                of differential equations.  If you do not provide a
      !                subroutine to evaluate it analytically (see the
      !                description of the item DJAC in the call list), it will
      !                be approximated by numerical differencing in this code.
      !                Although it is less trouble for you to have the code
      !                compute partial derivatives by numerical differencing,
      !                the solution will be more reliable if you provide the
      !                derivatives via DJAC.  Sometimes numerical differencing
      !                is cheaper than evaluating derivatives in DJAC and
      !                sometimes it is not - this depends on your problem.
      !
      !                If your problem is linear, i.e. has the form
      !                DU/DX = DF(X,U) = J(X)*U + G(X)   for some matrix J(X)
      !                and vector G(X), the Jacobian matrix  DF/DU = J(X).
      !                Since you must provide a subroutine to evaluate DF(X,U)
      !                analytically, it is little extra trouble to provide
      !                subroutine DJAC for evaluating J(X) analytically.
      !                Furthermore, in such cases, numerical differencing is
      !                much more expensive than analytic evaluation.
      !
      !             **** Do you want the code to evaluate the partial
      !                  derivatives automatically by numerical differences ...
      !                   YES -- Set INFO(5)=0
      !                    NO -- Set INFO(5)=1
      !                          and provide subroutine DJAC for evaluating the
      !                          Jacobian matrix ****
      !
      !         INFO(6) -- DDEBDF will perform much better if the Jacobian
      !                matrix is banded and the code is told this.  In this
      !                case, the storage needed will be greatly reduced,
      !                numerical differencing will be performed more cheaply,
      !                and a number of important algorithms will execute much
      !                faster.  The differential equation is said to have
      !                half-bandwidths ML (lower) and MU (upper) if equation I
      !                involves only unknowns Y(J) with
      !                               I-ML .LE. J .LE. I+MU
      !                for all I=1,2,...,NEQ.  Thus, ML and MU are the widths
      !                of the lower and upper parts of the band, respectively,
      !                with the main diagonal being excluded.  If you do not
      !                indicate that the equation has a banded Jacobian,
      !                the code works with a full matrix of NEQ**2 elements
      !                (stored in the conventional way).  Computations with
      !                banded matrices cost less time and storage than with
      !                full matrices if  2*ML+MU .LT. NEQ.  If you tell the
      !                code that the Jacobian matrix has a banded structure and
      !                you want to provide subroutine DJAC to compute the
      !                partial derivatives, then you must be careful to store
      !                the elements of the Jacobian matrix in the special form
      !                indicated in the description of DJAC.
      !
      !             **** Do you want to solve the problem using a full
      !                  (dense) Jacobian matrix (and not a special banded
      !                  structure) ...
      !                   YES -- Set INFO(6)=0
      !                    NO -- Set INFO(6)=1
      !                          and provide the lower (ML) and upper (MU)
      !                          bandwidths by setting
      !                          IWORK(1)=ML
      !                          IWORK(2)=MU ****
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
      !              (More specifically, a root-mean-square norm is used to
      !              measure the size of vectors, and the error test uses the
      !              magnitude of the solution at the beginning of the step.)
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
      !              Setting ATOL=0. results in a pure relative error test on
      !              that component.  Setting RTOL=0. results in a pure abso-
      !              lute error test on that component.  A mixed test with non-
      !              zero RTOL and ATOL corresponds roughly to a relative error
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
      !              (For some problems it may not be permissible to integrate
      !              past a point TSTOP because a discontinuity occurs there
      !              or the solution or its derivative is not defined beyond
      !              TSTOP.)
      !
      !       LRW -- Set it to the declared length of the RWORK array.
      !              You must have
      !                   LRW .GE. 250+10*NEQ+NEQ**2
      !              for the full (dense) Jacobian case (when INFO(6)=0),  or
      !                   LRW .GE. 250+10*NEQ+(2*ML+MU+1)*NEQ
      !              for the banded Jacobian case (when INFO(6)=1).
      !
      !       IWORK(*) -- Dimension this INTEGER work array of length LIW in
      !              your calling program.
      !
      !       IWORK(1), IWORK(2) -- If you have set INFO(6)=0, you can ignore
      !              these optional input parameters. Otherwise you must define
      !              the half-bandwidths ML (lower) and MU (upper) of the
      !              Jacobian matrix by setting    IWORK(1) = ML   and
      !              IWORK(2) = MU.  (The code will work with a full matrix
      !              of NEQ**2 elements unless it is told that the problem has
      !              a banded Jacobian, in which case the code will work with
      !              a matrix containing at most  (2*ML+MU+1)*NEQ  elements.)
      !
      !       LIW -- Set it to the declared length of the IWORK array.
      !              You must have LIW .GE. 56+NEQ.
      !
      !       RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
      !              INTEGER type, respectively. You can use them for
      !              communication between your program that calls DDEBDF and
      !              the  DF subroutine (and the DJAC subroutine). They are not
      !              used or altered by DDEBDF. If you do not need RPAR or
      !              IPAR, ignore these parameters by treating them as dummy
      !              arguments. If you do choose to use them, dimension them in
      !              your calling program and in DF (and in DJAC) as arrays of
      !              appropriate length.
      !
      !       DJAC -- If you have set INFO(5)=0, you can ignore this parameter
      !              by treating it as a dummy argument. (For some compilers
      !              you may have to write a dummy subroutine named  DJAC  in
      !              order to avoid problems associated with missing external
      !              routine names.)  Otherwise, you must provide a subroutine
      !              of the form
      !                           DJAC(X,U,PD,NROWPD,RPAR,IPAR)
      !              to define the Jacobian matrix of partial derivatives DF/DU
      !              of the system of differential equations   DU/DX = DF(X,U).
      !              For the given values of X and the vector
      !              U(*)=(U(1),U(2),...,U(NEQ)), the subroutine must evaluate
      !              the non-zero partial derivatives  DF(I)/DU(J)  for each
      !              differential equation I=1,...,NEQ and each solution
      !              component J=1,...,NEQ , and store these values in the
      !              matrix PD.  The elements of PD are set to zero before each
      !              call to DJAC so only non-zero elements need to be defined.
      !
      !              Subroutine DJAC must not alter X, U(*), or NROWPD. You
      !              must declare the name DJAC in an external statement in
      !              your program that calls DDEBDF. NROWPD is the row
      !              dimension of the PD matrix and is assigned by the code.
      !              Therefore you must dimension PD in DJAC according to
      !                               DIMENSION PD(NROWPD,1)
      !              You must also dimension U in DJAC.
      !
      !              The way you must store the elements into the PD matrix
      !              depends on the structure of the Jacobian which you
      !              indicated by INFO(6).
      !              *** INFO(6)=0 -- Full (Dense) Jacobian ***
      !                  When you evaluate the (non-zero) partial derivative
      !                  of equation I with respect to variable J, you must
      !                  store it in PD according to
      !                                 PD(I,J) = * DF(I)/DU(J) *
      !              *** INFO(6)=1 -- Banded Jacobian with ML Lower and MU
      !                  Upper Diagonal Bands (refer to INFO(6) description of
      !                  ML and MU) ***
      !                  When you evaluate the (non-zero) partial derivative
      !                  of equation I with respect to variable J, you must
      !                  store it in PD according to
      !                                 IROW = I - J + ML + MU + 1
      !                                 PD(IROW,J) = * DF(I)/DU(J) *
      !
      !              RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
      !              arrays which you can use for communication between your
      !              calling program and your Jacobian subroutine DJAC. They
      !              are not altered by DDEBDF. If you do not need RPAR or
      !              IPAR, ignore these parameters by treating them as dummy
      !              arguments.  If you do choose to use them, dimension them
      !              in your calling program and in DJAC as arrays of
      !              appropriate length.
      !
      !  **********************************************************************
      !  * OUTPUT -- After any return from DDEBDF *
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
      !              IDID = -4,-5  -- Not applicable for this code but used
      !                        by other members of DEPAC.
      !
      !              IDID = -6 -- DDEBDF had repeated convergence test failures
      !                        on the last attempted step.
      !
      !              IDID = -7 -- DDEBDF had repeated error test failures on
      !                        the last attempted step.
      !
      !              IDID = -8,..,-32  -- Not applicable for this code but
      !                        used by other members of DEPAC or possible
      !                        future extensions.
      !
      !                          *** Task Terminated ***
      !                    Reported by the value of IDID=-33
      !
      !              IDID = -33 -- The code has encountered trouble from which
      !                        it cannot recover.  A message is printed
      !                        explaining the trouble and control is returned
      !                        to the calling program.  For example, this
      !                        occurs when invalid input is detected.
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
      !              RWORK(12)--If the tolerances have been increased by the
      !                         code (IDID = -2) , they were multiplied by the
      !                         value in RWORK(12).
      !
      !              RWORK(13)--which contains the current value of the
      !                         independent variable, i.e. the farthest point
      !                         integration has reached.  This will be
      !                         different from T only when interpolation has
      !                         been performed (IDID=3).
      !
      !              RWORK(20+I)--which contains the approximate derivative
      !                         of the solution component Y(I).  In DDEBDF, it
      !                         is never obtained by calling subroutine DF to
      !                         evaluate the differential equation using T and
      !                         Y(*), except at the initial point of
      !                         integration.
      !
      !  **********************************************************************
      !  ** INPUT -- What To Do To Continue The Integration **
      !  **             (calls after the first)             **
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
      !         Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
      !         unless you are going to restart the code.
      !
      !         The parameter INFO(1) is used by the code to indicate the
      !         beginning of a new problem and to indicate whether integration
      !         is to be continued.  You must input the value  INFO(1) = 0
      !         when starting a new problem.  You must input the value
      !         INFO(1) = 1  if you wish to continue after an interrupted task.
      !         Do not set  INFO(1) = 0  on a continuation call unless you
      !         want the code to restart at the current T.
      !
      !                          *** Following a Completed Task ***
      !          If
      !              IDID = 1, call the code again to continue the integration
      !                      another step in the direction of TOUT.
      !
      !              IDID = 2 or 3, define a new TOUT and call the code again.
      !                      TOUT must be different from T.  You cannot change
      !                      the direction of integration without restarting.
      !
      !                          *** Following an Interrupted Task ***
      !                      To show the code that you realize the task was
      !                      interrupted and that you want to continue, you
      !                      must take appropriate action and reset INFO(1) = 1
      !          If
      !              IDID = -1, the code has attempted 500 steps.
      !                      If you want to continue, set INFO(1) = 1 and
      !                      call the code again.  An additional 500 steps
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
      !              IDID = -4,-5  --- cannot occur with this code but used
      !                      by other members of DEPAC.
      !
      !              IDID = -6, repeated convergence test failures occurred
      !                      on the last attempted step in DDEBDF.  An inaccu-
      !                      rate Jacobian may be the problem.  If you are
      !                      absolutely certain you want to continue, restart
      !                      the integration at the current T by setting
      !                      INFO(1)=0 and call the code again.
      !
      !              IDID = -7, repeated error test failures occurred on the
      !                      last attempted step in DDEBDF.  A singularity in
      !                      the solution may be present.  You should re-
      !                      examine the problem being solved.  If you are
      !                      absolutely certain you want to continue, restart
      !                      the integration at the current T by setting
      !                      INFO(1)=0 and call the code again.
      !
      !              IDID = -8,..,-32  --- cannot occur with this code but
      !                      used by other members of DDEPAC or possible future
      !                      extensions.
      !
      !                          *** Following a Terminated Task ***
      !          If
      !              IDID = -33, you cannot continue the solution of this
      !                      problem.  An attempt to do so will result in your
      !                      run being terminated.
      !
      !  **********************************************************************
      !
      !          ***** Warning *****
      !
      !      If DDEBDF is to be used in an overlay situation, you must save and
      !      restore certain items used internally by DDEBDF  (values in the
      !      common block DDEBD1).  This can be accomplished as follows.
      !
      !      To save the necessary values upon return from DDEBDF, simply call
      !         DSVCO(RWORK(22+NEQ),IWORK(21+NEQ)).
      !
      !      To restore the necessary values before the next call to DDEBDF,
      !      simply call    DRSCO(RWORK(22+NEQ),IWORK(21+NEQ)).
      !
      ! ***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
      !                  oriented package of ODE solvers, Report SAND79-2374,
      !                  Sandia Laboratories, 1979.
      ! ***ROUTINES CALLED  DLSOD, XERMSG
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    891024  Changed references from DVNORM to DHVNRM.  (WRB)
      !    891024  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900326  Removed duplicate information from DESCRIPTION section.
      !            (WRB)
      !    900510  Convert XERRWV calls to XERMSG calls, make Prologue comments
      !            consistent with DEBDF.  (RWC)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DDEBDF
      INTEGER IACOR, IBAND, IBEGIN, ICOMI, ICOMR, IDELSN, IDID, IER,&
            IEWT, IINOUT, IINTEG, IJAC, ILRW, INFO, INIT,&
            IOWNS, IPAR, IQUIT, ISAVF, ITOL, ITSTAR, ITSTOP, IWM,&
            IWORK, IYH, IYPOUT, JSTART, KFLAG, KSTEPS, L, LIW, LRW,&
            MAXORD, METH, MITER, ML, MU, N, NEQ, NFE, NJE, NQ, NQU,&
            NST
      DOUBLE PRECISION ATOL, EL0, H, HMIN, HMXI, HU, ROWNS, RPAR,&
            RTOL, RWORK, T, TN, TOLD, TOUT, UROUND, Y
      LOGICAL INTOUT
      CHARACTER*8 XERN1, XERN2
      CHARACTER*16 XERN3
      !
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),&
                RPAR(*),IPAR(*)
      !
      COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,&
                      IQUIT,INIT,IYH,IEWT,IACOR,ISAVF,IWM,KSTEPS,IBEGIN,&
                      ITOL,IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART,&
                      KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
      !
      EXTERNAL DF, DJAC
      !
      !         CHECK FOR AN APPARENT INFINITE LOOP
      !
      ! ***FIRST EXECUTABLE STATEMENT  DDEBDF
      IF (INFO(1) .EQ. 0) IWORK(LIW) = 0
      !
      IF (IWORK(LIW).GE. 5) THEN
         IF (T .EQ. RWORK(21+NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDEBDF',&
               'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //&
               'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //&
               ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //&
               'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //&
               'CODE, PARTICULARLY INFO(1).', 13, 2)
            RETURN
         ENDIF
      ENDIF
      !
      IDID = 0
      !
      !         CHECK VALIDITY OF INFO PARAMETERS
      !
      IF (INFO(1) .NE. 0 .AND. INFO(1) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(1) MUST BE SET TO 0 ' //&
            'FOR THE  START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' //&
            'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' //&
            'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE ' //&
            'CODE WITH  INFO(1) = ' // XERN1, 3, 1)
         IDID = -33
      ENDIF
      !
      IF (INFO(2) .NE. 0 .AND. INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(2) MUST BE 0 OR 1 ' //&
            'INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //&
            'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //&
            XERN1, 4, 1)
         IDID = -33
      ENDIF
      !
      IF (INFO(3) .NE. 0 .AND. INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(3) MUST BE 0 OR 1 ' //&
            'INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF ' //&
            'INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED THE CODE ' //&
            'WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID = -33
      ENDIF
      !
      IF (INFO(4) .NE. 0 .AND. INFO(4) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(4)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(4) MUST BE 0 OR 1 ' //&
            'INDICATING WHETHER OR NOT THE INTEGRATION INTERVAL IS ' //&
            'TO BE RESTRICTED BY A POINT TSTOP.  YOU HAVE CALLED ' //&
            'THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
         IDID = -33
      ENDIF
      !
      IF (INFO(5) .NE. 0 .AND. INFO(5) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(5)
         CALL XERMSG ('SLATEC',  'DDEBDF', 'INFO(5) MUST BE 0 OR 1 ' //&
            'INDICATING WHETHER THE CODE IS TOLD TO FORM THE ' //&
            'JACOBIAN MATRIX BY NUMERICAL DIFFERENCING OR YOU ' //&
            'PROVIDE A SUBROUTINE TO EVALUATE IT ANALYTICALLY.  ' //&
            'YOU HAVE CALLED THE CODE WITH INFO(5) = ' // XERN1, 15, 1)
         IDID = -33
      ENDIF
      !
      IF (INFO(6) .NE. 0 .AND. INFO(6) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(6)
         CALL XERMSG ('SLATEC', 'DDEBDF', 'INFO(6) MUST BE 0 OR 1 ' //&
            'INDICATING WHETHER THE CODE IS TOLD TO TREAT THE ' //&
            'JACOBIAN AS A FULL (DENSE) MATRIX OR AS HAVING A ' //&
            'SPECIAL BANDED STRUCTURE.  YOU HAVE CALLED THE CODE ' //&
            'WITH INFO(6) = ' // XERN1, 16, 1)
         IDID = -33
      ENDIF
      !
      ILRW = NEQ
      IF (INFO(6) .NE. 0) THEN
         !
         !         CHECK BANDWIDTH PARAMETERS
         !
         ML = IWORK(1)
         MU = IWORK(2)
         ILRW = 2*ML + MU + 1
         !
         IF (ML.LT.0 .OR. ML.GE.NEQ .OR. MU.LT.0 .OR. MU.GE.NEQ) THEN
            WRITE (XERN1, '(I8)') ML
            WRITE (XERN2, '(I8)') MU
            CALL XERMSG ('SLATEC', 'DDEBDF', 'YOU HAVE SET INFO(6) ' //&
               '= 1, TELLING THE CODE THAT THE JACOBIAN MATRIX HAS ' //&
               'A SPECIAL BANDED STRUCTURE.  HOWEVER, THE LOWER ' //&
               '(UPPER) BANDWIDTHS  ML (MU) VIOLATE THE CONSTRAINTS ' //&
               'ML,MU .GE. 0 AND  ML,MU .LT. NEQ.  YOU HAVE CALLED ' //&
               'THE CODE WITH ML = ' // XERN1 // ' AND MU = ' // XERN2,&
               17, 1)
            IDID = -33
         ENDIF
      ENDIF
      !
      !         CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
      !
      IF (LRW .LT. 250 + (10 + ILRW)*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         IF (INFO(6) .EQ. 0) THEN
            CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY RWORK ' //&
               'MUST BE AT LEAST 250 + 10*NEQ + NEQ*NEQ.$$' //&
               'YOU HAVE CALLED THE CODE WITH  LRW = ' // XERN1, 1, 1)
         ELSE
            CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY RWORK ' //&
               'MUST BE AT LEAST 250 + 10*NEQ + (2*ML+MU+1)*NEQ.$$' //&
               'YOU HAVE CALLED THE CODE WITH  LRW = ' // XERN1, 18, 1)
         ENDIF
         IDID = -33
      ENDIF
      !
      IF (LIW .LT. 56 + NEQ) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DDEBDF', 'LENGTH OF ARRAY IWORK ' //&
            'BE AT LEAST  56 + NEQ.  YOU HAVE CALLED THE CODE WITH ' //&
            'LIW = ' // XERN1, 2, 1)
         IDID = -33
      ENDIF
      !
      !         COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK
      !         ARRAY AND RESTORE COMMON BLOCK DATA
      !
      ICOMI = 21 + NEQ
      IINOUT = ICOMI + 33
      !
      IYPOUT = 21
      ITSTAR = 21 + NEQ
      ICOMR = 22 + NEQ
      !
      IF (INFO(1) .NE. 0) INTOUT = IWORK(IINOUT) .NE. (-1)
      !      CALL DRSCO(RWORK(ICOMR),IWORK(ICOMI))
      !
      IYH = ICOMR + 218
      IEWT = IYH + 6*NEQ
      ISAVF = IEWT + NEQ
      IACOR = ISAVF + NEQ
      IWM = IACOR + NEQ
      IDELSN = IWM + 2 + ILRW*NEQ
      !
      IBEGIN = INFO(1)
      ITOL = INFO(2)
      IINTEG = INFO(3)
      ITSTOP = INFO(4)
      IJAC = INFO(5)
      IBAND = INFO(6)
      RWORK(ITSTAR) = T
      !
      CALL DLSOD(DF,NEQ,T,Y,TOUT,RTOL,ATOL,IDID,RWORK(IYPOUT),&
                 RWORK(IYH),RWORK(IYH),RWORK(IEWT),RWORK(ISAVF),&
                 RWORK(IACOR),RWORK(IWM),IWORK(1),DJAC,INTOUT,&
                 RWORK(1),RWORK(12),RWORK(IDELSN),RPAR,IPAR)
      !
      IWORK(IINOUT) = -1
      IF (INTOUT) IWORK(IINOUT) = 1
      !
      IF (IDID .NE. (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T .NE. RWORK(ITSTAR)) IWORK(LIW) = 0
      !      CALL DSVCO(RWORK(ICOMR),IWORK(ICOMI))
      RWORK(11) = H
      RWORK(13) = TN
      INFO(1) = IBEGIN
      !
      RETURN
      END
      ! DECK DLSOD
      SUBROUTINE DLSOD (DF, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT,&
         YH, YH1, EWT, SAVF, ACOR, WM, IWM, DJAC, INTOUT, TSTOP, TOLFAC,&
         DELSGN, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DLSOD
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (LSOD-S, DLSOD-D)
      ! ***AUTHOR  (UNKNOWN)
      ! ***DESCRIPTION
      !
      !    DDEBDF  merely allocates storage for  DLSOD  to relieve the user of
      !    the inconvenience of a long call list.  Consequently  DLSOD  is used
      !    as described in the comments for  DDEBDF .
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  D1MACH, DHSTRT, DINTYD, DSTOD, DVNRMS, XERMSG
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    900510  Convert XERRWV calls to XERMSG calls.  (RWC)
      ! ***END PROLOGUE  DLSOD
      !
      INTEGER IBAND, IBEGIN, IDID, IER, IINTEG, IJAC, INIT, INTFLG,&
            IOWNS, IPAR, IQUIT, ITOL, ITSTOP, IWM, JSTART, K, KFLAG,&
            KSTEPS, L, LACOR, LDUM, LEWT, LSAVF, LTOL, LWM, LYH, MAXNUM,&
            MAXORD, METH, MITER, N, NATOLP, NEQ, NFE, NJE, NQ, NQU,&
            NRTOLP, NST
      DOUBLE PRECISION ABSDEL, ACOR, ATOL, BIG, D1MACH, DEL,&
            DELSGN, DT, DVNRMS, EL0, EWT,&
            H, HA, HMIN, HMXI, HU, ROWNS, RPAR, RTOL, SAVF, T, TOL,&
            TOLD, TOLFAC, TOUT, TSTOP, U, WM, X, Y, YH, YH1, YPOUT
      LOGICAL INTOUT
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
      !
      DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*),&
                ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      !
      !
      COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,X,U,IQUIT,INIT,&
                      LYH,LEWT,LACOR,LSAVF,LWM,KSTEPS,IBEGIN,ITOL,&
                      IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART,&
                      KFLAG,LDUM,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
      !
      EXTERNAL DF, DJAC
      !
      !      ..................................................................
      !
      !        THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
      !        NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE
      !        COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
      !        EXCESSIVE WORK.
      SAVE MAXNUM
      !
      DATA MAXNUM /500/
      !
      !      ..................................................................
      !
      ! ***FIRST EXECUTABLE STATEMENT  DLSOD
      IF (IBEGIN .EQ. 0) THEN
         !
         !         ON THE FIRST CALL , PERFORM INITIALIZATION --
         !         DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
         !         FUNCTION ROUTINE D1MACH. THE USER MUST MAKE SURE THAT THE
         !         VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
         !
         U = D1MACH(4)
         !                           -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER
         WM(1) = SQRT(U)
         !                           -- SET TERMINATION FLAG
         IQUIT = 0
         !                           -- SET INITIALIZATION INDICATOR
         INIT = 0
         !                           -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS = 0
         !                           -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT = .FALSE.
         !                           -- SET START INDICATOR FOR DSTOD CODE
         JSTART = 0
         !                           -- SET BDF METHOD INDICATOR
         METH = 2
         !                           -- SET MAXIMUM ORDER FOR BDF METHOD
         MAXORD = 5
         !                           -- SET ITERATION MATRIX INDICATOR
         !
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 0) MITER = 2
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 0) MITER = 1
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 1) MITER = 5
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 1) MITER = 4
         !
         !                           -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK
         N = NEQ
         NST = 0
         NJE = 0
         HMXI = 0.0D0
         NQ = 1
         H = 1.0D0
         !                           -- RESET IBEGIN FOR SUBSEQUENT CALLS
         IBEGIN = 1
      ENDIF
      !
      !      ..................................................................
      !
      !       CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
      !
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DLSOD',&
            'IN DDEBDF, THE NUMBER OF EQUATIONS MUST BE A ' //&
            'POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' //&
            XERN1, 6, 1)
         IDID=-33
      ENDIF
      !
      NRTOLP = 0
      NATOLP = 0
      DO 60 K = 1, NEQ
         IF (NRTOLP .LE. 0) THEN
            IF (RTOL(K) .LT. 0.) THEN
               WRITE (XERN1, '(I8)') K
               WRITE (XERN3, '(1PE15.6)') RTOL(K)
               CALL XERMSG ('SLATEC', 'DLSOD',&
                  'IN DDEBDF, THE RELATIVE ERROR TOLERANCES MUST ' //&
                  'BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH ' //&
                  'RTOL(' // XERN1 // ') = ' // XERN3 // '$$IN THE ' //&
                  'CASE OF VECTOR ERROR TOLERANCES, NO FURTHER ' //&
                  'CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
               IDID = -33
               IF (NATOLP .GT. 0) GO TO 70
               NRTOLP = 1
            ELSEIF (NATOLP .GT. 0) THEN
               GO TO 50
            ENDIF
         ENDIF
         !
         IF (ATOL(K) .LT. 0.) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DLSOD',&
               'IN DDEBDF, THE ABSOLUTE ERROR ' //&
               'TOLERANCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED ' //&
               'THE CODE WITH ATOL(' // XERN1 // ') = ' // XERN3 //&
               '$$IN THE CASE OF VECTOR ERROR TOLERANCES, NO FURTHER '&
               // 'CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID=-33
            IF (NRTOLP .GT. 0) GO TO 70
            NATOLP=1
         ENDIF
   50    IF (ITOL .EQ. 0) GO TO 70
   60 CONTINUE
   !
   70 IF (ITSTOP .EQ. 1) THEN
         IF (SIGN(1.0D0,TOUT-T) .NE. SIGN(1.0D0,TSTOP-T) .OR.&
            ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DLSOD',&
               'IN DDEBDF, YOU HAVE CALLED THE ' //&
               'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' //&
               'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' //&
               'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.$$' //&
               'THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
      !
      !         CHECK SOME CONTINUATION POSSIBILITIES
      !
      IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DLSOD',&
               'IN DDEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' //&
               XERN3 // '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.',&
               9, 1)
            IDID=-33
         ENDIF
         !
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DLSOD',&
               'IN DDEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' //&
               XERN3 // ' TO ' // XERN4 //&
               '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
         !
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.0D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DLSOD',&
                  'IN DDEBDF, BY CALLING THE CODE WITH TOUT = ' //&
                  XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' //&
                  'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' //&
                  'WITHOUT RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
      !
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
            !                        INVALID INPUT DETECTED
            IQUIT=-33
            IBEGIN=-1
         ELSE
            CALL XERMSG ('SLATEC', 'DLSOD',&
               'IN DDEBDF, INVALID INPUT WAS DETECTED ON ' //&
               'SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED ' //&
               'BECAUSE YOU HAVE NOT CORRECTED THE PROBLEM, ' //&
               'SO EXECUTION IS BEING TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
      !
      !         ...............................................................
      !
      !              RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED
      !              AS ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS
      !              CASE, THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE
      !              SMALLEST VALUE 100*U WHICH IS LIKELY TO BE REASONABLE FOR
      !              THIS METHOD AND MACHINE
      !
      DO 180 K = 1, NEQ
         IF (RTOL(K) + ATOL(K) .GT. 0.0D0) GO TO 170
            RTOL(K) = 100.0D0*U
            IDID = -2
  170    CONTINUE
         !      ...EXIT
         IF (ITOL .EQ. 0) GO TO 190
  180 CONTINUE
  190 CONTINUE
      !
      IF (IDID .NE. (-2)) GO TO 200
         !         RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
         !                                  SMALL POSITIVE VALUE
         IBEGIN = -1
      GO TO 460
  200 CONTINUE
                     !         BEGIN BLOCK PERMITTING ...EXITS TO 450
                     !            BEGIN BLOCK PERMITTING ...EXITS TO 430
                     !               BEGIN BLOCK PERMITTING ...EXITS TO 260
                     !                  BEGIN BLOCK PERMITTING ...EXITS TO 230
                     !
                     !                     BRANCH ON STATUS OF INITIALIZATION INDICATOR
                     !                            INIT=0 MEANS INITIAL DERIVATIVES AND
                     !                            NOMINAL STEP SIZE
                     !                                   AND DIRECTION NOT YET SET
                     !                            INIT=1 MEANS NOMINAL STEP SIZE AND
                     !                            DIRECTION NOT YET SET INIT=2 MEANS NO
                     !                            FURTHER INITIALIZATION REQUIRED
                     !
                     IF (INIT .EQ. 0) GO TO 210
                        !                  ......EXIT
                        IF (INIT .EQ. 1) GO TO 230
                        !               .........EXIT
                        GO TO 260
  210                CONTINUE
                     !
                     !                     ................................................
                     !
                     !                          MORE INITIALIZATION --
                     !                                              -- EVALUATE INITIAL
                     !                                              DERIVATIVES
                     !
                     INIT = 1
                     CALL DF(T,Y,YH(1,2),RPAR,IPAR)
                     NFE = 1
                     !                  ...EXIT
                     IF (T .NE. TOUT) GO TO 230
                     IDID = 2
                     DO 220 L = 1, NEQ
                        YPOUT(L) = YH(L,2)
  220                CONTINUE
                     TOLD = T
                     !         ............EXIT
                     GO TO 450
  230             CONTINUE
                  !
                  !                  -- COMPUTE INITIAL STEP SIZE
                  !                  -- SAVE SIGN OF INTEGRATION DIRECTION
                  !                  -- SET INDEPENDENT AND DEPENDENT VARIABLES
                  !                                       X AND YH(*) FOR DSTOD
                  !
                  LTOL = 1
                  DO 240 L = 1, NEQ
                     IF (ITOL .EQ. 1) LTOL = L
                     TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
                     IF (TOL .EQ. 0.0D0) GO TO 390
                     EWT(L) = TOL
  240             CONTINUE
                  !
                  BIG = SQRT(D1MACH(2))
                  CALL DHSTRT(DF,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG,&
                              YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR,&
                              IPAR,H)
                  !
                  DELSGN = SIGN(1.0D0,TOUT-T)
                  X = T
                  DO 250 L = 1, NEQ
                     YH(L,1) = Y(L)
                     YH(L,2) = H*YH(L,2)
  250             CONTINUE
                  INIT = 2
  260          CONTINUE
               !
               !               ......................................................
               !
               !                  ON EACH CALL SET INFORMATION WHICH DETERMINES THE
               !                  ALLOWED INTERVAL OF INTEGRATION BEFORE RETURNING
               !                  WITH AN ANSWER AT TOUT
               !
               DEL = TOUT - T
               ABSDEL = ABS(DEL)
  !
  !               ......................................................
  !
  !                  IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND
  !                  RETURN
  !
  270          CONTINUE
                        !                  BEGIN BLOCK PERMITTING ...EXITS TO 400
                        !                     BEGIN BLOCK PERMITTING ...EXITS TO 380
                        IF (ABS(X-T) .LT. ABSDEL) GO TO 290
                           CALL DINTYD(TOUT,0,YH,NEQ,Y,INTFLG)
                           CALL DINTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
                           IDID = 3
                           IF (X .NE. TOUT) GO TO 280
                              IDID = 2
                              INTOUT = .FALSE.
  280                      CONTINUE
                           T = TOUT
                           TOLD = T
                           !         ..................EXIT
                           GO TO 450
  290                   CONTINUE
                        !
                        !                        IF CANNOT GO PAST TSTOP AND SUFFICIENTLY
                        !                        CLOSE, EXTRAPOLATE AND RETURN
                        !
                        IF (ITSTOP .NE. 1) GO TO 310
                        IF (ABS(TSTOP-X) .GE. 100.0D0*U*ABS(X))&
                           GO TO 310
                           DT = TOUT - X
                           DO 300 L = 1, NEQ
                              Y(L) = YH(L,1) + (DT/H)*YH(L,2)
  300                      CONTINUE
                           CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
                           NFE = NFE + 1
                           IDID = 3
                           T = TOUT
                           TOLD = T
                           !         ..................EXIT
                           GO TO 450
  310                   CONTINUE
                        !
                        IF (IINTEG .EQ. 0 .OR. .NOT.INTOUT) GO TO 320
                           !
                           !                           INTERMEDIATE-OUTPUT MODE
                           !
                           IDID = 1
                        GO TO 370
  320                   CONTINUE
                        !
                        !                        .............................................
                        !
                        !                             MONITOR NUMBER OF STEPS ATTEMPTED
                        !
                        IF (KSTEPS .LE. MAXNUM) GO TO 330
                           !
                           !                           A SIGNIFICANT AMOUNT OF WORK HAS BEEN
                           !                           EXPENDED
                           IDID = -1
                           KSTEPS = 0
                           IBEGIN = -1
                        GO TO 370
  330                   CONTINUE
                           !
                           !                           ..........................................
                           !
                           !                              LIMIT STEP SIZE AND SET WEIGHT VECTOR
                           !
                           HMIN = 100.0D0*U*ABS(X)
                           HA = MAX(ABS(H),HMIN)
                           IF (ITSTOP .EQ. 1)&
                              HA = MIN(HA,ABS(TSTOP-X))
                           H = SIGN(HA,H)
                           LTOL = 1
                           DO 340 L = 1, NEQ
                              IF (ITOL .EQ. 1) LTOL = L
                              EWT(L) = RTOL(LTOL)*ABS(YH(L,1))&
                                       + ATOL(LTOL)
                              !                     .........EXIT
                              IF (EWT(L) .LE. 0.0D0) GO TO 380
  340                      CONTINUE
                           TOLFAC = U*DVNRMS(NEQ,YH,EWT)
                           !                  .........EXIT
                           IF (TOLFAC .LE. 1.0D0) GO TO 400
                           !
                           !                           TOLERANCES TOO SMALL
                           IDID = -2
                           TOLFAC = 2.0D0*TOLFAC
                           RTOL(1) = TOLFAC*RTOL(1)
                           ATOL(1) = TOLFAC*ATOL(1)
                           IF (ITOL .EQ. 0) GO TO 360
                              DO 350 L = 2, NEQ
                                 RTOL(L) = TOLFAC*RTOL(L)
                                 ATOL(L) = TOLFAC*ATOL(L)
  350                         CONTINUE
  360                      CONTINUE
                           IBEGIN = -1
  370                   CONTINUE
                        !            ............EXIT
                        GO TO 430
  380                CONTINUE
  !
  !                     RELATIVE ERROR CRITERION INAPPROPRIATE
  390                CONTINUE
                     IDID = -3
                     IBEGIN = -1
                     !            .........EXIT
                     GO TO 430
  400             CONTINUE
                  !
                  !                  ...................................................
                  !
                  !                       TAKE A STEP
                  !
                  CALL DSTOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM,&
                             DF,DJAC,RPAR,IPAR)
                  !
                  JSTART = -2
                  INTOUT = .TRUE.
               IF (KFLAG .EQ. 0) GO TO 270
               !
               !               ......................................................
               !
               IF (KFLAG .EQ. -1) GO TO 410
                  !
                  !                  REPEATED CORRECTOR CONVERGENCE FAILURES
                  IDID = -6
                  IBEGIN = -1
               GO TO 420
  410          CONTINUE
                  !
                  !                  REPEATED ERROR TEST FAILURES
                  IDID = -7
                  IBEGIN = -1
  420          CONTINUE
  430       CONTINUE
            !
            !            .........................................................
            !
            !                                   STORE VALUES BEFORE RETURNING TO
            !                                   DDEBDF
            DO 440 L = 1, NEQ
               Y(L) = YH(L,1)
               YPOUT(L) = YH(L,2)/H
  440       CONTINUE
            T = X
            TOLD = T
            INTOUT = .FALSE.
  450    CONTINUE
  460 CONTINUE
      RETURN
      END
      ! DECK DSTOD
      SUBROUTINE DSTOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM,&
         DF, DJAC, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DSTOD
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (STOD-S, DSTOD-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DSTOD integrates a system of first order odes over one step in the
      !    integrator package DDEBDF.
      !  ----------------------------------------------------------------------
      !  DSTOD  performs one step of the integration of an initial value
      !  problem for a system of ordinary differential equations.
      !  Note.. DSTOD  is independent of the value of the iteration method
      !  indicator MITER, when this is .NE. 0, and hence is independent
      !  of the type of chord method used, or the Jacobian structure.
      !  Communication with DSTOD  is done with the following variables..
      !
      !  Y      = An array of length .GE. N used as the Y argument in
      !           all calls to DF and DJAC.
      !  NEQ    = Integer array containing problem size in NEQ(1), and
      !           passed as the NEQ argument in all calls to DF and DJAC.
      !  YH     = An NYH by LMAX array containing the dependent variables
      !           and their approximate scaled derivatives, where
      !           LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
      !           J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
      !           (J = 0,1,...,NQ).  On entry for the first step, the first
      !           two columns of YH must be set from the initial values.
      !  NYH    = A constant integer .GE. N, the first dimension of YH.
      !  YH1    = A one-dimensional array occupying the same space as YH.
      !  EWT    = An array of N elements with which the estimated local
      !           errors in YH are compared.
      !  SAVF   = An array of working storage, of length N.
      !  ACOR   = A work array of length N, used for the accumulated
      !           corrections.  On a successful return, ACOR(I) contains
      !           the estimated one-step local error in Y(I).
      !  WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
      !           matrix operations in chord iteration (MITER .NE. 0).
      !  DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
      !           if a chord method is being used.
      !  DSLVS   = Name of routine to solve linear system in chord iteration.
      !  H      = The step size to be attempted on the next step.
      !           H is altered by the error control algorithm during the
      !           problem.  H can be either positive or negative, but its
      !           sign must remain constant throughout the problem.
      !  HMIN   = The minimum absolute value of the step size H to be used.
      !  HMXI   = Inverse of the maximum absolute value of H to be used.
      !           HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
      !           HMIN and HMXI may be changed at any time, but will not
      !           take effect until the next change of H is considered.
      !  TN     = The independent variable. TN is updated on each step taken.
      !  JSTART = An integer used for input only, with the following
      !           values and meanings..
      !                0  Perform the first step.
      !            .GT.0  Take a new step continuing from the last.
      !               -1  Take the next step with a new value of H, MAXORD,
      !                     N, METH, MITER, and/or matrix parameters.
      !               -2  Take the next step with a new value of H,
      !                     but with other inputs unchanged.
      !           On return, JSTART is set to 1 to facilitate continuation.
      !  KFLAG  = a completion code with the following meanings..
      !                0  The step was successful.
      !               -1  The requested error could not be achieved.
      !               -2  Corrector convergence could not be achieved.
      !           A return with KFLAG = -1 or -2 means either
      !           ABS(H) = HMIN or 10 consecutive failures occurred.
      !           On a return with KFLAG negative, the values of TN and
      !           the YH array are as of the beginning of the last
      !           step, and H is the last step size attempted.
      !  MAXORD = The maximum order of integration method to be allowed.
      !  METH/MITER = The method flags.  See description in driver.
      !  N      = The number of first-order differential equations.
      !  ----------------------------------------------------------------------
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  DCFOD, DPJAC, DSLVS, DVNRMS
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      !    920422  Changed DIMENSION statement.  (WRB)
      ! ***END PROLOGUE  DSTOD
      !
      INTEGER I, I1, IALTH, IER, IOD, IOWND, IPAR, IPUP, IREDO, IRET,&
            IWM, J, JB, JSTART, KFLAG, KSTEPS, L, LMAX, M, MAXORD,&
            MEO, METH, MITER, N, NCF, NEQ, NEWQ, NFE, NJE, NQ, NQNYH,&
            NQU, NST, NSTEPJ, NYH
      DOUBLE PRECISION ACOR, CONIT, CRATE, DCON, DDN,&
            DEL, DELP, DSM, DUP, DVNRMS, EL, EL0, ELCO,&
            EWT, EXDN, EXSM, EXUP, H, HMIN, HMXI, HOLD, HU, R, RC,&
            RH, RHDN, RHSM, RHUP, RMAX, ROWND, RPAR, SAVF, TESCO,&
            TN, TOLD, UROUND, WM, Y, YH, YH1
      EXTERNAL DF, DJAC
      !
      DIMENSION Y(*),YH(NYH,*),YH1(*),EWT(*),SAVF(*),ACOR(*),WM(*),&
                IWM(*),RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,CONIT,CRATE,EL(13),ELCO(13,12),HOLD,RC,RMAX,&
                      TESCO(3,12),EL0,H,HMIN,HMXI,HU,TN,UROUND,IOWND(7),&
                      KSTEPS,IOD(6),IALTH,IPUP,LMAX,MEO,NQNYH,NSTEPJ,&
                      IER,JSTART,KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE,&
                      NJE,NQU
            !
            !
            !      BEGIN BLOCK PERMITTING ...EXITS TO 690
            !         BEGIN BLOCK PERMITTING ...EXITS TO 60
            ! ***FIRST EXECUTABLE STATEMENT  DSTOD
            KFLAG = 0
            TOLD = TN
            NCF = 0
            IF (JSTART .GT. 0) GO TO 160
            IF (JSTART .EQ. -1) GO TO 10
               IF (JSTART .EQ. -2) GO TO 90
               !               ---------------------------------------------------------
               !                ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER
               !                VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY
               !                WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS
               !                INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H,
               !                BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE OCCURS
               !                (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT
               !                2 FOR THE NEXT INCREASE.
               !               ---------------------------------------------------------
               LMAX = MAXORD + 1
               NQ = 1
               L = 2
               IALTH = 2
               RMAX = 10000.0D0
               RC = 0.0D0
               EL0 = 1.0D0
               CRATE = 0.7D0
               DELP = 0.0D0
               HOLD = H
               MEO = METH
               NSTEPJ = 0
               IRET = 3
            GO TO 50
   10       CONTINUE
                  !               BEGIN BLOCK PERMITTING ...EXITS TO 30
                  !                  ------------------------------------------------------
                  !                   THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN
                  !                   JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX
                  !                   UPDATE.  IF AN ORDER INCREASE IS ABOUT TO BE
                  !                   CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO
                  !                   POSTPONE CONSIDERATION ONE MORE STEP.  IF THE CALLER
                  !                   HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE
                  !                   COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS
                  !                   CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
                  !                   ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN
                  !                   ACCORDINGLY.  IF H IS TO BE CHANGED, YH MUST BE
                  !                   RESCALED.  IF H OR METH IS BEING CHANGED, IALTH IS
                  !                   RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H
                  !                   FOR THAT MANY STEPS.
                  !                  ------------------------------------------------------
                  IPUP = MITER
                  LMAX = MAXORD + 1
                  IF (IALTH .EQ. 1) IALTH = 2
                  IF (METH .EQ. MEO) GO TO 20
                     CALL DCFOD(METH,ELCO,TESCO)
                     MEO = METH
                     !               ......EXIT
                     IF (NQ .GT. MAXORD) GO TO 30
                     IALTH = L
                     IRET = 1
                     !         ............EXIT
                     GO TO 60
   20             CONTINUE
                  IF (NQ .LE. MAXORD) GO TO 90
   30          CONTINUE
               NQ = MAXORD
               L = LMAX
               DO 40 I = 1, L
                  EL(I) = ELCO(I,NQ)
   40          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               DDN = DVNRMS(N,SAVF,EWT)/TESCO(1,L)
               EXDN = 1.0D0/L
               RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
               RH = MIN(RHDN,1.0D0)
               IREDO = 3
               IF (H .EQ. HOLD) GO TO 660
               RH = MIN(RH,ABS(H/HOLD))
               H = HOLD
               GO TO 100
   50       CONTINUE
            !            ------------------------------------------------------------
            !             DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS
            !             FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED
            !             CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT
            !             THE START OF THE PROBLEM.
            !            ------------------------------------------------------------
            CALL DCFOD(METH,ELCO,TESCO)
   60    CONTINUE
   70    CONTINUE
               !            BEGIN BLOCK PERMITTING ...EXITS TO 680
               DO 80 I = 1, L
                  EL(I) = ELCO(I,NQ)
   80          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               GO TO (90,660,160), IRET
   !               ---------------------------------------------------------
   !                IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
   !                RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH
   !                IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT
   !                MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST
   !                FAILURE.
   !               ---------------------------------------------------------
   90          CONTINUE
               IF (H .EQ. HOLD) GO TO 160
               RH = H/HOLD
               H = HOLD
               IREDO = 3
  100          CONTINUE
  110          CONTINUE
                  RH = MIN(RH,RMAX)
                  RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
                  R = 1.0D0
                  DO 130 J = 2, L
                     R = R*RH
                     DO 120 I = 1, N
                        YH(I,J) = YH(I,J)*R
  120                CONTINUE
  130             CONTINUE
                  H = H*RH
                  RC = RC*RH
                  IALTH = L
                  IF (IREDO .NE. 0) GO TO 150
                     RMAX = 10.0D0
                     R = 1.0D0/TESCO(2,NQU)
                     DO 140 I = 1, N
                        ACOR(I) = ACOR(I)*R
  140                CONTINUE
                     !      ...............EXIT
                     GO TO 690
  150             CONTINUE
  !                  ------------------------------------------------------
  !                   THIS SECTION COMPUTES THE PREDICTED VALUES BY
  !                   EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL
  !                   TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD
  !                   VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS
  !                   FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
  !                   TO FORCE DPJAC TO BE CALLED, IF A JACOBIAN IS
  !                   INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST
  !                   EVERY 20-TH STEP.
  !                  ------------------------------------------------------
  160             CONTINUE
  170             CONTINUE
                           !                     BEGIN BLOCK PERMITTING ...EXITS TO 610
                           !                        BEGIN BLOCK PERMITTING ...EXITS TO 490
                           IF (ABS(RC-1.0D0) .GT. 0.3D0) IPUP = MITER
                           IF (NST .GE. NSTEPJ + 20) IPUP = MITER
                           TN = TN + H
                           I1 = NQNYH + 1
                           DO 190 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 180 I = I1, NQNYH
                                 YH1(I) = YH1(I) + YH1(I+NYH)
  180                         CONTINUE
  190                      CONTINUE
                           KSTEPS = KSTEPS + 1
  !                           ---------------------------------------------
  !                            UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A
  !                            CONVERGENCE TEST IS MADE ON THE R.M.S. NORM
  !                            OF EACH CORRECTION, WEIGHTED BY THE ERROR
  !                            WEIGHT VECTOR EWT.  THE SUM OF THE
  !                            CORRECTIONS IS ACCUMULATED IN THE VECTOR
  !                            ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE
  !                            CORRECTOR LOOP.
  !                           ---------------------------------------------
  200                      CONTINUE
                              M = 0
                              DO 210 I = 1, N
                                 Y(I) = YH(I,1)
  210                         CONTINUE
                              CALL DF(TN,Y,SAVF,RPAR,IPAR)
                              NFE = NFE + 1
                              IF (IPUP .LE. 0) GO TO 220
                                 !                                 ---------------------------------------
                                 !                                  IF INDICATED, THE MATRIX P = I -
                                 !                                  H*EL(1)*J IS REEVALUATED AND
                                 !                                  PREPROCESSED BEFORE STARTING THE
                                 !                                  CORRECTOR ITERATION.  IPUP IS SET TO 0
                                 !                                  AS AN INDICATOR THAT THIS HAS BEEN
                                 !                                  DONE.
                                 !                                 ---------------------------------------
                                 IPUP = 0
                                 RC = 1.0D0
                                 NSTEPJ = NST
                                 CRATE = 0.7D0
                                 CALL DPJAC(NEQ,Y,YH,NYH,EWT,ACOR,SAVF,&
                                            WM,IWM,DF,DJAC,RPAR,IPAR)
                                 !                           ......EXIT
                                 IF (IER .NE. 0) GO TO 440
  220                         CONTINUE
                              DO 230 I = 1, N
                                 ACOR(I) = 0.0D0
  230                         CONTINUE
  240                         CONTINUE
                                 IF (MITER .NE. 0) GO TO 270
                                    !                                    ------------------------------------
                                    !                                     IN THE CASE OF FUNCTIONAL
                                    !                                     ITERATION, UPDATE Y DIRECTLY FROM
                                    !                                     THE RESULT OF THE LAST FUNCTION
                                    !                                     EVALUATION.
                                    !                                    ------------------------------------
                                    DO 250 I = 1, N
                                       SAVF(I) = H*SAVF(I) - YH(I,2)
                                       Y(I) = SAVF(I) - ACOR(I)
  250                               CONTINUE
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 260 I = 1, N
                                       Y(I) = YH(I,1) + EL(1)*SAVF(I)
                                       ACOR(I) = SAVF(I)
  260                               CONTINUE
                                 GO TO 300
  270                            CONTINUE
                                    !                                    ------------------------------------
                                    !                                     IN THE CASE OF THE CHORD METHOD,
                                    !                                     COMPUTE THE CORRECTOR ERROR, AND
                                    !                                     SOLVE THE LINEAR SYSTEM WITH THAT
                                    !                                     AS RIGHT-HAND SIDE AND P AS
                                    !                                     COEFFICIENT MATRIX.
                                    !                                    ------------------------------------
                                    DO 280 I = 1, N
                                       Y(I) = H*SAVF(I)&
                                              - (YH(I,2) + ACOR(I))
  280                               CONTINUE
                                    CALL DSLVS(WM,IWM,Y,SAVF)
                                    !                              ......EXIT
                                    IF (IER .NE. 0) GO TO 430
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 290 I = 1, N
                                       ACOR(I) = ACOR(I) + Y(I)
                                       Y(I) = YH(I,1) + EL(1)*ACOR(I)
  290                               CONTINUE
  300                            CONTINUE
                                 !                                 ---------------------------------------
                                 !                                  TEST FOR CONVERGENCE.  IF M.GT.0, AN
                                 !                                  ESTIMATE OF THE CONVERGENCE RATE
                                 !                                  CONSTANT IS STORED IN CRATE, AND THIS
                                 !                                  IS USED IN THE TEST.
                                 !                                 ---------------------------------------
                                 IF (M .NE. 0)&
                                    CRATE = MAX(0.2D0*CRATE,DEL/DELP)
                                 DCON = DEL*MIN(1.0D0,1.5D0*CRATE)&
                                        /(TESCO(2,NQ)*CONIT)
                                 IF (DCON .GT. 1.0D0) GO TO 420
                                    !                                    ------------------------------------
                                    !                                     THE CORRECTOR HAS CONVERGED.  IPUP
                                    !                                     IS SET TO -1 IF MITER .NE. 0, TO
                                    !                                     SIGNAL THAT THE JACOBIAN INVOLVED
                                    !                                     MAY NEED UPDATING LATER.  THE LOCAL
                                    !                                     ERROR TEST IS MADE AND CONTROL
                                    !                                     PASSES TO STATEMENT 500 IF IT
                                    !                                     FAILS.
                                    !                                    ------------------------------------
                                    IF (MITER .NE. 0) IPUP = -1
                                    IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
                                    IF (M .GT. 0)&
                                       DSM = DVNRMS(N,ACOR,EWT)&
                                             /TESCO(2,NQ)
                                    IF (DSM .GT. 1.0D0) GO TO 380
                                          !                                       BEGIN BLOCK
                                          !                                       PERMITTING ...EXITS TO 360
                                          !                                          ------------------------------
                                          !                                           AFTER A SUCCESSFUL STEP,
                                          !                                           UPDATE THE YH ARRAY.
                                          !                                           CONSIDER CHANGING H IF IALTH
                                          !                                           = 1.  OTHERWISE DECREASE
                                          !                                           IALTH BY 1.  IF IALTH IS THEN
                                          !                                           1 AND NQ .LT. MAXORD, THEN
                                          !                                           ACOR IS SAVED FOR USE IN A
                                          !                                           POSSIBLE ORDER INCREASE ON
                                          !                                           THE NEXT STEP.  IF A CHANGE
                                          !                                           IN H IS CONSIDERED, AN
                                          !                                           INCREASE OR DECREASE IN ORDER
                                          !                                           BY ONE IS CONSIDERED ALSO.  A
                                          !                                           CHANGE IN H IS MADE ONLY IF
                                          !                                           IT IS BY A FACTOR OF AT LEAST
                                          !                                           1.1.  IF NOT, IALTH IS SET TO
                                          !                                           3 TO PREVENT TESTING FOR THAT
                                          !                                           MANY STEPS.
                                          !                                          ------------------------------
                                          KFLAG = 0
                                          IREDO = 0
                                          NST = NST + 1
                                          HU = H
                                          NQU = NQ
                                          DO 320 J = 1, L
                                             DO 310 I = 1, N
                                                YH(I,J) = YH(I,J)&
                                                          + EL(J)&
                                                            *ACOR(I)
  310                                        CONTINUE
  320                                     CONTINUE
                                          IALTH = IALTH - 1
                                          IF (IALTH .NE. 0) GO TO 340
                                             !                                             ---------------------------
                                             !                                              REGARDLESS OF THE SUCCESS
                                             !                                              OR FAILURE OF THE STEP,
                                             !                                              FACTORS RHDN, RHSM, AND
                                             !                                              RHUP ARE COMPUTED, BY
                                             !                                              WHICH H COULD BE
                                             !                                              MULTIPLIED AT ORDER NQ -
                                             !                                              1, ORDER NQ, OR ORDER NQ +
                                             !                                              1, RESPECTIVELY.  IN THE
                                             !                                              CASE OF FAILURE, RHUP =
                                             !                                              0.0 TO AVOID AN ORDER
                                             !                                              INCREASE.  THE LARGEST OF
                                             !                                              THESE IS DETERMINED AND
                                             !                                              THE NEW ORDER CHOSEN
                                             !                                              ACCORDINGLY.  IF THE ORDER
                                             !                                              IS TO BE INCREASED, WE
                                             !                                              COMPUTE ONE ADDITIONAL
                                             !                                              SCALED DERIVATIVE.
                                             !                                             ---------------------------
                                             RHUP = 0.0D0
                                             !                        .....................EXIT
                                             IF (L .EQ. LMAX) GO TO 490
                                             DO 330 I = 1, N
                                                SAVF(I) = ACOR(I)&
                                                          - YH(I,LMAX)
  330                                        CONTINUE
                                             DUP = DVNRMS(N,SAVF,EWT)&
                                                   /TESCO(3,NQ)
                                             EXUP = 1.0D0/(L+1)
                                             RHUP = 1.0D0&
                                                    /(1.4D0*DUP**EXUP&
                                                      + 0.0000014D0)
                                             !                        .....................EXIT
                                             GO TO 490
  340                                     CONTINUE
                                          !                                       ...EXIT
                                          IF (IALTH .GT. 1) GO TO 360
                                          !                                       ...EXIT
                                          IF (L .EQ. LMAX) GO TO 360
                                          DO 350 I = 1, N
                                             YH(I,LMAX) = ACOR(I)
  350                                     CONTINUE
  360                                  CONTINUE
                                       R = 1.0D0/TESCO(2,NQU)
                                       DO 370 I = 1, N
                                          ACOR(I) = ACOR(I)*R
  370                                  CONTINUE
                                       !      .................................EXIT
                                       GO TO 690
  380                               CONTINUE
                                    !                                    ------------------------------------
                                    !                                     THE ERROR TEST FAILED.  KFLAG KEEPS
                                    !                                     TRACK OF MULTIPLE FAILURES.
                                    !                                     RESTORE TN AND THE YH ARRAY TO
                                    !                                     THEIR PREVIOUS VALUES, AND PREPARE
                                    !                                     TO TRY THE STEP AGAIN.  COMPUTE THE
                                    !                                     OPTIMUM STEP SIZE FOR THIS OR ONE
                                    !                                     LOWER ORDER.  AFTER 2 OR MORE
                                    !                                     FAILURES, H IS FORCED TO DECREASE
                                    !                                     BY A FACTOR OF 0.2 OR LESS.
                                    !                                    ------------------------------------
                                    KFLAG = KFLAG - 1
                                    TN = TOLD
                                    I1 = NQNYH + 1
                                    DO 400 JB = 1, NQ
                                       I1 = I1 - NYH
                                       DO 390 I = I1, NQNYH
                                          YH1(I) = YH1(I) - YH1(I+NYH)
  390                                  CONTINUE
  400                               CONTINUE
                                    RMAX = 2.0D0
                                    IF (ABS(H) .GT. HMIN*1.00001D0)&
                                       GO TO 410
                                       !                                       ---------------------------------
                                       !                                        ALL RETURNS ARE MADE THROUGH
                                       !                                        THIS SECTION.  H IS SAVED IN
                                       !                                        HOLD TO ALLOW THE CALLER TO
                                       !                                        CHANGE H ON THE NEXT STEP.
                                       !                                       ---------------------------------
                                       KFLAG = -1
                                       !      .................................EXIT
                                       GO TO 690
  410                               CONTINUE
                                    !                     ...............EXIT
                                    IF (KFLAG .LE. -3) GO TO 610
                                    IREDO = 2
                                    RHUP = 0.0D0
                                    !                        ............EXIT
                                    GO TO 490
  420                            CONTINUE
                                 M = M + 1
                                 !                              ...EXIT
                                 IF (M .EQ. 3) GO TO 430
                                 !                              ...EXIT
                                 IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP)&
                                    GO TO 430
                                 DELP = DEL
                                 CALL DF(TN,Y,SAVF,RPAR,IPAR)
                                 NFE = NFE + 1
                              GO TO 240
  430                         CONTINUE
                              !                              ------------------------------------------
                              !                               THE CORRECTOR ITERATION FAILED TO
                              !                               CONVERGE IN 3 TRIES.  IF MITER .NE. 0 AND
                              !                               THE JACOBIAN IS OUT OF DATE, DPJAC IS
                              !                               CALLED FOR THE NEXT TRY.  OTHERWISE THE
                              !                               YH ARRAY IS RETRACTED TO ITS VALUES
                              !                               BEFORE PREDICTION, AND H IS REDUCED, IF
                              !                               POSSIBLE.  IF H CANNOT BE REDUCED OR 10
                              !                               FAILURES HAVE OCCURRED, EXIT WITH KFLAG =
                              !                               -2.
                              !                              ------------------------------------------
                              !                           ...EXIT
                              IF (IPUP .EQ. 0) GO TO 440
                              IPUP = MITER
                           GO TO 200
  440                      CONTINUE
                           TN = TOLD
                           NCF = NCF + 1
                           RMAX = 2.0D0
                           I1 = NQNYH + 1
                           DO 460 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 450 I = I1, NQNYH
                                 YH1(I) = YH1(I) - YH1(I+NYH)
  450                         CONTINUE
  460                      CONTINUE
                           IF (ABS(H) .GT. HMIN*1.00001D0) GO TO 470
                              KFLAG = -2
                              !      ........................EXIT
                              GO TO 690
  470                      CONTINUE
                           IF (NCF .NE. 10) GO TO 480
                              KFLAG = -2
                              !      ........................EXIT
                              GO TO 690
  480                      CONTINUE
                           RH = 0.25D0
                           IPUP = MITER
                           IREDO = 1
                           !                  .........EXIT
                           GO TO 650
  490                   CONTINUE
                        EXSM = 1.0D0/L
                        RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
                        RHDN = 0.0D0
                        IF (NQ .EQ. 1) GO TO 500
                           DDN = DVNRMS(N,YH(1,L),EWT)/TESCO(1,NQ)
                           EXDN = 1.0D0/NQ
                           RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
  500                   CONTINUE
                        IF (RHSM .GE. RHUP) GO TO 550
                           IF (RHUP .LE. RHDN) GO TO 540
                              NEWQ = L
                              RH = RHUP
                              IF (RH .GE. 1.1D0) GO TO 520
                                 IALTH = 3
                                 R = 1.0D0/TESCO(2,NQU)
                                 DO 510 I = 1, N
                                    ACOR(I) = ACOR(I)*R
  510                            CONTINUE
                                 !      ...........................EXIT
                                 GO TO 690
  520                         CONTINUE
                              R = EL(L)/L
                              DO 530 I = 1, N
                                 YH(I,NEWQ+1) = ACOR(I)*R
  530                         CONTINUE
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
                              !            ..................EXIT
                              GO TO 680
  540                      CONTINUE
                        GO TO 580
  550                   CONTINUE
                        IF (RHSM .LT. RHDN) GO TO 580
                           NEWQ = NQ
                           RH = RHSM
                           IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0)&
                              GO TO 560
                              IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
                              !                              ------------------------------------------
                              !                               IF THERE IS A CHANGE OF ORDER, RESET NQ,
                              !                               L, AND THE COEFFICIENTS.  IN ANY CASE H
                              !                               IS RESET ACCORDING TO RH AND THE YH ARRAY
                              !                               IS RESCALED.  THEN EXIT FROM 680 IF THE
                              !                               STEP WAS OK, OR REDO THE STEP OTHERWISE.
                              !                              ------------------------------------------
                              !                  ............EXIT
                              IF (NEWQ .EQ. NQ) GO TO 650
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
                              !            ..................EXIT
                              GO TO 680
  560                      CONTINUE
                           IALTH = 3
                           R = 1.0D0/TESCO(2,NQU)
                           DO 570 I = 1, N
                              ACOR(I) = ACOR(I)*R
  570                      CONTINUE
                           !      .....................EXIT
                           GO TO 690
  580                   CONTINUE
                        NEWQ = NQ - 1
                        RH = RHDN
                        IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
                        IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 590
                           IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
                           !                           ---------------------------------------------
                           !                            IF THERE IS A CHANGE OF ORDER, RESET NQ, L,
                           !                            AND THE COEFFICIENTS.  IN ANY CASE H IS
                           !                            RESET ACCORDING TO RH AND THE YH ARRAY IS
                           !                            RESCALED.  THEN EXIT FROM 680 IF THE STEP
                           !                            WAS OK, OR REDO THE STEP OTHERWISE.
                           !                           ---------------------------------------------
                           !                  .........EXIT
                           IF (NEWQ .EQ. NQ) GO TO 650
                           NQ = NEWQ
                           L = NQ + 1
                           IRET = 2
                           !            ...............EXIT
                           GO TO 680
  590                   CONTINUE
                        IALTH = 3
                        R = 1.0D0/TESCO(2,NQU)
                        DO 600 I = 1, N
                           ACOR(I) = ACOR(I)*R
  600                   CONTINUE
                        !      ..................EXIT
                        GO TO 690
  610                CONTINUE
                     !                     ---------------------------------------------------
                     !                      CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES
                     !                      HAVE OCCURRED.  IF 10 FAILURES HAVE OCCURRED, EXIT
                     !                      WITH KFLAG = -1.  IT IS ASSUMED THAT THE
                     !                      DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY
                     !                      HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
                     !                      DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO
                     !                      1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE
                     !                      STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES
                     !                      HMIN.
                     !                     ---------------------------------------------------
                     IF (KFLAG .NE. -10) GO TO 620
                        !                        ------------------------------------------------
                        !                         ALL RETURNS ARE MADE THROUGH THIS SECTION.  H
                        !                         IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE
                        !                         H ON THE NEXT STEP.
                        !                        ------------------------------------------------
                        KFLAG = -1
                        !      ..................EXIT
                        GO TO 690
  620                CONTINUE
                     RH = 0.1D0
                     RH = MAX(HMIN/ABS(H),RH)
                     H = H*RH
                     DO 630 I = 1, N
                        Y(I) = YH(I,1)
  630                CONTINUE
                     CALL DF(TN,Y,SAVF,RPAR,IPAR)
                     NFE = NFE + 1
                     DO 640 I = 1, N
                        YH(I,2) = H*SAVF(I)
  640                CONTINUE
                     IPUP = MITER
                     IALTH = 5
                     !               ......EXIT
                     IF (NQ .NE. 1) GO TO 670
                  GO TO 170
  650             CONTINUE
  660             CONTINUE
                  RH = MAX(RH,HMIN/ABS(H))
               GO TO 110
  670          CONTINUE
               NQ = 1
               L = 2
               IRET = 3
  680       CONTINUE
         GO TO 70
  690 CONTINUE
      HOLD = H
      JSTART = 1
      RETURN
      !      ----------------------- END OF SUBROUTINE DSTOD
      !      -----------------------
      END
      ! DECK DCFOD
      SUBROUTINE DCFOD (METH, ELCO, TESCO)
      ! ***BEGIN PROLOGUE  DCFOD
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (CFOD-S, DCFOD-D)
      ! ***AUTHOR  (UNKNOWN)
      ! ***DESCRIPTION
      !
      !    DCFOD defines coefficients needed in the integrator package DDEBDF
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      ! ***END PROLOGUE  DCFOD
      !
      !
      INTEGER I, IB, METH, NQ, NQM1, NQP1
      DOUBLE PRECISION AGAMQ, ELCO, FNQ, FNQM1, PC, PINT, RAGQ,&
            RQ1FAC, RQFAC, TESCO, TSIGN, XPIN
      DIMENSION ELCO(13,12),TESCO(3,12)
      !      ------------------------------------------------------------------
      !       DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
      !       NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
      !       GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
      !       THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH =
      !       2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
      !       DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
      !       AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
      !
      !       THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
      !       THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF
      !       ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A
      !       GENERATING POLYNOMIAL, I.E.,
      !           L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
      !       FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
      !           DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) =
      !       0.  FOR THE BDF METHODS, L(X) IS GIVEN BY
      !           L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
      !       WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
      !
      !       THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
      !       LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
      !       AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
      !       SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
      !       NQ + 1 IF K = 3.
      !      ------------------------------------------------------------------
      DIMENSION PC(12)
      !
      ! ***FIRST EXECUTABLE STATEMENT  DCFOD
      GO TO (10,60), METH
   !
   10 CONTINUE
         ELCO(1,1) = 1.0D0
         ELCO(2,1) = 1.0D0
         TESCO(1,1) = 0.0D0
         TESCO(2,1) = 2.0D0
         TESCO(1,2) = 1.0D0
         TESCO(3,12) = 0.0D0
         PC(1) = 1.0D0
         RQFAC = 1.0D0
         DO 50 NQ = 2, 12
            !            ------------------------------------------------------------
            !             THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
            !                 POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1).
            !             INITIALLY, P(X) = 1.
            !            ------------------------------------------------------------
            RQ1FAC = RQFAC
            RQFAC = RQFAC/NQ
            NQM1 = NQ - 1
            FNQM1 = NQM1
            NQP1 = NQ + 1
            !            FORM COEFFICIENTS OF P(X)*(X+NQ-1).
            !            ----------------------------------
            PC(NQ) = 0.0D0
            DO 20 IB = 1, NQM1
               I = NQP1 - IB
               PC(I) = PC(I-1) + FNQM1*PC(I)
   20       CONTINUE
            PC(1) = FNQM1*PC(1)
            !            COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X).
            !            -----------------------
            PINT = PC(1)
            XPIN = PC(1)/2.0D0
            TSIGN = 1.0D0
            DO 30 I = 2, NQ
               TSIGN = -TSIGN
               PINT = PINT + TSIGN*PC(I)/I
               XPIN = XPIN + TSIGN*PC(I)/(I+1)
   30       CONTINUE
            !            STORE COEFFICIENTS IN ELCO AND TESCO.
            !            --------------------------------
            ELCO(1,NQ) = PINT*RQ1FAC
            ELCO(2,NQ) = 1.0D0
            DO 40 I = 2, NQ
               ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
   40       CONTINUE
            AGAMQ = RQFAC*XPIN
            RAGQ = 1.0D0/AGAMQ
            TESCO(2,NQ) = RAGQ
            IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
            TESCO(3,NQM1) = RAGQ
   50    CONTINUE
      GO TO 100
   !
   60 CONTINUE
         PC(1) = 1.0D0
         RQ1FAC = 1.0D0
         DO 90 NQ = 1, 5
            !            ------------------------------------------------------------
            !             THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
            !                 POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ).
            !             INITIALLY, P(X) = 1.
            !            ------------------------------------------------------------
            FNQ = NQ
            NQP1 = NQ + 1
            !            FORM COEFFICIENTS OF P(X)*(X+NQ).
            !            ------------------------------------
            PC(NQP1) = 0.0D0
            DO 70 IB = 1, NQ
               I = NQ + 2 - IB
               PC(I) = PC(I-1) + FNQ*PC(I)
   70       CONTINUE
            PC(1) = FNQ*PC(1)
            !            STORE COEFFICIENTS IN ELCO AND TESCO.
            !            --------------------------------
            DO 80 I = 1, NQP1
               ELCO(I,NQ) = PC(I)/PC(2)
   80       CONTINUE
            ELCO(2,NQ) = 1.0D0
            TESCO(1,NQ) = RQ1FAC
            TESCO(2,NQ) = NQP1/ELCO(1,NQ)
            TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
            RQ1FAC = RQ1FAC/FNQ
   90    CONTINUE
  100 CONTINUE
      RETURN
      !      ----------------------- END OF SUBROUTINE DCFOD
      !      -----------------------
      END
      ! DECK DVNRMS
      DOUBLE PRECISION FUNCTION DVNRMS (N, V, W)
      ! ***BEGIN PROLOGUE  DVNRMS
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D)
      ! ***AUTHOR  (UNKNOWN)
      ! ***DESCRIPTION
      !
      !    DVNRMS computes a weighted root-mean-square vector norm for the
      !    integrator package DDEBDF.
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  (NONE)
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      ! ***END PROLOGUE  DVNRMS
      INTEGER I, N
      DOUBLE PRECISION SUM, V, W
      DIMENSION V(*),W(*)
      ! ***FIRST EXECUTABLE STATEMENT  DVNRMS
      SUM = 0.0D0
      DO 10 I = 1, N
         SUM = SUM + (V(I)/W(I))**2
   10 CONTINUE
      DVNRMS = SQRT(SUM/N)
      RETURN
      !      ----------------------- END OF FUNCTION DVNRMS
      !      ------------------------
      END
      ! DECK DINTYD
      SUBROUTINE DINTYD (T, K, YH, NYH, DKY, IFLAG)
      ! ***BEGIN PROLOGUE  DINTYD
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DINTYD approximates the solution and derivatives at T by polynomial
      !    interpolation. Must be used in conjunction with the integrator
      !    package DDEBDF.
      !  ----------------------------------------------------------------------
      !  DINTYD computes interpolated values of the K-th derivative of the
      !  dependent variable vector Y, and stores it in DKY.
      !  This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may
      !  also be called by the user for any K up to the current order.
      !  (see detailed instructions in LSODE usage documentation.)
      !  ----------------------------------------------------------------------
      !  The computed values in DKY are gotten by interpolation using the
      !  Nordsieck history array YH.  This array corresponds uniquely to a
      !  vector-valued polynomial of degree NQCUR or less, and DKY is set
      !  to the K-th derivative of this polynomial at T.
      !  The formula for DKY is..
      !               Q
      !   DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
      !              J=K
      !  where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
      !  The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
      !  communicated by common.  The above sum is done in reverse order.
      !  IFLAG is returned negative if either K or T is out of bounds.
      !  ----------------------------------------------------------------------
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  (NONE)
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      ! ***END PROLOGUE  DINTYD
      !
      INTEGER I, IC, IER, IFLAG, IOWND, IOWNS, J, JB, JB2, JJ, JJ1,&
            JP1, JSTART, K, KFLAG, L, MAXORD, METH, MITER, N, NFE,&
            NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION C, DKY, EL0, H, HMIN, HMXI, HU, R, ROWND,&
            ROWNS, S, T, TN, TP, UROUND, YH
      DIMENSION YH(NYH,*),DKY(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,&
                      IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,&
                      MAXORD,N,NQ,NST,NFE,NJE,NQU
         !
         !      BEGIN BLOCK PERMITTING ...EXITS TO 130
         ! ***FIRST EXECUTABLE STATEMENT  DINTYD
         IFLAG = 0
         IF (K .LT. 0 .OR. K .GT. NQ) GO TO 110
            TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
            IF ((T - TP)*(T - TN) .LE. 0.0D0) GO TO 10
               IFLAG = -2
               !      .........EXIT
               GO TO 130
   10       CONTINUE
            !
            S = (T - TN)/H
            IC = 1
            IF (K .EQ. 0) GO TO 30
               JJ1 = L - K
               DO 20 JJ = JJ1, NQ
                  IC = IC*JJ
   20          CONTINUE
   30       CONTINUE
            C = IC
            DO 40 I = 1, N
               DKY(I) = C*YH(I,L)
   40       CONTINUE
            IF (K .EQ. NQ) GO TO 90
               JB2 = NQ - K
               DO 80 JB = 1, JB2
                  J = NQ - JB
                  JP1 = J + 1
                  IC = 1
                  IF (K .EQ. 0) GO TO 60
                     JJ1 = JP1 - K
                     DO 50 JJ = JJ1, J
                        IC = IC*JJ
   50                CONTINUE
   60             CONTINUE
                  C = IC
                  DO 70 I = 1, N
                     DKY(I) = C*YH(I,JP1) + S*DKY(I)
   70             CONTINUE
   80          CONTINUE
               !      .........EXIT
               IF (K .EQ. 0) GO TO 130
   90       CONTINUE
            R = H**(-K)
            DO 100 I = 1, N
               DKY(I) = R*DKY(I)
  100       CONTINUE
         GO TO 120
  110    CONTINUE
            !
            IFLAG = -1
  120    CONTINUE
  130 CONTINUE
      RETURN
      !      ----------------------- END OF SUBROUTINE DINTYD
      !      -----------------------
      END
      ! DECK DPJAC
      SUBROUTINE DPJAC (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, DF,&
         DJAC, RPAR, IPAR)
      ! ***BEGIN PROLOGUE  DPJAC
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (PJAC-S, DPJAC-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DPJAC sets up the iteration matrix (involving the Jacobian) for the
      !    integration package DDEBDF.
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  DGBFA, DGEFA, DVNRMS
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890911  Removed unnecessary intrinsics.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      !    920422  Changed DIMENSION statement.  (WRB)
      ! ***END PROLOGUE  DPJAC
      !
      INTEGER I, I1, I2, IER, II, IOWND, IOWNS, IPAR, IWM, J, J1,&
            JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND,&
            MEB1, MEBAND, METH, MITER, ML, ML3, MU, N, NEQ,&
            NFE, NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION CON, DI, DVNRMS, EL0, EWT,&
            FAC, FTEM, H, HL0, HMIN, HMXI, HU, R, R0, ROWND, ROWNS,&
            RPAR, SAVF, SRUR, TN, UROUND, WM, Y, YH, YI, YJ, YJJ
      EXTERNAL DF, DJAC
      DIMENSION Y(*),YH(NYH,*),EWT(*),FTEM(*),SAVF(*),WM(*),IWM(*),&
                RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,&
                      IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,&
                      MAXORD,N,NQ,NST,NFE,NJE,NQU
                  !      ------------------------------------------------------------------
                  !       DPJAC IS CALLED BY DSTOD  TO COMPUTE AND PROCESS THE MATRIX
                  !       P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
                  !       HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE DJAC IF
                  !       MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5.
                  !       IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
                  !       J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN
                  !       SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
                  !       OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
                  !       BY DGEFA IF MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5.
                  !
                  !       IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
                  !       WITH DPJAC USES THE FOLLOWING..
                  !       Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
                  !       FTEM = WORK ARRAY OF LENGTH N (ACOR IN DSTOD ).
                  !       SAVF = ARRAY CONTAINING DF EVALUATED AT PREDICTED Y.
                  !       WM   = DOUBLE PRECISION WORK SPACE FOR MATRICES.  ON OUTPUT IT
                  !       CONTAINS THE
                  !              INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU
                  !              DECOMPOSITION OF P IF MITER IS 1, 2 , 4, OR 5.
                  !              STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
                  !              WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
                  !              WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN
                  !              INCREMENTS.  WM(2) = H*EL0, SAVED FOR LATER USE IF MITER =
                  !              3.
                  !       IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
                  !              AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
                  !              THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER
                  !              IS 4 OR 5.
                  !       EL0  = EL(1) (INPUT).
                  !       IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF
                  !              P MATRIX FOUND TO BE SINGULAR.
                  !       THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
                  !       MITER, N, NFE, AND NJE.
                  ! -----------------------------------------------------------------------
                  !      BEGIN BLOCK PERMITTING ...EXITS TO 240
                  !         BEGIN BLOCK PERMITTING ...EXITS TO 220
                  !            BEGIN BLOCK PERMITTING ...EXITS TO 130
                  !               BEGIN BLOCK PERMITTING ...EXITS TO 70
                  ! ***FIRST EXECUTABLE STATEMENT  DPJAC
                  NJE = NJE + 1
                  HL0 = H*EL0
                  GO TO (10,40,90,140,170), MITER
   !                  IF MITER = 1, CALL DJAC AND MULTIPLY BY SCALAR.
   !                  -----------------------
   10             CONTINUE
                  LENP = N*N
                  DO 20 I = 1, LENP
                     WM(I+2) = 0.0D0
   20             CONTINUE
                  CALL DJAC(TN,Y,WM(3),N,RPAR,IPAR)
                  CON = -HL0
                  DO 30 I = 1, LENP
                     WM(I+2) = WM(I+2)*CON
   30             CONTINUE
                  !               ...EXIT
                  GO TO 70
   !                  IF MITER = 2, MAKE N CALLS TO DF TO APPROXIMATE J.
   !                  --------------------
   40             CONTINUE
                  FAC = DVNRMS(N,SAVF,EWT)
                  R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
                  IF (R0 .EQ. 0.0D0) R0 = 1.0D0
                  SRUR = WM(1)
                  J1 = 2
                  DO 60 J = 1, N
                     YJ = Y(J)
                     R = MAX(SRUR*ABS(YJ),R0*EWT(J))
                     Y(J) = Y(J) + R
                     FAC = -HL0/R
                     CALL DF(TN,Y,FTEM,RPAR,IPAR)
                     DO 50 I = 1, N
                        WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
   50                CONTINUE
                     Y(J) = YJ
                     J1 = J1 + N
   60             CONTINUE
                  NFE = NFE + N
   70          CONTINUE
               !               ADD IDENTITY MATRIX.
               !               -------------------------------------------------
               J = 3
               DO 80 I = 1, N
                  WM(J) = WM(J) + 1.0D0
                  J = J + (N + 1)
   80          CONTINUE
               !               DO LU DECOMPOSITION ON P.
               !               --------------------------------------------
               CALL DGEFA(WM(3),N,N,IWM(21),IER)
               !      .........EXIT
               GO TO 240
   !               IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND
   !               P. ---------
   90          CONTINUE
               WM(2) = HL0
               IER = 0
               R = EL0*0.1D0
               DO 100 I = 1, N
                  Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
  100          CONTINUE
               CALL DF(TN,Y,WM(3),RPAR,IPAR)
               NFE = NFE + 1
               DO 120 I = 1, N
                  R0 = H*SAVF(I) - YH(I,2)
                  DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
                  WM(I+2) = 1.0D0
                  IF (ABS(R0) .LT. UROUND*EWT(I)) GO TO 110
                     !            .........EXIT
                     IF (ABS(DI) .EQ. 0.0D0) GO TO 130
                     WM(I+2) = 0.1D0*R0/DI
  110             CONTINUE
  120          CONTINUE
               !      .........EXIT
               GO TO 240
  130       CONTINUE
            IER = -1
            !      ......EXIT
            GO TO 240
  !            IF MITER = 4, CALL DJAC AND MULTIPLY BY SCALAR.
  !            -----------------------
  140       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            ML3 = 3
            MBAND = ML + MU + 1
            MEBAND = MBAND + ML
            LENP = MEBAND*N
            DO 150 I = 1, LENP
               WM(I+2) = 0.0D0
  150       CONTINUE
            CALL DJAC(TN,Y,WM(ML3),MEBAND,RPAR,IPAR)
            CON = -HL0
            DO 160 I = 1, LENP
               WM(I+2) = WM(I+2)*CON
  160       CONTINUE
            !         ...EXIT
            GO TO 220
  !            IF MITER = 5, MAKE MBAND CALLS TO DF TO APPROXIMATE J.
  !            ----------------
  170       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            MBAND = ML + MU + 1
            MBA = MIN(MBAND,N)
            MEBAND = MBAND + ML
            MEB1 = MEBAND - 1
            SRUR = WM(1)
            FAC = DVNRMS(N,SAVF,EWT)
            R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
            IF (R0 .EQ. 0.0D0) R0 = 1.0D0
            DO 210 J = 1, MBA
               DO 180 I = J, N, MBAND
                  YI = Y(I)
                  R = MAX(SRUR*ABS(YI),R0*EWT(I))
                  Y(I) = Y(I) + R
  180          CONTINUE
               CALL DF(TN,Y,FTEM,RPAR,IPAR)
               DO 200 JJ = J, N, MBAND
                  Y(JJ) = YH(JJ,1)
                  YJJ = Y(JJ)
                  R = MAX(SRUR*ABS(YJJ),R0*EWT(JJ))
                  FAC = -HL0/R
                  I1 = MAX(JJ-MU,1)
                  I2 = MIN(JJ+ML,N)
                  II = JJ*MEB1 - ML + 2
                  DO 190 I = I1, I2
                     WM(II+I) = (FTEM(I) - SAVF(I))*FAC
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
            NFE = NFE + MBA
  220    CONTINUE
         !         ADD IDENTITY MATRIX.
         !         -------------------------------------------------
         II = MBAND + 2
         DO 230 I = 1, N
            WM(II) = WM(II) + 1.0D0
            II = II + MEBAND
  230    CONTINUE
         !         DO LU DECOMPOSITION OF P.
         !         --------------------------------------------
         CALL DGBFA(WM(3),MEBAND,N,ML,MU,IWM(21),IER)
  240 CONTINUE
      RETURN
      !      ----------------------- END OF SUBROUTINE DPJAC
      !      -----------------------
      END
      ! DECK DSLVS
      SUBROUTINE DSLVS (WM, IWM, X, TEM)
      ! ***BEGIN PROLOGUE  DSLVS
      ! ***SUBSIDIARY
      ! ***PURPOSE  Subsidiary to DDEBDF
      ! ***LIBRARY   SLATEC
      ! ***TYPE      DOUBLE PRECISION (SLVS-S, DSLVS-D)
      ! ***AUTHOR  Watts, H. A., (SNLA)
      ! ***DESCRIPTION
      !
      !    DSLVS solves the linear system in the iteration scheme for the
      !    integrator package DDEBDF.
      !
      ! ***SEE ALSO  DDEBDF
      ! ***ROUTINES CALLED  DGBSL, DGESL
      ! ***COMMON BLOCKS    DDEBD1
      ! ***REVISION HISTORY  (YYMMDD)
      !    820301  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900328  Added TYPE section.  (WRB)
      !    910722  Updated AUTHOR section.  (ALS)
      !    920422  Changed DIMENSION statement.  (WRB)
      ! ***END PROLOGUE  DSLVS
      !
      INTEGER I, IER, IOWND, IOWNS, IWM, JSTART, KFLAG, L, MAXORD,&
            MEBAND, METH, MITER, ML, MU, N, NFE, NJE, NQ, NQU, NST
      DOUBLE PRECISION DI, EL0, H, HL0, HMIN, HMXI, HU, PHL0,&
            R, ROWND, ROWNS, TEM, TN, UROUND, WM, X
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,&
                      IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,&
                      MAXORD,N,NQ,NST,NFE,NJE,NQU
            !      ------------------------------------------------------------------
            !       THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING
            !       FROM A CHORD ITERATION.  IT IS CALLED BY DSTOD  IF MITER .NE. 0.
            !       IF MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS.
            !       IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL
            !       MATRIX, AND THEN COMPUTES THE SOLUTION.
            !       IF MITER IS 4 OR 5, IT CALLS DGBSL.
            !       COMMUNICATION WITH DSLVS USES THE FOLLOWING VARIABLES..
            !       WM  = DOUBLE PRECISION WORK SPACE CONTAINING THE INVERSE DIAGONAL
            !       MATRIX IF MITER
            !             IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE.
            !             STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
            !             WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
            !             WM(1) = SQRT(UROUND) (NOT USED HERE),
            !             WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER =
            !             3.
            !       IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
            !             AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
            !             THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS
            !             4 OR 5.
            !       X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION
            !             VECTOR ON OUTPUT, OF LENGTH N.
            !       TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION.
            !       IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED.
            !             IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3.
            !       THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N.
            ! -----------------------------------------------------------------------
            !      BEGIN BLOCK PERMITTING ...EXITS TO 80
            !         BEGIN BLOCK PERMITTING ...EXITS TO 60
            ! ***FIRST EXECUTABLE STATEMENT  DSLVS
            IER = 0
            GO TO (10,10,20,70,70), MITER
   10       CONTINUE
            CALL DGESL(WM(3),N,N,IWM(21),X,0)
            !      ......EXIT
            GO TO 80
   !
   20       CONTINUE
            PHL0 = WM(2)
            HL0 = H*EL0
            WM(2) = HL0
            IF (HL0 .EQ. PHL0) GO TO 40
               R = HL0/PHL0
               DO 30 I = 1, N
                  DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
                  !         .........EXIT
                  IF (ABS(DI) .EQ. 0.0D0) GO TO 60
                  WM(I+2) = 1.0D0/DI
   30          CONTINUE
   40       CONTINUE
            DO 50 I = 1, N
               X(I) = WM(I+2)*X(I)
   50       CONTINUE
            !      ......EXIT
            GO TO 80
   60    CONTINUE
         IER = -1
         !      ...EXIT
         GO TO 80
   !
   70    CONTINUE
         ML = IWM(1)
         MU = IWM(2)
         MEBAND = 2*ML + MU + 1
         CALL DGBSL(WM(3),MEBAND,N,ML,MU,IWM(21),X,0)
   80 CONTINUE
      RETURN
      !      ----------------------- END OF SUBROUTINE DSLVS
      !      -----------------------
      END
      ! DECK DGBFA
      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
      ! ***BEGIN PROLOGUE  DGBFA
      ! ***PURPOSE  Factor a band matrix using Gaussian elimination.
      ! ***LIBRARY   SLATEC (LINPACK)
      ! ***CATEGORY  D2A2
      ! ***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
      ! ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
      ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
      ! ***DESCRIPTION
      !
      !      DGBFA factors a double precision band matrix by elimination.
      !
      !      DGBFA is usually called by DGBCO, but it can be called
      !      directly with a saving in time if  RCOND  is not needed.
      !
      !      On Entry
      !
      !         ABD     DOUBLE PRECISION(LDA, N)
      !                 contains the matrix in band storage.  The columns
      !                 of the matrix are stored in the columns of  ABD  and
      !                 the diagonals of the matrix are stored in rows
      !                 ML+1 through 2*ML+MU+1 of  ABD .
      !                 See the comments below for details.
      !
      !         LDA     INTEGER
      !                 the leading dimension of the array  ABD .
      !                 LDA must be .GE. 2*ML + MU + 1 .
      !
      !         N       INTEGER
      !                 the order of the original matrix.
      !
      !         ML      INTEGER
      !                 number of diagonals below the main diagonal.
      !                 0 .LE. ML .LT.  N .
      !
      !         MU      INTEGER
      !                 number of diagonals above the main diagonal.
      !                 0 .LE. MU .LT.  N .
      !                 More efficient if  ML .LE. MU .
      !      On Return
      !
      !         ABD     an upper triangular matrix in band storage and
      !                 the multipliers which were used to obtain it.
      !                 The factorization can be written  A = L*U  where
      !                 L  is a product of permutation and unit lower
      !                 triangular matrices and  U  is upper triangular.
      !
      !         IPVT    INTEGER(N)
      !                 an integer vector of pivot indices.
      !
      !         INFO    INTEGER
      !                 = 0  normal value.
      !                 = K  if  U(K,K) .EQ. 0.0 .  This is not an error
      !                      condition for this subroutine, but it does
      !                      indicate that DGBSL will divide by zero if
      !                      called.  Use  RCOND  in DGBCO for a reliable
      !                      indication of singularity.
      !
      !      Band Storage
      !
      !            If  A  is a band matrix, the following program segment
      !            will set up the input.
      !
      !                    ML = (band width below the diagonal)
      !                    MU = (band width above the diagonal)
      !                    M = ML + MU + 1
      !                    DO 20 J = 1, N
      !                       I1 = MAX(1, J-MU)
      !                       I2 = MIN(N, J+ML)
      !                       DO 10 I = I1, I2
      !                          K = I - J + M
      !                          ABD(K,J) = A(I,J)
      !                 10    CONTINUE
      !                 20 CONTINUE
      !
      !            This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
      !            In addition, the first  ML  rows in  ABD  are used for
      !            elements generated during the triangularization.
      !            The total number of rows needed in  ABD  is  2*ML+MU+1 .
      !            The  ML+MU by ML+MU  upper left triangle and the
      !            ML by ML  lower right triangle are not referenced.
      !
      ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
      !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
      ! ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
      ! ***REVISION HISTORY  (YYMMDD)
      !    780814  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900326  Removed duplicate information from DESCRIPTION section.
      !            (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
      !
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
      !
      ! ***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
      !
      !      ZERO INITIAL FILL-IN COLUMNS
      !
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
      !
      !      GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      !
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
         !
         !         ZERO NEXT FILL-IN COLUMN
         !
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
         !
         !         FIND L = PIVOT INDEX
         !
         LM = MIN(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
         !
         !         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
         !
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
            !
            !            INTERCHANGE IF NECESSARY
            !
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
            !
            !            COMPUTE MULTIPLIERS
            !
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
            !
            !            ROW ELIMINATION WITH COLUMN INDEXING
            !
            JU = MIN(MAX(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      ! DECK DGBSL
      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      ! ***BEGIN PROLOGUE  DGBSL
      ! ***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
      !             the factors computed by DGBCO or DGBFA.
      ! ***LIBRARY   SLATEC (LINPACK)
      ! ***CATEGORY  D2A2
      ! ***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
      ! ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
      ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
      ! ***DESCRIPTION
      !
      !      DGBSL solves the double precision band system
      !      A * X = B  or  TRANS(A) * X = B
      !      using the factors computed by DGBCO or DGBFA.
      !
      !      On Entry
      !
      !         ABD     DOUBLE PRECISION(LDA, N)
      !                 the output from DGBCO or DGBFA.
      !
      !         LDA     INTEGER
      !                 the leading dimension of the array  ABD .
      !
      !         N       INTEGER
      !                 the order of the original matrix.
      !
      !         ML      INTEGER
      !                 number of diagonals below the main diagonal.
      !
      !         MU      INTEGER
      !                 number of diagonals above the main diagonal.
      !
      !         IPVT    INTEGER(N)
      !                 the pivot vector from DGBCO or DGBFA.
      !
      !         B       DOUBLE PRECISION(N)
      !                 the right hand side vector.
      !
      !         JOB     INTEGER
      !                 = 0         to solve  A*X = B ,
      !                 = nonzero   to solve  TRANS(A)*X = B , where
      !                             TRANS(A)  is the transpose.
      !
      !      On Return
      !
      !         B       the solution vector  X .
      !
      !      Error Condition
      !
      !         A division by zero will occur if the input factor contains a
      !         zero on the diagonal.  Technically this indicates singularity
      !         but it is often caused by improper arguments or improper
      !         setting of LDA .  It will not occur if the subroutines are
      !         called correctly and if DGBCO has set RCOND .GT. 0.0
      !         or DGBFA has set INFO .EQ. 0 .
      !
      !      To compute  INVERSE(A) * C  where  C  is a matrix
      !      with  P  columns
      !            CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
      !            IF (RCOND is too small) GO TO ...
      !            DO 10 J = 1, P
      !               CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
      !         10 CONTINUE
      !
      ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
      !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
      ! ***ROUTINES CALLED  DAXPY, DDOT
      ! ***REVISION HISTORY  (YYMMDD)
      !    780814  DATE WRITTEN
      !    890531  Changed all specific intrinsics to generic.  (WRB)
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900326  Removed duplicate information from DESCRIPTION section.
      !            (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
      !
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
      ! ***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
         !
         !         JOB = 0 , SOLVE  A * X = B
         !         FIRST SOLVE L*Y = B
         !
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
         !
         !         NOW SOLVE  U*X = Y
         !
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
         !
         !         JOB = NONZERO, SOLVE  TRANS(A) * X = B
         !         FIRST SOLVE  TRANS(U)*Y = B
         !
         DO 60 K = 1, N
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
         !
         !         NOW SOLVE TRANS(L)*X = Y
         !
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      ! DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
      ! ***BEGIN PROLOGUE  DGEFA
      ! ***PURPOSE  Factor a matrix using Gaussian elimination.
      ! ***LIBRARY   SLATEC (LINPACK)
      ! ***CATEGORY  D2A1
      ! ***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
      ! ***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
      !              MATRIX FACTORIZATION
      ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
      ! ***DESCRIPTION
      !
      !      DGEFA factors a double precision matrix by Gaussian elimination.
      !
      !      DGEFA is usually called by DGECO, but it can be called
      !      directly with a saving in time if  RCOND  is not needed.
      !      (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
      !
      !      On Entry
      !
      !         A       DOUBLE PRECISION(LDA, N)
      !                 the matrix to be factored.
      !
      !         LDA     INTEGER
      !                 the leading dimension of the array  A .
      !
      !         N       INTEGER
      !                 the order of the matrix  A .
      !
      !      On Return
      !
      !         A       an upper triangular matrix and the multipliers
      !                 which were used to obtain it.
      !                 The factorization can be written  A = L*U  where
      !                 L  is a product of permutation and unit lower
      !                 triangular matrices and  U  is upper triangular.
      !
      !         IPVT    INTEGER(N)
      !                 an integer vector of pivot indices.
      !
      !         INFO    INTEGER
      !                 = 0  normal value.
      !                 = K  if  U(K,K) .EQ. 0.0 .  This is not an error
      !                      condition for this subroutine, but it does
      !                      indicate that DGESL or DGEDI will divide by zero
      !                      if called.  Use  RCOND  in DGECO for a reliable
      !                      indication of singularity.
      !
      ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
      !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
      ! ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
      ! ***REVISION HISTORY  (YYMMDD)
      !    780814  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900326  Removed duplicate information from DESCRIPTION section.
      !            (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
      !
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
      !
      !      GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      !
      ! ***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
         !
         !         FIND L = PIVOT INDEX
         !
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
         !
         !         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
         !
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
            !
            !            INTERCHANGE IF NECESSARY
            !
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
            !
            !            COMPUTE MULTIPLIERS
            !
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
            !
            !            ROW ELIMINATION WITH COLUMN INDEXING
            !
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      ! DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
      ! ***BEGIN PROLOGUE  DGESL
      ! ***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
      !             factors computed by DGECO or DGEFA.
      ! ***LIBRARY   SLATEC (LINPACK)
      ! ***CATEGORY  D2A1
      ! ***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
      ! ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
      ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
      ! ***DESCRIPTION
      !
      !      DGESL solves the double precision system
      !      A * X = B  or  TRANS(A) * X = B
      !      using the factors computed by DGECO or DGEFA.
      !
      !      On Entry
      !
      !         A       DOUBLE PRECISION(LDA, N)
      !                 the output from DGECO or DGEFA.
      !
      !         LDA     INTEGER
      !                 the leading dimension of the array  A .
      !
      !         N       INTEGER
      !                 the order of the matrix  A .
      !
      !         IPVT    INTEGER(N)
      !                 the pivot vector from DGECO or DGEFA.
      !
      !         B       DOUBLE PRECISION(N)
      !                 the right hand side vector.
      !
      !         JOB     INTEGER
      !                 = 0         to solve  A*X = B ,
      !                 = nonzero   to solve  TRANS(A)*X = B  where
      !                             TRANS(A)  is the transpose.
      !
      !      On Return
      !
      !         B       the solution vector  X .
      !
      !      Error Condition
      !
      !         A division by zero will occur if the input factor contains a
      !         zero on the diagonal.  Technically this indicates singularity
      !         but it is often caused by improper arguments or improper
      !         setting of LDA .  It will not occur if the subroutines are
      !         called correctly and if DGECO has set RCOND .GT. 0.0
      !         or DGEFA has set INFO .EQ. 0 .
      !
      !      To compute  INVERSE(A) * C  where  C  is a matrix
      !      with  P  columns
      !            CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
      !            IF (RCOND is too small) GO TO ...
      !            DO 10 J = 1, P
      !               CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
      !         10 CONTINUE
      !
      ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
      !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
      ! ***ROUTINES CALLED  DAXPY, DDOT
      ! ***REVISION HISTORY  (YYMMDD)
      !    780814  DATE WRITTEN
      !    890831  Modified array declarations.  (WRB)
      !    890831  REVISION DATE from Version 3.2
      !    891214  Prologue converted to Version 4.0 format.  (BAB)
      !    900326  Removed duplicate information from DESCRIPTION section.
      !            (WRB)
      !    920501  Reformatted the REFERENCES section.  (WRB)
      ! ***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
      !
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
      ! ***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
         !
         !         JOB = 0 , SOLVE  A * X = B
         !         FIRST SOLVE  L*Y = B
         !
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
         !
         !         NOW SOLVE  U*X = Y
         !
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
         !
         !         JOB = NONZERO, SOLVE  TRANS(A) * X = B
         !         FIRST SOLVE  TRANS(U)*Y = B
         !
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
         !
         !         NOW SOLVE TRANS(L)*X = Y
         !
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
