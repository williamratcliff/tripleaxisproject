!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: OPTIMIZATION_PROCEDURES
!!----   INFO: Module implementing several algorithms for global and local
!!----         optimization.
!!----
!!---- HISTORY
!!----    Update: February - 2005
!!----
!!----    April - 2004 Created by JRC
!!---- DEPENDENCIES
!!--++    use MATH_GEN, only: cp
!!----
!!---- VARIABLES
!!----    ERR_OPTIM
!!----    ERR_MESS_OPTIM
!!----    OPT_CONDITIONS_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       CG_QUASI_NEWTON
!!--++       CHOLA                 [Private]
!!----       CSENDES_GLOBAL
!!--++       FUN                   [Private]
!!----       INIT_ERR_OPTIM
!!----       INIT_OPT_CONDITIONS
!!--++       LOCAL                 [Private]
!!----       NELDER_MEAD_SIMPLEX
!!--++       PRINT_TRI_MATRIX      [Private]
!!--++       SYMINV                [Private]
!!--++       UPDATE                [Private]
!!--++       URDMN                 [Private]
!!----
!!
 Module Optimization_Procedures
    !---- Use Files ----!
    use Math_Gen, only: cp

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public :: Csendes_Global,Cg_Quasi_Newton,Init_Err_Optim, Init_Opt_Conditions, Nelder_Mead_Simplex

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private :: Chola, Local, Print_Tri_Matrix, Syminv, Update, Urdmn, Fun

    !---- Variable Definitions ----!

    !!----
    !!---- ERR_MESS_OPTIM
    !!----    character(len=150), public :: Err_Mess_Optim
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: Err_Mess_Optim

    !!----
    !!---- ERR_OPTIM
    !!----    logical, public :: Err_Optim
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: Err_Optim

    !!----
    !!---- TYPE :: OPT_CONDITIONS_TYPE
    !!--<<
    !!---- Type, public :: Opt_Conditions_Type
    !!----    integer      :: nmeth ! (in)  nmeth=0 => conjugate gradient
    !!----                          !       nmeth=1 => the BFGS method is used
    !!----    integer      :: mxfun ! (in)  Maximum number function calls
    !!----    integer      :: iout  ! (in)  Printing parameter,
    !!----                          !       if iout= 0 no printing for Quasi_Newton & Conjugate Gradient
    !!----                          !       if iout= 0 partial printing for Simplex (<0 no printing)
    !!----                          !       if iout>0 printing each iout iterations/evaluations.
    !!----    integer      :: loops ! (in)  Useful for SIMPLEX method: = 0
    !!----                          !
    !!----    integer      :: iquad ! (in)  For SIMPLEX, if iquad/= 0 fitting to a quadratic
    !!----                          !
    !!----    integer      :: nflag ! (out) Flag value states which condition caused the exit of
    !!----                          !       the optimization subroutine
    !!----                          !       If NFLAG=0, the algorithm has converged.
    !!----                          !       If NFLAG=1, the maximum number of function
    !!----                          !          evaluations have been used.
    !!----                          !       If NFLAG=2, the linear search has failed to
    !!----                          !          improve the function value. This is the
    !!----                          !          usual exit if either the function or the
    !!----                          !          gradient is incorrectly coded.
    !!----                          !       If NFLAG=3, The search vector was not
    !!----                          !          a descent direction. This can only be caused
    !!----                          !          by roundoff,and may suggest that the
    !!----                          !          convergence criterion is too strict.
    !!----    integer      :: ifun  ! (out) Total number of function and gradient evaluations
    !!----    integer      :: iter  ! (out) Total number of search directions used in the algorithm
    !!----    real(kind=cp):: eps   ! (in)  Convergence occurs when the norm of the gradient
    !!----                          !       is less than or equal to EPS times the maximum
    !!----                          !       of one and the norm of the vector X
    !!----                          !       Initialized to 1.0E-6
    !!----    real(kind=cp):: acc   ! (in)  ACC is a user supplied estimate of machine
    !!----                          !       accuracy. A linear search is unsuccessfully
    !!----                          !       terminated when the norm of the step size
    !!----                          !       becomes smaller than ACC. In practice,
    !!----                          !       ACC=10.0E-20 has proved satisfactory.
    !!----                          !       This is the default value.
    !!----                          !       For Simplex method the meaning is different
    !!----                          !       (see below) and this should be changed to 1.0e-6
    !!---- End Type Opt_Conditions_Type
    !!----
    !!----    This TYPE has been introduced to simplify the call to optimization
    !!----    procedures. It contains the optimization parameters useful for different
    !!----    algorithms.
    !!----    All integer components are initialized to zero and the real components
    !!----    are initilized as indicated below.
    !!----    A variable of this type should be defined by the user and all their
    !!----    input parameters (in) must be provided before calling the procedures.
    !!----    On output from the procedure the (out) items are provided for checking.
    !!-->>
    !!---- Update: February - 2005
    !!
    Type, public :: Opt_Conditions_Type
      integer        :: nmeth
      integer        :: mxfun
      integer        :: loops
      integer        :: iquad
      integer        :: iout
      integer        :: nflag
      integer        :: ifun
      integer        :: iter
      real(kind=cp)  :: eps
      real(kind=cp)  :: acc
    End type Opt_Conditions_Type

 Contains

    !!----
    !!---- Subroutine Cg_Quasi_Newton(Model_Functn,Nparm,X,F,G,C,Ipr)
    !!----    integer,                          intent(in)     :: n      ! The number of variables in the function to be
    !!----                                                               ! minimized.
    !!----    real(kind=cp),dimension(n),       intent(in out) :: x      ! The vector containing the current estimate to
    !!----                                                               ! the minimizer.
    !!----                                                               ! In -> Must contain an initial estimate supplied by the user.
    !!----                                                               ! Out-> X will hold the best estimate to the minimizer obtained
    !!----    real(kind=cp),                    intent(   out) :: f      ! Out-> F will contain the lowest value of the object function obtained.
    !!----    real(kind=cp),dimension(n),       intent(   out) :: g      ! Out-> G =(g(1),...g(n)) will contain the elements of the gradient of
    !!----                                                               !       F evaluated at the point contained in X=(x(1),...x(N))
    !!----    type(Opt_conditions),             intent(in out) :: C      ! Conditions for the algorithm. C is of type(Opt_conditions)
    !!----    integer, optional,                intent(in)     :: ipr   ! Logical unit for printing if the parameter C%IOUT /= 0.
    !!----
    !!--<<    Interface
    !!----       Subroutine Model_Functn(n,x,f,g)
    !!----          use Math_Gen,  only: cp
    !!----          integer,                    intent(in)     :: n
    !!----          real(kind=cp),dimension(:), intent(in)     :: x
    !!----          real(kind=cp),              intent(out)    :: f
    !!----          real(kind=cp),dimension(:), intent(out)    :: g
    !!----       End Subroutine Model_Functn
    !!-->>    End Interface
    !!----
    !!----    Minimization Of Unconstrained Multivariate Functions
    !!----    Subroutine CG_QUASI_NEWTON minimizes an unconstrained nonlinear
    !!----    scalar valued function of a vector variable X either by the
    !!----    BFGS variable metric algorithm or by a beale restarted conjugate
    !!----    gradient algorithm.
    !!----    (Authors: D.F. Shanno and K.H. Phua, The original name of the
    !!----    subroutine was CONMIN)
    !!----    ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE 6 (DECEMBER 1980), 618-622
    !!----
    !!--<<
    !!---- REMARKS:    In addition to the specified values in the above
    !!----             argument list, the user must supply a subroutine
    !!----             "Model_Functn" which calculates the function and gradient at
    !!----             X and places them in F and G(1),...,G(N) respectively.
    !!-->>
    !!--..             An example subroutine for the Rosenbrock function is:
    !!--..
    !!--..             Subroutine Model_Functn(n,x,f,g)
    !!--..                integer,parameter  :: cp=selected_real_kind(14, 300)
    !!--..                integer,                    intent(in) :: n
    !!--..                real(kind=cp),dimension(:), intent(in) :: x
    !!--..                real(kind=cp),              intent(out):: f
    !!--..                real(kind=cp),dimension(:), intent(out):: g
    !!--..                real(kind=cp) :: t1,t2
    !!--..
    !!--..                t1=x(2)-x(1)*x(1)
    !!--..                t2=1.0-x(1)
    !!--..                f=100.0*t1*t1+t2*t2
    !!--..                g(1)=-400.0*t1*x(1)-2.0*t2
    !!--..                g(2)=200.0*t1
    !!--..
    !!--..                return
    !!--..             End Subroutine Model_Functn
    !!--..
    !!--..    Code converted using TO_F90 by Alan Miller
    !!--..    Modified and adapted to F-language by Juan Rodriguez-Carvajal
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Cg_Quasi_Newton(Model_Functn, N, X, F, G, C, Ipr)
       !---- Arguments ----!
       integer,                    intent(in)       :: n
       real(kind=cp),dimension(n), intent(in out)   :: x
       real(kind=cp),              intent(   out)   :: f
       real(kind=cp),dimension(n), intent(   out)   :: g
       type(Opt_conditions_Type),  intent(in out)   :: C
       integer, optional,          intent(in)       :: Ipr

       Interface
          Subroutine Model_Functn(n,x,f,g)
             use Math_Gen,  only: cp
             integer,                    intent(in) :: n
             real(kind=cp),dimension(:), intent(in) :: x
             real(kind=cp),              intent(out):: f
             real(kind=cp),dimension(:), intent(out):: g
          End Subroutine Model_Functn
       End Interface

       !---- Local variables ----!
       logical       :: rsw
       integer       :: ioutk, nx
       integer       :: ng, nry, nrd, ncons, ncons1, ncons2, nrst, i, ncalls
       integer       :: nxpi, ngpi, nrdpi, nrypi, ij, j, ii, ngpj
       real(kind=cp) :: fp,fmin,alpha,at,ap,gsq,dg,dg1
       real(kind=cp) :: dp,step,dal,u1,u2,u3,u4
       real(kind=cp) :: xsq,rtst

       !---- W is a vector of working storagE.If C%NMETH=0,
       !---- W must be dimensioned 5*N+2. If C%NMETH=1,
       !---- W must be dimensioned N*(N+7)/2. In both cases,
       !---- W must be real.
       real(kind=cp),dimension(:),allocatable :: w

       !---- Initialize ITER,IFUN,NFLAG,and IOUTK,which counts output iterations.
       c%iter=0
       c%ifun=0
       ioutk=0
       c%nflag=0

       !---- Set parameters to extract vectors from W.
       !---- W(I) holds the search vector,W(NX+I) holds the best current
       !---- estimate to the minimizer,and W(NG+I) holds the gradient
       !---- at the best current estimate.
       nx=n
       ng=nx+n

       !---- Test which method is being used.
       !---- If C%NMETH=0, W(NRY+I) holds the restart Y vector and
       !---- W(NRD+I) holds the restart SEARCH vector.
       !---- If C%NMETH=1,W(NCONS+I) holds the approximate inverse Hessian.
       if (c%nmeth == 1) then
          ncons=3*n
          if (allocated(w)) deallocate(w)
          allocate(w(n*(n+7)/2))     !W has been removed from the argument
       else
          nry=ng+n
          nrd=nry+n
          ncons=5*n
          ncons1=ncons+1
          ncons2=ncons+2
          if (allocated(w)) deallocate(w)
          allocate(w(5*n+2))
       end if

       !---- Calculate the function and gradient at the initial
       !---- point and initialize NRST,which is used to determine
       !---- whether a beale restart is being done. NRST=N means that this
       !---- iteration is a restart iteration. initialize RSW,which indicates
       !---- that the current search direction is a gradient direction.
       do_20: do
          call Model_Functn(n,x,f,g)
          c%ifun=c%ifun+1
          nrst=n
          rsw=.true.

          !---- Calculate the initial search direction , the norm of X squared,
          !---- and the norm of G squared. DG1 is the current directional
          !---- derivative,while XSQ and GSQ are the squared norms.
          dg1=0.0
          xsq=0.0
          do i=1,n
             w(i)=-g(i)
             xsq=xsq+x(i)*x(i)
             dg1=dg1-g(i)*g(i)
          end do
          gsq=-dg1

          !---- Test if the initial point is the minimizer.
          if (gsq <= c%eps*c%eps*max(1.0_cp,xsq)) return

          !---- Begin the major iteration loop. NCALLS is used to guarantee that
          !---- at least two points have been tried when C%NMETH=0. FMIN is the
          !---- current function value.
          do_40: do
             fmin=f
             ncalls=c%ifun

             !---- If output is desired,test if this is the correct iteration
             !---- and if so, write output.
             if (present(ipr)) then
                if (c%iout  /= 0) then
                   if (ioutk == 0)then
                      write(unit=ipr,fmt="(a,i5,a,i6,2(a,f20.8))")"  Iteration #",c%iter, &
                            "    Num F-eval:",c%ifun,"    F-value = ",fmin,"    G-squared = ",gsq
                   end if
                   ioutk=ioutk+1
                   if (ioutk == c%iout) ioutk=0
                end if
             end if

             !---- Begin linear search. ALPHA is the steplength.
             !---- Set ALPHA to the nonrestart conjugate gradient ALPHA.
             alpha=alpha*dg/dg1

             !---- If C%NMETH=1 or a restart has been performed, set ALPHA=1.0.
             if (nrst == 1 .OR. c%nmeth == 1) alpha=1.0

             !---- If a gradient direction is used, set ALPHA=1.0/DSQRT(GSQ),
             !---- which scales the initial search vector to unity.
             if (rsw) alpha=1.0/sqrt(gsq)

             !---- The linear search fits a cubic to F and DAL, the function
             !---- and its derivative at ALPHA, and to FP and DP,the function
             !---- and derivative at the previous trial point AP.
             !---- Initialize AP ,FP,and DP.
             ap=0.0
             fp=fmin
             dp=dg1

             !---- Save the current derivative to scale the next search vector.
             dg=dg1

             !---- Update the iteration.
             c%iter=c%iter+1

             !---- Calculate the current steplength  and store the current X and G.
             step=0.0
             do i=1,n
                step=step+w(i)*w(i)
                nxpi=nx+i
                ngpi=ng+i
                w(nxpi)=x(i)
                w(ngpi)=g(i)
             end do
             step=SQRT(step)

             !---- Begin the linear search iteration.
             !---- Test for failure of the linear search.
             do_80: do
                if (alpha*step <= c%acc) then
                   !---- Test if direction is a gradient direction.
                   if (.not. rsw) cycle do_20
                   c%nflag=2
                   return
                end if

                !---- Calculate the trial point.
                do i=1,n
                   nxpi=nx+i
                   x(i)=w(nxpi)+alpha*w(i)
                end do

                !---- Evaluate the function at the trial point.
                call Model_Functn(n,x,f,g)

                !---- Test if the maximum number of function calls have been used.
                c%ifun=c%ifun+1
                if (c%ifun > c%mxfun) then
                   c%nflag=1
                   return
                end if

                !---- Compute the derivative of f at alpha.
                dal=0.0
                do i=1,n
                   dal=dal+g(i)*w(i)
                end do

                !---- Test whether the new point has a negative slope but a higher
                !---- function value than ALPHA=0. If this is the case,the search
                !---- has passed through a local maximun and is heading for a distant local
                !---- minimum.
                if (f > fmin .and. dal < 0.0) then
                   !---- A relative maximum has been passed. Reduce ALPHA and restart the search.
                   alpha=alpha/3.0
                   ap=0.0
                   fp=fmin
                   dp=dg
                   cycle do_80
                end if

                !---- If not, test whether the steplength criteria have been met.
                if (.not. (f > fmin+0.0001*alpha*dg .or. abs(dal/dg)  > 0.9)) then
                   !---- If they have been met, test if two points have been tried
                   !---- If NMETH=0 and if the true line minimum has not been found.
                   if(.not.((c%ifun-ncalls) <= 1 .and. abs(dal/dg) > c%eps .and. c%nmeth == 0)) exit do_80
                end if

                !---- A new point must be tried. Use cubic interpolation to find
                !---- the trial point AT.
                u1=dp+dal-3.0*(fp-f)/(ap-alpha)
                u2=u1*u1-dp*dal
                if (u2 < 0.0) u2=0.0
                u2=sqrt(u2)
                at=alpha-(alpha-ap)*(dal+u2-u1)/(dal-dp+2.0*u2)

                !---- Test whether the line minimum has been bracketed.
                if ( dal/dp > 0.0) then
                   !---- The minimum has not been bracketed. Test if both points are
                   !---- greater than the minimum and the trial point is sufficiently
                   !---- smaller than either.
                   if (.not.(dal > 0.0 .and. 0.0 < at .and. at < (0.99*min(ap,alpha)))) then
                      !---- Test if both points are less than the minimum and the trial point
                      !---- is sufficiently large.
                      if (.not.(dal <= 0.0 .and. at > (1.01*max(ap,alpha)))) then
                         !---- If the trial point is too small,double the largest prior point.
                         if (dal <= 0.0) at=2.0*max(ap,alpha)
                         !---- If the trial point is too large, halve the smallest prior point.
                         if (dal > 0.0) at=min(ap,alpha)/2.0
                      end if
                   end if
                else
                   !---- The minimum has been bracketed. Test whether the trial point lies
                   !---- sufficiently within the bracketed interval.
                   !---- If it does not, choose at as the midpoint of the interval.
                   if (at < (1.01*min(alpha,ap)) .or. at > (0.99*max(alpha,ap))) at=(alpha+ap)/2.0
                end if

                !---- Set AP=ALPHA, ALPHA=AT,and continue search.
                ap=alpha
                fp=f
                dp=dal
                alpha=at

             end do do_80

             !---- The line search has converged. Test for convergence of the algorithm.
             gsq=0.0
             xsq=0.0
             do i=1,n
                gsq=gsq+g(i)*g(i)
                xsq=xsq+x(i)*x(i)
             end do
             if (gsq <= c%eps*c%eps*max(1.0_cp,xsq)) return

             !---- Search continues. Set W(I)=ALPHA*W(I),the full step vector.
             w(1:n)=alpha*w(1:n)

             !---- Compute the new search vector. First test whether a
             !---- conjugate gradient or a variable metric vector is used.
             if (c%nmeth /= 1) then
                !---- Conjugate gradient update section.
                !---- Test if a Powell restart is indicated.
                rtst=0.0
                do i=1,n
                   ngpi=ng+i
                   rtst=rtst+g(i)*w(ngpi)
                end do
                if ( abs(rtst/gsq) > 0.2) nrst=n

                !---- If a restart is indicated, save the current D and Y
                !---- as the beale restart vectors and save D'Y and Y'Y
                !----- in W(NCONS+1) and W(NCONS+2).
                if (nrst == n) then
                   w(ncons+1)=0.0
                   w(ncons+2)=0.0
                   do i=1,n
                      nrdpi=nrd+i
                      nrypi=nry+i
                      ngpi=ng+i
                      w(nrypi)=g(i)-w(ngpi)
                      w(nrdpi)=w(i)
                      w(ncons1)=w(ncons1)+w(nrypi)*w(nrypi)
                      w(ncons2)=w(ncons2)+w(i)*w(nrypi)
                   end do
                end if

                !---- Calculate  the restart Hessian times the current gradient.
                u1=0.0
                u2=0.0
                do  i=1,n
                  nrdpi=nrd+i
                  nrypi=nry+i
                  u1=u1-w(nrdpi)*g(i)/w(ncons1)
                  u2=u2+w(nrdpi)*g(i)*2.0/w(ncons2)-w(nrypi)*g(i)/w(ncons1)
                end do
                u3=w(ncons2)/w(ncons1)
                do  i=1,n
                  nxpi=nx+i
                  nrdpi=nrd+i
                  nrypi=nry+i
                  w(nxpi)=-u3*g(i)-u1*w(nrypi)-u2*w(nrdpi)
                end do

                !---- If this is a restart iteration, W(NX+I) contains the new search
                !---- vector.
                if (nrst /= n) then
                   !---- Not a restart iteration. Calculate the restart Hessian
                   !---- times the current Y.
                   u1=0.0
                   u2=0.0
                   u3=0.0
                   u4=0.0
                   do i=1,n
                      ngpi=ng+i
                      nrdpi=nrd+i
                      nrypi=nry+i
                      u1=u1-(g(i)-w(ngpi))*w(nrdpi)/w(ncons1)
                      u2=u2-(g(i)-w(ngpi))*w(nrypi)/w(ncons1)+2.0*w(nrdpi)*(g(i)-w(ngpi))/w(ncons2)
                      u3=u3+w(i)*(g(i)-w(ngpi))
                   end do
                   step=0.0
                   do i=1,n
                      ngpi=ng+i
                      nrdpi=nrd+i
                      nrypi=nry+i
                      step=(w(ncons2)/w(ncons1))*(g(i)-w(ngpi)) +u1*w(nrypi)+u2*w(nrdpi)
                      u4=u4+step*(g(i)-w(ngpi))
                      w(ngpi)=step
                   end do

                   !---- Calculate the doubly updated Hessian times the current
                   !---- gradient to obtain the search vector.
                   u1=0.0
                   u2=0.0
                   do i=1,n
                      u1=u1-w(i)*g(i)/u3
                      ngpi=ng+i
                      u2=u2+(1.0+u4/u3)*w(i)*g(i)/u3-w(ngpi)*g(i)/u3
                   end do
                   do i=1,n
                      ngpi=ng+i
                      nxpi=nx+i
                      w(nxpi)=w(nxpi)-u1*w(ngpi)-u2*w(i)
                   end do
                end if

                !---- Calculate the derivative along the new search vector.
                dg1=0.0
                do i=1,n
                   nxpi=nx+i
                   w(i)=w(nxpi)
                   dg1=dg1+w(i)*g(i)
                end do

                !---- If the new direction is not a descent direction,stop.
                if (dg1 > 0.0) then
                   ! roundoff has produced a bad direction.
                   c%nflag=3
                   return
                end if

                !---- Update nrst to assure at least one restart every n iterations.
                if (nrst == n)nrst=0
                nrst=nrst+1
                rsw=.false.
                cycle do_40
             end if

             !---- A variable metric algoritm is being used. Calculate Y and D'Y.
             u1=0.0
             do i=1,n
                ngpi=ng+i
                w(ngpi)=g(i)-w(ngpi)
                u1=u1+w(i)*w(ngpi)
             end do

             !---- If RSW=.TRUE.,set up the initial scaled approximate Hessian.
             if (rsw) then
                !---- calculate y'y.
                u2=0.0
                do i=1,n
                   ngpi=ng+i
                   u2=u2+w(ngpi)*w(ngpi)
                end do

                !---- Calculate the initial Hessian as H=(P'Y/Y'Y)*I
                !---- and the initial U2=Y'HY and W(NX+I)=HY.
                ij=1
                u3=u1/u2
                do i=1,n
                   do j=i,n
                      ncons1=ncons+ij
                      w(ncons1)=0.0
                      if (i == j)w(ncons1)=u3
                      ij=ij+1
                   end do
                   nxpi=nx+i
                   ngpi=ng+i
                   w(nxpi)=u3*w(ngpi)
                end do
                u2=u3*u2

             else  !RSW
                !---- Calculate W(NX+I)=HY and U2=Y'HY.
                u2=0.0
                do i=1,n
                   u3=0.0
                   ij=i
                   if (i /= 1) then
                      ii=i-1
                      do j=1,ii
                         ngpj=ng+j
                         ncons1=ncons+ij
                         u3=u3+w(ncons1)*w(ngpj)
                         ij=ij+n-j
                      end do
                   end if
                   do j=i,n
                      ncons1=ncons+ij
                      ngpj=ng+j
                      u3=u3+w(ncons1)*w(ngpj)
                      ij=ij+1
                   end do
                   ngpi=ng+i
                   u2=u2+u3*w(ngpi)
                   nxpi=nx+i
                   w(nxpi)=u3
                end do
             end if !RSW

             !---- Calculate the updated approximate Hessian.
             u4=1.0+u2/u1
             do i=1,n
                nxpi=nx+i
                ngpi=ng+i
                w(ngpi)=u4*w(i)-w(nxpi)
             end do
             ij=1
             do i=1,n
                nxpi=nx+i
                u3=w(i)/u1
                u4=w(nxpi)/u1
                do j=i,n
                   ncons1=ncons+ij
                   ngpj=ng+j
                   w(ncons1)=w(ncons1)+u3*w(ngpj)-u4*w(j)
                   ij=ij+1
                end do
             end do

             !---- Calculate the new search direction W(I)=-HG and its derivative.
             dg1=0.0
             do i=1,n
                u3=0.0
                ij=i
                if (i /= 1) then
                   ii=i-1
                   do j=1,ii
                      ncons1=ncons+ij
                      u3=u3-w(ncons1)*g(j)
                      ij=ij+n-j
                   end do
                end if
                do j=i,n
                   ncons1=ncons+ij
                   u3=u3-w(ncons1)*g(j)
                   ij=ij+1
                end do
                dg1=dg1+u3*g(i)
                w(i)=u3
             end do

             !---- Test for a downhill direction.
             if (dg1 > 0.0) then
                !---- Roundoff has produced a bad direction.
                c%nflag=3
                return
             end if
             rsw=.false.
          end do do_40
          exit
       end do do_20

       return
    End Subroutine Cg_Quasi_Newton

    !!--++
    !!--++ Subroutine Chola(A, N, U, Nullty, Ifault, Rmax, R)
    !!--++    real (kind=cp),dimension(:), intent(in)   :: A        ! a +ve definite matrix stored in lower-triangular form
    !!--++    integer,                     intent(in)   :: N        ! The order of A
    !!--++    real (kind=cp),dimension(:), intent(out)  :: U        ! a lower triangular matrix such that U*U' = A.
    !!--++    integer,                     intent(out)  :: Nullty   ! the rank deficiency of A.
    !!--++    integer,                     intent(out)  :: Ifault   ! Error Code: 1 if N < 1, 2 If A Is not +ve semi-definite
    !!--++                                                          !             0 Otherwise
    !!--++    real (kind=cp),              intent(out)  :: Rmax     ! an estimate of the relative accuracy of the diagonal elements of U.
    !!--++    real (kind=cp),dimension(:), intent(out)  :: R        ! array containing bounds on the relative accuracy of each diagonal element of U.
    !!--++
    !!--++    ALGORITHM AS6, Applied statistics, VOL.17, 1968, with  modifications
    !!--++    by A.J.MILLER
    !!--<<
    !!--++    Note: Eta should be set equal to the smallest +ve value such that
    !!--++          1.0_cp + eta is calculated as being greater than 1.0_cp in
    !!--++          the accuracy
    !!-->>
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Chola(A, N, U, Nullty, Ifault, Rmax, R)
       !---- Arguments ----!
       real(kind=cp),dimension(:), intent(in)   :: a
       integer,                    intent(in)   :: n
       real(kind=cp),dimension(:), intent(out)  :: u
       integer,                    intent(out)  :: nullty
       integer,                    intent(out)  :: ifault
       real(kind=cp),              intent(out)  :: rmax
       real(kind=cp),dimension(:), intent(out)  :: r

       !---- Local variables ----!
       integer                   :: i, icol, irow, j, k, l, m
       real (kind=cp), parameter :: eta = epsilon(1.0_cp), zero = 0.0_cp
       real (kind=cp)            :: rsq, w

       ifault = 1
       if (n > 0) then
          ifault = 2
          nullty = 0
          rmax = eta
          r(1) = eta
          j = 1
          k = 0
          do icol = 1, n         !     Factorize column by column, ICOL = Column no.
             l = 0
             do irow = 1, icol    !     IROW = Row number within column ICOL
                k = k + 1
                w = a(k)
                if (irow == icol) rsq = (w*eta) ** 2
                m = j
                do i = 1, irow
                   l = l + 1
                   if (i == irow) exit
                   w = w - u(l) * u(m)
                   if (irow == icol) rsq = rsq + (u(l)**2*r(i)) ** 2
                   m = m + 1
                end do
                if (irow == icol) exit
                if (u(l) /= zero) then
                   u(k) = w / u(l)
                else
                   u(k) = zero
                   if (abs(w) > abs(rmax*a(k))) return
                end if
             end do

             !---- End of row, estimate relative accuracy of diagonal element.
             rsq = sqrt(rsq)
             if (abs(w) > 5.0*rsq) then
                if (w < zero) return
                u(k) = sqrt(w)
                r(i) = rsq / w
                if (r(i) > rmax) rmax = r(i)
             else
                u(k) = zero
                nullty = nullty + 1
             end if
             j = j + icol
          end do
          ifault = zero
       end if

       return
    End Subroutine Chola

    !!----
    !!---- Subroutine Csendes_Global(Model_Functn,Mini, Maxi, Nparm, N100, Ng0, Nsig, X0, Nc, F0, Ipr)
    !!----    Model_Functn : is the dummy name of the objective function to be optimized.
    !!----    real(kind=cp), dimension(:),   intent(in)    :: mini   ! Vector of Length Nparm Containing The Lower Bounds
    !!----    real(kind=cp), dimension(:),   intent(in)    :: maxi   ! Vector of Length Nparm Containing The Upper Bounds
    !!----    integer,                       intent(in)    :: nparm  ! Number of Parameters
    !!----    integer,                       intent(in out):: n100   ! Number of Sample Points To Be Drawn Uniformly In One Cycle < 10000
    !!----                                                           ! Suggested value is 100*Nparm
    !!----    integer,                       intent(in out):: ng0    ! Number of Best Points Selected From The Actual Sample
    !!----                                                           ! The Suggested Value Is Twice The Expected Number Of Local Minima
    !!----    integer,                       intent(in)    :: nsig   ! Convergence Criterion
    !!----    real(kind=cp), dimension(:,:), intent(in out):: x0     ! Output (up to Npar x Nc) Matrix Containing Nc Local Minimizers Found.
    !!----    integer,                       intent(out)   :: nc     ! Number of Different Local Minimizers Found
    !!----    real(kind=cp), dimension(:),   intent(in out):: f0     ! Output Vector Of Nc  Objective Function Values,
    !!----                                                           ! F0(I) Belongs To The Parameters X0(1,I), X0(2,I), ..., X0(Nparm,I)
    !!----    integer,                       intent(in)    :: ipr    ! Printing Information
    !!----
    !!--<<    Interface
    !!----       Subroutine Model_Functn(Nparm,X, F)
    !!----          real(kind=cp),dimension(:), intent(in)  :: x
    !!----          real(kind=cp),              intent(out) :: f
    !!----          integer,                    intent(in)  :: nparm
    !!----       End Subroutine Model_Functn
    !!----    End Interface
    !!-->>
    !!----
    !!----    Global optimization using the Boender-Timmer-Rinnoy Kan algorithm
    !!----
    !!--..    GLOBAL MINIMUM OF FUNCTION OF N VARIABLES USING A LOCAL SEARCH METHOD
    !!--..
    !!--..    Information
    !!--..
    !!--..    Convergence criterium: The accuracy required in the parameter estimates.
    !!--..    this convergence criterion is satisfied if on two successive iterations
    !!--..    the parameter estimates agree, component by component, to nsig digits.
    !!--..    the suggested value is 6.
    !!--..
    !!--..        SHORT DESCRIPTION OF THE GLOBAL OPTIMIZATION SUBROUTINE
    !!--..        -------------------------------------------------------
    !!--..        (Based in the original explanations by Tibor Csendes)
    !!--..
    !!--..   Global optimization is a part of nonlinear optimization, it deals with
    !!--..   problems with (possibly) several local minima. The presented method is
    !!--..   stochastic (i.e. not deterministic). The framework procedure, the GLOBAL
    !!--..   routine gives a computational evidence, that the best local minimum found is
    !!--..   with high probability the global minimum.  This routine calls a local search
    !!--..   routine, and a routine for generating random numbers.
    !!--..
    !!--..   Let F(X) be a real function of NPARM parameters and we are looking for
    !!--..   parameter values X(I) from the given intervals [MIN(I), MAX(I)] for each
    !!--..   I = 1, 2, ..., NPARM.  The problem is to determine such a point X*, that the
    !!--..   function value F(X) is greater than or equal to F(X*) for every X in the
    !!--..   NPARM-dimensional interval specified by MIN(I)'s and MAX(I)'s.
    !!--..
    !!--..   The given version allows 15 parameters to be optimized.  For modifying the
    !!--..   program to accept larger problems change the numbers 15 everywhere in the
    !!--..   source lines to the new value N, and change 120 in a declaration to N(N+1)/2.
    !!--..
    !!--..   I. THE ALGORITHM
    !!--..   ----------------
    !!--..
    !!--..   The algorithm consists of the following steps:
    !!--..
    !!--..   0. step: initialization of the parameters of the algorithm.  X0 is the set of
    !!--..            the local minima found, F0 contains corresponding values of the
    !!--..            objective function.  The local minima are ordered increasingly
    !!--..            according to the function values.  X1 is the set of points starting
    !!--..            from which the local search procedure led to a local minimum.  These
    !!--..            points are called seed points, and as such they are used as seed
    !!--..            points in the clustering phase.  At the start the mentioned sets are
    !!--..            empty.  The subroutine checks the parameter bounds, and gives error
    !!--..            message, and stops if they are not correct, or if they do not meet
    !!--..            the requirements.
    !!--..
    !!--..   1. step: generate sample points with uniform distribution, and add them to the
    !!--..            sample.  The number of generated points is specified by NSAMPL.
    !!--..
    !!--..   2. step: generate the actual sample, as the best 100 * NSEL/NSAMPL percentage
    !!--..            of the sample points generated so far (according to the objective
    !!--..            function values).  This actual sample contains in general NSEL more
    !!--..            points than the previous one.
    !!--..
    !!--..   3. step: form clusters in the actual sample by Single Linkage method, by
    !!--..            growing the clusters first around the seed points (elements of the
    !!--..            sets X0 and X1).  A new point will join a cluster, if there is a point
    !!--..            in the cluster with which the distance is less than a critical
    !!--..            distance calculated automatically by the program for the given
    !!--..            problem, and if the point in the cluster has a smaller objective
    !!--..            function value than that of the considered one.  The critical distance
    !!--..            depends on the number of points in the whole sample, and on the
    !!--..            dimension of the problem NPARM.  If all points of the actual sample
    !!--..            are successfully ordered to some of the existing clusters, then comes
    !!--..            the 5th step.
    !!--..
    !!--..   4. step: start local search from the actual sample points not yet clustered in
    !!--..            ascending order by the values of the respective function values. If
    !!--..            the result of the local search is close to an element of the sets X0
    !!--..            and X1, then the starting point will be added to the set X1.  If the
    !!--..            result of the local search cannot be clustered to any of the existing
    !!--..            clusters then the point is regarded as a new local minimizer, and is
    !!--..            added to X0.  Choose now this result point of the local search to be a
    !!--..            seed point, and try to find (not clustered) points in the current
    !!--..            sample that can be clustered to this one.  If a new local minimum was
    !!--..            found in step 4, go to the step 1.
    !!--..
    !!--..   5. step: determine the element of the set X0 with the smallest function value.
    !!--..            This is the candidate of the program for the global minimizer.  Stop.
    !!--..
    !!--..
    !!--..   The presented program is a modification of the algorithm by Boender et al.
    !!--..   (see [1]).  The following changes were made.
    !!--..
    !!--..   1. The local search procedure is an algorithm of Quasi-Newton type which uses
    !!--..      the so called DFP (Davidon-Fletcher-Powell) update formula.  The comparison
    !!--..      of this local search method with others can be found in [2].  For smooth
    !!--..      objective functions this seems to be a good choice, for problems with
    !!--..      discontinuous objective function or derivatives the robust random search
    !!--..      method UNIRANDI (more details in [5]) can be recommended.
    !!--..
    !!--..   2. The subroutine GLOBAL will automatically scale the parameters to be
    !!--..
    !!--..      optimized (as in [3]) in order to keep all the starting points of the
    !!--..      parameters between -1 and 1 for the local minimizer procedure.  The user
    !!--..      does not have to care about the scaling, because the results are
    !!--..      transformed back to the original interval before each output, and before
    !!--..      giving the control back to the calling program.
    !!--..
    !!--..   3. We have also made use of later results [6] on this method.  For example, the
    !!--..      condition used at clustering has been changed to another one.
    !!--..
    !!--..   4. Instead of the Euclidean distance the greatest difference in absolute
    !!--..      values is used in our GLOBAL routine.  This change leads to a simpler and
    !!--..      quicker version of the algorithm.
    !!--..
    !!--..   II. HOW TO CALL THE SUBROUTINE GLOBAL
    !!--..   -------------------------------------
    !!--..
    !!--..       CALL GLOBAL (Model_Functn, MIN, MAX, NPARM, NSAMPL, NSEL, IPR, NSIG, X0, NC, F0)
    !!--..
    !!--..   MIN and MAX are real vectors (input) of NPARM elements.  They specify an
    !!--..   NPARM-dimensional interval that contains the starting points for the local
    !!--..   searches, i.e. the Ith coordinate of a random sample point, X(I) is always in
    !!--..   [MIN(I), MAX(I)].  The values of this vector will not be changed during the run
    !!--..   of the GLOBAL subroutine.
    !!--..
    !!--..   NPARM is an integer constant (input) containing the number of the parameters
    !!--..   to be optimized.  The value will not be changed during the run of GLOBAL.
    !!--..
    !!--..   NSAMPL is an integer constant (input) containing the number of the sampling
    !!--..   points to be generated in the parameter space with uniform distribution in one
    !!--..   loop.  The value will not be changed during the run of the GLOBAL subroutine.
    !!--..
    !!--..   NSEL is an integer constant (input).  In one sampling the number of points
    !!--..   selected for starting points of local search routine is NSEL.  These points are
    !!--..   those with the smallest objective function value.  The value will not be
    !!--..   changed by the GLOBAL routine.
    !!--..
    !!--..   IPR is an integer constant (input), the FORTRAN logical file number for the
    !!--..   output file. The value will not be changed during the GLOBAL run.  The user
    !!--..   should take care on the corresponding OPEN statement in the main program.  In
    !!--..   the IBM-PC version of the GLOBAL routine each output of the routine will be
    !!--..   repeated for the default unit (the screen).
    !!--..
    !!--..   NSIG is an integer constant (input), a parameter for the stopping criterion of
    !!--..   the local search algorithm.  If the value of the objective function is the same
    !!--..   in the last two steps in the first NSIG significant digits, then the local
    !!--..   search will be cancelled.  The value of NSIG will not be changed during the run
    !!--..   of GLOBAL.
    !!--..
    !!--..   X0 is a two-dimensional array of 15 * 20 real elements (output).  After the run
    !!--..   of GLOBAL this matrix contains the parameters of the local minimizer points
    !!--..   found: the J-th local minimizer vector is stored in X0(I,J) I=1,2,...,NPARM.
    !!--..   These minimizers are ordered in ascending order according to their objective
    !!--..   function values.  X0 will be set by GLOBAL.
    !!--..
    !!--..   NC is an integer constant (output) containing the number of the local minima
    !!--..   found. This constant will be set by the GLOBAL routine.
    !!--..
    !!--..   F0 is a real array of 20 elements (output) that contains the objective
    !!--..   function values of the respective local minimizer points.  Thus F(X0(.,J)) =
    !!--..   F0(J).  F0 will be set by GLOBAL.
    !!--..
    !!--..   III. THE FOLLOWING SUBROUTINES ARE CALLED BY THE SUBROUTINE GLOBAL:
    !!--..   -------------------------------------------------------------------
    !!--..
    !!--..   URDMN, FUN, LOCAL, UPDATE, Model_Functn, TIMER
    !!--..
    !!--..   The user has to provide only the routine Model_Functn (see later) out of the
    !!--..   mentioned.  URDMN generates the uniformly distributed pseudorandom numbers.  If
    !!--..   necessary, it can be replaced with any other routine (in case the latter
    !!--..   fulfills the calling requirements).
    !!--..
    !!--..   The starting value for the random number generation is given by the TIMER
    !!--..   routine.  This starting value is based on the computer clock, hence it is
    !!--..   expected to be different for all runs.  The routine TIMER is written in the
    !!--..   ASSEMBLER language of the IBM-PC (MASM).  It can be changed to another routine
    !!--..   that asks the user what should be the starting number - if it is preferred.
    !!--..
    !!--..   FUN transforms the scaled variables back before calling the Model_Functn routine.  It
    !!--..   is provided together with the GLOBAL routine.  Thus, the user can write the
    !!--..   Model_Functn routine according to the usual parameter space.
    !!--..
    !!--..   The subroutines LOCAL and UPDATE contain the Quasi-Newton local search method.
    !!--..   If a random search algorithm is preferred, the routine UNIRANDI should be
    !!--..   linked, it has the same calling form as LOCAL.
    !!--..
    !!--..   IV. THE SUBROUTINE Model_Functn FOR COMPUTING THE OBJECTIVE FUNCTION
    !!--..   --------------------------------------------------------------------
    !!--..
    !!--..   The calling sequence is:
    !!--..
    !!--..                     CALL Model_Functn (NPARM,X, F)
    !!--..
    !!--..   where the parameters should be understood as follows:
    !!--..
    !!--..   X is a real vector of NPARM elements (input) containing the values of the
    !!--..   parameters, i.e. the co-ordinates of the point in the parameter space.  The
    !!--..   variables X(I) are now not scaled.  (They are transformed back to the interval
    !!--..   given by the user.)  The vector X must not be changed in the routine Model_Functn.
    !!--..
    !!--..   F is a real constant (output).  This variable has to contain the computed
    !!--..   objective function value after returning from the routine FUNCT.
    !!--..
    !!--..   If there are more necessary information besides the values of X, and NPARM
    !!--..   for the routine Model_Functn, then they can be passed through global module variables.
    !!--..
    !!--..
    !!--..   V. LIMITATIONS
    !!--..   --------------
    !!--..
    !!--..   parameter           limitation               otherwise
    !!--..
    !!--..   ------------------------------------------------------------------------------
    !!--..
    !!--..   MIN(I)              MIN(I) = MAX(I)          error message, interrupt
    !!--..
    !!--..   MAX(I)              MIN(I) = MAX(I)          error message, interrupt
    !!--..
    !!--..   NPARM               1 <= NPARM <= 15         error message, interrupt
    !!--..
    !!--..   NSAMPL              20 <= NSAMPL <= 10000    the value of NSAMPL is changed to
    !!--..                                                20 or 10000 whichever is closer
    !!--..                                                to the previous value
    !!--..
    !!--..   NSEL                1 <= NSEL <= 20          the value of NSEL is changed to 1
    !!--..                                                or 20 whichever is closer to the
    !!--..                                                previous value
    !!--..
    !!--..   VI. INPUT - OUTPUT
    !!--..   ------------------
    !!--..
    !!--..   Input is a reponsability of the user. The program writes to the
    !!--..   given logical device (specified by IPR).  This output gives a document of
    !!--..   the run, and provides a list of the local minima found.
    !!--..
    !!--..
    !!--..   VII. THE SUGGESTED VALUES FOR THE PARAMETERS OF THE ALGORITHM
    !!--..   -------------------------------------------------------------
    !!--..
    !!--..   The filling out of the arrays MIN and MAX does not mean generally any problem.
    !!--..   The program accepts the input parameters even if MAX(I) is less than MIN(I).
    !!--..   However, if MIN(I) = MAX(I) the program halts with an error message.  According
    !!--..   to this description of the algorithm, the sampling points are generated in the
    !!--..   NPARM-dimensional interval given by vectors MIN and MAX.  However, the local
    !!--..   search may leave this interval.  If the given limits are ment as strict limits,
    !!--..   the objective function should be modified to incorporate some form of penalty
    !!--..   for points outside the interval.
    !!--..
    !!--..   NSAMPL and NSEL can affect the reliability of the GLOBAL subroutine.  If NSAMPL
    !!--..   = NSEL, then the subroutine corresponds to the so-called Multiply Starts
    !!--..   method (usually it is not sufficiently efficient).  It is worth to choose
    !!--..   NSAMPL to be at least as large as the number of function evaluations used by a
    !!--..   single local search. The smaller the NSEL/NSAMPL ratio, the smaller can be the
    !!--..   region of attraction of the global minimum.  Thus, decreasing this ratio will
    !!--..   cause the GLOBAL routine to be more reliable.
    !!--..
    !!--..   It is not worth to give small value for NSIG.  Take into account that the local
    !!--..   search procedure is capable to determine the location of the local minimum
    !!--..   only to that extent what is allowed by the objective function.  As a minimum, 6
    !!--..   is suggested for a value of NSIG. The clustering phase of the GLOBAL routine
    !!--..   is more effective if the local minima are well approximated by the local
    !!--..   search procedure.
    !!--..
    !!--..
    !!--..   VIII. SAMPLE PROGRAM TO SHOW THE USAGE OF THE GLOBAL SUBROUTINE ON PC-S
    !!--..   -----------------------------------------------------------------------
    !!--..
    !!--..         REAL X0(15,20), F0(20), MIN(15), MAX(15)
    !!--..
    !!--..         M = 1
    !!--..         NPARM = 1
    !!--..         NSAMPL = 100
    !!--..         NSEL = 2
    !!--..         IPR = 6
    !!--..
    !!--..         OPEN(6, FILE='OUTPUT')
    !!--..
    !!--..         NSIG = 6
    !!--..         MIN(1) = -100.0
    !!--..         MAX(1) =  100.0
    !!--..
    !!--..         CALL GLOBAL(MIN, MAX, NPARM, NSAMPL, NSEL, IPR, NSIG, X0, NC, F0)
    !!--..
    !!--..         END
    !!--..
    !!--..
    !!--..
    !!--..         SUBROUTINE Model_Functn(NPARM, X, VALUE)
    !!--..
    !!--..         DIMENSION X(:)
    !!--..
    !!--..         VALUE = 1.0 - COS(X(1)) + (X(1)/100.0)**2
    !!--..
    !!--..         RETURN
    !!--..         SUBROUTINE Model_Functn
    !!--..
    !!--..   The given test example is a one-dimensional problem, it has several local
    !!--..   minima.  The global minimum is 0, and the objective function reaches this value
    !!--..   at the point x(1) = 0.0. The test program found 5 local minima in 48 sec. -
    !!--..   among them the global one - with 523 function evaluations.  Note, that if the
    !!--..   NSEL/NSAMPL ratio is too small, then the local minima with relatively large
    !!--..   objective function values cannot be found by the program.  This feature is
    !!--..   useful if we are interested mainly only in the global minimum.
    !!--..
    !!--..   Dr. Tibor Csendes,
    !!--..   Institute of Informatics, Jozsef Attila University
    !!--..   H-6701 Szeged, Pf. 652, Hungary
    !!--..   E-mail: csendes@inf.u-szeged.hu
    !!--..
    !!--..   URL: http://www.inf.u-szeged.hu/~csendes/
    !!--..
    !!--..   You can find additional important information in the comments of the source codes.
    !!--..
    !!--..
    !!--..
    !!--..   REFERENCES
    !!--..   ----------
    !!--..
    !!--..   1. Boender, C.G.E., A.H.G. Rinnooy Kan, G.T. Timmer, L. Stougie: A stochastic
    !!--..      method for global optimization, Mathematical Programming 22(1982) 125-140.
    !!--..
    !!--..   2. Gill, P.E., W. Murray, M.H. Wright: Practical Optimization (Academic Press,
    !!--..      London, 1981).
    !!--..
    !!--..   3. Csendes, T., B. Daroczy, Z. Hantos: Nonlinear parameter estimation by
    !!--..      global optimization: comparison of local search methods in respiratory
    !!--..      system modelling, in: System Modelling and Optimization (Springer-Verlag,
    !!--..      Berlin, 1986) 188-192.
    !!--..
    !!--..   4. Csendes, T.: Nonlinear parameter estimation by global optimization -
    !!--..      Efficiency and reliability, Acta Cybernetica 8(1988) 361-370.
    !!--..
    !!--..   5. Jarvi, T.: A random search optimizer with an application to a maxmin
    !!--..      problem, Publications of the Inst. of Appl. Math., Univ. of Turku, No. 3,
    !!--..      1973.
    !!--..
    !!--..   6. Timmer, G.T.: Global optimization: a stochastic approach, (Ph.D. Thesis,
    !!--..      Erasmus University, Rotterdam, 1984).
    !!--..
    !!---- Update: February - 2005
    !!
    Subroutine Csendes_Global(Model_Functn, Mini, Maxi, Nparm, N100, Ng0, Nsig, X0, Nc, F0, Ipr)
       !---- Arguments ----!
       real(kind=cp), dimension(:),  intent(in)      :: mini
       real(kind=cp), dimension(:),  intent(in)      :: maxi
       integer,                      intent(in)      :: nparm
       integer,                      intent(in out)  :: n100
       integer,                      intent(in out)  :: ng0
       integer,                      intent(in)      :: nsig
       real(kind=cp), dimension(:,:),intent(in out)  :: x0
       integer,                      intent(out)     :: nc
       real(kind=cp), dimension(:),  intent(in out)  :: f0
       integer,                      intent(in)      :: ipr

       Interface
          Subroutine Model_Functn(Nparm,X, F)
             !---- Arguments ----!
             Use Math_gen, only : cp
             real(kind=cp),dimension(:), intent(in)  :: x
             real(kind=cp),              intent(out) :: f
             integer,                    intent(in)  :: nparm
          End Subroutine Model_Functn
       End Interface

       !---- Local Variables ----!
       logical                  :: new_locmin=.false.
       integer                  :: i, i1, icc, icj, ig, ii, iii, im, in1, inum, inum1, inum2, it, iv, &
                                   j, jj, l1, maxfn, n, n0, n1, ncp, nfe, nfe1, ng, ng10, nm, nn100, ns

       integer, dimension(100)       :: ic
       integer, dimension(size(f0))  :: ic1

       real(kind=cp)                             :: a,  b  , alfa, b1, bb, fc, ff, fm, relcon
       real(kind=cp), dimension(nparm,100)       :: x, xcl
       real(kind=cp), dimension(nparm,size(f0))  :: x1
       real(kind=cp), dimension(100,nparm)       :: r
       real(kind=cp), dimension(nparm)           :: w, y, mmin, mmax
       real(kind=cp), dimension(100)             :: f,  fcl
       real(kind=cp), dimension(size(f0))        :: f1
       real(kind=cp), parameter                  :: zero = 0.0_cp, one = 1.0_cp, two = 2.0_cp, ten = 10.0_cp

       if (nparm <= 0) then
          Err_Optim=.true.
          Err_Mess_Optim="  ERROR: Negative number of parameters! "
          write(unit=ipr,fmt="(a)")   Err_Mess_Optim
          return
       end if

       do i = 1, nparm
          mmin(i) = mini(i)
          mmax(i) = maxi(i)
          if (mmin(i) == mmax(i)) then
             Err_Optim=.true.
             write(unit=Err_Mess_Optim,fmt="(a,i3)")"  ERROR: Minimum=Maximum, for parameter #",i
             write(unit=ipr,fmt="(a)")   Err_Mess_Optim
             return
          end if
       end do

       b1 = one / real(nparm)
       if ( ng0 <  1)     ng0 = 1
       if ( ng0 > 20)     ng0 = 20
       if (n100 < 20)    n100 = 20
       if (n100 > 10000) n100 = 10000
       if (n100 < 100) then
          nn100 = n100
          n = 1
       else
          nn100 = 100
          n = n100 / 100
          n100 = n * 100
       end if
       ng10 = 100
       do i = 1, ng10
          f(i)  = 9.9e10
          ic(i) = 0
          fcl(i)= 0.0
       end do
       do i = 1, nparm
          mmax(i) = (mmax(i)-mmin(i)) / two
          mmin(i) = mmin(i) + mmax(i)
       end do
       alfa = 0.01
       nfe = 0
       ng = 0
       ns = 0
       nc = 0
       ncp = 1
       n0 = 0
       n1 = 0
       im = 1
       ig = 0
       fm = 9.9E10
       maxfn = 500 * nparm
       relcon = ten ** (-nsig)

       !----  SAMPLING ----!
       do_main: DO      !Main infinite loop controlled by variable "it"
          n0 = n0 + n100
          nm = n0 - 1
          ng = ng + ng0
          ns = ns + 1
          if (ns*ng0 > 100) then
             write(unit=ipr,fmt="(a)") " ***   Too many sampling"
             exit do_main
          end if

          b = (one - alfa**(one/real(nm))) ** b1
          bb = 0.1 * b
          do i1 = 1, n
             call urdmn(r,100*nparm)
             do j = 1, nn100
                do i = 1, nparm
                   y(i) = two * r(j,i) - one
                end do
                call fun(y, fc, nparm, mmin, mmax,Model_Functn)
                if (fc < fm) then
                   f(im) = fc
                   do i = 1, nparm
                      x(i,im) = y(i)
                   end do
                   if (im <= ng .and. ic(im) > 0) ig = ig - 1
                   ic(im) = 0
                   im = 1
                   fm = f(1)
                   do i = 2, ng10
                      if (f(i) >= fm) then
                         im = i
                         fm = f(i)
                      end if
                   end do
                end if
             end do
          end do
          nfe = nfe + n100
          write(unit=ipr,fmt="(/,a,i5,a)") " ", n100," Function Evaluations Used For Sampling"

          !---- SORTING ----!
          inum = ng10 - 1
          do i = 1, inum
             im = i
             fm = f(i)
             inum1 = i + 1
             do j = inum1, ng10
                if (f(j) < fm) then
                   im = j
                   fm = f(j)
                end if
             end do
             if (im > i) then
                a = fm
                do j = 1, nparm
                   y(j) = x(j,im)
                end do
                if (i <= ng .and. im > ng) then
                   if (ic(ng) == 0 .and. ic(im)  > 0) ig = ig + 1
                   if (ic(ng)  > 0 .and. ic(im) == 0) ig = ig - 1
                end if
                icc = ic(im)
                inum1 = im - i
                do j = 1, inum1
                   inum2 = im - j
                   f(inum2+1) = f(inum2)
                   ic(inum2+1) = ic(inum2)
                   do jj = 1, nparm
                      x(jj,inum2+1) = x(jj,inum2)
                   end do
                end do
                f(i) = a
                do j = 1, nparm
                   x(j,i) = y(j)
                end do
                ic(i) = icc
             end if
          end do

          if (nc > 0) then
             !---- CLUSTERING TO    X*
             do iii = 1, nc
                i = 1
                in1 = i
                fcl(i) = f0(iii)
                do j = 1, nparm
                   xcl(j,i) = x0(j,iii)
                end do
                do j = 1, ng
                   if (ic(j) == iii) then
                      in1 = in1 + 1
                      xcl(1:nparm,in1) = x(1:nparm,j)
                   end if
                end do

                do    ! while i <= in1
                   do j = 1, ng
                      if (ic(j) == 0) then
                         if (fcl(i) < f(j)) then
                            do l1 = 1, nparm
                               w(l1) = abs(xcl(l1,i)-x(l1,j))
                            end do
                            a = zero
                            do l1 = 1, nparm
                               if (w(l1) > a) a = w(l1)
                            end do
                            if (a < b) then
                               write(unit=ipr,fmt="(a,i2)") " SAMPLE POINT ADDED TO THE CLUSTER NO. ",iii
                               do ii = 1, nparm
                                  w(ii) = x(ii,j) * mmax(ii) + mmin(ii)
                               end do
                               write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") f(j), (w(ii),ii = 1,nparm)
                               ig = ig + 1
                               if (ig >= ng) exit do_main
                               in1 = in1 + 1
                               fcl(in1) = f(j)
                               do ii = 1, nparm
                                  xcl(ii,in1) = x(ii,j)
                               end do
                               ic(j) = iii
                            end if
                         end if
                      end if
                   end do
                   i = i + 1
                   if (i > in1) exit
                end do
             end do

             if (n1 > 0) then
                !---- CLUSTERING TO    X1
                do iii = 1, n1
                   i = 1
                   in1 = i
                   fcl(i) = f1(iii)
                   do j = 1, nparm
                      xcl(j,i) = x1(j,iii)
                   end do
                   do   ! while i <= in1
                      do j = 1, ng
                         if (ic(j) == 0) then
                            if (fcl(i) < f(j)) then
                               do l1 = 1, nparm
                                  w(l1) = abs(xcl(l1,i)-x(l1,j))
                               end do
                               a = zero
                               do l1 = 1, nparm
                                  if (w(l1) > a) a = w(l1)
                               end do
                               if (a < b) then
                                  write(unit=ipr,fmt="(a,i2)") " SAMPLE POINT ADDED TO THE CLUSTER NO. ",ic1(iii)
                                  do ii = 1, nparm
                                     w(ii) = x(ii,j) * mmax(ii) + mmin(ii)
                                  end do
                                  write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") f(j), w(1:nparm)
                                  ig = ig + 1
                                  if (ig >= ng) exit do_main
                                  in1 = in1 + 1
                                  fcl(in1) = f(j)
                                  do ii = 1, nparm
                                     xcl(ii,in1) = x(ii,j)
                                  end do
                                  ic(j) = ic1(iii)
                               end if
                            end if
                         end if
                      end do
                      i = i + 1
                      if (i > in1) exit
                   end do
                end do
             end if
          end if

          !---- LOCAL SEARCH ----!
          it = 0
          do i1 = 1, ng
             if (ic(i1) == 0) then
                y(1:nparm) = x(1:nparm,i1)
                ff = f(i1)
                call local(size(mini), nparm, relcon, maxfn, y, ff, nfe1, r(:,1), mmin, mmax,Model_Functn)

                new_locmin=.true.
                if (nc > 0) then
                   do iv = 1, nc
                      do l1 = 1, nparm
                         w(l1) = abs(x0(l1,iv) - y(l1))
                      end do
                      a = zero
                      do l1 = 1, nparm
                         if (w(l1) > a) a = w(l1)
                      end do
                      if (a < bb) then
                         !---- new seed-point
                         n1 = n1 + 1
                         write(unit=ipr,fmt="(a,i2,a,i5)") " New Seed Point Added To The Cluster No. ",iv,", NFEV=", nfe1
                         do ii = 1, nparm
                            w(ii) = x(ii,i1) * mmax(ii) + mmin(ii)
                         end do
                         write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") ff, (w(ii),ii = 1,nparm)
                         if (ff < f0(iv)) then
                            write(unit=ipr,fmt="(a,i2,2(a,f16.8))")   &
                                 " *** Improvement On The Local Minimum No. ",iv,":", f0(iv), " For ", ff
                            do ii = 1, nparm
                               w(ii) = y(ii) * mmax(ii) + mmin(ii)
                            end do
                            write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") ff, (w(ii),ii = 1,nparm)
                            f0(iv) = ff
                            do ii = 1, nparm
                               x0(ii,iv) = y(ii)
                            end do
                         end if
                         if (n1 > 20) then
                            write(unit=ipr,fmt="(a)") " ***   Too many new seed points"
                            exit do_main
                         end if
                         do ii = 1, nparm
                            x1(ii,n1) = x(ii,i1)
                            xcl(ii,1) = x(ii,i1)
                         end do
                         f1(n1) = f(i1)
                         fcl(1) = f(i1)
                         ic1(n1) = iv
                         icj = iv
                         new_locmin=.false.
                         exit
                      end if
                   end do
                end if

                if (new_locmin) then
                   !---- NEW LOCAL MINIMUM
                   nc = nc + 1
                   ncp = ncp + 1
                   write(unit=ipr,fmt="(a,i2,a,f14.8,a,i5)") " *** The Local Minimum No. ",nc,": ", ff,", NFEV=", nfe1
                   do ii = 1, nparm
                      w(ii) = y(ii) * mmax(ii) + mmin(ii)
                   end do
                   write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") ff, (w(ii),ii = 1,nparm)
                   do ii = 1, nparm
                      x0(ii,nc) = y(ii)
                      xcl(ii,1) = y(ii)
                   end do
                   fcl(1) = ff
                   f0(nc) = ff
                   if (nc >= 20) then
                      write(unit=ipr,fmt="(a)")  " ***   Too many clusters"
                      exit do_main
                   end if
                   it = 1
                   icj = nc
                end if

                !----  CLUSTERING TO THE NEW POINT
                nfe = nfe + nfe1
                ic(i1) = icj
                ig = ig + 1
                if (ig >= ng) exit
                i = 1
                in1 = i

                do   ! while i < in1
                   do j = 1, ng
                      if (ic(j) == 0) then
                         if (fcl(i) < f(j)) then
                            do l1 = 1, nparm
                               w(l1) = abs(xcl(l1,i)-x(l1,j))
                            end do
                            a = zero
                            do l1 = 1, nparm
                               if (w(l1) > a) a = w(l1)
                            end do
                            if (a < b) then
                               in1 = in1 + 1
                               do ii = 1, nparm
                                  xcl(ii,in1) = x(ii,j)
                               end do
                               fcl(in1) = f(j)
                               ic(j) = icj
                               write(unit=ipr,fmt="(a,i2)") " Sample Point Added To The Cluster No. ",icj
                               do ii = 1, nparm
                                  w(ii) = x(ii,j) * mmax(ii) + mmin(ii)
                               end do
                               write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") f(j), (w(ii),ii = 1,nparm)
                               ig = ig + 1
                               if (ig >= ng) exit
                            end if
                         end if
                      end if
                   end do
                   i = i + 1
                   if (i >= in1) exit
                end do
             end if
          end do

          if (it == 0) exit
       end do do_main      !End of main infinite loop controlled by variable "it"

       !----  PRINT RESULTS ----!
       write(unit=ipr,fmt="(/,/,/,a,/,/)") " Local Minima Found:"
       if (nc > 1) then
          inum = nc - 1
          do i = 1, inum
             im = i
             fm = f0(i)
             inum1 = i + 1
             do j = inum1, nc
                if (f0(j) < fm) then
                   im = j
                   fm = f0(j)
                end if
             end do
             if (im > i) then
                a = fm
                do j = 1, nparm
                   y(j) = x0(j,im)
                end do
                inum1 = im - i
                do j = 1, inum1
                   inum2 = im - j
                   f0(inum2+1) = f0(inum2)
                   do jj = 1, nparm
                      x0(jj,inum2+1) = x0(jj,inum2)
                   end do
                end do
                f0(i) = a
                do j = 1, nparm
                   x0(j,i) = y(j)
                end do
             end if
          end do
       end if

       if (nc > 0) then
          do i = 1, nc
             do ii = 1, nparm
                x0(ii,i) = x0(ii,i) * mmax(ii) + mmin(ii)
             end do
             write(unit=ipr,fmt="(tr1,f14.8,3(/,t5, 5(f14.8,tr1)))") f0(i), (x0(ii,i),ii = 1,nparm)
          end do
       end if
       write(unit=ipr,fmt="(a,i5,a)") " Normal Termination After ",Nfe," Function Evaluations"

       return
    End Subroutine Csendes_Global

    !!--++
    !!--++ Subroutine Fun(R, F, Nparm, M, Mini, Maxi, Model_Functn)
    !!--++
    !!--++    (Private)
    !!--++    Subroutine used by Csendes Subroutine
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Fun(R, F, Nparm, Mini, Maxi, Model_Functn)
       !---- Arguments ----!
       real(kind=cp),dimension(:),  intent(in)      :: r
       real(kind=cp),               intent(in out)  :: f
       integer,                     intent(in)      :: nparm
       real(kind=cp),dimension(:),  intent(in)      :: mini
       real(kind=cp),dimension(:),  intent(in)      :: maxi

       Interface
          Subroutine Model_Functn(Nparm,X, F)
             Use Math_gen, only : cp
             real(kind=cp),dimension(:), intent(in)  :: x
             real(kind=cp),              intent(out) :: f
             integer,                    intent(in)  :: nparm
          End Subroutine Model_Functn
       End Interface

       !---- Local Variables ----!
       integer                            :: i
       real(kind=cp), dimension(nparm)    :: x

       do i = 1, nparm
          x(i) = maxi(i) * r(i) + mini(i)
       end do
       call Model_Functn(nparm,x, f)

       return
    End Subroutine Fun

    !!----
    !!---- Subroutine Init_Err_Optim()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Optim()

       err_optim=.false.
       err_mess_optim=" "

       return
    End Subroutine Init_Err_Optim

    !!----
    !!---- Subroutine Init_Opt_Contitions(Opt)
    !!----    type(Opt_Conditions_Type), intent(out) :: Opt
    !!----
    !!----    Initialize the variable Opt
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Opt_Conditions(Opt)
       !---- Arguments ----!
       type(Opt_Conditions_Type), intent(in out) :: Opt
       Opt=Opt_Conditions_Type(0,0,0,0,0,0,0,0,1.0e-6_cp,1.0e-20_cp)

       return
    End Subroutine Init_Opt_Conditions

    !!--++
    !!--++ Subroutine Local (Siz,N,Eps,Maxfn,X,F,Nfev,W,Mmin,Mmax)
    !!--++  integer,            intent(in)      :: siz     ! Physical size of the parameter vector
    !!--++  integer,            intent(in)      :: n       ! The Number Of Parameters
    !!--++  real,               intent(in)      :: eps     ! Convergence Criterion
    !!--++  integer,            intent(in out)  :: maxfn   ! Maximum Number Of Function Evaluations
    !!--++  real, dimension(n), intent(in out)  :: x       ! Vector Of Length N Containing Parameter Values
    !!--++                                                 ! In -> Must Contain The Initial Parameter Estimates
    !!--++                                                 ! Out-> The Final Parameter Estimates As Determined By Local
    !!--++  real,               intent(out   )  :: f       ! The Value Of The Function At The Final Parameter Estimates
    !!--++  integer,            intent(out)     :: nfev    ! The Number Of Function Evaluations
    !!--++  real, dimension(:), intent(out)     :: w       ! A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
    !!--++  real, dimension(:), intent(in)      :: mmin    ! Lower bounds of the parameters
    !!--++  real, dimension(:), intent(in)      :: mmax    ! Upper bounds of the parameters
    !!--++
    !!--++  Interface
    !!--++     Subroutine Model_Functn(nparm,x, f)
    !!--++        Use Math_gen, only : cp
    !!--++        real(kind=cp),dimension(:), intent(in)  :: x
    !!--++        real(kind=cp),              intent(out) :: f
    !!--++        integer,                    intent(in)  :: nparm
    !!--++     End Subroutine Model_Functn
    !!--++  End Interface
    !!--++
    !!--++
    !!--++  MINIMUM OF A FUNCTION OF N VARIABLES USING A QUASI-NEWTON METHOD
    !!--++
    !!--..  Information
    !!--..
    !!--..  Convergence criterion: this convergence condition is satisfied if
    !!--..  on two successive iterations, the parameter estimates (i.e.,
    !!--..  x(i), i=1,...,n) differs, component by component, by at most eps.
    !!--..
    !!--..  Searching Limits: x(i) is searched in the interval (mmin(i),mmax(i))
    !!--..
    !!--..  Defining a Subroutine FUN
    !!--..
    !!--..  a user supplied subroutine which calculates the function f for given
    !!--..  parameter values x(1),x(2),...,x(n).
    !!--..  the calling sequence has the following form
    !!--..
    !!--..      CALL FUN(X, F, N, MMIN, MMAX,Model_Functn)
    !!--..
    !!--..  where x is a vector of length n. Model_Functn must have an explicit
    !!--..  interface in the calling program. Model_Functn must not alter the values
    !!--..  of x(i),i=1,...,n or n.
    !!--..
    !!--++ Update: February - 2005
    !!
    Subroutine Local(Siz, N, Eps, Maxfn, X, F, Nfev, W, Mmin, Mmax, Model_Functn)
       !---- Arguments ----!
       integer,                     intent(in)      :: siz
       integer,                     intent(in)      :: n
       real(kind=cp),               intent(in)      :: eps
       integer,                     intent(in out)  :: maxfn
       real(kind=cp), dimension(n), intent(in out)  :: x
       real(kind=cp),               intent(out)     :: f
       integer,                     intent(out)     :: nfev
       real(kind=cp), dimension(:), intent(out)     :: w
       real(kind=cp), dimension(:), intent(in)      :: mmin
       real(kind=cp), dimension(:), intent(in)      :: mmax

       Interface
          Subroutine Model_Functn(nparm,x, f)
             Use Math_gen, only : cp
             real(kind=cp),dimension(:), intent(in)  :: x
             real(kind=cp),              intent(out) :: f
             integer,                    intent(in)  :: nparm
          End Subroutine Model_Functn
       End Interface

       !---- Local Variables ----!
       logical              :: from290=.false.
       integer              :: ig, igg, is, idiff, ir, ij, i, iopt, j, nm1, jj, jp1, l, kj, k,  &
                               link, itn, ii, im1, jnt, np1, jb, nj, ier
       real(kind=cp)        :: hh, hjj, v, df, relx, gs0, diff, aeps, alpha, ff, tot, f1, f2,  &
                               z, gys, dgs, sig, zz, hhh, ghh
       real(kind=cp), dimension(siz)           :: g
       real(kind=cp), dimension(siz*(siz+1)/2) :: h
       real(kind=cp), parameter      :: reps = 1.1921e-07_cp, zero = 0.0_cp, one = 1.0_cp, half = 0.5_cp, &
                                        seven = 7.0_cp, five = 5.0_cp, twelve = 12.0_cp, p1 = 0.1_cp

       !---- Initialization ----!

       !---- IOPT: OPTIONS SELECTOR.
       !---- = 0 Causes LOCAL to initialize the hessian matrix H to the identity matrix.
       !---- = 1 Indicates that H has been initialized by the user to a positive definite matrix.
       !---- = 2 Causes LOCAL to compute the diagonal values of the hessian matrix and set H to
       !----     a diagonal matrix containing these values.
       !---- = 3 Causes LOCAL to compute an estimate of the hessian in H.
       iopt = 0
       ier = 0
       hh = sqrt(reps)
       ig = n
       igg = n + n
       is = igg
       idiff = 1
       ir = n
       w(1) = -one
       w(2) = zero
       w(3) = zero
       h=zero

       !----  Evaluate function at starting point
       g(1:n) = x(1:n)
       call fun(g, f, n, mmin, mmax,Model_Functn)
       nfev = 1
       do
          if (iopt /= 1) then
             !---- set off-diagonal elements of h to 0.0
             if (n /= 1) then
                ij = 2
                do i = 2, n
                   do j = 2, i
                      h(ij) = zero
                      ij = ij + 1
                   end do
                   ij = ij + 1
                end do

                if (iopt == 0) then
                   !---- set diagonal elements of h to one
                   ij = 0
                   do i = 1, n
                      ij = ij + i
                      h(ij) = one
                   end do
                   exit     !from infinite loop
                end if
             end if

             !---- get diagonal elements of hessian
             im1 = 1
             nm1 = 1
             np1 = n + 1
             do i = 2, np1
                hhh = hh * abs(x(im1))
                g(im1) = x(im1) + hhh
                call fun(g, f2, n, mmin, mmax,Model_Functn)
                g(im1) = g(im1) + hhh
                call fun(g, ff, n, mmin, mmax,Model_Functn)
                h(nm1) = (ff-f2+f-f2) / (hhh*hhh)
                g(im1) = x(im1)
                im1 = i
                nm1 = i + nm1
             end do
             nfev = nfev + n + n

             if (iopt == 3 .and. n /= 1) then
                !---- get the rest of the hessian
                jj = 1
                ii = 2
                do i = 2, n
                   ghh = hh * abs(x(i))
                   g(i) = x(i) + ghh
                   call fun(g, f2, n, mmin, mmax,Model_Functn)
                   do j = 1, jj
                      hhh = hh * abs(x(j))
                      g(j) = x(j) + hhh
                      call fun(g, ff, n, mmin, mmax,Model_Functn)
                      g(i) = x(i)
                      call fun(g, f1, n, mmin, mmax,Model_Functn)
                      h(ii) = (ff-f1-f2+f) / (hhh*ghh)
                      ii = ii + 1
                      g(j) = x(j)
                   end do
                   jj = jj + 1
                   ii = ii + 1
                end do
                nfev = nfev + ((n*n-n)/2)
             end if
          end if

          !---- Factor H to L*D*L-TRANSPOSE
          ir = n
          if (n <= 1) then
             if (h(1) > zero) exit  !exit from infinite loop
             h(1) = zero
             ir = 0
          else
             nm1 = n - 1
             jj = 0
             do j = 1, n
                jp1 = j + 1
                jj = jj + j
                hjj = h(jj)
                if (hjj <= zero) then
                   h(jj) = zero
                   ir = ir - 1
                else
                   if (j /= n) then
                      ij = jj
                      l = 0
                      do i = jp1, n
                         l = l + 1
                         ij = ij + i - 1
                         v = h(ij) / hjj
                         kj = ij
                         do k = i, n
                            h(kj+l) = h(kj+l) - h(kj) * v
                            kj = kj + k
                         end do
                         h(ij) = v
                      end do
                   end if
                end if
             end do
          end if

          if (ir /= n) then
             ier = 129
             return
          end if
          exit
       end do

       itn = 0
       df = -one
       do_120: DO
          !---- Evaluate gradienT W(IG+I),I=1,...,N
          link = 1

          !---- Evaluate gradient
          do_410: DO
             If (idiff /= 2) then
                !---- forward differences
                !---- gradient =    w(ig+i), i=1,...,n
                do i = 1, n
                   z = hh * abs(x(i))
                   zz = x(i)
                   x(i) = zz + z
                   call fun(x, f1, n, mmin, mmax,Model_Functn)
                   w(ig+i) = (f1-f) / z
                   x(i) = zz
                end do
                nfev = nfev + n

             else ! now idiff == 2
                !---- central differences
                !---- gradient =    w(ig+i), i=1,...,n
                do i = 1, n
                   z = hh * abs(x(i))
                   zz = x(i)
                   x(i) = zz + z
                   call fun(x, f1, n, mmin, mmax,Model_Functn)
                   x(i) = zz - z
                   call fun(x, f2, n, mmin, mmax,Model_Functn)
                   w(ig+i) = (f1-f2) / (z+z)
                   x(i) = zz
                end do
                nfev = nfev + n + n
             end if

             if ( link == 2 ) then
                do  ! 1-iteration loop
                    if (nfev < maxfn) then
                       gys = zero
                       do i = 1, n
                          gys = gys + w(ig+i) * w(is+i)
                          w(igg+i) = w(i)
                       end do
                       df = ff - f
                       dgs = gys - gs0
                       if (dgs <= zero) exit ! 1-iteration loop
                       if (dgs+alpha*gs0 <= zero) then
                          !---- update hessian h using complementary dfp formula
                          sig = one / gs0
                          ir = -ir
                          call update(h, n, w, sig, g, ir, 0, zero)
                          do i = 1, n
                             g(i) = w(ig+i) - w(igg+i)
                          end do
                          sig = one / (alpha*dgs)
                          ir = -ir
                          call update(h, n, g, sig, w, ir, 0, zero)
                          exit ! 1-iteration loop
                       end if

                       !---- update hessian using dfp formula
                       zz = alpha / (dgs-alpha*gs0)
                       sig = -zz
                       call update(h, n, w, sig, g, ir, 0, reps)
                       z = dgs * zz - one
                       do i = 1, n
                          g(i) = w(ig+i) + z * w(igg+i)
                       end do
                       sig = one / (zz*dgs*dgs)
                       call update(h, n, g, sig, w, ir, 0, zero)
                    end if
                    exit
                end do !1-iteration loop
             end if

             !---- Begin iteration loop
             if (nfev >= maxfn) exit do_410
             itn = itn + 1
             do i = 1, n
                w(i) = -w(ig+i)
             end do

             !---- Determine search direction W by solving    H*W = -G
             !---- Where H = L*D*L-TRANSPOSE
             if (ir >= n) then
                !---- n .eq. 1
                g(1) = w(1)
                if (n <= 1) then
                   w(1) = w(1) / h(1)
                else
                   !---- n .gt. 1
                   ii = 1
                   !---- solve l*w = -g
                   do i = 2, n
                      ij = ii
                      ii = ii + i
                      v = w(i)
                      im1 = i - 1
                      do j = 1, im1
                         ij = ij + 1
                         v = v - h(ij) * w(j)
                      end do
                      g(i) = v
                      w(i) = v
                   end do
                   !---- solve (d*lt)*z = w where lt = l-transpose
                   w(n) = w(n) / h(ii)
                   jj = ii
                   nm1 = n - 1
                   do nj = 1, nm1
                      !---- j = n-1,n-2,...,1
                      j = n - nj
                      jp1 = j + 1
                      jj = jj - jp1
                      v = w(j) / h(jj)
                      ij = jj
                      do i = jp1, n
                         ij = ij + i - 1
                         v = v - h(ij) * w(i)
                      end do
                      w(j) = v
                   end do
                end if
             end if

             !---- Determine step length ALPHA
             relx = zero
             gs0 = zero
             do i = 1, n
                w(is+i) = w(i)
                diff = abs(w(i)) / abs(x(i))
                relx = max(relx, diff)
                gs0 = gs0 + w(ig+i) * w(i)
             end do
             if (relx == zero) then
                if (idiff /= 2) then
                   !---- change to central differences
                   idiff = 2
                   cycle do_120
                else
                   exit do_410
                end if
             end if
             aeps = eps / relx
             ier = 130

             if (gs0 >= zero) then
                if (idiff /= 2) then
                   !---- change to central differences
                   idiff = 2
                   cycle do_120
                else
                   exit do_410
                end if
             end if
             if (df == zero) then
                if (idiff /= 2) then
                   !---- change to central differences
                   idiff = 2
                   cycle do_120
                else
                   exit do_410
                end if
             end if

             ier = 0
             alpha = (-df-df) / gs0
             if (alpha <= zero) alpha = one
             alpha = min(alpha, one)
             if (idiff == 2) alpha = max(p1,alpha)
             ff = f
             tot = zero
             jnt = 0

             !---- Search along    X + ALPHA*W
             do_200: DO
                if (nfev >= maxfn) exit do_410
                do i = 1, n
                   w(i) = x(i) + alpha * w(is+i)
                end do
                call fun(w, f1, n, mmin, mmax,Model_Functn)
                nfev = nfev + 1
                if (f1 < f) then
                   f2 = f
                   tot = tot + alpha
                   do
                      ier = 0
                      f = f1
                      do i = 1, n
                         x(i) = w(i)
                      end do
                      if (jnt-1 < 0) then
                         if (nfev >= maxfn) exit do_410
                         do i = 1, n
                            w(i) = x(i) + alpha * w(is+i)
                         end do
                         call fun(w, f1, n, mmin, mmax,Model_Functn)
                         nfev = nfev + 1
                         if (f1 >= f) then
                            from290=.true.
                            exit do_200
                         end if
                         if (f1+f2 >= f+f .and. seven*f1+five*f2 > twelve*f) jnt = 2
                         tot = tot + alpha
                         alpha = alpha + alpha
                      else if (jnt-1 == 0) then
                         exit do_200
                      else
                         from290=.true.
                         exit do_200
                      end if
                   end do
                end if

                if (f == ff .and. idiff == 2 .and. relx > eps) ier = 130
                if (alpha < aeps) then
                   if (idiff /= 2) then
                      !---- change to central differences
                      idiff = 2
                      cycle do_120
                   else
                      exit do_410
                   end if
                end if
                if (nfev >= maxfn) exit do_410
                alpha = half * alpha
                do i = 1, n
                   w(i) = x(i) + alpha * w(is+i)
                end do
                call fun(w, f2, n, mmin, mmax,Model_Functn)
                nfev = nfev + 1
                if (f2 < f) then
                   tot = tot + alpha
                   ier = 0
                   f = f2
                   do i = 1, n
                      x(i) = w(i)
                   end do
                else
                   z = p1
                   if (f1+f > f2+f2) z = one + half * (f-f1) / (f+f1-f2-f2)
                   z = max(p1,z)
                   alpha = z * alpha
                   jnt = 1
                   cycle do_200
                end if
                exit
             end do do_200

             if (.not. from290) then
                if (tot < aeps) then
                   if (idiff /= 2) then
                      !---- change to central differences
                      idiff = 2
                      cycle do_120
                   else
                      exit do_410
                   end if
                end if
             end if
             from290=.false.
             alpha = tot

             !---- Save old gradient
             do i = 1, n
                w(i) = w(ig+i)
             end do

             !---- Evaluate gradient W(IG+I), I=1,...,N
             link = 2
          end do do_410

          !---- MAXFN function evaluations
          if (relx > eps .and. ier == 0) cycle do_120
          exit
       end do do_120

       !---- Compute H = L*D*L-TRANSPOSE and output
       if (n == 1) return
       np1 = n + 1
       nm1 = n - 1
       jj = (n*(np1)) / 2
       do jb = 1, nm1
          jp1 = np1 - jb
          jj = jj - jp1
          hjj = h(jj)
          ij = jj
          l = 0
          do i = jp1, n
             l = l + 1
             ij = ij + i - 1
             v = h(ij) * hjj
             kj = ij
             do k = i, n
                h(kj+l) = h(kj+l) + h(kj) * v
                kj = kj + k
             end do
             h(ij) = v
          end do
          hjj = h(jj)
       end do

       return
    End Subroutine Local

    !!----
    !!---- Subroutine Nelder_Mead_Simplex(Model_Functn, Nop, P, Step, Var, Func, C, Ipr)
    !!----    integer,                      intent(in)      :: nop   ! In -> no. of parameters, incl. any to be held fixed
    !!----    real (kind=cp), dimension(:), intent(in out)  :: p     ! In -> starting values of parameters
    !!----                                                           ! Out-> final values of parameters
    !!----    real (kind=cp), dimension(:), intent(in out)  :: step  ! In -> initial step sizes
    !!----    real (kind=cp), dimension(:), intent(out)     :: var   ! Out-> contains the diagonal elements of the inverse of
    !!----                                                           !       the information matrix.
    !!----    real (kind=cp),               intent(out)     :: func  ! Out-> the function value corresponding to the final
    !!----                                                           !       parameter values.
    !!----    type(opt_conditions),         intent(in out)  :: c     ! Optimization Conditions
    !!----    integer, optional,            intent(in)      :: Ipr   ! Logical unit for printing if the parameter C%IOUT /= 0.
    !!----
    !!--<<    Interface
    !!----       Subroutine Model_Functn(n,p, func)                  ! name of the user's subroutine - arguments (P,FUNC)
    !!----          use math_gen, only: cp                           ! which returns the function value for a given set of
    !!----          integer,                      intent(in)  :: n   ! Number of parameters
    !!----          real (kind=cp), dimension(:), intent(in)  :: p   ! parameter values in array P.
    !!----          real (kind=cp),               intent(out) :: func
    !!----       End Subroutine Model_Functn
    !!-->>    End Interface
    !!----
    !!----    A subroutine for function minimization using the SIMPLEX method.
    !!----
    !!--<<
    !!----    Optimization Conditions type with the following components:
    !!----       C%MXFUN = Input, the maximum no. of function evaluations allowed.
    !!----                   Say, 20 times the number of parameters, NOP.
    !!----       C%IOUT  = Input, print control parameter
    !!----                   < 0 No printing
    !!----                   = 0 Printing of parameter values and the function
    !!----                       value after initial evidence of convergence.
    !!----                   > 0 As for C%IOUT = 0 plus progress reports after every
    !!----                       C%IOUT evaluations, plus printing for the initial simplex.
    !!----       C%EPS   = Input, stopping criterion.
    !!----                 The criterion is applied to the standard deviation of
    !!----                 the values of FUNC at the points of the simplex.
    !!----       C%Loops = Input, the stopping rule is applied after every NLOOP
    !!----                 function evaluations.   Normally NLOOP should be slightly
    !!----                 greater than NOP, say NLOOP = 2*NOP.
    !!----       C%Iquad = Input, = 1 If fitting of a quadratic surface is required
    !!----                        = 0 If not
    !!----                 N.B. The fitting of a quadratic surface is strongly
    !!----                 recommended, provided that the fitted function is
    !!----                 continuous in the vicinity of the minimum.   It is often
    !!----                 a good indicator of whether a premature termination of
    !!----                 the search has occurred.
    !!----       C%ACC   = Input, criterion for expanding the simplex to overcome
    !!----                 rounding errors before fitting the quadratic surface.
    !!----                 The simplex is expanded so that the function values at
    !!----                 the points of the simplex exceed those at the supposed
    !!----                 minimum by at least an amount SIMP.
    !!----
    !!----       C%NFLAG = Output, = 0 for successful termination
    !!----                   = 1 If maximum no. of function evaluations exceeded
    !!----                   = 2 If information matrix is not +ve semi-definite
    !!----                   = 3 If NOP < 1
    !!----                   = 4 If C%Loops < 1
    !!----
    !!----       N.B. P, STEP and VAR (If C%Iquad = 1) must have dimension at least NOP
    !!----            in the calling program.
    !!-->>
    !!--..    For details, see Nelder & Mead, The Computer JournaL, January 1965
    !!--..    Programmed by D.E.Shaw,
    !!--..    CSIRO, Division of Mathematics & Statistics
    !!--..    P.O. BOX 218, Lindfield, N.S.W. 2070
    !!--..
    !!--..    With amendments by R.W.M.WEDDERBURN
    !!--..    Rothamsted Experimental Station
    !!--..    Harpenden, Hertfordshire, ENGLAND
    !!--..
    !!--..    Further amended by Alan Miller
    !!--..    CSIRO Division of Mathematical & Information Sciences
    !!--..    Private Bag 10, CLAYTON, VIC. 3169
    !!--..
    !!--..    Fortran 90 conversion by Alan Miller, June 1995
    !!--..    Alan.Miller @ vic.cmis.csiro.au
    !!--..    Latest revision - 5 December 1999
    !!--..
    !!--..    Conversion to F-language and more modifications
    !!--..    by Juan Rodriguez-Carvajal (LLB-CEA)- 7 April 2004
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Nelder_Mead_Simplex(Model_Functn, Nop, P, Step, Var, Func, C, Ipr)
       !---- Arguments ----!
       integer,                      intent(in)      :: nop
       real(kind=cp), dimension(:),  intent(in out)  :: p, step
       real(kind=cp), dimension(:),  intent(out)     :: var
       real(kind=cp),                intent(out)     :: func
       type(opt_conditions_Type),    intent(in out)  :: c
       integer, optional,            intent(in)      :: Ipr

       Interface
          Subroutine Model_Functn(n,P,Func)
             use math_gen, only: cp
             integer,                      intent(in)  :: n
             real (kind=cp), dimension(:), intent(in)  :: p
             real (kind=cp),               intent(out) :: func
          End Subroutine Model_Functn
       End Interface

       !---- Local variables ----!
       real(kind=cp), dimension(nop+1,nop)     :: g
       real(kind=cp), dimension(nop+1)         :: h
       real(kind=cp), dimension(nop)           :: pbar, pstar, pstst,aval, pmin
       real(kind=cp), dimension(nop*(nop+1)/2) :: bmat, vc

       real(kind=cp)                           :: ymin, rmax, hstst, a0, hmin, test, hmean, &
                                                  hstd, hstar, hmax, savemn

       real(kind=cp), parameter                :: zero = 0.0_cp, half = 0.5_cp, one = 1.0_cp, two = 2.0_cp

       integer                                 :: i, i1, i2, iflag, ii, ij, imax, imin, irank, irow, j, j1, jj, &
                                                  k, l, loop, nap, neval, nmore, np1, nullty

       !---- A = Reflection coefficient, B = Contraction coefficient, and
       !---- C = Expansion coefficient.
       real(kind=cp), parameter :: a = 1.0_cp, b = 0.5_cp, cc = 2.0_cp

       !---- If progress reports have been requested, print heading
       if (c%iout > 0) then
          if (present(ipr)) then
             write(unit=ipr,fmt="(a,i4,a,/,a)") " Progress Report every",c%iout, &
                  " function evaluations"," EVAL.   FUNC.VALUE.          PARAMETER VALUES"
          end if
       end if

       !---- Check input arguments
       c%nflag = 0
       if (nop <= 0) c%nflag = 3
       if (c%loops <= 0) c%nflag = 4
       if (c%nflag /= 0) return

       !---- SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0
       nap = count(step /= zero)
       neval = 0
       loop  = 0
       iflag = 0

       !---- IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN
       if (nap <= 0) then
          call Model_Functn(Nop,p,func)
          return
       end if

       !---- SET UP THE INITIAL SIMPLEX
       do          ! main infinite loop setting the initial simplex
          g(1,:) = p
          irow = 2
          do i = 1, nop
             if (step(i) /= zero) then
                g(irow,:) = p
                g(irow,i) = p(i) + step(i)
                irow = irow + 1
             end if
          end do

          np1 = nap + 1
          do i = 1, np1
             p = g(i,:)
             call Model_Functn(Nop,p,h(i))
             neval = neval + 1
             if (c%iout > 0 .and. present(ipr) ) then
                write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") " ",neval, "  ", h(i),"  ", p
             end if
          end do

          !---- START OF MAIN CYCLE.
          !---- FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).
          Main_loop: do
             loop = loop + 1
             imax = 1
             imin = 1
             hmax = h(1)
             hmin = h(1)
             do i = 2, np1
                if (h(i) > hmax) then
                   imax = i
                   hmax = h(i)
                else
                   if (h(i) < hmin) then
                      imin = i
                      hmin = h(i)
                   end if
                end if
             end do

             !---- FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)
             pbar = zero
             do i = 1, np1
                if (i /= imax) then
                   pbar = pbar + g(i,:)
                end if
             end do
             pbar = pbar / nap

             !---- REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
             !---- HSTAR = FUNCTION VALUE AT PSTAR.
             pstar = a * (pbar - g(imax,:)) + pbar
             call Model_Functn(Nop,pstar,hstar)

             neval = neval + 1
             if (c%iout > 0 .and. present(ipr)) then
                if (modulo(neval,c%iout) == 0) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                   " ",neval, "  ", hstar,"  ", pstar
             end if

             !---- IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
             !---- HSTST = FUNCTION VALUE AT PSTST.
             if (hstar < hmin) then
                pstst = cc * (pstar - pbar) + pbar
                call Model_Functn(Nop,pstst,hstst)
                neval = neval + 1
                if (c%iout > 0 .and. present(ipr)) then
                   if (modulo(neval,c%iout) == 0) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                      " ",neval, "  ", hstst, "  ", pstst
                end if

                !---- IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
                !---- HMAX BY HSTST, THEN TEST FOR CONVERGENCE.
                if (hstst >= hmin) then   ! replace maximum point by pstar & h(imax) by hstar.
                   g(imax,:) = pstar
                   h(imax) = hstar
                else
                   g(imax,:) = pstst
                   h(imax) = hstst
                end if

             else   ! (hstar < hmin)
                !---- HSTAR is not < HMIN.
                !---- Test whether it is < function value at some point other than
                !---- P(IMAX).   If it is replace P(IMAX) by PSTAR & HMAX by HSTAR.
                do_250: do
                   do i = 1, np1
                      if (i /= imax) then
                         if (hstar < h(i)) then  ! replace maximum point by pstar & h(imax) by hstar.
                            g(imax,:) = pstar
                            h(imax) = hstar
                            exit do_250
                         end if
                      end if
                   end do

                   !---- HSTAR > all function values except possibly HMAX.
                   !---- If HSTAR <= HMAX, replace P(IMAX) by PSTAR & HMAX by HSTAR.
                   if (hstar <= hmax) then
                      g(imax,:) = pstar
                      hmax = hstar
                      h(imax) = hstar
                   end if

                   !---- Contracted STEP to the point PSTST,
                   !---- HSTST = Function value at PSTST.
                   pstst = b * g(imax,:) + (one-b) * pbar
                   call Model_Functn(Nop,pstst,hstst)
                   neval = neval + 1
                   if (c%iout > 0 .and. present(ipr)) then
                      if (modulo(neval,c%iout) == 0) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                         " ",neval, "  ", hstst, "  ", pstst
                   end if

                   !---- IF HSTST < HMAX replace P(IMAX) by PSTST & HMAX by HSTST.
                   if (hstst <= hmax) then
                      g(imax,:) = pstst
                      h(imax) = hstst
                   else     !(hstst <= hmax)
                      !---- HSTST > HMAX.
                      !---- Shrink the simplex by replacing each point, other than the current
                      !---- minimum, by a point mid-way between its current position and the
                      !---- minimum.
                      do i = 1, np1
                         if (i /= imin) then
                            do j = 1, nop
                               if (step(j) /= zero) g(i,j) = (g(i,j) + g(imin,j)) * half
                               p(j) = g(i,j)
                            end do
                            call Model_Functn(Nop,p,h(i))
                            neval = neval + 1
                            if (c%iout > 0 .and. present(ipr)) then
                               if (modulo(neval,c%iout) == 0) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                                  " ",neval, "  ", h(i), "  ", p
                            end if
                         end if
                      end do
                   end if   !(hstst <= hmax)
                   exit
                end do do_250
             end if  ! (hstar < hmin)

             !---- If LOOP = NLOOP test for convergence, otherwise repeat main cycle.
             if (loop < c%loops) cycle main_loop

             !---- Calculate mean & standard deviation of function values for the
             !---- current simplex.
             hmean = SUM( h(1:np1) ) / np1
             hstd  = SUM( (h(1:np1) - hmean) ** 2 )
             hstd  = SQRT(hstd / np1)

             !---- If the RMS > STOPCR, set IFLAG & LOOP to zero and go to the
             !---- start of the main cycle again.
             if (hstd > c%eps .and. neval <= c%mxfun) then
                iflag = 0
                loop = 0
                cycle main_loop
             end if

             !---- Find the centroid of the current simplex and the function value there.
             do i = 1, nop
                if (step(i) /= zero) then
                   p(i) = sum( g(1:np1,i) ) / np1
                end if
             end do
             call Model_Functn(Nop,p,func)
             neval = neval + 1
             if (c%iout > 0 .and. present(ipr)) then
                if (modulo(neval,c%iout) == 0) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                   " ",neval, "  ", func, "  ", p
             end if

             !---- Test whether the no. of function values allowed, c%mxfun, has been
             !---- overrun; if so, exit with C%NFLAG = 1.
             if (neval > c%mxfun) then
                c%nflag = 1
                if (c%iout < 0) return
                if (present(ipr)) then
                   Write(unit=ipr,fmt="(a,i5)") "  No. of function evaluations > ",c%mxfun
                   Write(unit=ipr,fmt="(a, f14.6)") "  RMS of function values of last simplex =",hstd
                   Write(unit=ipr,fmt="(a,4(/,6f13.5))") "  Centroid of last simplex =",p
                   Write(unit=ipr,fmt="(a, f14.6)") "  Function value at centroid =",func
                end if
                return
             end if

             !---- Convergence criterion satisfied.
             !---- If IFLAG = 0, set IFLAG & save HMEAN.
             !---- If IFLAG = 1 & change in HMEAN <= STOPCR then search is complete.
             if (c%iout >= 0 .and. present(ipr)) then
                Write(unit=ipr,fmt="(/,a)") "  EVIDENCE OF CONVERGENCE"
                Write(unit=ipr,fmt="(a,4(/,6f13.5))") "  Centroid of last simplex =",p
                Write(unit=ipr,fmt="(a, f14.6)") "  Function value at centroid =",func
             end if

             if (iflag == 0 .or. abs(savemn-hmean) >= c%eps) then
                iflag = 1
                savemn = hmean
                loop = 0
             else
                exit main_loop
             end if
          end do main_loop

          if (c%iout >= 0 .and. present(ipr)) then
             Write(unit=ipr,fmt="(/,a,i5,a)") "  Minimum found after ",neval," function evaluations"
             Write(unit=ipr,fmt="(a,4(/,6f13.6))") "  Minimum at",p
             Write(unit=ipr,fmt="(a, f14.6)") "  Function value at minimum =",func
          end if
          if (c%iquad <= 0) return

          !----
          !---- Quadratic surface fitting
          !----
          if (c%iout >= 0 .and. present(ipr)) &
             Write(unit=ipr,fmt="(/,a,/)")  "  Fitting quadratic surface about supposed minimum"

          !---- Expand the final simplex, if necessary, to overcome rounding errors.
          hmin = func
          nmore = 0
          do i = 1, np1
             do
                test = abs(h(i)-func)
                if (test < c%acc) then
                   do j = 1, nop
                      if (step(j) /= zero) g(i,j) = (g(i,j)-p(j)) + g(i,j)
                      pstst(j) = g(i,j)
                   end do
                   call Model_Functn(Nop,pstst,h(i))
                   nmore = nmore + 1
                   neval = neval + 1
                   if (h(i) >= hmin) cycle
                   hmin = h(i)
                   if (c%iout >= 0 .and. present(ipr)) write(unit=ipr,fmt="(a,i4,a,f12.5,a, 5f11.4, 3(/,t22, 5f11.4))") &
                                                                            " ",neval, "  ", hmin, "  ", pstst
                else
                   exit
                end if
             end do
          end do

          !---- Function values are calculated at an additional NAP points.
          do i = 1, nap
             i1 = i + 1
             pstar = (g(1,:) + g(i1,:)) * half
             call Model_Functn(Nop,pstar,aval(i))
             nmore = nmore + 1
             neval = neval + 1
          end do

          !---- The matrix of estimated second derivatives is calculated and its
          !---- lower triangle stored in BMAT.
          a0 = h(1)
          do i = 1, nap
             i1 = i - 1
             i2 = i + 1
             do j = 1, i1
                j1 = j + 1
                pstst = (g(i2,:) + g(j1,:)) * half
                call Model_Functn(Nop,pstst,hstst)
                nmore = nmore + 1
                neval = neval + 1
                l = i * (i-1) / 2 + j
                bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
             end do
          end do

          l = 0
          do i = 1, nap
             i1 = i + 1
             l = l + i
             bmat(l) = two * (h(i1) + a0 - two*aval(i))
          end do

          !---- The vector of estimated first derivatives is calculated and
          !---- stored in aval.
          do i = 1, nap
             i1 = i + 1
             aval(i) = two * aval(i) - (h(i1) + 3.0_cp*a0) * half
          end do

          !---- The matrix Q of Nelder & Mead is calculated and stored in G.
          pmin = g(1,:)
          do i = 1, nap
             i1 = i + 1
             g(i1,:) = g(i1,:) - g(1,:)
          end do

          do i = 1, nap
             i1 = i + 1
             g(i,:) = g(i1,:)
          end do

          !---- invert bmat
          call syminv(bmat, nap, bmat, nullty, c%nflag, rmax)
          if (c%nflag == 0) then
             irank = nap - nullty
             exit !quit the infinite loop
          else                                 ! BMAT not +ve definite
                                               ! Resume search for the minimum
             if (c%iout >= 0 .and. present(ipr)) write(unit=ipr,fmt="(/,a,/,a,/)")    &
                "  MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.", &
                "  MINIMUM PROBABLY NOT FOUND"
             c%nflag = 2
             if (neval > c%mxfun) return
             if (present(ipr)) write(unit=ipr,fmt="(/,t11,a,/)")   "Search restarting"
             step = half * step
             cycle    !
          end if
       end do ! Main infinite loop setting the initial simplex

       !---- BMAT*A/2 IS CALCULATED AND STORED IN H.
       do i = 1, nap
          h(i) = zero
          do j = 1, nap
             if (j <= i) then
                l = i * (i-1) / 2 + j
             else
                l = j * (j-1) / 2 + i
             end if
             h(i) = h(i) + bmat(l) * aval(j)
          end do
       end do

       !---- Find the position, PMIN, & value, YMIN, of the minimum of the
       !---- quadratic.
       ymin = dot_product( h(1:nap), aval(1:nap) )
       ymin = a0 - ymin
       do i = 1, nop
          pstst(i) = dot_product( h(1:nap), g(1:nap,i) )
       end do
       pmin = pmin - pstst
       if (c%iout >= 0 .and. present(ipr)) then
          write(unit=ipr,fmt="(a,f14.6,a,4(/,6f13.5))") "  Minimum of quadratic surface =",ymin," at", pmin
          write(unit=ipr,fmt="(a,/,a,/)")                                        &
               "  If this differs by much from the minimum estimated from the minimization,",          &
               "  The minimum may be false &/or the information matrix may be inaccurate"
       end if

       !---- Q*BMAT*Q'/2 is calculated & its lower triangle stored in VC
       do i = 1, nop
          do j = 1, nap
             h(j) = zero
             do k = 1, nap
                if (k <= j) then
                   l = j * (j-1) / 2 + k
                else
                   l = k * (k-1) / 2 + j
                end if
                h(j) = h(j) + bmat(l) * g(k,i) * half
             end do
          end do

          do j = i, nop
             l = j * (j-1) / 2 + i
             vc(l) = dot_product( h(1:nap), g(1:nap,j) )
          end do
       end do

       !---- The diagonal elements of VC are copied into VAR.
       j = 0
       do i = 1, nop
          j = j + i
          var(i) = vc(j)
       end do
       if (c%iout < 0) return
       if (present(ipr)) then
          write(unit=ipr,fmt="(a,i3,/,a)") "  Rank of information matrix =",irank, &
                                            "  Inverse of information matrix:-"
          call print_tri_matrix(vc, nop, ipr)

          write(unit=ipr,fmt="(5(/,a),/)")                           &
               "  If the function minimized was -LOG(LIKELIHOOD),"   ,&
               "  this is the covariance matrix of the parameters."  ,&
               "  If the function was a sum of squares of residuals,",&
               "  this matrix must be multiplied by twice the estimated residual variance", &
               "  to obtain the covariance matrix."
       end if
       call syminv(vc, nap, bmat, nullty, c%nflag, rmax)

       !---- BMAT NOW CONTAINS THE INFORMATION MATRIX
       if (present(ipr)) then
          write(unit=ipr,fmt="(a,/)")    " INFORMATION MATRIX:-"
          call print_tri_matrix(bmat, nop, ipr)
       end if
       ii = 0
       ij = 0
       do i = 1, nop
          ii = ii + i
          if (vc(ii) > zero) then
             vc(ii) = one / sqrt(vc(ii))
          else
             vc(ii) = zero
          end if
          jj = 0
          do j = 1, i - 1
             jj = jj + j
             ij = ij + 1
             vc(ij) = vc(ij) * vc(ii) * vc(jj)
          end do
          ij = ij + 1
       end do

       if (present(ipr)) then
          write(unit=ipr,fmt="(/,a)")   " CORRELATION MATRIX:-"
          ii = 0
          do i = 1, nop
             ii = ii + i
             if (vc(ii) /= zero) vc(ii) = one
          end do
          call print_tri_matrix(vc, nop, ipr)

          !---- Exit, on successful termination.
          write(unit=ipr,fmt="(/,a,i4,a,/)") " A further", nmore, &
                                               " function evaluations have been used"
       end if

       return
    End Subroutine Nelder_Mead_Simplex

    !!--++
    !!--++ Subroutine Print_Tri_Matrix(A, N, Iunit)
    !!--++    real (kind=cp),dimension(:), intent(in)   :: A       ! Matrix
    !!--++    integer,                     intent(in)   :: N       ! The order of A
    !!--++    integer,                     intent(in)   :: Iunit   ! Output unit
    !!--++
    !!--++    (Private)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Print_Tri_Matrix(A, N, Iunit)
       !---- Arguments ----!
       real(kind=cp),dimension(:),  intent(in)  :: a
       integer,                     intent(in)  :: n
       integer,                     intent(in)  :: iunit

       !---- Local variables ----!
       integer  :: i, ii, i1, i2, l

       l = 1
       do l = 1, n, 6
          ii = l * (l-1) / 2
          do i = l, n
             i1 = ii + l
             ii = ii + i
             i2 = min(ii,i1+5)
             write(unit=iunit,fmt="(tr1,6f13.5)") a(i1:i2)
          end do
       end do

       return
    End Subroutine Print_Tri_Matrix

    !!--++
    !!--++ Subroutine Syminv(A, N, C, Nullty, Ifault, Rmax)
    !!--++   real(kind=cp),dimension(:), intent(in)    :: A       ! the symmetric matrix to be inverted, stored in lower triangular form
    !!--++   integer,                    intent(in)    :: N       ! The order of A
    !!--++   real(kind=cp),dimension(:), intent(out)   :: C       ! the inverse of A (A generalized inverse if C is singular), also stored in lower triangular form.
    !!--++   integer,                    intent(out)   :: Nullty  ! the rank deficiency of A.
    !!--++   integer,                    intent(out)   :: Ifault  ! Error Code: 1 if N < 1, 2 If A Is not +ve semi-definite
    !!--++                                                             !             0 Otherwise
    !!--++   real (kind=cp),             intent(out)   :: Rmax        ! an estimate of the relative accuracy of the diagonal elements of C.
    !!--++
    !!--++    ALGORITHM AS7, Applied statistics, VOL.17, 1968, with  modifications
    !!--++    by A.J.MILLER
    !!--<<
    !!--++    Note: If RMAX = 1.E-04 then the diagonal elements of C will be
    !!--++          accurate to about 4 dec. digits.
    !!-->>
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Syminv(A, N, C, Nullty, Ifault, Rmax)
       !---- Arguments ----!
       real(kind=cp), dimension(:),  intent(in)     :: A
       integer,                      intent(in)     :: N
       real(kind=cp), dimension(:),  intent(out)    :: C
       integer,                      intent(out)    :: Nullty
       integer,                      intent(out)    :: Ifault
       real(kind=cp),                intent(out)    :: Rmax

       !---- Local variables ----!
       integer                          :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow
       real(kind=cp), parameter         :: zero = 0.0_cp, one = 1.0_cp
       real(kind=cp)                    :: x
       real(kind=cp), dimension(size(a)):: w  !workspace array that was in the argument

       w=zero
       nrow = n
       ifault = 1
       if (nrow > 0) then
          ifault = 0
          !---- Cholesky factorization of A, result in C
          call chola(a, nrow, c, nullty, ifault, rmax, w)
          if (ifault == 0) then
             !---- Invert C & form the product (CINV)'*CINV, where CINV is the inverse
             !---- of C, row by row starting with the last row.
             !---- IROW = The row number, NDIAG = Location of last element in the row.
             nn = nrow * (nrow+1) / 2
             irow = nrow
             ndiag = nn
             do
                if (c(ndiag) /= zero) then
                   l = ndiag
                   do i = irow, nrow
                      w(i) = c(l)
                      l = l + i
                   end do
                   icol = nrow
                   jcol = nn
                   mdiag = nn
                   do
                      l = jcol
                      x = zero
                      if (icol == irow) x = one / w(irow)
                      k = nrow
                      do
                         if (k /= irow) then
                            x = x - w(k) * c(l)
                            k = k - 1
                            l = l - 1
                            if (l > mdiag) l = l - k + 1
                            cycle
                         else
                            exit
                         end if
                      end do

                      c(l) = x / w(irow)
                      if (icol == irow) exit
                      mdiag = mdiag - icol
                      icol = icol - 1
                      jcol = jcol - 1
                   end do

                else
                   l = ndiag
                   do j = irow, nrow
                      c(l) = zero
                      l = l + j
                   end do
                end if ! (c(ndiag) /= zero)

                ndiag = ndiag - irow
                irow = irow - 1
                if (irow /= 0) cycle
                exit
             end do
          end if
       end if

       return
    End Subroutine Syminv

    !!--++
    !!--++ Subroutine Update(A, N, Z, Sig, W, Ir, Mk, Eps)
    !!--++    real,dimension(:), intent(out)      :: A       !
    !!--++    integer,           intent(in)       :: N       !
    !!--++    real,dimension(n), intent(in out)   :: Z       !
    !!--++    real,              intent(in)       :: Sig
    !!--++    real,dimension(n), intent(in out)   :: W
    !!--++    integer,           intent(   out)   :: Ir
    !!--++    integer,           intent(   out)   :: mk
    !!--++    real,              intent(in)       :: eps
    !!--++
    !!--++    Routine called by LOCAL Subroutine
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Update(a, n, z, sig, w, ir, mk, eps)
       !---- Arguments ----!
       real(kind=cp), dimension(:),intent(in out)   :: A       !
       integer,                    intent(in)       :: N       !
       real(kind=cp),dimension(n), intent(in out)   :: Z       !
       real(kind=cp),              intent(in)       :: Sig
       real(kind=cp),dimension(n), intent(in out)   :: W
       integer,                    intent(in out)   :: Ir
       integer,                    intent(in)       :: mk
       real(kind=cp),              intent(in)       :: eps

       !---- Local Variables ----!
       integer                   :: j, jj, ij, jp1, i, ii, mm
       real(kind=cp)             :: ti, v, tim, al, r, b, gm, y
       real(kind=cp), parameter  :: zero = 0.0_cp, one = 1.0_cp, four = 4.0_cp

       !---- UPDATE FACTORS GIVEN IN A  SIG*Z*Z-TRANSPOSE IS ADDED
       !---- FIRST EXECUTABLE STATEMENT
       if (n <= 1) then
          a(1) = a(1) + sig * z(1) * z(1)
          ir = 1
          if (a(1) > zero) return
          a(1) = zero
          ir = 0
       else
          do    !one-iteration loop
             if (sig <= zero) then
                if (sig == zero .or. ir == 0) return
                ti = one / sig
                jj = 0
                if (mk /= 0) then  !   l*w = z on input
                   do j = 1, n
                      jj = jj + j
                      if (a(jj) /= zero) ti = ti + (w(j)*w(j)) / a(jj)
                   end do
                else
                   !---- solve l*w = z
                   w(1:n) = z(1:n)
                   do j = 1, n
                      jj = jj + j
                      v = w(j)
                      if (a(jj) <= zero) then
                         w(j) = zero
                      else
                         ti = ti + (v*v) / a(jj)
                         if (j /= n) then
                            ij = jj
                            jp1 = j + 1
                            do i = jp1, n
                               ij = ij + i - 1
                               w(i) = w(i) - v * a(ij)
                            end do
                         end if
                      end if
                   end do
                end if

                !---- SET    TI, TIM    AND W
                if (ir > 0) then
                   if (ti > zero) then
                      ti = eps / sig
                      if (eps == zero) ir = ir - 1
                   else
                      if (mk-1 <= 0) then
                         mm = 0
                         tim = one / sig
                         exit   !one-iteration loop
                      end if
                   end if
                else
                   ti = zero
                   ir = -ir - 1
                end if

                tim = ti
                ii = jj
                i = n
                do j = 1, n
                   if (a(ii) /= zero) tim = ti - (w(i)*w(i)) / a(ii)
                   w(i) = ti
                   ti = tim
                   ii = ii - i
                   i = i - 1
                end do
                mm = 1
                jj=0
                exit   !one-iteration loop
             else
                mm = 0
                tim = one / sig
                jj = 0
             end if
             exit !one-iteration loop
          end do   !one-iteration loop

          !---- UPDATE A
          do j = 1, n
             jj = jj + j
             ij = jj
             jp1 = j + 1
             !---- update a(j,j)
             v = z(j)
             if (a(jj) <= zero) then
                !----  a(j,j) .eq. zero
                if (ir <= 0 .and. sig >= zero .and. v /= zero) then
                   ir = 1 - ir
                   a(jj) = (v*v) / tim
                   if (j == n) return
                   do i = jp1, n
                      ij = ij + i - 1
                      a(ij) = z(i) / v
                   end do
                   return
                end if
                ti = tim
             else
                !---- a(j,j) .gt. zero
                al = v / a(jj)
                ti = w(j)
                if (mm == 0) ti = tim + v * al
                r = ti / tim
                a(jj) = r * a(jj)
                if (r == zero) exit
                if (j == n) exit

                !---- update remainder of column j
                b = al / ti
                if (r <= four) then
                   do i = jp1, n
                      ij = ij + i - 1
                      z(i) = z(i) - v * a(ij)
                      a(ij) = a(ij) + b * z(i)
                   end do
                else
                   gm = tim / ti
                   do i = jp1, n
                      ij = ij + i - 1
                      y = a(ij)
                      a(ij) = b * z(i) + y * gm
                      z(i) = z(i) - v * y
                   end do
                end if
                tim = ti
             end if
          end do
          if (ir < 0) ir = -ir
       end if

       return
    End Subroutine Update

    !!--++
    !!--++ Subroutine Urdmn(x, n)
    !!--++    integer,           intent(in)    :: N
    !!--++    real,dimension(n), intent(out)   :: X
    !!--++
    !!--++    Generates n random numbers
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Urdmn(X, N)
       !---- Arguments ----!
       integer,                    intent(in)   :: n
       real(kind=cp),dimension(n), intent(out)  :: x

       call random_number(x)

       return
    End Subroutine Urdmn

 End Module Optimization_Procedures
