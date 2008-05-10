!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: OPTIMIZATION_LSQ
!!----   INFO: Module implementing several algorithms for non-linear least-squares.
!!----         At present only the Levenberg-Marquardt method is implemented.
!!----
!!---- HISTORY
!!----    Update: January - 2004
!!----
!!----    October - 1997 Created by JRC
!!---- DEPENDENCIES
!!--++    Use Math_Gen,     only : Sp, Dp, Invert_Matrix
!!----
!!---- VARIABLES
!!----    ERR_LSQ
!!----    ERR_MESS_LSQ
!!----    NPAR
!!----    LSQ_CONDITIONS_TYPE
!!----    LSQ_STATE_VECTOR_TYPE
!!--++    CH                      [Private]
!!--++    CORREL                  [Private]
!!--++    CURV_MAT                [Private]
!!--++    NAMFREE                 [Private]
!!--++    PN                      [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!--++       FCHISQ              [Private]
!!----
!!----    Subroutines:
!!--++       BOX_CONSTRAINTS     [Private]
!!--++       CURFIT              [Private]
!!----       MARQUARDT_FIT
!!--++       OUTPUT_CYC          [Private]
!!--++       OUTPUT_FINAL        [Private]
!!----
!!
 Module Optimization_LSQ
    !---- Use Files ----!
    Use Math_Gen, only: Sp, Dp, Invert_Matrix

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public  :: marquardt_fit

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: fchisq

    !---- List of private subroutines ----!
    private :: curfit,box_constraints,output_cyc,output_final

    !---- Variable Definitions ----!

    !!----
    !!---- ERR_LSQ
    !!----    logical, public  :: ERR_LSQ
    !!----
    !!----    Logical variable. The vaule .true. indicates that an error condition occurs
    !!----
    !!---- Update: February - 2005
    !!
    logical, public   :: Err_lsq =.false.

    !!----
    !!---- ERR_MESS_LSQ
    !!----    Character(len=132), public  :: ERR_MESS_LSQ
    !!----
    !!----    Character variable containing the error message associated to the
    !!----    ocurrence of an error condition
    !!----
    !!---- Update: February - 2005
    !!
    Character(len=132), public  :: ERR_MESS_LSQ

    !!----
    !!---- NPAR
    !!----    integer, parameter, public  :: npar
    !!----
    !!----    Maximum number of free parameters (500)
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public   :: npar=500   !Maximum number of free parameters

    !!----
    !!----  Type, public :: LSQ_Conditions_type
    !!----     logical          :: constr=.false.  ! if true box constraint of percent% are applied to parameters
    !!----     logical          :: reached =.false.! if true convergence was reached in the algorithm
    !!----     integer          :: corrmax=50      ! value of correlation in % to output
    !!----     integer          :: icyc            ! number of cycles of refinement
    !!----     integer          :: npvar           !number of effective free parameters of the model
    !!----     integer          :: iw              ! indicator for weighting scheme (if iw=1 => w=1/yc)
    !!----     real(kind=sp)    :: percent         ! %value of maximum variation of a parameter w.r.t.
    !!----                                         ! the intial value before fixing it
    !!----  End Type LSQ_Conditions_type
    !!----
    !!----  Derived type encapsulating all necessary conditions for running the LSQ algorithm
    !!----
    !!---- Update: March - 2005
    !!
    Type, public :: LSQ_Conditions_type
       logical          :: constr          ! if true box constraint of percent% are applied to parameters
       logical          :: reached         ! if true convergence was reached in the algorithm
       integer          :: corrmax         ! value of correlation in % to output
       integer          :: icyc            ! number of cycles of refinement
       integer          :: npvar           ! number of effective free parameters of the model
       integer          :: iw              ! indicator for weighting scheme (if iw=1 => w=1/yc)
       real(kind=sp)    :: percent         ! %value of maximum variation of a parameter w.r.t.
                                           ! the intial value before fixing it
    End Type LSQ_Conditions_type

    !!----
    !!----  Type, public :: LSQ_State_Vector_type
    !!----     integer                            :: np         !total number of model parameters <= npar
    !!----     real(kind=sp),     dimension(npar) :: pv         !Vector of parameters
    !!----     real(kind=sp),     dimension(npar) :: spv        !Vector of standard deviations
    !!----     integer,           dimension(npar) :: code       !pointer for selecting variable parameters
    !!----     character(len=15), dimension(npar) :: nampar=" " !Names of parameters
    !!----  End Type LSQ_State_Vector_type
    !!----
    !!----  Derived type encapsulating the vector state defining a set of parameter
    !!----  for calculating the model function and running the LSQ algorithm.
    !!----
    !!---- Update: March - 2005
    !!

    Type, public :: LSQ_State_Vector_type
       integer                            :: np         !total number of model parameters <= npar
       real(kind=sp),     dimension(npar) :: pv         !Vector of parameters
       real(kind=sp),     dimension(npar) :: spv        !Vector of standard deviations
       integer,           dimension(npar) :: code       !pointer for selecting variable parameters
       character(len=15), dimension(npar) :: nampar     !Names of parameters
    End Type LSQ_State_Vector_type

    !!--++
    !!--++ CH
    !!--++    real(kind=sp),     dimension(npar), private   :: ch
    !!--++
    !!--++    (PRIVATE)
    !!--++    Vector holding the change in the values of parameters (ch = pn - pv)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(npar), private   :: ch         ! ch = pn - pv

    !!--++
    !!--++ CORREL
    !!--++    real(kind=sp), dimension(npar,npar), public  :: correl
    !!--++
    !!--++    Variance/covariance/correlation matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(npar,npar), private  :: correl     !Variance/covariance/correlation matrix


    !!--++
    !!--++ CURV_MAT
    !!--++    real(kind=sp), dimension(npar,npar), public  :: curv_mat
    !!--++
    !!--++    Curvature matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(npar,npar), private  :: curv_mat   !Curvature matrix

    !!--++
    !!--++ NAMFREE
    !!--++    character(len=15), dimension(npar), private   :: namfree
    !!--++
    !!--++  (PRIVATE)
    !!--++  Names of refined parameters
    !!--++
    !!--++ Update: February - 2005
    !!
    character(len=15), dimension(npar), private   :: namfree    !Names of refined parameters

    !!--++
    !!--++ PN
    !!--++    real(kind=sp), dimension(npar), public   :: pn
    !!--++
    !!--++    Vector with new values of parameters
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(npar), private   :: pn   !Vector with new values of parameters


 Contains

    !---- Functions ----!

    !!--++
    !!--++ Function Fchisq(Nfr,Nobs,Y,W,Yc) Result(Chisq)
    !!--++    integer, intent(in)              :: Nfr
    !!--++    integer, intent(in)              :: Nobs
    !!--++    real,    intent(in),dimension(:) :: y
    !!--++    real,    intent(in),dimension(:) :: w
    !!--++    real,    intent(in),dimension(:) :: yc
    !!--++
    !!--++    (PRIVATE)
    !!--++    Evaluate reduced chi2
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fchisq(Nfr,Nobs,Y,W,Yc) Result(Chisq)
       !---- Arguments ----!
       integer, intent(in)              :: nfr,nobs
       real,    intent(in),dimension(:) :: y
       real,    intent(in),dimension(:) :: w
       real,    intent(in),dimension(:) :: yc
       real(kind=sp)                    :: chisq

       !---- Local variables ----!
       integer :: i

       chisq=0.0
       do i=1,nobs
          chisq=chisq+w(i)*(y(i)-yc(i))**2
       end do
       chisq=chisq/real(nfr)

       return
    End Function Fchisq

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Box_Constraints(A,Sa,Fixed,c,vs)
    !!--++    real, dimension (npar), intent(in out) :: a
    !!--++    real, dimension (npar), intent(in out) :: sa
    !!--++    logical               , intent(   out) :: fixed
    !!--++    Type(LSQ_Conditions_type),  intent(in) :: c     !conditions for refinement
    !!--++    Type(LSQ_State_Vector_type),intent(in) :: vs    !State vector
    !!--++
    !!--++    (PRIVATE)
    !!--++    This subroutine avoid a peak-position parameter undergoing a change
    !!--++    greater than percent% w.r.t. the initial value. This is probably enough
    !!--++    to avoid lack of convergence when a peak makes an excursion outside
    !!--++    the current angular range. The parameter is fixed for the next cycle
    !!--++    The intensity parameters are also constrained to be strictly positive
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Box_Constraints(A,Sa,Fixed,c,vs)
       !---- Arguments ----!
       real, dimension (npar), intent(in out) :: a
       real, dimension (npar), intent(in out) :: sa
       logical               , intent(   out) :: fixed
       Type(LSQ_Conditions_type),  intent(in) :: c     !conditions for refinement
       Type(LSQ_State_Vector_type),intent(in) :: vs    !State vector

       !---- Local variables ----!
       integer :: i,ncount
       real    :: per
       logical :: ifi

       fixed=.false.
       ncount=0
       per=c%percent*0.01
       do i=1,vs%np
          ifi=.false.
          if (vs%code(i) /= 0) then
             ncount=ncount+1
             if (abs(vs%pv(i)) > 0.0000001) then
                if (abs(a(ncount)-vs%pv(i))/vs%pv(i) > per )  ifi=.true.
                if (ifi) then
                   fixed=.true.
                   a(ncount) = vs%pv(i)
                   namfree(ncount) = trim(namfree(ncount))//"-*"
                   sa(ncount) = 0.0
                   pn(i) = vs%pv(i)
                   ch(i) = 0.0
                end if
             end if
          end if
       end do

       return
    End Subroutine Box_Constraints

    !!--++
    !!--++  Subroutine Curfit(Model_Functn, X, Y, W, Nobs, c,vs, A, Sa, Fl, Yc, Chir, Ifail)
    !!--++     real,    dimension(:),      intent(in)      :: x     !vector with abcisae
    !!--++     real,    dimension(:),      intent(in)      :: y     !Observed values
    !!--++     real,    dimension(:),      intent(in out)  :: w     !weight of observations
    !!--++     integer,                    intent(in)      :: nobs  !number of observations
    !!--++     Type(LSQ_Conditions_type),  intent(in)      :: c     !conditions for refinement
    !!--++     Type(LSQ_State_Vector_type),intent(in out)  :: vs    !State vector
    !!--++     real,dimension(:),          intent(in out)  :: a     !vector of parameter
    !!--++     real,dimension(:),          intent(in out)  :: sa    !estimated standard deviations
    !!--++     real,                       intent(in out)  :: fl    !Marquardt LAMBDA value
    !!--++     real,dimension(:),          intent(out)     :: yc    !Calculated
    !!--++     real,                       intent(out)     :: chir
    !!--++     integer,                    intent(out)     :: ifail
    !!--++
    !!--++     Interface
    !!--++      Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
    !!--++           integer,                    intent(in) :: iv
    !!--++           real,                       intent(in) :: xv
    !!--++           real,dimension(:),          intent(in) :: aa
    !!--++           real,                       intent(out):: ycalc
    !!--++           real,dimension(:),optional, intent(out):: der
    !!--++      End Subroutine Model_Functn
    !!--++     End Interface
    !!--++
    !!--++ Update: February - 2003
    !!
    Subroutine Curfit(Model_Functn, X, Y, W, Nobs, c, vs, A, Sa, Fl, Yc, Chir, Ifail)
       !---- Arguments ----!
       real,    dimension(:),      intent(in)      :: x     !vector with abcisae
       real,    dimension(:),      intent(in)      :: y     !Observed values
       real,    dimension(:),      intent(in out)  :: w     !weight of observations
       integer,                    intent(in)      :: nobs  !number of observations
       Type(LSQ_Conditions_type),  intent(in)      :: c     !conditions for refinement
       Type(LSQ_State_Vector_type),intent(in out)  :: vs    !State vector
       real,dimension(:),          intent(in out)  :: a     !vector of parameter
       real,dimension(:),          intent(in out)  :: sa    !estimated standard deviations
       real,                       intent(in out)  :: fl    !Marquardt LAMBDA value
       real,dimension(:),          intent(out)     :: yc    !Calculated
       real,                       intent(out)     :: chir
       integer,                    intent(out)     :: ifail

       Interface
        Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
             integer,                    intent(in) :: iv
             real,                       intent(in) :: xv
             real,dimension(:),          intent(in) :: aa
             real,                       intent(out):: ycalc
             real,dimension(:),optional, intent(out):: der
        End Subroutine Model_Functn
       End Interface

       !---- Local variables ----!
       real                     :: chisq1
       integer                  :: ntrials
       integer                  :: nfr,j,i,k,ntr
       logical                  :: change_par
       real,  dimension(c%npvar):: beta, b, der
       logical                  :: singular


       ntrials=100

       !---- Optimization with MARQUARDT ALGORITHM ----!
       nfr=nobs-c%npvar
       ifail=0
       change_par=.false.
       if (nobs > 2000) ntrials= 30

       !---- Evaluate ALFA and BETA matrices       ----!
       !---- (BETA=DP*ALFA)--- DP=BETA*ALFA(-1)    ----!
       beta(:) =0.0
       curv_mat(:,:)=0.0


       do i=1,nobs
          call model_functn(i,x(i),yc(i),a,der)
          if (c%iw == 1 .and. yc(i) /= 0.0) w(i)=1.0/yc(i)
          do j=1,c%npvar
             beta(j)=beta(j)+w(i)*(y(i)-yc(i))*der(j)
             do k=1,j
                curv_mat(j,k)=curv_mat(j,k)+w(i)*der(j)*der(k)
             end do
          end do
       end do

       do i=1,c%npvar
          if (curv_mat(i,i) < 1.0e-30) then
             Err_lsq =.true.
             write(unit=err_mess_lsq,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",i," -> ",trim(namfree(i))
             ifail=2
             return
          end if
       end do

       do j=1,c%npvar
          do k=1,j
             curv_mat(k,j)=curv_mat(j,k)
          end do
       end do

       !---- Evaluate CHI2 at starting point ----!
       chisq1=fchisq(nfr,nobs,y,w,yc)  !chi2 before inverting and applying shifts

       !---- Invert modified curvature matrix to find new parameters
       do ntr=1,ntrials
          do j=1,c%npvar
             do k=1,c%npvar
                correl(j,k)=curv_mat(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
             correl(j,j)=1.0+fl
          end do

          call Invert_Matrix(correl(1:c%npvar,1:c%npvar),correl(1:c%npvar,1:c%npvar),singular)

          if (singular) then
             Err_lsq =.true.
             do i=1,c%npvar
                if (correl(i,i) < 1.0e-20) then
                   j=i !iperm(i)
                   exit
                end if
             end do
             write(unit=err_mess_lsq,fmt="(a,i5,a,a)")  &
                  " => Singular matrix!!, problem with parameter no.:",j," -> ",trim(namfree(j))
             ifail=2
             return
          end if

          if (fl < 1.0e-20)  then
             chir=chisq1
             exit
          end if

          do j=1,c%npvar
             b(j)=a(j)
             do k=1,c%npvar
                b(j)=b(j)+beta(k)*correl(j,k)/sqrt(curv_mat(j,j)*curv_mat(k,k))
             end do
          end do


          do i=1,nobs
             call model_functn(i,x(i),yc(i),b)
          end do

          chir=fchisq(nfr,nobs,y,w,yc)

          if (ntrials > 1) then
             if (ntr == ntrials .or. fl > 1.0e+20 ) then    !no improvement
                ifail=1
                change_par=.false.
                exit
             end if
          end if

          !---- If Chi2 increase, increase fl and try again ----!
          if (chisq1-chir < 0.0) then  !Chi2 increases
             fl=10.0*fl                !Increase fl
             cycle
          end if
          change_par=.true.
          exit
       end do  ! ntrials

       !---- Evaluate parameters and uncertainties ----!
       do j=1,c%npvar
          sa(j)=sqrt(abs(correl(j,j)/curv_mat(j,j)))
       end do
       if (change_par) then
          do j=1,c%npvar
             a(j)=b(j)
          end do
       end if
       fl=fl/10.0

       return
    End Subroutine Curfit

    !!----
    !!---- Subroutine Marquardt_Fit(Model_Functn, X, Y, W, Yc, Nobs, c, vs, Ipr, Chi2, scroll_lines)
    !!----    real(kind=sp), dimension(:),intent(in)     :: x      !Vector of x-values
    !!----    real(kind=sp), dimension(:),intent(in)     :: y      !Vector of observed y-values
    !!----    real(kind=sp), dimension(:),intent(in out) :: w      !Vector of weights-values (1/variance)
    !!----    real(kind=sp), dimension(:),intent(   out) :: yc     !Vector of calculated y-values
    !!----    integer                    ,intent(in)     :: nobs   !Number of effective components of x,y,w,yc
    !!----    type(LSQ_conditions_type),  intent(in out) :: c      !Conditions for the algorithm
    !!----    type(LSQ_State_Vector_type),intent(in out) :: vs     !State vector for the model calculation
    !!----    integer                    ,intent(in)     :: Ipr    !Logical unit for writing
    !!----    real(kind=sp),              intent(out)    :: chi2   !Reduced Chi-2
    !!----    character(len=*),dimension(:), intent(out), optional  :: scroll_lines  !If present, part of the output is stored
    !!----                                                                           !in this text for treatment in the calling program
    !!----
    !!----    Model_functn                                            !Name of the subroutine calculating yc(i) for point x(i)
    !!----    Interface                                               !Interface for the Model_Functn subroutine
    !!----       Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
    !!----            integer,                    intent(in) :: iv     !Number of the component: "i" in x(i)
    !!----            real,                       intent(in) :: xv     !Value of x(i)
    !!----            real,                       intent(out):: ycalc  !Value of yc at point x(i) => ycalc=yc(i)
    !!----            real,dimension(:),          intent(in) :: aa     !Vector of parameters
    !!----            real,dimension(:),optional, intent(out):: der    !Derivatives of the function w.r.t. free parameters
    !!----       End Subroutine Model_Functn
    !!----    End Interface
    !!----
    !!----
    !!----    Subroutine for applying the Levenberg-Marquardt method for Least-Squares.
    !!----    The user must provide a model function according to the interface above.
    !!----    The model function should use at least some of the public variables of the
    !!----    present module in order to set the derivatives with respect to the model
    !!----    parameters. Examples of using this module are given in the program
    !!----    templates CW_fit and TOF_fit.
    !!----
    !!--..    INFORMATION
    !!--..        Author: Juan Rodriguez-Carvajal (based in text descriptions of the literature)
    !!--..                Translated, extracted and modified from the Fortran 77 code: XRFIT (JRC 1986)
    !!--..        Module implementing the Levenberg-Marquardt method for non-linear least-squares.
    !!--..        For using this module, the user must provide:
    !!--..            1: Number of cycles (c%icyc), type of weighting scheme (c%iw), number of
    !!--..               model parameters (vs%np), constraint conditions (c%constr, c%percent)
    !!--..            2: A character(len=15) name for all possible parameters of the
    !!--..               model stored in the array vs%nampar(:).
    !!--..            3: A set of initial value for all parameters stored in vs%pv(:)
    !!--..            4: A set of flags values to refine or fix the parameters vs%code(:)
    !!--..            5: The model function (subroutine model_functn) must be provided by the user
    !!--..               according to the interface described below. The actual name of the model function
    !!--..               is arbitrary and it is passed to the only public procedure "marquardt_fit" as a
    !!--..               dummy argument.
    !!--..
    !!--..       The values of all possible refinable parameters are stored in the array vs%pv(:).
    !!--..       The derivatives must be calculated within model_functn, by using the array vs%der(:)
    !!--..       The actual refined parameters a(:) are selected from the larger vs%pv(:) array by
    !!--..       the integer array of flags: vs%code(:). A value vs%code(j)=1 means that the j-th
    !!--..       parameter is to be varied. A value vs%code(k)=0 means that the k-th parameter is
    !!--..       kept fixed through the refinement cycles.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Marquardt_Fit(Model_Functn,X,Y,W,Yc,Nobs,c,vs,Ipr,Chi2,scroll_lines)
       !---- Arguments ----!
       real(kind=sp),   dimension(:), intent(in)             :: x,y
       real(kind=sp),   dimension(:), intent(in out)         :: w
       real(kind=sp),   dimension(:), intent(   out)         :: yc
       integer,                       intent(in)             :: nobs,Ipr
       type(LSQ_conditions_type),     intent(in out)         :: c
       type(LSQ_State_Vector_type),   intent(in out)         :: vs
       real(kind=sp),                 intent(   out)         :: chi2
       character(len=*),dimension(:), intent(out), optional  :: scroll_lines

       Interface
        Subroutine Model_Functn(iv,Xv,ycalc,aa,der)
             integer,                    intent(in) :: iv
             real,                       intent(in) :: xv
             real,dimension(:),          intent(in) :: aa
             real,                       intent(out):: ycalc
             real,dimension(:),optional, intent(out):: der
        End Subroutine Model_Functn
       End Interface

       !---- local variables ----!
       real(kind=sp), dimension(npar) :: a,sa
       integer                        :: ifail
       integer                        :: i,ncount,li,ntex
       real(kind=sp)                  :: fl,chi1
       logical                        :: fixed, write_cyc
       character(len=132)             :: line

       fixed = .false.
       write_cyc=.true.
       if (c%icyc > 50) write_cyc=.false.

       !---- Beginning the iterations ----!
       fl=0.001
       chi1=1.0e+30
       c%reached =.false.
       write(unit=ipr,fmt="(/,a)")   " -------------------------------------------------"
       write(unit=ipr,fmt="(a,i3,a)")" => Marquardt Least Squares Fitting for ",c%icyc," cycles"
       write(unit=ipr,fmt="(a,/)")   " -------------------------------------------------"

       write(unit=ipr,fmt="(a,i6)")        " => Number of observations   : ",nobs
       ncount=0
       do i=1,vs%np
          if (vs%code(i) == 0) cycle
          ncount=ncount+1
       end do
       write(unit=ipr,fmt="(a,i6)")        " => Number of free parameters: ",ncount
       write(unit=ipr,fmt="(2(a,f12.4),a)")" => Range of variable X      : (",x(1),",",x(nobs),")"

       ntex=0
       if(present(scroll_lines)) scroll_lines=" "  !Clear the stored text
       do li=1,c%icyc                      !Loop over icyc cycles
          ncount=0
          do i=1,vs%np
             if (vs%code(i) == 0) cycle
             ncount=ncount+1
             a(ncount)=vs%pv(i)
             namfree(ncount)=vs%nampar(i)
          end do

          c%npvar=ncount
          call curfit(Model_Functn,x,y,w,nobs,c,vs,a,sa,fl,yc,chi2,ifail)

          if (ifail /= 0) then
             line= " "
             write(unit=line,fmt="(a,i3)") " => IFAIL /= 0 : ",ifail
             if(present(scroll_lines)) then
             	ntex=ntex+1
                scroll_lines(ntex)=line
             end if
             if (ABS(chi2-chi1)*100/Chi2 <= 0.001 .and. ifail==1) then
                c%reached=.true.
                line= " "
                write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                if(present(scroll_lines)) then
                    ntex=ntex+1
                    scroll_lines(ntex)=line
                end if
                write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.00001 * Chi2(old)"
                exit
             end if

             write(unit=ipr,fmt="(/,a)")" => Marquardt convergence not reached (F-lambda increased too much)."
             write(unit=ipr,fmt="(a)" ) "    If Chi2 is reasonably low, you may have reached the minimum"
             write(unit=ipr,fmt="(a)" ) "    within the precision of the machine.                       "
             write(unit=ipr,fmt="(a)" ) "    Otherwise, look at the output file to search for the reason "
             write(unit=ipr,fmt="(a,/)")"    of anomalous behaviour and try to start from new parameters"

             write(unit=ipr,fmt="(/,a,/)")" =>  Diagonal elements of LSQ-matrix before inversion"
             do i=1,c%npvar
                write(unit=ipr,fmt="(a,i3,a,i3,a,f16.6)") "    "//namfree(i)//"  (",i,",",i,")", curv_mat(i,i)
             end do

             exit
          end if  !ifail==1

          if (abs(chi2-chi1)*100.0/chi2 >= 0.0001 ) then
             c%reached=.false.
             chi1=chi2
             ncount=0
             do i=1,vs%np
                if (vs%code(i) == 0) then
                   vs%spv(i)=0.0
                   pn(i)=vs%pv(i)
                   ch(i)=0.0
                else
                   ncount=ncount+1
                   pn(i)=a(ncount)
                   ch(i)=pn(i)-vs%pv(i)
                   vs%spv(i)=sa(ncount)
                   namfree(ncount)=vs%nampar(i)
                end if
             end do

             if (c%constr)  then
                call box_constraints(a,sa,fixed,c,vs)
                if (fixed) then
                   write(unit=ipr,fmt="(/,a,/,a,i2,a)") " => Some parameters have been restored to their initial value ", &
                                                        "    because of local divergence (change > ",nint(c%percent),"%)"
                   write(unit=ipr,fmt="(a,/)")          "    Look at the above list for parameters having a '*'"
                end if  !fixed
             end if

             line= " "
             write(unit=line,fmt="(a,i5,2(a,f14.6))")" => Cycle No.:",li,"  Chi2 =",chi2, "  Marquardt-Lambda: ",fl
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             end if

             if (write_cyc) call output_cyc(li,Ipr,chi2,vs)
             vs%pv(1:vs%np)=pn(1:vs%np)
             cycle
          else
             line= " "
             write(unit=line,fmt="(a)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             if(present(scroll_lines)) then
                 ntex=ntex+1
                 scroll_lines(ntex)=line
             end if
             write(unit=Ipr,fmt="(/,a,/)") " => Convergence reached! : Chi2(old)-Chi2(new) < 0.000001 * Chi2(old)"
             c%reached=.true.
             exit
          end if
       end do  !li=1,c%icyc

       if (c%reached) then
          if (fl <= 10.0) then
             fl=0.0
          else
             fl=0.5e-20
          end if
          call curfit(Model_Functn,x,y,w,nobs,c,vs,a,sa,fl,yc,chi2,ifail)
          ncount=0
          do i=1,vs%np
             if (vs%code(i) == 0) then
                vs%spv(i)=0.0
                pn(i)=vs%pv(i)
                ch(i)=0.0
             else
                ncount=ncount+1
                namfree(ncount)=vs%nampar(i)
                pn(i)=a(ncount)
                ch(i)=pn(i)-vs%pv(i)
                vs%spv(i)=sa(ncount)
                vs%pv(i)=a(ncount)
             end if
          end do
       end if

       call output_final(chi2,FL,nobs,x,y,yc,w,Ipr,c,vs)

       return
    End Subroutine Marquardt_Fit



    !!--++
    !!--++  Subroutine Output_Cyc(Ic,Lun,Chi2,vs)
    !!--++   integer,                    intent(in) :: ic  !cycle number
    !!--++   integer,                    intent(in) :: lun !logical number of the output file
    !!--++   real,                       intent(in) :: chi2
    !!--++   type(LSQ_State_Vector_type),intent(in) :: vs
    !!--++
    !!--++  (PRIVATE)
    !!--++  Subroutine for output information on each cycle
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Output_Cyc(Ic,Lun,Chi2,vs)
       !---- Arguments ----!
       integer,                    intent(in) :: ic  !cycle number
       integer,                    intent(in) :: lun !logical number of the output file
       real,                       intent(in) :: chi2
       type(LSQ_State_Vector_type),intent(in) :: vs
       !---- Local variables ----!
       integer                          :: i,j
       real                             :: rat

       !---- Writing during cycles ----!
       write(unit=lun,fmt="(/,/,a,i5,a,f14.6)")" => Cycle No.:",ic,"  Chi2 =",chi2
       write(unit=lun,fmt="(/,/,a,/)") &
            "    Name-Par        No.      Old-Value          Change        New-Value         Sigma        Change/Sigma"
       j=0
       do i=1,vs%np
          if (vs%code(i)/=0) then
             j=j+1
             if (vs%spv(i) > 1.0e-36) then
                rat=ch(i)/vs%spv(i)
             else
                rat=0.0
             end if
             write(unit=lun,fmt="(a,i6,a,5f16.6)") " "//namfree(j),i," ",vs%pv(i),ch(i),pn(i),vs%spv(i),rat
          end if
       end do

       return
    End Subroutine Output_Cyc

    !!--++
    !!--++  Subroutine Output_Final(Chi2,Nobs,X,Y,Yc,W,Lun,c,vs)
    !!--++   real,                       intent(in)     :: chi2
    !!--++   real,                       intent(in)     :: FL
    !!--++   integer,                    intent(in)     :: nobs
    !!--++   real,dimension(:),          intent(in)     :: x
    !!--++   real,dimension(:),          intent(in)     :: y
    !!--++   real,dimension(:),          intent(in)     :: yc
    !!--++   real,dimension(:),          intent(in)     :: w
    !!--++   integer,                    intent(in)     :: lun
    !!--++   type(LSQ_conditions_type),  intent(in)     :: c
    !!--++   type(LSQ_State_Vector_type),intent(in)     :: vs
    !!--++
    !!--++  (PRIVATE)
    !!--++  Subroutine for output information at the end of refinement
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Output_Final(Chi2,FL,Nobs,X,Y,Yc,W,Lun,c,vs)
       !---- Arguments ----!
       real,                       intent(in)     :: chi2
       real,                       intent(in)     :: FL
       integer,                    intent(in)     :: nobs
       real,dimension(:),          intent(in)     :: x
       real,dimension(:),          intent(in)     :: y
       real,dimension(:),          intent(in)     :: yc
       real,dimension(:),          intent(in)     :: w
       integer,                    intent(in)     :: lun
       type(LSQ_conditions_type),  intent(in)     :: c
       type(LSQ_State_Vector_type),intent(in)     :: vs

       !---- Local variables ----!
       real(kind=dp)                    :: rfact,rwfact,riobs,rex
       integer                          :: i,j,inum
       real                             :: del,g2

       !---- Final calculation and writings R-Factors calculations ----!
       rfact=0.0
       rwfact=0.0
       riobs=0.0
       do i=1,nobs
          riobs=riobs+y(i)
          del=y(i)-yc(i)
          rfact=rfact+abs(del)
          rwfact=rwfact+w(i)*del*del
       end do
       rfact=rfact/riobs*100.0
       rwfact=sqrt(rwfact/riobs)*100.0
       rex=sqrt(real(nobs-c%npvar)/riobs)*100.0
       write(unit=lun,fmt="(/,(3(a,f8.3)))") "  Rfact= ",rfact,"   Rwfact= ",rwfact,"   Rex= ",rex
       write(unit=lun,fmt="(/,a,F16.3)") "  Final value of Marquardt F-Lambda = ",FL

       !---- Correlation matrix ----!
       write(unit=lun,fmt="(/,a,/)")   " => Correlation Matrix: "

       do i=1,c%npvar
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/sqrt(curv_mat(i,i)*curv_mat(j,j))
             correl(j,i)=correl(i,j)
          end do
       end do
       do i=1,c%npvar
         g2=sqrt(correl(i,i))
          do j=i,c%npvar
             correl(i,j)=correl(i,j)/g2/sqrt(correl(j,j))*100.0
             correl(j,i)=correl(i,j)
          end do
       end do

       inum=0
       do i=1,c%npvar-1
          do j=i+1,c%npvar
             if (correl(i,j) > real(c%corrmax) ) then
                write(unit=lun,fmt="(a,i4,a,i2,a,a15,a,a15)") "    Correlation:",nint(correl(i,j)),  &
                     " > ",c%corrmax,"% for parameters:   ", adjustr(namfree(i))," & ", namfree(j)
                inum=inum+1
             end if
          end do
       end do
       if (inum == 0) then
          write(unit=lun,fmt="(/,a,i2,a)") " => There is no correlations greater than ",c%corrmax,"% "
       else
          write(unit=lun,fmt="(/,a,i3,a,i2,a)") " => There are ",inum," values of Correlation > ",c%corrmax,"%"
       end if

       write(unit=lun,fmt="(/,/,a,/,a,/)") "    FINAL LIST OF REFINED PARAMETERS AND STANDARD DEVIATIONS",&
                                           "    --------------------------------------------------------"
       write(unit=lun,fmt="(/,a,/)")       "       Name-Par      No.          Final-Value    Standard Deviation "
       do i=1,vs%np
          if (vs%code(i)/=0) write(unit=lun,fmt="(a,i6,a,2f20.6)") "   "//vs%nampar(i),i," ",vs%pv(i),vs%spv(i)
       end do
       write(unit=lun,fmt="(/,a,f10.5)") " => Final value of Chi2: ",chi2

     ! !---- Output of a file with the observed and calculated curves ----!
     ! open(unit=22,file="marquart_fit.xy",status="replace", action="write")
     ! write(unit=22,fmt="(a)") "!  X   Y-obs   Y-calc  Sigma"
     ! do i=1,nobs
     !    write(unit=22,fmt="(4f16.4)") x(i),y(i),yc(i),sqrt(1.0/w(i))
     ! end do
     ! close(unit=22)

       return
    End Subroutine Output_Final

 End Module Optimization_LSQ
