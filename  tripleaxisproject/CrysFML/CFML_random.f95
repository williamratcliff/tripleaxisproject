!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: RANDOM_GENER
!!----   INFO: A module for random number generation for differents
!!----         distributions.
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----
!!--..    Distribution                    Function/subroutine name
!!--..    --------------------------------------------------------
!!--..    Normal (Gaussian)               random_normal
!!--..    Gamma                           random_gamma
!!--..    Chi-squared                     random_chisq
!!--..    Exponential                     random_exponential
!!--..    Weibull                         random_Weibull
!!--..    Beta                            random_beta
!!--..    t                               random_t
!!--..    Multivariate normal             random_mvnorm
!!--..    Generalized inverse Gaussian    random_inv_gauss
!!--..    Poisson                         random_Poisson
!!--..    Binomial                        random_binomial1   *
!!--..                                    random_binomial2   *
!!--..    Negative binomial               random_neg_binomial
!!--..    von Mises                       random_von_Mises
!!--..    Cauchy                          random_Cauchy
!!--..
!!--..    Generate a random ordering of the integers 1 .. N
!!--..    random_order
!!--..
!!--..    Initialize (seed) the uniform random number generator
!!--..    for ANY compiler seed_random_number
!!--..
!!--..    Two functions are provided for the binomial distribution.
!!--..    If the parameter values remain constant, it is recommended that the
!!--..    first subroutine is used (random_binomial1).   If one or both of the
!!--..    parameters change, use the second subroutine (random_binomial2).
!!--..
!!--..    The compilers own random number generator,
!!--..    SUBROUTINE RANDOM_NUMBER(r), is used to provide a source of
!!--..    uniformly distributed random numbers.
!!--..
!!--..    At this stage, only one random number is generated at each call to
!!--..    one of the functions above.
!!--..
!!--..    The module uses the following functions which are included here:
!!--..    bin_prob to calculate a single binomial probability lngamma
!!--..    to calculate the logarithm to base e of the gamma subroutine
!!--..
!!--..    Some of the code is adapted from Dagpunar"s book:
!!--..        Dagpunar, J. "Principles of random variate generation"
!!--..        Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!!--..
!!--..    In most of Dagpunar"s routines, there is a test to see whether
!!--..    the value of one or two floating-point parameters has changed
!!--..    since the last call.
!!--..    These tests have been replaced by using a logical variable FIRST.
!!--..    This should be set to .TRUE. on the first call using new values
!!--..    of the parameters, and .FALSE. if the parameter values are the
!!--..    same as for the previous call.
!!--..
!!--..    Version 1.11, 4 January 1999
!!--..    Changes from version 1.01
!!--..    1. The random_order, random_Poisson & random_binomial routines
!!--..       have been replaced with more efficient routines.
!!--..    2. A routine, seed_random_number, has been added to seed the
!!--..       uniform random number generator. This requires input of the
!!--..       required number of seeds for the particular compiler from a
!!--..       specified I/O unit such as a keyboard.
!!--..    3. Made compatible with Lahey"s ELF90.
!!--..
!!--..       Author: Alan Miller
!!--..       CSIRO Division of Mathematical & Information Sciences
!!--..       Private Bag 10, Clayton South MDC
!!--..       Clayton 3169, Victoria, Australia
!!--..       Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
!!--..       e-mail: Alan.Miller @ vic.cmis.csiro.au
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    ERR_MESS_RANDOM
!!----    ERR_RANDOM
!!--++    DP                        [Private]
!!--++    SP                        [Private]
!!--++    HALF                      [Private]
!!--++    ONE                       [Private]
!!--++    TWO                       [Private]
!!--++    VLARGE                    [Private]
!!--++    VSMALL                    [Private]
!!--++    ZERO                      [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       BIN_PROB
!!--++       GG_F                   [Private]
!!--++       GPG_F                  [Private]
!!--++       GPP_F                  [Private]
!!----       INIT_ERR_RANDOM
!!----       INTEGRAL
!!----       LNGAMMA
!!----       RANDOM_BETA
!!----       RANDOM_BINOMIAL1
!!----       RANDOM_BINOMIAL2
!!----       RANDOM_CAUCHY
!!----       RANDOM_CHISQ
!!----       RANDOM_EXPONENTIAL
!!----       RANDOM_GAMMA
!!----       RANDOM_GAMMA1
!!----       RANDOM_GAMMA2
!!----       RANDOM_INV_GAUSS
!!----       RANDOM_MVNORM
!!----       RANDOM_NEG_BINOMIAL
!!----       RANDOM_NORMAL
!!----       RANDOM_ORDER
!!----       RANDOM_POISSON
!!----       RANDOM_T
!!----       RANDOM_VON_MISES
!!----       RANDOM_WEIBULL
!!----       SEED_RANDOM_NUMBER
!!----
!!
 Module Random_Gener

    !---- Use Modules ----!

    implicit none

    !---- List of public subroutines ----!
    public  :: bin_prob, lngamma, random_beta, random_binomial1, random_binomial2, random_cauchy,   &
               random_chisq, random_exponential, random_gamma, random_gamma1, random_gamma2,        &
               random_inv_gauss, random_neg_binomial, random_normal, random_poisson, random_t,      &
               random_von_mises, random_weibull, init_err_random, integral, random_mvnorm,          &
               random_order, seed_random_number

    !---- List of private subroutines ----!
    private :: gpp_f, gpg_f, gg_f

    !---- Definitions ----!

    !---- Local Variables ----!

    !!----
    !!---- ERR_MESS_RANDOM
    !!----    character(len=150), public :: Err_Mess_Random
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: Err_Mess_Random = " "

    !!----
    !!---- ERR_RANDOM
    !!----    logical, public :: err_random
    !!----
    !!----    Logical Variable indicating an error in RANDOM_GENER module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_random

    !!--++
    !!--++ DP
    !!--++    integer, parameter, private  :: dp
    !!--++
    !!--++    (PRIVATE)
    !!--++    DP: Double precision
    !!--++
    !!--++ Update: February - 2005
    !!
    integer, parameter, private  :: dp = selected_real_kind(12, 60)

    !!--++
    !!--++ SP
    !!--++    integer, parameter, private  :: sp
    !!--++
    !!--++    (PRIVATE)
    !!--++    SP: Simple Precision
    !!--++
    !!--++ Update: February - 2005
    !!
    integer, parameter, private  :: sp = selected_real_kind(6, 30)

    !!--++
    !!--++ HALF
    !!--++    real, parameters, private  :: half
    !!--++
    !!--++    (PRIVATE)
    !!--++    Half=0.5
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: half = 0.5

    !!--++
    !!--++ ONE
    !!--++    real, parameters, private  :: one
    !!--++
    !!--++    (PRIVATE)
    !!--++    One=1.0
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: one = 1.0

    !!--++
    !!--++ TWO
    !!--++    real, parameters, private  :: two
    !!--++
    !!--++    (PRIVATE)
    !!--++    Two=2.0
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: two = 2.0

    !!--++
    !!--++ VLARGE
    !!--++    real, parameters, private  :: vlarge
    !!--++
    !!--++    (PRIVATE)
    !!--++    VLarge=huge(1.0)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: vlarge = huge(1.0)

    !!--++
    !!--++ VSMALL
    !!--++    real, parameters, private  :: vsmall
    !!--++
    !!--++    (PRIVATE)
    !!--++    VSmall=tiny(1.0)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: vsmall = tiny(1.0)

    !!--++
    !!--++ ZERO
    !!--++    real, parameters, private  :: zero
    !!--++
    !!--++    (PRIVATE)
    !!--++    Zero=0.0
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private  :: zero = 0.0

 Contains

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Bin_Prob(N, P, R,Fn_Val)
    !!----    integer,          intent(in) :: n         !  In ->
    !!----    real(kind=sp),    intent(in) :: p         !  In ->
    !!----    integer,          intent(in) :: r         !  In ->
    !!----    real(kind=sp)    ,intent(out):: fn_val    ! Out ->
    !!----
    !!----    Calculate a binomial probability
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Bin_Prob(N, P, R, Fn_Val)
       !---- Arguments ----!
       integer, intent(in)         :: n, r
       real(kind=sp), intent(in)   :: p
       real(kind=sp),intent(out)   :: fn_val

       !---- Local variable ----!
       real(kind=dp)               :: n1,r1,nr1,rn1,rr1,rnr1

       n1=real(n+1)
       r1=real(r+1)
       nr1=real(n-r+1)
       call lngamma(n1,rn1)
       call lngamma(r1,rr1)
       call lngamma(nr1,rnr1)
       fn_val = exp( rn1 - rr1 - rnr1 &
                    + r*log(p) + (n-r)*log(one - p) )

       return
    End Subroutine Bin_Prob

    !!--++
    !!--++ Subroutine Gg_F(Methodeg,Gg)
    !!--++    integer,          intent(in) :: methodeg
    !!--++    real(kind=sp)    ,intent(out):: gg
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Gg_F(Methodeg,Gg)
       !---- Arguments ----!
       integer,       intent(in) :: methodeg
       real(kind=sp), intent(out):: gg

       !---- Local variables ----!
       real(kind=sp), parameter :: dpi  =6.283185308       ! 2*pi
       real(kind=sp), parameter :: sqdse=0.857763885       ! sqrt(2/e)
       real(kind=sp)            :: u,v,ran1,ran2
       logical                  :: sifin

       if (methodeg == 2)then
          call random_number(ran1)
          call random_number(ran2)
          gg=sqrt(-2*log(ran1))*cos(dpi*ran2)
       else if (methodeg == 3) then
          sifin=.false.
          do
             if (sifin) exit
             call random_number(ran1)
             call random_number(ran2)
             u=ran1
             v=(2*ran2-1)*sqdse
             sifin = ((0.25*v*v/(u*u)) <= -log(u))
          end do
          gg=v/u
       end if

       return
    End Subroutine Gg_F

    !!--++
    !!--++ Subroutine Gpg_F(Mt,Methodeg,Gpstab,Gpg)
    !!--++    real(kind=sp),    intent(in) :: mt
    !!--++    integer,          intent(in) :: methodeg
    !!--++    integer,          intent(in) :: gpstab
    !!--++    integer,          intent(out):: gpg
    !!--++
    !!--++    (PRIVATE)
    !!--++    Poisson distribution by approx.gauss. (N(0,1) -> P(MT))
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Gpg_F(Mt,Methodeg,Gpstab,Gpg)
       !---- Arguments ----!
       real(kind=sp), intent(in) :: mt
       integer,       intent(in) :: methodeg,gpstab
       integer      , intent(out):: gpg

       !---- Local variables ----!
       integer       :: x
       real(kind=sp) :: gg

       if (gpstab <= 0) then
          call gg_f(methodeg,gg)
          x=nint(mt+sqrt(mt)*gg)
       else if (gpstab == 1) then
          call gg_f(methodeg,gg)
          x=nint((0.5*gg + sqrt(mt))**2 - 0.33)
       end if
       if (x < 0) x=0
       gpg=x

       return
    End Subroutine Gpg_F

    !!--++
    !!--++ Subroutine Gpp_F(Mt,Gpp)
    !!--++    real(kind=sp),    intent(in) :: mt
    !!--++    integer,          intent(out):: gpp
    !!--++
    !!--++    (PRIVATE)
    !!--++    Poisson distribution by the product method
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Gpp_F(Mt,Gpp)
       !---- Arguments ----!
       real(kind=sp), intent(in) :: mt
       integer     , intent(out) :: gpp

       !---- Local variables ----!
       real(kind=sp) :: p,r,ran
       integer       :: k

       call random_number(ran)
       p=ran
       k=0
       r=exp(-mt)
       do
          if (p < r) exit
          call random_number(ran)
          p=p*ran
          k=k+1
       end do
       gpp=k

       return
    End Subroutine Gpp_F

    !!----
    !!---- Subroutine Init_Err_Random()
    !!----
    !!----    Initialize the errors flags in Random_Gener
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Random()

       err_random=.false.
       err_mess_random=" "

       return
    End Subroutine Init_Err_Random

    !!----
    !!---- Subroutine Integral(A, B, Resultt, Dk)
    !!----    real(kind=sp), intent(in)  :: a
    !!----    real(kind=sp), intent(in)  :: b
    !!----    real(kind=sp), intent(out) :: resultt
    !!----    real(kind=sp), intent(in)  :: dk
    !!----
    !!----    Gaussian integration of exp(k.cosx) from a to b.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Integral(A, B, Resultt, Dk)
       !---- Arguments ----!
       real(kind=dp), intent(in)      :: dk
       real(kind=sp), intent(in)      :: a, b
       real(kind=sp), intent(out)     :: resultt

       !---- Local variables ----!
       real(kind=dp)             :: xmid, rangee, x1, x2
       real(kind=dp),dimension(3):: x = (/0.238619186083197_dp, 0.661209386466265_dp, 0.932469514203152_dp/)
       real(kind=dp),dimension(3):: w = (/0.467913934572691_dp, 0.360761573048139_dp, 0.171324492379170_dp/)
       integer                   :: i

       xmid = (a + b)/2.0_dp
       rangee = (b - a)/2.0_dp
       resultt = 0.0_dp
       do i = 1, 3
          x1 = xmid + x(i)*rangee
          x2 = xmid - x(i)*rangee
          resultt = resultt + w(i)*(exp(dk*cos(x1)) + exp(dk*cos(x2)))
       end do

       resultt = resultt * rangee

       return
    End Subroutine Integral

    !!----
    !!---- Subroutine Lngamma(X,Fn_Val)
    !!----    real(kind=sp) (kind=dp), intent(in) :: x         !  In ->
    !!----    real(kind=sp) (kind=dp)             :: fn_val    ! Out ->
    !!----
    !!----    Logarithm to base e of the gamma subroutine.
    !!----     Accurate to about 1.e-14. (Alan Miller)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Lngamma(X,Fn_Val)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp),intent(out) :: fn_val

       !---- Local variables ----!
       real(kind=dp) :: a1 = -4.166666666554424e-02_dp, a2 = 2.430554511376954e-03_dp,  &
                        a3 = -7.685928044064347e-04_dp, a4 = 5.660478426014386e-04_dp,  &
                        temp, arg, productt, lnrt2pi = 9.189385332046727e-1_dp,      &
                        pi = 3.141592653589793_dp, eps=tiny(1.0_dp)
       logical       :: reflect

       !---- lngamma is not defined if x = 0 or a negative integer.
       if (.not.(x > 0.0_dp .or. abs(x-real(int(x),dp)) > eps) ) then
          fn_val = 0.0_dp
          return
       end if

       !---- If x < 0, use the reflection formula:
       !----  gamma(x) * gamma(1-x) = pi * cosec(pi.x)
       reflect = (x < 0.0_dp)
       if (reflect) then
          arg = 1.0_dp - x
       else
          arg = x
       end if

       !---- Increase the argument, if necessary, to make it > 10.
       productt = 1.0_dp
       do
          if (arg <= 10.0_dp) then
             productt = productt * arg
             arg = arg + 1.0_dp
             cycle
          end if
          exit
       end do

       !---- Use a polynomial approximation to Stirling"s formula.
       !---- N.B. The real(kind=sp) Stirling"s formula is used here, not
       !---- the simpler, but less accurate formula given by De Moivre
       !---- in a letter to Stirling, which is the one usually quoted.
       arg = arg - 0.5_dp
       temp = 1.0_dp/arg**2
       fn_val = lnrt2pi + arg * (log(arg) - 1.0_dp  + &
                (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - log(productt)
       if (reflect) then
          temp = sin(pi * x)
          fn_val = log(pi/temp) - fn_val
       end if

       return
    End Subroutine Lngamma

    !!----
    !!---- Subroutine Random_Beta(Aa, Bb, First,Fn_Val)
    !!----    real(kind=sp), intent(in)    :: aa          !  In -> shape parameter from distribution (0 < real)
    !!----    real(kind=sp), intent(in)    :: bb          !  In -> shape parameter from distribution (0 < real)
    !!----    logical, intent(in)          :: first       !  In ->
    !!----    real(kind=sp)                :: fn_val      ! Out ->
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate in [0,1] from a beta distribution with
    !!----    density proportional to beta**(aa-1) * (1-beta)**(bb-1) using cheng"s log
    !!----    logistic method.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Beta(Aa, Bb, First,Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)    :: aa, bb
       logical, intent(in)          :: first
       real(kind=sp),intent(out)    :: fn_val

       !---- Local variables ----!
       real(kind=sp), parameter  :: aln4 = 1.3862944_sp
       real(kind=sp)             :: a, b, g, r, s, x, y, z
       real(kind=sp), save       :: d, f, h, t, c
       logical,       save       :: swap

       call init_err_random()
       if (aa <= zero .or. bb <= zero) then
          err_random=.true.
          err_mess_random="Impermissible Shape Parameter Value(s)"
          return
       end if

       if (first) then                        ! initialization, if necessary
          a = aa
          b = bb
          swap = b > a
          if (swap) then
             g = b
             b = a
             a = g
          end if
          d = a/b
          f = a+b
          if (b > one) then
             h = sqrt((two*a*b - f)/(f - two))
             t = one
          else
             h = b
             t = one/(one + (a/(vlarge*b))**b)
          end if
          c = a+h
       end if

       do
          call random_number(r)
          call random_number(x)
          s = r*r*x
          if (r < vsmall .or. s <= zero) cycle
          if (r < t) then
             x = log(r/(one - r))/h
             y = d*exp(x)
             z = c*x + f*log((one + d)/(one + y)) - aln4
             if (s - one > z) then
                if (s - s*z > one) cycle
                if (log(s) > z) cycle
             end if
             fn_val = y/(one + y)
          else
             if (4.0*s > (one + one/d)**f) cycle
             fn_val = one
          end if
          exit
       end do

       if (swap) fn_val = one - fn_val

       return
    End Subroutine Random_Beta

    !!----
    !!---- Subroutine Random_Binomial1(N, P, First,Ival)
    !!----    integer,       intent(in) :: n          !  In -> Number Of Bernoulli Trials       (1 <= Integer)
    !!----    real(kind=sp), intent(in) :: p          !  In -> Bernoulli Success Probability    (0 <= Real <= 1)
    !!----    logical,       intent(in) :: first      !  In -> .TRUE. for the first call using the current parameter values
    !!----                                                     .FALSE. if the values of (n,p) are unchanged from last call
    !!----    integer,       intent(out):: ival       ! Out ->
    !!----
    !!----    Generates A Random Binomial Variate Using C.D.Kemp"s method.
    !!----    This algorithm is suitable when many random variates are
    !!----    required with the SAME parameter values for n & p.
    !!----
    !!--..    Reference: Kemp, C.D. (1986). `A modal method for generating
    !!--..    binomial variables", Commun. Statist. - Theor. Meth. 15(3),
    !!--..    805-813.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Binomial1(N, P, First, Ival)
       !---- Arguments ----!
       integer,       intent(in)  :: n
       real(kind=sp), intent(in)  :: p
       logical,       intent(in)  :: first
       integer,       intent(out) :: ival

       !---- Local variables ----!
       integer                  :: ru, rd
       integer, save            :: r0
       real(kind=sp)            :: u, pd, pu
       real(kind=sp), save      :: odds_ratio, p_r

       if (first) then
          r0 = (n+1)*p
          call bin_prob(n, p, r0, p_r)
          odds_ratio = p / (one - p)
       end if

       call random_number(u)
       u = u - p_r
       if (u < zero) then
          ival = r0
          return
       end if

       pu = p_r
       ru = r0
       pd = p_r
       rd = r0
       do
          rd = rd - 1
          if (rd >= 0) then
             pd = pd * (rd+1) / (odds_ratio * (n-rd))
             u = u - pd
             if (u < zero) then
                ival = rd
                return
             end if
          end if

          ru = ru + 1
          if (ru <= n) then
             pu = pu * (n-ru+1) * odds_ratio / ru
             u = u - pu
             if (u < zero) then
                ival = ru
                return
             end if
          end if
       end do

       !---- This point should not be reached, but just in case:
       ival = r0

       return
    End Subroutine Random_Binomial1

    !!----
    !!---- Subroutine Random_Binomial2(N, Pp, First, Ival)
    !!----    integer,       intent(in) :: n      !  In -> The number of trials in the binomial distribution
    !!----                                           from which a random deviate is to be generated.
    !!----    real(kind=sp), intent(in) :: pp     !  In -> The probability of an event in each trial of the
    !!----                                           binomial distribution from which a random deviate
    !!----                                           is to be generated.
    !!----    logical,       intent(in) :: first  !  In -> .TRUE. for the first call to perform initialization
    !!----                                           the set FIRST = .FALSE. for further calls using the
    !!----                                           same pair of parameter values (N, P)
    !!----    integer,       intent(out):: ival
    !!----
    !!----    Generates a single random deviate from a binomial
    !!----    distribution whose number of trials is N and whose
    !!----    probability of an event in each trial is P.
    !!----    Random_binomial2 <-- A random deviate yielding the number
    !!----    of events from N independent trials, each of which has a
    !!----    probability of event P.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Binomial2(N, Pp, First, Ival)
       !---- Arguments ----!
       real(kind=sp), intent(in)    :: pp
       integer, intent(in)          :: n
       logical, intent(in)          :: first
       integer, intent(out)         :: ival

       !---- local variables ----!
       real(kind=sp)            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
       integer                  :: i, ix, ix1, k, mp
       integer, save            :: m
       real(kind=sp), save      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
                                   xlr, p2, p3, p4, qn, r, g

       !---- setup, perform only when parameters change
       if (first) then
          p = min(pp, one-pp)
          q = one - p
          xnp = n * p
       end if

       if (xnp > 30.0) then
          if (first) then
             ffm = xnp + p
             m = ffm
             fm = m
             xnpq = xnp * q
             p1 = int(2.195*sqrt(xnpq) - 4.6*q) + half
             xm = fm + half
             xl = xm - p1
             xr = xm + p1
             c = 0.134 + 20.5 / (15.3 + fm)
             al = (ffm-xl) / (ffm - xl*p)
             xll = al * (one + half*al)
             al = (xr - ffm) / (xr*q)
             xlr = al * (one + half*al)
             p2 = p1 * (one + c + c)
             p3 = p2 + c / xll
             p4 = p3 + c / xlr
          end if

          !---- generate variate, binomial mean at least 30.
          do

             call random_number(u)
             u = u * p4
             call random_number(v)

             !---- triangular region
             if (u <= p1) then
                ix = xm - p1 * v + u
                if (pp > half) ix = n - ix
                ival = ix
                return
             end if

             !---- parallelogram region
             if (u <= p2) then
                x = xl + (u-p1) / c
                v = v * c + one - abs(xm-x) / p1
                if (v > one .or. v <= zero) cycle
                ix = x
             else
                !---- left tail
                if (u <= p3) then
                   ix = xl + log(v) / xll
                   if (ix < 0) cycle
                   v = v * (u-p2) * xll
                else
                   !---- right tail
                   ix = xr - log(v) / xlr
                   if (ix > n) cycle
                   v = v * (u-p3) * xlr
                end if
             end if

             !---- determine appropriate way to perform accept/reject test
             k = abs(ix-m)
             if (k <= 20 .or. k >= xnpq/2-1) then
                !---- explicit evaluation
                f = one
                r = p / q
                g = (n+1) * r
                if (m < ix) then
                   mp = m + 1
                   do i = mp, ix
                      f = f * (g/i-r)
                   end do
                else if (m > ix) then
                   ix1 = ix + 1
                   do i = ix1, m
                      f = f / (g/i-r)
                   end do
                end if

                if (v > f) then
                   cycle
                else
                   if (pp > half) ix = n - ix
                   ival = ix
                   return
                end if
             end if

             !---- squeezing using upper and lower bounds on log(f(x))
             amaxp = (k/xnpq) * ((k*(k/3.0 + 0.625) + 0.1666666666666)/xnpq + half)
             ynorm = -k * k / (2.0*xnpq)
             alv = log(v)
             if (alv<ynorm - amaxp) then
                if (pp > half) ix = n - ix
                ival = ix
                return
             end if
             if (alv>ynorm + amaxp) cycle

             !---- stirling"s (actually de moivre"s) formula to machine accuracy
             !---- for the final acceptance/rejection test
             x1 = ix + 1
             f1 = fm + one
             z = n + 1 - fm
             w = n - ix + one
             z2 = z * z
             x2 = x1 * x1
             f2 = f1 * f1
             w2 = w * w
             if (alv - (xm*log(f1/x1) + (n-m+half)*log(z/w) + (ix-m)*log(w*p/(x1*q)) +    &
                (13860.0-(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0 +               &
                (13860.0-(462.0-(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0 +                &
                (13860.0-(462.0-(132.0-(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0 +               &
                (13860.0-(462.0-(132.0-(99.0-140.0/w2)/w2)/w2)/w2)/w/166320.0) > zero) then
                cycle
             else
                if (pp > half) ix = n - ix
                ival = ix
                return
             end if
             exit
          end do

       else
          !---- inverse cdf logic for mean less than 30
          if (first) then
             qn = q ** n
             r = p / q
             g = r * (n+1)
          end if

          do
             ix = 0
             f = qn
             call random_number(u)
             do
                if (u >= f) then
                   if (ix > 110) exit
                   u = u - f
                   ix = ix + 1
                   f = f * (g/ix - r)
                   cycle
                end if
                exit
             end do
             if (ix > 110) cycle
             exit
          end do
       end if

       if (pp > half) ix = n - ix
       ival = ix

       return
    End Subroutine Random_Binomial2

    !!----
    !!---- Subroutine Random_Cauchy(Fn_Val)
    !!----    real(kind=sp) ,intent(out):: fn_val
    !!----
    !!----    Generate a random deviate from the standard Cauchy distribution
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Cauchy(Fn_Val)
       !---- Arguments ----!
       real(kind=sp),intent(out):: fn_val

       !---- Local variables ----!
       real(kind=sp),dimension(2)     :: v

       do
          call random_number(v)
          v = two*(v - half)
          if (abs(v(2)) < vsmall) cycle               ! test for zero
          if (v(1)**2 + v(2)**2 < one) exit
       end do

       fn_val = v(1) / v(2)

       return
    End Subroutine Random_Cauchy

    !!----
    !!---- Subroutine Random_Chisq(Ndf, First,Fn_Val)
    !!----    integer,       intent(in) :: ndf      !  In ->
    !!----    logical,       intent(in) :: first    !  In ->
    !!----    real(kind=sp) ,intent(out):: fn_val   ! Out ->
    !!----
    !!----    Generates a random variate from the chi-squared
    !!----    distribution with ndf degrees of freedom
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Chisq(Ndf, First, Fn_Val)
       !---- Arguments ----!
       integer, intent(in) :: ndf
       logical, intent(in) :: first
       real(kind=sp) ,intent(out)  :: fn_val

       call random_gamma(half*ndf, first, fn_val)
       fn_val = two * fn_val

       return
    End Subroutine Random_Chisq

    !!----
    !!---- Subroutine Random_Exponential(Fn_Val)
    !!----    real(kind=sp) ,intent(out)   :: fn_val
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate in [0,infinity) from a negative exponential
    !!----    distribution wlth density proportional to exp(-random_exponential), using inversion.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Exponential(Fn_Val)
       !---- Arguments ----!
       real(kind=sp),intent(out)  :: fn_val

       !---- Local variable ----!
       real(kind=sp)  :: r

       do
          call random_number(r)
          if (r > zero) exit
       end do

       fn_val = -log(r)

       return
    End Subroutine Random_Exponential

    !!----
    !!---- Subroutine Random_Gamma(S, First,Fn_Val)
    !!----    real(kind=sp), intent(in)  :: s       !  In -> shape parameter of distribution (0.0 < real)
    !!----    logical,       intent(in)  :: first   !  In ->
    !!----    real(kind=sp) ,intent(out) :: fn_val  ! Out ->
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random gamma variate.
    !!--..
    !!--..    calls either random_gamma1 (S > 1.0)
    !!--..    or random_exponential (S = 1.0)
    !!--..    or random_gamma2 (S < 1.0).
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Gamma(S, First,Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)  :: s
       logical, intent(in)        :: first
       real(kind=sp) ,intent(out) :: fn_val

       call init_err_random()
       if (s <= zero) then
          err_random=.true.
          err_mess_random="Shape Parameter Value Must Be Positive"
          return
       end if

       if (s > one) then
          call random_gamma1(s, first,fn_val)
       else if (s < one) then
          call random_gamma2(s, first,fn_val)
       else
          call random_exponential(fn_val)
       end if

       return
    End Subroutine Random_Gamma

    !!----
    !!---- Subroutine Random_Gamma1(S, First,Fn_Val)
    !!----    real(kind=sp), intent(in) :: s       !  In -> shape parameter of distribution (1.0 < real)
    !!----    logical,       intent(in) :: first
    !!----    real(kind=sp) ,intent(out):: fn_val
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate in [0,infinity) from a gamma distribution
    !!----    with density proportional to gamma**(s-1)*exp(-gamma), based upon best"s t
    !!----    distribution method
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Gamma1(S, First,Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)   :: s
       logical,       intent(in)   :: first
       real(kind=sp) ,intent(out)  :: fn_val

       !---- Local variables ----!
       real(kind=sp)            :: d, r, g, f, x
       real(kind=sp), save      :: b, h
       real(kind=sp), parameter :: sixty4 = 64.0, three = 3.0, pt75 = 0.75

       call init_err_random()
       if (s <= one) then
          err_random=.true.
          err_mess_random="Impermissible Shape Parameter Value"
          return
       end if

       if (first) then                        ! initialization, if necessary
          b = s - one
          h = sqrt(three*s - pt75)
       end if

       do
          call random_number(r)
          g = r - r*r
          if (g <= zero) cycle
          f = (r - half)*h/sqrt(g)
          x = b + f
          if (x <= zero) cycle
          call random_number(r)
          d = sixty4*g*(r*g)**2
          if (d <= zero) exit
          if (d*x < x - two*f*f) exit
          if (log(d) < two*(b*log(x/b) - f)) exit
       end do
       fn_val = x

       return
    End Subroutine Random_Gamma1

    !!----
    !!---- Subroutine Random_Gamma2(S, First,Fn_Val)
    !!----    real(kind=sp), intent(in)  :: s        !  In -> shape parameter of distribution (1.0 < real)
    !!----    logical,       intent(in)  :: first
    !!----    real(kind=sp) ,intent(out) :: fn_val
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate in [0,infinity) from
    !!----    a gamma distribution with density proportional to
    !!----    gamma2**(s-1) * exp(-gamma2), using a switching method.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Gamma2(S, First, Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)  :: s
       logical,       intent(in)  :: first
       real(kind=sp) ,intent(out) :: fn_val

       !---- Local variables ----!
       real(kind=sp)       :: r, x, w
       real(kind=sp), save :: a, p, c, uf, vr, d

       call init_err_random()
       if (s <= zero .or. s >= one) then
          err_random=.true.
          err_mess_random="Shape Parameter Value Outside Permitted Range"
          return
       end if

       if (first) then                        ! initialization, if necessary
          a = one - s
          p = a/(a + s*exp(-a))
          if (s < vsmall) then
             err_random=.true.
             err_mess_random="Shape Parameter Value Too Small"
             return
          end if
          c = one/s
          uf = p*(vsmall/a)**s
          vr = one - vsmall
          d = a*log(a)
       end if

       do
          call random_number(r)
          if (r >= vr) then
             cycle
          else if (r > p) then
             x = a - log((one - r)/(one - p))
             w = a*log(x)-d
          else if (r > uf) then
             x = a*(r/p)**c
             w = x
          else
             fn_val = zero
             return
          end if

          call random_number(r)
          if (one-r <= w .and. r > zero) then
             if (r*(w + one) >= one) cycle
             if (-log(r) <= w) cycle
          end if
          exit
       end do
       fn_val = x

       return
    End Subroutine Random_Gamma2

    !!----
    !!---- Subroutine Random_Inv_Gauss(H, B, First, Fn_Val)
    !!----    real(kind=sp), intent(in)    :: h       !  In -> parameter of distribution (0 <= real)
    !!----    real(kind=sp), intent(in)    :: b       !  In -> parameter of distribution (0 < real)
    !!----    logical,       intent(in)    :: first
    !!----    real(kind=sp), intent(out)   :: fn_val
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate in [0,infinity] from
    !!----    a reparameterised generalised inverse gaussian (gig) distribution
    !!----    with density proportional to  gig**(h-1) * exp(-0.5*b*(gig+1/gig))
    !!----    using a ratio method.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Inv_Gauss(H, B, First, Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)  :: h, b
       logical,       intent(in)  :: first
       real(kind=sp), intent(out) :: fn_val

       !---- Local variables ----!
       real(kind=sp)            :: ym, xm, r, w, r1, r2, x
       real(kind=sp), save      :: a, c, d, e
       real(kind=sp), parameter :: quart = 0.25

        call init_err_random()
       if (h < zero .or. b <= zero) then
           err_random=.true.
           err_mess_random="Impermissible Distribution Parameter Values"
          return
       end if

       if (first) then                        ! initialization, if necessary
          if (h > quart*b*sqrt(vlarge)) then
             err_random=.true.               !Not possible in F (pure Functions!)
             err_mess_random="The Ratio H:B Is Too Small"
             return
          end if
          e = b*b
          d = h + one
          ym = (-d + sqrt(d*d + e))/b
          if (ym < vsmall) then
             err_random=.true.
             err_mess_random="The Value Of B Is Too Small"
             return
          end if

          d = h - one
          xm = (d + sqrt(d*d + e))/b
          d = half*d
          e = -quart*b
          r = xm + one/xm
          w = xm*ym
          a = w**(-half*h) * sqrt(xm/ym) * exp(-e*(r - ym - one/ym))
          if (a < vsmall) then
             err_random=.true.
             err_mess_random="The Value Of H Is Too Large"
             return
          end if
          c = -d*log(xm) - e*r
       end if

       do
          call random_number(r1)
          if (r1 <= zero) cycle
          call random_number(r2)
          x = a*r2/r1
          if (x <= zero) cycle
          if (log(r1) < d*log(x) + e*(x + one/x) + c) exit
       end do

       fn_val = x

       return
    End Subroutine Random_Inv_Gauss

    !!----
    !!---- Subroutine Random_Mvnorm(N, H, D, F, First, X, Ier)
    !!----    integer,       intent(in)  :: n            !  In -> number of variates in vector (input,integer >= 1)
    !!----    real(kind=sp), intent(in)  :: h(n)         !  In -> vector of means
    !!----    real(kind=sp), intent(in)  :: d(n*(n+1)/2) !  In -> variance matrix (j> = i)
    !!----    real(kind=sp), intent(out) :: f(n*(n+1)/2) ! Out -> lower triangular decomposition of variance matrix (j <= i)
    !!----    real(kind=sp), intent(out) :: x(n)         ! Out -> delivered vector
    !!----    logical,       intent(in)  :: first        !  In -> .true. if this is the first call of the routine
    !!----                                                         or if the distribution has changed since the last
    !!----                                                         call of the routine. otherwise set to .false.
    !!----    integer,       intent(out) :: ier          ! Out ->  = 1 if the input covariance matrix is not +ve definite
    !!----                                                         = 0 otherwise
    !!----
    !!----    Generates an n variate random normal vector using
    !!----    a cholesky decomposition.
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!--..
    !!---- Update: February - 2005
    !!
    Subroutine Random_Mvnorm(N, H, D, F, First, X, Ier)
       !---- Arguments ----!
       integer, intent(in)  :: n
       real(kind=sp), dimension(n),        intent(in) :: h
       real(kind=sp), dimension(n*(n+1)/2),intent(in) :: d
       real(kind=sp), dimension(n),        intent(out):: x
       real(kind=sp), dimension(n*(n+1)/2),intent(out):: f
       logical,                            intent(in) :: first
       integer,                            intent(out):: ier

       !---- Local variables ----!
       integer       :: j, i, m
       real(kind=sp) :: y, v
       integer, save :: n2

        call init_err_random()
       if (n < 1) then
          err_random=.true.
          err_mess_random="Size Of Vector Is Non Positive"
          return
       end if
       ier = 0
       if (first) then                        ! initialization, if necessary
          n2 = 2*n
          if (d(1) < zero) then
             ier = 1
             return
          end if
          f(1) = sqrt(d(1))
          y = one/f(1)
          do j = 2,n
             f(j) = d(1+j*(j-1)/2) * y
          end do

          do i = 2,n
             v = d(i*(i-1)/2+i)
             do m = 1,i-1
                v = v - f((m-1)*(n2-m)/2+i)**2
             end do
             if (v < zero) then
                ier = 1
                return
             end if
             v = sqrt(v)
             y = one/v
             f((i-1)*(n2-i)/2+i) = v
             do j = i+1,n
                v = d(j*(j-1)/2+i)
                do m = 1,i-1
                   v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
                end do ! m = 1,i-1
                f((i-1)*(n2-i)/2 + j) = v*y
             end do ! j = i+1,n
          end do ! i = 2,n
       end if

       x(1:n) = h(1:n)
       do j = 1,n
           call random_normal(y)
          do i = j,n
             x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
          end do ! i = j,n
       end do ! j = 1,n

       return
    End Subroutine Random_Mvnorm

    !!----
    !!---- Subroutine Random_Neg_Binomial(Sk, P, Ival)
    !!----    real(kind=sp), intent(in)    :: sk     !  In -> Number of failures required (dagpunar's words!)
    !!----                                                    the "power" parameter of the negative binomial  (0 < real)
    !!----    real(kind=sp), intent(in)    :: p      !  In -> bernoulli success probability  (0 < real(kind=sp) < 1)
    !!----    integer,       intent(out)   :: ival
    !!----
    !!----    Generates a random negative binomial variate using unstored
    !!----    inversion and/or the reproductive property.
    !!----
    !!--..    the parameter h is set so that unstored inversion only is
    !!--..    used when p <= h, otherwise a combination of unstored
    !!--..    inversion and the reproductive property is used.
    !!----
    !!--..    adapted from fortran 77 code from the book:
    !!--..    dagpunar, j. "principles of random variate generation"
    !!--..    clarendon press, oxford, 1988.   isbn 0-19-852202-9
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Neg_Binomial(Sk, P, Ival)
       !---- Arguments ----!
       real(kind=sp), intent(in)   :: sk, p
       integer , intent(out)       :: ival

       !---- Local variables ----!
       !---- the parameter uln = -log(machine"s smallest real(kind=sp) number).
       real(kind=sp), parameter    :: h = 0.7
       real(kind=sp)               :: q, x, st, uln, v, r, s, y, g
       integer            :: k, i, n

       call init_err_random()
       if (sk <= zero .or. p <= zero .or. p >= one) then
          err_random=.true.
          err_mess_random="impermissible distribution parameter values"
          return
       end if
       q = one - p
       x = zero
       st = sk
       if (p > h) then
          v = one/log(p)
          k = st
          do i = 1,k
             do
                call random_number(r)
                if (r > zero) exit
             end do
             n = v*log(r)
             x = x + n
          end do
          st = st - k
       end if
       s = zero
       uln = -LOG(vsmall)
       if (st > -uln/log(q)) then
          err_random=.true.
          err_mess_random=" P Is Too Large For This Value Of Sk"
          return
       end if
       y = q**st
       g = st
       call random_number(r)
       do
          if (y > r) exit
          r = r - y
          s = s + one
          y = y*p*g/s
          g = g + one
       end do
       ival = x + s + half

       return
    End Subroutine Random_Neg_Binomial

    !!----
    !!---- Subroutine Random_Normal(Fn_Val)
    !!----    real(kind=sp)  :: fn_val
    !!----
    !!--..    Adapted from the following Fortran 77 code ALGORITHM 712,
    !!--..    Collected Algorithms From Acm.
    !!--..    This Work Published In Transactions On Mathematical Software,
    !!--..    Vol. 18, No. 4, December, 1992, Pp. 434-435.
    !!----
    !!----    The subroutine random_normal() returns a normally distributed
    !!----    pseudo-random number with zero mean and unit variance.
    !!----    The algorithm uses the ratio of uniforms method of A.J. Kinderman
    !!----    and J.F. Monahan augmented with quadratic bounding curves.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Normal(Fn_Val)
       !---- Arguments ----!
       Real(kind=sp),intent(out) :: fn_val

       !---- Local variables ----!
       Real :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
              r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

       !---- Generate P = (u,v) uniform in rectangle enclosing ----!
       !---- acceptance region                                 ----!
       do
          call random_number(u)
          call random_number(v)
          v = 1.7156 * (v - half)

          !---- Evaluate the quadratic form ----!
          x = u - s
          y = abs(v) - t
          q = x**2 + y*(a*y - b*x)

          !---- Accept P if inside inner ellipse ----!
          if (q < r1) exit

          !---- Reject P if outside outer ellipse ----!
          if (q > r2) cycle

          !---- Reject P if outside acceptance region ----!
          if (v**2 < -4.0*log(u)*u**2) exit
       end do

       !---- Return ratio of P"s coordinates as the normal deviate ----!
       fn_val = v/u

       return
    End Subroutine Random_Normal

    !!----
    !!---- Subroutine Random_Order(Order, N)
    !!----    integer, intent(in)  :: n
    !!----    integer, intent(out) :: order(n)
    !!----
    !!----    Generate a random ordering of the integers 1 ... n.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Order(Order, N)
       !---- Arguments ----!
       integer,              intent(in)  :: n
       integer, dimension(n),intent(out) :: order

       !---- Local variables ----!
       integer :: i, j, k
       real(kind=sp)    :: wk

       do i = 1, n
          order(i) = i
       end do

       !---- starting at the end, swap the current last indicator with one
       !---- randomly chosen from those preceeding it.
       do i = n, 2, -1
          call random_number(wk)
          j = 1 + i * wk
          if (j < i) then
             k = order(i)
             order(i) = order(j)
             order(j) = k
          end if
       end do

       return
    End Subroutine Random_Order

    !!----
    !!---- Subroutine Random_Poisson(Mu,Genpoi)
    !!----    real(kind=sp), intent(in)    :: mu
    !!----    integer , intent(out)        :: genpoi
    !!----
    !!----    Generates a single random deviate from a Poisson distribution
    !!----    with mean mu.
    !!--..    Based on: J.Berruyer, A.Antoniadis, A.Filhol  (1985)  Rapport ILL 85AN19T
    !!----
    !!----    Method for generation of a random number according to Poisson distribution
    !!----        METHODE  = 1/2/3 if (centred limit )/(box-muller)/(reject)
    !!----        METHODEG = 1/2/3 if (FR-1)/(product)/(gauss)
    !!----        GPSTAB   = 0/1/2 stabilisation of the variance
    !!----                   (non)/(0.33)/(0.375) cf. EFRON)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Poisson(mt, Genpoi)
       !---- Arguments ----!
       real(kind=sp), intent(in)    :: mt
       integer ,      intent(out)   :: Genpoi

       !---- Local Variables ----!
       integer ::  methode,methodeg,gpstab

       if (mt < 0.0) then
          genpoi=mt
          return
       end if

       if (mt < 30.0) then
          methode=2
       else if ((mt >= 30.0) .and. (mt < 87.0)) then
          methode=3
          methodeg=2
          gpstab=1
       else if ((mt >= 87.0) .and. (mt < 500.0)) then
          methode=3
          methodeg=3
          gpstab=1
       else if (mt >= 500.0) then
          methode=3
          methodeg=2
          gpstab=0
       end if

       if (methode == 2) then
          call gpp_f(mt, genpoi)
       else if (methode == 3) then
          call gpg_f(mt,methodeg,gpstab, genpoi)
       else
          genpoi=mt
       end if

       return
    End Subroutine Random_Poisson

    !!----
    !!---- Subroutine Random_T(M, Fn_Val)
    !!----    integer,      intent(in) :: m       !  In -> degrees of freedom of distribution (1 <= 1nteger)
    !!----    real(kind=sp),intent(out):: fn_val
    !!----
    !!--..    Adapted from Fortran 77 code from the book:
    !!--..    Dagpunar, J. "Principles of random variate generation"
    !!--..    Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !!----
    !!----    Subroutine generates a random variate from a
    !!----    t distribution using kinderman and monahan"s ratio method.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_T(M, Fn_Val)
       !---- Arguments ----!
       integer,      intent(in)  :: m
       real(kind=sp),intent(out) :: fn_val

       !---- Local variables ----!
       real(kind=sp), save      :: s, c, a, f, g
       real(kind=sp)            :: r, x, v
       real(kind=sp), parameter :: three = 3.0, four = 4.0, quart = 0.25,   &
                                   five = 5.0, sixteen = 16.0
       integer                  :: mm = 0

       call init_err_random()
       if (m < 1) then
          err_random=.true.
          err_mess_random="Impermissible Degrees Of Freedom"
          return
       end if

       if (m /= mm) then                    ! initialization, if necessary
          s = m
          c = -quart*(s + one)
          a = four/(one + one/s)**c
          f = sixteen/a
          if (m > 1) then
             g = s - one
             g = ((s + one)/g)**c*sqrt((s+s)/g)
          else
             g = one
          end if
          mm = m
       end if

       do
          call random_number(r)
          if (r <= zero) cycle
          call random_number(v)
          x = (two*v - one)*g/r
          v = x*x
          if (v > five - a*r) then
             if (m >= 1 .and. r*(v + three) > f) cycle
             if (r > (one + v/s)**c) cycle
          end if
          exit
       end do

       fn_val = x

       return
    End Subroutine Random_T

    !!----
    !!---- Subroutine Random_Von_Mises(K, First, Fn_Val)
    !!----    real(kind=sp), intent(in)  :: k        !  In -> Parameter of the von Mises distribution
    !!----    logical,       intent(in)  :: first    !  In -> set to .TRUE. the first time that the subroutine
    !!----                                                    is called
    !!----    real(kind=sp) ,intent(out) :: fn_val
    !!----
    !!--..    When first = .TRUE., the subroutine sets up starting
    !!--..    values and may be very much slower.
    !!----
    !!----    Von Mises Distribution
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Von_Mises(K, First, Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)  :: k
       logical,       intent(in)  :: first
       real(kind=sp) ,intent(out) :: fn_val

       !---- Local variables ----!
       integer                             :: j, n
       integer, save                       :: nk
       real(kind=sp), parameter            :: pi = 3.14159265
       real(kind=sp), save,dimension(20)   :: p
       real(kind=sp), save,dimension(0:20) :: theta
       real(kind=sp)                       :: sump, r, th, lambda, rlast
       real(kind=dp)                       :: dk

       call init_err_random()
       if (first) then                        ! initialization, if necessary
          if (k < zero) then
             err_random=.true.
             err_mess_random="Error in argument k for random_von_Mises"
             return
          end if
          nk = k + k + one
          if (nk > 20) then
             err_random=.true.
             err_mess_random="Error in argument k for random_von_Mises"
             return
          end if
          dk = k
          theta(0) = zero
          if (k > half) then
             !---- set up array p of probabilities.
             sump = zero
             do j = 1, nk
                if (j < nk) then
                   theta(j) = acos(one - j/k)
                else
                   theta(nk) = pi
                end if

                !---- numerical integration of e^[k.cos(x)] from theta(j-1)
                !---- to theta(j)
                call integral(theta(j-1), theta(j), p(j), dk)
                sump = sump + p(j)
             end do ! j = 1, nk
             p(1:nk) = p(1:nk) / sump
          else
             p(1) = one
             theta(1) = pi
          end if                         ! if k > 0.5
       end if                           ! if first
       call random_number(r)
       do j = 1, nk
          r = r - p(j)
          if (r < zero) exit
       end do
       r = -r/p(j)

       do
          th = theta(j-1) + r*(theta(j) - theta(j-1))
          lambda = k - j + one - k*cos(th)
          n = 1
          rlast = lambda

          do
             call random_number(r)
             if (r > rlast) exit
             n = n + 1
             rlast = r
          end do
          if (n /= 2*(n/2)) exit         ! is n even?
          call random_number(r)
       end do

       fn_val = sign(th, (r - rlast)/(one - rlast) - half)

       return
    End Subroutine Random_Von_Mises

    !!----
    !!---- Subroutine Random_Weibull(A, Fn_Val)
    !!----    real(kind=sp), intent(in)  :: a
    !!----    real(kind=sp), intent(out) :: fn_val
    !!----                                a
    !!----                         a-1  -x
    !!----               f(x) = a.x    e
    !!----
    !!----    Generates a random variate from the Weibull distribution with
    !!----    probability density as shown before.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Random_Weibull(a, Fn_Val)
       !---- Arguments ----!
       real(kind=sp), intent(in)  :: a
       real(kind=sp) ,intent(out) :: fn_val

       !---- For speed, there is no checking that a is
       !---- not zero or very small.
       call random_exponential(fn_val)
       fn_val = fn_val** (one/a)

       return
    End Subroutine Random_Weibull

    !!----
    !!---- Subroutine Seed_Random_Number(I_input,I_output)
    !!----    integer, optional, intent(in)  :: I_input
    !!----    integer, optional, intent(in)  :: I_output
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Seed_Random_Number(I_input,I_output)
       !---- Arguments ----!
       integer, optional, intent(in)  :: I_input
       integer, optional, intent(in)  :: I_output

       !---- Local variables ----!
       integer                            :: k,lun1,lun2
       integer, dimension(:), allocatable :: seed

       call random_seed(size=k)
       allocate( seed(k) )

       lun1=5
       lun2=6
       if (present(i_input))  lun1=i_input
       if (present(i_output)) lun2=i_output

       write(unit=*, fmt= "(a, i2, a)")" Enter ", k, " integers for random no. seeds: "
       read(unit=lun1, fmt=*) seed

       write(unit=lun2,fmt="(a, (7i10))") " Random no. seeds: ", seed
       call random_seed(put=seed)

       return
    End Subroutine Seed_Random_Number

 End Module Random_Gener
