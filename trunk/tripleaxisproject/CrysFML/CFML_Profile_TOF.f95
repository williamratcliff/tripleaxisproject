!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: TOF_PROFILES
!!----   INFO: This module contains variables and procedures used by programs aiming
!!----         to handle T.O.F. powder diffraction patterns. The public procedures
!!----         have only one (erfc, erfcp) two (tof_Jorgensen, tof_Jorgensen_VonDreele)
!!----         or three arguments (tof_Ikeda_Carpenter).
!!----         The rest of needed variables are passed to or from the module as global
!!----         variables.
!!----
!!----         Example:
!!----             yi= IntegInt * tof_Jorgensen(dt,ider)
!!----         Calculates the intensity in channel i, with dt=TOFi-TOF(Bragg), ider=1 if
!!----         derivatives are to be calculated.
!!----
!!----         In FullProf, Omega is the general name given to the peak shape function
!!----
!!---- HISTORY
!!----    Update:  September  - 2004
!!----    Author : Juan Rodriguez-Carvajal (CEA/DSM/DRECAM/LLB)
!!----             Laurent C. Chapon (ISIS/RAL) => Ikeda-Carperter function
!!----
!!---- VARIABLES
!!--++    DP
!!--++    SP
!!--++    INV_8LN2
!!--++    TWO_OVER_PI
!!----    LORCOMP
!!----    DERIV_TOF
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       ERFC
!!--++       DERFCC            [Overloaded]
!!--++       ERFCC             [Overloaded]
!!--++       EXPI_E1           [Private]
!!----
!!----    Subroutines:
!!----       TOF_CARPENTER
!!----       TOF_JORGENSEN
!!----       TOF_JORGENSEN_VONDREELE
!!----
!!
 Module Tof_Profiles
    !---- Use Modules ----!

    !---- Variables ----!
    implicit none

    !---- List of public functions ----!
    public :: Tof_Jorgensen, Tof_Jorgensen_Vondreele, Tof_Carpenter, Erfc, Erfcp

    !---- List of private functions ----!
    private:: Expi_E1, Erfcc, Derfccp, Erfccp, Derfcc


    !---- Definitions ----!

    !!--++
    !!--++ DP
    !!--++    integer, parameter,private   :: dp = selected_real_kind(14, 300)
    !!--++
    !!--++    Private variable
    !!--++
    !!--++ Update: October - 2005
    !!
    integer, parameter,private   :: dp = selected_real_kind(14, 300)

    !!--++
    !!--++ SP
    !!--++    integer, parameter,private   :: sp = selected_real_kind(6, 30)
    !!--++
    !!--++    Private variable
    !!--++
    !!--++ Update: October - 2005
    !!
    integer, parameter,private   :: sp = selected_real_kind(6, 30)

    !!--++
    !!--++ INV_8LN2
    !!--++    real(kind=dp),private, parameter ::    inv_8ln2=0.18033688011112042591999058512524_dp
    !!--++
    !!--++    Private variable
    !!--++
    !!--++ Update: October - 2005
    !!
    real(kind=dp),private, parameter ::    inv_8ln2=0.18033688011112042591999058512524_dp

    !!--++
    !!--++ TWO_OVER_PI
    !!--++    real(kind=dp),private, parameter ::  two_over_pi=0.63661977236758134307553505349006_dp
    !!--++
    !!--++    Private variable
    !!--++
    !!--++ Update: October - 2005
    !!
    real(kind=dp),private, parameter :: two_over_pi=0.63661977236758134307553505349006_dp

    !!----
    !!---- LORCOMP
    !!----    logical, public  :: lorcomp
    !!----
    !!----    .true. if there are Lorentzian components
    !!----
    !!---- Update: October - 2005
    !!
    logical, public  :: lorcomp


    !!--.. The following public type contains the derivatives with respect to different
    !!--.. parameters of the peak shape functions
    !!--.. The calling program should define a variable of type deriv_TOF, if needed, and
    !!--.. invoke the functions with the appropriate arguments

    !!----
    !!---- TYPE :: DERIV_TOF
    !!--..
    !!---- Type, public :: deriv_TOF
    !!----    real(kind=sp) :: alfa     ! omega_a  DOmega/Dalpha
    !!----    real(kind=sp) :: beta     ! omega_b  DOmega/Dbeta
    !!----    real(kind=sp) :: dt       ! omega_t  DOmega/Ddt      (dt=TOFi-TOF(Bragg))
    !!----    real(kind=sp) :: sigma    ! omega_s  DOmega/Dsigma   (for tof_Jorgensen function)
    !!----    real(kind=sp) :: gamma    ! omega_g  DOmega/Dgamma   (for tof_Jorgensen_VonDreele function)
    !!----    real(kind=sp) :: eta      ! omega_e  DOmega/Deta                     "
    !!----    real(kind=sp) :: kappa    ! omega_e  DOmega/kappa    (for tof_Carpenter function)
    !!---- End Type deriv_TOF
    !!----
    !!---- Update: October - 2005
    !!
    Type, Public :: Deriv_TOF
       real(kind=sp) :: alfa     ! omega_a  DOmega/Dalpha
       real(kind=sp) :: beta     ! omega_b  DOmega/Dbeta
       real(kind=sp) :: dt       ! omega_t  DOmega/Ddt      (dt=TOFi-TOF(Bragg))
       real(kind=sp) :: sigma    ! omega_s  DOmega/Dsigma   (for tof_Jorgensen function)
       real(kind=sp) :: gamma    ! omega_g  DOmega/Dgamma   (for tof_Jorgensen_VonDreele function)
       real(kind=sp) :: eta      ! omega_e  DOmega/Deta                     "
       real(kind=sp) :: kappa    ! omega_e  DOmega/kappa    (for tof_Carpenter function)
    End Type Deriv_TOF

    !---- Interfaces ----!
    Interface  Erfc
       module procedure erfcc
       module procedure derfcc
    End Interface

    Interface  Erfcp
       module procedure erfccp
       module procedure derfccp
    End Interface


 Contains
    !!----
    !!---- Function Erfc(X) Result(Cerrf)
    !!----    real(kind=dp/sp), intent(in)  :: x
    !!----    real(kind=dp/sp)              :: cerrf
    !!----
    !!----    Complementary error function
    !!----
    !!---- Update: October - 2005
    !!

    !!--++
    !!--++ Function Derfcc(X) Result(Cerrf)
    !!--++    real(kind=dp), intent(in)  :: x
    !!--++    real(kind=dp)              :: cerrf
    !!--++
    !!--++    Complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Function Derfcc(X) Result(Cerrf)
       !---- Argument ----!
       real(kind=dp), intent(in)  :: x
       real(kind=dp)              :: cerrf

       !---- Local variables ----!
       real(kind=dp) :: t,z,y

       z=abs(x)
       t=1.0_dp/(1.0_dp+0.5_dp*z)
       y= -z*z-1.26551223_dp+t*(1.00002368_dp+t*(0.37409196_dp+  &
          t*(0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+t*(-1.13520398_dp+  &
          t*(1.48851587_dp+t*(-0.82215223_dp+t*0.17087277_dp))))))))
       cerrf=t*exp(y)
       if (x < 0.0_dp) cerrf=2.0_dp-cerrf

       return
    End Function Derfcc

    !!--++
    !!--++ Function Erfcc(X) Result(Cerrf)
    !!--++    real(kind=dp), intent(in)  :: x
    !!--++    real(kind=dp)              :: cerrf
    !!--++
    !!--++    Complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Function Erfcc(X) Result(Cerrf)
       !---- Argument ----!
       real(kind=sp), intent(in)  :: x
       real(kind=sp)              :: cerrf

       !---- Local Variables ----!
       real(kind=sp) :: t,z,y

       z=abs(x)
       t=1.0/(1.0+0.5*z)
       y= -z*z-1.26551223+t*(1.00002368+t*(0.37409196+  &
          t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+  &
          t*(1.48851587+t*(-0.82215223+t*0.17087277))))))))
       cerrf=t*exp(y)
       if (x < 0.0) cerrf=2.0-cerrf

       return
    End Function Erfcc

    !!----
    !!---- Function Erfcp(X) Result(Der)
    !!----    real(kind=dp/sp), intent(in)  :: x
    !!----    real(kind=dp/sp)              :: der
    !!----
    !!----    Derivative of the complementary error function
    !!----
    !!---- Update: October - 2005
    !!

    !!--++
    !!--++ Function Derfccp(X) Result(Der)
    !!--++    real(kind=dp), intent(in)  :: x
    !!--++    real(kind=dp)              :: der
    !!--++
    !!--++    Derivative of the complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Function Derfccp(X) Result(Der)
       !---- Argument ----!
       real(kind=dp), intent(in)    :: x
       real(kind=dp)                :: der

       der = - 1.128379167096_dp * exp(-x*x)   !from Abramovitz

       return
    End Function Derfccp

    !!--++
    !!--++ Function Erfccp(X) Result(Der)
    !!--++    real(kind=dp), intent(in)  :: x
    !!--++    real(kind=dp)              :: der
    !!--++
    !!--++    Derivative of the complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Function Erfccp(X)  Result(Der)
       !---- Argument ----!
       real(kind=sp), intent(in)    :: x
       real(kind=sp)                :: der

       der = - 1.128379167096 * exp(-x*x)

       return
    End Function Erfccp

    !!--++
    !!--++ Function Expi_e1(z) Result(Exp_e1)
    !!--++    complex(kind=dp), intent (in) ::  z
    !!--++    complex(kind=dp)              ::  exp_e1
    !!--++
    !!--++    Derivative of the complementary error function
    !!--++
    !!--++ Update: October - 2005
    !!
    Function Expi_E1(Z)  Result(Exp_E1)
       !---- Argument ----!
       complex(kind=dp), intent (in) ::  z
       complex(kind=dp)              ::  exp_e1

       !---- Local variables ----!
       integer                  :: k
       real(kind=dp), parameter :: pi=3.1415926535897932_dp
       real(kind=dp), parameter :: el=0.5772156649015328_dp
       real(kind=dp)            :: a0,x
       complex(kind=dp)         :: cr,ct0

       x=real(z)
       a0=abs(z)
       if (a0 <= 0.0_dp) then
          exp_e1=(1.0e+300_dp,0.0_dp)
       else if (a0 <= 10.0 .or. x <= 0.0 .and. a0 < 20.0) then
          exp_e1=(1.0_dp,0.0_dp)
          cr=(1.0_dp,0.0_dp)
          do k=1,150
             cr=-cr*k*z/(k+1.0_dp)**2
             exp_e1=exp_e1+cr
             if ( abs(cr) <=  abs(exp_e1)*1.0e-15_dp) exit
          end do
          exp_e1=(-el- log(z)+z*exp_e1)*exp(z)
       else
          ct0=(0.0_dp,0.0_dp)
          do k=120,1,-1
             ct0=k/(1.0_dp+k/(z+ct0))
          end do
          exp_e1=1.0_dp/(z+ct0)
          if (x <= 0.0_dp .and. abs(aimag(z)) <= 1.0e-9_dp)  &
                          exp_e1=exp_e1-pi*(0.0_dp,1.0_dp)*exp(z)
       end if

       return
    End Function Expi_E1

    !!----
    !!---- Subroutine Tof_Carpenter(Dt,D,Alfa,Beta,Gamma,Eta,Kappa,Tof_Theta,Tof_Peak,Deriv)
    !!----    real(kind=sp),             intent( in) :: dt        ! dt = TOF(channel i) -TOF(Bragg position)
    !!----    real(kind=sp),             intent( in) :: d         ! d-spacing of the peak in A
    !!----    real(kind=sp),             intent( in) :: alfa      !  alfa  : units microsecs-1
    !!----    real(kind=sp),             intent( in) :: beta      !  beta  : units microsecs-1
    !!----    real(kind=sp),             intent( in) :: gamma     !  gamma : units microsecs
    !!----    real(kind=sp),             intent( in) :: eta       !  eta   : mixing coefficient calculated using TCH
    !!----    real(kind=sp),             intent( in) :: kappa     ! Mixing coeficient of the Ikeda-Carpenter function
    !!----    real(kind=sp),             intent( in) :: tof_theta ! This is the value of 2sin(theta)
    !!----    real(kind=sp),             intent(out) :: tof_peak
    !!----    type(deriv_TOF), optional, intent(out) :: deriv     ! present if derivatives are to be calculated
    !!----
    !!---- Laurent C Chapon
    !!----
    !!---- Update: October - 2005
    !!
    Subroutine Tof_Carpenter(Dt,D,Alfa,Beta,Gamma,Eta,Kappa,Tof_Theta,Tof_Peak,Deriv)
       !---- Arguments ----!
       real(kind=sp),             intent( in) :: dt        ! dt = TOF(channel i) -TOF(Bragg position)
       real(kind=sp),             intent( in) :: d         ! d-spacing of the peak in A
       real(kind=sp),             intent( in) :: alfa      !  alfa  : units microsecs-1
       real(kind=sp),             intent( in) :: beta      !  beta  : units microsecs-1
       real(kind=sp),             intent( in) :: gamma     !  gamma : units microsecs
       real(kind=sp),             intent( in) :: eta       !  eta   : mixing coefficient calculated using TCH
       real(kind=sp),             intent( in) :: kappa     ! Mixing coeficient of the Ikeda-Carpenter function
       real(kind=sp),             intent( in) :: tof_theta ! This is the value of 2sin(theta)
       real(kind=sp),             intent(out) :: tof_peak
       type(deriv_TOF), optional, intent(out) :: deriv     ! present if derivatives are to be calculated

       !---- local variables ----!
       integer          :: i,udiv,vdiv,sdiv,rdiv
       complex(kind=dp) :: zu,zv,zs,zr,fzu,fzv,fzs,fzr
       real(kind=dp)    :: lambda,sigma,R,Nu,Nv,Ns,Nr,deno,denoinv,ki,alfa_min,alfa_plu, &
                           xik,yik,zik,norm,g_u,g_v,g_s,g_r,y_u,y_v,y_s,y_r, &
                           expu,expv,exps,expr,erfu,erfv,erfs,erfr,omg, &
                           oml_u,oml_v,oml_s,oml_r,one_e,oml, &
                           omg_u,omg_v,omg_s,omg_r,erfpu,erfpv,erfps,erfpr,a_u,a_v,a_s,a_r, &
                           domg_t, domg_a, domg_b, domg_g, domg_k, doml_t, doml_a, &
                           doml_b, doml_g, doml_k, dnuda,dnudb,dnvda,dnvdb,dnsda,dnsdb,dnrda,dnrdb, &
                           expru,exprv,exprs,exprr,im_const,re_const

       ! Definitions of all parameters for the Gaussian part :
       sigma=gamma*gamma*inv_8ln2   ! Define sigma with respect to gamma (sigma is the variance of sigma)
       lambda=d*tof_theta
       R=exp(-81.799_dp/(kappa*lambda*lambda))
       deno=SQRT(2.0_dp*sigma)
       denoinv=1.0_dp/deno
       ki=0.05_dp
       alfa_min=alfa*(1.0_dp-ki)
       alfa_plu=alfa*(1.0_dp+ki)
       xik=alfa_min-beta
       yik=alfa-beta
       zik=alfa_plu-beta
       norm=0.25_dp*alfa*(1-ki*ki)/ki/ki
       Nu=1.0_dp-(R*alfa_min/xik)
       Nv=1.0_dp-(R*alfa_plu/zik)
       Ns=-2.0_dp*(1-R*alfa/yik)
       Nr=2.0_dp*R*alfa*alfa*beta*ki*ki/xik/yik/zik
       y_u=(alfa_min*sigma-dt)*denoinv
       y_v=(alfa_plu*sigma-dt)*denoinv
       y_s=(alfa*sigma-dt)*denoinv
       y_r=(beta*sigma-dt)*denoinv
       g_u=0.5_dp*alfa_min*(alfa_min*sigma-2.0*dt)
       g_v=0.5_dp*alfa_plu*(alfa_plu*sigma-2.0*dt)
       g_s=0.5_dp*alfa*(alfa*sigma-2.0*dt)
       g_r=0.5_dp*beta*(beta*sigma-2.0*dt)
       udiv=max(1,nint(0.003*g_u))
       vdiv=max(1,nint(0.003*g_v))
       sdiv=max(1,nint(0.003*g_s))
       rdiv=max(1,nint(0.003*g_r))
       g_u=g_u/udiv
       g_v=g_v/vdiv
       g_s=g_s/sdiv
       g_r=g_r/rdiv
       ! End of definitions

       expu=exp(g_u)
       expv=exp(g_v)
       exps=exp(g_s)
       expr=exp(g_r)
       erfu=erfc(y_u)
       erfv=erfc(y_v)
       erfs=erfc(y_s)
       erfr=erfc(y_r)

       if (y_u < 27) then      ! LCC , november 2003: Modifications to the original code
          omg_u=expu*erfu      ! to solve underflow problems occuring
          do i=1, udiv-1       ! when erfc(y_*) becomes too small and exp(g_*) too large,
             omg_u=omg_u*expu  ! --> Use an approximation at third order to erfc.
          end do
       else
          omg_u=exp(g_u*udiv-y_u**2)/(sqrt(2/two_over_pi)*y_u)* &
          (1-1.0/(2.0*y_u**2)+3.0/(2.0*y_u**2)**2+15.0/(2.0*y_u**2)**3)
       endif

       if (y_v < 27) then
          omg_v=expv*erfv
          do i=1, vdiv-1
             omg_v=omg_v*expv
          end do
       else
          omg_v=exp(g_v*vdiv-y_v**2)/(sqrt(2/two_over_pi)*y_v)* &
          (1-1.0/(2.0*y_v**2)+3.0/(2.0*y_v**2)**2+15.0/(2.0*y_v**2)**3)
       endif

       if (y_s < 27) then
          omg_s=exps*erfs
          do i=1, sdiv-1
             omg_s=omg_s*exps
          end do
       else
          omg_s=exp(g_s*sdiv-y_s**2)/(sqrt(2.0/two_over_pi)*y_s)* &
          (1.0-1.0/(2.0*y_s**2)+3.0/(2.0*y_s**2)**2+15.0/(2.0*y_s**2)**3)
       endif

       if (y_r < 27) then
          omg_r=expr*erfr
          do i=1, rdiv-1
             omg_r=omg_r*expr
          end do
       else
          omg_r=exp(g_r*rdiv-y_r**2)/(sqrt(2.0/two_over_pi)*y_r)* &
          (1.0-1.0/(2.0*y_r**2)+3.0/(2.0*y_r**2)**2+15.0/(2.0*y_r**2)**3)
       endif

       omg=(Nu*omg_u+Nv*omg_v+Ns*omg_s+Nr*omg_r) ! End of gaussian part of the function

       if (lorcomp) then
          zs=CMPLX(-alfa*dt,0.5*alfa*gamma,kind=dp)
          zu=(1.0-ki)*zs
          zv=(1.0+ki)*zs
          zr=CMPLX(-beta*dt,0.5*beta*gamma,kind=dp)
          fzu=expi_e1(zu)
          fzv=expi_e1(zv)
          fzs=expi_e1(zs)
          fzr=expi_e1(zr)
          oml_u=-AIMAG(fzu)*two_over_pi
          oml_v=-AIMAG(fzv)*two_over_pi
          oml_s=-AIMAG(fzs)*two_over_pi
          oml_r=-AIMAG(fzr)*two_over_pi
          oml=Nu*oml_u+Nv*oml_v+Ns*oml_s+Nr*oml_r   ! End of Lorentzian part of the function
          one_e=1.0-eta
          tof_peak=norm*(one_e*omg+eta*oml)         ! Total function
       else
          tof_peak=norm*omg
       end if

       ! Derivatives

       if(.not. present(deriv)) return

       ! Derivatives of Omega(Gaussian)

       if(omg <= 1.0E-35) then
          domg_t= 0.0           ! DOmG/Ddt
          domg_a= 0.0           ! DOmG/Dalfa
          domg_b= 0.0           ! DOmG/Dbeta
          domg_g= 0.0           ! DOmG/Dgamma
          domg_k= 0.0           ! DOmG/Dkappa
       else
          dnuda=R*beta/xik/xik*(1-ki)  ! Partial derivatives of Nu,Nv,Ns and Nr /dalpha, dbeta
          dnudb=-R*alfa_min/xik/xik
          dnvda=R*beta/zik/zik*(1+ki)
          dnvdb=-R*alfa_plu/zik/zik
          dnsda=-2*R*beta/yik/yik
          dnsdb=2*R*alfa/yik/yik
          dnrda=-2*R*alfa*beta*ki*ki*(alfa**3-3.0*alfa*beta**2-alfa**3*ki**2+2.0*beta**3)/xik**2/yik**2/zik**2
          dnrdb=-dnrda*alfa/beta
          erfpu=erfcp(y_u)                         ! d(erfc(y_u))/dy_u
          erfpv=erfcp(y_v)                         ! d(erfc(y_v))/dy_v
          erfps=erfcp(y_s)                         ! d(erfc(y_s))/dy_s
          erfpr=erfcp(y_r)                         ! d(erfc(y_r))/dy_r
          a_u=expu*erfpu
          a_v=expv*erfpv
          a_s=exps*erfps
          a_r=expr*erfpr
          do i=1, udiv-1
            a_u=a_u*expu
          end do
          do i=1, vdiv-1
            a_v=a_v*expv
          end do
          do i=1, sdiv-1
            a_s=a_s*exps
          end do
          do i=1, rdiv-1
            a_r=a_r*expr
          end do
          domg_t=-Nu*alfa_min*omg_u-Nv*alfa_plu*omg_v-Ns*alfa*omg_s-Nr*beta*omg_r -  &
                  (Nu*a_u+Nv*a_v+Ns*a_s+Nr*a_r)*denoinv
          domg_a= dnuda*omg_u+dnvda*omg_v+dnsda*omg_s+dnrda*omg_r + &
                  deno*(Nu*(1-ki)*(y_u*omg_u+0.5_dp*a_u)+Nv*(1+ki)*(y_v*omg_v+0.5_dp*a_v)+Ns*(y_s*omg_s+0.5_dp*a_s))
          domg_b= dnudb*omg_u+dnvdb*omg_v+dnsdb*omg_s+dnrdb*omg_r+deno*Nr*(y_r*omg_r+0.5_dp*a_r)
          domg_g= inv_8ln2*gamma*(Nu*alfa_min**2*omg_u+Nv*alfa_plu**2*omg_v+Ns*alfa**2*omg_s+Nr*beta**2*omg_r) + &
                  2.0*inv_8ln2*gamma*denoinv**3*(Nu*(alfa_min*sigma+dt)*a_u+Nv*(alfa_plu*sigma+dt)*a_v+          &
                  Ns*(alfa*sigma+dt)*a_s+Nr*(beta*sigma+dt)*a_r)
          domg_k= 81.799_dp/kappa**2/lambda**2*R*(-alfa_min/xik*omg_u-alfa_plu/zik*omg_v+2.0*alfa/yik*omg_s+     &
                  2.0*alfa**2*beta*ki**2/xik/yik/zik*omg_r)
       end if

       if(lorcomp) then

        ! Derivatives of Omega(Lorentzian)
          im_const=AIMAG(alfa/zs)
          re_const=REAL(alfa/zs)
          expru=REAL(fzu)
          exprv=REAL(fzv)
          exprs=REAL(fzs)
          exprr=REAL(fzr)
          doml_t=-Nu*alfa_min*oml_u-Nv*alfa_plu*oml_v-Ns*alfa*oml_s-Nr*beta*oml_r-two_over_pi*im_const*(Nu+Nv+Ns+Nr)
          doml_a=dnuda*oml_u+dnvda*oml_v+dnsda*oml_s+dnrda*oml_r+Nu*(1.0-ki)*(-dt*oml_u-0.5*two_over_pi*gamma*expru)+ &
                 Nv*(1.0+ki)*(-dt*oml_v-0.5*two_over_pi*gamma*exprv)+Ns*(-dt*oml_s-0.5*two_over_pi*gamma*exprs)
          doml_b=dnudb*oml_u+dnvdb*oml_v+dnsdb*oml_s+dnrdb*oml_r+Nr*(-dt*oml_r-0.5*two_over_pi*gamma*exprr)
          doml_g=0.5_dp*two_over_pi*(-Nu*alfa_plu*expru-Nv*alfa_min*exprv-Ns*alfa*exprs-Nr*beta*exprr+(Nu+Nv+Ns+Nr)*re_const)
          doml_k=81.799_dp/kappa**2/lambda**2*R*(-alfa_min/xik*oml_u-alfa_plu/zik*oml_v+2.0*alfa/yik*oml_s+  &
                 2.0*alfa**2*beta*ki**2/xik/yik/zik*oml_r)

        ! Total derivatives

          deriv%dt    = norm*(one_e*domg_t+eta*doml_t)
          deriv%alfa  = tof_peak/alfa+norm*(one_e*domg_a+eta*doml_a)
          deriv%beta  = norm*(one_e*domg_b+eta*doml_b)
          deriv%gamma = norm*(one_e*domg_g+eta*doml_g)
          deriv%kappa = norm*(one_e*domg_k+eta*doml_k)
          deriv%eta   = norm*(oml-omg)

       else

          deriv%dt    = norm*domg_t
          deriv%alfa  = tof_peak/alfa+norm*domg_a
          deriv%beta  = norm*domg_b
          deriv%gamma = norm*domg_g
          deriv%kappa = norm*domg_k
          deriv%eta   = 0.0
       end if
       deriv%sigma = deriv%gamma/(2.0*inv_8ln2*gamma)

       return
    End Subroutine Tof_Carpenter

    !!----
    !!---- Subroutine Tof_Jorgensen(Dt,Alfa,Beta,Sigma,Tof_Peak,Deriv)
    !!----    real(kind=sp),             intent( in)  :: dt       !  dt = TOF(channel i) -TOF(Bragg position): units microsecs
    !!----    real(kind=sp),             intent( in)  :: alfa     !  alfa  : units microsecs-1
    !!----    real(kind=sp),             intent( in)  :: beta     !  beta  : units microsecs-1
    !!----    real(kind=sp),             intent( in)  :: sigma    !  sigma : units microsecs**2
    !!----    real(kind=sp),             intent(out)  :: tof_peak
    !!----    type(deriv_TOF), optional, intent(out)  :: deriv    ! present if derivatives are to be calculated
    !!----
    !!---- Laurent C Chapon
    !!----
    !!---- Update: October - 2005
    !!
    Subroutine Tof_Jorgensen(Dt,Alfa,Beta,Sigma,Tof_Peak,Deriv)
       !---- Arguments ----!
       real(kind=sp),             intent( in)  :: dt       !  dt = TOF(channel i) -TOF(Bragg position): units microsecs
       real(kind=sp),             intent( in)  :: alfa     !  alfa  : units microsecs-1
       real(kind=sp),             intent( in)  :: beta     !  beta  : units microsecs-1
       real(kind=sp),             intent( in)  :: sigma    !  sigma : units microsecs**2
       real(kind=sp),             intent(out)  :: tof_peak
       type(deriv_TOF), optional, intent(out)  :: deriv    ! present if derivatives are to be calculated

       !---- Local Variables ----!
       integer        :: i, udiv, vdiv
       real(kind=dp)  :: u,v,y,z,expu,expv,erfy,erfz,a,b,omega,omegb
       real(kind=dp)  :: norm, deno, denoinv,ca,cb,a2,b2,d3,omeg
       real(kind=dp)  :: afpbet,erfyp,erfzp

       u=0.5*alfa*(alfa*sigma+2.0*dt)
       v=0.5*beta*(beta*sigma-2.0*dt)
       ! To avoid pathological behaviour EXP(U) is calculated as
       ! {EXP(U/N)}^N using mutiplicative loops
       udiv=max(1,nint(0.003*u))   !udiv=max(1,int(2.0*u))
       vdiv=max(1,nint(0.003*v))   !vdiv=max(1,int(2.0*v))
       u=u/udiv
       v=v/vdiv
       afpbet=alfa+beta
       norm=0.5*alfa*beta/afpbet
       deno=SQRT(2.0_dp*sigma)
       denoinv=1.0_dp/deno
       y=(alfa*sigma+dt)*denoinv
       z=(beta*sigma-dt)*denoinv

       expu=exp(u)
       erfy=erfc(y)              !error function of y
       omega=expu*erfy
       do i=1,udiv-1
          omega=omega*expu
       end do
       expv=exp(v)
       erfz=erfc(z)             !error function of z
       omegb=expv*erfz
       do i=1,vdiv-1
          omegb=omegb*expv
       end do
       omeg=norm*(omega+omegb)
       tof_peak=omeg

       ! Derivatives
       if(.not. present(deriv)) return

       ! Derivatives of Omega(Gaussian)
       if(omeg <= 1.0e-30) then
          deriv%dt   = 0.0          ! DOmeG/Ddt
          deriv%alfa = 0.0          ! DOmeG/Dalfa
          deriv%beta = 0.0          ! DOmeG/Dbeta
          deriv%sigma= 0.0          ! DOmeG/Dsigma
       else
          erfyp=erfcp(y)         ! erfcc'(y)
          erfzp=erfcp(z)         ! erfcc'(z)
          a=expu*erfyp
          b=expv*erfzp
          do i=1,udiv-1
           a=a*expu              !a=exp(u)*erfcc'(y)
          end do
          do i=1,vdiv-1
           b=b*expv              !b=exp(v)*erfcc'(z)
          end do
          a2=alfa*alfa
          b2=beta*beta
          d3= 2.0_dp*(denoinv*denoinv*denoinv)
          ca=d3*(alfa*sigma-dt)
          cb=d3*(beta*sigma+dt)
          deriv%dt    = norm*(alfa*omega-beta*omegb+(a-b)*denoinv)
          deriv%alfa  = norm*(2.0_dp*omeg/a2+deno*(y*omega+0.5_dp*a))
          deriv%beta  = norm*(2.0_dp*omeg/b2+deno*(z*omegb+0.5_dp*b))
          deriv%sigma = 0.5_dp*norm*(a2*omega+b2*omegb+a*ca+b*cb)
          deriv%gamma = 2.0*deriv%sigma*sqrt(sigma*inv_8ln2)
       end if

       return
    End Subroutine Tof_Jorgensen

    !!----
    !!---- Subroutine Tof_Jorgensen(Dt,Alfa,Beta,Sigma,Tof_Peak,Deriv)
    !!----    real(kind=sp),             intent( in) :: dt       ! dt = TOF(channel i) -TOF(Bragg position)
    !!----    real(kind=sp),             intent( in) :: alfa     !  alfa  : units microsecs-1
    !!----    real(kind=sp),             intent( in) :: beta     !  beta  : units microsecs-1
    !!----    real(kind=sp),             intent( in) :: gamma    !  gamma : units microsecs
    !!----    real(kind=sp),             intent( in) :: eta      !  eta   : mixing coefficient calculated using TCH
    !!----    real(kind=sp),             intent(out) :: tof_peak
    !!----    type(deriv_TOF), optional, intent(out) :: deriv    ! present if derivatives are to be calculated
    !!----
    !!---- Laurent C Chapon
    !!----
    !!---- Update: October - 2005
    !!
    Subroutine Tof_Jorgensen_Vondreele(Dt,Alfa,Beta,Gamma,Eta,Tof_Peak,Deriv)
       !---- Arguments ----!
       real(kind=sp),             intent( in) :: dt       ! dt = TOF(channel i) -TOF(Bragg position)
       real(kind=sp),             intent( in) :: alfa     !  alfa  : units microsecs-1
       real(kind=sp),             intent( in) :: beta     !  beta  : units microsecs-1
       real(kind=sp),             intent( in) :: gamma    !  gamma : units microsecs
       real(kind=sp),             intent( in) :: eta      !  eta   : mixing coefficient calculated using TCH
       real(kind=sp),             intent(out) :: tof_peak
       type(deriv_TOF), optional, intent(out) :: deriv    ! present if derivatives are to be calculated

       !---- local variables ----!
       complex(kind=dp):: z1,z2,fz1,fz2
       real(kind=dp)   :: u,v,expu,expv,a, b, erfy, erfz, y,z
       real(kind=dp)   :: norm, deno,denoinv,ca,cb,a2,b2,d3, sigma
       real(kind=dp)   :: omeg, omel, afpbet,oml_a,oml_b,one_e, &
                          erfyp,erfzp,domg_t,domg_a,domg_b,domg_g,exper1,exper2,  &
                          doml_t,doml_a,doml_b,doml_g,omega,omegb
       integer         :: i, udiv, vdiv

       sigma=gamma*gamma*inv_8ln2
       u=0.5*alfa*(alfa*sigma+2.0*dt)
       v=0.5*beta*(beta*sigma-2.0*dt)
       ! To avoid pathological behaviour EXP(U) is calculated as
       ! {EXP(U/N)}^N using mutiplicative loops
       udiv=max(1,nint(0.003*u))
       vdiv=max(1,nint(0.003*v))
       u=u/udiv
       v=v/vdiv
       afpbet=alfa+beta
       norm=0.5*alfa*beta/afpbet
       deno=SQRT(2.0_dp*sigma)
       denoinv=1.0_dp/deno
       y=(alfa*sigma+dt)*denoinv
       z=(beta*sigma-dt)*denoinv

       expu=exp(u)
       erfy=erfc(y)              !error function of y
       omega=expu*erfy
       do i=1,udiv-1
          omega=omega*expu
       end do
       expv=exp(v)
       erfz=erfc(z)             !error function of z
       omegb=expv*erfz
       do i=1,vdiv-1
          omegb=omegb*expv
       end do

       omeg=norm*(omega+omegb)    !Gaussian contribution

       if(lorcomp) then
          z1=CMPLX( alfa*dt,0.5*alfa*gamma,kind=dp)
          z2=CMPLX(-beta*dt,0.5*beta*gamma,kind=dp)
          fz1=expi_e1(z1)                  ! exp(p).E1(p)
          fz2=expi_e1(z2)                  ! exp(q).E1(q)
          oml_a=-AIMAG(fz1)*two_over_pi    ! OmL,alfa
          oml_b=-AIMAG(fz2)*two_over_pi    ! OmL,beta
          omel= norm*(oml_a+oml_b)         ! OmL = OmL,alfa + OmL,beta (Lorentzian contribution)
          one_e=1.0-eta
          tof_peak = one_e * omeg + eta * omel
       else
          tof_peak = omeg
       end if

       ! Derivatives
       if(.not. present(deriv)) return

       ! Derivatives of Omega(Gaussian)
       a2=alfa*alfa
       b2=beta*beta
       d3= 2.0_dp*(denoinv*denoinv*denoinv)
       ca=d3*(alfa*sigma-dt)
       cb=d3*(beta*sigma+dt)

       if(omeg <= 1.e-30) then
          domg_t= 0.0           ! DOmG/Ddt
          domg_a= 0.0           ! DOmG/Dalfa
          domg_b= 0.0           ! DOmG/Dbeta
          domg_g= 0.0           ! DOmG/Dgamma
       else
          erfyp=erfcp(y)         ! erfcc'(y)
          erfzp=erfcp(z)         ! erfcc'(z)
          a=expu*erfyp
          b=expv*erfzp
          do i=1,udiv-1
           a=a*expu              !a=exp(u)*erfcc'(y)
          end do
          do i=1,vdiv-1
           b=b*expv              !b=exp(v)*erfcc'(z)
          end do
          domg_t = norm*(alfa*omega-beta*omegb+(a-b)*denoinv)
          domg_a = norm*(2.0_dp*omeg/a2+deno*(y*omega+0.5_dp*a))
          domg_b = norm*(2.0_dp*omeg/b2+deno*(z*omegb+0.5_dp*b))
          !Multiply by Dsigma/Dgamma=gamma/4ln2
          domg_g = inv_8ln2 * norm * gamma * (a2*omega+b2*omegb+a*ca+b*cb)
       end if

       if(lorcomp) then

        ! Derivatives of Omega(Lorentzian)

          exper1=-REAL(fz1)*two_over_pi
          exper2=-REAL(fz2)*two_over_pi
          doml_t= norm*(alfa* oml_a - beta* oml_b)                               ! DOmL/Ddt
          doml_a= norm*( 2.0_dp * omel/a2 + dt* oml_a + 0.5_dp * gamma * exper1) ! DOmL/Dalfa
          doml_b= norm*( 2.0_dp * omel/b2 - dt* oml_b + 0.5_dp * gamma * exper2) ! DOmL/Dbeta
          doml_g= 0.5 * norm * ( alfa*exper1 + beta*exper2 )                     ! DOmL/Dgamma

        ! Total derivatives

          deriv%dt    = one_e * domg_t + eta * doml_t
          deriv%alfa  = one_e * domg_a + eta * doml_a
          deriv%beta  = one_e * domg_b + eta * doml_b
          deriv%gamma = one_e * domg_g + eta * doml_g
          deriv%eta   = omel-omeg
       else
          deriv%dt    = domg_t
          deriv%alfa  = domg_a
          deriv%beta  = domg_b
          deriv%gamma = domg_g
          deriv%eta   = 0.0
       end if
       deriv%sigma = deriv%gamma/(2.0*inv_8ln2*gamma)

       return
    End Subroutine Tof_Jorgensen_Vondreele

 End Module  Tof_Profiles
