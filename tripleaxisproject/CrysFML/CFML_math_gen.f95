!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: MATH_GEN
!!----   INFO: Mathematic general utilities for use in Crystallography and
!!----         Solid State Physics and Chemistry.
!!----
!!---- HISTORY
!!----    Update: February - 2005  JGP
!!----            August   - 1990  Based in public codes. Created by JRC
!!----
!!---- DEPENDENCIES
!!----
!!----    MOD_FUN    !To be commented for non-F compilers
!!----
!!---- VARIABLES
!!----    DP
!!----    SP
!!----    CP
!!----    DEPS
!!----    EPS
!!--++    EPSS                         [Private]
!!----    ERR_MATH_GEN
!!----    ERR_MESS_MATH_GEN
!!----    PI
!!----    TO_DEG
!!----    TO_RAD
!!----    TPI
!!----
!!---- PROCEDURES
!!----    Functions:
!!--..
!!--..    Trigonometric Functions
!!----       ACOSD
!!--++       ACOSD_dp                  [Overloaded]
!!--++       ACOSD_sp                  [Overloaded]
!!----       ASIND
!!--++       ASIND_dp                  [Overloaded]
!!--++       ASIND_sp                  [Overloaded]
!!----       ATAN2D
!!--++       ATAN2D_dp                 [Overloaded]
!!--++       ATAN2D_sp                 [Overloaded]
!!----       ATAND
!!--++       ATAND_dp                  [Overloaded]
!!--++       ATAND_sp                  [Overloaded]
!!----       COSD
!!--++       COSD_dp                   [Overloaded]
!!--++       COSD_sp                   [Overloaded]
!!----       SIND
!!--++       SIND_dp                   [Overloaded]
!!--++       SIND_sp                   [Overloaded]
!!----       TAND
!!--++       TAND_dp                   [Overloaded]
!!--++       TAND_sp                   [Overloaded]
!!--..
!!--..    Special Functions
!!----       BESSJ0
!!----       BESSJ1
!!--..
!!--..    Scalar Functions
!!----       FACTORIAL
!!----       NEGLIGIBLE
!!--++       NEGLIGIBLEC               [Overloaded]
!!--++       NEGLIGIBLER               [Overloaded]
!!----       PGCD
!!----       PPCM
!!----       PYTHAG
!!--++       PYTHAG_dp                 [Overloaded]
!!--++       PYTHAG_sp                 [Overloaded]
!!--..
!!--..    Arrays and Vectors Functions
!!----       CO_LINEAR
!!--++       CO_LINEAR_C               [Overloaded]
!!--++       CO_LINEAR_I               [Overloaded]
!!--++       CO_LINEAR_R               [Overloaded]
!!----       EQUAL_MATRIX
!!--++       EQUAL_MATRIX_I            [Overloaded]
!!--++       EQUAL_MATRIX_R            [Overloaded]
!!----       EQUAL_VECTOR
!!--++       EQUAL_VECTOR_I            [Overloaded]
!!--++       EQUAL_VECTOR_R            [Overloaded]
!!----       IMAXLOC
!!--++       IMAXLOC_I                 [Overloaded]
!!--++       IMAXLOC_R                 [OVerloaded]
!!----       IMINLOC
!!--++       IMINLOC_I                 [Overloaded]
!!--++       IMINLOC_R                 [OVerloaded]
!!----       LOCATE
!!--++       LOCATE_I                  [Overloaded]
!!--++       LOCATE_R                  [Overloaded]
!!----       MODULO_LAT
!!----       OUTERPROD
!!--++       OUTERPROD_dp              [Overloaded]
!!--++       OUTERPROD_sp              [Overloaded]
!!----       TRAZA
!!--++       TRAZA_C                   [Overloaded]
!!--++       TRAZA_I                   [Overloaded]
!!--++       TRAZA_R                   [Overloaded]
!!----       ZBELONG
!!--++       ZBELONGM                  [Overloaded]
!!--++       ZBELONGN                  [Overloaded]
!!--++       ZBELONGV                  [Overloaded]

!!----
!!----    Subroutines:
!!--..
!!--..    Init Routine
!!----       INIT_ERR_MATHGEN
!!----       SET_EPSG
!!----       SET_EPSG_DEFAULT
!!--..
!!--..    Trigonometric Subroutines
!!----       RTAN
!!--++       RTAN_dp                   [Overloaded]
!!--++       RTAN_sp                   [Overloaded]
!!--..
!!--..    Arrays and Vectors Functions
!!----       DETERMINANT
!!--++       DETERMINANT_C             [Overloaded]
!!--++       DETERMINANT_R             [Overloaded]
!!----       DIAGONALIZE_SH
!!--++       DIAGONALIZE_HERM          [Overloaded]
!!--++       DIAGONALIZE_SYMM          [Overloaded]
!!--++       EIGSRT                    [Private]
!!----       INVERT_MATRIX
!!----       LINEAR_DEPENDENT
!!--++       LINEAR_DEPENDENTC         [Overloaded]
!!--++       LINEAR_DEPENDENTI         [Overloaded]
!!--++       LINEAR_DEPENDENTR         [Overloaded]
!!----       LU_BACKSUB
!!----       LU_DECOMP
!!----       MATINV
!!--++       PARTITION                 [Private]
!!----       RANK
!!--++       RANK_dp                   [Overloaded]
!!--++       RANK_sp                   [Overloaded]
!!----       SORT
!!--++       SORT_I                    [Overloaded]
!!--++       SORT_R                    [Overloaded]
!!----       SORT_STRINGS
!!----       SPLINE
!!----       SPLINT
!!----       SVDCMP
!!--++       SVDCMP_dp                 [Overloaded]
!!--++       SVDCMP_sp                 [Overloaded]
!!----       SWAP
!!--++       SWAP_C                    [Overloaded]
!!--++       SWAP_CM                   [Overloaded]
!!--++       SWAP_CV                   [Overloaded]
!!--++       SWAP_I                    [Overloaded]
!!--++       SWAP_IM                   [Overloaded]
!!--++       SWAP_IV                   [Overloaded]
!!--++       SWAP_R                    [Overloaded]
!!--++       SWAP_RM                   [Overloaded]
!!--++       SWAP_RV                   [Overloaded]
!!--++       MASKED_SWAP_R             [Overloaded]
!!--++       MASKED_SWAP_RM            [Overloaded]
!!--++       MASKED_SWAP_RV            [Overloaded]
!!--++       TQLI1                     [Private]
!!--++       TQLI2                     [Private]
!!--++       TRED1                     [Private]
!!--++       TRED2                     [Private]
!!--++
!!
 Module Math_Gen

    !---- Use Modules ----!
    !Use Mod_fun    !To be commented for non-F compilers

    !---- Variables ----!
    implicit none
    private

    !---- List of public functions ----!
    public :: Bessj0, Bessj1, Factorial, Pgcd, Ppcm, Modulo_Lat

    !---- List of public overloaded procedures: functions ----!
    public :: Acosd, Asind, Atan2d, Atand, Cosd, Sind, Tand, Negligible, Pythag,  &
              Co_Linear, Equal_Matrix, Equal_Vector, Locate, Outerprod, Traza,     &
              Zbelong, Imaxloc, Iminloc

    !---- List of private functions ----!
    private :: Acosd_dp, Acosd_sp, Asind_dp, Asind_sp, Atan2d_dp, Atan2d_sp,      &
               Atand_dp, Atand_sp, Cosd_dp, Cosd_sp, Sind_dp, Sind_sp, Tand_dp,   &
               Tand_sp, Negligiblec, Negligibler, Pythag_dp, Pythag_sp,           &
               Co_linear_C, Co_linear_I, Co_linear_R, Equal_Matrix_I,              &
               Equal_Matrix_R, Equal_Vector_I, Equal_Vector_R, Locate_I, Locate_R, &
               Outerprod_dp, Outerprod_sp, Traza_C, Traza_I, Traza_R, ZbelongM,   &
               ZbelongN, ZbelongV, Imaxloc_I, Imaxloc_R, Iminloc_R, Iminloc_I

    !---- List of public subroutines ----!
    public ::  Init_Err_Mathgen, Invert_Matrix, LU_Decomp, LU_Backsub, Matinv,     &
               Sort_Strings, Spline, Splint, Set_Epsg, Set_Epsg_Default

    !---- List of public overloaded procedures: subroutines ----!
    public ::  RTan, Determinant, Diagonalize_Sh, Linear_Dependent, Rank, Sort,   &
               Svdcmp, Swap

    !---- List of private subroutines ----!
    private :: RTan_dp, RTan_sp, Determinant_C,Determinant_R, Diagonalize_Herm,   &
               Diagonalize_Symm, Eigsrt, Linear_DependentC, Linear_DependentI,    &
               Linear_DependentR, Rank_dp, Rank_sp, Sort_I, Sort_R, Svdcmp_dp,    &
               Svdcmp_sp, Swap_C, Swap_Cm, Swap_Cv, Swap_I, Swap_Im, Swap_Iv,     &
               Swap_R, Swap_Rm, Swap_Rv, Masked_Swap_R, Masked_Swap_Rm,           &
               Masked_Swap_Rv, Tqli1, Tqli2, Tred1, Tred2, Partition

    !---- Definitions ----!
    !!----
    !!---- DP
    !!----    DP: Double precision ( dp = selected_real_kind(14,80) )
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: dp = selected_real_kind(14,80)

    !!----
    !!---- SP
    !!----    SP: Single precision ( sp = selected_real_kind(6,30) )
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: sp = selected_real_kind(6,30)

    !!----
    !!---- CP
    !!----    CP: Current precision
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: cp = sp

    !!----
    !!---- DEPS
    !!----    real(kind=dp), parameter :: deps=0.00000001_dp
    !!----
    !!----    Epsilon value
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=dp), parameter, public :: deps=0.00000001_dp

    !!----
    !!----  EPS
    !!----     real(kind=sp), public ::  eps=0.00001_sp
    !!----
    !!----     Epsilon value
    !!----
    !!----  Update: February - 2005
    !!
    real(kind=sp),  parameter, public  ::  eps=0.00001_sp

    !!--++
    !!--++ EPSS
    !!--++    real(kind=sp)  :: epss=1.0E-5_sp
    !!--++
    !!--++    Internal epsilon value used for comparing reals to integers
    !!--++    in crystallographic applications where the maximum precision in the
    !!--++    measured values is of the order of 10^-5.
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=sp),   private :: epss=1.0E-5_sp

    !!--++
    !!--++ EP_SS
    !!--++    real(kind=sp), parameter, private  :: ep_ss=1.0E-12_sp
    !!--++
    !!--++    Internal epsilon value used for comparison in matrix operations
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private :: ep_ss=1.0E-12_sp
    !!----
    !!---- ERR_MATH_GEN
    !!----    logical :: err_math_gen
    !!----
    !!----    Logical Variable indicating an error in MATH_GEN module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_math_gen

    !!----
    !!---- ERR_MESS_MATH_GEN
    !!----    character(len=150) :: err_mess_math_gen
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public:: err_mess_math_gen

    !!----
    !!---- PI
    !!----    real(kind=dp), parameter ::  pi = 3.141592653589793238463_dp
    !!----
    !!----    Pi value
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=dp), parameter, public    ::  pi = 3.141592653589793238463_dp

    !!----
    !!---- TO_DEG
    !!----    real(kind=dp), parameter ::  to_DEG = 180.0_dp/pi
    !!----
    !!----    Conversion from Rad to Degree
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=dp), parameter, public    ::  to_DEG  = 180.0_dp/pi

    !!----
    !!---- TO_RAD
    !!----    real(kind=dp), parameter ::  to_RAD  = pi/180.0_dp
    !!----
    !!----    Conversion from Degree to Rad
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=dp), parameter, public    ::  to_RAD  = pi/180.0_dp

    !!----
    !!---- TPI
    !!----  real(kind=dp), parameter ::  tpi = 6.283185307179586476925_dp
    !!----
    !!----  2.0*Pi value
    !!----
    !!---- Update: February - 2005
    !!
    real(kind=dp), parameter, public ::  tpi = 6.283185307179586476925_dp


    !---- Interfaces - Overloaded ----!
    Interface  Acosd
       Module Procedure Acosd_dp
       Module Procedure Acosd_sp
    End Interface

    Interface  Asind
       Module Procedure Asind_dp
       Module Procedure Asind_sp
    End Interface

    Interface  Atan2d
       Module Procedure Atan2d_dp
       Module Procedure Atan2d_sp
    End Interface

    Interface  Atand
       Module Procedure Atand_dp
       Module Procedure Atand_sp
    End Interface

    Interface  Cosd
       Module Procedure Cosd_dp
       Module Procedure Cosd_sp
    End Interface

    Interface  Sind
       Module Procedure Sind_dp
       Module Procedure Sind_sp
    End Interface

    Interface  Tand
       Module Procedure Tand_dp
       Module Procedure Tand_sp
    End Interface

    Interface  Negligible
       Module Procedure Negligibler
       Module Procedure Negligiblec
    End Interface

    Interface  Pythag
       Module Procedure Pythag_dp
       Module Procedure Pythag_sp
    End Interface

    Interface  Co_Linear
       Module Procedure Co_linear_C
       Module Procedure Co_linear_I
       Module Procedure Co_linear_R
    End Interface

    Interface  Equal_Matrix
       Module Procedure Equal_Matrix_I
       Module Procedure Equal_Matrix_R
    End Interface

    Interface  Equal_Vector
       Module Procedure Equal_Vector_I
       Module Procedure Equal_Vector_R
    End Interface

    Interface  IMaxloc
       Module Procedure IMaxloc_I
       Module Procedure IMaxloc_R
    End Interface

    Interface  IMinloc
       Module Procedure IMinloc_I
       Module Procedure IMinloc_R
    End Interface

    Interface  Locate
       Module Procedure Locate_I
       Module Procedure Locate_R
    End Interface

    Interface  Outerprod
       Module Procedure Outerprod_dp
       Module Procedure Outerprod_sp
    End Interface

    Interface  Traza
       Module Procedure Traza_C
       Module Procedure Traza_I
       Module Procedure Traza_R
    End Interface

    Interface  Zbelong
       Module Procedure ZbelongM
       Module Procedure ZbelongN
       Module Procedure ZbelongV
    End Interface

    Interface  Rtan
       Module Procedure Rtan_dp
       Module Procedure Rtan_sp
    End Interface

    Interface  Determinant
       Module Procedure Determinant_c
       Module Procedure Determinant_r
    End Interface

    Interface  Diagonalize_SH
       Module Procedure Diagonalize_HERM
       Module Procedure Diagonalize_SYMM
    End Interface

    Interface  Linear_Dependent
       Module Procedure Linear_Dependentc
       Module Procedure Linear_Dependenti
       Module Procedure Linear_Dependentr
    End Interface

    Interface  Rank
       Module Procedure Rank_dp
       Module Procedure Rank_sp
    End Interface

    Interface  Sort
       Module Procedure Sort_I
       Module Procedure Sort_R
    End Interface

    Interface  Svdcmp
       Module Procedure Svdcmp_dp
       Module Procedure Svdcmp_sp
    End Interface

    Interface Swap
        Module Procedure swap_c
        Module Procedure swap_cm
        Module Procedure swap_cv
        Module Procedure swap_i
        Module Procedure swap_im
        Module Procedure swap_iv
        Module Procedure swap_r
        Module Procedure swap_rm
        Module Procedure swap_rv
        Module Procedure masked_swap_r
        Module Procedure masked_swap_rm
        Module Procedure masked_swap_rv
    End interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Elemental Function Acosd(x) Result(arc_cos)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_cos
    !!----
    !!----    Inverse cosine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Acosd_dp(x) Result(arc_cos)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_dp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_cos

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp)  then
             arc_cos=0.0_dp
          else
             arc_cos=180.0_dp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_dp

    !!--++
    !!--++ Elemental Function Acosd_sp(x) Result(arc_cos)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_sp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_cos

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp)  then
             arc_cos=0.0_sp
          else
             arc_cos=180.0_sp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_sp

    !!----
    !!---- Function Asind(x) Result(arc_sin)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_sin
    !!----
    !!----    Inverse sine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Asind_dp(x) result(arc_sin)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_dp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_sin

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp) then
             arc_sin=90.0_dp
          else
             arc_sin=-90.0_dp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_dp

    !!--++
    !!--++ Elemental Function Asind_sp(x) result(arc_sin)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_sp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_sin

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp) then
             arc_sin=90.0_sp
          else
             arc_sin=-90.0_sp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_sp

    !!----
    !!---- Elemental Function Atan2d(y,x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: y,x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function of y/x
    !!----    y,x have the same units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atan2d_dp(y,x) Result(atande)
    !!--++    real(kind=dp), intent(in) :: y,x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_dp(y,x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: y,x
       real(kind=dp)             :: atand

       atand=atan2(y,x)*to_DEG

       return
    End Function Atan2d_dp

    !!--++
    !!--++ Elemental Function Atan2d_sp(y,x) Result(atande)
    !!--++    real(kind=sp), intent(in) :: y,x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_sp(y,x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: y,x
       real(kind=sp)             :: atande

       atande=atan2(y,x)*to_DEG

       return
    End Function Atan2d_sp

    !!----
    !!---- Elemental Function Atand(x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function, X no units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atand_dp(x) result(atande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_dp(x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: atand

       atand=atan(x)*to_DEG

       return
    End Function Atand_dp

    !!--++
    !!--++ Function Atand_sp(x) result(atande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_sp(x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: atande

       atande=atan(x)*to_DEG

       return
    End Function Atand_sp

    !!----
    !!---- Elemental Function Cosd(x) Result(cosine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: cosine
    !!----
    !!----    Cosine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Cosd_dp(x) Result(cosine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_dp(x) Result(cosine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_dp

    !!--++
    !!--++ Elemental Function Cosd_sp(x) Result(cosine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_sp(x) Result(cosine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_sp

    !!----
    !!---- Elemental Function Sind(x) Result(sine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: sine
    !!----
    !!----    Sine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Sind_dp(x) Result(sine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_dp(x) Result(sine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_dp

    !!--++
    !!--++ Elemental Function Sind_sp(x) Result(sine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_sp(x) Result(sine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_sp

    !!----
    !!---- Elemental Function Tand(x) Result(tande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: tande
    !!----
    !!----    Tangent function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Tand_dp(x) Result(tande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_dp(x) Result(tand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: tand

       tand=tan(to_RAD*x)

       return
    End Function Tand_dp

    !!--++
    !!--++ Elemental Function Tand_sp(x) Result(tande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_sp(x) Result(tande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: tande

       tande=tan(to_RAD*x)

       return
    End Function Tand_sp

    !!----
    !!---- Elemental Function Bessj0(x) Result(bessj_0)
    !!----    real(kind=sp), intent(in) :: x
    !!----    real(kind=sp)             :: bessj_0
    !!----
    !!----    Bessel Fuction J0(x)
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function Bessj0(x) Result(bessj_0)
       !---- Arguments ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: bessj_0

       !---- Local variables ----!
       real(kind=dp), parameter :: p1=   1.0_dp
       real(kind=dp), parameter :: p2=  -0.1098628627e-2_dp
       real(kind=dp), parameter :: p3=   0.2734510407e-4_dp
       real(kind=dp), parameter :: p4=  -0.2073370639e-5_dp
       real(kind=dp), parameter :: p5=   0.2093887211e-6_dp
       real(kind=dp), parameter :: q1=  -0.1562499995e-1_dp
       real(kind=dp), parameter :: q2=   0.1430488765e-3_dp
       real(kind=dp), parameter :: q3=  -0.6911147651e-5_dp
       real(kind=dp), parameter :: q4=   0.7621095161e-6_dp
       real(kind=dp), parameter :: q5=  -0.934945152e-7
       real(kind=dp), parameter :: r1=   57568490574.0_dp
       real(kind=dp), parameter :: r2=   57568490574.0_dp
       real(kind=dp), parameter :: r3=     651619640.7_dp
       real(kind=dp), parameter :: r4=     -11214424.18_dp
       real(kind=dp), parameter :: r5=         77392.33017_dp
       real(kind=dp), parameter :: r6=          -184.9052456_dp
       real(kind=dp), parameter :: s1=   57568490411.0_dp
       real(kind=dp), parameter :: s2=    1029532985.0_dp
       real(kind=dp), parameter :: s3=       9494680.718_dp
       real(kind=dp), parameter :: s4=         59272.64853_dp
       real(kind=dp), parameter :: s5=           267.8532712_dp
       real(kind=dp), parameter :: s6=             1.0_dp

       real(kind=dp)            :: y
       real(kind=sp)            :: ax, xx, z

       if (abs(x) < 1.0e-05) then
          bessj_0=1.0
          return
       end if
       if (abs(x) < 8.0)then
          y=x**2
          bessj_0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*  &
                  (s4+y*(s5+y*s6)))))
       else
          ax=abs(x)
          z=8.0/ax
          y=z**2
          xx=ax-0.785398164
          bessj_0=sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*  &
                  p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
       end if

       return
    End Function Bessj0

    !!----
    !!---- Elemental Function Bessj1(x) Result(bessj_1)
    !!----    real(kind=sp), intent(in) : x
    !!----    real(kind=sp)             : bessj_1
    !!----
    !!----    Bessel Fuction J1(x)
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function Bessj1(x) Result(bessj_1)
      !---- Arguments ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: bessj_1

       !---- Local variales ----!
       real(kind=dp), parameter :: p1= 1.0_dp
       real(kind=dp), parameter :: p2=  0.183105e-2_dp
       real(kind=dp), parameter :: p3= -0.3516396496e-4_dp
       real(kind=dp), parameter :: p4=  0.2457520174e-5_dp
       real(kind=dp), parameter :: p5= -0.240337019e-6_dp
       real(kind=dp), parameter :: q1=  0.04687499995_dp
       real(kind=dp), parameter :: q2= -0.2002690873e-3_dp
       real(kind=dp), parameter :: q3=  0.8449199096e-5_dp
       real(kind=dp), parameter :: q4= -0.88228987e-6_dp
       real(kind=dp), parameter :: q5=  0.105787412e-6_dp
       real(kind=dp), parameter :: r1=  72362614232.0_dp
       real(kind=dp), parameter :: r2=  -7895059235.0_dp
       real(kind=dp), parameter :: r3=    242396853.1_dp
       real(kind=dp), parameter :: r4=     -2972611.439_dp
       real(kind=dp), parameter :: r5=        15704.48260_dp
       real(kind=dp), parameter :: r6=          -30.16036606_dp
       real(kind=dp), parameter :: s1= 144725228442.0_dp
       real(kind=dp), parameter :: s2=   2300535178.0_dp
       real(kind=dp), parameter :: s3=     18583304.74_dp
       real(kind=dp), parameter :: s4=        99447.43394_dp
       real(kind=dp), parameter :: s5=          376.9991397_dp
       real(kind=dp), parameter :: s6=            1.0_dp

       real(kind=dp)            :: y
       real(kind=sp)            :: ax,xx,z

       if (abs(x) < 1.0e-05) then
          bessj_1=0.0
          return
       end if
       if (abs(x) < 8.0)then
          y=x**2
          bessj_1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+  &
                  y*(s4+y*(s5+y*s6)))))
       else
          ax=abs(x)
          z=8.0/ax
          y=z**2
          xx=ax-2.356194491
          bessj_1=sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))  &
                       -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.0,x)
       end if

       return
    End Function Bessj1

    !!----
    !!---- Elemental Function Factorial(n) Result(fact)
    !!----    integer, intent(in) : n
    !!----
    !!----    Factorial of N
    !!----
    !!---- Update: February - 2005
    !!
    Elemental Function Factorial(n) Result(fact)
       !---- Argument ----!
       integer, intent(in) :: n
       integer             :: fact

       !---- Local variables ----!
       integer   :: nt, np

       if (n ==0) then
          fact=1
       else
          nt=1
          np=abs(n)
          do
             nt=nt*np
             np=np-1
             if(np == 1) exit
          end do
          fact=nt
       end if

       return
    End Function Factorial

    !!----
    !!---- Elemental Function Negligible(v)
    !!----    complex/real(kind=sp),    intent( in) :: v
    !!----
    !!----    Provides the value .TRUE. if the real(kind=sp) (or complex)
    !!----    number V is less than EPS
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Negligiblec(v)
    !!--++    complex, intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate if a complex number is negligible
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligiblec(v) Result(Neglig)
       !---- Argument ----!
       complex, intent( in) :: v
       logical              :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligiblec

    !!--++
    !!--++ Elemental Function Negligibler(v)
    !!--++    real(kind=sp), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is negligible (abs < EPSS)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Negligibler(v) Result(neglig)
       !---- Argument ----!
       real(kind=sp), intent( in) :: v
       logical                    :: Neglig

       Neglig=.false.
       if (abs(v) > epss) return
       Neglig=.true.

       return
    End Function Negligibler

    !!----
    !!---- Function Pgcd(a,b) Result(mcd)
    !!----    integer, intent(in) :: a
    !!----    integer, intent(in) :: b
    !!----    integer             :: mcd
    !!----
    !!----    Function calculating the maximum common divisor of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Pgcd(a,b) Result(mcd)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcd

       !---- Local variables ----!
       integer  :: u,v,m

       u=max(a,b)
       v=min(a,b)
       m=0
       do
          if (m == 1) exit
          m=mod(u,v)
          u=v
          v=m
       end do
       mcd=u

       return
    End Function Pgcd

    !!----
    !!---- Function Ppcm(a,b) result(mcm)
    !!----    integer, intent(in) :: a
    !!----    integer, intent(in) :: b
    !!----    integer             :: mcm
    !!----
    !!----    Function calculating the minimum common multiple of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Function Ppcm(a,b) result(mcm)
       !---- Arguments ----!
       integer, intent(in) :: a,b
       integer             :: mcm

       !---- Local variables ----!
       integer :: u,v,w,i

       u=max(a,b)
       v=min(a,b)
       mcm=1
       if (v <= 1) then
          mcm=u
          return
       end if
       w=int(sqrt(real(u)))+1
       do i=2,w
          do
             if(.not. ((mod(u,i)==0) .or. (mod(v,i)==0)) ) exit
             mcm=mcm*i
             if (modulo(u,i) == 0) u=u/i
             if (modulo(v,i) == 0) v=v/i
          end do
       end do

       return
    End Function Ppcm

    !!----
    !!---- Function Pythag(a,b) Result (c)
    !!----    real(sp/dp),intent(in):: a,b
    !!----    real(sp/dp)           :: c
    !!--<<
    !!----    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!----    Adapted from Numerical Recipes.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Pythag_dp(a,b) Result (c)
    !!--++    real(dp),intent(in):: a,b
    !!--++    real(dp)           :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_dp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=dp),intent(in):: a,b
       real(kind=dp)           :: c

       !---- Local variables ----!
       real(kind=dp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa >absb)then
          c=absa*sqrt(1.0_dp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_dp))then
             c=0.0
          else
             c=absb*sqrt(1.0_dp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_dp

    !!--++
    !!--++ Function Pythag_sp(a,b) result (c)
    !!--++    real(sp),intent(in):: a,b
    !!--++    real(sp)           :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes c=sqrt(a^2 +b^2 ) without destructive underflow or overflow.
    !!--++    Adapted from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Pythag_sp(a,b) Result (c)
       !---- Arguments ----!
       real(kind=sp),intent(in):: a,b
       real(kind=sp)           :: c

       !---- Local variables ----!
       real(kind=sp)           :: absa,absb

       absa=abs(a)
       absb=abs(b)
       if (absa > absb) then
          c=absa*sqrt(1.0_sp+(absb/absa)**2)
       else
          if (absb < tiny(1.0_sp)) then
             c=0.0
          else
             c=absb*sqrt(1.0_sp+(absa/absb)**2)
          end if
       end if

       return
    End Function Pythag_sp

    !!----
    !!---- Logical Function Co_Linear(A,B,N)
    !!----    complex/integer/real(kind=sp), dimension(:), intent(in)  :: a
    !!----    complex/integer/real(kind=sp), dimension(:), intent(in)  :: b
    !!----    integer,                                     intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the vectors A and B are co-linear
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Co_Linear_C(A, B, N)
    !!--++    complex, dimension(:), intent(in)  :: a
    !!--++    complex, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two complex vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_C(a,b,n) Result(co_linear)
       !---- Argument ----!
       complex, dimension(:), intent(in) :: a,b
       integer,               intent(in) :: n
       logical                           :: co_linear

       !---- Local variables ----!
       integer :: i,ia,ib
       complex :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_C

    !!--++
    !!--++ Logical Function Co_Linear_I(A, B, N)
    !!--++    integer, dimension(:), intent(in)  :: a
    !!--++    integer, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_I(a,b,n) Result(co_linear)
       !---- Argument ----!
       integer, dimension(:), intent(in) :: a,b
       integer,               intent(in) :: n
       logical                           :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=sp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > 0) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > 0) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=real(a(ia))/real(b(ib))
          do i=1,n
             if (abs( a(i)-nint(c*real(b(i))) ) > 0) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_I

    !!--++
    !!--++ Logical Function Co_Linear_R(A, B, N)
    !!--++    real(kind=sp), dimension(:), intent(in)  :: a
    !!--++    real(kind=sp), dimension(:), intent(in)  :: b
    !!--++    integer,                     intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Co_linear_R(a,b,n) Result(co_linear)
       !---- Argument ----!
       real(kind=sp), dimension(:), intent(in) :: a,b
       integer,                     intent(in) :: n
       logical                                 :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib
       real(kind=sp) :: c

       co_linear=.true.
       do i=1,n
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do
       do i=1,n
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,n
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_R

    !!----
    !!---- Logical Function Equal_Matrix(A,B,N)
    !!----    integer/real(kind=sp), dimension(:,:), intent(in)  :: a,b
    !!----    integer,                               intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the array A is equal to array B
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Equal_Matrix_I(A, B, N)
    !!--++    integer, dimension(:,:), intent(in)  :: a
    !!--++    integer, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a,b
       integer                , intent(in) :: n
       logical                             :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (a(i,j) /= b(i,j)) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_I

    !!--++
    !!--++ Logical Function Equal_Matrix_R(A, B, N)
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: a
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: b
    !!--++    integer,                       intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Matrix_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=sp), dimension(:,:)   , intent(in) :: a,b
       integer,                          intent(in) :: n
       logical                                      :: info

       !---- Local variables ----!
       integer :: i,j

       info=.false.
       do i=1,n
          do j=1,n
             if (abs(a(i,j) - b(i,j)) > epss ) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_R

    !!----
    !!---- Logical Function Equal_Vector(A,B,N)
    !!----    integer/real(kind=sp), dimension(:),   intent(in)  :: a,b
    !!----    integer,                               intent(in)  :: n
    !!----
    !!----    Provides the value .TRUE. if the vector A is equal to vector B
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function Equal_Vector_I(A, B, N)
    !!--++    integer, dimension(:), intent(in)  :: a
    !!--++    integer, dimension(:), intent(in)  :: b
    !!--++    integer,               intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_I(a,b,n) result(info)
       !---- Argument ----!
       integer, dimension(:),   intent(in) :: a,b
       integer                , intent(in) :: n
       logical                             :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (a(i) /= b(i)) return
       end do
       info=.true.

       return
    End Function Equal_Vector_I

    !!--++
    !!--++ Logical Function Equal_Vector_R(A, B, N)
    !!--++    real(kind=sp), dimension(:), intent(in)  :: a
    !!--++    real(kind=sp), dimension(:), intent(in)  :: b
    !!--++    integer,                     intent(in)  :: n
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Equal_Vector_R(a,b,n) result(info)
       !---- Argument ----!
       real(kind=sp), dimension(:)   ,   intent(in) :: a,b
       integer,                          intent(in) :: n
       logical                                      :: info

       !---- Local variables ----!
       integer :: i

       info=.false.
       do i=1,n
          if (abs(a(i) - b(i)) > epss ) return
       end do
       info=.true.

       return
    End Function Equal_Vector_R

    !!----
    !!---- Function Imaxloc(arr) Result(mav)
    !!----  real(kind=sp)/integer, dimension(:), intent(in) :: arr
    !!----  integer                                         :: mav
    !!----
    !!----   Index of maxloc on an array
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Imaxloc_I(arr) Result(mav)
    !!--++  integer, dimension(:), intent(in) :: arr
    !!--++  integer                           :: mav
    !!--++
    !!--++   Index of maxloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Imaxloc_I(iarr) Result(mav)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: iarr
       integer                           :: mav

       !---- Local variables ----!
       integer, dimension(1) :: imax

       imax=maxloc(iarr(:))
       mav=imax(1)

       return
    End Function Imaxloc_I

    !!--++
    !!--++ Function Imaxloc_R(arr) Result(mav)
    !!--++  real(kind=sp), dimension(:), intent(in) :: arr
    !!--++  integer                                 :: mav
    !!--++
    !!--++   Index of maxloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Imaxloc_R(arr) Result(mav)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in) :: arr
       integer                                 :: mav

       !---- Local variables ----!
       integer, dimension(1) :: imax

       imax=maxloc(arr(:))
       mav=imax(1)

       return
    End Function Imaxloc_R

    !!----
    !!---- Function Iminloc(arr)  Result(miv)
    !!----  real(kind=sp)/integer, dimension(:), intent(in) :: arr
    !!----  integer                                         :: miv
    !!----
    !!----   Index of minloc on an array  (from Numerical Recipes)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Iminloc_I(arr)  Result(miv)
    !!--++  integer, dimension(:), intent(in) :: arr
    !!--++  integer                           :: miv
    !!--++
    !!--++   Index of minloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Iminloc_I(arr)  Result(miv)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: arr
       integer                           :: miv

       !---- Local variables ----!
       integer, dimension(1) :: imin

       imin=minloc(arr(:))
       miv=imin(1)

       return
    End Function Iminloc_I

    !!--++
    !!--++ Function Iminloc_R(arr)  Result(miv)
    !!--++  real(kind=sp), dimension(:), intent(in) :: arr
    !!--++  integer                                 :: miv
    !!--++
    !!--++   Index of minloc on an array (from Numerical Recipes)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Iminloc_R(arr)  Result(miv)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in) :: arr
       integer                                 :: miv

       !---- Local variables ----!
       integer, dimension(1) :: imin

       imin=minloc(arr(:))
       miv=imin(1)

       return
    End Function Iminloc_R

    !!----
    !!---- Function Locate(xx, n, x) Result(j)
    !!----    integer/real(kind=sp), dimension(n),intent(in)  :: xx
    !!----    integer ,                           intent(in)  :: n
    !!----    integer/real(kind=sp),              intent(in)  :: x
    !!----    integer ,                           intent(out) :: j
    !!----
    !!----    Subroutine for locating the index J of an array XX(N)
    !!----    satisfying:
    !!--<<
    !!----               XX(J) <= X < XX(J+1)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Locate_I(xx, n, x) Result(j)
    !!--++    integer, dimension(:),intent(in)  :: xx
    !!--++    integer ,             intent(in)  :: n
    !!--++    integer,              intent(in)  :: x
    !!--++    integer ,             intent(out) :: j
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Locate_I(xx,n,x) Result(j)
       !---- Argument ----!
       integer, dimension(:), intent(in):: xx
       integer ,              intent(in):: n
       integer,               intent(in):: x
       integer                          :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_I

    !!--++
    !!--++ Function Locate_R(xx, n, x) Result(j)
    !!--++    real(kind=sp), dimension(:),intent(in)  :: xx
    !!--++    integer ,                   intent(in)  :: n
    !!--++    real(kind=sp),              intent(in)  :: x
    !!--++    integer ,                   intent(out) :: j
    !!--++
    !!--++    Subroutine for locating the index J of an array XX(N)
    !!--++    satisfying:
    !!--++
    !!--++               XX(J) <= X < XX(J+1)
    !!--++
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Locate_R(xx,n,x) Result(j)
       !---- Argument ----!
       real(kind=sp), dimension(:), intent(in):: xx
       integer ,                    intent(in):: n
       real(kind=sp),               intent(in):: x
       integer                                :: j

       !---- Local Variables ----!
       integer :: jl, ju, jm

       jl=0
       ju=n+1
       do
          if(ju-jl <= 1) exit
          jm=(ju+jl)/2
          if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
             jl=jm
          else
             ju=jm
          end if
       end do
       j=jl

       return
    End Function Locate_R

    !!----
    !!---- Function Modulo_Lat(U)
    !!----    real(kind=sp), dimension(:), intent(in) :: u
    !!----
    !!----    Reduces a real vector to another with components in
    !!----    the interval [0,1)
    !!----
    !!---- Update: February - 2005
    !!
    Function Modulo_Lat(u) result(v)
       !---- Argument ----!
       real(kind=sp), dimension(:), intent( in) :: u
       real(kind=sp), dimension(1:size(u))      :: v

       v=mod(u+10.0_sp,1.0_sp)

       return
    End Function  Modulo_Lat

    !!----
    !!---- Function Outerprod(a,b) Result(c)
    !!----    real(sp/dp),dimension(:),intent(in)    :: a,b
    !!----    real(sp/dp),dimension(size(a),size(b)) :: c
    !!----
    !!----    Computes the outer product (tensorial product) of two
    !!----    vectors to give a tensor (matrix) as the result:
    !!--<<
    !!----                   c(i,j) = a(i)*b(j).
    !!-->>
    !!----    It uses the intrinsic Fortran 90 function SPREAD.
    !!----    Function adapted from Numerical Recipes.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Outerprod_dp(a,b) Result(c)
    !!--++    real(dp),dimension(:),intent(in)    :: a,b
    !!--++    real(dp),dimension(size(a),size(b)) :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_dp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=dp),dimension(:),intent(in)    :: a,b
       real(kind=dp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_dp

    !!--++
    !!--++ Function Outerprod_sp(a,b) Result(c)
    !!--++    real(sp),dimension(:),intent(in)    :: a,b
    !!--++    real(sp),dimension(size(a),size(b)) :: c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the outer product (tensorial product) of two
    !!--++    vectors to give a tensor (matrix) as the result:
    !!--++                   c(i,j) = a(i)*b(j).
    !!--++
    !!--++    It uses the intrinsic Fortran 90 function SPREAD.
    !!--++    Taken from Numerical Recipes.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Outerprod_sp(a,b)  Result(c)
       !---- Arguments ----!
       real(kind=sp),dimension(:),intent(in)    :: a,b
       real(kind=sp),dimension(size(a),size(b)) :: c

       c =spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

       return
    End Function Outerprod_sp

    !!----
    !!---- Function Traza(A)
    !!----    complex/integer/real(kind=sp), dimension(:,:), intent(in)  :: a
    !!----
    !!----    Provides the trace of a complex/real or integer matrix
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Traza_C(A)
    !!--++    complex, dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a complex nxn array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Traza_C(a) Result(b)
       !---- Argument ----!
       complex, dimension(:,:), intent(in) :: a
       complex                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=(0.0,0.0)
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Traza_C

    !!--++
    !!--++ Function Traza_I(A)
    !!--++    integer, dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of an integer 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Traza_I(a) Result(b)
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a
       integer                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Traza_I

    !!--++
    !!--++ Function Traza_R(A)
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a real(kind=sp) 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Traza_R(a) Result(b)
       !---- Argument ----!
       real(kind=sp), dimension(:,:), intent(in) :: a
       real(kind=sp)                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0.0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Traza_R

    !!----
    !!---- Logical Function Zbelong(V)
    !!----    real(kind=sp),   dimension(:,:), intent( in) :: v
    !!----                      or
    !!----    real(kind=sp),   dimension(:),   intent( in) :: v
    !!----                      or
    !!----    real(kind=sp),                   intent( in) :: v
    !!----
    !!----    Provides the value .TRUE. if the real(kind=sp) number (or array) V is close enough
    !!----    (whithin EPSS) to an integer.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Logical Function ZbelongM(V)
    !!--++    real(kind=sp),   dimension(:,:), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real(kind=sp) array is an Integer matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongM(v) Result(belong)
       !---- Argument ----!
       real(kind=sp),   dimension(:,:), intent( in) :: v
       logical                                      :: belong

       !---- Local variables ----!
       real(kind=sp),   dimension(size(v,1),size(v,2)) :: vec

       vec= abs(real(nint (v))-v)
       belong=.not. ANY(vec > epss)

       return
    End Function ZbelongM

    !!--++
    !!--++ Logical Function ZbelongN(A)
    !!--++    real(kind=sp),  intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real(kind=sp) number is an Integer
    !!--++
    !!--++ Update: February - 2005
   !!
    Function ZbelongN(a) Result(belong)
       !---- Argument ----!
       real(kind=sp), intent( in) :: a
       logical                    :: belong

       belong=.false.
       if (abs(real(nint (a))-a) > epss) return
       belong=.true.

       return
    End Function ZbelongN

    !!--++
    !!--++ Logical Function ZbelongV(V)
    !!--++    real(kind=sp),   dimension(:), intent( in) :: v
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real(kind=sp) vector is an Integer vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Function ZbelongV(v) Result(belong)
       !---- Argument ----!
       real(kind=sp),   dimension(:), intent( in) :: v
       logical                                    :: belong

       !---- Local variables ----!
       integer                             :: i
       real(kind=sp),   dimension(size(v)) :: vec

       belong=.false.
       vec= abs(real(nint (v))-v)
       do i=1,size(v)
          if (vec(i) > epss) return
       end do
       belong=.true.

       return
    End Function ZbelongV

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Mathgen()
    !!----
    !!----    Initialize the errors flags in Math_Gen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_MathGen()

       err_math_gen=.false.
       err_mess_math_gen=" "

       return
    End Subroutine Init_Err_MathGen

    !!----
    !!---- Subroutine Set_Epsg(Neweps)
    !!----    real(kind=sp), intent( in) :: neweps
    !!----
    !!----    Sets global EPSS to the value "neweps"
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg(Neweps)
       !---- Arguments ----!
       real(kind=sp), intent( in) :: neweps

       epss=neweps

       return
    End Subroutine Set_Epsg

    !!----
    !!---- Subroutine Set_Epsg_Default()
    !!----
    !!----    Sets global EPSS to the default value: epss=1.0E-5_sp
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Set_Epsg_Default()

       epss=1.0E-5_sp

       return
    End Subroutine Set_Epsg_Default

    !!----
    !!---- Subroutine Rtan(x,y,ang,deg)
    !!----    real(sp/dp),               intent( in) :: x,y
    !!----    real(sp/dp),               intent(out) :: ang
    !!----    character(len=*),optional, intent( in) :: deg
    !!----
    !!----    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!----    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Rtan_dp(x,y,ang,deg)
    !!--++    real(dp),                  intent( in) :: x,y
    !!--++    real(dp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_dp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=dp),              Intent( In)   :: x,y
       real(kind=dp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- Local variables ----!
       real(kind=dp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_dp
          return
       else if(abx < eps) then
          ang = pi/2.0_dp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_dp) ang = pi-ang
       else
          ang = pi/2.0_dp - atan(abx/aby)
          if(x < 0.0_dp) ang = pi-ang
       end if
       if (y < 0.0_dp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_dp

    !!--++
    !!--++ Subroutine Rtan_sp(x,y,ang,deg)
    !!--++    real(sp),                  intent( in) :: x,y
    !!--++    real(sp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_sp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=sp),              Intent( In)   :: x,y
       real(kind=sp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- local variables ----!
       real(kind=sp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < eps) .and. (aby < eps)) then
          ang = 0.0_sp
          return
       else if(abx < eps) then
          ang = pi/2.0_sp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_sp) ang = pi-ang
       else
          ang = pi/2.0_sp - atan(abx/aby)
          if(x < 0.0_sp) ang = pi-ang
       end if
       if(y < 0.0_sp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_sp

    !!----
    !!---- Subroutine Determinant(A,n,determ)
    !!----    complex/real(sp), dimension(:,:), intent( in) :: A      !input square matrix (n,n)
    !!----    integer,                          intent( in) :: n      !actual dimension of A
    !!----    real(kind=sp),                    intent(out) :: determ !det(A) if real
    !!----                                                             det(AR)^2 + det(AI)^2 if complex
    !!----
    !!----    Calculates the determinant of a real(kind=sp) square matrix.
    !!----    Calculates the pseudo-determinant of a complex square matrix.
    !!----    The calculated value is only useful for linear dependency purposes.
    !!----    It tell us if the complex matrix is singular or not.
    !!--..
    !!--..    Calculates the determinant of a complex square matrix selected from a rectangular
    !!--..    matrix A, n x m, where m >= n. determ=determinant_of_A(1:n,icol:icol+n-1)
    !!--..    If icol is absent, the calculation is performed as if icol=1.
    !!--..    If icol+n-1 > m, or m < n, determ is set to 0.0 and an error message is generated.
    !!----
    !!--..    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Determinant_C(A,n,determ)
    !!--++    complex,          dimension(:,:), intent( in) :: A      !input square matrix (n,n)
    !!--++    integer,                          intent( in) :: n      !actual dimension of A
    !!--++    real(kind=sp),                    intent(out) :: determ !det(A) if real
    !!--++                                                             det(AR)^2 + det(AI)^2 if complex
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real(kind=sp) square matrix.
    !!--++    Calculates the pseudo-determinant of a complex square matrix.
    !!--++    The calculated value is only useful for linear dependency purposes.
    !!--++    It tell us if the complex matrix is singular or not.
    !!--++
    !!--++    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Determinant_C(A,n,determ)
       !---- Arguments ----!
       complex, dimension(:,:), intent( in) :: A
       integer,                 intent( in) :: n
       real(kind=sp),           intent(out) :: determ

       !---- local variables ----!
       real(kind=sp),    dimension(2*n,2*n) :: AC   !real(kind=sp) square matrix
       real(kind=sp)                        :: d
       integer                              :: i,nn
       logical                              :: singular

       nn=2*n
       AC(  1:n ,  1:n ) =  real(A(1:n ,1:n))
       AC(n+1:nn,  1:n ) = aimag(A(1:n ,1:n))
       AC(n+1:nn,n+1:nn) =    AC(  1:n ,1:n)
       AC(  1:n ,n+1:nn) =   -AC(n+1:nn,1:n)

       call lu_decomp(ac(1:nn,1:nn),d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,nn
             d=d*sign(1.0,ac(i,i))
             determ=determ+ log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_C

    !!--++
    !!--++ Subroutine Determinant_R(A,n,determ)
    !!--++    real(kind=sp), dimension(:,:),intent( in) :: A   (input square matrix (n,n))
    !!--++    integer,                      intent( in) :: n   (actual dimension of A)
    !!--++    real(kind=sp),                intent(out) :: determ  (determinant )
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real(kind=sp) square matrix.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Determinant_R(A,n,determ)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent( in) :: A
       integer,                       intent( in) :: n
       real(kind=sp),                 intent(out) :: determ

       !---- local variables ----!
       real(kind=sp),    dimension(n,n)  :: AC
       real(kind=sp)                     :: d
       integer                           :: i
       logical                           :: singular

       ac=A(1:n,1:n)
       call lu_decomp(ac,d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,n
             d=d*sign(1.0,ac(i,i))
             determ=determ + log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_R

    !!----
    !!---- Subroutine Diagonalize_SH(A,n,e_val,e_vect)
    !!----    complex/real,      dimension(:,:), intent( in)  :: A
    !!----    integer,                           intent( in)  :: n
    !!----    real(kind=sp),     dimension(:),   intent(out)  :: E_val
    !!----    complex, optional, dimension(:,:), intent(out)  :: E_vect
    !!----
    !!----    Diagonalize Symmetric/Hermitian matrices.
    !!----    The eigen_values E_val are sorted in descending order. The columns
    !!----    of E_vect are the corresponding eigenvectors.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Diagonalize_Herm(a,n,e_val,e_vect)
    !!--++    complex,           dimension(:,:), intent( in)  :: A
    !!--++    integer,                           intent( in)  :: n
    !!--++    real(kind=sp),     dimension(:),   intent(out)  :: E_val
    !!--++    complex, optional, dimension(:,:), intent(out)  :: E_vect
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize Hermitian matrices.
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Herm(a,n,e_val,e_vect)
       !---- Arguments ----!
       complex,           dimension(:,:), intent( in)  :: A
       integer,                           intent( in)  :: n
       real(kind=sp),     dimension(:),   intent(out)  :: E_val
       complex, optional, dimension(:,:), intent(out)  :: E_vect

       !---- Local variables ----!
       real(kind=sp),        dimension(2*n,2*n)   :: aux
       real(kind=sp),        dimension(2*n)       :: e,d
       integer :: nn

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          err_math_gen=.true.
          err_mess_math_gen=" Diagonalize_HERM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       nn=2*n
       aux(  1:n ,  1:n ) =  real(a(1:n ,1:n))   !      (  U   V )
       aux(n+1:nn,n+1:nn) =  real(a(1:n ,1:n))   !   M=(          ),   A = U + i V
       aux(n+1:nn,  1:n ) = aimag(a(1:n ,1:n))   !      ( -V   U )
       aux(  1:n ,n+1:nn) =-aimag(a(1:n ,1:n))   !

       if (present(E_vect)) then
          call tred2(aux,nn,d,e)
          call tqli2(d,e,nn,aux)
          call eigsrt(d,aux,nn,1)
          e_vect(1:n,1:n)=cmplx(aux(1:n,1:nn:2),aux(n+1:nn,1:nn:2))
       else
          call tred1(aux,nn,d,e)
          call tqli1(d,e,nn)
          call eigsrt(d,aux,nn,0)
       end if
       e_val(1:n)=d(1:nn:2)

       return
    End Subroutine Diagonalize_Herm

    !!--++
    !!--++ Subroutine Diagonalize_Symm(a,n,e_val,e_vect)
    !!--++    real(kind=sp)            dimension(:,:),intent( in)  :: A      (input matrix with)
    !!--++    integer,                                intent( in)  :: n      (actual dimension)
    !!--++    real(kind=sp),           dimension(:),  intent(out)  :: E_val  (eigenvalues)
    !!--++    real(kind=sp), optional, dimension(:,:),intent(out)  :: E_vect (eigenvectors)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Diagonalize symmetric matrices
    !!--++    The eigen_values E_val are sorted in descending order. The columns
    !!--++    of E_vect are the corresponding eigenvectors.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Diagonalize_Symm(A,n,E_Val,E_vect)
       !---- Arguments ----!
       real(kind=sp),           dimension(:,:), intent( in)  :: A
       integer,                                 intent( in)  :: n
       real(kind=sp),           dimension(:),   intent(out)  :: E_val
       real(kind=sp), optional, dimension(:,:), intent(out)  :: E_vect

       !---- Local variables ----!
       real(kind=sp),        dimension(n,n)   :: aux
       real(kind=sp),        dimension(n)     :: e

       e_val=0.0
       call init_err_mathgen()
       if (n > size(A,1) .or. n > size(A,2)) then
          err_math_gen=.true.
          err_mess_math_gen=" Diagonalize_SYMM: Error in dimension of input matrix: A(m,m) with m < n "
          return
       end if

       aux=a(1:n,1:n)
       if (present(E_vect)) then
          call tred2(aux,n,E_val,e)
          call tqli2(E_val,e,n,aux)
          call eigsrt(E_val,aux,n,1)
          e_vect(1:n,1:n)=aux
       else
          call tred1(aux,n,E_val,e)
          call tqli1(E_val,e,n)
          call eigsrt(E_val,aux,n,0)
       end if

       return
    End Subroutine Diagonalize_Symm

    !!--++
    !!--++ Subroutine Eigsrt(d,v,n,io)
    !!--++    real(kind=sp), dimension(:),   intent(in out) :: d
    !!--++    real(kind=sp), dimension(:,:), intent(in out) :: v
    !!--++    integer,                       intent (in)    :: n
    !!--++    integer,                       intent (in)    :: io
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for sorting eigenvalues in d(n) and eigenvectors
    !!--++    in columns of v(n,n). Sorts d(n) in descending order and
    !!--++    rearranges v(n,n) correspondingly. The method is the straight
    !!--++    insertion. If io=0 order  only the eigenvalues are treated.
    !!--++    Adapted from Numerical Recipes. Valid for hermitian matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Eigsrt(d,v,n,io)
       !---- Arguments ----!
       real(kind=sp), dimension(:),   intent(in out) :: d
       real(kind=sp), dimension(:,:), intent(in out) :: v
       integer,                       intent(in)     :: n
       integer,                       intent(in)     :: io

       !---- Local Variables ----!
       integer          :: i,j,k
       real(kind=sp)    :: p

       do i=1,n-1
          k=i
          p=d(i)
          do j=i+1,n
             if (d(j) >= p) then
                k=j
                p=d(j)
             end if
          end do
          if (k /= i) then
             d(k)=d(i)
             d(i)=p
             if (io == 1) then
                do j=1,n
                   p=v(j,i)
                   v(j,i)=v(j,k)
                   v(j,k)=p
                end do
             end if
          end if
       end do

       return
    End Subroutine Eigsrt

    !!----
    !!---- Subroutine Invert_Matrix(a,b,singular,perm)
    !!----    real(kind=sp), dimension(:,:),  intent( in) :: a
    !!----    real(kind=sp), dimension(:,:),  intent(out) :: b
    !!----    LOGICAL,                        intent(out) :: singular
    !!----    integer, dimension(:),optional, intent(out) :: perm
    !!--<<
    !!----    Subroutine to invert a real matrix using LU decomposition.
    !!----    In case of singular matrix (singular=.true.) instead of the inverse
    !!----    matrix, the subroutine provides the LU decomposed matrix as used
    !!----    in Numerical Recipes.
    !!----    The input matrix is preserved and its inverse (or its LU decomposition)
    !!----    is provided in "b". The optional argument "perm" holds the row permutation
    !!----    performed to obtain the LU decomposition.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Invert_Matrix(a,b,singular,perm)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),  intent(in ) :: a
       real(kind=sp), dimension(:,:),  intent(out) :: b
       logical,                        intent(out) :: singular
       integer, dimension(:),optional, intent(out) :: perm

       !---- Local variables ----!
       integer                                       :: i,n
       integer,       dimension(size(a,1))           :: indx
       real(kind=sp)                                 :: d, det
       real(kind=sp), dimension(size(a,1),size(a,1)) :: lu

       n=size(a,1)
       lu=a(1:n,1:n)

       call LU_Decomp(lu,d,singular,indx)
       if (present(perm)) perm(1:n)=indx(1:n)

       if (singular) then
          b=lu
          return
       else
          det=0.0
          do i=1,n
             d=d*sign(1.0,lu(i,i))
             det=det + log(abs(lu(i,i)))
          end do
          det=d*exp(det)
          if (abs(det) <= 1.0e-30) then
             singular=.true.
             b=lu
             return
          end if
       end if

       b=0.0
       do i=1,n
          b(i,i)=1.0
          call LU_backsub(lu,indx,b(:,i))
       end do

       return
    End Subroutine Invert_Matrix

    !!----
    !!---- Subroutine Linear_Dependent(a,na,b,nb,mb,linear_dependent)
    !!----    complex/integer/real(kind=sp), dimension(:),   intent(in)  :: a
    !!----    complex/integer/real(kind=sp), dimension(:,:), intent(in)  :: b
    !!----    integer,                                       intent(in)  :: na,nb,mb
    !!----    logical,                                       intent(out) :: Linear_Dependent
    !!--<<
    !!----    Provides the value .TRUE. if the vector A is linear dependent of the
    !!----    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!----    are the number of rows and columns of B to be considered. The actual
    !!----    dimension of vector a should be na=max(nb,mb).
    !!----    The problem is equivalent to determine the rank (in algebraic sense)
    !!----    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!----    case it is supposed that na = mb and in the second na = nb.
    !!----    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!----    is generated. The function uses floating arithmetic for all types.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Linear_DependentC(a,na,b,nb,mb,linear_dependent)
    !!--++    complex, dimension(:),   intent(in)  :: a
    !!--++    complex, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: na,nb,mb
    !!--++    logical,                 intent(out) :: Linear_Dependent
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++
    !!--++    For the case of complex vectors in Cn the problem can be reduced to real vectors
    !!--++    of dimension R2n. Each complex vector contributes as two real vectors of dimension
    !!--++    2n: (R,I) and (-I,R). A complex vector V is linearly dependent on n complex vectors
    !!--++    if V can be written as: V = Sigma{j=1,n}(Cj.Vj), with Cj complex numbers and Vj
    !!--++    having n complex components. One may write:
    !!--++
    !!--++     V = Sigma{j=1,n}(Cj.Vj)
    !!--++     (R,I) = Sigma{j=1,n} (Cjr Vj + i Cji Vj) = Sigma{j=1,n} (Cjr (Rj,Ij) +  Cji (-Ij,Rj) )
    !!--++     (R,I) = Sigma{j=1,n} (aj (Rj,Ij) + bj (-Ij,Rj) )  = Sigma{j=1,2n} (Aj.Uj)
    !!--++     Were Uj=(Rj,Ij) and U(j+1)= (-Ij,Rj)
    !!--++
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentC(A,na,B,nb,mb,Linear_Dependent)
       !---- Arguments ----!
       complex, dimension(:),   intent(in)  :: a
       complex, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: Linear_Dependent

       !---- Local variables ----!
       integer                                                     :: r,n1
       real(kind=dp), parameter                                    :: tol= 100.0_dp*deps
       real(kind=dp), dimension(2*max(nb+1,mb+1),2*max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       Linear_Dependent=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentC: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=2*nb+1
          if(n1+1 > 2*mb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(1:nb,     mb+1:mb+na) = aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,      1:mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:2*nb,mb+1:mb+na) =  real(b(1:nb,1:mb))
          c(n1,             1:mb) =  real(a(1:na))
          c(n1,      mb+1:mb+na ) = aimag(a(1:na))
          c(n1+1,           1:mb) =-aimag(a(1:na))
          c(n1+1,    mb+1:mb+na ) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*mb)) Linear_Dependent=.false.
       else if( na == nb) then
          n1=2*mb+1
          if(n1+1 > 2*nb) return !the vector is linear dependent
          c(1:nb,           1:mb) =  real(b(1:nb,1:mb))
          c(nb+1:nb+na,     1:mb) = aimag(b(1:nb,1:mb))
          c(1:nb,      mb+1:2*mb) =-aimag(b(1:nb,1:mb))
          c(nb+1:nb+na,mb+1:2*mb) =  real(b(1:nb,1:mb))
          c(1:na,             n1) =  real(a(1:na))
          c(nb+1:nb+na,       n1) = aimag(a(1:na))
          c(1:na,           1+n1) =-aimag(a(1:na))
          c(nb+1:nb+na,     1+n1) =  real(a(1:na))
          call rank(c,tol,r)
          if(r == min(n1+1,2*nb)) Linear_Dependent=.false.
       else
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentC: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentC

    !!--++
    !!--++ Subroutine Linear_DependentI(a,na,b,nb,mb,linear_dependent)
    !!--++    integer, dimension(:),   intent(in)  :: a
    !!--++    integer, dimension(:,:), intent(in)  :: b
    !!--++    integer,                 intent(in)  :: na,nb,mb
    !!--++    logical,                 intent(out) :: Linear_Dependent
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentI(A,na,B,nb,mb,Linear_Dependent)
       !---- Arguments ----!
       integer, dimension(:),   intent(in)  :: a
       integer, dimension(:,:), intent(in)  :: b
       integer,                 intent(in)  :: na,nb,mb
       logical,                 intent(out) :: Linear_Dependent

       !---- Local variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       Linear_Dependent=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentI: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(n1,  1:mb)=real(a(1:na))      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) Linear_Dependent=.false.
       else if( na == nb) then
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=real(b(1:nb,1:mb))
          c(1:nb,  n1)=real(a(1:na))     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) Linear_Dependent=.false.
       else
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentI: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentI

    !!--++
    !!--++ Subroutine Linear_DependentR(a,na,b,nb,mb,linear_dependent)
    !!--++    real(kind=sp), dimension(:),   intent(in)  :: a
    !!--++    real(kind=sp), dimension(:,:), intent(in)  :: b
    !!--++    integer,                       intent(in)  :: na,nb,mb
    !!--++    logical,                                       intent(out) :: Linear_Dependent
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Provides the value .TRUE. if the vector A is linear dependent of the
    !!--++    vectors constituting the rows (columns) of the matrix B. In input nb & mb
    !!--++    are the number of rows and columns of B to be considered. The actual
    !!--++    dimension of vector a should be na=max(nb,mb).
    !!--++    The problem is equivalent to determine the rank (in algebraic sense)
    !!--++    of the composite matrix C(nb+1,mb)=(B/A) or C(nb,mb+1)=(B|A). In the first
    !!--++    case it is supposed that na = mb and in the second na = nb.
    !!--++    and the rank of B is min(nb, mb). If na /= nb and na /= mb an error condition
    !!--++    is generated
    !!--++    The function uses floating arithmetic for all types.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Linear_DependentR(A,na,B,nb,mb,Linear_Dependent)
       !---- Arguments ----!
       real(kind=sp), dimension(:),   intent(in)  :: a
       real(kind=sp), dimension(:,:), intent(in)  :: b
       integer,                       intent(in)  :: na,nb,mb
       logical,                       intent(out) :: Linear_Dependent

       !---- Local Variables ----!
       integer                                                 :: r,n1
       real(kind=dp), parameter                                :: tol= 100.0_dp*deps
       real(kind=dp), dimension(max(nb+1,mb+1),max(nb+1,mb+1)) :: c

       c=0.0
       call init_err_mathgen()
       Linear_Dependent=.true.
       if (nb > size(b,1) .or. mb > size(b,2) .or. na > size(a) ) then
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentR: Error in dimension of input matrix or vector"
          return
       end if

       if ( na == mb) then    !Vector added as an additional row
          n1=nb+1
          if(n1 > mb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(n1,  1:mb)=a(1:na)      !C(nb+1,mb)
          call rank(c,tol,r)
          if(r == min(n1,mb)) Linear_Dependent=.false.
       else if( na == nb) then   !Vector added as an additional column
          n1=mb+1
          if(n1 > nb) return !the vector is linear dependent
          c(1:nb,1:mb)=b(1:nb,1:mb)
          c(1:nb,  n1)=a(1:na)     !C(nb,mb+1)
          call rank(c,tol,r)
          if(r == min(n1,nb)) Linear_Dependent=.false.
       else
          err_math_gen=.true.
          err_mess_math_gen=" Linear_DependentR: input dimension of vector incompatible with matrix"
       end if

       return
    End Subroutine Linear_DependentR

    !!----
    !!---- Subroutine LU_Backsub(a,indx,b)
    !!----    real(kind=sp),    dimension(:,:),intent(in)     :: a
    !!----    integer,          dimension(:),  intent(in)     :: indx
    !!----    real(kind=sp),    dimension(:),  intent(in out) :: b
    !!--<<
    !!----    Adapted from Numerical Recipes.
    !!----    Solves the set of N linear equations A � X = B. Here the N � N matrix A is input,
    !!----    not as the original matrix A, but rather as its LU decomposition, determined
    !!----    by the routine LU_DECOMP. INDX is input as the permutation vector of length N
    !!----    returned by LU_DECOMP. B is input as the right-hand-side vector B,
    !!----    also of length N, and returns with the solution vector X.
    !!----    A and INDX are not modified by this routine and can be left in place for successive calls
    !!----    with different right-hand sides B. This routine takes into account the possibility that B will
    !!----    begin with many zero elements, so it is efficient for use in matrix inversion.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Backsub(a,indx,b)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in)     :: a
       integer,         dimension(:), intent(in)     :: indx
       real(kind=sp),   dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer       :: i,ii,ll,n
       real(kind=sp) :: summ

       n=size(a,1)
       ii=0              !When ii is set to a positive value, it will become the index
       do i=1,n          !of the first nonvanishing element of b. We now do
          ll=indx(i)     !the forward substitution. The only new wrinkle is to
          summ=b(ll)     !unscramble the permutation as we go.
          b(ll)=b(i)
          if (ii /= 0) then
             summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
          else if(summ /= 0.0) then   !A nonzero element was encountered, so from now on
             ii=i                       !we will have to do the dot product above.
          end if
          b(i)=summ
       end do

       do i=n,1,-1       !Now we do the backsubstitution
          b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do

       return
    End Subroutine LU_Backsub

    !!----
    !!---- Subroutine LU_Decomp(a,d,singular,indx)
    !!----    real(kind=sp),    dimension(:,:),intent(in out) :: a
    !!----    real(kind=sp),                   intent(out)    :: d
    !!----    logical,                         intent(out)    :: singular
    !!----    integer, dimension(:), optional, intent(out)    :: indx
    !!--<<
    !!----    Subroutine to make the LU decomposition of an input matrix A.
    !!----    The input matrix is destroyed and replaced by a matrix containing
    !!----    in its upper triangular part (plus diagonal) the matrix U. The
    !!----    lower triangular part contains the nontrivial part (Lii=1) of matrix L.
    !!----    The output is rowwise permutation of the initial matrix. The vector INDX
    !!----    recording the row permutation. D is output as +/-1 depending on whether
    !!----    the number of row interchanges was even or odd, respectively.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine LU_Decomp(a,d,singular,indx)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a
       real(kind=sp),                 intent(out)    :: d
       logical,                       intent(out)    :: singular
       integer,  dimension(:), intent(out), optional :: indx

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1)):: vv  !vv stores the implicit scaling of each row.
       real(kind=sp), parameter           :: vtiny = 1.0e-20_sp !A small number.
       integer                            :: j,imax,n

       singular=.false.
       n=size(a,1)
       d=1.0                      !No row interchanges yet.
       vv=maxval(abs(a),dim=2)    !Loop over rows to get the implicit scaling information.
       if (any(abs(vv) <= vtiny)) then   !There is a row of zeros.
          singular=.true.
          return
       end if
       vv=1.0_sp/vv     !Save the scaling.
       do j=1,n
          imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))   !Find the pivot row.
          if (j /= imax) then                         !Do we need to interchange rows?
             call swap(a(imax,:),a(j,:))              !Yes, do so...
             d=-d                                     !...and change the parity of d.
             vv(imax)=vv(j)                           !Also interchange the scale factor.
          end if
          if (present(indx)) indx(j)=imax
          if (abs(a(j,j)) <= vtiny) then !If the pivot element is zero the matrix is singular.
             a(j,j)=vtiny                !(at least to the precision of the algorithm)
             singular=.true.             !For some applications on singular matrices,
             return                      !it is desirable to substitute vtiny for zero.
          end if                         ! This is actually the present case
          a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                    !Divide by the pivot element.
          a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))  !Reduce remaining submatrix.
       end do

       return
    End Subroutine LU_Decomp

    !!----
    !!---- Subroutine Matinv(a,n)
    !!----    real(kind=sp), dimension(:,:),intent(in out) :: a
    !!----    integer     ,                 intent(in)     :: n
    !!----
    !!----  Subroutine for inverting a real square matrix.
    !!----  The input matrix is replaced in output with its inverse.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matinv(a,n)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a
       integer     ,                  intent(in)     :: n

       !---- Local variables ----!
       real(kind=sp)                 :: amax,savec
       integer, dimension(size(a,1)) :: ik,jk
       integer                       :: i,j,k,l

       !---- Subroutine to invert a real matrix ----!
       do k=1,n
          amax=0.0
          do
             do
                do i=k,n
                   do j=k,n
                      if (abs(amax)-abs(a(i,j)) > 0.0) cycle
                      amax=a(i,j)
                      ik(k)=i
                      jk(k)=j
                   end do
                end do
                i=ik(k)
                if (i-k < 0) cycle
                exit
             end do

             if (i-k /= 0) then
                do j=1,n
                   savec=a(k,j)
                   a(k,j)=a(i,j)
                   a(i,j)=-savec
                end do
             end if

             j=jk(k)
             if (j-k < 0) cycle
             exit
          end do

          if (j-k /= 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=a(i,j)
                a(i,j)=-savec
             end do
          end if

          do i=1,n
             if (i-k /= 0)  then
                a(i,k)=-a(i,k)/amax
             end if
          end do
          do i=1,n
             do j=1,n
                if (i-k == 0 .or. j-k == 0) cycle
                a(i,j)=a(i,j)+a(i,k)*a(k,j)
             end do
          end do
          do j=1,n
             if (j-k == 0)   cycle
             a(k,j)=a(k,j)/amax
          end do
          a(k,k)=1.0/amax
       end do     !k

       do l=1,n
          k=n-l+1
          j=ik(k)
          if (j-k > 0) then
             do i=1,n
                savec=a(i,k)
                a(i,k)=-a(i,j)
                a(i,j)=savec
             end do
          end if
          i=jk(k)
          if (i-k > 0) then
             do j=1,n
                savec=a(k,j)
                a(k,j)=-a(i,j)
                a(i,j)=savec
             end do
          end if
       end do

       return
    End Subroutine Matinv

    !!--++
    !!--++ Subroutine Partition(A, marker)
    !!--++    character(len=*), dimension(:), intent(in out) :: A
    !!--++    integer,                        intent(out)    :: marker
    !!--++
    !!--++    (Private)
    !!--++    Utilised by Sort_Strings.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Partition(A, Marker)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in out) :: A
       integer,                        intent(   out) :: marker

       !---- Local variables ----!
       integer                  :: i, j
       character(len=len(A(1))) :: temp
       character(len=len(A(1))) :: x      ! pivot point

       x = A(1)
       i= 0
       j= size(A) + 1

       do
          j = j-1
          do
             if (A(j) <= x) exit
             j = j-1
          end do
          i = i+1
          do
             if (A(i) >= x) exit
             i = i+1
          end do
          if (i < j) then
             !---- exchange A(i) and A(j)
             temp = A(i)
             A(i) = A(j)
             A(j) = temp
          else if (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          end if
       end do

       return
    End Subroutine Partition

    !!----
    !!---- Subroutine Rank(a,tol,r)
    !!----    real(sp/dp), dimension(:,:), intent( in) :: a
    !!----    real(sp/dp),                 intent( in) :: tol
    !!----    integer,                     intent(out) :: r
    !!----
    !!----    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Rank_dp(a,tol,r)
    !!--++    real(dp), dimension(:,:), intent( in) :: a
    !!--++    real(dp),                 intent( in) :: tol
    !!--++    integer,                  intent(out) :: r
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_dp(a,tol,r)
       !---- Arguments ----!
       real(kind=dp), dimension(:,:),intent( in)      :: a
       real(kind=dp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Arguments ----!
       real(kind=dp), dimension(size(a,1),size(a,2))  :: u
       real(kind=dp), dimension(size(a,2))            :: w
       real(kind=dp), dimension(size(a,2),size(a,2))  :: v
       integer                                        :: i

       u=a
       call svdcmp(u,w,v)
       if (err_math_gen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_dp

    !!--++
    !!--++ Subroutine Rank_sp(a,tol,r)
    !!--++    real(sp), dimension(:,:), intent( in) :: a
    !!--++    real(sp),                 intent( in) :: tol
    !!--++    integer,                  intent(out) :: r
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rank_sp(a,tol,r)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),intent( in)      :: a
       real(kind=sp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2))  :: u
       real(kind=sp), dimension(size(a,2))            :: w
       real(kind=sp), dimension(size(a,2),size(a,2))  :: v
       integer :: i

       u=a
       call svdcmp(u,w,v)
       if (err_math_gen) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_sp

    !!---
    !!---- Subroutine Sort(arr,n,indx)
    !!----    integer/real(kind=sp)  dimension(:), intent( in) :: arr
    !!----    integer,                             intent( in) :: n
    !!----    integer,               dimension(:), intent(out) :: indx
    !!----
    !!----    Sort an array such the arr(indx(j)) is in ascending
    !!----    order for j=1,2,...,N.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Sort_I(Arr,N,Indx)
    !!--++    integer, dimension(:), intent( in) :: arr
    !!--++    integer,               intent( in) :: n
    !!--++    integer, dimension(:), intent(out) :: indx
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_I(arr,n,indx)
       !---- Arguments ----!
       integer, dimension(:), intent(in ) :: arr
       integer              , intent(in ) :: n
       integer, dimension(:), intent(out) :: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer                      :: i,indxt,ir,itemp,j,jstack,k,l
       integer                      :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                err_math_gen=.true.
                err_mess_Math_Gen=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_I

    !!--++
    !!--++ Subroutine Sort_R(arr,n,indx)
    !!--++    real(kind=sp),dimension(:), intent( in) :: arr
    !!--++    integer,                    intent( in) :: n
    !!--++    integer,      dimension(:), intent(out) :: indx
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sort an array such the arr(indx(j)) is in ascending
    !!--++    order for j=1,2,...,N.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sort_R(arr,n,indx)
       !---- Arguments ----!
       real(kind=sp),dimension(:), intent(in) :: arr
       integer,                    intent(in) :: n
       integer,      dimension(:), intent(out):: indx

       !---- Local Variables ----!
       integer, parameter           :: m=7
       integer, parameter           :: nstack=50  !nstack=2log2(n)
       integer, dimension(nstack)   :: istack
       integer :: i,indxt,ir,itemp,j,jstack,k,l
       real(kind=sp)    :: a

       call init_Err_MathGen()
       do j=1,n
          indx(j)=j
       end do

       jstack=0
       l=1
       ir=n
       do
          if (ir-l < m) then
             doext: do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if (arr(indx(i)) <= a)  then
                      indx(i+1)=indxt
                      cycle doext
                   end if
                   indx(i+1)=indx(i)
                end do
                i=0
                indx(i+1)=indxt
             end do doext

             if (jstack == 0) exit
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if (arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             end if
             if (arr(indx(l+1)) > arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             end if
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
             do
                i=i+1
                if (arr(indx(i)) < a)  cycle
                do
                   j=j-1
                   if (arr(indx(j)) > a) cycle
                   exit
                end do
                if (j < i) exit
                itemp=indx(i)
                indx(i)=indx(j)
                indx(j)=itemp
             end do
             indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if (jstack > nstack) then
                err_math_gen=.true.
                err_mess_Math_Gen=" NSTACK too small in SORT"
                return
             end if
             if (ir-i+1 >= j-l) then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do

       return
    End Subroutine Sort_R

    !!---
    !!---- Subroutine Sort_Strings(arr)
    !!----    character(len=*) dimension(:), intent( in out) :: arr
    !!----
    !!----    Sort an array of string
    !!----
    !!---- Update: March - 2005
    !!
    Recursive Subroutine Sort_Strings(Arr)
       !---- Argument ----!
       character(len=*), intent(in out), dimension(:) :: Arr

       !---- Local variables ----!
       integer :: iq

       if (size(Arr) > 1) then
          call Partition(Arr, iq)
          call Sort_Strings(Arr(:iq-1))
          call Sort_Strings(Arr(iq:))
       end if

       return
    End Subroutine Sort_Strings

    !!----
    !!---- Subroutine Spline(x, y, n, yp1, ypn, y2)
    !!----    real(kind=sp),    intent(in),     dimension(n) :: x     !  In -> Array X
    !!----    real(kind=sp),    intent(in),     dimension(n) :: y     !  In -> Array Yi=F(Xi)
    !!----    integer ,         intent(in)                   :: n     !  In -> Dimension of X, Y
    !!----    real(kind=sp),    intent(in)                   :: yp1   !  In -> Derivate of Point 1
    !!----    real(kind=sp),    intent(in)                   :: ypn   !  In -> Derivate of Point N
    !!----    real(kind=sp),    intent(out),    dimension(n) :: y2    ! Out -> array containing second derivatives
    !!----                                                                     at the given points
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spline(x,y,n,yp1,ypn,y2)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in)  :: x
       real(kind=sp), dimension(:), intent(in)  :: y
       integer ,                    intent(in)  :: n
       real(kind=sp),               intent(in)  :: yp1
       real(kind=sp),               intent(in)  :: ypn
       real(kind=sp), dimension(:), intent(out) :: y2

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=sp), dimension(n) :: u
       real(kind=sp)               :: sig, p, qn, un

       if (yp1 > 0.99e30) then
          y2(1)=0.0
          u(1)=0.0
       else
          y2(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       end if

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.0
          y2(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
       if (ypn > 0.99e30) then
          qn=0.0
          un=0.0
       else
          qn=0.5
          un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       end if
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
       do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
       end do

       return
    End Subroutine Spline

    !!----
    !!---- Subroutine Splint(xa, ya, y2a, n, x, y)
    !!----    real(kind=sp),    intent(in), dimension(n) :: xa  !  In -> Array X
    !!----    real(kind=sp),    intent(in), dimension(n) :: ya  !  In -> Array Y=F(X)
    !!----    real(kind=sp),    intent(in), dimension(n) :: y2a !  In -> Array Second Derivatives in X
    !!----    integer ,         intent(in)               :: n   !  In -> Dimension of XA,YA,Y2A
    !!----    real(kind=sp),    intent(in)               :: x   !  In -> Point to evaluate
    !!----    real(kind=sp),    intent(out),             :: y   ! Out -> Interpoled value
    !!----
    !!----    Spline Interpolation
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Splint(xa,ya,y2a,n,x,y)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in)  :: xa
       real(kind=sp), dimension(:), intent(in)  :: ya
       real(kind=sp), dimension(:), intent(in)  :: y2a
       integer ,                    intent(in)  :: n
       real(kind=sp),               intent(in)  :: x
       real(kind=sp),               intent(out) :: y

       !---- Local Variables ----!
       INTEGER          :: klo, khi, k
       real(kind=sp)    :: h, a, b

       klo=1
       khi=n
       do
          if (khi-klo > 1) then
             k=(khi+klo)/2
             if (xa(k) > x) then
                khi=k
             else
                klo=k
             end if
             cycle
          else
             exit
          end if
       end do

       h=xa(khi)-xa(klo)
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)* y2a(khi))*(h**2)/6.0

       return
    End Subroutine Splint

    !!----
    !!---- Subroutine Svdcmp(a,w,v)
    !!----    real(sp/dp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!----    real(sp/dp),dimension(:),  intent(   out) :: w  !W(n)
    !!----    real(sp/dp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--<<
    !!----    Given an M�N matrix A ,this routine computes its singular value decomposition,
    !!----    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!----    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!----    (not the transpose VT )is output as v .
    !!----    Adapted from Numerical Recipes. Valid for arbitrary real matrices
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Svdcmp_dp(a,w,v)
    !!--++    real(dp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!--++    real(dp),dimension(:),  intent(   out) :: w  !W(n)
    !!--++    real(dp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real(kind=sp) matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_dp(a,w,v)
       !---- Arguments ----!
       real(kind=dp),dimension(:,:),intent(in out) ::a
       real(kind=dp),dimension(:),  intent(   out) ::w
       real(kind=dp),dimension(:,:),intent(   out) ::v

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=dp)                               ::anorm,c,f,g,h,s,scal,x,y,z
       real(kind=dp),dimension(size(a,1))          ::tempm
       real(kind=dp),dimension(size(a,2))          ::rv1,tempn

       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          err_math_gen = .true.
          err_mess_math_gen = " => Physical dimensions of arguments in SVDcmp_dp are not compatible "
          return
       end if
       g=0.0_dp
       scal=0.0_dp
       do i=1,n
          l=i+1
          rv1(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if (i <=m)then
             scal=sum(abs(a(i:m,i)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i:m,i)=a(i:m,i)/scal
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scal*a(i:m,i)
             end if
          end if
          w(i)=scal*g
          g=0.0_dp
          scal=0.0_dp
          if ((i <=m).and.(i /=n))then
             scal=sum(abs(a(i,l:n)))
             if ( abs(scal) > tiny(1.0_dp) ) then
                a(i,l:n)=a(i,l:n)/scal
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scal*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1
          if (i <n) then
             if ( abs(g) > tiny(1.0_dp) ) then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0_dp
             v(l:n,i)=0.0_dp
          end if
          v(i,i)=1.0_dp
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1
          l=i+1
          g=w(i)
          a(i,l:n)=0.0_dp
          if ( abs(g) > tiny(1.0_dp) ) then
             g=1.0_dp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0_dp
          end if
          a(i,i)=a(i,i)+1.0_dp
       end do
       do k=n,1,-1
          do its=1,num_its
             do l=k,1,-1
                nm=l-1
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0_dp
                   s=1.0_dp
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_dp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k)then
                if (z <0.0_dp)then
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                err_math_gen = .true.
                err_mess_math_gen = " => SVDcmp_dp: convergence not reached ! "
                return
             end if
             x=w(l)
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
             g=pythag(f,1.0_dp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0_dp
             s=1.0_dp
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z
                if ( abs(z) > tiny(1.0_dp) ) then
                   z=1.0_dp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0_dp
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_dp

    !!--++
    !!--++ Subroutine Svdcmp_sp(a,w,v)
    !!--++    real(sp),dimension(:,:),intent(in out) :: a  !A(m,n)
    !!--++    real(sp),dimension(:),  intent(   out) :: w  !W(n)
    !!--++    real(sp),dimension(:,:),intent(   out) :: v  !V(n,n)
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Given an M �N matrix A ,this routine computes its singular value decomposition,
    !!--++    A = U �W �VT . The matrix U replaces A on output. The diagonal matrix of
    !!--++    singular values W is output as the N-dimensional vector w. The N�N matrix V
    !!--++    (not the transpose VT )is output as v .
    !!--++    Adapted from Numerical Recipes. Valid for arbitrary real(kind=sp) matrices
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Svdcmp_sp(a,w,v)
       !---- Arguments ----!
       real(kind=sp),dimension(:,:),intent(in out) :: a
       real(kind=sp),dimension(:),  intent(   out) :: w
       real(kind=sp),dimension(:,:),intent(   out) :: v

       !---- Local variables ----!
       integer, parameter                          :: num_its=500
       integer                                     ::i,its,j,k,l,m,n,nm
       real(kind=sp)                               ::anorm,c,f,g,h,s,scala,x,y,z
       real(kind=sp),dimension(size(a,1))          ::tempm
       real(kind=sp),dimension(size(a,2))          ::rv1,tempn


       m=size(a,1)
       n=size(a,2)
       call init_err_mathgen()
       if ( .not. (size(v,1) == n .and. size(v,2) == n .and. size(w) == n)) then
          err_math_gen = .true.
          err_mess_math_gen = " => Physical dimensions of arguments in SVDcmp_sp are not compatible "
          return
       end if
       g=0.0
       scala=0.0
       do i=1,n                        !Householder reduction to bidiagonal form.
          l=i+1
          rv1(i)=scala*g
          g=0.0
          scala=0.0
          if (i <=m)then
             scala=sum(abs(a(i:m,i)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i:m,i)=a(i:m,i)/scala
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scala*a(i:m,i)
             end if
          end if
          w(i)=scala*g
          g=0.0
          scala=0.0
          if ((i <=m).and.(i /=n))then
             scala=sum(abs(a(i,l:n)))
             if (abs(scala) > tiny(1.0_sp))then
                a(i,l:n)=a(i,l:n)/scala
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scala*a(i,l:n)
             end if
          end if
       end do
       anorm=maxval(abs(w)+abs(rv1))
       do i=n,1,-1                    ! Accumulation of right-hand transformations.
          if (i <n)then
             if (abs(g) > tiny(1.0_sp))then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g   !Double division to avoid possible underflow.
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
             end if
             v(i,l:n)=0.0
             v(l:n,i)=0.0
          end if
          v(i,i)=1.0
          g=rv1(i)
          l=i
       end do
       do i=min(m,n),1,-1  !Accumulation of left-hand transformations.
          l=i+1
          g=w(i)
          a(i,l:n)=0.0
          if (abs(g) > tiny(1.0_sp))then
             g=1.0_sp/g
             tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=a(i:m,i)*g
          else
             a(i:m,i)=0.0
          end if
          a(i,i)=a(i,i)+1.0_sp
       end do
       do k=n,1,-1         !Diagonalization of the idiagonal form:Loop over
          do its=1,num_its   !singular values,and over allowed iterations.
             do l=k,1,-1     !Test for splitting.
                nm=l-1        !Note that rv1(1)is always zero,so can never fall through bottom of loop.
                if ((abs(rv1(l))+anorm)==anorm) exit
                if ((abs(w(nm))+anorm)==anorm) then
                   c=0.0       ! Cancellation of rv1(l),if l >1 .
                   s=1.0
                   do i=l,k
                      f=s*rv1(i)
                      rv1(i)=c*rv1(i)
                      if ((abs(f)+anorm)==anorm)exit
                      g=w(i)
                      h=pythag(f,g)
                      w(i)=h
                      h=1.0_sp/h
                      c=(g*h)
                      s=-(f*h)
                      tempm(1:m)=a(1:m,nm)
                      a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                   end do
                   exit
                end if
             end do
             z=w(k)
             if (l ==k) then    !Convergence.
                if (z <0.0)then !Singular value is made nonnegative.
                   w(k)=-z
                   v(1:n,k)=-v(1:n,k)
                end if
                exit
             end if
             if (its == num_its) then
                err_math_gen = .true.
                err_mess_math_gen = " => SVDcmp_sp: convergence not reached ! "
                return
             end if
             x=w(l)             !Shift from ottom 2-y-2 minor.
             nm=k-1
             y=w(nm)
             g=rv1(nm)
             h=rv1(k)
             f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
             g=pythag(f,1.0_sp)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.0  ! Next QR transformation:
             s=1.0
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z                 !Rotation can e arbitrary if z =0 .
                if (abs(z) > tiny(1.0_sp) )then
                   z=1.0_sp/z
                   c=f*z
                   s=h*z
                end if
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
             end do
             rv1(l)=0.0
             rv1(k)=f
             w(k)=x
          end do
       end do

       return
    End Subroutine Svdcmp_sp

    !!----
    !!---- Subroutine Swap(a,b) or Swap(a,b,mask)
    !!----    integer,real(sp),complex, intent( in out) :: a, b
    !!----      or
    !!----    integer,real(sp),complex, dimension(:), intent( in out) :: a, b
    !!----      or
    !!----    integer,real(sp),complex, dimension(:,:), intent( in out) :: a, b
    !!----      or
    !!----    real(kind=sp),  intent(in out) :: a,b
    !!----    logical,        intent(in)     :: mask
    !!----      or
    !!----    real(kind=sp), dimension(:), intent(in out) :: a,b
    !!----    logical,       dimension(:), intent(in)     :: mask
    !!----      or
    !!----    real(kind=sp), dimension(:,:), intent(in out) :: a,b
    !!----    logical,       dimension(:,:), intent(in)     :: mask
    !!----
    !!----    Swap the contents of a and b, when mask (if given) is true.
    !!----
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Swap_C(a,b)
    !!--++    complex, intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_C(a,b)
       !---- Arguments ----!
       complex, intent(in out) :: a
       complex, intent(in out) :: b

       !---- Local variables ----!
       complex :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_C

    !!--++
    !!--++ Subroutine Swap_Cm(A,B)
    !!--++    complex, dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cm(a,b)
       !---- Arguments ----!
       complex, dimension(:,:), intent(in out) :: a
       complex, dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cm

    !!--++
    !!--++ Subroutine Swap_Cv(a,b)
    !!--++    complex, dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Cv(a,b)
       !---- Arguments ----!
       complex, dimension(:), intent(in out) :: a
       complex, dimension(:), intent(in out) :: b

       !---- Local variables ----!
       complex, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Cv

    !!--++
    !!--++ Subroutine Swap_I(A,B)
    !!--++    integer , intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_I(A,B)
       !---- Arguments ----!
       integer , intent(in out) :: a
       integer , intent(in out) :: b

       !---- Local variables ----!
       integer  :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_I

    !!--++
    !!--++ Subroutine Swap_Im(A,B)
    !!--++    integer, dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Im(A,B)
       !---- Arguments ----!
       integer, dimension(:,:), intent(in out) :: a
       integer, dimension(:,:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Im

    !!--++
    !!--++ Subroutine Swap_Iv(A,B)
    !!--++    integer, dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Iv(A,B)
       !---- Arguments ----!
       integer, dimension(:), intent(in out) :: a
       integer, dimension(:), intent(in out) :: b

       !---- Local Variables ----!
       integer, dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Iv

    !!--++
    !!--++ Subroutine Swap_R(A,B)
    !!--++    real(kind=sp) , intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_R(A,B)
       !---- Arguments ----!
       real(kind=sp), intent(in out) :: a
       real(kind=sp), intent(in out) :: b

       !---- Local variables ----!
       real(kind=sp) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_R

    !!--++
    !!--++ Subroutine Swap_Rm(A,B)
    !!--++    real(kind=sp), dimension(:,:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rm(A,B)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a
       real(kind=sp), dimension(:,:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rm

    !!--++
    !!--++ Subroutine Swap_Rv(A,B)
    !!--++    real(kind=sp), dimension(:), intent(in out) :: a,b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Swap_Rv(A,B)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in out) :: a
       real(kind=sp), dimension(:), intent(in out) :: b

       !---- Local variables ----!
       real(kind=sp), dimension(size(a)) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_Rv

    !!--++
    !!--++ Subroutine Masked_Swap_R(A,B,Mask)
    !!--++    real(kind=sp), intent(in out) :: a,b
    !!--++    logical,           intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b if mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_R(A,B,Mask)
       !---- Arguments ----!
       real(kind=sp), intent(in out) :: a
       real(kind=sp), intent(in out) :: b
       logical,           intent(in) :: mask

       !---- Local Variables ----!
       real(kind=sp) :: swp

       if (mask) then
          swp=a
          a=b
          b=swp
       end if

       return
    End Subroutine Masked_Swap_R

    !!--++
    !!--++ Subroutine Masked_Swap_Rm(A,B,Mask)
    !!--++    real(kind=sp), dimension(:,:),intent(in out) :: a,b
    !!--++    logical,       dimension(:,:),    intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rm(A,B,Mask)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a
       real(kind=sp), dimension(:,:), intent(in out) :: b
       logical,       dimension(:,:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2)) :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rm

    !!--++
    !!--++ Subroutine Masked_Swap_Rv(A,B,Mask)
    !!--++    real(kind=sp), dimension(:),intent(in out) :: a,b
    !!--++    logical,       dimension(:),    intent(in) :: mask
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b where mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Masked_Swap_Rv(A,B,Mask)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in out) :: a
       real(kind=sp), dimension(:), intent(in out) :: b
       logical,       dimension(:), intent(in)     :: mask

       !---- Local variables ----!
       real(kind=sp), dimension(size(a))           :: swp

       where (mask)
          swp=a
          a=b
          b=swp
       end where

       return
    End Subroutine Masked_Swap_Rv

    !!--++
    !!--++ Subroutine Tqli1(d,e,n)
    !!--++    real(kind=sp), dimension(:), intent (in out):: d
    !!--++    real(kind=sp), dimension(:), intent (in out):: e
    !!--++    integer,                     intent (in)    :: n
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real(kind=sp) tridiagonal symmetric matrix, or of
    !!--++    a real(kind=sp) symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    In TLQ1 only the eigenvalues are calculated
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tqli1(d,e,n)
       !---- Arguments ----!
       real(kind=sp), dimension(:), intent(in out):: d, e ! d(np),e(np)
       integer,                     intent(in )   :: n

       !---- Local variables ----!
       integer      :: i, iter, l, m, mv
       real(kind=sp):: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do
       e(n)=0.0
       do l=1,n
          iter=0
          do_g : do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv

             if (m /= l) then
                if (iter == 40) then
                   err_math_gen=.true.
                   err_mess_Math_Gen=" Too many iterations in TQLI1"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r)  <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli1

    !!--++
    !!--++ Subroutine Tqli2(d,e,n,z)
    !!--++    real(kind=sp), dimension(:)  , intent (in out):: d
    !!--++    real(kind=sp), dimension(:)  , intent (in out):: e
    !!--++    integer,                       intent (in)    :: n
    !!--++    real(kind=sp), dimension(:,:), intent (in out):: z
    !!--++
    !!--++    (PRIVATE)
    !!--++    QL-algorithm with implicit shifts, to determine the eigenvalues
    !!--++    and eigenvectors of a real(kind=sp) tridiagonal symmetric matrix, or of
    !!--++    a real(kind=sp) symmetric matrix previously reduced by tred. D is a vector
    !!--++    with the diagonal elements of the tridiagonal matrix. on output
    !!--++    it returns the eigenvalues. the vector e inputs the subdiagonal
    !!--++    elements of the tridiagonal matrix, with E(1) arbitrary. on
    !!--++    output e is destroyed.
    !!--++    The eigenvectors of the tridiagonal matrix are calculated in TLQ2
    !!--++    by providing the matrix Z  as the identity matrix on input. if the
    !!--++    eigenvectors of the matrix reduced by tred are required, then Z
    !!--++    is input as the matrix output of tred. in either cased, the k-th
    !!--++    column of Z returns the mormalized eigenvector corresponding to
    !!--++    D(k).
    !!--++
    !!--++  Update: February - 2005
    !!
    Subroutine Tqli2(d,e,n,z)
       !---- Arguments ----!
       real(kind=sp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)
       integer,                       intent(in )    :: n
       real(kind=sp), dimension(:,:), intent(in out) :: z    ! z(np,np)

       !---- Local Variables ----!
       integer       :: i, iter, k, l, m, mv
       real(kind=sp) :: b, c, dd, f, g, p, r, s, comp

       call init_Err_MathGen()
       do i=2,n
          e(i-1)=e(i)
       end do

       e(n)=0.0
       do l=1,n
          iter=0
          do_g: do
             mv=n
             do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                comp= abs(e(m))+dd
                if (abs(comp-dd) <= ep_ss) then
                   mv=m
                   exit
                end if
             end do
             m=mv
             if (m /= l) then
                if (iter == 40) then
                   err_math_gen=.true.
                   err_mess_Math_Gen=" Too many iterations in TQLI2"
                   exit
                end if

                iter=iter+1
                g=(d(l+1)-d(l))/(2.0*e(l))
                r=sqrt(g*g+1.0)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                   f=s*e(i)
                   b=c*e(i)
                   r=sqrt(f*f+g*g)
                   e(i+1)=r
                   if (abs(r) <= ep_ss) then
                      d(i+1)=d(i+1)-p
                      e(m)=0.0
                      cycle do_g
                   end if
                   s=f/r
                   c=g/r
                   g=d(i+1)-p
                   r=(d(i)-g)*s+2.0*c*b
                   p=s*r
                   d(i+1)=g+p
                   g=c*r-b

                   !---- omit lines from here ...
                   do k=1,n
                      f=z(k,i+1)
                      z(k,i+1)=s*z(k,i)+c*f
                      z(k,i)=c*z(k,i)-s*f
                   end do

                   !---- ... to here when finding only eigenvalues.
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
                cycle do_g
             end if
             exit
          end do do_g
       end do

       return
    End Subroutine Tqli2

    !!--++
    !!--++ Subroutine Tred1(a,n,d,e)
    !!--++    real(kind=sp), dimension(:,:), intent (in out):: a
    !!--++    integer,                       intent (in)    :: n
    !!--++    real(kind=sp), dimension(:)  , intent (in out):: d
    !!--++    real(kind=sp), dimension(:)  , intent (in out):: e
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find only eigenvalues
    !!--++    Householder reduction of a real(kind=sp) symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++    In tred1 several lines have been deleted and A contains no
    !!--++    useful information on output.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred1(a,n,d,e)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=sp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local Variables ----!
       integer :: i, j, k, l
       real(kind=sp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
                hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       e(1)=0.0
       do i=1,n
          d(i)=a(i,i)
       end do

       return
    End Subroutine Tred1

    !!--++
    !!--++ Subroutine Tred2(a,n,d,e)
    !!--++    real(kind=sp), dimension(:,:), intent (in out) :: a
    !!--++    integer,                       intent (in)     :: n
    !!--++    real(kind=sp), dimension(:)  , intent (in out) :: d
    !!--++    real(kind=sp), dimension(:)  , intent (in out) :: e
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine for preparing the matrix to find the complete set
    !!--++    of eigenvectors.
    !!--++    Householder reduction of a real(kind=sp) symetric nxn matrix A.
    !!--++    On output A is replaced by the orthogonal matrix Q effecting
    !!--++    the transformation. D returns the diagonal elements of the tri-
    !!--++    diagonal matrix and E the off-diagonal elements with E(1)=0.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Tred2(a,n,d,e)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:), intent(in out) :: a    ! a(np,np)
       integer,                       intent(in)     :: n
       real(kind=sp), dimension(:),   intent(in out) :: d, e ! d(np),e(np)

       !---- Local variables ----!
       integer :: i, j, k, l
       real(kind=sp)    :: f, g, h, hh, scala

       do i=n,2,-1
          l=i-1
          h=0.0
          scala=0.0
          if (l > 1)then
             do k=1,l
                scala=scala+abs(a(i,k))
             end do
             if (abs(scala) <= ep_ss) then
                e(i)=a(i,l)
             else
                do k=1,l
                   a(i,k)=a(i,k)/scala
                   h=h+a(i,k)**2
                end do
                f=a(i,l)
                g=-sign(sqrt(h),f)
                e(i)=scala*g
                h=h-f*g
                a(i,l)=f-g
                f=0.0
                do j=1,l
                   !---- omit following line if finding only eigenvalues
                   a(j,i)=a(i,j)/h
                   g=0.0
                   do k=1,j
                      g=g+a(j,k)*a(i,k)
                   end do
                   do k=j+1,l
                      g=g+a(k,j)*a(i,k)
                   end do
                   e(j)=g/h
                   f=f+e(j)*a(i,j)
                end do
               hh=f/(h+h)
                do j=1,l
                   f=a(i,j)
                   g=e(j)-hh*f
                   e(j)=g
                   do k=1,j
                      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                   end do
                end do
             end if
          else
             e(i)=a(i,l)
          end if
          d(i)=h
       end do

       !---- omit following line if finding only eigenvalues.
       d(1)=0.0
       e(1)=0.0
       do i=1,n
          !---- delete lines from here ...
          l=i-1
          if (abs(d(i)) > ep_ss)then
             do j=1,l
                g=0.0
                do k=1,l
                   g=g+a(i,k)*a(k,j)
                end do
                do k=1,l
                   a(k,j)=a(k,j)-g*a(k,i)
                end do
             end do
          end if
          !---- ... to here when finding only eigenvalues.
          d(i)=a(i,i)
          !---- also delete lines from here ...
          a(i,i)=1.0
          do j=1,l
             a(i,j)=0.0
             a(j,i)=0.0
          end do
          !---- ... to here when finding only eigenvalues.
       end do

       return
    End Subroutine Tred2

 End Module Math_Gen

