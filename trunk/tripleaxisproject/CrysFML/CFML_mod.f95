!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: MOD_FUN
!!----   INFO:  Module defining the function MOD for F-compilers.
!!----          This public  function overrides the intrinsic Fortran MOD
!!----
!!---- HISTORY
!!----    Update:January  - 2005
!!----           October  - 2002  Created by JRC
!!----
!!---- VARIABLES
!!--++    DP        [Private]
!!--++    SP        [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       MOD               [Overrides intrinsic MOD]
!!--++       MOD_DP            [Overloaded]
!!--++       MOD_INT           [Overloaded]
!!--++       MOD_SP            [Overloaded]
!!----
!!
 Module Mod_Fun

    !---- Variables ----!
    implicit none

    !---- List of public functions ----!
    public :: Mod

    !---- List of private functions ----!
    private :: Mod_Dp, Mod_Int, Mod_Sp

    !---- Definitions ----!

    !!--++
    !!--++ DP
    !!--++    integer, private, parameter :: dp
    !!--++
    !!--++    (PRIVATE)
    !!--++    Double precision ( kind(1.0d0) )
    !!--++
    !!--++ Update: February - 2003
    !!
    integer, private, parameter :: dp = selected_real_kind(14,80)

    !!--++
    !!--++ SP
    !!--++    integer, private, parameter :: sp
    !!--++
    !!--++    (PRIVATE)
    !!--++    Single precision ( kind(1.0) )
    !!--++
    !!--++ Update: February - 2003
    !!
    integer, private, parameter :: sp = selected_real_kind(6,30)


    !---- Interfaces - Overlapp ----!

    Interface  Mod     !Redefinition of the intrinsic MOD of Fortran (to be used with F)
       Module Procedure Mod_Dp
       Module Procedure Mod_Int
       Module Procedure Mod_Sp
    End Interface

 Contains

    !!----
    !!---- Function Mod(x,y) Result(md)
    !!----    real(sp/dp)/integer, intent(in) : x,y
    !!----    real(sp/dp)/integer             : md
    !!----
    !!----    Overrides intrinsic MOD function of Fortran
    !!----    (Not to be used for non-F compilers)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Mod_Dp(x,y) Result(md)
    !!--++    real(dp), intent(in) : x,y
    !!--++    real(dp)             : md
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Overrides intrinsic MOD function of Fortran
    !!--++    (Not to be used for non-F compilers)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Mod_Dp(x,y) Result(md)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x,y
       real(kind=dp)             :: md

       md=modulo(abs(x),abs(y))*sign(1.0_dp,x)

       return
    End Function Mod_Dp

    !!--++
    !!--++ Function Mod_Int(x,y) Result(md)
    !!--++    integer, intent(in) : x,y
    !!--++    integer             : md
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Overrides intrinsic MOD function of Fortran
    !!--++    (Not to be used for non-F compilers)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Mod_Int(x,y) Result(md)
       !---- Arguments ----!
       integer, intent(in) :: x,y
       integer             :: md

       md=modulo(abs(x),abs(y))*sign(1,x)

       return
    End Function Mod_Int

    !!--++
    !!--++ Function Mod_Sp(x,y) Result(md)
    !!--++    real(sp), intent(in) : x,y
    !!--++    real(sp)             : md
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Overrides intrinsic MOD function of Fortran
    !!--++    (Not to be used for non-F compilers)
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Mod_sp(x,y) result(md)
       !---- Arguments ----!
       real(kind=sp), intent(in) :: x,y
       real(kind=sp)             :: md

       md=modulo(abs(x),abs(y))*sign(1.0_sp,x)

       return
    End Function Mod_Sp

 End Module Mod_fun
