!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!----
!!---- MODULE: CRYSTAL_TYPES
!!----   INFO: Module to define crystallographic types and to provide
!!----         automatic crystallographic operations.
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----
!!----       July - 2000 Re-organised and updated by JGP and JRC.
!!----    October - 1997 Created by JRC
!!----
!!--.. INFORMATION
!!--..
!!--..    List Of Matrix Relationships For Crystallographic Applications
!!--..
!!--..    Small "t" is for transpose, inv(F) is the inverse of matrix F
!!--..
!!--..    Basis vectors as symbolic matrices
!!--..       At = (a,b,c)  At'=(a',b',c') ;  At* = (a*,b*,c*)  At*'=(a*',b*',c*')
!!--..
!!--..    Direct and reciprocal metric tensors: G, G*=inv(G)
!!--..    X  column vector in     direct space, referred to basis A
!!--..    X* column vector in reciprocal space, referred to basis A*
!!--..
!!--..       A'  = M  A           X'  = inv(Mt) X
!!--..       A*  = G* A           X*  =   G     X
!!--..       A*' = inv(Mt) A*     X*' =   M     X*
!!--..
!!--..       G' = M G Mt          G*' = inv(Mt) G* inv(M)
!!--..
!!--..   Symmetry operator defined in bases: A, A', A*, A*'
!!--..       C = (R,T), C'= (R',T'), C*= (R*,T*), C*'= (R*',T*')
!!--..
!!--..       R'  = inv(Mt) R Mt  ; T' = inv(Mt) T
!!--..       R*' =  M  R* inv(M) ; T*' = M T*
!!--..       R*  = G R G*  = inv(Rt)
!!--..
!!--..   If a change of origin is performed the positions are changed
!!--..   Ot=(o1,o2,o3) origin of the new basis A' w.r.t. old basis A
!!--..
!!--..       X' = inv(Mt) (X-O)
!!--..
!!--..   Changing just the origin   Xn  = C  X  = R  X  + T
!!--..                              Xn' = C' X' = R' X' + T'
!!--..          R=R'                X'  = X -O
!!--..                              Xn' = Xn-O
!!--..                  Xn-O = R' (X-O) + T' = R X + T - O
!!--..                   R X - R O + T' = R X + T - O
!!--..                               T' = T - (O - R O) = T - (E-R)O
!!--..
!!--..   Changing the basis (A,o) -> (A',o')
!!--..                  Xn  = C  X  = R  X  + T
!!--..                  Xn' = C' X' = R' X' + T'
!!--..                  X'= inv(Mt) (X-O), Xn' = inv(Mt) (Xn-O)
!!--..
!!--..            inv(Mt) (Xn-O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) (R  X  + T -O) = R' inv(Mt) (X-O) + T'
!!--..            inv(Mt) R X + inv(Mt)(T-O) = R' inv(Mt) X - R' inv(Mt) O + T'
!!--..            inv(Mt) R = R' inv(Mt)  => R' = inv(Mt) R Mt
!!--..            inv(Mt) (T-O)  = - R' inv(Mt) O + T'
!!--..            T' = R' inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R Mt inv(Mt) O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) (T-O)
!!--..            T' = inv(Mt) R  O + inv(Mt) T - inv(Mt) O
!!--..            T' = inv(Mt)( R  O + T -  O) = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..                       R' = inv(Mt) R Mt
!!--..
!!--..                       T' = inv(Mt) (T -(E-R)O)
!!--..
!!--..
!!--..   A symmetry operator does not change the modulus of vectors and
!!--..   the angles between vectors (dot product is invariant):
!!--..
!!--..      X' = R X ,  Y' = R Y  =>  Xt' = Xt Rt,  Yt' = Yt Rt
!!--..
!!--..      Xt' G Y' = Xt Rt G R Y = Xt G Y  =>  G = Rt G R
!!--..
!!--..
!!--..   Second rank tensor Q and Q* defined in bases A and A*.
!!--..
!!--..      Q' = M Q Mt      Q* = G* Q G*     Q*'= inv(Mt) Q* inv(M)
!!--..
!!--..   A symmetry operator R is equivalent to a transformation
!!--..   M = inv(Rt) acting on basis vectors => G' = inv(Rt) G inv(R) = G
!!--..   The anisotropic temperature factors Beta is defined in reciprocal
!!--..   space: is a tensor like Q*, the transformation of beta under a
!!--..   symmetry operator is then :
!!--..
!!--..           Beta' = Inv(Mt) Beta inv(M) = R Beta Rt
!!--..
!!----
!!---- DEPENDENCIES
!!--++    Use MATH_GEN: Sp, Eps, Cosd, Sind, Acosd, Pi
!!--++    Use MATH_3D : Matrix_Inverse
!!----
!!---- VARIABLES
!!----    CRYSTAL_CELL_TYPE
!!----    ERR_CRYS
!!----    ERR_MESS_CRYS
!!--++    TPI2                           [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CART_U_VECTOR
!!----       CART_VECTOR
!!----       CONVERT_B_BETAS
!!----       CONVERT_B_U
!!----       CONVERT_BETAS_B
!!----       CONVERT_BETAS_U
!!----       CONVERT_U_B
!!----       CONVERT_U_BETAS
!!--++       METRICS                     [Private]
!!----       ROT_MATRIX
!!----       U_EQUIV
!!----
!!----    Subroutines:
!!----       CHANGE_SETTING_CELL
!!----       GET_CRYST_FAMILY
!!--++       GET_CRYST_ORTHOG_MATRIX     [Private]
!!----       INIT_ERR_CRYS
!!--++       RECIP                       [Private]
!!----       SET_CRYSTAL_CELL
!!----       WRITE_CRYSTAL_CELL
!!----
!!
 Module Crystal_Types

    !---- Use files ----!
    use Math_Gen , only: cosd, sind, acosd, eps, sp, pi
    use Math_3D,   only: matrix_inverse

    implicit none

    private

    !---- List of public variables ----!

    !---- List of public functions ----!
    public :: Cart_u_vector, Cart_vector, Convert_B_Betas, Convert_B_U, &
              Convert_Betas_B, Convert_Betas_U, Convert_U_B,            &
              Convert_U_Betas, Rot_matrix, U_Equiv

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public :: Init_err_crys, Change_Setting_Cell,Set_Crystal_Cell,   &
              Get_Cryst_Family, Write_Crystal_Cell, get_deriv_Orth_cell

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: metrics

    !---- List of private Subroutines ----!
    private :: Recip, Get_Cryst_Orthog_Matrix

    !---- Definitions ----!

    !!----
    !!----  TYPE :: CRYSTAL_CELL_TYPE
    !!--..
    !!----  Type, public :: Crystal_Cell_Type
    !!----     real(kind=sp),dimension(3)   :: cell, ang          ! Direct cell parameters
    !!----     real(kind=sp),dimension(3)   :: cell_std, ang_std  ! Standar deviations cell parameters
    !!----     real(kind=sp),dimension(3)   :: rcell,rang         ! Reciprocal cell parameters
    !!----     real(kind=sp),dimension(3,3) :: GD,GR              ! Direct and reciprocal Metric Tensors
    !!----     real(kind=sp),dimension(3,3) :: Cr_Orth_cel        ! P-Matrix transforming Orthonormal
    !!----                                                        ! basis to direct Crytal cell (as I.T.)
    !!----                                                        ! (or crystallographic components to
    !!----                                                        !  cartesian components)
    !!----     real(kind=sp),dimension(3,3) :: Orth_Cr_cel        ! Inv(Cr_Orth_cel) -> Cartesian to cryst. components
    !!----     real(kind=sp)                :: CellVol            ! Direct and Reciprocal
    !!----     real(kind=sp)                :: RCellVol           ! Cell volumes
    !!----     Character (len=1)            :: CartType           ! Cartesian Frame type: if CartType='A'
    !!----                                                        ! the Cartesian Frame has x // a.
    !!----  End Type Crystal_Cell_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Crystal_Cell_Type
       real(kind=sp),dimension(3)   :: cell, ang
       real(kind=sp),dimension(3)   :: cell_std, ang_std
       real(kind=sp),dimension(3)   :: rcell, rang
       real(kind=sp),dimension(3,3) :: GD,GR
       real(kind=sp),dimension(3,3) :: Cr_Orth_cel
       real(kind=sp),dimension(3,3) :: Orth_Cr_cel
       real(kind=sp)                :: CellVol
       real(kind=sp)                :: RCellVol
       character (len=1)            :: CartType
    End Type Crystal_Cell_Type

    !!----
    !!---- ERR_CRYS
    !!----    logical, public :: err_crys
    !!----
    !!----    Logical Variable indicating an error in CRYSTAL_TYPES module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public          :: Err_Crys

    !!----
    !!---- ERR_MESS_CRYS
    !!----    character(len=150), public :: Err_Mess_Crys
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: Err_Mess_Crys

    !!--++
    !!--++ TPI2
    !!--++    real(kind=sp), parameter :: tpi2=2.0*pi*pi
    !!--++
    !!--++    (PRIVATE)
    !!--++    Two times PI squared
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private :: tpi2=2.0*pi*pi

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Function Cart_U_Vector(Code,V,Celda) Result(Vc)
    !!----    character (len=*),             intent(in) :: code    !  In -> D: Direct, R: Reciprocal
    !!----    real(kind=sp), dimension(3),   intent(in) :: v       !  In -> Vector
    !!----    Type (Crystal_Cell_Type),      intent(in) :: Celda   !  In -> Cell Variable
    !!----    real(kind=sp),dimension(3)                :: vc      ! Out ->
    !!----
    !!----    Convert a vector in crystal space to unitary cartesian components
    !!----
    !!---- Update: February - 2005
    !!
    Function Cart_U_Vector(Code,V,Celda) Result(Vc)
       !---- Arguments ----!
       character (len=*),           intent(in) :: code
       real(kind=sp), dimension(3), intent(in) :: v
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=sp), dimension(3)             :: vc

       !---- Local Variables ----!
       real(kind=sp) :: vmod

       vc=cart_vector(code,v,celda)
       vmod=sqrt(dot_product(vc,vc))
       if (vmod > 0.0) then
          vc=vc/vmod
       end if

       return
    End Function Cart_U_Vector

    !!----
    !!---- Function Cart_Vector(Code,V,Celda) Result(Vc)
    !!----    character (len=*),             intent(in) :: code     !  In -> D: Direct, R: Reciprocal
    !!----    real(kind=sp), dimension(3),   intent(in) :: v        !  In -> Vector
    !!----    Type (Crystal_Cell_Type),      intent(in) :: Celda    !  In -> Cell variable
    !!----    real(kind=sp) dimension(3)                :: vc       ! Out ->
    !!----
    !!----    Convert a vector in crystal space to cartesian components
    !!----
    !!---- Update: February - 2005
    !!
    Function Cart_Vector(Code,V,Celda) Result(Vc)
       !---- Arguments ----!
       character (len=*),           intent(in) :: code
       real(kind=sp), dimension(3), intent(in) :: v
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=sp), dimension(3)             :: vc

       select case (code(1:1))
          case("d","D")
             vc = matmul(celda%Cr_Orth_cel,v)

          case ("r","R")
             vc = matmul(celda%GR,v)
             vc = matmul(celda%Cr_Orth_cel,vc)
       end select

       return
    End Function Cart_Vector

    !!--..
    !!--.. Betas are defined by the following expression of the temperature factor:
    !!--.. Taniso= exp( -(beta11 h^2 + beta22 k^2 + beta33 l^2 + 2 (beta12 h k + beta13 h l + beta23 k l)) )
    !!--.. Taniso= exp( -(bet(1) h^2 + bet(2) k^2 + bet(3) l^2 + 2 (bet(4) h k + bet(5) h l + bet(6) k l)) )
    !!--..
    !!--.. Us are defined by the following expression of the temperature factor:
    !!--.. Taniso= exp( -2pi^2 (h^2 (a*)^2 U11+ k^2 (b*)^2 U22+ l^2 (c*)^2 U33+
    !!--..                2 (h k (a*) (b*) U12+ h l (a*) (c*) U13+  k l (b*) (c*) U23)) )
    !!--..

    !!----
    !!---- Function Convert_B_Betas(B,Cell) Result(Beta)
    !!----    real(kind=sp),dimension(6), intent(in)  :: B
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=sp),dimension(6)              :: Beta
    !!----
    !!----    Convert Thermal factors from B to Betas
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_B_Betas(B,Cell) Result(Beta)
       !---- Arguments ----!
       real(kind=sp),dimension(6), intent(in)  :: B
       type (Crystal_cell_Type),   intent(in)  :: Cell
       real(kind=sp),dimension(6)              :: Beta

       beta(1)=0.25*b(1)*cell%gr(1,1)                ! beta11
       beta(2)=0.25*b(2)*cell%gr(2,2)                ! beta22
       beta(3)=0.25*b(3)*cell%gr(3,3)                ! beta33
       beta(4)=0.25*b(4)*cell%rcell(1)*cell%rcell(2) ! beta12
       beta(5)=0.25*b(5)*cell%rcell(1)*cell%rcell(3) ! beta13
       beta(6)=0.25*b(6)*cell%rcell(2)*cell%rcell(3) ! beta23

       return
    End Function Convert_B_Betas

    !!----
    !!---- Function Convert_B_U(B) Result(U)
    !!----    real(kind=sp),dimension(6), intent(in)  :: B
    !!----    real(kind=sp),dimension(6)              :: U
    !!----
    !!----    Convert Thermal factors from B to U
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_B_U(B) Result(U)
       !---- Arguments ----!
       real(kind=sp),dimension(6),  intent(in)  :: B
       real(kind=sp),dimension(6)               :: U

       u=b/(4.0*tpi2)

       return
    End Function Convert_B_U

    !!----
    !!---- Function Convert_Betas_B(Beta,Cell) Result(B)
    !!----    real(kind=sp),dimension(6), intent(in)  :: Beta
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=sp),dimension(6)              :: B
    !!----
    !!----    Convert Thermal factors from Betas to B
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_Betas_B(Beta,Cell) Result(B)
       !---- Arguments ----!
       real(kind=sp),dimension(6), intent(in)  :: Beta
       type (Crystal_cell_Type),   intent(in)  :: Cell
       real(kind=sp),dimension(6)              :: B

       b(1)=4.0*beta(1)/cell%gr(1,1)                  ! B11
       b(2)=4.0*beta(2)/cell%gr(2,2)                  ! B22
       b(3)=4.0*beta(3)/cell%gr(3,3)                  ! B33
       b(4)=4.0*beta(4)/(cell%rcell(1)*cell%rcell(2)) ! B12
       b(5)=4.0*beta(5)/(cell%rcell(1)*cell%rcell(3)) ! B13
       b(6)=4.0*beta(6)/(cell%rcell(2)*cell%rcell(3)) ! B23

       return
    End Function Convert_Betas_B

    !!----
    !!---- Function Convert_Betas_U(Beta,Cell) Result(U)
    !!----    real(kind=sp),dimension(6), intent(in)  :: Beta
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=sp),dimension(6)              :: U
    !!----
    !!----    Convert Thermal factors from Betas to U
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_Betas_U(Beta,Cell) Result(U)
       !---- Arguments ----!
       real(kind=sp),dimension(6),intent(in)  :: Beta
       type (Crystal_cell_Type),  intent(in)  :: Cell
       real(kind=sp),dimension(6)             :: U

       u(1)=beta(1)/(tpi2*cell%gr(1,1))                ! U11
       u(2)=beta(2)/(tpi2*cell%gr(2,2))                ! U22
       u(3)=beta(3)/(tpi2*cell%gr(3,3))                ! U33
       u(4)=beta(4)/(tpi2*cell%rcell(1)*cell%rcell(2)) ! U12
       u(5)=beta(5)/(tpi2*cell%rcell(1)*cell%rcell(3)) ! U13
       u(6)=beta(6)/(tpi2*cell%rcell(2)*cell%rcell(3)) ! U23

       return
    End Function Convert_Betas_U

    !!----
    !!---- Function Convert_U_B(U) Result(B)
    !!----    real(kind=sp),dimension(6), intent(in)  :: U
    !!----    real(kind=sp),dimension(6)              :: B
    !!----
    !!----    Convert Thermal factors from U to B
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_U_B(U) Result(B)
       !---- Arguments ----!
       real(kind=sp),dimension(6),        intent(in)  :: U
       real(kind=sp),dimension(6)                     :: B

       b=4.0*tpi2*u

       return
    End Function Convert_U_B

    !!----
    !!---- Function Convert_U_Betas(U,Cell) Result(Beta)
    !!----    real(kind=sp),dimension(6), intent(in)  :: U
    !!----    type (Crystal_cell_Type),   intent(in)  :: Cell
    !!----    real(kind=sp),dimension(6)              :: Beta
    !!----
    !!----    Convert Thermal factors from U to Betas
    !!----
    !!---- Update: February - 2003
    !!
    Function Convert_U_Betas(U,Cell) Result(Beta)
       !---- Arguments ----!
       real(kind=sp),dimension(6),intent(in)  :: U
       type (Crystal_cell_Type),  intent(in)  :: Cell
       real(kind=sp),dimension(6)             :: Beta

       beta(1)=tpi2*u(1)*cell%gr(1,1)                ! beta11
       beta(2)=tpi2*u(2)*cell%gr(2,2)                ! beta22
       beta(3)=tpi2*u(3)*cell%gr(3,3)                ! beta33
       beta(4)=tpi2*u(4)*cell%rcell(1)*cell%rcell(2) ! beta12
       beta(5)=tpi2*u(5)*cell%rcell(1)*cell%rcell(3) ! beta13
       beta(6)=tpi2*u(6)*cell%rcell(2)*cell%rcell(3) ! beta23

       return
    End Function Convert_U_Betas

    !!--++
    !!--++ Function Metrics(A,B) Result(G)
    !!--++    real(kind=sp), dimension(3)  , intent(in ) :: a   !  In -> Cell Parameters
    !!--++    real(kind=sp), dimension(3)  , intent(in ) :: b   !  In -> Ang Parameters
    !!--++    real(kind=sp), dimension(3,3)              :: g   ! Out -> Metrics array
    !!--++
    !!--++    (PRIVATE)
    !!--++    Constructs the metric tensor
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Metrics(A,B) Result(G)
       !---- Arguments ----!
       real(kind=sp), dimension(3)  , intent(in ) :: a
       real(kind=sp), dimension(3)  , intent(in ) :: b
       real(kind=sp), dimension(3,3)              :: g

       !---- Local Variables ----!
       integer :: i

       G(1,2)= a(1)*a(2)*cosd(b(3))
       G(1,3)= a(1)*a(3)*cosd(b(2))
       G(2,3)= a(2)*a(3)*cosd(b(1))

       do i=1,3
          G(i,i)= a(i)*a(i)
       end do

       G(2,1)=G(1,2)
       G(3,1)=G(1,3)
       G(3,2)=G(2,3)

       return
    End Function Metrics

    !!----
    !!---- Function Rot_Matrix(U, Phi, Celda)
    !!----    real(kind=sp), dimension(3),        intent(in) :: U
    !!----    real(kind=sp),                      intent(in) :: Phi
    !!----    type (Crystal_Cell_Type), optional, intent(in) :: Celda
    !!----    real(kind=sp), dimension(3,3)                  :: Rm
    !!----
    !!----    Returns the matrix (Gibbs matrix) of the active rotation of "phi" degrees
    !!----    along the "U" direction: R v = v', the vector v is tranformed to vector v'
    !!----    keeping the reference frame unchanged.
    !!----
    !!----    If one wants to calculate the components of the vector "v" in a rotated
    !!----    reference frame it suffices to invoke the function using "-phi".
    !!----    If "Celda" is present, "U" is in "Celda" coordinates,
    !!----    if not "U" is in cartesian coordinates.
    !!----
    !!----
    !!---- Update: February - 2005
    !!
    Function Rot_Matrix(U,Phi,Celda) Result(Rm)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent(in)        :: U
       real(kind=sp), intent(in)                      :: phi
       type (Crystal_Cell_Type), optional, intent(in) :: Celda
       real(kind=sp), dimension(3,3)                  :: RM

       !---- Local variables ----!
       real(kind=sp)               :: c, s, umc, umod
       real(kind=sp), dimension(3) :: UU


       if (present(celda)) then
          uu= matmul(celda%cr_orth_cel,u)
          umod=sqrt(uu(1)**2+uu(2)**2+uu(3)**2)
          if (umod < tiny(1.0)) then
             uu=(/0.0,0.0,1.0/)
          else
             uu= uu/umod
          end if
       else
          if (umod < tiny(1.0)) then
             uu=(/0.0,0.0,1.0/)
          else
             umod=sqrt(u(1)**2+u(2)**2+u(3)**2)
             uu= u/umod
          end if
       end if
       c= cosd(phi)
       s= sind(phi)
       umc = 1.0-c
       rm(1,1)= c+ umc*uu(1)**2
       rm(1,2)= umc*uu(1)*uu(2)- s*uu(3)
       rm(1,3)= umc*uu(1)*uu(3)+ s*uu(2)

       rm(2,1)= umc*uu(2)*uu(1)+ s*uu(3)
       rm(2,2)= c+ umc*uu(2)**2
       rm(2,3)= umc*uu(2)*uu(3)- s*uu(1)

       rm(3,1)= umc*uu(3)*uu(1)- s*uu(2)
       rm(3,2)= umc*uu(3)*uu(2)+ s*uu(1)
       rm(3,3)= c+ umc*uu(3)**2

       return
    End Function Rot_Matrix

    !!----
    !!---- Function U_Equiv(Cell, Th_U) Result(Uequi)
    !!----    type(Crystal_Cell_Type),    intent(in)     :: Cell    !  In -> Cell variable
    !!----    real(kind=sp), dimension(6),intent(in)     :: Th_U    !  In -> U parameters
    !!----
    !!----    Subroutine to obtain the U equiv from U11 U22 U33 U23 U13 U12
    !!----
    !!---- Update: February - 2005
    !!
    Function U_Equiv(Cell, Th_U) Result(Uequi)
       !---- Arguments ----!
       type (Crystal_cell_Type),    intent(in)  :: Cell
       real(kind=sp), dimension(6), intent(in)  :: Th_U
       real(kind=sp)                            :: Uequi

       !---- Local variables ----!
       real(kind=sp)    :: a, b, c, as, bs, cs, cosa, cosb, cosg
       real(kind=sp)    :: u11, u22, u33, u23, u13, u12

       a  =cell%cell(1)
       b  =cell%cell(2)
       c  =cell%cell(3)
       as =cell%rcell(1)
       bs =cell%rcell(2)
       cs =cell%rcell(3)
       cosa=cosd(cell%ang(1))
       cosb=cosd(cell%ang(2))
       cosg=cosd(cell%ang(3))

       u11=Th_u(1)
       u22=Th_u(2)
       u33=Th_u(3)
       u12=Th_u(4)
       u13=Th_u(5)
       u23=Th_u(6)
       uequi= (1.0/3.0) * (u11 * a * a * as * as + &
                           u22 * b * b * bs * bs + &
                           u33 * c * c * cs * cs + &
                           2.0*u12 * a * b * as * bs * cosg + &
                           2.0*u13 * a * c * as * cs * cosb + &
                           2.0*u23 * b * c * bs * cs * cosa )

       return
    End Function U_Equiv

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
    !!----    type (Crystal_Cell_Type),      intent( in)    :: Cell
    !!----    real(kind=sp), dimension (3,3),intent( in)    :: Mat
    !!----    type (Crystal_Cell_Type),      intent(out)    :: Celln
    !!----    character (len=*), optional,   intent (in)    :: matkind
    !!----
    !!----    Calculates a new cell
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Change_Setting_Cell(Cell,Mat,Celln,Matkind)
       !---- Arguments ----!
       type (Crystal_Cell_Type),      intent( in)    :: Cell
       real(kind=sp), dimension (3,3),intent( in)    :: Mat
       type (Crystal_Cell_Type),      intent(out)    :: Celln
       character(len=*),  optional,   intent (in)    :: matkind

       !--- Local variables ---!
       integer                       :: i
       real(kind=sp), dimension(3)   :: cellv,angl
       real(kind=sp), dimension(3,3) :: S,Gn,ST

       if (present(matkind)) then
          if (matkind(1:2) == "it" .or. matkind(1:2) == "IT" ) then
             S=Mat
            ST=transpose(Mat)
          else
             S=transpose(Mat)
            ST=Mat
          end if
       else
          S=transpose(Mat)
         ST=Mat
       end if

       !---- Get the new metric tensor
       !---- GDN= Mat GD MatT  or GDN= ST GD S
       gn=matmul(ST,matmul(Cell%GD,S))

       !---- Calculate new cell parameters from the new metric tensor
       do i=1,3
          Cellv(i)=sqrt(gn(i,i))
       end do
       angl(1)=acosd(Gn(2,3)/(cellv(2)*cellv(3)))
       angl(2)=acosd(Gn(1,3)/(cellv(1)*cellv(3)))
       angl(3)=acosd(Gn(1,2)/(cellv(1)*cellv(2)))

       !---- Construct the new cell
       call Set_Crystal_Cell(cellv,angl,Celln)

       return
    End Subroutine Change_Setting_Cell

    !!----
    !!---- Subroutine Get_Cryst_Family(Cell,Car_Family,Car_Symbol,Car_System)
    !!----    type(Crystal_Cell_type),         intent(in ) :: Cell
    !!----    character(len=*),                intent(out) :: Car_Family
    !!----    character(len=*),                intent(out) :: Car_Symbol
    !!----    character(len=*),                intent(out) :: Car_System
    !!----
    !!---- Obtain the Crystal Family, Symbol and System from cell parameters
    !!----
    !!---- Update: May - 2005
    !!----
    Subroutine Get_Cryst_Family(Cell,Car_Family,Car_Symbol,Car_System)
       !---- Arguments ----!
       type(Crystal_Cell_type),         intent(in ) :: Cell
       character(len=*),                intent(out) :: Car_Family
       character(len=*),                intent(out) :: Car_Symbol
       character(len=*),                intent(out) :: Car_System

       !---- Local variables ----!
       integer, dimension(3) :: icodp, icoda
       integer               :: n1,n2

       Car_Family=" "
       Car_Symbol=" "
       Car_System=" "

       icodp=0
       icoda=0

       !---- Cell Parameters ----!

       !---- a ----!
       icodp(1)=1

       !---- b ----!
       if (abs(cell%cell(2)-cell%cell(1)) <= 0.0001) then
          icodp(2)=icodp(1)
       else
          icodp(2)=2
       end if

       !---- c ----!
       if (abs(cell%cell(3)-cell%cell(1)) <= 0.0001) then
          icodp(3)=icodp(1)
       else
          icodp(3)=3
       end if

       !---- Angles Parameters ----!

       !---- alpha ----!
       icoda(1)=1

       !---- beta ----!
       if (abs(cell%ang(2)-cell%ang(1)) <= 0.0001) then
          icoda(2)=icoda(1)
       else
          icoda(2)=2
       end if

       !---- gamma ----!
       if (abs(cell%ang(3)-cell%ang(1)) <= 0.0001) then
          icoda(3)=icoda(1)
       else
          icoda(3)=3
       end if


       n1=count(icoda==icoda(1))
       n2=count(icodp==icodp(1))
       select case (n1)
          case (1) ! all are differents
                 if (n2 ==1) then
                          Car_Family="Triclinic"
                          Car_Symbol ="a"
                          Car_System ="Triclinic"
                 else
                          err_crys=.true.
                          err_mess_crys=" Error obtaining Crystal Familiy"
                 end if

          case (2) ! two angles are equal
                 if (icoda(1) == icoda(2)) then
                          if (abs(cell%ang(3)-120.0) <= 0.0001) then
                                 if (icodp(1)==icodp(2)) then
                                          !---- Hexagonal ----!
                                          Car_Family="Hexagonal"
                                Car_Symbol ="h"
                                Car_System ="Hexagonal"
                                 else
                                          err_crys=.true.
                                err_mess_crys=" Error obtaining Crystal Familiy"
                                 end if
                          else
                   !---- Monoclinic ----!
                   Car_Family="Monoclinic"
                             Car_Symbol ="m"
                             Car_System ="Monoclinic"
                          end if

                 else
                          !---- Monoclic b-unique setting ----!
                          if (abs(cell%ang(1)-90.0) <= 0.0001) then
                                 Car_Family="Monoclinic"
                             Car_Symbol ="m"
                             Car_System ="Monoclinic"
                          else
                                 err_crys=.true.
                             err_mess_crys=" Error obtaining Crystal Familiy"
                          end if
                 end if

          case (3) ! all are the same angle
                 if (abs(cell%ang(1) - 90.000) <= 0.0001) then
                          select case (n2)
                                 case (1)
                                          !---- Orthorhombic ----!
                                          Car_Family="Orthorhombic"
                                          Car_Symbol ="o"
                                          Car_System ="Orthorhombic"

                                 case (2)
                                          !---- Tetragonal ----!
                                          if (icodp(1)==icodp(2)) then
                                                 Car_Family="Tetragonal"
                                             Car_Symbol ="t"
                                             Car_System ="Tetragonal"
                                          else
                                                 err_crys=.true.
                                   err_mess_crys=" Error obtaining Crystal Familiy"
                                          end if

                                 case (3)
                                          !---- Cubic ----!
                                          Car_Family="Cubic"
                                          Car_Symbol ="o"
                                          Car_System ="Cubic"

                          end select
                 else
                          if (n2 == 3) then
                                 !---- Hexagonal with rhombohedral axes ----!
                                 Car_Family="Hexagonal"
                                 Car_Symbol ="h"
                                 Car_System ="Trigonal"
                          else
                                 err_crys=.true.
                             err_mess_crys=" Error obtaining Crystal Familiy"
                          end if
                 end if

       end select ! n

       return
    End Subroutine Get_Cryst_Family

    !!--++
    !!--++ Subroutine Get_Cryst_Orthog_Matrix(Cellv,Ang, Crystort,Cartype)
    !!--++    real(kind=sp), dimension(3  ), intent (in ) :: cellv           !  In ->  a,b,c parameters
    !!--++    real(kind=sp), dimension(3  ), intent (in ) :: ang             !  In ->  angles parameters of cell
    !!--++    real(kind=sp), dimension(3,3), intent (out) :: CrystOrt        ! Out ->  Conversion matrix (a) = (e) CrystOrt
    !!--++    character (len=1), optional,   intent (in)  :: CarType         !  In ->  Type of Cartesian axes
    !!--++
    !!--++    (PRIVATE)
    !!--++    Obtains the matrix giving the crystallographic basis in
    !!--++    direct space in terms of a cartesian basis. The output matrix
    !!--++    can be directly used for transforming crystallographic components
    !!--++    to Cartesian components of the components of a vector considered
    !!--++    as a column vector:   XC = CrystOrt X.
    !!--++
    !!--++    If CartType is not present, or if it is not equal to 'A',
    !!--++    the cartesian system is defined as:
    !!--++          z // c; y is in the bc-plane; x is y ^ z
    !!--++    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
    !!--++    b = (         0         ,     b sinalpha      , b cosalpha)
    !!--++    c = (         0         ,         0           , c         )
    !!--++
    !!--++    If CartType = 'A', the the cartesian system is defined as:
    !!--++         x // a; y is in the ab-plane; z is x ^ z
    !!--++    a = (       a   ,         0           ,       0             )
    !!--++    b = ( b cosgamma,    b singamma       ,       0             )
    !!--++    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
    !!--++
    !!--++    The output matrix is the tranposed of the above one(s) so that the
    !!--++    matrix can directly be used for transforming "components" given
    !!--++    in a crystallographic basis to "components" in cartesian basis
    !!--++    when the components are used as "column" vectors.
    !!--++
    !!--++      [a] = C [e] , In [a],[e] basis vectors are in column form
    !!--++      (a) = (e) CT, In (a),(e) basis vectors are in row form
    !!--++      CrystOrt = CT  => (a) = (e) CystOrt
    !!--++
    !!--++    Remember that CT.C = C.CT = GD (direct cell metrics)
    !!--++
    !!--++
    !!--++      Xc = CrystOrt X (Xc Cartesian componets, X crystallographic components)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cryst_Orthog_Matrix(Cellv,Ang, Crystort,CarType)
       !---- Arguments ----!
       real(kind=sp), dimension(3  ), intent (in ) :: cellv,ang
       real(kind=sp), dimension(3,3), intent (out) :: CrystOrt
       character (len=1), optional,   intent (in ) :: CarType

       !---- Local Variables ----!
       real(kind=sp) :: cosgas, singas

       if (present(CarType)) then
          if (CarType == "A" .or. CarType == "a" ) then  ! x//a
             !  Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )
             cosgas =(cosd(ang(3))*cosd(ang(2))-cosd(ang(1)))/(sind(ang(3))*sind(ang(2)))
             singas = sqrt(1.0-cosgas**2)
             CrystOrt(1,1) = cellv(1)
             CrystOrt(1,2) = cellv(2)*cosd(ang(3))
             CrystOrt(1,3) = cellv(3)*cosd(ang(2))
             CrystOrt(2,1) = 0.0
             CrystOrt(2,2) = cellv(2)*sind(ang(3))
             CrystOrt(2,3) =-cellv(3)*sind(ang(2))*cosgas
             CrystOrt(3,1) = 0.0
             CrystOrt(3,2) = 0.0
             CrystOrt(3,3) = cellv(3)*sind(ang(2))*singas
             return
          end if
       end if

       !
       !  By default, the cartesian frame is such as z//c
       !  Transponse of the following matrix:
       !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
       !    b = (         0         ,     b sinalpha      , b cosalpha)
       !    c = (         0         ,         0           , c         )
       cosgas =(cosd(ang(1))*cosd(ang(2))-cosd(ang(3)))/(sind(ang(1))*sind(ang(2)))
       singas = sqrt(1.0-cosgas**2)
       CrystOrt(1,1) = cellv(1)*sind(ang(2))*singas
       CrystOrt(1,2) = 0.0
       CrystOrt(1,3) = 0.0
       CrystOrt(2,1) =-cellv(1)*sind(ang(2))*cosgas
       CrystOrt(2,2) = cellv(2)*sind(ang(1))
       CrystOrt(2,3) = 0.0
       CrystOrt(3,1) = cellv(1)*cosd(ang(2))
       CrystOrt(3,2) = cellv(2)*cosd(ang(1))
       CrystOrt(3,3) = cellv(3)

       return
    End Subroutine Get_Cryst_Orthog_Matrix

    !!----
    !!---- Subroutine Get_Deriv_Orth_Cell(Cellp,De_Orthcell,Cartype)
    !!----    type(Crystal_Cell_type),         intent(in ) :: cellp
    !!----    real(kind=sp), dimension(3,3,6), intent(out) :: de_Orthcell
    !!----    character (len=1), optional,    intent (in ) :: CarType
    !!----
    !!----    Subroutine to get derivative matrix of the transformation matrix
    !!----    to orthogonal frame. Useful for calculations of standard deviations
    !!----    of distances and angles. The specialised subroutine calculating
    !!----    sigmas of distances "distance_and_sigma" is in Atom_mod.
    !!----    The output matrices "de_Orthcell" are the derivatives of, with
    !!----    respect to a(1),b(2),c(3),alpha(4),beta(5) and gamma(6) of the
    !!----    matrix   "Cellp%Cr_Orth_cel".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Deriv_Orth_Cell(Cellp,De_Orthcell,Cartype)
       !---- Arguments ----!
       type(Crystal_Cell_type),         intent(in ) :: cellp
       real(kind=sp), dimension(3,3,6), intent(out) :: de_Orthcell
       character (len=1), optional,    intent (in ) :: CarType

       !---- Local Variables ----!
       real(kind=sp) ::  ca,cb,cg,sa,sb,sg,f,g, fa,fb,fc,ga,gb,gc

       de_Orthcell=0.0
       ca=cosd(cellp%ang(1))
       cb=cosd(cellp%ang(2))
       cg=cosd(cellp%ang(3))
       sa=sind(cellp%ang(1))
       sb=sind(cellp%ang(2))
       sg=sind(cellp%ang(3))

       if (present(CarType)) then
          if (CarType == "A" .or. CarType == "a" ) then  ! x//a

             f=(ca-cb*cg)/sg    !-cosgas*sinbeta
             g=SQRT(sb*sb-f*f)  ! singas*sinbeta
             fa=-sa/sg          ! df/dalpha
             fb=sb*cg/sg        ! df/dbeta
             fc=cb/sg**2        ! df/dgamma
             ga=-f*fa/g         ! dg/dalpha
             gb=(sb*cb-f*fb)/g  ! dg/dbeta
             gc=f/g*fc          ! dg/dgamma

             ! M: Transponse of the following matrix:
             !    a = (       a   ,         0           ,       0             )
             !    b = ( b cosgamma,    b singamma       ,       0             )
             !    c = (  c cosbeta, -c sinbeta cosalpha*, c sinbeta sinalpha* )

             !
             !        (   a         b*cg        c*cb )
             !    M = (   0         b*sg        c*f  )
             !        (   0          0          c*g  )
             !
             !           (   1      0      0 )
             !  dM_da =  (   0      0      0 )
             !           (   0      0      0 )
             de_Orthcell(1,1,1) = 1.0
             !           (   0      cg     0 )
             !  dM_db =  (   0      sg     0 )
             !           (   0      0      0 )
             de_Orthcell(1,2,2) = cg
             de_Orthcell(2,2,2) = sg
             !
             !            (   0          0          cb )
             !  dM_dc =   (   0          0          f  )
             !            (   0          0          g  )
             de_Orthcell(1,3,3) = cb
             de_Orthcell(2,3,3) = f
             de_Orthcell(3,3,3) = g
             !
             !             (   0          0           0   )
             ! dM_dalpha=  (   0          0          c*fa )
             !             (   0          0          c*ga )
             !
             de_Orthcell(2,3,4) = cellp%cell(3)*fa
             de_Orthcell(3,3,4) = cellp%cell(3)*ga
             !
             !             (   0          0         -c*sb )
             ! dM_dbeta =  (   0          0          c*fb )
             !             (   0          0          c*gb )
             !
             de_Orthcell(1,3,5) = -cellp%cell(3)*sb
             de_Orthcell(2,3,5) =  cellp%cell(3)*fb
             de_Orthcell(3,3,5) =  cellp%cell(3)*gb
             !
             !              (   0        -b*sg         0   )
             ! dM_dgamma =  (   0         b*cg        c*fc )
             !              (   0          0          c*gc )
             !
             de_Orthcell(1,2,6) = -cellp%cell(2)*sg
             de_Orthcell(2,2,6) =  cellp%cell(2)*cg
             de_Orthcell(2,3,6) =  cellp%cell(3)*fc
             de_Orthcell(3,3,6) =  cellp%cell(3)*gc
             return
          end if
       end if

       !
       !  By default, the cartesian frame is such as z//c
       !  Transponse of the following matrix:
       !    a = (a sinbeta singamma*, -a sinbeta cosgamma*, a cosbeta )
       !    b = (         0         ,     b sinalpha      , b cosalpha)
       !    c = (         0         ,         0           , c         )

       !         ( a sinbeta singamma*          0             0 )
       !    M =  (-a sinbeta cosgamma*      b sinalpha        0 )
       !         ( a cosbeta                b cosalpha        c )

       f=(cg-ca*cb)/sa    !-sinbeta . cosgamma*
       g=SQRT(sb*sb-f*f)  ! sinbeta . singamma*
       fa= cb/sa**2       ! df/dalpha
       fb=sb*ca/sa        ! df/dbeta
       fc=-sb/sa          ! df/dgamma
       ga=-f*fa/g         ! dg/dalpha
       gb=(sb*cb-f*fb)/g  ! dg/dbeta
       gc=f/g*fc          ! dg/dgamma

       !         ( a*g        0         0 )
       !    M =  ( a*f      b*sa        0 )
       !         ( a*cb     b*ca        c )

       !
       !           (   g       0      0 )
       !  dM_da =  (   f       0      0 )
       !           (   cb      0      0 )
       de_Orthcell(1,1,1) = g
       de_Orthcell(1,2,1) = f
       de_Orthcell(1,3,1) = cb

       !           (   0      0      0 )
       !  dM_db =  (   0      sa     0 )
       !           (   0      ca     0 )
       de_Orthcell(1,2,2) = sa
       de_Orthcell(3,2,2) = ca
       !
       !            (   0      0      0  )
       !  dM_dc =   (   0      0      0  )
       !            (   0      0      1  )
       de_Orthcell(3,3,3) = 1
       !
       !             ( a*ga         0          0 )
       ! dM_dalpha=  ( a*fa       -b*ca        0 )
       !             (   0         b*sa        0 )
       !
       de_Orthcell(1,1,4) = cellp%cell(1)*ga
       de_Orthcell(2,1,4) = cellp%cell(1)*fa
       de_Orthcell(2,2,4) =-cellp%cell(2)*ca
       de_Orthcell(3,2,4) = cellp%cell(2)*sa
       !
       !             (  a*gb        0         0 )
       ! dM_dbeta =  (  a*fb        0         0 )
       !             ( -a*sb        0         0 )
       !
       de_Orthcell(1,1,5) = cellp%cell(1)*gb
       de_Orthcell(2,1,5) = cellp%cell(1)*fb
       de_Orthcell(3,1,5) =-cellp%cell(1)*sb
       !
       !              (  a*gc     0      0   )
       ! dM_dgamma =  (  a*fc     0      0   )
       !              (   0       0      0   )
       !
       de_Orthcell(1,1,6) = cellp%cell(1)*gc
       de_Orthcell(2,1,6) = cellp%cell(1)*fc

       return
    End Subroutine Get_Deriv_Orth_Cell

    !!----
    !!---- SUBROUTINE INIT_ERR_CRYS()
    !!----
    !!----    Initialize Flags of Errors in this module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Crys()

       err_crys=.false.
       err_mess_crys=" "

       return
    End Subroutine Init_Err_Crys

    !!--++
    !!--++ Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
    !!--++    real(kind=sp), dimension(3), intent(in ) :: a        !  In -> a,b,c
    !!--++    real(kind=sp), dimension(3), intent(in ) :: ang      !  In -> alpha,beta,gamma
    !!--++    real(kind=sp), dimension(3), intent(out) :: ar       !  In -> a*,b*,c*
    !!--++    real(kind=sp), dimension(3), intent(out) :: angr     !  In -> alpha*,beta*,gamma*
    !!--++    real(kind=sp),               intent(out) :: vol      ! Out -> Vol
    !!--++    real(kind=sp),               intent(out) :: volr     ! Out -> Vol*
    !!--++
    !!--++    (PRIVATE)
    !!--++    Calculates the reciprocal lattice vectors and cell volume
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Recip(A,Ang,Ar,Angr,Vol,Volr)
       !---- Arguments ----!
       real(kind=sp), dimension(3), intent(in ) :: a,ang
       real(kind=sp), dimension(3), intent(out) :: ar,angr
       real(kind=sp),               intent(out) :: vol,volr

       !---- Local Variables ----!
       integer        :: i
       real(kind=sp)  :: s,p,cose

       p=1.0
       s=1.0
       do i=1,3
          cose=cosd(ang(i))
          p=p*cose
          s=s-cose*cose
       end do
       vol=sqrt(abs(s+2.0*p))

       do i=1,3
          vol=vol*a(i)
       end do
       volr=1.0/vol

       ar(1)=a(2)*a(3)*sind(ang(1))/vol
       ar(2)=a(3)*a(1)*sind(ang(2))/vol
       ar(3)=a(1)*a(2)*sind(ang(3))/vol
       angr(1)=(cosd(ang(2))*cosd(ang(3))-cosd(ang(1)))/(sind(ang(2))*sind(ang(3)))
       angr(2)=(cosd(ang(1))*cosd(ang(3))-cosd(ang(2)))/(sind(ang(1))*sind(ang(3)))
       angr(3)=(cosd(ang(2))*cosd(ang(1))-cosd(ang(3)))/(sind(ang(2))*sind(ang(1)))
       do i=1,3
          angr(i)=acosd(angr(i))
       end do

       return
    End Subroutine Recip

    !!----
    !!---- Subroutine Set_Crystal_Cell(Cellv,Angl,Celda,Cartype,Scell,Sangl)
    !!----    real(kind=sp), dimension (3),        intent(in ) :: cellv   !  In -> a,b,c
    !!----    real(kind=sp), dimension (3),        intent(in ) :: angl    !  In -> angles of cell parameters
    !!----    Type (Crystal_Cell_Type),            intent(out) :: Celda   !  Out-> Celda components
    !!----    character (len=1),          optional,intent(in ) :: CarType !  In -> Type of Cartesian Frame
    !!----    real(kind=sp), dimension(3),optional,intent(in ) :: scell,sangl
    !!----
    !!----    Constructs the object "Celda" of type Crystal_Cell. Control for error
    !!----    is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Crystal_Cell(cellv,angl,Celda,CarType,scell,sangl)
       !---- Arguments ----!
       real(kind=sp), dimension (3),        intent(in ) :: cellv, angl
       Type (Crystal_Cell_Type),            intent(out) :: Celda
       character (len=1),          optional,intent(in ) :: CarType
       real(kind=sp), dimension(3),optional,intent(in ) :: scell,sangl

       !---- Local Variables ----!
       integer :: ifail

       call init_err_crys()

       if (present(scell) .and. present(sangl)) then
         Celda%cell_std=scell
         Celda%ang_std=sangl
       else
         Celda%cell_std=0.0
         Celda%ang_std=0.0
       end if

       Celda%cell=cellv
       Celda%ang=angl
       where(Celda%ang < eps) Celda%ang =90.0
       call recip(cellv,angl,Celda%rcell,Celda%rang,Celda%CellVol,Celda%RCellVol)
       if (present(CarType)) then
          call Get_Cryst_Orthog_matrix(cellv,angl,Celda%Cr_Orth_cel,CarType)
          Celda%CartType=CarType
       else
          call Get_Cryst_Orthog_matrix(cellv,angl,Celda%Cr_Orth_cel)
          Celda%CartType="C"
       end if
       call matrix_inverse(Celda%Cr_Orth_cel,Celda%Orth_Cr_cel,ifail)

       if (ifail /= 0) then
          err_crys=.true.
          err_mess_crys=" Bad cell parameters "
          return
       end if

       Celda%GD=Metrics(cellv,angl)
       Celda%GR=Metrics(Celda%rcell,Celda%rang)

       return
    End Subroutine Set_Crystal_Cell

    !!----
    !!---- Subroutine Write_Crystal_Cell(Celda,Lun)
    !!----    Type (Crystal_Cell_Type),  intent(in)  :: Celda   !  In -> Cell variable
    !!----    Integer,optional           intent(in)  :: lun     !  In -> Unit to write
    !!----
    !!----    Writes the cell characteristics in a file associated to the
    !!----    logical unit lun
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Crystal_Cell(Celda,Lun)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),  intent(in) :: Celda
       Integer,intent(in) ,  optional        :: lun

       !---- Local variables ----!
       integer            :: iunit
       integer            :: i,j

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")    "        Metric information:"
       Write(unit=iunit,fmt="(a,/)")    "        -------------------"
       Write(unit=iunit,fmt="(a,/)")    " => Direct cell parameters:"
       Write(unit=iunit,fmt="(3(a,f12.4))")"         a = ", Celda%cell(1),"      b = ", Celda%cell(2), "      c = ", Celda%cell(3)
       Write(unit=iunit,fmt="(3(a,f12.3))")"     alpha = ", Celda%ang(1) ,"   beta = ", Celda%ang(2) , "  gamma = ", Celda%ang(3)
       Write(unit=iunit,fmt="(a,f12.4)")   "                        Direct Cell Volume = ",Celda%CellVol
       Write(unit=iunit,fmt="(/,a,/)")     " => Reciprocal cell parameters:"
       Write(unit=iunit,fmt="(3(a,f12.6))")"         a*= ", Celda%rcell(1),"      b*= ",Celda%rcell(2),"      c*= ", Celda%rcell(3)
       Write(unit=iunit,fmt="(3(a,f12.3))")"     alpha*= ", Celda%rang(1) ,"   beta*= ",Celda%rang(2) ,"  gamma*= ", Celda%rang(3)
       Write(unit=iunit,fmt="(a,f12.8)")   "                    Reciprocal Cell Volume = ",Celda%RCellVol
       Write(unit=iunit,fmt="(/,a,/)")     " => Direct and Reciprocal Metric Tensors:"
       Write(unit=iunit,fmt="(a)")         "                   GD                                       GR"

       do i=1,3
          Write(unit=iunit,fmt="(3f12.4,a,3f12.6)") (Celda%GD(i,j),j=1,3),"      ", (Celda%GR(i,j),j=1,3)
       end do
       If(Celda%CartType == "A") then
         Write(unit=iunit,fmt="(/,a,/)") " =>  Cartesian frame: x // a; y is in the ab-plane; z is x ^ y   "
       else
         Write(unit=iunit,fmt="(/,a,/)") " =>  Cartesian frame: z // c; y is in the bc-plane; x is y ^ z   "
       end if
       Write(unit=iunit,fmt="(a)")       "     Crystal_to_Orthonormal_Matrix              Orthonormal_to_Crystal Matrix"
       Write(unit=iunit,fmt="(a)")       "              Cr_Orth_cel                               Orth_Cr_cel  "
       do i=1,3
          Write(unit=iunit,fmt="(3f12.4,a,3f12.6)") (Celda%Cr_Orth_cel(i,j),j=1,3),"      ", (Celda%Orth_Cr_cel(i,j),j=1,3)
       end do

       return
    End Subroutine Write_Crystal_Cell

 End Module Crystal_Types
