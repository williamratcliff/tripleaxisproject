!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: MATH_3D
!!----   INFO: Simple mathematics general utilities for 3D Systems
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----            November - 2000  Updated by JGP
!!----            October  - 1996  Routines created by JRC
!!----
!!---- DEPENDENCIES
!!--++    Use MATH_GEN,   only: sp, dp, pi, cosd, sind, to_rad, to_deg
!!----
!!---- VARIABLES
!!--++    EPS                          [Private]
!!----    ERR_MATH_3D
!!----    ERR_MESS_MATH_3D
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CROSS_PRODUCT
!!--++       CROSS_PRODUCT_dp          [Overloaded]
!!--++       CROSS_PRODUCT_sp          [Overloaded]
!!----       DETERM_A
!!--++       DETERM_A_I                [Overloaded]
!!--++       DETERM_A_R                [Overloaded]
!!----       DETERM_V
!!--++       DETERM_V_I                [Overloaded]
!!--++       DETERM_V_R                [Overloaded]
!!----       INVERT_A
!!--++       INVERT_DP                 [Overloaded]
!!--++       INVERT_SP                 [Overloaded]
!!----       ROTATE_OX
!!----       ROTATE_OY
!!----       ROTATE_OZ
!!----       VECLENGTH
!!----
!!----    Subroutines:
!--..
!!--..    Init Routine
!!----       INIT_ERR_MATH3D
!!----       SET_EPS
!!----       SET_EPS_DEFAULT
!--..
!!--..    Matrix and Vectors Subroutines
!!----       GET_CART_FROM_CYLIN
!!--++       GET_CART_FROM_CYLIN_DP    [Overloaded]
!!--++       GET_CART_FROM_CYLIN_SP    [Overloaded]
!!----       GET_CYLINDR_COORD
!!--++       GET_CYLINDR_COORD_DP      [Overloaded]
!!--++       GET_CYLINDR_COORD_SP      [Overloaded]
!!----       GET_CART_FROM_SPHER
!!--++       GET_CART_FROM_SPHER_DP    [Overloaded]
!!--++       GET_CART_FROM_SPHER_SP    [Overloaded]
!!----       GET_PLANE_FROM_POINTS
!!----       GET_SPHERIC_COORD
!!--++       GET_SPHERIC_COORD_DP      [Overloaded]
!!--++       GET_SPHERIC_COORD_SP      [Overloaded]
!!----       MATRIX_DIAGEIGEN
!!----       MATRIX_INVERSE
!!----       RESOLV_SIST_1X2
!!----       RESOLV_SIST_1X3
!!----       RESOLV_SIST_2X2
!!----       RESOLV_SIST_2X3
!!----       RESOLV_SIST_3X3
!!----
!!
 Module Math_3D

    !---- Use Modules ----!
    Use Math_gen, only: sp, dp, pi, cosd, sind, to_rad, to_deg

    implicit none

    private

    !---- List of public functions ----!
    public :: Rotate_OX, Rotate_OY, Rotate_OZ, Veclength

    !---- List of public overloaded procedures: functions ----!
    public :: Cross_Product, Determ_A, Determ_V, Invert_A

    !---- List of public subroutines ----!
    public :: Init_Err_Math3D, Set_Eps, Set_Eps_Default, Matrix_DiagEigen, Matrix_Inverse, &
              Resolv_Sist_1X2, Resolv_Sist_1X3, Resolv_Sist_2X2, Resolv_Sist_2X3,          &
              Resolv_Sist_3X3, Get_Plane_from_Points

    !---- List of public overloaded procedures: subroutines ----!
    public :: Get_Cart_From_Cylin, Get_Cylindr_Coord, Get_Cart_From_Spher, Get_Spheric_Coord

    !----  Make private the overloaded procedures ----!
    private :: Cross_Product_dp, Cross_Product_sp, Determ_A_I, Determ_A_R, Determ_V_I,  &
               Determ_V_R, Invert_dp, Invert_sp, Get_Cart_From_Cylin_dp,                &
               Get_Cart_From_Cylin_sp, Get_Cylindr_Coord_dp, Get_Cylindr_Coord_sp,      &
               Get_Cart_From_Spher_dp, Get_Cart_From_Spher_sp, Get_Spheric_Coord_dp,    &
               Get_Spheric_Coord_sp

    !---- Definitions ----!

    !!--++
    !!--++ EPS
    !!--++    real(kind=sp), private :: eps
    !!--++
    !!--++    (PRIVATE)
    !!--++    epsilon value may be changed by calling to the public procedure
    !!--++    Set_Eps
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), private ::  eps = 0.00001

    !!----
    !!---- ERR_MATH_3D
    !!----    logical :: err_math_3d
    !!----
    !!----    Logical Variable indicating an error in MATH_3D module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: err_math_3d

    !!----
    !!---- ERR_MESS_MATH_3D
    !!----    character(len=150) :: err_mess_math_3d
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_math_3d

    !---- Interfaces - Overlapp ----!
    Interface  Cross_Product
       Module Procedure Cross_product_sp
       Module Procedure Cross_product_dp
    End Interface

    Interface  Determ_A
       Module Procedure Determ_A_I
       Module Procedure Determ_A_R
    End Interface

    Interface  Determ_V
       Module Procedure Determ_V_I
       Module Procedure Determ_V_R
    End Interface

    Interface  Invert_A
       Module Procedure Invert_sp
       Module Procedure Invert_dp
    End Interface

    Interface  Get_Cart_from_Cylin
       Module Procedure Get_Cart_from_Cylin_dp
       Module Procedure Get_Cart_from_Cylin_sp
    End Interface

    Interface  Get_Cylindr_Coord
       Module Procedure Get_Cylindr_Coord_dp
       Module Procedure Get_Cylindr_Coord_sp
    End Interface

    Interface  Get_Cart_from_Spher
       Module Procedure Get_Cart_from_Spher_dp
       Module Procedure Get_Cart_from_Spher_sp
    End Interface

    Interface  Get_Spheric_Coord
       Module Procedure Get_Spheric_Coord_dp
       Module Procedure Get_Spheric_Coord_sp
    End Interface

 Contains

    !!----
    !!---- Function  Cross_Product(U,V) Result(W)
    !!----    real(kind=sp/dp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!----    real(kind=sp/dp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!----    real(kind=sp/dp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!----
    !!----    Calculates the cross product of vectors u and v
    !!----    Vectors, w= u x v, are given in cartesian components.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function  Cross_Product_dp(U,V) Result(W)
    !!--++    real(kind=dp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    real(kind=dp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    real(kind=dp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_dp(u,v) Result(w)
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: u,v
       real(kind=dp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)
       w(2)=u(3)*v(1)-u(1)*v(3)
       w(3)=u(1)*v(2)-u(2)*v(1)

       return
    End Function Cross_Product_dp

    !!--++
    !!--++ Function  Cross_Product_sp(U,V) Result(W)
    !!--++    real(kind=sp), dimension(3), intent( in) :: u   !  In -> Vector 1
    !!--++    real(kind=sp), dimension(3), intent( in) :: v   !  In -> Vector 2
    !!--++    real(kind=sp), dimension(3)              :: w   ! Out -> Vector 1 x vector 2
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the cross product of vectors u and v
    !!--++    Vectors, w= u x v, are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Cross_Product_sp(u,v) Result(w)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: u,v
       real(kind=sp), dimension(3)              :: w

       w(1)=u(2)*v(3)-u(3)*v(2)  ! i  j   k !
       w(2)=u(3)*v(1)-u(1)*v(3)  !u1  u2  u3! = (u2.v3 - u3.v2)i + (v1.u3 - u1.v3)j + (u1.v2-u2.v1)k
       w(3)=u(1)*v(2)-u(2)*v(1)  !v1  v2  v3!

       return
    End Function Cross_Product_sp

    !!----
    !!---- Function Determ_A(A)
    !!----    integer/real(kind=sp), dimension(3,3), intent(in)  :: a
    !!----
    !!----    Calculates the determinant of an integer/real 3x3 matrix
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Determ_A_I(A)
    !!--++    integer, dimension(3,3), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of an integer 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_I(A) Result(determ)
       !---- Argument ----!
       integer, dimension(3,3), intent(in) :: A
       integer                             :: determ

       determ=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(1,2)*A(2,3)*A(3,1) &
             -A(1,3)*A(2,2)*A(3,1)-A(1,1)*A(3,2)*A(2,3)-A(1,2)*A(2,1)*A(3,3)

       return
    End Function Determ_A_I

    !!--++
    !!--++ Function Determ_A_R(A)
    !!--++    real(kind=sp), dimension(3,3), intent(in)  :: a
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_A_R(A) Result (determ)
       !---- Argument ----!
       real(kind=sp), dimension(3,3), intent(in) :: A
       real(kind=sp)                             :: determ

       determ=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(1,2)*A(2,3)*A(3,1) &
             -A(1,3)*A(2,2)*A(3,1)-A(1,1)*A(3,2)*A(2,3)-A(1,2)*A(2,1)*A(3,3)

       return
    End Function Determ_A_R

    !!----
    !!---- Function  Determ_V(a,b,c)
    !!----    integer/real(kind=sp), dimension(3), intent(in) :: a,b,c
    !!----
    !!----    Calculates the determinant of the components of three vectors
    !!----
    !!----  Update: February - 2005
    !!

    !!--++
    !!--++ Function Determ_V_I(A,B,C)
    !!--++    integer, dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_I(a,b,c) Result(det)
       !---- Arguments ----!
       integer, dimension(3), intent(in) :: a,b,c
       integer                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+a(i)*(b(j)*c(k)-b(k)*c(j))
       end do

       return
    End Function Determ_V_I

    !!--++
    !!--++ Function Determ_V_R(A,B,C)
    !!--++    real(kin=sp), dimension(3), intent(in) :: a,b,c
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of the components of three vectors
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Determ_V_R(a,b,c) Result(det)
       !---- Arguments ----!
       real(kind=sp), dimension(3), intent(in) :: a,b,c
       real(kind=sp)                           :: det

       !---- Local variables ----!
       integer :: i,j,k

       det = 0.0
       do i = 1,3
          j = i+1
          if (j == 4) j = 1
          k = 6-i-j
          det = det+a(i)*(b(j)*c(k)-b(k)*c(j))
       end do

       return
    End Function Determ_V_R

    !!----
    !!---- Funcion Invert_A(A) Result(b)
    !!----    real(kind=sp/dp), dimension(3,3), intent(in) :: a
    !!----    real(Kind=sp/dp), dimension(3,3)             :: b
    !!----
    !!----    Calculate de inverse of a real 3x3 matrix. If the routine fails,
    !!----    then a 0.0 matrix is returned.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Funcion Invert_Dp(A) Result(b)
    !!--++    real(kind=dp), dimension(3,3), intent(in) :: a
    !!--++    real(Kind=dp), dimension(3,3)             :: b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Dp(a) Result(b)
       !---- Arguments ----!
       real(kind=dp),dimension(3,3), intent(in) :: a
       real(kind=dp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=dp)  :: dmat

       b(1,1) =   a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) =   a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) =   a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) =   a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) =   a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1) !determinant of A

       if (abs(dmat) > tiny(dmat)) then
          b= b/dmat
       else
          b=0.0_dp
       end if

       return
    End Function Invert_Dp

    !!--++
    !!--++ Funcion Invert_Sp(A) Result(b)
    !!--++    real(kind=sp), dimension(3,3), intent(in) :: a
    !!--++    real(Kind=sp), dimension(3,3)             :: b
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate de inverse of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Invert_Sp(a) Result(b)
       !---- Arguments ----!
       real(kind=sp),dimension(3,3), intent(in) :: a
       real(kind=sp),dimension(3,3)             :: b

       !---- Local variables ----!
       real(kind=sp)  :: dmat

       b(1,1) =   a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) =   a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) =   a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) =   a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) =   a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1) !determinant of A

       if (abs(dmat) > tiny(dmat)) then
          b= b/dmat
       else
          b=0.0
       end if

       return
    End Function Invert_Sp

    !!----
    !!---- Function Rotate_OX(X,Angle) Result (Vec)
    !!----    real, dimension(3), intent(in) :: x       !  In -> Vector
    !!----    real,               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real, dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    X Rotation. Positive rotation is unclockwise
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OX(X,Angle) Result(vec)
       !---- Arguments ----!
       real, dimension(3), intent(in) :: x
       real,               intent(in) :: angle
       real, dimension(3)             :: vec

       !---- Variables locales ----!
       real, dimension(3,3)           :: rot

       rot(1,1)=  1.0
       rot(2,1)=  0.0
       rot(3,1)=  0.0

       rot(1,2)=  0.0
       rot(2,2)=  cosd(angle)
       rot(3,2)=  sind(angle)

       rot(1,3)=  0.0
       rot(2,3)=  -sind(angle)
       rot(3,3)=  cosd(angle)

       vec=matmul(rot,x)

       return
    End Function Rotate_OX

    !!----
    !!---- Function Rotate_OY(Y,Angle) Result (Vec)
    !!----    real, dimension(3), intent(in) :: y       !  In -> Vector
    !!----    real,               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real, dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    Y Rotation.
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OY(Y,Angle) Result(vec)
       !---- Arguments ----!
       real, dimension(3), intent(in) :: y
       real,               intent(in) :: angle     ! Angulo en grados
       real, dimension(3)             :: vec

       !---- Variables locales ----!
       real, dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  0.0
       rot(3,1)=  -sind(angle)

       rot(1,2)=  0.0
       rot(2,2)=  1.0
       rot(3,2)=  0.0

       rot(1,3)= sind(angle)
       rot(2,3)= 0.0
       rot(3,3)= cosd(angle)

      vec=matmul(rot,y)

       return
    End Function Rotate_OY

    !!----
    !!---- Function Rotate_OZ(Z,Angle) Result (Vec)
    !!----    real, dimension(3), intent(in) :: z       !  In -> Vector
    !!----    real,               intent(in) :: angle   !  In -> Angle (Degrees)
    !!----    real, dimension(3)             :: vec     ! Out -> Vector
    !!----
    !!----    Z Rotation
    !!----
    !!---- Update: February - 2005
    !!
    Function Rotate_OZ(Z,Angle) Result(vec)
       !---- Arguments ----!
       real, dimension(3), intent(in) :: z
       real,               intent(in) :: angle
       real, dimension(3)             :: vec

       !---- Variables locales ----!
       real, dimension(3,3)           :: rot

       rot(1,1)=  cosd(angle)
       rot(2,1)=  sind(angle)
       rot(3,1)=  0.0

       rot(1,2)=  -sind(angle)
       rot(2,2)=  cosd(angle)
       rot(3,2)=  0.0

       rot(1,3)=  0.0
       rot(2,3)=  0.0
       rot(3,3)=  1.0

       vec=matmul(rot,z)

       return
    End Function Rotate_OZ

    !!----
    !!---- Function Veclength(A,B) Result(c)
    !!----    real(kind=sp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=sp), dimension(3),   intent(in)  :: b
    !!----    real(kind=sp),                             :: c
    !!----
    !!----    Length of vector B when A is the Crystallographic
    !!----    to orthogonal matrix length=c
    !!----
    !!---- Update: February - 2005
    !!
    Function Veclength(a,b) Result(c)
       !---- Arguments ----!
       real(kind=sp), intent(in)  , dimension(3,3)       :: a
       real(kind=sp), intent(in)  , dimension(3  )       :: b
       real(kind=sp)                                     :: c

       !---- Local variables ----!
       integer                     :: i,j
       real(kind=sp), dimension(3) :: v

       v=0.0
       do i = 1,3
          do j = 1,3
             v(i) = v(i)+a(i,j)*b(j)
          end do
       end do

       c = sqrt(v(1)**2+v(2)**2+v(3)**2)

       return
    End Function Veclength

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Math3D()
    !!----
    !!----    Initialize the errors flags in Math_3d
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Math3D()

       err_math_3d=.false.
       err_mess_math_3d=" "

       return
    End Subroutine Init_Err_Math3D

    !!----
    !!---- Subroutine Set_Eps(Neweps)
    !!----    real(kind=sp), intent( in) :: neweps
    !!----
    !!----    Sets global EPS to the value "neweps"
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Eps(Neweps)
       !---- Arguments ----!
       real(kind=sp), intent( in) :: neweps

       eps=neweps

       return
    End Subroutine Set_Eps

    !!----
    !!---- Subroutine Set_Eps_Default()
    !!----
    !!----    Sets global EPS to the default value: eps=0.00001
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Eps_Default()

       eps=0.00001

       return
    End Subroutine Set_Eps_Default

    !!----
    !!---- Subroutine Get_Cart_from_Cylin(rho,Phi,zeta,Xo,Mode)
    !!----    real(kind=sp/dp),              intent( in)           :: rho
    !!----    real(kind=sp/dp),              intent( in)           :: phi
    !!----    real(kind=sp/dp),              intent( in)           :: zeta
    !!----    real(kind=sp/dp), dimension(3),intent(out)           :: xo
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the Cartesian coordinates from cylindrical coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_dp(rho,Phi,zeta,Xo,Mode)
    !!--++    real(kind=dp),              intent( in)           ::  rho
    !!--++    real(kind=dp),              intent( in)           ::  phi
    !!--++    real(kind=dp),              intent( in)           ::  zeta
    !!--++    real(kind=dp), dimension(3),intent(out)           ::  xo
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine  Get_Cart_from_Cylin_dp(rho,Phi,zeta,Xo,Mode)
       !---- Arguments ----!
       real(kind=dp),              intent( in)           ::  rho
       real(kind=dp),              intent( in)           ::  phi
       real(kind=dp),              intent( in)           ::  zeta
       real(kind=dp), dimension(3),intent(out)           ::  xo
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       real(kind=dp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=zeta

       return
    End Subroutine Get_Cart_from_Cylin_dp

    !!--++
    !!--++ Subroutine  Get_Cart_from_Cylin_sp(rho,Phi,zeta,Xo,Mode)
    !!--++    real(kind=sp),              intent( in)           ::  rho
    !!--++    real(kind=sp),              intent( in)           ::  phi
    !!--++    real(kind=sp),              intent( in)           ::  zeta
    !!--++    real(kind=sp), dimension(3),intent(out)           ::  xo
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from cylindrical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine  Get_Cart_from_Cylin_sp(rho,Phi,zeta,Xo,Mode)
       real(kind=sp),              intent( in)           ::  rho
       real(kind=sp),              intent( in)           ::  phi
       real(kind=sp),              intent( in)           ::  zeta
       real(kind=sp), dimension(3),intent(out)           ::  xo
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       real(kind=sp) :: ph

       ph=phi
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") ph=phi*to_rad
       end if
       xo(1)=rho*cos(ph)
       xo(2)=rho*sin(ph)
       xo(3)=zeta

       return
    End Subroutine Get_Cart_from_Cylin_sp

    !!----
    !!---- Subroutine Get_Cylindr_Coord(Xo,rho,Phi,zeta,Mode)
    !!----    real(kind=sp/dp), dimension(3),intent( in)           :: xo
    !!----    real(kind=sp/dp),              intent(out)           :: rho
    !!----    real(kind=sp/dp),              intent(out)           :: phi
    !!----    real(kind=sp/dp),              intent(out)           :: zeta
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the cylindrical coordinates from Cartesian coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine  Get_Cylindr_Coord_dp(Xo,rho,Phi,zeta,Mode)
    !!--++    real(kind=dp), dimension(3),intent( in)           ::  xo
    !!--++    real(kind=dp),              intent(out)           ::  rho
    !!--++    real(kind=dp),              intent(out)           ::  phi
    !!--++    real(kind=dp),              intent(out)           ::  zeta
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine  Get_Cylindr_Coord_dp(Xo,rho,Phi,zeta,Mode)
       !---- Arguments ----!
       real(kind=dp), dimension(3),intent( in)           ::  xo
       real(kind=dp),              intent(out)           ::  rho
       real(kind=dp),              intent(out)           ::  phi
       real(kind=dp),              intent(out)           ::  zeta
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       zeta=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_dp
       end if
       rho=0.0_dp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylindr_Coord_dp

    !!--++
    !!--++ Subroutine  Get_Cylindr_Coord_sp(Xo,rho,Phi,zeta,Mode)
    !!--++    real(kind=sp), dimension(3),intent( in)           ::  xo
    !!--++    real(kind=sp),              intent(out)           ::  rho
    !!--++    real(kind=sp),              intent(out)           ::  phi
    !!--++    real(kind=sp),              intent(out)           ::  zeta
    !!--++    character(len=*),           intent( in), optional ::  mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the cylindrical coordinates from Cartesian coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine  Get_Cylindr_Coord_sp(Xo,rho,Phi,zeta,Mode)
       !---- Arguments ----!
       real(kind=sp), dimension(3),intent( in)           ::  xo
       real(kind=sp),              intent(out)           ::  rho
       real(kind=sp),              intent(out)           ::  phi
       real(kind=sp),              intent(out)           ::  zeta
       character(len=*),           intent( in), optional ::  mode

       !---- Local Variables ----!
       integer :: j

       zeta=xo(3)
       if( abs(xo(2)) > eps .or. abs(xo(1)) > eps) then
          phi=atan2(xo(2),xo(1))
       else
          phi= 0.0_sp
       end if
       rho=0.0_sp
       do j=1,2
          rho=rho+xo(j)*xo(j)
       end do
       rho=sqrt(rho)

       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") phi=phi*to_deg
       end if

       return
    End Subroutine Get_Cylindr_Coord_sp

    !!----
    !!---- Subroutine Get_Cart_from_Spher(r,Theta,Phi,Xo,Mode)
    !!----    real(kind=sp/dp),              intent( in)           :: r
    !!----    real(kind=sp/dp),              intent( in)           :: Theta
    !!----    real(kind=sp/dp),              intent( in)           :: Phi
    !!----    real(kind=sp/dp), dimension(3),intent(out)           :: xo
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the Cartesian coordinates from spherical coordinates.
    !!----    If Mode='D' the angle phi is provided in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
    !!--++    real(kind=dp),              intent( in)           :: r
    !!--++    real(kind=dp),              intent( in)           :: Theta
    !!--++    real(kind=dp),              intent( in)           :: Phi
    !!--++    real(kind=dp), dimension(3),intent(out)           :: xo
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_dp(r,Theta,Phi,Xo,Mode)
       !---- Arguments ----!
       real(kind=dp),              intent( in)           :: r
       real(kind=dp),              intent( in)           :: Theta
       real(kind=dp),              intent( in)           :: phi
       real(kind=dp), dimension(3),intent(out)           :: xo
       character(len=*),           intent( in), optional :: mode

       !---- Local Variables ----!
       real(kind=dp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_dp

    !!--++
    !!--++ Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
    !!--++    real(kind=sp),              intent( in)           :: r
    !!--++    real(kind=sp),              intent( in)           :: Theta
    !!--++    real(kind=sp),              intent( in)           :: Phi
    !!--++    real(kind=sp), dimension(3),intent(out)           :: xo
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the Cartesian coordinates from spherical coordinates.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cart_from_Spher_sp(r,Theta,Phi,Xo,Mode)
       !---- Arguments ----!
       real(kind=sp),              intent( in)           :: r
       real(kind=sp),              intent( in)           :: Theta
       real(kind=sp),              intent( in)           :: phi
       real(kind=sp), dimension(3),intent(out)           :: xo
       character(len=*),           intent( in), optional :: mode

       !---- Local Variables ----!
       real(kind=sp) :: ph,th

       ph=Phi
       th=Theta
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             ph=Phi*to_rad
             th=Theta*to_rad
          end if
       end if
       xo(1)=r*cos(ph)*sin(th)
       xo(2)=r*sin(ph)*sin(th)
       xo(3)=r*cos(th)

       return
    End Subroutine Get_Cart_from_Spher_sp

    !!----
    !!---- Subroutine Get_Plane_from_Points(P1,P2,P3,A,B,C,D)
    !!----    real(kind=sp), dimension(3), intent(in) :: P1
    !!----    real(kind=sp), dimension(3), intent(in) :: P2
    !!----    real(kind=sp), dimension(3), intent(in) :: P3
    !!----    real(kind=sp),               intent(out):: A
    !!----    real(kind=sp),               intent(out):: B
    !!----    real(kind=sp),               intent(out):: C
    !!----    real(kind=sp),               intent(out):: D
    !!----
    !!----    Caculate the implicit form of a Plane in 3D as
    !!----    A * X + B * Y + C * Z + D = 0
    !!----
    !!---- Update: July - 2005
    !!
    Subroutine Get_Plane_from_Points(P1, P2, P3, A, B, C, D)
       !---- Arguments ----!
       real(kind=sp), dimension(3), intent(in) :: P1
       real(kind=sp), dimension(3), intent(in) :: P2
       real(kind=sp), dimension(3), intent(in) :: P3
       real(kind=sp),               intent(out):: A
       real(kind=sp),               intent(out):: B
       real(kind=sp),               intent(out):: C
       real(kind=sp),               intent(out):: D

       a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
           - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

       b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
           - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

       c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
           - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

       d = - p2(1) * a - p2(2) * b - p2(3) * c

       return
    End Subroutine Get_Plane_from_Points

    !!----
    !!---- Subroutine Get_Spheric_Coord(Xo,Ss,Theta,Phi,Mode)
    !!----    real(kind=sp/dp), dimension(3),intent( in)           :: xo
    !!----    real(kind=sp/dp),              intent(out)           :: ss
    !!----    real(kind=sp/dp),              intent(out)           :: theta
    !!----    real(kind=sp/dp),              intent(out)           :: phi
    !!----    character(len=*),              intent( in), optional :: mode
    !!----
    !!----    Determine the spheric coordinates from rectangular coordinates.
    !!----    If Mode='D' the angles will be done in Degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Get_Spheric_Coord_dp(Xo,Ss,Theta,Phi,Mode)
    !!--++    real(kind=dp), dimension(3),intent( in)           :: xo
    !!--++    real(kind=dp),              intent(out)           :: ss
    !!--++    real(kind=dp),              intent(out)           :: theta
    !!--++    real(kind=dp),              intent(out)           :: phi
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spheric_Coord_dp(xo,ss,theta,phi,mode)
       !---- Arguments ----!
       real(kind=dp), intent( in), dimension(3)   :: xo
       real(kind=dp), intent(out)                 :: ss
       real(kind=dp), intent(out)                 :: theta
       real(kind=dp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       ss=0.0_dp
       do j=1,3
          ss=ss+xo(j)*xo(j)
       end do
       ss=sqrt(ss)
       if (ss > 0.0_dp) then
          theta=xo(3)/ss
          if (abs(theta) > 1.0_dp) then
             theta=sign(1.0_dp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_dp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_dp
          phi=0.0_dp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spheric_Coord_dp

    !!--++
    !!--++ Subroutine Get_Spheric_Coord_sp(Xo,Ss,Theta,Phi,Mode)
    !!--++    real(kind=sp), dimension(3),intent( in)           :: xo
    !!--++    real(kind=sp),              intent(out)           :: ss
    !!--++    real(kind=sp),              intent(out)           :: theta
    !!--++    real(kind=sp),              intent(out)           :: phi
    !!--++    character(len=*),           intent( in), optional :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determine the spheric coordinates from rectangular coordinates
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Spheric_Coord_sp(xo,ss,theta,phi,mode)
       !---- Arguments ----!
       real(kind=sp), intent( in), dimension(3)   :: xo
       real(kind=sp), intent(out)                 :: ss
       real(kind=sp), intent(out)                 :: theta
       real(kind=sp), intent(out)                 :: phi
       character(len=*), intent(in), optional     :: mode

       !---- Local Variables ----!
       integer :: j

       ss=0.0_sp
       do j=1,3
          ss=ss+xo(j)*xo(j)
       end do
       ss=sqrt(ss)
       if (ss > 0.0_sp) then
          theta=xo(3)/ss
          if (abs(theta) > 1.0_sp) then
             theta=sign(1.0_sp,theta)
          end if
          theta=acos(theta)
          if (abs(theta) < eps .or. abs(theta-pi) < eps) then
             phi=0.0_sp
          else
             phi=atan2(xo(2),xo(1))
          end if
       else
          theta=0.0_sp
          phi=0.0_sp
       end if
       if (present(mode)) then
          if (mode(1:1) == "D" .or. mode(1:1) == "d") then
             theta=theta*to_deg
             phi=phi*to_deg
          end if
       end if

       return
    End Subroutine Get_Spheric_Coord_sp

    !!----
    !!---- Subroutine Matrix_DiagEigen(A, V, C)
    !!----    real(kind=sp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=sp), dimension(3),   intent(out) :: v
    !!----    real(kind=sp), dimension(3,3), intent(out) :: c
    !!----
    !!----    Diagonalize the matrix A, put eigenvalues in V and
    !!----    eigenvectors in C
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matrix_DiagEigen(a,v,c)
       !---- Arguments ----!
       real(kind=sp), intent(in)  , dimension(3,3)    :: a
       real(kind=sp), intent(out) , dimension(3)      :: v
       real(kind=sp), intent(out) , dimension(3,3)    :: c

       !---- Local Variables ----!
       integer, parameter            :: n=3
       integer                       :: i, j, k, itmax, nm1, ip1, iter
       real(kind=sp), dimension(3)   :: u
       real(kind=sp), dimension(3,3) :: e
       real(kind=sp), parameter      :: eps1=1.e-7 , eps2=1.e-7 , eps3=1.e-7
       real(kind=sp)                 :: sigma1, offdsq, p, q, spq, csa, sna
       real(kind=sp)                 :: holdik, holdki, sigma2

       call init_err_math3d()
       nm1=n-1
       itmax=50
       do i=1,n
          do j=1,n
             e(i,j)=a(i,j)
             c(i,j)=0.0
             if (j < i) e(i,j)=0.0
          end do
       end do
       sigma1=0.0
       offdsq=0.0

       do i=1,n
          sigma1=sigma1+e(i,i)**2
          c(i,i)=1.0
          ip1=i+1
          if (i >= n) exit
          do j=ip1,n
             offdsq=offdsq+e(i,j)**2
          end do
       end do

       do iter=1,itmax
          do i=1,nm1
             ip1=i+1
             do j=ip1,n
                q=abs(e(i,i)-e(j,j))
                if (q <= eps1) then
                   csa=1.0/sqrt(2.0)
                   sna=csa
                else
                   if (abs(e(i,j)) <= eps2) then
                      e(i,j)=0.0
                      cycle
                   end if
                   p=2.0*e(i,j)*q/(e(i,i)-e(j,j))
                   spq=sqrt(p*p+q*q)
                   csa=sqrt((1.0+q/spq)/2.0)
                   sna=p/(2.0*csa*spq)
                end if
                do k=1,n
                   holdki=c(k,i)
                   c(k,i)=holdki*csa+c(k,j)*sna
                   c(k,j)=holdki*sna-c(k,j)*csa
                end do
                do k=i,n
                   if (k > j) then
                      holdik=e(i,k)
                      e(i,k)=csa*holdik+sna*e(j,k)
                      e(j,k)=sna*holdik-csa*e(j,k)
                   else
                      u(k)=e(i,k)
                      e(i,k)=csa*u(k)+sna*e(k,j)
                      if (k /= j) cycle
                      e(j,k)=sna*u(k)-csa*e(j,k)
                   end if
                end do
                u(j)=sna*u(i)-csa*u(j)
                do k=1,j
                   if (k <= i)  then
                      holdki=e(k,i)
                      e(k,i)=csa*holdki+sna*e(k,j)
                      e(k,j)=sna*holdki-csa*e(k,j)
                   else
                      e(k,j)=sna*u(k)-csa*e(k,j)
                   end if
                end do
                e(i,j)=0.0
             end do
          end do
          sigma2=0.0
          do i=1,n
             v(i)=e(i,i)
             sigma2=sigma2+v(i)*v(i)
          end do
          if (1.0-sigma1/sigma2 <= eps3) return
          sigma1=sigma2
       end do

       Err_Math_3D =.true.
       Err_Mess_Math_3D=" Convergence not reached in diagonalization "

       return
    End Subroutine Matrix_DiagEigen

    !!----
    !!---- Subroutine Matrix_Inverse(A, B, Ifail)
    !!----    real(kind=sp), dimension(3,3), intent(in)  :: a
    !!----    real(kind=sp), dimension(3,3), intent(out) :: b
    !!----    integer                      , intent(out) :: ifail
    !!----                                                  0 = OK; 1 = Fail
    !!----
    !!----    Inverts a 3x3 Matrix
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Matrix_Inverse(a,b,ifail)
       !---- Argument ----!
       real(kind=sp), dimension(3,3), intent(in)  :: a
       real(kind=sp), dimension(3,3), intent(out) :: b
       integer                      , intent(out) :: ifail

       !---- Local variables ----!
       real(kind=sp), parameter :: epso=1.0e-20
       real(kind=sp)            :: dmat

       ifail=0
       call init_err_math3d()

       b(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
       b(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
       b(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
       b(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
       b(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
       b(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
       b(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
       b(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
       b(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)
       dmat = a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)

       if (abs(dmat) < epso) then
          ifail=1
          Err_Math_3D =.true.
          Err_Mess_Math_3D="Singular Matrix: inversion impossible"
          return
       end if

       b = b/dmat

       return
    End Subroutine Matrix_Inverse

    !!----
    !!---- Subroutine Resolv_Sist_1X2(W,T,Ts,X,Ix)
    !!----    integer,       dimension(2),      intent(in) :: w     !  In -> Input vector
    !!----    real(kind=sp),                    intent(in) :: t     !  In -> Input value
    !!----    real(kind=sp), dimension(2),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=sp), dimension(2),      intent(out):: x     ! Out -> Fixed value for x,y
    !!----    integer, dimension(2),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----              w11 x1 + w12 x2  = t1
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2), intent( in) :: w
       real(kind=sp),                 intent( in) :: t
       real(kind=sp), dimension(2),   intent(out) :: ts
       real(kind=sp), dimension(2),   intent(out) :: x
       integer,dimension(2), intent(out) :: ix

       !---- Initialize ----!
       ts = 0.0
       x  = 1.0
       ix = 0
       call init_err_math3d()

       !---- Both are zeros ----!
       if ( all(w == 0)) then
          if (abs(t) < eps) then
             ix(1)=1
             ix(2)=2
          else
             err_math_3d=.true.
             err_mess_math_3d="Inconsistent solution (1x2)"
          end if
          return
       end if

       !---- Any is zero ----!
       if (any(w == 0)) then
          if ( w(1) == 0 ) then
             ix(1)=1
             ts(2)=t/real(w(2))
              x(2)=0.0
          else
             ts(1)=t/real(w(1))
              x(1)=0.0
             ix(2)=2
          end if
       else
          ix(1)=1
          ts(2)=t/real(w(2))
           x(2)=-real(w(1))/real(w(2))
          ix(2)=1
       end if

       return
    End Subroutine Resolv_Sist_1x2

    !!----
    !!---- Subroutine Resolv_Sist_1X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(3),            intent(in) :: w     !  In -> Input vector
    !!----    real(kind=sp),                    intent(in) :: t     !  In -> Input value
    !!----    real(kind=sp), dimension(3),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=sp), dimension(3),      intent(out):: x     ! Out -> Fixed value for x,y,z
    !!----    integer, dimension(3),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_1x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(3), intent( in) :: w
       real(kind=sp),                 intent( in) :: t
       real(kind=sp), dimension(3),   intent(out) :: ts
       real(kind=sp), dimension(3),   intent(out) :: x
       integer,dimension(3), intent(out) :: ix

       !---- Local Variables ----!
       integer               :: i, zeros
       integer, dimension(2) :: w1
       integer, dimension(2) :: ix1
       real(kind=sp), dimension(2)    :: ts1
       real(kind=sp), dimension(2)    :: x1

       !---- Initialize ----!
       ts = 0.0
       x  = 1.0
       ix = 0
       call init_err_math3d()

       !---- Are there zeros? ----!
       zeros=0
       do i=1,3
          if (w(i) == 0) zeros=zeros+1
       end do
       select case (zeros)
          case (3)
             if (abs(t) < eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                err_math_3d=.true.
                err_mess_math_3d="Inconsistent solution (1 x 3)"
             end if

          case (2)
             do i=1,3
                if (w(i) /= 0) then
                   ts(i)=t/real(w(i))
                   x(i) =0.0
                else
                   ix(i)=i
                end if
             end do

          case (1)
             do i=1,3
                if (w(i) == 0) exit
             end do
             select case (i)
                case (1)
                   w1=w(2:3)

                case (2)
                   w1(1)=w(1)
                   w1(2)=w(3)

                case (3)
                   w1=w(1:2)
             end select
             call resolv_sist_1x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x(2:3) = x1
                   if (ix1(1)==1) ix(2)=2
                   if (ix1(1)==2) ix(2)=3
                   if (ix1(2)==1) ix(3)=2
                   if (ix1(2)==2) ix(3)=3

                  case (2)
                     ix(2)= 2
                     ts(1)= ts1(1)
                     ts(3)= ts1(2)
                     x(1) = x1(1)
                     x(3) = x1(2)
                     if (ix1(1)==1) ix(1)=1
                     if (ix1(1)==2) ix(1)=3
                     if (ix1(2)==1) ix(3)=1
                     if (ix1(2)==2) ix(3)=3

                  case (3)
                     ix(3)  = 3
                     ts(1:2)= ts1
                     x(1:2) = x1
                     ix(1:2)= ix1
               end select

          case (0)
             err_math_3d=.true.
             err_mess_math_3d="Inconsistent case ax+by+cz=t (1x3)"
       end select

       return
    End Subroutine Resolv_Sist_1x3

    !!----
    !!---- Subroutine Resolv_Sist_2X2(W,T,Ts,X,Ix)
    !!----    integer, dimension(2,2),          intent(in) :: w     !  In -> Input vector
    !!----    real(kind=sp), dimension(2),      intent(in) :: t     !  In -> Input value
    !!----    real(kind=sp), dimension(2),      intent(out):: ts    ! Out -> Fixed value of solution
    !!----    real(kind=sp), dimension(2),      intent(out):: x     ! Out -> Fixed value for x,y
    !!----    integer, dimension(2),            intent(out):: ix    ! Out -> determine if solution
    !!----                                                                   1: x, 2: y, 3: z
    !!--<<
    !!----                 w11 x1 + w12 x2  = t1
    !!----                 w21 x1 + w22 x2  = t2
    !!----                 x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x2(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2,2), intent( in) :: w
       real(kind=sp),dimension(2),      intent( in) :: t
       real(kind=sp),dimension(2),      intent(out) :: ts
       real(kind=sp),dimension(2),      intent(out) :: x
       integer,dimension(2),   intent(out) :: ix

       !---- Local Variables ----!
       integer                 :: i,deter
       integer, dimension(2)   :: zeros,colum
       real(kind=sp)           :: rden, rnum

       !---- Initialize ----!
       ts    = 0.0
       x     = 1.0
       ix    = 0
       call init_err_math3d()

       deter = w(1,1)*w(2,2) - w(1,2)*w(2,1)
       rden=real(deter)
       if (deter /= 0) then
          !---- X(1) ----!
          rnum=t(1)*w(2,2) - w(1,2)*t(2)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rnum=w(1,1)*t(2) - t(1)*w(2,1)
          ts(2)=rnum/rden

          x =0.0

       else                        ! Singular Matrix
          !---- Are there zero rows? ----!
          zeros=0
          do i=1,2
             if (w(i,1) == 0 .and. w(i,2) == 0 )  zeros(i)=1
          end do
          select case (sum(zeros))
             case (2)
                if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                   ix(1)=1
                   ix(2)=2
                else
                   err_math_3d=.true.
                   err_mess_math_3d="Inconsistent solution (2x2)"
                end if

             case (1)
                do i=1,2
                   if (zeros(i) == 0) exit
                end do
                call resolv_sist_1x2(w(i,:),t(i),ts,x,ix)

             case (0)
                !---- Are there zero columns? ----!
                colum=0
                do i=1,2
                   if (w(1,i) == 0 .and. w(2,i) == 0 ) colum(i)=1
                end do
                select case (sum(colum))
                   case (1)
                      do i=1,2
                         if (colum(i) == 0) exit
                      end do
                      if (w(1,i) /= 0) then
                         ts(i)=t(1)/real(w(1,i))
                      else
                         ts(i)=t(2)/real(w(2,i))
                      end if
                      x(i)=0.0
                      if (i == 1) then
                         ix(2)=2
                      else
                         ix(1)=1
                      end if

                   case (0)
                      call resolv_sist_1x2(w(1,:),t(1),ts,x,ix)

                end select
          end select
       end if

       return
    End Subroutine Resolv_Sist_2x2

    !!----
    !!---- Subroutine Resolv_Sist_2X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(2,3),          intent(in) :: w      !  In -> Input vector
    !!----    real(kind=sp), dimension(2),      intent(in) :: t      !  In -> Input value
    !!----    real(kind=sp), dimension(3),      intent(out):: ts     ! Out -> Fixed value of solution
    !!----    real(kind=sp), dimension(3),      intent(out):: x      ! Out -> Fixed value for x,y
    !!----    integer, dimension(3),            intent(out):: ix     ! Out -> determine if solution
    !!----                                                                    1: x, 2: y, 3: z
    !!----               w11 x1 + w12 x2 + w13 x3 = t1
    !!----               w21 x1 + w22 x2 + w23 x3 = t2
    !!----               x_sol(i)= ts(i) + x(i) ix(i)
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Resolv_Sist_2x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer,dimension(2,3),          intent( in) :: w
       real(kind=sp),dimension(2),      intent( in) :: t
       real(kind=sp),dimension(3),      intent(out) :: ts
       real(kind=sp),dimension(3),      intent(out) :: x
       integer,dimension(3),            intent(out) :: ix

       !---- Local Variables ----!
       integer                 :: i, j
       integer, dimension(2)   :: fila
       integer, dimension(2)   :: ix1
       integer, dimension(3)   :: colum
       integer, dimension(2,2) :: w1
       integer, dimension(2,3) :: wm
       integer, dimension(2)   :: wc
       real(kind=sp)                    :: tc
       real(kind=sp), dimension(2)      :: tm
       real(kind=sp), dimension(2)      :: ts1, x1

       !---- Initialize ----!
       ts    = 0.0
       x     = 1.0
       ix    = 0
       call init_err_math3d()

       !---- Are there zero columns? ----!
       colum=0
       do i=1,3
            if (all(w(:,i) == 0)) colum(i)=1
       end do
       select case (sum(colum))
          case (3)
             if (abs(t(1)) <= eps .and. abs(t(2)) <= eps) then
                do i=1,3
                   ix(i)=i
                end do
             else
                err_math_3d=.true.
                err_mess_math_3d="Inconsistent solution in (2x3)"
             end if

          case (2)
             do i=1,3
                if (colum(i) == 0) exit
             end do
             if (w(1,i) /= 0) then
                ts(i)=t(1)/real(w(1,i))
             else
                ts(i)=t(2)/real(w(2,i))
             end if
             x(i)=0.0
             select case (i)
                case (1)
                   ix(2)=2
                   ix(3)=3

                case (2)
                   ix(1)=1
                   ix(3)=3

                case (3)
                   ix(1)=1
                   ix(2)=2
             end select

          case (1)
             do i=1,3
                if (colum(i) == 1) exit
             end do
             select case (i)
                case (1)
                   w1=w(:,2:3)

                case (2)
                   w1(1,1)=w(1,1)
                   w1(1,2)=w(1,3)
                   w1(2,1)=w(2,1)
                   w1(2,2)=w(2,3)

                case (3)
                   w1=w(:,1:2)
             end select
             call resolv_sist_2x2(w1,t,ts1,x1,ix1)
             select case (i)
                case (1)
                   ix(1)  = 1
                   ts(2:3)= ts1
                   x (2:3)= x1
                   if (ix1(1) == 1) ix(2)=2
                   if (ix1(1) == 2) ix(2)=3
                   if (ix1(2) == 1) ix(3)=2
                   if (ix1(2) == 2) ix(3)=3

                case (2)
                   ix(2)=2
                   ts(1)=ts1(1)
                   ts(3)=ts1(2)
                   x(1) = x1(1)
                   x(3) = x1(2)
                   if (ix1(1) == 1) ix(1)=1
                   if (ix1(1) == 2) ix(1)=3
                   if (ix1(2) == 1) ix(3)=1
                   if (ix1(2) == 2) ix(3)=3

                case (3)
                   ix(3)  = 3
                   ts(1:2)= ts1
                   x (1:2)= x1
                   ix(1:2)= ix1
             end select

          case (0)
             !---- Are there zeros in any element of rows? ----!
             fila = 0
             do i=1,2
                if (all(w(i,:)==0)) fila(i)=1
             end do
             select case (sum(fila))
                case (1)
                   if (w(1,1) /= 0) then
                      call resolv_sist_1x3(w(1,:),t(1),ts,x,ix)
                   else
                      call resolv_sist_1x3(w(2,:),t(2),ts,x,ix)
                   end if

                case (0)
                   fila = 0
                   wm   = w
                   tm   = t
                   !---- Are there zeros in any element of rows? ----!
                   do i=1,2
                      do j=1,3
                         if (w(i,j)==0) fila(i)=fila(i)+1
                      end do
                   end do
                   if ( fila(2) > fila(1) ) then
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(1,:)
                      tm(1)  =t(2)
                      tm(2)  =t(1)
                          j  =fila(1)
                      fila(1)=fila(2)
                      fila(2)=j
                   end if
                   select case (fila(1))
                      case (2)
                         do i=1,3
                            if (wm(1,i) /= 0) exit
                         end do
                         ts(i)=tm(1)/real(wm(1,i))
                         x(i)=0.0
                         select case (i)
                            case (1)
                               wc(1)=wm(2,2)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,1)*ts(i))

                            case (2)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,3)
                               tc=tm(2)-(wm(2,2)*ts(i))

                            case (3)
                               wc(1)=wm(2,1)
                               wc(2)=wm(2,2)
                               tc=tm(2)-(wm(2,3)*ts(i))
                         end select
                         call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                         select case(i)
                            case (1)
                               ts(2:3)=ts1
                                x(2:3)=x1
                                if (ix1(1)==1) ix(2)=2
                                if (ix1(1)==2) ix(2)=3
                                if (ix1(2)==1) ix(3)=2
                                if (ix1(2)==2) ix(3)=3

                            case (2)
                               ts(1)=ts1(1)
                               ts(3)=ts1(2)
                                x(1)=x1(1)
                                x(3)=x1(2)
                                if (ix1(1)==1) ix(1)=1
                                if (ix1(1)==2) ix(1)=3
                                if (ix1(2)==1) ix(3)=1
                                if (ix1(2)==2) ix(3)=3

                            case (3)
                               ts(1:2)=ts1
                                x(1:2)=x1
                               ix(1:2)=ix1
                         end select

                      case (1)
                         do i=1,3
                            if (wm(1,i) == 0) exit
                         end do
                         select case (fila(2))
                            case (1)
                               do j=1,3
                                  if (wm(2,j) == 0) exit
                               end do
                               select case (i)
                                  case (1)             ! 0 en w(1,1)
                                     select case (j)
                                        case (2)
                                           wc(1)=-wm(2,1)/wm(2,3)
                                           wc(2)= wm(1,2)/wm(1,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(2)/real(wm(2,3)) - ts(1)*wm(2,1)/real(wm(2,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(1)/real(wm(1,3)) - ts(2)*wm(1,2)/real(wm(1,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3)=-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) =-real(wc(1))/real(wc(2))
                                                 ix(2)=1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,1)/wm(2,2)
                                           wc(2)= wm(1,3)/wm(1,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(2)/real(wm(2,2)) - ts(1)*wm(2,1)/real(wm(2,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(1)/real(wm(1,2)) - ts(3)*wm(1,3)/real(wm(1,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(3)=-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (2)             ! 0 en w(1,2)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,3)
                                           wc(2)=-wm(2,2)/wm(2,3)
                                           tc=tm(1)/real(wm(1,3)) - tm(2)/real(wm(2,3))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1:2)=ts1
                                           x(1:2) =x1
                                           ix(1:2)=ix1
                                           if (ix(1) == 0) then
                                              ts(3)=tm(1)/real(wm(1,3)) - ts(1)*wm(1,1)/real(wm(1,3))
                                              x(3)=0.0
                                           else
                                              if (ix(2) == 0) then
                                                 ts(3)=tm(2)/real(wm(2,3)) - ts(2)*wm(2,2)/real(wm(2,3))
                                                 x(3)=0.0
                                              else
                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3)=-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1

                                                 ts(2)=tc/real(wc(2))
                                                 x(2) = -real(wc(1))/real(wc(2))
                                                 ix(2)= 1
                                              end if
                                           end if

                                        case (3)
                                           wc(1)=-wm(2,2)/wm(2,1)
                                           wc(2)= wm(1,3)/wm(1,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(2)/real(wm(2,1)) - ts(2)*wm(2,2)/real(wm(2,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(1)/real(wm(1,1)) - ts(3)*wm(1,3)/real(wm(1,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(2)/real(wm(2,2))
                                                 x(2) =-real(wm(2,1))/real(wm(2,2))
                                                 ix(2)=1

                                                 ts(3)=tm(1)/real(wm(1,3))
                                                 x(3) =-real(wm(1,1))/real(wm(1,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select

                                  case (3)             ! 0 en w(1,3)
                                     select case (j)
                                        case (1)
                                           wc(1)= wm(1,1)/wm(1,2)
                                           wc(2)=-wm(2,3)/wm(2,2)
                                           tc=tm(1)/real(wm(1,2)) - tm(2)/real(wm(2,2))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(1)=ts1(1)
                                           ts(3)=ts1(2)
                                           x(1) =x1(1)
                                           x(3) =x1(2)
                                           if (ix1(1) == 1) ix(1)=1
                                           if (ix1(1) == 2) ix(1)=3
                                           if (ix1(2) == 1) ix(3)=1
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(1) == 0) then
                                              ts(2)=tm(1)/real(wm(1,2)) - ts(1)*wm(1,1)/real(wm(1,2))
                                              x(2)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(2)=tm(2)/real(wm(2,2)) - ts(3)*wm(2,3)/real(wm(2,2))
                                                 x(2)=0.0
                                              else
                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tc/real(wc(2))
                                                 x(3) =-real(wc(1))/real(wc(2))
                                                 ix(3)=1
                                              end if
                                           end if

                                        case (2)
                                           wc(1)= wm(1,2)/wm(1,1)
                                           wc(2)=-wm(2,3)/wm(2,1)
                                           tc=tm(1)/real(wm(1,1)) - tm(2)/real(wm(2,1))
                                           call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                           ts(2:3)=ts1
                                           x(2:3) =x1
                                           if (ix1(1) == 1) ix(2)=2
                                           if (ix1(1) == 2) ix(2)=3
                                           if (ix1(2) == 1) ix(3)=2
                                           if (ix1(2) == 2) ix(3)=3
                                           if (ix(2) == 0) then
                                              ts(1)=tm(1)/real(wm(1,1)) - ts(2)*wm(1,2)/real(wm(1,1))
                                              x(1)=0.0
                                           else
                                              if (ix(3) == 0) then
                                                 ts(1)=tm(2)/real(wm(2,1)) - ts(3)*wm(2,3)/real(wm(2,1))
                                                 x(1)=0.0
                                              else
                                                 ix(1)=1

                                                 ts(2)=tm(1)/real(wm(1,2))
                                                 x(2) =-real(wm(1,1))/real(wm(1,2))
                                                 ix(2)=1

                                                 ts(3)=tm(2)/real(wm(2,3))
                                                 x(3) =-real(wm(2,1))/real(wm(2,3))
                                                 ix(3)=1
                                              end if
                                           end if
                                     end select
                               end select

                            case (0)
                               select case (i)
                                  case (1)
                                     wc(1)=wm(2,1)
                                     wc(2)=wm(2,2)- wm(2,3)*wm(1,2)/wm(1,3)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                     ts(1:2)=ts1
                                     x(1:2)=x1
                                     ix(1:2)=ix1
                                     if (ix(2) == 0) then
                                        ts(3)=tm(1)/real(wm(1,3)) - ts(2)*real(wm(1,2))/real(wm(1,3))
                                        x(3)=0.0
                                     else
                                        ix(1)=1

                                        ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3))) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        x(2) =-real(wm(2,1)) / &
                                              (real(wm(2,2)) - real(wm(2,3)*wm(1,2))/real(wm(1,3)) )
                                        ix(2)=1

                                        ts(3)= tm(1)/real(wm(1,3)) - (real(wm(1,2))/real(wm(1,3)))*ts(2)
                                        x(3) =- (real(wm(1,2))/real(wm(1,3)))*x(2)
                                        ix(3)=1
                                     end if

                                  case (2)
                                     wc(1)=wm(2,1)-wm(2,3)*wm(1,1)/wm(1,3)
                                     wc(2)=wm(2,2)
                                     tc=tm(2)-tm(1)*wm(2,3)/real(wm(1,3))
                                     call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1:2)=ts1
                                    x(1:2)=x1
                                    ix(1:2)=ix1
                                    if (ix(1) == 0) then
                                       ts(3)=tm(1)/real(wm(1,3)) - ts(1)*real(wm(1,1))/real(wm(1,3))
                                       x(3)=0.0
                                    else
                                       ix(1)=1

                                       ts(2)=(tm(2) - tm(1)*wm(2,3)/real(wm(1,3)))/real(wm(2,2))
                                       x(2) =(real(wm(1,1)*wm(2,3))/real(wm(1,3)) - real(wm(2,1)))/real(wm(2,2))
                                       ix(2)=1

                                       ts(3)=tm(1)/real(wm(1,3))
                                       x(3) =-real(wm(1,1))/real(wm(1,3))
                                       ix(3)=1
                                    end if

                                 case (3)
                                    wc(1)=wm(2,1)-wm(1,1)*wm(2,2)/wm(1,2)
                                    wc(2)=wm(2,3)
                                    tc=tm(2)-tm(1)*wm(2,2)/real(wm(1,2))
                                    call resolv_sist_1x2(wc,tc,ts1,x1,ix1)
                                    ts(1)=ts1(1)
                                    ts(3)=ts1(2)
                                    x(1)=x1(1)
                                    x(3)=x1(2)
                                    if (ix1(1) == 1) ix(1)=1
                                    if (ix1(1) == 2) ix(1)=3
                                    if (ix1(2) == 1) ix(3)=1
                                    if (ix1(2) == 2) ix(3)=3
                                    if (ix(1) == 0) then
                                       ts(2)=tm(1)/real(wm(1,2)) - ts(1)*real(wm(1,1))/real(wm(1,2))
                                       x(2)=0.0
                                    else
                                       ix(1) =1

                                       ts(2)=tm(1)/real(wm(1,2))
                                       x(2) =-real(wm(1,1))/real(wm(1,2))
                                       ix(2)=1

                                       ts(3)=(tm(2) - tm(1)*wm(2,2)/real(wm(1,2)))/real(wm(2,3))
                                       x(3) =(real(wm(1,1)*wm(2,2))/real(wm(1,2)) - real(wm(2,1)))/real(wm(2,3))
                                       ix(3)=1
                                    end if
                               end select
                         end select

                      case (0)
                         call resolv_sist_1x3(wm(1,:),tm(1),ts,x,ix)
                   end select

             end select
       end select

       return
    End Subroutine Resolv_Sist_2x3

    !!----
    !!---- Subroutine Resolv_Sist_3X3(W,T,Ts,X,Ix)
    !!----    integer, dimension(3,3),          intent(in) :: w      !  In -> Input vector
    !!----    real(kind=sp), dimension(3),      intent(in) :: t      !  In -> Input value
    !!----    real(kind=sp), dimension(3),      intent(out):: ts     ! Out -> Fixed value of solution
    !!----    real(kind=sp), dimension(3),      intent(out):: x      ! Out -> Fixed value for x,y
    !!----    integer, dimension(3),            intent(out):: ix     ! Out -> determine if solution
    !!----                                                                     1: x, 2: y, 3: z
    !!--<<
    !!----              w11 x1 + w12 x2 + w13 x3 = t1
    !!----              w21 x1 + w22 x2 + w23 x3 = t2
    !!----              w31 x1 + w32 x2 + w33 x3 = t3
    !!----              x_sol(i)= ts(i) + x(i) ix(i)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Resolv_Sist_3x3(w,t,ts,x,ix)
       !---- Arguments ----!
       integer, dimension(3,3),          intent(in) :: w
       real(kind=sp), dimension(3),      intent(in) :: t
       real(kind=sp), dimension(3),      intent(out):: ts
       real(kind=sp), dimension(3),      intent(out):: x
       integer, dimension(3),            intent(out):: ix

       !---- Local variables ----!
       integer                 :: i,j,deter
       integer, dimension(3)   :: fila
       integer, dimension(3,3) :: w1
       integer, dimension(2,3) :: wm
       real(kind=sp)                    :: rnum, rden
       real(kind=sp), dimension(3)      :: t1
       real(kind=sp), dimension(2)      :: tm
       real(kind=sp),dimension(3,3)     :: rw

       !---- Initialize ----!
       ts  = 0.0
       x   = 1.0
       ix  = 0
       call init_err_math3d()

       deter=determ_a(w)
       rden=real(deter)

       if (deter /= 0) then
          !---- X(1) ----!
          rw=real(w)
          rw(:,1)=t
          rnum=determ_a(rw)
          ts(1)=rnum/rden

          !---- X(2) ----!
          rw=real(w)
          rw(:,2)=t
          rnum=determ_a(rw)
          ts(2)=rnum/rden

          !---- X(3) ----!
          rw=real(w)
          rw(:,3)=t
          rnum=determ_a(rw)
          ts(3)=rnum/rden

          x=0.0

       else                     !  Singular Matrix
          !---- Are there zero rows? ----!
          fila=0
          do i=1,3
             if (all(w(i,:) == 0)) fila(i)=1
          end do
          select case (sum(fila))
             !---- All values are zeros ----!
             case (3)
                if (all(abs(t) < eps)) then
                   do i=1,3
                      ix(i)=i
                   end do
                else
                   err_math_3d=.true.
                   err_mess_math_3d="Inconsistent system (3 x 3)"
                end if

             !---- Two rows with zeroes ----!
             case (2)
                do i=1,3
                   if (fila(i) == 0) exit
                end do
                call resolv_sist_1x3(w(i,:),t(i),ts,x,ix)

             !---- One row with zeroes ----!
             case (1)
                do i=1,3
                   if (fila(i) == 1) exit
                end do
                select case(i)
                   case (1)
                      wm(1,:)=w(2,:)
                      wm(2,:)=w(3,:)
                      tm=t(2:3)

                   case (2)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(3,:)
                      tm(1)=t(1)
                      tm(2)=t(3)

                   case (3)
                      wm(1,:)=w(1,:)
                      wm(2,:)=w(2,:)
                      tm=t(1:2)

                end select
                call resolv_sist_2x3(wm,tm,ts,x,ix)

             !---- Non zero rows ----!
             case (0)
                w1=w
                t1=t

                !---- Are there 2 rows proportional? ----!
                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(2,i)) ) then
                      if (w1(2,i) /= 0) then
                         j=w1(1,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(1,1) .and. j*w1(2,2) == w1(1,2) .and. &
                             j*w1(2,3) == w1(1,3) ) then
                            w1(1,:)=w1(2,:)
                            t1(1)  =t1(2)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(2,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(2,1) .and. j*w1(1,2) == w1(2,2) .and. &
                             j*w1(1,3) == w1(2,3) ) then
                            w1(2,:)=w1(1,:)
                            t1(2)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(1,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(1,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(1,1) .and. j*w1(3,2) == w1(1,2) .and. &
                             j*w1(3,3) == w1(1,3) ) then
                            w1(1,:)=w1(3,:)
                            t1(1)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(1,i) /= 0) then
                         j=w1(3,i)/w1(1,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(1,1) == w1(3,1) .and. j*w1(1,2) == w1(3,2) .and. &
                             j*w1(1,3) == w1(3,3) ) then
                            w1(3,:)=w1(1,:)
                            t1(3)  =t1(1)
                            exit
                         end if
                      end if
                   end if
                end do

                do i=1,3
                   if ( abs(w1(2,i)) > abs(w1(3,i)) ) then
                      if (w1(3,i) /= 0) then
                         j=w1(2,i)/w1(3,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(3,1) == w1(2,1) .and. j*w1(3,2) == w1(2,2) .and. &
                             j*w1(3,3) == w1(2,3) ) then
                            w1(2,:)=w1(3,:)
                            t1(2)  =t1(3)
                            exit
                         end if
                      end if
                   else
                      if (w1(2,i) /= 0) then
                         j=w1(3,i)/w1(2,i)
                      else
                         j=0
                      end if
                      if (j /= 0) then
                         if (j*w1(2,1) == w1(3,1) .and. j*w1(2,2) == w1(3,2) .and. &
                             j*w1(2,3) == w1(3,3) ) then
                            w1(3,:)=w1(2,:)
                            t1(3)  =t1(2)
                            exit
                         end if
                      end if
                   end if
                end do

                !---- Are there 3 rows equal? ----!
                if ( (w1(1,1) == w1(2,1)) .and. (w1(1,1) == w1(3,1)) .and. &
                     (w1(1,2) == w1(2,2)) .and. (w1(1,2) == w1(3,2)) .and. &
                     (w1(1,3) == w1(2,3)) .and. (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_1x3(w1(1,:),t1(1),ts,x,ix)

                !---- Are there 2 rows equal? ----!
                elseif( (w1(1,1) == w1(2,1)) .and. (w1(1,2) == w1(2,2)) .and. &
                        (w1(1,3) == w1(2,3)) ) then

                   call resolv_sist_2x3(w1(2:3,:),t1(2:3),ts,x,ix)

                elseif( (w1(1,1) == w1(3,1)) .and. (w1(1,2) == w1(3,2)) .and. &
                        (w1(1,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                elseif( (w1(2,1) == w1(3,1)) .and. (w1(2,2) == w1(3,2)) .and. &
                        (w1(2,3) == w1(3,3)) ) then

                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                !---- Are linear combinations? ----!
                else
                   call resolv_sist_2x3(w1(1:2,:),t1(1:2),ts,x,ix)

                end if

          end select
       end if

       return
    End Subroutine Resolv_Sist_3x3

 End Module Math_3D
