!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: GEOM_CALCULATIONS
!!----   INFO: Routines for Geometry Calculations
!!----
!!---- HISTORY
!!----    Update: January - 2004
!!----
!!----    November - 2000 Updated by JGP and JRC
!!----
!!----
!!---- DEPENDENCIES
!!--++    Math_3D:  Cross_Product
!!--++    Math_Gen: Eps, Pi, Acosd, Cp, Sp, Cosd, Sind, To_Rad, To_Deg
!!--++    Crystal_Types: Crystal_Cell_Type
!!----
!!---- VARIABLES
!!----    COORDINATION_TYPE
!!----    COORD_INFO
!!--++    EPSI
!!----    ERR_GEOM
!!----    ERR_MESS_GEOM
!!----    POINT_LIST_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       ANGLE_DIHEDRAL
!!--++       ANGLE_DIHEDRAL_IJKN            [Overloaded]
!!--++       ANGLE_DIHEDRAL_UVW             [Overloaded]
!!----       ANGLE_MOD
!!--++       ANGLE_MODN                     [Overloaded]
!!--++       ANGLE_MODV                     [Overloaded]
!!----       ANGLE_UV
!!--++       ANGLE_UVI                      [Overloaded]
!!--++       ANGLE_UVR                      [Overloaded]
!!----       COORD_MOD
!!--++       COORD_MODN                     [Overloaded]
!!--++       COORD_MODV                     [Overloaded]
!!----       DISTANCE
!!--++       DISTANCE_FR                    [Overloaded]
!!--++       DISTANCE_SC                    [Overloaded]
!!----       MATRIX_PHITHECHI
!!----       MATRIX_RX
!!----       MATRIX_RY
!!----       MATRIX_RZ
!!----
!!----    Subroutines:
!!----       ALLOCATE_COORDINATION_TYPE
!!----       ALLOCATE_POINT_LIST
!!----       CALC_DIST_ANGLE
!!----       CALC_DIST_ANGLE_SIGMA
!!----       DEALLOCATE_COORDINATION_TYPE
!!----       DEALLOCATE_POINT_LIST
!!----       DISTANCE_AND_SIGMA
!!----       GET_EULER_FROM_FRACT
!!----       GET_PHITHECHI
!!----       GET_TRANSF_LIST
!!----       INIT_ERR_GEOM
!!----       P1_DIST
!!----       PRINT_DISTANCES
!!----       SET_ORBITS_INLIST
!!----       SET_TDIST_COORDINATION
!!----
!!
 Module Geom_Calculations

    !---- Use Modules ----!
    !Use Mod_fun    !To be commented for non-F compilers
    use Math_Gen,                   only: Sp, Cp, eps, pi, acosd, cosd, sind, to_rad, to_deg, Modulo_Lat
    use Math_3D,                    only: Cross_Product, Matrix_Inverse
    use String_Utilities,           only: Frac_Trans_1Dig, L_Case,U_Case,pack_string,setnum_std, get_logunit
    use Crystal_Types,              only: Crystal_Cell_Type, Get_Deriv_Orth_Cell
    use Atom_Module,                only: atom_list_type,Atoms_Cell_Type,Equiv_Atm, Wrt_Lab
    use Crystallographic_Symmetry,  only: Space_Group_Type, ApplySo, Lattice_Trans, Get_Multip_Pos, &
                                          searchop, Read_SymTrans_Code, Write_SymTrans_Code

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!
    public :: Angle_Dihedral, Angle_Mod, Angle_Uv, Coord_Mod, Distance, Matrix_PhiTheChi, Matrix_Rx, &
              Matrix_Ry, Matrix_Rz

    !---- List of public subroutines ----!
    public :: Allocate_Coordination_Type, Allocate_Point_List, Calc_Dist_Angle, Calc_Dist_Angle_Sigma, &
              Deallocate_Coordination_Type, Deallocate_Point_List, Distance_and_Sigma, Get_Euler_From_Fract, &
              Get_PhiTheChi, init_err_geom, P1_Dist, Print_Distances, Set_Orbits_InList, Set_TDist_Coordination, &
              Get_Transf_List

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: Angle_Dihedral_Uvw,  Angle_Dihedral_Ijkn, Angle_Uvi, Angle_Uvr, Angle_Modn, Angle_Modv, &
               Coord_Modn, Coord_Modv, Distance_fr, Distance_sc

    !---- Definitions ----!

    !!----
    !!---- TYPE :: COORDINATION_TYPE
    !!--..
    !!---- Type, public :: Coordination_Type
    !!----    integer                                      :: Natoms    ! number of atoms
    !!----    integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
    !!----    integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
    !!----    integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
    !!----                                                              ! atom to the atom given by the first index
    !!----    integer,       dimension(:,:),   allocatable :: N_Sym     !
    !!----    real(kind=sp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
    !!----    real(kind=sp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
    !!----    real(kind=sp), dimension(:,:,:), allocatable :: Tr_coo    !
    !!---- End type Coordination_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Coordination_Type
       integer                                      :: Natoms    ! number of atoms
       integer                                      :: Max_Coor  ! Maximum number of connected atoms to a given one
       integer,       dimension(:),     allocatable :: Coord_Num ! Counter of distances connected to the current atom
       integer,       dimension(:,:),   allocatable :: N_Cooatm  ! Pointer to the ordinal number in the list of the attached
                                                                 ! atom to the atom given by the first index
       integer,       dimension(:,:),   allocatable :: N_Sym     ! Number of symmetry operator to apply to N_Cooatm
       real(kind=sp), dimension(:,:),   allocatable :: Dist      ! List of distances related to an atom
       real(kind=sp), dimension(:,:),   allocatable :: S_Dist    ! List of Sigma(distances)
       real(kind=sp), dimension(:,:,:), allocatable :: Tr_coo
    End type Coordination_Type

    !!----
    !!---- COORD_INFO
    !!----    type(Coordination_Type), public :: coord_info
    !!----
    !!----    Coordination Information
    !!----
    !!---- Update: March - 2005
    !!
    type(Coordination_Type), public :: coord_info

    !!--++
    !!--++ EPSI
    !!--++    real(kind=sp), parameter :: epsi=0.001
    !!--++
    !!--++    (PRIVATE)
    !!--++    Epsilon for roughly comparing distances
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), parameter, private :: epsi=0.001

    !!----
    !!---- ERR_GEOM
    !!----    logical, public  :: err_geom
    !!----
    !!----    Logical Variable indicating an error in GEOM_CALCULATIONS module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: err_geom

    !!----
    !!---- ERR_MESS_GEOM
    !!----    character(len=150), public :: err_mess_geom
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_geom


    !!----
    !!---- TYPE :: POINT_LIST_TYPE
    !!--..
    !!---- Type, public :: Point_List_Type
    !!----    integer                                       :: np   !number of points in list
    !!----    character(len=12), dimension(:),  allocatable :: nam  !name/label associated to each point
    !!----    integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
    !!----    real,              dimension(:,:),allocatable :: x    !fractional coordinates of points
    !!---- End type Point_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: point_list_type
       integer                                       :: np   !number of points in list
       character(len=12), dimension(:),  allocatable :: nam  !name/label associated to each point
       integer,           dimension(:),  allocatable :: p    !integer pointer for various purposes
       real,              dimension(:,:),allocatable :: x    !fractional coordinates of points
    End type point_list_type


    !---- Interfaces - Overlapp ----!
    Interface  Angle_Dihedral
       Module Procedure Angle_Dihedral_Ijkn
       Module Procedure Angle_Dihedral_Uvw
    End Interface

    Interface  Angle_Uv
       Module Procedure Angle_UvI
       Module Procedure Angle_UvR
    End Interface

    Interface  Angle_Mod
       Module Procedure Angle_ModN
       Module Procedure Angle_ModV
    End Interface

    Interface  Coord_Mod
       Module Procedure Coord_ModN
       Module Procedure Coord_ModV
    End Interface

    Interface  Distance
       Module Procedure Distance_FR
       Module Procedure Distance_SC
    End Interface

 Contains

    !---- Functions ----!

    !!----
    !!---- Function Angle_Dihedral(U,V,W) Or (Ri,Rj,Rk,Rn)   Result(Angle)
    !!----    real(kind=cp), dimension(3), intent( in) :: u       !  In -> Vector 1
    !!----    real(kind=cp), dimension(3), intent( in) :: v       !  In -> Vector 2
    !!----    real(kind=cp), dimension(3), intent( in) :: w       !  In -> Vector 3
    !!----    or
    !!----    real(kind=cp), dimension(3), intent( in) :: ri      !  In -> Vector position ri
    !!----    real(kind=cp), dimension(3), intent( in) :: rj      !  In -> Vector position rj
    !!----    real(kind=cp), dimension(3), intent( in) :: rk      !  In -> Vector position rk
    !!----    real(kind=cp), dimension(3), intent( in) :: rl      !  In -> Vector position rn
    !!----    real(kind=cp)                            :: angle   ! Out -> Dihedral angle
    !!----
    !!----    Calculates the dihedral angle between planes "u-v" and "v-w", where vectors U,V,W
    !!----    are given in cartesian components.
    !!----    Calculates the dihedral angle corresponding to the four points (ri,rj,rk,rn)
    !!----    given in cartesian components. The definition used for the dihedral angle
    !!----    is the following:
    !!--<<
    !!----    Phi(i,j,k,n) = acos { (rij x rjk) (rjk x rkn) / |rij x rjk| / |rjk x rkn| }
    !!----
    !!----    with this definition the sign of Phi is positive if the vector product
    !!----    (rij x rjk) x (rjk x rkn) is in the same direction as rjk, and negative if
    !!----    the direction is opposite.
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_Dihedral_Ijkn(Ri,Rj,Rk,Rn) Result(Angle)
    !!--++    real(kind=cp), dimension(3), intent( in) :: ri       !  In -> Vector position ri
    !!--++    real(kind=cp), dimension(3), intent( in) :: rj       !  In -> Vector position rj
    !!--++    real(kind=cp), dimension(3), intent( in) :: rk       !  In -> Vector position rk
    !!--++    real(kind=cp), dimension(3), intent( in) :: rl       !  In -> Vector position rn
    !!--++    real(kind=cp)                            :: angle    ! Out -> Dihedral angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the dihedral angle corresponding to the four points (ri,rj,rk,rn)
    !!--++    given in cartesian components. The definition used for the dihedral angle
    !!--++    is the following:
    !!--++
    !!--++    Phi(i,j,k,n) = acos { (rij x rjk) (rjk x rkn) / |rij x rjk| / |rjk x rkn| }
    !!--++
    !!--++    with this definition the sign of Phi is positive if the vector product
    !!--++    (rij x rjk) x (rjk x rkn) is in the same direction as rjk, and negative if
    !!--++    the direction is opposite.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_Dihedral_Ijkn(ri,rj,rk,rn) result(angle)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent( in) :: ri,rj,rk,rn
       real(kind=cp)                            :: angle

       angle=Angle_Dihedral_Uvw(rj-ri ,rk-rj, rn-rk )

       return
    End Function Angle_Dihedral_Ijkn

    !!--++
    !!--++ Function Angle_Dihedral_Uvw(U,V,W) Result(Angle)
    !!--++    real(kind=cp), dimension(3), intent( in) :: u       !  In -> Vector 1
    !!--++    real(kind=cp), dimension(3), intent( in) :: v       !  In -> Vector 2
    !!--++    real(kind=cp), dimension(3), intent( in) :: w       !  In -> Vector 3
    !!--++    real(kind=cp)                            :: angle   ! Out -> Dihedral angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the dihedral angle between planes u-v and v-w
    !!--++    Vectors u,v,w are given in cartesian components.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_Dihedral_Uvw(u,v,w) result(angle)
       !---- Argument ----!
       real(kind=cp), dimension(3), intent( in) :: u,v,w
       real(kind=cp)                            :: angle

       !---- Local variables ----!
       real(kind=cp)               :: uvmod, vwmod, sig
       real(kind=cp), dimension(3) :: uv,vw

       angle=0.0

       uv=cross_product(u,v)
       vw=cross_product(v,w)
       sig = -sign(1.0_cp, dot_product(cross_product(uv,vw),v))
       uvmod=sqrt(dot_product(uv,uv))
       vwmod=sqrt(dot_product(vw,vw))
       if (uvmod < eps .or. vwmod < eps) return
       angle=acosd(dot_product(uv,vw)/uvmod/vwmod)*sig

       return
    End Function Angle_Dihedral_Uvw

    !!----
    !!---- Function Angle_Mod(X) Result (Y)
    !!----     real(kind=cp),               intent(in) :: x
    !!----                  or
    !!----     real(kind=cp), dimension(:), intent(in) :: x
    !!----
    !!----     Calculates the angle [-pi,pi)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_Modn(Angle) Result(AngMod)
    !!--++    real(kind=cp), intent(in) :: Angle    !  In/Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms angle in radians between -pi and +pi
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_ModN(Angle) Result(AngMod)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: Angle
       real(kind=cp)             :: AngMod

       AngMod=mod(angle+6.0*pi,2.0*pi)
       if (angmod > pi) angmod=angmod-2.0*pi

       return
    End Function Angle_ModN

    !!--++
    !!--++ Function Angle_Modv(V_Angle) Result(VAngMod)
    !!--++    real(kind=cp), dimension(:), intent(in) :: V_Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms angles in radians between -pi and +pi
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_ModV(V_Angle) Result(VAngMod)
       !---- Arguments ----!
       real(kind=cp), dimension(:),intent(in) :: V_Angle
       real(kind=cp), dimension(size(V_Angle)):: VAngMod

       !---- Local Variables ----!
       integer :: i

       VAngMod=mod(V_Angle+6.0*pi,2.0*pi)
       do i=1,size(V_Angle)
          if (VAngMod(i) > pi) VAngMod(i)=VAngMod(i)-2.0*pi
       end do

       return
    End Function Angle_ModV

    !!----
    !!---- Function Angle_Uv(U,V,G) Result(Angle)
    !!----    integer/real(kind=cp), dimension(:), intent( in)     :: u      !  In -> Vector 1
    !!----    integer/real(kind=cp), dimension(:), intent( in)     :: v      !  In -> Vector 2
    !!----    real(kind=cp), dimension(:,:), intent( in), optional :: g      !  In -> Metric tensor
    !!----    real(kind=cp)                                        :: angle  ! Out -> Angle
    !!----
    !!----    Calculates the angle between vectors u and v given in cartesian
    !!----    components. If g is not given cartesian components are assumed.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Angle_UvI(Ui,Vi,G) Result(Angle)
    !!--++    integer, dimension(:),                   intent(in) :: ui      !  In -> Vector 1
    !!--++    integer, dimension(:),                   intent(in) :: vi      !  In -> Vector 2
    !!--++    real(kind=cp), dimension(:,:), optional, intent(in) :: g       !  In -> Metric tensor
    !!--++    real(kind=cp)                                       :: angle   ! Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the angle between vectors u and v given in cartesians
    !!--++    or fractional components. If g is not given cartesian components
    !!--++    are assumed.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_UvI(Ui,Vi,G) Result(Angle)
       !---- Argument ----!
       integer, dimension(:),   intent( in)                 :: ui
       integer, dimension(:),   intent( in)                 :: vi
       real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
       real(kind=cp)                                        :: angle

       !---- Local variables ----!
       real(kind=cp)                      :: umod, vmod
       real(kind=cp), dimension(size(ui)) :: u
       real(kind=cp), dimension(size(vi)) :: v

       angle=0.0

       u=real(ui)
       v=real(vi)

       if (present(g)) then
          umod = sqrt(dot_product(u,matmul(g,u)))
          vmod = sqrt(dot_product(v,matmul(g,v)))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,matmul(g,v))/umod/vmod)
       else
          umod=sqrt(dot_product(u,u))
          vmod=sqrt(dot_product(v,v))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,v)/umod/vmod)
       end if

       return
    End Function Angle_uvi

    !!--++
    !!--++ Function Angle_Uvr(U,V,G) Result(Angle)
    !!--++    real(kind=cp), dimension(:), intent( in)             :: u      !  In -> Vector 1
    !!--++    real(kind=cp), dimension(:), intent( in)             :: v      !  In -> Vector 2
    !!--++    real(kind=cp), dimension(:,:), intent( in), optional :: g      !  In -> Metric tensor
    !!--++    real(kind=cp)                                        :: angle  ! Out -> Angle
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the angle between vectors u and v given in cartesian
    !!--++    or fractional components. If g is not given cartesian components
    !!--++    are assumed.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Angle_UvR(u,v,g) result(angle)
       !---- Argument ----!
       real(kind=cp), dimension(:),   intent( in)           :: u
       real(kind=cp), dimension(:),   intent( in)           :: v
       real(kind=cp), dimension(:,:), intent( in), optional :: g   !metric tensor
       real(kind=cp)                                        :: angle

       !---- Local variables ----!
       real(kind=cp)   :: umod, vmod

       angle=0.0

       if (present(g)) then
          umod = sqrt(dot_product(u,matmul(g,u)))
          vmod = sqrt(dot_product(v,matmul(g,v)))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,matmul(g,v))/umod/vmod)
       else
          umod=sqrt(dot_product(u,u))
          vmod=sqrt(dot_product(v,v))
          if (umod < eps .or. vmod < eps) return
          angle=acosd(dot_product(u,v)/umod/vmod)
       end if

       return
    End Function Angle_uvr

    !!----
    !!---- Function Coord_Mod(X) Result (Y)
    !!----    Real(Kind=Cp),               intent(in) :: x
    !!----                  or
    !!----    real(kind=cp), dimension(:), intent(in) :: x
    !!----
    !!----    Calculates the coordinates between [0,1)
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Coord_Modn(X) Result (XMod)
    !!--++    real(kind=cp), intent(in) :: x
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms reduced the value between 0 and 1
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Coord_ModN(x) result(Xmod)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: x
       real(kind=cp)             :: xmod

       xmod=mod(x+10.0_sp,1.0_sp)

       return
    End Function Coord_ModN

    !!--++
    !!--++ Function Coord_Modv(X) Result(XMod)
    !!--++    real(kind=cp), dimension(:), intent(in) :: x
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Transforms reduced coordinate between 0 and 1
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Coord_ModV(x) Result(Xmod)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in) :: x
       real(kind=cp), dimension(size(x))       :: xmod

       xmod=mod(x+10.0_sp,1.0_sp)

       return
    End Function Coord_ModV

    !!----
    !!---- Function Distance(X0,X1,Cell or Code) Result(D)
    !!----    real(kind=sp), dimension(3),        intent(in) :: x0     !  In -> Point 1
    !!----    real(kind=sp), dimension(3),        intent(in) :: x1     !  In -> Point 2
    !!----    Type (Crystal_Cell_Type),           intent(in) :: Cell   !  In -> Cell parameters
    !!----    Or
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=sp)                                  :: d      ! Out -> Distance
    !!----
    !!----    Calculate distance between two points.
    !!----       Fractional Coordinates: Use Cell
    !!----       Cartesian Coordiantes: Code="C" or Code=" "
    !!----       Spherical Coordinates: Code="S"
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Distance_Fr(X0,X1,Celda) Result(D)
    !!--++    real(kind=sp), dimension(3),  intent(in) :: x0     !  In -> Point 1
    !!--++    real(kind=sp), dimension(3),  intent(in) :: x1     !  In -> Point 2
    !!--++    Type (Crystal_Cell_Type),     intent(in) :: Celda  !  In -> Cell parameters
    !!--++    real(kind=sp)                                  :: d      ! Put -> Distance
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate distance between two points in Fractional
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Distance_Fr(X0,X1,Celda) Result(Dis)
       !---- Arguments ----!
       real(kind=sp), dimension(3), intent(in) :: x0,x1
       type (Crystal_Cell_Type),    intent(in) :: Celda
       real(kind=sp)                           :: dis

       !---- Local Variables ----!
       real(kind=sp), dimension(3) :: xr

       xr = matmul(celda%Cr_Orth_cel,x1-x0)
       dis=sqrt(dot_product(xr,xr))

       return
    End Function Distance_Fr

    !!--++
    !!--++ Function Distance_SC(X0,X1,Code) Result(D)
    !!--++    real(kind=sp), dimension(3),        intent(in) :: x0     !  In -> Point 1
    !!--++    real(kind=sp), dimension(3),        intent(in) :: x1     !  In -> Point 2
    !!--++    character(len=*), optional,         intent(in) :: Code
    !!--++    real(kind=sp)                                  :: d      ! Put -> Distance
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculate distance between two points in Cartesian or Spherical
    !!--++    If Code =="C" or Blank or not present then the coordinates are Cartesian.
    !!--++    If Code =="S" then the coordinates are spherical (R, Theta, Phi).
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Distance_SC(X0,X1,Code) Result(Dis)
       !---- Arguments ----!
       real(kind=sp), dimension(3), intent(in) :: x0,x1
       character(len=*), optional,  intent(in) :: Code
       real(kind=sp)                           :: dis

       !---- Local Variables ----!
       real(kind=sp), dimension(3) :: xr,xi,xj

       xr=0.0
       if (present(code)) then
          select case (code(1:1))
             case("S","s") ! Spherical
                xi(1)=x0(1)*cosd(x0(3))*sind(x0(2))  ! R * cos(Phi) * sin(Theta)
                xi(2)=x0(1)*sind(x0(3))*sind(x0(2))  ! R * sin(Phi) * sin(Theta)
                xi(3)=x0(1)*cosd(x0(2))              ! R * cos(Theta)

                xj(1)=x1(1)*cosd(x1(3))*sind(x1(2))  ! R * cos(Phi) * sin(Theta)
                xj(2)=x1(1)*sind(x1(3))*sind(x1(2))  ! R * sin(Phi) * sin(Theta)
                xj(3)=x1(1)*cosd(x1(2))              ! R * cos(Theta)

                xr=xi-xj
             case("C","c") ! Cartesian
                xr=x1-x0
          end select
       else
          !---- Cartesian ----!
          xr=x1-x0
       end if
       dis=sqrt(dot_product(xr,xr))

       return
    End Function Distance_SC

    !!----
    !!---- Function Matrix_Phithechi(Phi,Theta,Chi,Code) Result(M)
    !!----    real(kind=sp),                intent(in) :: Phi
    !!----    real(kind=sp),                intent(in) :: Theta
    !!----    real(kind=sp),                intent(in) :: Chi
    !!----    character(len=*), optional,   intent(in) :: Code
    !!----    real(kind=sp), dimension(3,3)            :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the composition
    !!----    of a rotation around z of angle Chi, followed by a rotation of angle Theta
    !!----    around the y-axis and a subsequent rotation of angle Phi around z.
    !!----    The matrix is M = Rz(Chi) . Ry(Theta) . Rz(Chi)
    !!----    The colums represent the components of the unitary vectors {u,v,w} that
    !!----    may be considered a an alternative orthonormal frame to the canonical {i,j,k}.
    !!----    Apply the matrix M to a point in {i,j,k} gives another point in {i,j,k} obtained
    !!----    by the successive application of the three rotations given above. The transpose
    !!----    (inverse) of the M-matrix, when applied to a point in {i,j,k}, gives the coordinates
    !!----    of the same point referred to the frame {u,v,w}.
    !!----    If Code =="R" or Blank or not present then the input angles are given in radians.
    !!----    If Code =="D" then the input angles are given in degrees (Phi, Theta, Chi).
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Phithechi(Phi,Theta,Chi,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=sp),                intent(in) :: Phi
       real(kind=sp),                intent(in) :: Theta
       real(kind=sp),                intent(in) :: Chi
       character(len=*), optional,   intent(in) :: Code
       real(kind=sp), dimension(3,3)            :: Mt

       !---- Local Variables ----!
       real(kind=sp) :: p,t,c

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Phi*to_rad
                t=Theta*to_rad
                c=Chi*to_rad
             case default ! radians
                p=Phi
                t=Theta
                c=Chi
          end select
       else
          !---- radians ----!
          p=Phi
          t=Theta
          c=Chi
       end if
       Mt(1,1)= cos(p)*cos(t)*cos(c)-sin(p)*sin(c)    !
       Mt(2,1)= sin(p)*cos(t)*cos(c)+cos(p)*sin(c)    !  u
       Mt(3,1)=-sin(t)*cos(c)                         !
       Mt(1,2)=-cos(p)*cos(t)*sin(c)-sin(p)*cos(c)    !
       Mt(2,2)=-sin(p)*cos(t)*sin(c)+cos(p)*cos(c)    !  v
       Mt(3,2)= sin(t)*sin(c)                         !
       Mt(1,3)= cos(p)*sin(t)                         !
       Mt(2,3)= sin(p)*sin(t)                         !  w
       Mt(3,3)= cos(t)                                !

       return
    End Function Matrix_Phithechi

    !!----
    !!---- Function Matrix_Rx(Phi,Code) Result(M)
    !!----    real(kind=sp),                      intent(in) :: Phi
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=sp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the x-axis.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Rx(Phi,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=sp),               intent(in) :: Phi
       character(len=*), optional,  intent(in) :: Code
       real(kind=sp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=sp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Phi*to_rad
             case default ! radians
                p=Phi
          end select
       else
          !---- radians ----!
          p=Phi
       end if
       Mt(1,1)= 1.0        !
       Mt(2,1)= 0.0        !  u
       Mt(3,1)= 0.0        !
       Mt(1,2)= 0.0        !
       Mt(2,2)= cos(p)     !  v
       Mt(3,2)= sin(p)     !
       Mt(1,3)= 0.0        !
       Mt(2,3)=-sin(p)     !  w
       Mt(3,3)= cos(p)     !

       return
    End Function Matrix_Rx

    !!----
    !!---- Function Matrix_Ry(Phi,Code) Result(M)
    !!----    real(kind=sp),                      intent(in) :: Phi
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=sp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the y-axis.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Ry(Phi,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=sp),               intent(in) :: Phi
       character(len=*), optional,  intent(in) :: Code
       real(kind=sp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=sp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Phi*to_rad
             case default ! radians
                p=Phi
          end select
       else
          !---- radians ----!
          p=Phi
       end if
       Mt(1,1)= cos(p)  !
       Mt(2,1)= 0.0     !  u
       Mt(3,1)=-sin(p)  !
       Mt(1,2)= 0.0     !
       Mt(2,2)= 1.0     !  v
       Mt(3,2)= 0.0     !
       Mt(1,3)= sin(p)  !
       Mt(2,3)= 0.0     !  w
       Mt(3,3)= cos(p)  !

       return
    End Function Matrix_Ry

    !!----
    !!---- Function Matrix_Rz(Phi,Code) Result(M)
    !!----    real(kind=sp),                      intent(in) :: Phi
    !!----    character(len=*), optional,         intent(in) :: Code
    !!----    real(kind=sp), dimension(3,3)                  :: M    ! Put -> Active Rotation Matrix
    !!----
    !!----    Calculate the active rotation matrix corresponding to the positive rotation
    !!----    of an angle Phi around the z-axis.
    !!----    If Code =="R" or Blank or not present then the input angle is given in radians.
    !!----    If Code =="D" then the input angle is given in degrees.
    !!----
    !!---- Update: February - 2005
    !!
    Function Matrix_Rz(Phi,Code) Result(Mt)
       !---- Arguments ----!
       real(kind=sp),               intent(in) :: Phi
       character(len=*), optional,  intent(in) :: Code
       real(kind=sp), dimension(3,3)           :: Mt

       !---- Local Variables ----!
       real(kind=sp) :: p

       if (present(code)) then
          select case (code(1:1))
             case("D","d") ! degrees
                p=Phi*to_rad
             case default ! radians
                p=Phi
          end select
       else
          !---- radians ----!
          p=Phi
       end if
       Mt(1,1)= cos(p)  !
       Mt(2,1)= sin(p)  !  u
       Mt(3,1)= 0.0     !
       Mt(1,2)=-sin(p)  !
       Mt(2,2)= cos(p)  !  v
       Mt(3,2)= 0.0     !
       Mt(1,3)= 0.0     !
       Mt(2,3)= 0.0     !  w
       Mt(3,3)= 1.0     !

       return
    End Function Matrix_Rz

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Allocate_Coordination_Type(nasu,numops,dmax,Max_Coor)
    !!----    integer,       intent(in) :: nasu      !  In -> Number of atoms in asymmetric unit
    !!----    integer,       intent(in) :: numops    !  In -> Number of S.O. excluding lattice centerings
    !!----    real(kind=sp), intent(in) :: dmax      !  In -> Maximun distance to be calculated
    !!----    integer,      intent(out) :: Max_Coor  !  Maximum coordination allowed
    !!----
    !!----    Allocation of Coordination_Type.
    !!----    Should be called before using this module.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_Coordination_Type(nasu,numops,dmax,Max_Coor)
       !---- Arguments ----!
       integer,       intent(in) :: nasu
       integer,       intent(in) :: numops
       real(kind=sp), intent(in) :: dmax
       integer,      intent(out) :: Max_Coor

       !---- local variables ----!
       real, parameter :: r_atom=0.4 !Radius of a typical atom

       if (allocated(Coord_Info%Coord_Num)) deallocate(Coord_Info%Coord_Num)
       if (allocated(Coord_Info%N_Cooatm))  deallocate(Coord_Info%N_Cooatm)
       if (allocated(Coord_Info%N_Sym))     deallocate(Coord_Info%N_Sym)
       if (allocated(Coord_Info%Dist))      deallocate(Coord_Info%Dist)
       if (allocated(Coord_Info%S_Dist))    deallocate(Coord_Info%S_Dist)
       if (allocated(Coord_Info%Tr_Coo))    deallocate(Coord_Info%Tr_Coo)


       max_coor= (dmax/r_atom)**3
       max_coor=max(max_coor,nasu*numops)

       !---- Assigninmg the new values ----!
       Coord_Info%Natoms=nasu
       Coord_Info%Max_Coor= max_coor

       allocate (Coord_Info%Coord_Num(nasu))
       allocate (Coord_Info%N_Cooatm(nasu,max_coor))
       allocate (Coord_Info%N_Sym(nasu,max_coor))
       allocate (Coord_Info%Dist(nasu,max_coor))
       allocate (Coord_Info%S_Dist(nasu,max_coor))
       allocate (Coord_Info%Tr_Coo(3,nasu,max_coor))

       Coord_Info%Coord_Num=0
       Coord_Info%N_Cooatm =0
       Coord_Info%N_Sym    =0
       Coord_Info%Dist     =0.0
       Coord_Info%S_Dist   =0.0
       Coord_Info%Tr_Coo   =0.0

       return
    End Subroutine Allocate_Coordination_Type

    !!----
    !!---- Subroutine Allocate_Point_List(N,Pl,Ier)
    !!----    integer,               intent(in)     :: n      !  In -> Dimension for allocating components of the type
    !!----    type(point_list_type), intent(in out) :: pl     !  In Out-> Type with allocatable components
    !!----    integer,               intent(out)    :: ier    !  Out -> if ier /= 0 an error occurred.
    !!----
    !!----    Allocation of an objet of type Point_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Allocate_Point_List(N,Pl,Ier)
       !---- Arguments ----!
       integer,               intent(in)     :: n
       type(point_list_type), intent(in out) :: pl
       integer,               intent(out)    :: ier

       ier=0
       if (n <= 0) then
          ier=1
          return
       end if

       if ( .not. allocated(pl%nam) ) allocate(pl%nam(n),stat=ier)
       if ( .not. allocated(pl%p) )   allocate(pl%p(n),stat=ier)
       if ( .not. allocated(pl%x) )   allocate(pl%x(3,n),stat=ier)

       pl%nam= " "
       pl%np=0
       pl%p=0
       pl%x=0.0

       return
    End subroutine Allocate_Point_List

    !!----
    !!---- Subroutine Calc_Dist_Angle(Dmax, Dangl, Cell, Spg, A, Lun)
    !!----    real(kind=sp),            intent(in)             :: dmax   !  In -> Max. Distance to calculate
    !!----    real(kind=sp),            intent(in)             :: dangl  !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)             :: Cell   !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)             :: SpG    !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)             :: A      !  In -> Object of atom_list_type
    !!----    integer,                  optional, intent(in)   :: lun    !  In -> Logical Unit for writing
    !!----
    !!----    Subroutine to calculate distances and angles, below the prescribed distances
    !!----    "dmax" and "dangl" (angles of triplets at distance below "dangl" to an atom),
    !!----    without standard deviations. If dangl=0.0, no angle calculations are done.
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (or type atom_list, that should be allocated in the calling program).
    !!----    Writes results in file (unit=lun) if lun is present
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_Dist_Angle(Dmax, Dangl, Cell, Spg, A, Lun)
       !---- Arguments ----!
       real(kind=sp),            intent(in)   :: Dmax, Dangl
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),   intent(in)   :: A
       integer, optional,        intent(in)   :: lun

       !---- Local Variables ----!
       logical                            :: iprin
       integer                            :: i,j,k,lk,i1,i2,i3,jl,npeq,nn,L,nlines, max_coor
       character(len= 80), dimension(12)  :: texto = " "
       character(len=  5)                 :: nam,nam1,nam2
       character(len= 16)                 :: transla
       character(len=160)                 :: form3
       character(len= 90)                 :: form2= &
                                             "("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",a,""  "",3F8.4)"
       integer, dimension(3)              :: ic1,ic2
       real(kind=sp),    dimension(3)     :: xx,x1,xo,Tn,xr, QD
       real(kind=sp)                      :: T,dd, da1,da2,da12,cang12,ang12,cang1,ang2,ang1

       real(kind=sp), allocatable,dimension(:,:) :: uu
       real(kind=sp), allocatable,dimension(:,:) :: bcoo

       iprin=.false.
       if (present(lun)) then
          if (lun > 0) iprin=.true.
       end if

       call init_err_geom()

       call allocate_coordination_type(A%natoms,Spg%multip,Dmax,Max_coor)
       if(allocated(uu)) deallocate(uu)
       allocate(uu(3,Max_coor))
       if(allocated(bcoo)) deallocate(bcoo)
       allocate(bcoo(3,Max_coor))

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+1.0)
       ic1(:)=-ic2(:)
       npeq=spg%numops
       if (dangl > epsi .and. iprin ) then
          form3="(""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",f8.3/"
          form3=trim(form3)//"""         ("",a,"") :"",3f8.4,""  ("",a,"") :"",3f8.4)"
       end if

       if (spg%centred == 2) then
          npeq=2*npeq
          if (iprin) then
             write(unit=lun,fmt="(/,a)")" => Symmetry operators combined with inversion centre:"
             nlines=1
             do i=SpG%NumOps+1,npeq
                if (mod(i,2) == 0) then
                   write(unit=texto(nlines)(36:70),fmt="(a,i2,a,a)") &
                               " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                   nlines=nlines+1
                else
                   write(unit=texto(nlines)( 1:34),fmt="(a,i2,a,a)")  &
                               " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                end if
             end do
             do i=1,min(nlines,12)
                write(unit=lun,fmt="(a)") texto(i)
             end do
          end if
       end if

       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)
          nam=a%atom(i)%lab
          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") &
                       " Orig. extr. p.equiv.           Distance     tx   ty   tz       x_ext   y_ext   z_ext"
          end if
          Coord_Info%Coord_Num(i)=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             do j=1,npeq
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do i3=ic1(3),ic2(3)
                         do_jl:do jl=1,Spg%NumLat
                            Tn(1)=real(i1)+Spg%Latt_trans(1,jl)
                            Tn(2)=real(i2)+Spg%Latt_trans(2,jl)
                            Tn(3)=real(i3)+Spg%Latt_trans(3,jl)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle do_jl
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle do_jl
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            Coord_Info%Coord_Num(i)=Coord_Info%Coord_Num(i)+1

                            if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                               err_geom=.true.
                               err_mess_geom=" => Too many distances around atom: "//nam
                               return
                            end if

                            lk=lk+1
                            uu(:,lk)=x1(:)
                            Coord_Info%Dist(i,Coord_Info%Coord_Num(i))=dd
                            Coord_Info%N_Cooatm(i,Coord_Info%Coord_Num(i))=k
                            bcoo(:,Coord_Info%Coord_Num(i))=x1(:)
                            if (iprin) then
                               call Frac_Trans_1Dig(tn,transla)
                               write(unit=lun,fmt=form2) i,k,j,nam,nam1,dd,transla,x1(:)
                            end if
                         end do do_jl
                      end do !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k

          if (dangl <= epsi) cycle     !loop on "i" still running

          !---- Angle calculations for bonded atoms at distance lower than DANGL

          if (iprin) then
                write(unit=lun,fmt="(/,/,a)")       "   -------------------------------------------------------"
                write(unit=lun,fmt="(a,a,3f8.4)")   "   -  Angles around atom: ",nam, xo
                write(unit=lun,fmt="(a,/)")         "   -------------------------------------------------------"
          end if
          do j=1,Coord_Info%Coord_Num(i)
             if (Coord_Info%dist(i,j) < epsi .or. Coord_Info%dist(i,j) > dangl) cycle
             da1=Coord_Info%dist(i,j)
             i1=Coord_Info%N_Cooatm(i,j)
             nam1=a%atom(i1)%lab
             do k=j+1,Coord_Info%Coord_Num(i)
                if (Coord_Info%dist(i,k) < epsi .OR. Coord_Info%dist(i,k) > dangl) cycle
                da2=Coord_Info%dist(i,k)
                i2=Coord_Info%N_Cooatm(i,k)
                nam2=a%atom(i2)%lab
                xx(:)=bcoo(:,k)-bcoo(:,j)
                xr = matmul(Cell%Cr_Orth_cel,xx)
                da12=sqrt(dot_product(xr,xr))
                cang12=0.5_sp*(da1/da2+da2/da1-da12*da12/da1/da2)
                ang12=acosd(cang12)
                cang1=0.5_sp*(da12/da2+da2/da12-da1*da1/da12/da2)
                ang1=acosd(cang1)
                ang2=180.0_sp-ang1-ang12

                if (iprin) then
                    write(unit=lun,fmt="(/,3(a,f8.4))")  &
                         "     Atm-1   Atm-2   Atm-3            d12 =",da1,"  d23 =",da2,"   d13 =",da12
                    write(unit=lun,fmt=form3)  nam1,nam,nam2,ang12,   &
                         nam,nam2,nam1,ang1, nam,nam1,nam2,ang2,  &
                         nam1,bcoo(:,j),nam2, bcoo(:,k)
                end if
             end do !k
          end do !j
       end do !i

       return
    End Subroutine Calc_Dist_Angle

    !!----
    !!---- Subroutine Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spg, A, Lun, Lun_cons)
    !!----    real(kind=sp),            intent(in)   :: dmax     !  In -> Max. Distance to calculate
    !!----    real(kind=sp),            intent(in)   :: dangl    !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)   :: Cell     !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)   :: SpG      !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),    intent(in)   :: A        !  In -> Object of atom_list_type
    !!----    integer, optional,        intent(in)   :: lun      !  In -> Logical Unit for writing
    !!----    integer, optional,        intent(in)   :: lun_cons !  In -> Logical unit for writing restraints
    !!----
    !!----    Subroutine to calculate distances and angles, below the prescribed distances
    !!----    "dmax" and "dangl" (angles of triplets at distance below "dangl" to an atom),
    !!----    with standard deviations. If dangl=0.0, no angle calculations are done.
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (or type atom_list, that should be allocated in the calling program).
    !!----    Writes results in file (unit=lun) if iprin=.true.
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_Dist_Angle_Sigma(Dmax, Dangl, Cell, Spg, A, Lun, Lun_cons)
       !---- Arguments ----!
       real(kind=sp),            intent(in)   :: dmax, dangl
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A
       integer, optional,        intent(in)   :: lun
       integer, optional,        intent(in)   :: lun_cons

       !---- Local Variables ----!
       logical                            :: iprin
       integer,parameter                  :: nconst=500
       integer                            :: i,j,k,lk,i1,i2,i3,jl,nn,L,&
                                             itnum1,itnum2,num_const, max_coor,num_angc
       character(len=  5)                 :: nam,nam1,nam2
       character(len= 16)                 :: transla
       character(len= 20)                 :: text,tex,texton
       character(len=132)                 :: line
       character(len=160)                 :: form3
       character(len= 90)                 :: form2= &
                                             "("" "",3I4,""  ("",a ,"")-("",a ,""):"",a12,""   "",a,""  "",3F9.5)"
       integer, dimension(3)              :: ic1,ic2
       integer, dimension(192)            :: itnum
       real(kind=sp),dimension(3,3,6)     :: DerM
       real(kind=sp),    dimension(3)     :: xx,x1,xo,Tn, QD,so,ss,s1,s2,x2,tr1,tr2
       real(kind=sp)                      :: T,dd, da1,da2,da12,cang12,ang12,cang1,ang2,ang1
       real(kind=sp)                      :: sdd,sda1,sda2,sda12,sang12,sang2,sang1,srel1,srel2,srel12

       real(kind=sp), allocatable, dimension(:,:) :: uu
       real(kind=sp), allocatable, dimension(:,:) :: bcoo
       real(kind=sp), allocatable, dimension(:,:) :: sbcoo
       real(kind=sp), allocatable, dimension(:,:) :: trcoo

       character(len=132), dimension(:), allocatable  :: const_text
       character(len=132), dimension(:), allocatable  :: dist_text
       character(len=132), dimension(:), allocatable  :: angl_text
       character(len=8) :: codesym
       logical :: esta

       iprin=.false.
       if (present(lun)) then
          if (lun > 0) iprin=.true.
       end if
       call init_err_geom()
       call Allocate_Coordination_Type(A%natoms,Spg%Multip,Dmax,max_coor)

       if(allocated(uu)) deallocate(uu)
       allocate(uu(3,max_coor))
       if(allocated(bcoo)) deallocate(bcoo)
       allocate(bcoo(3,max_coor))
       if(allocated(sbcoo)) deallocate(sbcoo)
       allocate(sbcoo(3,max_coor))
       if(allocated(trcoo)) deallocate(trcoo)
       allocate(trcoo(3,max_coor))


       call get_deriv_Orth_cell(cell,DerM,"A")

       if (present(lun_cons)) then
          num_angc=0
          num_const=0
          open (unit=lun_cons, file="CFML_Restraints.tpcr", status="unknown")
          write(unit=lun_cons,fmt="(a)") " FILE with lines for soft distance and angle constraints (restraints)."
          write(unit=lun_cons,fmt="(a)") " It is intended to help editing PCR files with restraints by pasting, "
          write(unit=lun_cons,fmt="(a)") " after correcting the values as wished, to the appropriate lines.  "
          write(unit=lun_cons,fmt="(a)") " Lines with repeated identical distances have been excluded because symmetry "
          write(unit=lun_cons,fmt="(a)") " already force a hard constraint."
          write(unit=lun_cons,fmt="(a)") " Accidental coincidences have also been excluded, check that in list of distances! "
          write(unit=lun_cons,fmt="(/,a)")   " Warning! "
          write(unit=lun_cons,fmt="(a,/,a/)") " Symmetry constrained angles have not been eliminated,",&
                                            " this has to be performed by hand!"

          !---- Set ITnum ----!
          i=0
          i1=1
          i2=24
          if (spg%hexa) then
             i1=25
             i2=36
          end if
          do j=1,Spg%multip
             call searchop(SpG%Symop(j)%Rot(:,:),i1,i2,i)
             Itnum(j)=i
          end do
          if (allocated(const_text)) deallocate(const_text)
          allocate(const_text(nconst)) !Maximum number of restraints
          const_text(:)(1:132)=" "
          if (allocated(dist_text)) deallocate(dist_text)
          allocate(dist_text(nconst)) !Maximum number of restraints
          dist_text(:)(1:132)=" "
          if (allocated(angl_text)) deallocate(angl_text)
          allocate(angl_text(nconst)) !Maximum number of restraints
          angl_text(:)(1:132)=" "
       end if

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+1.5)
       ic1(:)=-ic2(:)
       if (dangl > epsi .and. iprin ) then
          form3=            "(""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""    ("",a,"")-("",a,"")-("",a,""):"",a12/"
          form3=trim(form3)//"""         ("",a,"") :"",3f9.5,""  ("",a,"") :"",3f9.5)"
       end if
       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)
          so(:)=a%atom(i)%x_std(:)
          nam=a%atom(i)%lab
          Select Case (len_trim(nam))
             case(1)
                nam="  "//trim(nam)
             case(2,3)
                nam=" "//trim(nam)
          End Select
          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") &
                  " Orig. extr. p.equiv.           Distance       tx   ty   tz        x_ext    y_ext    z_ext"
          end if
          Coord_Info%Coord_Num(i)=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             Select Case (len_trim(nam1))
               case(1)
                  nam1="  "//trim(nam1)
               case(2,3)
                  nam1=" "//trim(nam1)
             End Select
             ss(:)=A%atom(k)%x_std(:)
             do j=1,Spg%Multip
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do_i3:do i3=ic1(3),ic2(3)

                            Tn(1)=real(i1)
                            Tn(2)=real(i2)
                            Tn(3)=real(i3)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle  do_i3
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                            end do
                            call distance_and_sigma(Cell,DerM,xo,x1,so,ss,dd,sdd)
                            if (dd > dmax .or. dd < 0.001) cycle
                            Coord_Info%Coord_Num(i)=Coord_Info%Coord_Num(i)+1
                            if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                               err_geom=.true.
                               err_mess_geom=" => Too many distances around atom: "//nam
                               return
                            end if
                            lk=lk+1
                            uu(:,lk)=x1(:)

                            Coord_Info%Dist(i,Coord_Info%Coord_Num(i))=dd
                            Coord_Info%S_Dist(i,Coord_Info%Coord_Num(i))=sdd
                            Coord_Info%N_Cooatm(i,Coord_Info%Coord_Num(i))=k
                            Coord_Info%N_sym(i,Coord_Info%Coord_Num(i))=j

                            bcoo(:,Coord_Info%Coord_Num(i))=x1(:)
                            sbcoo(:,Coord_Info%Coord_Num(i))=ss(:)
                            trcoo(:,Coord_Info%Coord_Num(i))=Tn(:)
                            if (iprin) then
                               call Frac_Trans_1Dig(tn,transla)
                               call setnum_std(dd,sdd,text)
                               write(unit=lun,fmt=form2) i,k,j,nam,nam1,text,transla,x1(:)
                            end if

                            if(present(lun_cons)) then
                              esta=.false.
                              write(unit=line,fmt="(a4,tr2,a4,i5,3f10.5,tr5,2f7.4)") A%atom(i)%lab ,A%atom(k)%lab ,&
                                     Itnum(j), tn(:)+SpG%Symop(j)%tr(:) ,dd, sdd
                              if(num_const == 0) then
                                const_text(1)=line(1:132)
                                num_const=1
                                write(unit=dist_text(1),fmt="(a,2f9.5,a)") "DFIX ",dd,sdd, &
                                                                           "  "//trim(A%atom(i)%lab)//"  "//trim(A%atom(k)%lab)
                                call Write_SymTrans_Code(j,tn,codesym)
                                dist_text(1)=trim(dist_text(1))//codesym
                              else
                                do l=num_const,1,-1
                                 if( (line(1:4) == const_text(l)(1:4) .and. line(7:10) == const_text(l)(7:10)) .or. &
                                     (line(1:4) == const_text(l)(7:10) .and. line(7:10) == const_text(l)(1:4)) ) then
                                   if(line(51:132) == const_text(l)(51:132)) then
                                        esta=.true.
                                        exit
                                   end if
                                 end if
                                end do
                                if(.not. esta) then
                                  num_const=num_const+1
                                  if(num_const > NCONST) then
                                     num_const=num_const-1
                                  end if
                                  const_text(num_const)=line(1:132)
                                  write(unit=dist_text(num_const),fmt="(a,2f9.5,a)") "DFIX ",dd,sdd,&
                                        "  "//trim(A%atom(i)%lab)//"  "//trim(A%atom(k)%lab)
                                  call Write_SymTrans_Code(j,tn,codesym)
                                  dist_text(num_const)=trim(dist_text(num_const))//trim(codesym)
                                end if
                              end if
                            end if

                      end do do_i3 !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k

          if (dangl <= epsi) cycle     !loop on "i" still running

          !---- Angle calculations for bonded atoms at distance lower than DANGL
          if (present(lun_cons)) write(unit=lun_cons,fmt="(a,a)")"=> Help for possible angle restraints around atom ",A%atom(i)%lab

          if (iprin) then
             write(unit=lun,fmt="(/,/,a)")       "   -------------------------------------------------------"
             write(unit=lun,fmt="(a,a,3f8.4)")   "   -  Angles around atom: ",nam, xo
             write(unit=lun,fmt="(a,/)")         "   -------------------------------------------------------"
          end if
          do j=1,Coord_Info%Coord_Num(i)
             if (Coord_Info%Dist(i,j) < epsi .or. Coord_Info%Dist(i,j) > dangl) cycle
             da1=Coord_Info%Dist(i,j)
             sda1=Coord_Info%S_Dist(i,j)
             i1=Coord_Info%N_Cooatm(i,j)
             nam1=a%atom(i1)%lab
             Select Case (len_trim(nam1))
               case(1)
                  nam1="  "//trim(nam1)
               case(2,3)
                  nam1=" "//trim(nam1)
             End Select
             if (present(lun_cons)) then
               itnum1=itnum(Coord_Info%N_sym(i,j))
               tr1(:)=trcoo(:,j)+SpG%Symop(Coord_Info%N_sym(i,j))%tr(:)
             end if
             do k=j+1,Coord_Info%Coord_Num(i)
                if (Coord_Info%Dist(i,k) < epsi .OR. Coord_Info%Dist(i,k) > dangl) cycle
                da2=Coord_Info%Dist(i,k)
                sda2=Coord_Info%S_Dist(i,k)
                i2=Coord_Info%N_Cooatm(i,k)
                nam2=a%atom(i2)%lab
                Select Case (len_trim(nam2))
                  case(1)
                     nam2="  "//trim(nam2)
                  case(2,3)
                     nam2=" "//trim(nam2)
                End Select
                if (present(lun_cons)) then
                  itnum2=itnum(Coord_Info%N_sym(i,k))
                  tr2(:)=trcoo(:,k)+SpG%Symop(Coord_Info%N_sym(i,k))%tr(:)
                end if
                x1(:)=bcoo(:,k)
                x2(:)=bcoo(:,j)
                s1(:)=sbcoo(:,k)
                s2(:)=sbcoo(:,j)
                call distance_and_sigma(Cell,derM,x1,x2,s1,s2,da12,sda12)
                if( da12 < 0.0001) cycle

                cang12=0.5_sp*(da1/da2+da2/da1-da12*da12/da1/da2)
                ang12=ACOSd(cang12)
                cang1=0.5_sp*(da12/da2+da2/da12-da1*da1/da12/da2)
                ang1=ACOSd(cang1)
                ang2=180.0_sp-ang1-ang12

                !---- Alternative calculation of angle"s sigmas ----!
                srel1=(sda1/da1)**2
                srel12=(sda12/da12)**2
                srel2=(sda2/da2)**2
                sang12=SQRT(srel1+srel2+(sda12*da12/da1/da2)**2)*to_deg
                sang1=SQRT(srel12+srel2+(sda1*da1/da2/da12)**2)*to_deg
                sang2=SQRT(srel12+srel1+(sda2*da2/da1/da12)**2)*to_deg

                if (iprin) then
                   call setnum_std(da1,sda1,tex)
                   call setnum_std(da2,sda2,text)
                   call setnum_std(da12,sda12,texton)
                   write(unit=lun,fmt="(/,a,3a21)")  &
                        "     Atm-1   Atm-2   Atm-3           "," d12 ="//tex,"  d23 ="//text,"   d13 ="//texton
                   call setnum_std(ang12,sang12,tex)
                   call setnum_std(ang1,sang1,text)
                   call setnum_std(ang2,sang2,texton)
                   write(unit=lun,fmt=form3)  nam1,nam,nam2,tex,   &
                        nam,nam2,nam1,text, nam,nam1,nam2,texton,    &
                        nam1,bcoo(:,j),nam2, bcoo(:,k)
                end if

                if (present(lun_cons)) then

                  write(unit=lun_cons,fmt="(3(a4,tr1),i3,i4,tr1,3f8.4,tr1,3f8.4,2f7.2)") &
                  A%atom(i)%lab ,nam1 ,nam2 ,itnum1,itnum2,tr1(:),tr2(:),ang2,sang2

                  if(num_angc == 0) then
                    num_angc=1
                    line=" "
                    write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang2,sang2,&
                                                       "  "//trim(A%atom(i)%lab)//" "//trim(nam1)
                    call Write_SymTrans_Code(Coord_Info%N_sym(i,j),trcoo(:,j),codesym)
                    line=trim(line)//trim(codesym)//" "//trim(nam2)
                    call Write_SymTrans_Code(Coord_Info%N_sym(i,k),trcoo(:,k),codesym)
                    line=trim(line)//trim(codesym)
                    angl_text(1)=line(1:132)

                  else

                    line=" "
                    write(unit=line,fmt="(a,2f9.3,a)") "AFIX ",ang2,sang2,&
                                                       "  "//trim(A%atom(i)%lab)//" "//trim(nam1)
                    call Write_SymTrans_Code(Coord_Info%N_sym(i,j),trcoo(:,j),codesym)
                    line=trim(line)//trim(codesym)//" "//trim(nam2)
                    call Write_SymTrans_Code(Coord_Info%N_sym(i,k),trcoo(:,k),codesym)
                    line=trim(line)//trim(codesym)
                    esta=.false.
                    jl=index(line,"_")
                    if(jl == 0) jl=len_trim(line)
                    do l=num_angc,1,-1
                     if( line(1:jl) == angl_text(l)(1:jl)) then
                         esta=.true.
                         exit
                     end if
                    end do
                    if(.not. esta) then
                      num_angc=num_angc+1
                      if(num_angc > NCONST) num_angc=NCONST
                      angl_text(num_angc)=line(1:132)
                    end if
                  end if

                end if !present

             end do !k
          end do !j
       end do !i

       if (present(lun_cons)) then
          write(unit=lun_cons,fmt="(/,a,i5)")"=> Total number of independent distances: ",num_const
          write(unit=lun_cons,fmt="(a,/)")   "   List of possible restraints: "
          write(unit=lun_cons,fmt="(a)")" At1   At2  ITnum     T1        T2        T3          DIST   SIGMA"
          do i=1,num_const
             write(unit=lun_cons,fmt="(2x,a)") const_text(i)
          end do

          write(unit=lun_cons,fmt="(/,a)")   "   ========================================= "
          write(unit=lun_cons,fmt="(a  )")   "   List of possible restraints in CFL format "
          write(unit=lun_cons,fmt="(a,/)")   "   ========================================= "


          write(unit=lun_cons,fmt="(/a,i5)")"=> Total number of independent distance restraints: ",num_const
          do i=1,num_const
             write(unit=lun_cons,fmt="(a)") dist_text(i)
          end do
          write(unit=lun_cons,fmt="(/a,i5)")"=> Total number of possible angle restraints: ",num_angc
          do i=1,num_angc
             write(unit=lun_cons,fmt="(a)") angl_text(i)
          end do
          close(unit=lun_cons)
       end if

       return
    End Subroutine Calc_Dist_Angle_Sigma

    !!----
    !!---- Subroutine Deallocate_Coordination_Type()
    !!----
    !!----    Deallocation of Coordination_Type.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_Coordination_Type()

       if (allocated(Coord_Info%Coord_Num)) deallocate(Coord_Info%Coord_Num)
       if (allocated(Coord_Info%N_Cooatm))  deallocate(Coord_Info%N_Cooatm)
       if (allocated(Coord_Info%N_Sym))     deallocate(Coord_Info%N_Sym)
      if (allocated(Coord_Info%Dist))      deallocate(Coord_Info%Dist)
       if (allocated(Coord_Info%S_Dist))    deallocate(Coord_Info%S_Dist)
       if (allocated(Coord_Info%Tr_Coo))    deallocate(Coord_Info%Tr_Coo)

       !---- Assigninmg the new values ----!
       Coord_Info%Natoms=0
       Coord_Info%Max_Coor= 0

       return
    End Subroutine Deallocate_Coordination_Type

    !!----
    !!---- Subroutine Deallocate_Point_List(Pl)
    !!----    type(point_list_type), intent(in out) :: pl  !  In Out-> Type with allocatable components
    !!----
    !!----     De-allocation of an objet of type point_list_type
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Deallocate_Point_List(Pl)
       !---- Arguments ----!
       type(point_list_type), intent(in out) :: pl

       if (allocated(pl%nam) ) deallocate(pl%nam)
       if (allocated(pl%p) )   deallocate(pl%p)
       if (allocated(pl%x) )   deallocate(pl%x)

       return
    End Subroutine Deallocate_Point_List

    !!----
    !!---- Subroutine Distance_and_Sigma(Cellp,DerM,x0,x1,s0,s1,dis,s)
    !!----    Type(Crystal_Cell_Type),         intent(in)  :: Cellp         ! Cell object
    !!----    real(kind=sp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
    !!----    real(kind=sp), dimension(3),     intent(in)  :: x0,x1,s0,s1   ! Two points in fractional coordinates and sigmas
    !!----    real(kind=sp),                   intent(out) :: dis,s         ! Distance and sigma
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Distance_and_Sigma(Cellp,DerM,x0,x1,s0,s1,dis,s)
       !---- Arguments ----!
       Type(Crystal_Cell_Type),         intent(in)  :: Cellp         ! Cell object
       real(kind=sp), dimension(3,3,6), intent(in)  :: DerM          ! Matrix of derivatives of Cellp%Cr_Orth_cel
       real(kind=sp), dimension(3),     intent(in)  :: x0,x1,s0,s1   ! Two points in fractional coordinates and sigmas
       real(kind=sp),                   intent(out) :: dis,s         ! Distance and sigma

       !---- Local variables ----!
       integer                     :: i
       real(kind=sp), dimension(3) :: xc,xf
       real(kind=sp), dimension(6) :: dc,df

       xf=x1-x0
       xc = matmul(cellp%Cr_Orth_cel,xf)
       dis=sqrt(dot_product(xc,xc))
       do i=1,6
          dc(i) = dot_product(xc,matmul(DerM(:,:,i),xf))
       end do
       do i=1,3
          df(i) = dot_product(xc,Cellp%Cr_Orth_cel(:,i))
       end do
       df(4:6) =-df(1:3)
       s=0.0
       do i=1,3
          s = s + (dc(i)*Cellp%cell_std(i))**2
          s = s + (dc(i+3)*Cellp%ang_std(i)*to_rad)**2
          s = s + (df(i)*s1(i))**2 + (df(i+3)*s0(i))**2
       end do
       s=sqrt(s)/dis

       return
    End Subroutine Distance_and_Sigma

    !!----
    !!----  Subroutine Get_Euler_From_Fract(X1,X2,X3,Mt,Phi,Theta,Chi,Eum,Code)
    !!----    real(kind=sp),           dimension(3),   intent (in) :: x1,x2,x3
    !!----    real(kind=sp),           dimension(3,3), intent (in) :: M !Matrix transforming to Cartesian coordinates
    !!----    real(kind=sp),                           intent(out) :: theta,phi,chi
    !!----    real(kind=sp), optional, dimension(3,3), intent(out) :: EuM
    !!----    character(len=*), optional,              intent (in) :: Code
    !!----
    !!----  Subroutine to obtain the Euler angles (2nd setting) of a Cartesian frame having
    !!----  as origin the point x3, the z-axis along x1-x3 and the "xz" plane coincident with
    !!----  the plane generated by the two vectors (x2-x3,x1-x3). The
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Euler_From_Fract(X1,X2,X3,Mt,Phi,Theta,Chi,Eum,Code)
       !---- Arguments ----!
       real(kind=sp),           dimension(3),   intent (in) :: x1,x2,x3
       real(kind=sp),           dimension(3,3), intent (in) :: Mt
       real(kind=sp),                           intent(out) :: theta,phi,chi
       real(kind=sp), optional, dimension(3,3), intent(out) :: EuM
       character(len=*), optional,              intent (in) :: Code

       !---- Local variables ----!
       real(kind=sp), dimension(3)   :: u,v,w
       real(kind=sp), dimension(3,3) :: rot

!  U = ( cosPhi cosTheta cosChi - sinPhi sinChi,   sinPhi cosTheta cosChi+cosPhi sinChi,  -sinTheta cosChi)
!  V = (-sinPhi cosChi   - cosPhi cosTheta sinChi, cosPhi cosChi -sinPhi cosTheta sinChi,  sinTheta sinChi)
!  W = ( cosPhi sinTheta, sinPhi sinTheta,  cosTheta)
!
!     This corresponds to Euler angles defined in the following way:
!
!     In the starting position the cartesian frame (u,v,w) coincides with the crystallographic
!     cartesian frame (e1//a, e2 in the a-b plane and e3= e1 x e2). First a rotation Chi around
!     the e3 axis is applied, then a rotation Theta around the e2 axis and finally a rotation Phi
!     around e3. The total rotation matrix is
!
!          R(Phi,Theta,Chi) = R(e3,Phi) R(e2,Theta) R(e3,Chi) = [[ u, v, w]]
!
!     The columns of the active rotation matrix are the components of the unitary vectors u,v,w.

       w=matmul(Mt,x1-x3)
       w=w/sqrt(dot_product(w,w))
       u=matmul(Mt,x2-x3)
       u=u/sqrt(dot_product(u,u))
       v=cross_product(w,u)
       v=v/sqrt(dot_product(v,v))
       u=cross_product(v,w) !already normalized
       rot(:,1)=u; rot(:,2)=v;  rot(:,3)=w  !Matrix Rot ([u,v,w] columns)
       if (present(EuM)) EuM=rot
       if (present(Code)) then
          call get_PhiTheChi(rot,phi,theta,chi,Code)
       else
          call get_PhiTheChi(rot,phi,theta,chi)
       end if

       return
    End Subroutine Get_Euler_From_Fract

    !!----
    !!---- SUBROUTINE Get_PhiTheChi(Mt,Phi,Theta,Chi,Code)
    !!----    real(kind=sp), dimension(3,3),intent(in)  :: Mt
    !!----    real(kind=sp),                intent(out) :: Phi
    !!----    real(kind=sp),                intent(out) :: Theta
    !!----    real(kind=sp),                intent(out) :: Chi
    !!----    character(len=*), optional,   intent(in)  :: Code
    !!----
    !!----    Calculate the Euler Angles corresponding to an orthogonal matrix
    !!----    The definition of the Euler angles in this case correspond to the
    !!----    active rotation matrix obtained from the composition of a rotation
    !!----    around z of angle Chi, followed by a rotation of angle Theta
    !!----    around the y-axis and a subsequent rotation of angle Phi around z.
    !!----    The matrix is supposed to be of the form: M = Rz(Chi).Ry(Theta).Rz(Chi)
    !!----    If Code =="R" or not present then the output angles are provided in radians.
    !!----    If Code =="D" then the output angles are provided in degrees.
    !!----    A checking of the input matrix is given before calculating the angles.
    !!----    The user must check the logical variable "err_geom" after calling this
    !!----    subroutine. If err_geom=.true. it means that the input matrix is not orthogonal.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_PhiTheChi(Mt,Phi,Theta,Chi,Code)
       !---- Arguments ----!
       real(kind=sp), dimension(3,3),intent(in)  :: Mt
       real(kind=sp),                intent(out) :: Phi
       real(kind=sp),                intent(out) :: Theta
       real(kind=sp),                intent(out) :: Chi
       character(len=*), optional,   intent(in)  :: Code

       !---- Local Variables ----!
       real(kind=sp), dimension(3,3):: MTT
       real(kind=sp), parameter, dimension(3,3) :: &
                      identity = reshape ( (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),(/3,3/))

       MTT=transpose(Mt)
       MTT=matmul(MTT,Mt)-identity
       if (sum(abs(MTT)) > 5.0*eps) then
          err_geom=.true.
          err_mess_geom=" Error in Get_PhiTheChi ... the input matrix is not orthogonal! "
          return
       end if
       if (abs(Mt(3,3)-1.0) < eps) then  !M(3,3)=cos(Theta)
          Theta=0.0
          Phi=0.0
          Chi=acos(Mt(1,1))               !M(1,1)=cos(Phi)cos(Theta)cos(Chi)-sin(Phi)sin(Chi)
       else if(abs(Mt(3,3)+1.0) < eps) then
          Theta=pi
          Phi=0.0
          Chi=acos(-Mt(1,1))
       else
          Theta=acos(Mt(3,3))
          Phi=atan2(Mt(2,3),Mt(1,3))     !M(1,3)=cos(Phi)sin(Theta)  M(2,3)=sin(phi)sin(Theta)
          Chi=atan2(Mt(3,2),-Mt(3,1))    !M(3,1)= -sin(Theta)cos(Chi)   M(3,2)= sin(Theta)sin(Chi)
       end if
       if (present(Code)) then
          if (code(1:1)=="D" .or. code(1:1)=="d") then
             Phi=Phi*to_deg
             Theta=Theta*to_deg
             Chi=Chi*to_deg
          end if
       end if

       return
    End Subroutine Get_PhiTheChi

    !!----
    !!---- Subroutine Get_Transf_List(Trans,Ox,Pl,Npl,Ifail)
    !!----   real,            dimension(3,3), intent(in) :: trans   !Matrix transforming the basis
    !!----   real,            dimension(3  ), intent(in) :: ox      !Coordinates of origin of the new basis
    !!----   type(point_list_type),           intent(in) :: pl      !Input List of points
    !!----   type(point_list_type),       intent(in out) :: npl     !Output list of transformed points
    !!----   integer,                         intent(out):: ifail   !If ifail/=0 matrix inversion failed
    !!----
    !!----  Subroutine to get the fractional coordinates of the points of the input list "pl" in the
    !!----  new transformed cell ( a'= trans a) displaced to the new origing "ox". The coordinates
    !!----  are generated using only lattice translations. All coordinates are reduced to be
    !!----  between 0.0 and 1.0, so that  0.0 <= x,y,z < 1.0
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Transf_List(trans,ox,pl,npl,ifail)
       !---- Arguments ----!
       real,            dimension(3,3), intent(in) :: trans
       real,            dimension(3  ), intent(in) :: ox
       type(point_list_type),           intent(in) :: pl
       type(point_list_type),       intent(in out) :: npl
       integer,                         intent(out):: ifail

       !---- local variables ----!
       integer                 :: i,j,ia,ib,ic,nat,mm
       integer, dimension(3)   :: mini,maxi
       real,    dimension(7,3) :: vecpar
       real,    dimension(3,3) :: si
       real,    dimension(3  ) :: xx, xxn,v

       ifail=0
       call matrix_inverse(trans,si,ifail)
       if (ifail == 1) return

       !----  Construction of the 7 vertices of the new cell
       !----  1:a, 2:b, 3:c, 4:a+b, 5:a+c, 6:b+c 7:a+b+c
       do j=1,3
          do i=1,3
             vecpar(i,j)=trans(i,j)
          end do
          vecpar(4,j)=trans(1,j)+trans(2,j)
          vecpar(5,j)=trans(1,j)+trans(3,j)
          vecpar(6,j)=trans(2,j)+trans(3,j)
          vecpar(7,j)=trans(1,j)+trans(2,j)+trans(3,j)
       end do

       !---- Exploration of the vertex matrix
       mini(:)=1000
       maxi(:)=-1000
       do j=1,3
          do i=1,7
             if (vecpar(i,j) < mini(j)) mini(j)=nint(min(vecpar(i,j),0.0))
             if (vecpar(i,j) > maxi(j)) maxi(j)=nint(max(vecpar(i,j),1.0))
          end do
       end do

       !
       !   Explore the region  a-> min(1)---max(1)  where atoms will be generated
       !                       b-> min(2)---max(2)
       !                       c-> min(3)---max(3)
       !   and select those belonging to the interior of the new cell before
       !   translation to the new origin.
       !   set the translation to the new origin, put the atoms inside the new
       !   unit cell and, finally, print atoms coordinates
       !
       nat=0
       do mm=1,pl%np
          do ia=mini(1),maxi(1)
             xx(1)=pl%x(1,mm)+real(ia)
             do ib=mini(2),maxi(2)
                xx(2)=pl%x(2,mm)+real(ib)
                do_ic: do ic=mini(3),maxi(3)
                   xx(3)=pl%x(3,mm)+real(ic)
                   xxn=matmul(xx-ox,si)
                   xxn=Modulo_Lat(xxn)
                   do i=nat,1,-1
                      v=npl%x(:,i)-xxn(:)
                      if (Lattice_trans(v,"P") ) cycle do_ic
                   end do
                   nat=nat+1
                   npl%x(:,nat)= xxn
                   if ( nat < 10) then
                      write(unit=npl%nam(nat),fmt="(a,i1)") trim(pl%nam(mm))//"_",nat
                   else if( nat < 100) then
                      write(unit=npl%nam(nat),fmt="(a,i2)") trim(pl%nam(mm))//"_",nat
                   else
                      write(unit=npl%nam(nat),fmt="(a,i3)") trim(pl%nam(mm))//"_",nat
                   end if
                end do do_ic
             end do
          end do
       end do
       npl%np=nat

       return
    End Subroutine Get_Transf_List

    !!----
    !!---- SUBROUTINE INIT_ERR_GEOM()
    !!----
    !!----    Initialize the errors flags in Geom_Calculations
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Geom()

       err_geom=.false.
       err_mess_geom=" "

       return
    End Subroutine Init_Err_Geom

    !!----
    !!---- Subroutine P1_Dist(Dmax, Cell, Spg, Ac, Lun)
    !!----    real(kind=sp),            intent(in)    :: dmax      !  In -> Max. Distance to be calculated
    !!----    type (Crystal_cell_Type), intent(in)    :: Cell      !  In -> Object of Crystal_cell_Type
    !!----    type (Space_Group_Type),  intent(in)    :: SpG       !  In -> Object of Space_Group_Type
    !!----    type (Atoms_Cell_Type),   intent(in out):: Ac        !  In -> Object of Atoms_Cell_Type
    !!----                                                           Out -> Updated Object of Atoms_Cell_Type
    !!----    integer,optional,         intent(in)    :: lun       !  In -> Logical Unit for writing
    !!----
    !!----    Subroutine calculate distances, below the prescribed distances "dmax",
    !!----    without standard deviations. No symmetry is applied: only lattice translations.
    !!----    Need as input the objects "Cell" (of type Crystal_cell_type), "SpG" (of type Space_Group_Type)
    !!----    and "Ac" (or type Atoms_Cell). Complete the construction of Ac.
    !!----    Control for error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine P1_Dist(Dmax, Cell, Spg, Ac, Lun)
       !---- Arguments ----!
       real(kind=sp),            intent(in)       :: dmax
       type (Crystal_cell_Type), intent(in)       :: Cell
       type (Space_Group_Type),  intent(in)       :: SpG
       type (Atoms_Cell_Type),   intent(in out)   :: Ac
       integer, optional,        intent(in)       :: lun

       !---- Local Variables ----!
       logical                                :: iprint
       character(len=6 )                      :: nam,nam1
       character(len=16)                      :: transla
       character(len=90)                      :: form2= &
                       "("" "",2I4,""   ("",a,"")-("",a,""):"",f10.4,""   "",a,""  "",3F8.4)"
       integer                                :: i,k,lk,i1,i2,i3,jl,nn,L,inew,ne,id
       integer, dimension(3)                  :: ic1,ic2
       integer, dimension(Ac%nat,Ac%nat)      :: mn  !neighbouring matrix
       real(kind=sp)                          :: T,dd
       real(kind=sp), dimension(3)            :: xx,x1,xo,Tn,xr, QD
       real(kind=sp), dimension(3,Ac%nat*Ac%nat*spg%multip) :: u

       iprint=.false.
       if (present(lun)) then
          if (lun > 0) iprint=.true.
       end if
       call init_err_geom()
       id=3*nint(0.74048*(dmax/1.1)**3)

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+3.0)
       ic1(:)=-ic2(:)
       mn(:,:) = 0
       inew=0
       do i=1,ac%nat
          xo(:)=Ac%xyz(:,i)
          nam= ac%noms(i)
          if (iprint) then
             write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                       "    Distances less than",dmax,"  to atom: ",nam, xo(:)
             write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
             write(unit=lun,fmt="(/,/,a,/,/)") &
                       " Orig. extr.                    Distance     tx   ty   tz       x_ext   y_ext   z_ext"
          end if
          ne=0
          do k=1,Ac%nat
             lk=1
             u(:,lk)=xo(:)
             xx(:)=Ac%xyz(:,k)
             nam1= Ac%noms(k)
             do i1=ic1(1),ic2(1)
                do i2=ic1(2),ic2(2)
                   do i3=ic1(3),ic2(3)
                      do_jl:do jl=1,Spg%NumLat
                         Tn(1)=real(i1)+Spg%Latt_trans(1,jl)
                         Tn(2)=real(i2)+Spg%Latt_trans(2,jl)
                         Tn(3)=real(i3)+Spg%Latt_trans(3,jl)
                         x1(:)=xx(:)+tn(:)
                         do l=1,3
                            t=abs(x1(l)-xo(l))*qd(l)
                            if (t > dmax) cycle do_jl
                         end do
                         do nn=1,lk
                            if (sum(abs(u(:,nn)-x1(:)))  <= epsi) cycle do_jl
                         end do
                         xr = matmul(cell%cr_orth_cel,x1-xo)
                         dd=sqrt(dot_product(xr,xr))
                         if (dd > dmax .or. dd < 0.001) cycle
                         lk=lk+1
                         u(:,lk)=x1(:)
                         call Frac_Trans_1Dig(tn,transla)
                         if (iprint) write(unit=lun,fmt=form2)i,k,nam,nam1,dd,transla,x1(:)
                         mn(i,k)=mn(i,k)+1
                         ne=ne+1
                         IF (ne > id) THEN
                            err_geom=.true.
                            err_mess_geom="Too many connected atoms! in sub. P1_dist"
                            return
                         END IF
                         Ac%neighb_atom(i,ne)=k    !Pointer to the number of atom connected to i
                         Ac%distance   (i,ne)=dd   !Corresponding distance
                         Ac%trans(:,i,ne)=tn(:)    !corresponding lattice translation
                         do nn=1,inew
                            if (abs(dd-ac%ddist(nn)) <= epsi) then
                               if (equiv_atm(nam,nam1,ac%ddlab(nn)))  cycle do_jl
                            end if
                         end do
                         inew=inew+1
                         Ac%ddist(inew)=dd
                         Ac%ddlab(inew)=wrt_lab(nam,nam1)
                      end do do_jl
                   end do !i3
                end do !i2
             end do !i1
          end do !k
          Ac%neighb(i)=ne
       end do !i
       Ac%ndist=inew
       if (iprint) then
          write(unit=lun,fmt="(/,/,a)") " -------------------"
          write(unit=lun,fmt="(a)"  ) " Neighbouring matrix"
          write(unit=lun,fmt="(a)")   " -------------------"
          write(unit=lun,fmt="(a)")
          write(unit=lun,fmt="(24i3)")"     ",(i,i=1,Ac%nat)
          write(unit=lun,fmt="(a)")
          do i=1,ac%nat
             write(unit=lun,fmt="(i3,a,24i3)")i,"  ",(mn(i,k),k=1,Ac%nat)
          end do
          write(unit=lun,fmt="(a,/,/,/)")
       end if

       return
    End Subroutine P1_Dist

    !!----
    !!---- Subroutine Print_Distances(Lun, Dmax, Cell, Spg, A)
    !!----    integer,                  intent(in)   :: lun    !  In -> Logical Unit for writing
    !!----    real(kind=sp),            intent(in)   :: dmax   !  In -> Max. Distance to be calculated
    !!----    type (Crystal_cell_Type), intent(in)   :: Cell   !  In -> Object of Crystal_cell_Type
    !!----    type (Space_Group_Type),  intent(in)   :: SpG    !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)   :: A      !  In -> Object of atom_list_type
    !!----
    !!----    Subroutine to print distances, below the prescribed distances
    !!----    "dmax", without standard deviations.
    !!----    Need as input the objects "Cell" (of type Crystal_cell_type), "SpG"
    !!----    (of type Space_Group_type) and "A" (or type atom_list_type, that should be
    !!----    allocated in the calling program).
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Print_Distances(Lun, Dmax, Cell, Spg, A)
       !-- Arguments --!
       integer,                  intent(in)   :: lun
       real(kind=sp),            intent(in)   :: dmax
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),   intent(in)   :: A

       !---- Local Variables ----!
       integer                           :: i,j,k,lk,i1,i2,i3,jl,npeq,nn,L,nlines
       character(len=80), dimension(12)  :: texto=" "
       character(len=5 )                 :: nam,nam1
       character(len=16)                 :: transla
       character(len=54)                 :: form2= &
                                            "("" "",3I4,""  ("",a,"")-("",a,""):"",f9.4,""   "",a,""  "",3F8.4)"
       integer,          dimension(3)    :: ic1,ic2
       real(kind=sp),    dimension(3)    :: xx,x1,xo,Tn,xr, QD
       real(kind=sp)                     :: T,dd
       real(kind=sp), dimension(3,A%Natoms*Spg%multip) :: uu

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+1.0)
       ic1(:)=-ic2(:)
       npeq=spg%numops

       if (Spg%Centred == 2) then
          npeq=2*npeq
          write(unit=lun,fmt="(a)")" => Symmetry operators combined with inversion centre:"
          nlines=1
          do i=SpG%NumOps+1,npeq
             if (mod(i,2) == 0) then
                write(unit=texto(nlines)(36:70),fmt="(a,i2,a,a)") &
                                           " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:34),fmt="(a,i2,a,a)")  &
                                           " => SYMM(",i,"): ",trim(SpG%SymopSymb(i))
             end if
          end do
          do i=1,min(nlines,12)
             write(unit=lun,fmt="(a)") texto(i)
          end do
       end if

       do i=1,A%natoms
          nam=a%atom(i)%lab
          xo(:)=a%atom(i)%x(:)
          write(unit=lun,fmt="(/,/,a)")"    -------------------------------------------------------------------"
          write(unit=lun,fmt="(a,f8.4,a,a,3f8.4)")   &
                    "    Distances less than",dmax,"  to atom: ",nam, xo(:)
          write(unit=lun,fmt="(a,/,/)")"    -------------------------------------------------------------------"
          write(unit=lun,fmt="(/,/,a,/,/)") &
                    " Orig. extr. p.equiv.           Distance     tx   ty   tz       x_ext   y_ext   z_ext"
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             nam1=a%atom(k)%lab
             do j=1,npeq
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do i3=ic1(3),ic2(3)
                         do_jl:do jl=1,Spg%NumLat
                            Tn(1)=real(i1)+Spg%Latt_trans(1,jl)
                            Tn(2)=real(i2)+Spg%Latt_trans(2,jl)
                            Tn(3)=real(i3)+Spg%Latt_trans(3,jl)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle do_jl
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle do_jl
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            lk=lk+1
                            uu(:,lk)=x1(:)
                            call Frac_Trans_1Dig(tn,transla)
                            write(unit=lun,fmt=form2)i,k,j,nam ,nam1,dd,transla,x1(:)
                         end do do_jl
                      end do !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k
       end do !i

       return
    End Subroutine Print_Distances

    !!----
    !!---- Subroutine Set_Orbits_Inlist(Spg,Pl)
    !!----    type(space_group_type), intent(in)     :: SpG     !  In -> Space group
    !!----    type(point_list_type),  intent(in out) :: pl      !  In -> list of points
    !!----
    !!----    Set up of the integer pointer "pl%p" in the object "pl" of type point_list_type.
    !!----    Each point is associated with the number of an orbit. This pointer is useful
    !!----    to get the asymmetric unit with respect to the input space group of an arbitrary
    !!----    list of points (atom coordinates).
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Orbits_Inlist(Spg,Pl)
       !---- Arguments ----!
       type(space_group_type), intent(in)     :: SpG
       type(point_list_type),  intent(in out) :: pl

       !--- Local variables ---!
       integer               :: i,j,norb,nt
       real,    dimension(3) :: x,xx,v

       norb=0
       pl%p=0
       do i=1,pl%np
          if (pl%p(i) == 0) then
             norb=norb+1
             pl%p(i)=norb
             x=pl%x(:,i)
             do j=1,Spg%multip
                xx=ApplySO(Spg%SymOp(j),x)
                xx=modulo_lat(xx)
                do nt=1,pl%np
                   if (pl%p(nt) /= 0) cycle
                   v=pl%x(:,nt)-xx(:)
                   if (Lattice_trans(v,Spg%spg_lat)) pl%p(nt)=norb
                end do
             end do
          end if
       end do

       return
    End Subroutine Set_Orbits_Inlist

    !!----
    !!---- Subroutine Set_TDist_Coordination(Max_coor,Dmax, Cell, Spg, A)
    !!----    integer,                  intent(in)   :: max_coor !  Maximum expected coordination
    !!----    real(kind=sp),            intent(in)   :: dmax     !  In -> Max. Distance to calculate
    !!----    real(kind=sp),            intent(in)   :: dangl    !  In -> Max. distance for angle calculations
    !!----    type (Crystal_cell_type), intent(in)   :: Cell     !  In -> Object of Crytal_Cell_Type
    !!----    type (Space_Group_type),  intent(in)   :: SpG      !  In -> Object of Space_Group_Type
    !!----    type (atom_list_type),   intent(in)    :: A        !  In -> Object of atom_list_type
    !!----
    !!----    Subroutine to calculate distances, below the prescribed distance "dmax"
    !!----    Sets up the coordination type: Coord_Info for each atom in the asymmetric unit
    !!----    Needs as input the objects Cell (of type Crystal_cell), SpG (of type Space_Group)
    !!----    and A (or type atom_list, that should be allocated in the calling program).
    !!----    The input argument Max_Coor is obtained, before calling the present procedure,
    !!----    by a call to Allocate_Coordination_Type with arguments:(A%natoms,Spg%Multip,Dmax,max_coor)
    !!----    Further calls to this routine do not need a previous call to Allocate_Coordination_Type.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_TDist_Coordination(max_coor,Dmax, Cell, Spg, A)
       !---- Arguments ----!
       integer,                  intent(in)   :: max_coor
       real(kind=sp),            intent(in)   :: dmax
       type (Crystal_cell_Type), intent(in)   :: Cell
       type (Space_Group_Type),  intent(in)   :: SpG
       type (atom_list_type),    intent(in)   :: A

       !---- Local Variables ----!
       integer                              :: i,j,k,lk,i1,i2,i3,nn,L
       integer,       dimension(3)          :: ic1,ic2
       real(kind=sp), dimension(3)          :: xx,x1,xo,Tn,xr, QD
       real(kind=sp)                        :: T,dd
       real(kind=sp), dimension(3,max_coor) :: uu

      ! call init_err_geom()  !Control of error

       qd(:)=1.0/cell%rcell(:)
       ic2(:)= nint(dmax/cell%cell(:)+1.5)
       ic1(:)=-ic2(:)
       do i=1,a%natoms
          xo(:)=a%atom(i)%x(:)
          Coord_Info%Coord_Num(i)=0
          do k=1,a%natoms
             lk=1
             uu(:,lk)=xo(:)
             do j=1,Spg%Multip
                xx=ApplySO(Spg%SymOp(j),a%atom(k)%x)
                do i1=ic1(1),ic2(1)
                   do i2=ic1(2),ic2(2)
                      do_i3:do i3=ic1(3),ic2(3)
                            Tn(1)=real(i1); Tn(2)=real(i2); Tn(3)=real(i3)
                            x1(:)=xx(:)+tn(:)
                            do l=1,3
                               t=abs(x1(l)-xo(l))*qd(l)
                               if (t > dmax) cycle  do_i3
                            end do
                            do nn=1,lk
                               if (sum(abs(uu(:,nn)-x1(:)))  <= epsi) cycle  do_i3
                            end do
                            xr = matmul(cell%cr_orth_cel,x1-xo)
                            dd=sqrt(dot_product(xr,xr))
                            if (dd > dmax .or. dd < 0.001) cycle
                            Coord_Info%Coord_Num(i)=Coord_Info%Coord_Num(i)+1
                           ! Control not performed ... it is supposed that max_coor is large enough
                           !if (Coord_Info%Coord_Num(i) > Coord_Info%Max_Coor) then
                           !   err_geom=.true.
                           !   err_mess_geom=" => Too many distances around an atom"
                           !   return
                           !end if
                            lk=lk+1
                            uu(:,lk)=x1(:)
                            Coord_Info%Dist(i,Coord_Info%Coord_Num(i))=dd
                            Coord_Info%N_Cooatm(i,Coord_Info%Coord_Num(i))=k
                            Coord_Info%N_sym(i,Coord_Info%Coord_Num(i))=j
                      end do do_i3 !i3
                   end do !i2
                end do !i1
             end do !j
          end do !k
       end do !i

       return
    End Subroutine Set_TDist_Coordination

 End Module Geom_Calculations
