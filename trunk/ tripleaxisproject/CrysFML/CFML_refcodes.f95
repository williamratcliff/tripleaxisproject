!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: REFINEMENT_CODES
!!----   INFO: Refinable Codes for Parameters
!!----
!!---- HISTORY
!!----    Update: March - 2005
!!----
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    ERR_REFCODES
!!----    ERR_MESS_REFCODES
!!--++    NCODE                        [Private]
!!--++    CODE_NAM                     [Private]
!!--++    NKEY                         [Private]
!!--++    KEY_CODE                     [Private]
!!----    NP_CONS
!!----    NP_MAX
!!----    NP_REFI
!!----    NP_REST
!!----    V_BCON
!!----    V_BOUNDS
!!----    V_LIST
!!----    V_NAME
!!----    V_VEC
!!----    V_SHIFT
!!----
!!---- PUBLIC PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATING RESTPARAM
!!----       ALLOCATING V_PARAM
!!--++       DELETE_REFCODES             [Private]
!!--++       DELETE_REFCODES_FATOM       [Overloaded]
!!--++       DELETE_REFCODES_MOLCRYS     [Overloaded]
!!--++       DELETE_REFCODES_MOLEC       [Overloaded]
!!--++       FILL_REFCODES               [Private]
!!--++       FILL_REFCODES_FATOM         [Overloaded]
!!--++       FILL_REFCODES_MOLCRYS       [Overloaded]
!!--++       FILL_REFCODES_MOLEC         [Overloaded]
!!--++       GET_ATOMBET_CTR             [Private]
!!--++       GET_ATOMPOS_CTR             [Private]
!!--++       GET_CONCODES_LINE           [Private]
!!--++       GET_CONCODES_LINE_FATOM     [Overloaded]
!!--++       GET_CONCODES_LINE_MOLCRYS   [Overloaded]
!!--++       GET_CONCODES_LINE_MOLEC     [Overloaded]
!!--++       GET_REFCODES_LINE           [Private]
!!--++       GET_REFCODES_LINE_FATOM     [Overloaded]
!!--++       GET_REFCODES_LINE_MOLCRYS   [Overloaded]
!!--++       GET_REFCODES_LINE_MOLEC     [Overloaded]
!!----       INIT_ERR_REFCODES
!!----       INIT_REFCODES
!!--++       INIT_REFCODES_FATOM         [Overloaded]
!!--++       INIT_REFCODES_MOLCRYS       [Overloaded]
!!--++       INIT_REFCODES_MOLEC         [Overloaded]
!!----       READ_REFCODES_FILE
!!--++       READ_REFCODES_FILE_FATOM    [Overloaded]
!!--++       READ_REFCODES_FILE_MOLCRYS  [Overloaded]
!!--++       READ_REFCODES_FILE_MOLEC    [Overloaded]
!!--++       SPLIT_OPERATIONS            [Private]
!!----       VSTATE_TO_ATOMSPAR
!!--++       VSTATE_TO_ATOMSPAR_FATOM    [Overloaded]
!!--++       VSTATE_TO_ATOMSPAR_MOLCRYS  [Overloaded]
!!--++       VSTATE_TO_ATOMSPAR_MOLEC    [Overloaded]
!!----       WRITE_INFO_REFCODES
!!--++       WRITE_INFO_REFCODES_FATOM   [Overloaded]
!!--++       WRITE_INFO_REFCODES_MOLCRYS [Overloaded]
!!--++       WRITE_INFO_REFCODES_MOLEC   [Overloaded]
!!----       WRITE_RESTRAINTS_OBSCALC
!!--++
!!----
!!
 Module Refinement_Codes
    !---- Modules ----!
    Use Math_Gen,                  only: Sp, Sort
    Use String_Utilities,          only: Cutst, U_Case, L_Case, Getword, GetNum
    Use Crystallographic_Symmetry, only: Space_Group_Type, Get_Stabilizer, Symmetry_Symbol,   &
                                         Sym_B_Relations, Read_SymTrans_Code, Get_SymSymb
    Use Atom_Module,               only: Atom_list_Type  !, Atom_Type
    Use Molecular_Crystals,        only: Molecule_Type, Molecular_Crystal_Type
    Use IO_Formats,                only: File_List_Type

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Allocate_VParam, Init_RefCodes, Read_RefCodes_File, VState_to_AtomsPar,  &
              Write_Info_RefCodes, Get_RestAng_Line, Get_RestDis_Line, Get_RestTor_Line, &
              Allocate_RestParam, Write_Restraints_ObsCalc, Init_Err_RefCodes

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private :: Delete_RefCodes,   &
               Delete_RefCodes_FAtom, Delete_RefCodes_Molcrys, Delete_RefCodes_Molec,          &
               Fill_RefCodes,     &
               Fill_RefCodes_FAtom, Fill_RefCodes_Molcrys, Fill_RefCodes_Molec,                &
               Get_AtomBet_Ctr, Get_Atompos_Ctr,                                               &
               Get_ConCodes_Line, &
               Get_ConCodes_Line_FAtom, Get_ConCodes_Line_Molcrys, Get_ConCodes_Line_Molec,    &
               Get_RefCodes_Line, &
               Get_RefCodes_Line_FAtom, Get_RefCodes_Line_Molcrys, Get_RefCodes_Line_Molec,    &
               Init_RefCodes_FAtom, Init_RefCodes_Molcrys, Init_RefCodes_Molec,                &
               Read_RefCodes_File_FAtom, Read_RefCodes_File_Molcrys, Read_RefCodes_File_Molec, &
               Split_Operations,  &
               VState_to_AtomsPar_FAtom, VState_to_AtomsPar_Molcrys, VState_to_AtomsPar_Molec, &
               Write_Info_RefCodes_FAtom, Write_Info_RefCodes_Molcrys, Write_Info_RefCodes_Molec

    !---- Definitions ----!

    !!----
    !!---- TYPE :: ANGLE_RESTRAINT_TYPE
    !!----
    !!---- Type, public :: Angle_Restraint_Type
    !!----    real(kind=sp)                 :: AObs
    !!----    real(kind=sp)                 :: ACalc
    !!----    real(kind=sp)                 :: Sigma
    !!----    integer,dimension(3)          :: P
    !!----    character(len=8),dimension(2) :: STCode
    !!---- End Type Angle_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Angle_Restraint_Type
      real(kind=sp)                 :: AObs
      real(kind=sp)                 :: ACalc
      real(kind=sp)                 :: Sigma
      integer,dimension(3)          :: P
      character(len=8),dimension(2) :: STCode
    End Type Angle_Restraint_Type

    !!----
    !!---- ANG_REST
    !!----    type(Angle_Restraint_Type), public, dimension(:), allocatable :: Ang_rest
    !!----
    !!---- Relations for Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Angle_Restraint_Type), public, dimension(:), allocatable :: Ang_Rest

    !!----
    !!---- TYPE :: DISTANCE_RESTRAINT_TYPE
    !!----
    !!---- Type, public :: Distance_Restraint_Type
    !!----    real(kind=sp)        :: DObs
    !!----    real(kind=sp)        :: DCalc
    !!----    real(kind=sp)        :: Sigma
    !!----    integer,dimension(2) :: P
    !!----    character(len=8)     :: STCode    ! _N.ABC
    !!---- End Type Distance_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Distance_Restraint_Type
      real(kind=sp)        :: DObs
      real(kind=sp)        :: DCalc
      real(kind=sp)        :: Sigma
      integer,dimension(2) :: P
      character(len=8)     :: STCode
    End Type Distance_Restraint_Type

    !!----
    !!---- DIS_REST
    !!----    type(Distance_Restraint_Type), public, dimension(:), allocatable :: Dis_Rest
    !!----
    !!---- Relations for Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Distance_Restraint_Type), public, dimension(:), allocatable :: Dis_Rest

    !!----
    !!---- TYPE :: TORSION_RESTRAINT_TYPE
    !!----
    !!---- Type, public :: Torsion_Restraint_Type
    !!----    real(kind=sp)                 :: TObs
    !!----    real(kind=sp)                 :: TCalc
    !!----    real(kind=sp)                 :: Sigma
    !!----    integer,dimension(4)          :: P
    !!----    character(len=8),dimension(3) :: STCode
    !!---- End Type Torsion_Restraint_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Torsion_Restraint_Type
      real(kind=sp)                 :: TObs
      real(kind=sp)                 :: TCalc
      real(kind=sp)                 :: Sigma
      integer,dimension(4)          :: P
      character(len=8),dimension(3) :: STCode
    End Type Torsion_Restraint_Type

    !!----
    !!---- TOR_REST
    !!----    type(Torsion_Restraint_Type), public, dimension(:), allocatable :: Tor_Rest
    !!----
    !!---- Relations for Torsion Angle Restraints
    !!----
    !!---- Update: March - 2005
    !!
    type(Torsion_Restraint_Type), public, dimension(:), allocatable :: Tor_Rest


    !!----
    !!---- ERR_REFCODES
    !!----    logical, public :: Err_RefCodes
    !!----
    !!----    Error variable
    !!----
    !!---- Update: March - 2005
    !!
    logical, public :: Err_RefCodes = .false.

    !!----
    !!---- ERR_MESS_REFCODES
    !!----    character(len=150), public :: Err_Mess_RefCodes
    !!----
    !!----    Error variable messages
    !!----
    !!---- Update: March - 2005
    !!
    character(len=150), public :: Err_Mess_RefCodes = " "

    !!--++
    !!--++ NCODE
    !!--++    integer, private :: NCode
    !!--++
    !!--++    Number of Code variables
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private, parameter :: NCode=21

    !!----
    !!---- CODE_NAM
    !!----    character(len=6), dimension(NCode), public, parameter :: Code_Nam
    !!----
    !!----    Variable for treatement Codes
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), dimension(NCode), public, parameter :: Code_Nam=(/ "X_    ","Y_    ","Z_    ",       &
                                                                         "Biso_ ","Occ_  ","B11_  ",       &
                                                                         "B22_  ","B33_  ","B12_  ",       &
                                                                         "B13_  ","B23_  ","Banis_",       &
                                                                         "Xc_   ","Yc_   ","Zc_   ",       &
                                                                         "Theta_","Phi_  ","Chi_  ",       &
                                                                         "Th_L_ ","Th_T_ ","Th_S_ "/)
    !!--++
    !!--++ NKEY
    !!--++    integer, public :: NKey
    !!--++
    !!--++    Number of Keys variables
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private, parameter :: NKey=8

    !!----
    !!---- KEY_CODE
    !!----    character(len=3), dimension(NKey), public, parameter :: key_Code
    !!----
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), dimension(Nkey), public, parameter :: Key_Code=(/ "XYZ","OCC","BIS","BAN","ALL", &
                                                                        "CEN","ORI","THE"/)
    !!----
    !!---- NP_CONS
    !!----    integer, public :: NP_Cons
    !!----
    !!----    Number of Constraints relations
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Cons

    !!----
    !!---- NP_MAX
    !!----    integer, public :: NP_Max
    !!----
    !!----    Number of Maximum Parameters to Refine
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Max

    !!----
    !!---- NP_REFI
    !!----    integer, public :: NP_Refi
    !!----
    !!----    Number of Refinable Parameters
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Refi

    !!----
    !!---- NP_Rest_Ang
    !!----    integer, public :: NP_Rest_Ang
    !!----
    !!----    Number of Angle Restraints relations
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Rest_Ang
    !!----
    !!---- NP_Rest_Dis
    !!----    integer, public :: NP_Rest_Dis
    !!----
    !!----    Number of Distance Restraints relations
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Rest_Dis

    !!----
    !!---- NP_Rest_Tor
    !!----    integer, public :: NP_Rest_Tor
    !!----
    !!----    Number of Torsion Restraints relations
    !!----
    !!---- Update: March - 2005
    !!
    integer, public :: NP_Rest_Tor

    !!----
    !!---- V_BCon
    !!----    integer, public, dimension(:), allocatable :: V_BCon
    !!----
    !!----    Vector of Boundary Conditions
    !!----
    !!---- Update: March - 2005
    !!
    integer, public, dimension(:), allocatable :: V_BCon

    !!----
    !!---- V_Bounds
    !!----    real, public, dimension(:,:), allocatable :: V_Bounds
    !!----
    !!----    Vector of Lower, Upper limits and Step for Parameters
    !!----
    !!---- Update: March - 2005
    !!
    real, public, dimension(:,:),  allocatable :: V_Bounds

    !!----
    !!---- V_List
    !!----    integer, public, dimension(:), allocatable :: V_List
    !!----
    !!----    Vector of Index point the atom order
    !!----
    !!---- Update: March - 2005
    !!
    integer, public, dimension(:),  allocatable :: V_List

    !!----
    !!---- V_Name
    !!----    character(len=20), public, dimension(:), allocatable :: V_Name
    !!----
    !!----    Vector of  Name of Refinable Parameters
    !!----
    !!---- Update: March - 2005
    !!
    character(len=20), public, dimension(:), allocatable :: V_Name

    !!----
    !!---- V_Vector
    !!----    real, public, dimension(:), allocatable :: V_Vec
    !!----
    !!----    Vector of  Parameters
    !!----
    !!---- Update: March - 2005
    !!
    real, public, dimension(:),    allocatable :: V_Vec

    !!----
    !!---- V_Shift
    !!----    real, public, dimension(:), allocatable :: V_Shift
    !!----
    !!----    Vector of holding the shift of parameters
    !!----
    !!---- Update: March - 2005
    !!
    real, public, dimension(:),    allocatable :: V_Shift

    !---- Interfaces - Overloaded ----!
    Interface Delete_RefCodes
       Module Procedure Delete_RefCodes_FAtom
       Module Procedure Delete_RefCodes_Molcrys
       Module Procedure Delete_RefCodes_Molec
    End Interface

    Interface Fill_RefCodes
       Module Procedure Fill_RefCodes_FAtom
       Module Procedure Fill_RefCodes_Molcrys
       Module Procedure Fill_RefCodes_Molec
    End Interface

    Interface Get_ConCodes_Line
       Module Procedure Get_ConCodes_Line_FAtom
       Module Procedure Get_ConCodes_Line_Molcrys
       Module Procedure Get_ConCodes_Line_Molec
    End Interface

    Interface Get_RefCodes_Line
       Module Procedure Get_RefCodes_Line_FAtom
       Module Procedure Get_RefCodes_Line_Molcrys
       Module Procedure Get_RefCodes_Line_Molec
    End Interface

    Interface Init_RefCodes
       Module Procedure Init_RefCodes_FAtom
       Module Procedure Init_RefCodes_Molcrys
       Module Procedure Init_RefCodes_Molec
    End Interface

    Interface Read_RefCodes_File
       Module Procedure Read_RefCodes_File_FAtom
       Module Procedure Read_RefCodes_File_Molcrys
       Module Procedure Read_RefCodes_File_Molec
    End Interface

    Interface VState_to_AtomsPar
       Module Procedure VState_to_AtomsPar_FAtom
       Module Procedure VState_to_AtomsPar_Molcrys
       Module Procedure VState_to_AtomsPar_Molec
    End Interface

    Interface Write_Info_RefCodes
       Module Procedure Write_Info_RefCodes_FAtom
       Module Procedure Write_Info_RefCodes_Molcrys
       Module Procedure Write_Info_RefCodes_Molec
    End Interface

 Contains
    !!----
    !!---- Subroutine Allocate_VParam(N)
    !!----    integer, intent(in) :: N
    !!----
    !!---- Allocate vectors Vec, Bounds, Boundary_Conditions
    !!---- If N is equal zero the deallocating vectors
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_VParam(N)
       !---- Arguments ----!
       integer, intent(in) :: N

       if (allocated(V_Vec))    deallocate(V_Vec)
       if (allocated(V_Name))   deallocate(V_Name)
       if (allocated(V_Bounds)) deallocate(V_Bounds)
       if (allocated(V_BCon))   deallocate(V_BCon)
       if (allocated(V_Shift))  deallocate(V_Shift)
       if (allocated(V_List))   deallocate(V_List)

       if (N > 0) then
          allocate(V_Vec(n))
          V_Vec=0.0
          allocate(V_Name(n))
          V_Name=" "
          allocate(V_Bounds(3,n))
          V_Bounds=0.0
          allocate(V_BCon(n))
          V_BCon=0
          allocate(V_Shift(n))
          V_Shift=0.0
          allocate(V_List(n))
          V_list=0
          np_max=n
       else
          np_max=0
       end if

       return
    End Subroutine Allocate_VParam

    !!----
    !!---- Subroutine Allocate_RestParam(file_dat)
    !!----    Type(file_list_type),     intent( in)    :: file_dat
    !!----
    !!---- Allocate vectors Ang_Rest, Dist_Rest, Tor_Rest
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_RestParam(file_dat)
       !---- Arguments ----!
       Type(file_list_type),     intent( in)    :: file_dat

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=15),dimension(40) :: car
       integer                         :: i,nc,nr

       if (allocated(Ang_Rest)) deallocate(Ang_Rest)
       if (allocated(Dis_Rest)) deallocate(Dis_Rest)
       if (allocated(Tor_Rest)) deallocate(Tor_Rest)

       NP_Rest_Ang=0
       NP_Rest_Dis=0
       NP_Rest_Tor=0

       !---- Dimension for AFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "AFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/3
       end do
       if (nr >0) then
          allocate(Ang_Rest(nr))
          ang_rest%aobs =0.0
          ang_rest%acalc=0.0
          ang_rest%sigma=0.0
          ang_rest%p(1) = 0
          ang_rest%p(2) = 0
          ang_rest%STCode(1)=" "
          ang_rest%STCode(2)=" "
       end if

       !---- Dimension for DFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "DFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/2
          if (modulo(nc,2) == 0) nr=nr-1
       end do
       if (nr >0) then
          allocate(Dis_Rest(nr))
          dis_rest%dobs =0.0
          dis_rest%dcalc=0.0
          dis_rest%sigma=0.0
          dis_rest%p(1) = 0
          dis_rest%p(2) = 0
          dis_rest%STCode=" "
       end if

       !---- Dimension for TFIX ----!
       nr=0
       do i=1,file_dat%nlines
          line=adjustl(file_dat%line(i))
          if (u_case(line(1:4)) /= "TFIX") cycle
          call cutst(line)
          call getword(line,car,nc)
          nr=nr+nc/4
       end do
       if (nr >0) then
          allocate(Tor_Rest(nr))
          tor_rest%tobs=0.0
          tor_rest%tcalc=0.0
          tor_rest%sigma=0.0
          tor_rest%p(1) = 0
          tor_rest%p(2) = 0
          tor_rest%STCode(1)=" "
          tor_rest%STCode(2)=" "
          tor_rest%STCode(3)=" "
       end if

       return
    End Subroutine Allocate_RestParam

    !!--++
    !!--++ Subroutine Delete_RefCodes_FAtom(N, FAtom)
    !!--++    integer,              intent(in)     :: N
    !!--++    type(Atom_List_Type), intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_FAtom(N, FAtom)
       !---- Arguments ----!
       integer,              intent(in)     :: N
       type(Atom_List_Type), intent(in out) :: FAtom

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j

       deleted=.false.

       !---- Eliminate N Parameter ----!
       do i=1,FAtom%natoms
          do j=1,3
             if (FAtom%atom(i)%lx(j) == N) then
                 FAtom%atom(i)%lx(j)=0
                 FAtom%atom(i)%mx(j)=0.0
                 deleted=.true.
             end if
          end do

          if (FAtom%atom(i)%lbiso == N) then
              FAtom%atom(i)%lbiso=0
              FAtom%atom(i)%mbiso=0.0
              deleted=.true.
          end if

          if (FAtom%atom(i)%locc == N) then
              FAtom%atom(i)%locc=0
              FAtom%atom(i)%mocc=0.0
              deleted=.true.
          end if

          do j=1,6
             if (FAtom%atom(i)%lu(j) == N) then
                 FAtom%atom(i)%lu(j)=0
                 FAtom%atom(i)%mu(j)=0.0
                 deleted=.true.
             end if
          end do
       end do

       !---- Updating Variables ----!
       do i=1,FAtom%natoms
          do j=1,3
             if (FAtom%atom(i)%lx(j) > N) then
                FAtom%atom(i)%lx(j)=FAtom%atom(i)%lx(j)-1
             end if
          end do

          if (FAtom%atom(i)%lbiso > N) then
             FAtom%atom(i)%lbiso=FAtom%atom(i)%lbiso-1
          end if

          if (FAtom%atom(i)%locc > N) then
             FAtom%atom(i)%locc=FAtom%atom(i)%locc-1
          end if

          do j=1,6
             if (FAtom%atom(i)%lu(j) > N) then
                FAtom%atom(i)%lu(j)=FAtom%atom(i)%lu(j)-1
             end if
          end do
       end do

       !---- Updating V_Vectors ----!
       if (deleted) then
          do i=N+1,Np_refi
             V_Vec(i-1)=V_Vec(i)
             V_Name(i-1)=V_Name(i)
             V_Bounds(:,i-1)=V_Bounds(:,i)
             V_BCon(i-1)=V_BCon(i)
             V_List(i-1)=V_List(i)
          end do
          V_Vec(np_refi)=0.0
          V_Name(np_refi)=" "
          V_Bounds(:,np_refi)=0.0
          V_BCon(np_refi)=0
          V_List(np_refi)=0

          np_refi=np_refi-1
       end if

       return
    End Subroutine Delete_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Delete_RefCodes_MolCrys(N,MolCrys)
    !!--++    integer,                      intent(in)     :: N
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_MolCrys(N,Molcrys)
       !---- Arguments ----!
       integer,                      intent(in)     :: N
       type(molecular_Crystal_type), intent(in out) :: MolCrys

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j,k

       deleted=.false.

       if (MolCrys%N_Free > 0 ) then
          do i=1,MolCrys%N_Free
             do j=1,3
                if (MolCrys%Atm(i)%lx(j) == N) then
                    MolCrys%Atm(i)%lx(j)=0
                    MolCrys%Atm(i)%mx(j)=0.0
                    deleted=.true.
                end if
             end do

             if (MolCrys%Atm(i)%lbiso == N) then
                 MolCrys%Atm(i)%lbiso=0
                 MolCrys%Atm(i)%mbiso=0.0
                 deleted=.true.
             end if

             if (MolCrys%Atm(i)%locc == N) then
                 MolCrys%Atm(i)%locc=0
                 MolCrys%Atm(i)%mocc=0.0
                 deleted=.true.
             end if

             do j=1,6
                if (MolCrys%Atm(i)%lu(j) == N) then
                    MolCrys%Atm(i)%lu(j)=0
                    MolCrys%Atm(i)%mu(j)=0.0
                    deleted=.true.
                end if
             end do
          end do

          do i=1,MolCrys%N_Free
             do j=1,3
                if (MolCrys%Atm(i)%lx(j) > N) then
                    MolCrys%Atm(i)%lx(j)=MolCrys%Atm(i)%lx(j)-1
                end if
             end do

             if (MolCrys%Atm(i)%lbiso > N) then
                MolCrys%Atm(i)%lbiso=MolCrys%Atm(i)%lbiso-1
             end if

             if (MolCrys%Atm(i)%locc > N) then
                MolCrys%Atm(i)%locc=MolCrys%Atm(i)%locc-1
             end if

             do j=1,6
                if (MolCrys%Atm(i)%lu(j) > N) then
                   MolCrys%Atm(i)%lu(j)=MolCrys%Atm(i)%lu(j)-1
                end if
             end do
          end do
       end if

       if (MolCrys%N_Mol > 0 ) then

          do k=1,MolCrys%N_Mol
             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) == N) then
                   Molcrys%Mol(k)%lxcentre(j)=0
                   Molcrys%Mol(k)%mxcentre(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) == N) then
                   Molcrys%Mol(k)%lOrient(j)=0
                   Molcrys%Mol(k)%mOrient(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) == N) then
                   Molcrys%Mol(k)%lT_TLS(j)=0
                   Molcrys%Mol(k)%mT_TLS(j)=0.0
                   deleted=.true.
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) == N) then
                   Molcrys%Mol(k)%lL_TLS(j)=0
                   Molcrys%Mol(k)%mL_TLS(j)=0.0
                   deleted=.true.
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) == N) then
                      Molcrys%Mol(k)%lS_TLS(i,j)=0
                      Molcrys%Mol(k)%mS_TLS(i,j)=0.0
                      deleted=.true.
                   end if
                end do
             end do

             !---- Updating ----!
             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) > N) then
                   Molcrys%Mol(k)%lxcentre(j)=Molcrys%Mol(k)%lxcentre(j)-1
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) > N) then
                   Molcrys%Mol(k)%lOrient(j)=Molcrys%Mol(k)%lOrient(j)-1
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) > N) then
                   Molcrys%Mol(k)%lT_TLS(j)=Molcrys%Mol(k)%lT_TLS(j)-1
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) > N) then
                   Molcrys%Mol(k)%lL_TLS(j)=Molcrys%Mol(k)%lL_TLS(j)-1
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) > N) then
                      Molcrys%Mol(k)%lS_TLS(i,j)=Molcrys%Mol(k)%lS_TLS(i,j)-1
                   end if
                end do
             end do

             if (Molcrys%Mol(k)%natoms <=0) cycle

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (MolCrys%Mol(k)%lI_coor(j,i) == N) then
                      MolCrys%Mol(k)%lI_coor(j,i)=0
                      MolCrys%Mol(k)%mI_coor(j,i)=0.0
                      deleted=.true.
                   end if
                end do

                if (MolCrys%Mol(k)%lbiso(i) == N) then
                   MolCrys%Mol(k)%lbiso(i)=0
                   MolCrys%Mol(k)%mbiso(i)=0.0
                   deleted=.true.
                end if

                if (MolCrys%Mol(k)%locc(i) == N) then
                   MolCrys%Mol(k)%locc(i)=0
                   MolCrys%Mol(k)%mocc(i)=0.0
                   deleted=.true.
                end if
             end do

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (MolCrys%Mol(k)%lI_coor(j,i) > N) then
                      MolCrys%Mol(k)%lI_coor(j,i)=MolCrys%Mol(k)%lI_coor(j,i)-1
                   end if
                end do

                if (MolCrys%Mol(k)%lbiso(i) > N) then
                   MolCrys%Mol(k)%lbiso(i)=MolCrys%Mol(k)%lbiso(i)-1
                end if

                if (MolCrys%Mol(k)%locc(i) > N) then
                   MolCrys%Mol(k)%locc(i)=MolCrys%Mol(k)%locc(i)-1
                end if
             end do

          end do
       end if

       !---- Updating V_Vectors ----!
       if (deleted) then
          do i=N+1,Np_refi
             V_Vec(i-1)=V_Vec(i)
             V_Name(i-1)=V_Name(i)
             V_Bounds(:,i-1)=V_Bounds(:,i)
             V_BCon(i-1)=V_BCon(i)
             V_List(i-1)=V_List(i)
          end do
          V_Vec(np_refi)=0.0
          V_Name(np_refi)=" "
          V_Bounds(:,np_refi)=0.0
          V_BCon(np_refi)=0
          V_List(np_refi)=0

          np_refi=np_refi-1
       end if

       return
    End Subroutine Delete_RefCodes_MolCrys

    !!--++
    !!--++ Subroutine Delete_RefCodes_Molec(N,Molec)
    !!--++    integer,             intent(in)     :: N
    !!--++    type(molecule_type), intent(in out) :: Molec
    !!--++
    !!--++ Overloaded
    !!--++ Delete the number of Refinable Parameter (N) on the list
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Delete_RefCodes_Molec(N,Molec)
       !---- Arguments ----!
       integer,             intent(in)     :: N
       type(molecule_type), intent(in out) :: Molec

       !---- Local Variables ----!
       logical :: deleted
       integer :: i,j

       deleted=.false.

       do j=1,3
          if (Molec%lxcentre(j) == N) then
             Molec%lxcentre(j)=0
             Molec%mxcentre(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,3
          if (Molec%lOrient(j) == N) then
             Molec%lOrient(j)=0
             Molec%mOrient(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,6
          if (Molec%lT_TLS(j) == N) then
             Molec%lT_TLS(j)=0
             Molec%mT_TLS(j)=0.0
             deleted=.true.
          end if
       end do

       do j=1,6
          if (Molec%lL_TLS(j) == N) then
             Molec%lL_TLS(j)=0
             Molec%mL_TLS(j)=0.0
             deleted=.true.
          end if
       end do

       do i=1,3
          do j=1,3
             if (Molec%lS_TLS(i,j) == N) then
                Molec%lS_TLS(i,j)=0
                Molec%mS_TLS(i,j)=0.0
                deleted=.true.
             end if
          end do
       end do

       !---- Updating ----!
       do j=1,3
          if (Molec%lxcentre(j) > N) then
             Molec%lxcentre(j)=Molec%lxcentre(j)-1
          end if
       end do

       do j=1,3
          if (Molec%lOrient(j) > N) then
             Molec%lOrient(j)=Molec%lOrient(j)-1
          end if
       end do

       do j=1,6
          if (Molec%lT_TLS(j) > N) then
             Molec%lT_TLS(j)=Molec%lT_TLS(j)-1
          end if
       end do

       do j=1,6
          if (Molec%lL_TLS(j) > N) then
             Molec%lL_TLS(j)=Molec%lL_TLS(j)-1
          end if
       end do

       do i=1,3
          do j=1,3
             if (Molec%lS_TLS(i,j) > N) then
                Molec%lS_TLS(i,j)=Molec%lS_TLS(i,j)-1
             end if
          end do
       end do

       if (molec%natoms <=0) return

       do i=1,Molec%Natoms
          do j=1,3
             if (Molec%lI_coor(j,i) == N) then
                Molec%lI_coor(j,i)=0
                Molec%mI_coor(j,i)=0.0
                deleted=.true.
             end if
          end do

          if (Molec%lbiso(i) == N) then
             Molec%lbiso(i)=0
             Molec%mbiso(i)=0.0
             deleted=.true.
          end if

          if (Molec%locc(i) == N) then
             Molec%locc(i)=0
             Molec%mocc(i)=0.0
             deleted=.true.
          end if
       end do

       !---- Updating ----!
       do i=1,Molec%Natoms
          do j=1,3
             if (Molec%lI_coor(j,i) > N) then
                Molec%lI_coor(j,i)=Molec%lI_coor(j,i)-1
             end if
          end do

          if (Molec%lbiso(i) > N) then
             Molec%lbiso(i)=Molec%lbiso(i)-1
          end if

          if (Molec%locc(i) > N) then
             Molec%locc(i)=Molec%locc(i)-1
          end if
       end do

       !---- Updating V_Vectors ----!
       if (deleted) then
          do i=N+1,Np_refi
             V_Vec(i-1)=V_Vec(i)
             V_Name(i-1)=V_Name(i)
             V_Bounds(:,i-1)=V_Bounds(:,i)
             V_BCon(i-1)=V_BCon(i)
             V_List(i-1)=V_List(i)
          end do
          V_Vec(np_refi)=0.0
          V_Name(np_refi)=" "
          V_Bounds(:,np_refi)=0.0
          V_BCon(np_refi)=0
          V_List(np_refi)=0

          np_refi=np_refi-1
       end if

       return
    End Subroutine Delete_RefCodes_Molec

    !!--++
    !!--++ Subroutine Fill_RefCodes_FAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FAtom,Spg)
    !!--++    integer,                       intent(in)     :: Key
    !!--++    character(len=*),              intent(in)     :: Dire
    !!--++    integer,                       intent(in)     :: Na
    !!--++    integer,                       intent(in)     :: Nb
    !!--++    real,                          intent(in)     :: Xl
    !!--++    real,                          intent(in)     :: Xu
    !!--++    real,                          intent(in)     :: Xs
    !!--++    integer,                       intent(in)     :: Ic
    !!--++    type(Atom_List_Type),          intent(in out) :: FAtom
    !!--++    type(space_group_type),        intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Free Atoms
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_FAtom(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,FAtom,Spg)
       !---- Arguments ----!
       integer,                       intent(in)     :: Key
       character(len=*),              intent(in)     :: Dire
       integer,                       intent(in)     :: Na
       integer,                       intent(in)     :: Nb
       real,                          intent(in)     :: Xl
       real,                          intent(in)     :: Xu
       real,                          intent(in)     :: Xs
       integer,                       intent(in)     :: Ic
       type(Atom_List_Type),          intent(in out) :: FAtom
       type(space_group_type),        intent(in)     :: Spg

       !---- Local variables ----!
       integer           :: j,k,nc
       character(len=20) :: car

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)
                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (FAtom%atom(na)%lx(nb) /=0) then
                            nc=FAtom%atom(na)%lx(nb)
                            call Delete_RefCodes(nc,FAtom)
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (FAtom%atom(na)%lbiso /=0) then
                            nc=FAtom%atom(na)%lbiso
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (FAtom%atom(na)%locc /=0) then
                            nc=FAtom%atom(na)%locc
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case ( 6:11)
                         !---- B11_,...,B23_ ----!
                         if (FAtom%atom(na)%lu(nb-5) /=0) then
                            nc=FAtom%atom(na)%lu(nb-5)
                            call Delete_RefCodes(nc,fatom)
                         end if

                      case (12)
                         !---- Banis_ ----!
                         do j=1,6
                            if (FAtom%atom(na)%lu(j) /=0) then
                               nc=FAtom%atom(na)%lu(j)
                               call Delete_RefCodes(nc,fatom)
                            end if
                         end do

                      case (13:)
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined for this type of variable "
                         return
                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) /=0) then
                         nc=FAtom%atom(na)%lx(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (FAtom%atom(na)%locc /=0) then
                      nc=FAtom%atom(na)%locc
                      call Delete_RefCodes(nc,fatom)
                   end if

                case (3)
                   !---- BIS ----!
                   if (FAtom%atom(na)%lbiso /=0) then
                      nc=FAtom%atom(na)%lbiso
                      call Delete_RefCodes(nc,fatom)
                   end if

               case (4)
                  !---- BAN ----!
                  do j=1,6
                     if (FAtom%atom(na)%lu(j) /=0) then
                        nc=FAtom%atom(na)%lu(j)
                        call Delete_RefCodes(nc,fatom)
                     end if
                  end do

                case (5)
                   !---- ALL ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) /=0) then
                         nc=FAtom%atom(na)%lx(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do
                   if (FAtom%atom(na)%locc /=0) then
                      nc=FAtom%atom(na)%locc
                      call Delete_RefCodes(nc,fatom)
                   end if
                   if (FAtom%atom(na)%lbiso /=0) then
                      nc=FAtom%atom(na)%lbiso
                      call Delete_RefCodes(nc,fatom)
                   end if
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) /=0) then
                         nc=FAtom%atom(na)%lu(j)
                         call Delete_RefCodes(nc,fatom)
                      end if
                   end do

                case (6:)
                   err_refcodes=.true.
                   err_mess_refcodes="Incompatible information for this type of variable "
                   return
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)

                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (FAtom%atom(na)%lx(nb) ==0) then
                            FAtom%atom(na)%mx(nb)=1.0
                            call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                                 FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                            if (FAtom%atom(na)%lx(nb) == np_refi) then
                               V_Vec(np_refi)=FAtom%atom(na)%x(nb)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (FAtom%atom(na)%lbiso ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=FAtom%atom(na)%biso
                            V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                            FAtom%atom(na)%mbiso=1.0
                            FAtom%atom(na)%lbiso=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (FAtom%atom(na)%locc ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=FAtom%atom(na)%occ
                            V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                            FAtom%atom(na)%mocc=1.0
                            FAtom%atom(na)%locc=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 6:11)
                         !---- B11_,...,B23_ ----!
                         if (FAtom%atom(na)%lu(nb-5) ==0) then
                            FAtom%atom(na)%mu(nb-5)=1.0
                            call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                                 np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                            if (FAtom%atom(na)%lu(nb-5) == np_refi) then
                               V_Vec(np_refi)=FAtom%atom(na)%u(nb-5)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(FAtom%atom(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case (12)
                         !---- Banis_ ----!
                         do j=1,6
                            if (FAtom%atom(na)%lu(j) ==0) then
                               FAtom%atom(na)%mu(j)=1.0
                               call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                                    np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                               if (FAtom%atom(na)%lu(j) == np_refi) then
                                  V_Vec(np_refi)=FAtom%atom(na)%u(j)
                                  V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         end do

                      case (13:)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined by this type of variables"
                         return

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) ==0) then
                         FAtom%atom(na)%mx(j)=1.0
                         call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                              FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                         if (FAtom%atom(na)%lx(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%x(j)
                            V_Name(np_refi)=trim(code_nam(j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (FAtom%atom(na)%locc ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%occ
                      V_Name(np_refi)=trim(code_nam(5))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mocc=1.0
                      FAtom%atom(na)%locc=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (3)
                   !---- BIS ----!
                   if (FAtom%atom(na)%lbiso ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%biso
                      V_Name(np_refi)=trim(code_nam(4))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mbiso=1.0
                      FAtom%atom(na)%lbiso=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (4)
                   !---- BAN ----!
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) ==0) then
                         FAtom%atom(na)%mu(j)=1.0
                         call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                              np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                         if (FAtom%atom(na)%lu(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%u(j)
                            V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (5)
                   !---- ALL ----!
                   do j=1,3
                      if (FAtom%atom(na)%lx(j) ==0) then
                         FAtom%atom(na)%mx(j)=1.0
                         call get_atompos_ctr(FAtom%atom(na)%x, Spg, np_refi,   &
                                              FAtom%atom(na)%lx, FAtom%atom(na)%mx)
                         if (FAtom%atom(na)%lx(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%x(j)
                            V_Name(np_refi)=trim(code_nam(j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do
                   if (FAtom%atom(na)%locc ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%occ
                      V_Name(np_refi)=trim(code_nam(5))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mocc=1.0
                      FAtom%atom(na)%locc=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   if (FAtom%atom(na)%lbiso ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=FAtom%atom(na)%biso
                      V_Name(np_refi)=trim(code_nam(4))//trim(FAtom%atom(na)%lab)
                      FAtom%atom(na)%mbiso=1.0
                      FAtom%atom(na)%lbiso=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   do j=1,6
                      if (FAtom%atom(na)%lu(j) ==0) then
                         FAtom%atom(na)%mu(j)=1.0
                         call get_atombet_ctr(FAtom%atom(na)%x,FAtom%atom(na)%u,Spg, &
                                              np_refi,FAtom%atom(na)%lu,FAtom%atom(na)%mu)
                         if (FAtom%atom(na)%lu(j) == np_refi) then
                            V_Vec(np_refi)=FAtom%atom(na)%u(j)
                            V_Name(np_refi)=trim(code_nam(5+j))//trim(FAtom%atom(na)%lab)
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case(6:)
                   err_refcodes=.true.
                   err_mess_refcodes="Option Not defined by this type of variables"
                   return
             end select
       end select

       return
    End Subroutine Fill_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Fill_RefCodes_Molcrys(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molcrys,NMol)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    integer,                      intent(in)     :: Na
    !!--++    integer,                      intent(in)     :: Nb
    !!--++    real,                         intent(in)     :: Xl
    !!--++    real,                         intent(in)     :: Xu
    !!--++    real,                         intent(in)     :: Xs
    !!--++    integer,                      intent(in)     :: Ic
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    integer,                      intent(in)     :: NMol
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Free Atoms
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_Molcrys(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molcrys,Nmol)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       integer,                      intent(in)     :: Na
       integer,                      intent(in)     :: Nb
       real,                         intent(in)     :: Xl
       real,                         intent(in)     :: Xu
       real,                         intent(in)     :: Xs
       integer,                      intent(in)     :: Ic
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       integer,                      intent(in)     :: NMol

       !---- Local variables ----!
       character(len=2) :: car
       integer          :: i, j, k, nc, naa

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)
                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lx(nb) /=0) then
                               nc=molcrys%atm(na)%lx(nb)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lI_coor(nb,naa) /=0) then
                                  nc=molcrys%mol(i)%lI_coor(nb,naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lbiso /=0) then
                               nc=molcrys%atm(na)%lbiso
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lbiso(naa) /=0) then
                                  nc=molcrys%mol(i)%lbiso(naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%locc /=0) then
                               nc=molcrys%atm(na)%locc
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%locc(naa) /=0) then
                                  nc=molcrys%mol(i)%locc(naa)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         end if

                      case ( 6:11)
                         !---- B11_, ..., B23_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lu(nb-5) /=0) then
                               nc=molcrys%atm(na)%lu(nb-5)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         else
                            err_refcodes=.true.
                            err_mess_refcodes="Option Not defined"
                            return
                         end if

                      case (12)
                         !---- Banis_ ----!
                         if (na <= molcrys%n_free) then
                            do j=1,6
                               if (molcrys%atm(na)%lu(j) /=0) then
                                  nc=molcrys%atm(na)%lu(j)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                            end do
                         else
                            err_refcodes=.true.
                            err_mess_refcodes="Option Not defined"
                            return
                         end if

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               err_mess_refcodes="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  if (molcrys%mol(i)%lxcentre(nb-12) /=0) then
                                     nc=molcrys%mol(i)%lxcentre(nb-12)
                                     call Delete_RefCodes(nc,MolCrys)
                                  end if
                               end do

                            case (1:)
                               if (molcrys%mol(nmol)%lxcentre(nb-12) /=0) then
                                  nc=molcrys%mol(nmol)%lxcentre(nb-12)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                         end select

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               err_mess_refcodes="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  if (molcrys%mol(i)%lOrient(nb-15) /=0) then
                                     nc=molcrys%mol(i)%lOrient(nb-15)
                                     call Delete_RefCodes(nc,MolCrys)
                                  end if
                               end do

                            case (1:)
                               if (molcrys%mol(nmol)%lOrient(nb-15) /=0) then
                                  nc=molcrys%mol(nmol)%lOrient(nb-15)
                                  call Delete_RefCodes(nc,MolCrys)
                               end if
                         end select

                      case (19:21)
                         !!! Not yet implement !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) /=0) then
                            nc=molcrys%atm(na)%lx(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) /=0) then
                               nc=molcrys%mol(i)%lI_coor(j,naa)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         end do
                      end do
                   end if

                case (2)
                   !---- OCC ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc /=0) then
                         nc=molcrys%atm(na)%locc
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) /=0) then
                            nc=molcrys%mol(i)%locc(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                case (3)
                   !---- BIS ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso /=0) then
                         nc=molcrys%atm(na)%lbiso
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) /=0) then
                            nc=molcrys%mol(i)%lbiso(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                case (4)
                   !---- BAN ----!
                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) /=0) then
                            nc=molcrys%atm(na)%lu(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      err_refcodes=.true.
                      err_mess_refcodes="Option Not defined"
                      return
                   end if

                case (5)
                   !---- ALL ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) /=0) then
                            nc=molcrys%atm(na)%lx(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) /=0) then
                               nc=molcrys%mol(i)%lI_coor(j,naa)
                               call Delete_RefCodes(nc,MolCrys)
                            end if
                         end do
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc /=0) then
                         nc=molcrys%atm(na)%locc
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) /=0) then
                            nc=molcrys%mol(i)%locc(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso /=0) then
                         nc=molcrys%atm(na)%lbiso
                         call Delete_RefCodes(nc,MolCrys)
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) /=0) then
                            nc=molcrys%mol(i)%lbiso(naa)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) /=0) then
                            nc=molcrys%atm(na)%lu(j)
                            call Delete_RefCodes(nc,MolCrys)
                         end if
                      end do
                   end if

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) /=0) then
                                  nc=molcrys%mol(i)%lxcentre(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do

                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) /=0) then
                                  nc=molcrys%mol(i)%lOrient(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) /=0) then
                               nc=molcrys%mol(nmol)%lxcentre(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do

                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) /=0) then
                               nc=molcrys%mol(nmol)%lOrient(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                   !!! THE Not Implemented !!!

                case (6)
                   !---- CEN ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) /=0) then
                                  nc=molcrys%mol(i)%lxcentre(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) /=0) then
                               nc=molcrys%mol(nmol)%lxcentre(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                case (7)
                   !---- ORI ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) /=0) then
                                  nc=molcrys%mol(i)%lOrient(j)
                                  call Delete_RefCodes(nc,molcrys)
                               end if
                            end do
                         end do

                      case (1:)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) /=0) then
                               nc=molcrys%mol(nmol)%lOrient(j)
                               call Delete_RefCodes(nc,molcrys)
                            end if
                         end do
                   end select

                case (8)
                   !---- THE ----!
                   !!! Not yet Implemented !!!!
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)
                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lx(nb) ==0) then
                               molcrys%atm(na)%mx(nb)=1.0
                               call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                    molcrys%Spg,np_refi, &
                                                    molcrys%atm(na)%lx,  &
                                                    molcrys%atm(na)%mx)
                               if (molcrys%atm(na)%lx(nb) == np_refi) then
                                  V_Vec(np_refi)=molcrys%atm(na)%x(nb)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if

                               if (molcrys%mol(i)%lI_coor(nb,naa) ==0) then
                                  molcrys%mol(i)%mI_coor(nb,naa)=1.0
                                  call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa),  &
                                                       molcrys%Spg, np_refi,          &
                                                       molcrys%mol(i)%lI_coor(:,naa), &
                                                       molcrys%mol(i)%mI_coor(:,naa))
                                  if (molcrys%mol(i)%lI_coor(nb,naa) == np_refi) then
                                     V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                     V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=na
                                  else
                                     np_refi=np_refi-1
                                  end if
                               end if
                            end do
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lbiso ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%atm(na)%biso
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                               molcrys%atm(na)%mbiso=1.0
                               molcrys%atm(na)%lbiso=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%lbiso(naa) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                  molcrys%mol(i)%mbiso(naa)=1.0
                                  molcrys%mol(i)%lbiso(naa)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               end if
                            end do
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%locc ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%atm(na)%occ
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                               molcrys%atm(na)%mocc=1.0
                               molcrys%atm(na)%locc=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            end if
                         else
                            naa=na-molcrys%n_free
                            do i=1,molcrys%n_mol
                               if (naa > molcrys%mol(i)%natoms) then
                                  naa=naa-molcrys%mol(i)%natoms
                                  cycle
                               end if
                               if (molcrys%mol(i)%locc(naa) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%occ(naa)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%mol(i)%AtName(naa))
                                  molcrys%mol(i)%mocc(naa)=1.0
                                  molcrys%mol(i)%locc(naa)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               end if
                            end do
                         end if

                      case ( 6:11)
                         !---- B11_, ..., B23_ ----!
                         if (na <= molcrys%n_free) then
                            if (molcrys%atm(na)%lu(nb-5) ==0) then
                               molcrys%atm(na)%mu(nb-5)=1.0
                               call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                   np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                               if (molcrys%atm(na)%lu(nb-5) == np_refi) then
                                  V_Vec(np_refi)=molcrys%atm(na)%u(nb-5)
                                  V_Name(np_refi)=trim(code_nam(nb))//trim(molcrys%atm(na)%lab)
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         else
                            err_refcodes=.true.
                            err_mess_refcodes="Option Not defined"
                            return
                         end if

                      case (12)
                         !---- Banis_ ----!
                         if (na <= molcrys%n_free) then
                            do j=1,6
                               if (molcrys%atm(na)%lu(j) ==0) then
                                  molcrys%atm(na)%mu(j)=1.0
                                  call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                       np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                                  if (molcrys%atm(na)%lu(j) == np_refi) then
                                     V_Vec(np_refi)=molcrys%atm(na)%u(j)
                                     V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=na
                                  else
                                     np_refi=np_refi-1
                                  end if
                               end if
                            end do
                         else
                            err_refcodes=.true.
                            err_mess_refcodes="Option Not defined"
                            return
                         end if

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               err_mess_refcodes="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  write(unit=car,fmt="(i2)") i
                                  car=adjustl(car)
                                  if (molcrys%mol(i)%lxcentre(nb-12) ==0) then
                                     np_refi=np_refi+1
                                     V_Vec(np_refi)=molcrys%mol(i)%xcentre(nb-12)
                                     V_Name(np_refi)=trim(code_nam(nb))//"Mol"//trim(car)
                                     molcrys%mol(i)%mxcentre(nb-12)=1.0
                                     molcrys%mol(i)%lxcentre(nb-12)=np_refi
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=-i
                                  end if
                               end do

                            case (1:)
                               write(unit=car,fmt="(i2)") nmol
                               car=adjustl(car)
                               if (molcrys%mol(nmol)%lxcentre(nb-12) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(nb-12)
                                  V_Name(np_refi)=trim(code_nam(nb))//"entre_Mol"//trim(car)
                                  molcrys%mol(nmol)%mxcentre(nb-12)=1.0
                                  molcrys%mol(nmol)%lxcentre(nb-12)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-nmol
                               end if
                         end select

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         select case (nmol)
                            case (-1)
                               err_refcodes=.true.
                               err_mess_refcodes="Option Not defined"
                               return

                            case (0)
                               do i=1,molcrys%n_mol
                                  write(unit=car,fmt="(i2)") i
                                  car=adjustl(car)
                                  if (molcrys%mol(i)%lOrient(nb-15) ==0) then
                                     np_refi=np_refi+1
                                     V_Vec(np_refi)=molcrys%mol(i)%Orient(nb-15)
                                     V_Name(np_refi)=trim(code_nam(nb))//"Orient_Mol"//trim(car)
                                     molcrys%mol(i)%mOrient(nb-15)=1.0
                                     molcrys%mol(i)%lOrient(nb-15)=np_refi
                                     V_Bounds(1,np_refi)=xl
                                     V_Bounds(2,np_refi)=xu
                                     V_Bounds(3,np_refi)=xs
                                     V_BCon(np_refi)=ic
                                     V_list(np_refi)=-i
                                  end if
                               end do

                            case (1:)
                               write(unit=car,fmt="(i2)") nmol
                               car=adjustl(car)
                               if (molcrys%mol(nmol)%lOrient(nb-15) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(nmol)%Orient(nb-15)
                                  V_Name(np_refi)=trim(code_nam(nb))//"Orient_Mol"//trim(car)
                                  molcrys%mol(nmol)%mOrient(nb-15)=1.0
                                  molcrys%mol(nmol)%lOrient(nb-15)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-nmol
                               end if
                         end select

                      case (19:21)
                         !!! Not Yet Implemented !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) ==0) then
                            molcrys%atm(na)%mx(j)=1.0
                            call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                 molcrys%Spg,np_refi, &
                                                 molcrys%atm(na)%lx,  &
                                                 molcrys%atm(na)%mx)
                            if (molcrys%atm(na)%lx(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%x(j)
                               V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                             naa=naa-molcrys%mol(i)%natoms
                             cycle
                         end if
                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) ==0) then
                               molcrys%mol(i)%mI_coor(nb,naa)=1.0
                               call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa), &
                                                    molcrys%Spg,np_refi,   &
                                                    molcrys%mol(i)%lI_coor(:,naa), &
                                                    molcrys%mol(i)%mI_coor(:,naa))
                               if (molcrys%mol(i)%lI_coor(j,naa) == np_refi) then
                                  V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                  V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%mol(i)%AtName(naa))
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         end do
                      end do
                   end if

                case (2)
                   !---- OCC ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%occ
                         V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mocc=1.0
                         molcrys%atm(na)%locc=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if

                         if (molcrys%mol(i)%locc(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%Occ(naa)
                            V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mocc(naa)=1.0
                            molcrys%mol(i)%locc(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                case (3)
                   !---- BIS ----!
                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%biso
                         V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mbiso=1.0
                         molcrys%atm(na)%lbiso=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%lbiso(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                            V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mbiso(naa)=1.0
                            molcrys%mol(i)%lbiso(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                case (4)
                   !---- BAN ----!
                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) ==0) then
                            molcrys%atm(na)%mu(j)=1.0
                            call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                 np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                            if (molcrys%atm(na)%lu(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%u(j)
                               V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      err_refcodes=.true.
                      err_mess_refcodes="Option Not defined"
                      return
                   end if

                case (5)
                   !---- ALL ----!
                   if (na <= molcrys%n_free) then
                      do j=1,3
                         if (molcrys%atm(na)%lx(j) ==0) then
                            molcrys%atm(na)%mx(j)=1.0
                            call get_atompos_ctr(molcrys%atm(na)%x,   &
                                                 molcrys%Spg,np_refi, &
                                                 molcrys%atm(na)%lx,  &
                                                 molcrys%atm(na)%mx)
                            if (molcrys%atm(na)%lx(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%x(j)
                               V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                             naa=naa-molcrys%mol(i)%natoms
                             cycle
                         end if

                         do j=1,3
                            if (molcrys%mol(i)%lI_coor(j,naa) ==0) then
                               molcrys%mol(i)%mI_coor(nb,naa)=1.0
                               call get_atompos_ctr(molcrys%mol(i)%I_Coor(:,naa), &
                                                    molcrys%Spg,np_refi,   &
                                                    molcrys%mol(i)%lI_coor(:,naa), &
                                                    molcrys%mol(i)%mI_coor(:,naa))
                               if (molcrys%mol(i)%lI_coor(j,naa) == np_refi) then
                                  V_Vec(np_refi)=molcrys%mol(i)%I_Coor(nb,naa)
                                  V_Name(np_refi)=trim(code_nam(j))//trim(molcrys%mol(i)%AtName(naa))
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=na
                               else
                                  np_refi=np_refi-1
                               end if
                            end if
                         end do
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%locc ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%occ
                         V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mocc=1.0
                         molcrys%atm(na)%locc=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if
                         if (molcrys%mol(i)%locc(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%Occ(naa)
                            V_Name(np_refi)=trim(code_nam(5))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mocc(naa)=1.0
                            molcrys%mol(i)%locc(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      if (molcrys%atm(na)%lbiso ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molcrys%atm(na)%biso
                         V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%atm(na)%lab)
                         molcrys%atm(na)%mbiso=1.0
                         molcrys%atm(na)%lbiso=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=na
                      end if
                   else
                      naa=na-molcrys%n_free
                      do i=1,molcrys%n_mol
                         if (naa > molcrys%mol(i)%natoms) then
                            naa=naa-molcrys%mol(i)%natoms
                            cycle
                         end if

                         if (molcrys%mol(i)%lbiso(naa) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molcrys%mol(i)%biso(naa)
                            V_Name(np_refi)=trim(code_nam(4))//trim(molcrys%mol(i)%AtName(naa))
                            molcrys%mol(i)%mbiso(naa)=1.0
                            molcrys%mol(i)%lbiso(naa)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if
                      end do
                   end if

                   if (na <= molcrys%n_free) then
                      do j=1,6
                         if (molcrys%atm(na)%lu(j) ==0) then
                            molcrys%atm(na)%mu(j)=1.0
                            call get_atombet_ctr(molcrys%atm(na)%x,molcrys%atm(na)%u,molcrys%Spg, &
                                                 np_refi,molcrys%atm(na)%lu,molcrys%atm(na)%mu)
                            if (molcrys%atm(na)%lu(j) == np_refi) then
                               V_Vec(np_refi)=molcrys%atm(na)%u(j)
                               V_Name(np_refi)=trim(code_nam(5+j))//trim(molcrys%atm(na)%lab)
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if
                      end do
                   end if

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%xcentre(j)
                                  V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                                  molcrys%mol(i)%mxcentre(j)=1.0
                                  molcrys%mol(i)%lxcentre(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(j)
                               V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                               molcrys%mol(nmol)%mxcentre(j)=1.0
                               molcrys%mol(nmol)%lxcentre(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%Orient(j)
                                  V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                                  molcrys%mol(i)%mOrient(j)=1.0
                                  molcrys%mol(i)%lOrient(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%Orient(j)
                               V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                               molcrys%mol(nmol)%mOrient(j)=1.0
                               molcrys%mol(nmol)%lOrient(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (6)
                   !---- CEN ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lxcentre(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%xcentre(j)
                                  V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                                  molcrys%mol(i)%mxcentre(j)=1.0
                                  molcrys%mol(i)%lxcentre(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lxcentre(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%xcentre(j)
                               V_Name(np_refi)=trim(code_nam(12+j))//"entre_Mol"//trim(car)
                               molcrys%mol(nmol)%mxcentre(j)=1.0
                               molcrys%mol(nmol)%lxcentre(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (7)
                   !---- ORI  ----!
                   select case (nmol)
                      case (-1)
                         err_refcodes=.true.
                         err_mess_refcodes="Option Not defined"
                         return

                      case (0)
                         do i=1,molcrys%n_mol
                            write(unit=car,fmt="(i2)") i
                            car=adjustl(car)
                            do j=1,3
                               if (molcrys%mol(i)%lOrient(j) ==0) then
                                  np_refi=np_refi+1
                                  V_Vec(np_refi)=molcrys%mol(i)%Orient(j)
                                  V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                                  molcrys%mol(i)%mOrient(j)=1.0
                                  molcrys%mol(i)%lOrient(j)=np_refi
                                  V_Bounds(1,np_refi)=xl
                                  V_Bounds(2,np_refi)=xu
                                  V_Bounds(3,np_refi)=xs
                                  V_BCon(np_refi)=ic
                                  V_list(np_refi)=-i
                               end if
                            end do
                         end do

                      case (1:)
                         write(unit=car,fmt="(i2)") nmol
                         car=adjustl(car)
                         do j=1,3
                            if (molcrys%mol(nmol)%lOrient(j) ==0) then
                               np_refi=np_refi+1
                               V_Vec(np_refi)=molcrys%mol(nmol)%Orient(j)
                               V_Name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"//trim(car)
                               molcrys%mol(nmol)%mOrient(j)=1.0
                               molcrys%mol(nmol)%lOrient(j)=np_refi
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=-nmol
                            end if
                         end do
                   end select

                case (8)
                   !!! Not yet implemented !!!

             end select
       end select

       return
    End Subroutine Fill_RefCodes_Molcrys

    !!--++
    !!--++ Subroutine Fill_RefCodes_Molec(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molec,Spg)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    integer,                      intent(in)     :: Na
    !!--++    integer,                      intent(in)     :: Nb
    !!--++    real,                         intent(in)     :: Xl
    !!--++    real,                         intent(in)     :: Xu
    !!--++    real,                         intent(in)     :: Xs
    !!--++    integer,                      intent(in)     :: Ic
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    type(space_group_type),       intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Write on Vectors the Information for Molecule_Type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Fill_RefCodes_Molec(Key,Dire,Na,Nb,Xl,Xu,Xs,Ic,Molec,Spg)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       integer,                      intent(in)     :: Na
       integer,                      intent(in)     :: Nb
       real,                         intent(in)     :: Xl
       real,                         intent(in)     :: Xu
       real,                         intent(in)     :: Xs
       integer,                      intent(in)     :: Ic
       type(molecule_type),          intent(in out) :: Molec
       type(space_group_type),       intent(in)     :: Spg

       !---- Local variables ----!
       integer :: j, k, nc

       call init_err_refcodes()
       if (Na <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Number of atom no defined"
          return
       end if

       select case (dire)
          !---- FIX Directive ----!
          case ("fix")

             select case (key)
                case (0)

                   !---- nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case ( 1:3)
                         !---- X_, Y_, Z_----!
                         if (molec%lI_coor(nb,na) /=0) then
                            nc=molec%lI_coor(nb,na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (molec%lbiso(na) /=0) then
                            nc=molec%lbiso(na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (molec%locc(na) /=0) then
                            nc=molec%locc(na)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case ( 6:12)
                         !---- B11_, ..., B23_ ----!
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         if (molec%lxcentre(nb-12) /=0) then
                            nc=molec%lxcentre(nb-12)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         if (molec%lOrient(nb-15) /=0) then
                            nc=molec%lOrient(nb-15)
                            call Delete_RefCodes(nc,molec)
                         end if

                      case (19:21)
                         !!! Not yet implement !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) /=0) then
                         nc=molec%lI_coor(j,na)
                         call Delete_RefCodes(nc,molec)
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (molec%locc(na) /=0) then
                      nc=molec%locc(na)
                      call Delete_RefCodes(nc,molec)
                   end if

                case (3)
                   !---- BIS ----!
                   if (molec%lbiso(na) /=0) then
                      nc=molec%lbiso(na)
                      call Delete_RefCodes(nc,molec)
                   end if

                case (4)
                   !---- BAN ----!
                   err_refcodes=.true.
                   err_mess_refcodes="Option not defined"
                   return

                case (5)
                  !---- ALL ----!
                  do j=1,3
                     if (molec%lI_coor(j,na) /=0) then
                        nc=molec%lI_coor(j,na)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  if (molec%locc(na) /=0) then
                     nc=molec%locc(na)
                     call Delete_RefCodes(nc,molec)
                  end if
                  if (molec%lbiso(na) /=0) then
                     nc=molec%lbiso(na)
                     call Delete_RefCodes(nc,molec)
                  end if
                  do j=1,3
                     if (molec%lxcentre(j) /=0) then
                        nc=molec%lxcentre(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  do j=1,3
                     if (molec%lOrient(j) /=0) then
                        nc=molec%lOrient(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do
                  !!! Falta THermal TLS !!!

               case (6)
                  !---- CEN ----!
                  do j=1,3
                     if (molec%lxcentre(j) /=0) then
                        nc=molec%lxcentre(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do

               case (7)
                  !---- ORI ----!
                  do j=1,3
                     if (molec%lOrient(j) /=0) then
                        nc=molec%lOrient(j)
                        call Delete_RefCodes(nc,molec)
                     end if
                  end do

               case (8)
                  !---- THE ----!
                  !!! Not Yet Implemented !!!
             end select

          !---- VARY Directive ----!
          case ("var")

             select case (key)
                case (0)

                   !---- Nb must be different zero ----!
                   select case (nb)
                      case (0)
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case ( 1:3)
                         !--- X_, Y_, Z_ ----!
                         if (molec%lI_coor(nb,na) ==0) then
                            molec%mI_coor(nb,na)=1.0
                            call get_atompos_ctr(molec%I_Coor(:,na),  &
                                                 Spg, np_refi,  &
                                                 molec%lI_coor(:,na), &
                                                 molec%mI_coor(:,na))
                            if (molec%lI_coor(nb,na) == np_refi) then
                               V_Vec(np_refi)=molec%I_Coor(nb,na)
                               V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                               V_Bounds(1,np_refi)=xl
                               V_Bounds(2,np_refi)=xu
                               V_Bounds(3,np_refi)=xs
                               V_BCon(np_refi)=ic
                               V_list(np_refi)=na
                            else
                               np_refi=np_refi-1
                            end if
                         end if

                      case ( 4)
                         !---- Biso_ ----!
                         if (molec%lbiso(na) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%biso(na)
                            V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                            molec%mbiso(na)=1.0
                            molec%lbiso(na)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 5)
                         !---- Occ_ ----!
                         if (molec%locc(na) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%occ(na)
                            V_Name(np_refi)=trim(code_nam(nb))//trim(molec%AtName(na))
                            molec%mocc(na)=1.0
                            molec%locc(na)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         end if

                      case ( 6:12)
                         !---- B11_, ..., B23_ ----!
                         err_refcodes=.true.
                         err_mess_refcodes="Option not defined"
                         return

                      case (13:15)
                         !---- Xc_, Yc_, Zc_ ----!
                         if (molec%lxcentre(nb-12) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%xcentre(nb-12)
                            V_Name(np_refi)=trim(code_nam(nb))//"Mol"
                            molec%mxcentre(nb-12)=1.0
                            molec%lxcentre(nb-12)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=0
                         end if

                      case (16:18)
                         !---- Theta_, Phi_, Chi_ ----!
                         if (molec%lOrient(nb-15) ==0) then
                            np_refi=np_refi+1
                            V_Vec(np_refi)=molec%orient(nb-15)
                            V_Name(np_refi)=trim(code_nam(nb))//"Mol"
                            molec%mOrient(nb-15)=1.0
                            molec%lOrient(nb-15)=np_refi
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=0
                            V_list(np_refi)=0
                         end if

                      case (19:21)
                         !!! Not yet implement !!!

                   end select ! nb

                case (1)
                   !---- XYZ ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) ==0) then
                         molec%mI_coor(j,na)=1.0
                         call get_atompos_ctr(molec%I_Coor(:,na),  &
                                              Spg, np_refi,  &
                                              molec%lI_coor(:,na), &
                                              molec%mI_coor(:,na))
                         if (molec%lI_coor(j,na) == np_refi) then
                            V_Vec(np_refi)=molec%I_Coor(j,na)
                            V_Name(np_refi)=trim(code_nam(j))//trim(molec%AtName(na))
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do

                case (2)
                   !---- OCC ----!
                   if (molec%locc(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%occ(na)
                      V_Name(np_refi)=trim(code_nam(5))//trim(molec%AtName(na))
                      molec%mocc(na)=1.0
                      molec%locc(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (3)
                   !---- BIS ----!
                   if (molec%lbiso(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%biso(na)
                      V_Name(np_refi)=trim(code_nam(4))//trim(molec%AtName(na))
                      molec%mbiso(na)=1.0
                      molec%lbiso(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if

                case (4)
                   !---- BAN ----!
                   err_refcodes=.true.
                   err_mess_refcodes="Option not defined"
                   return

                case (5)
                   !---- ALL ----!
                   do j=1,3
                      if (molec%lI_coor(j,na) ==0) then
                         molec%mI_coor(j,na)=1.0
                         call get_atompos_ctr(molec%I_Coor(:,na),  &
                                              Spg, np_refi,  &
                                              molec%lI_coor(:,na), &
                                              molec%mI_coor(:,na))
                         if (molec%lI_coor(j,na) == np_refi) then
                            V_Vec(np_refi)=molec%I_Coor(j,na)
                            V_Name(np_refi)=trim(code_nam(j))//trim(molec%AtName(na))
                            V_Bounds(1,np_refi)=xl
                            V_Bounds(2,np_refi)=xu
                            V_Bounds(3,np_refi)=xs
                            V_BCon(np_refi)=ic
                            V_list(np_refi)=na
                         else
                            np_refi=np_refi-1
                         end if
                      end if
                   end do
                   if (molec%locc(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%occ(na)
                      V_Name(np_refi)=trim(code_nam(5))//trim(molec%AtName(na))
                      molec%mocc(na)=1.0
                      molec%locc(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   if (molec%lbiso(na) ==0) then
                      np_refi=np_refi+1
                      V_Vec(np_refi)=molec%biso(na)
                      V_Name(np_refi)=trim(code_nam(4))//trim(molec%AtName(na))
                      molec%mbiso(na)=1.0
                      molec%lbiso(na)=np_refi
                      V_Bounds(1,np_refi)=xl
                      V_Bounds(2,np_refi)=xu
                      V_Bounds(3,np_refi)=xs
                      V_BCon(np_refi)=ic
                      V_list(np_refi)=na
                   end if
                   do j=1,3
                      if (molec%lxcentre(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%xcentre(j)
                         V_name(np_refi)=trim(code_nam(12+j))//"entre_Mol"
                         molec%mxcentre(j)=1.0
                         molec%lxcentre(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do
                   do j=1,3
                      if (molec%lOrient(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%orient(j)
                         V_name(np_refi)=trim(code_nam(15+j))//"Orient_Mol"
                         molec%mOrient(j)=1.0
                         molec%lOrient(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=0
                         V_list(np_refi)=0
                      end if
                   end do

                   !!! Falta THE !!!

                case (6)
                   !---- CEN ----!
                   do j=1,3
                      if (molec%lxcentre(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%xcentre(j)
                         V_name(np_refi)=trim(code_nam(12+j))//"Mol"
                         molec%mxcentre(j)=1.0
                         molec%lxcentre(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do

                case (7)
                   !---- ORI ----!
                   do j=1,3
                      if (molec%lOrient(j) ==0) then
                         np_refi=np_refi+1
                         V_Vec(np_refi)=molec%orient(j)
                         V_name(np_refi)=trim(code_nam(15+j))//"Mol"
                         molec%mOrient(j)=1.0
                         molec%lOrient(j)=np_refi
                         V_Bounds(1,np_refi)=xl
                         V_Bounds(2,np_refi)=xu
                         V_Bounds(3,np_refi)=xs
                         V_BCon(np_refi)=ic
                         V_list(np_refi)=0
                      end if
                   end do

                case (8)
                   !---- THE ----!

                   !!! Not yet implemented !!!
             end select
       end select

       return
    End Subroutine Fill_RefCodes_Molec

    !!----
    !!----  Subroutine Get_Atombet_Ctr(X,Betas,Spgr,Codini,Codes,Ord,Ss,Debug)
    !!----     real(kind=sp), dimension(3),     intent(in    ) :: x      ! Atom position (fractional coordinates)
    !!----     real(kind=sp), dimension(6),     intent(in out) :: betas  !Anisotropic temperature factors
    !!----     type(Space_Group_type),          intent(in    ) :: Spgr   !Space Group
    !!----     Integer,                         intent(in out) :: codini !Last attributed parameter
    !!----     real(kind=sp), dimension(6),     intent(in out) :: codes  !codewords for positions
    !!----     integer,               optional, intent(in    ) :: ord    !Order of the stabilizer
    !!----     integer, dimension(:), optional, intent(in    ) :: ss     !Pointer to SymmOp. of stabilizer
    !!----     integer,               optional, intent(in    ) :: debug  !Debug variable
    !!----
    !!----  Subroutine to get the appropriate constraints in the refinement codes of
    !!----  anisotropic atomic displacement(thermal) parameters.
    !!----
    !!----  Uses: Space_Group_type, get_stabilizer, sym_b_relations
    !!----
    !!----  Updated: March - 2005
    !!----
    !!
    Subroutine Get_Atombet_Ctr(X,Betas,Spgr,Codini,ICodes,Multip,Ord,Ss,Ipr)
       !---- Arguments ----!
       real(kind=sp), dimension(3),     intent(in    ) :: x
       real(kind=sp), dimension(6),     intent(in out) :: betas
       type(Space_Group_type),          intent(in    ) :: Spgr
       integer,                         intent(in out) :: codini
       integer,       dimension(6),     intent(in out) :: Icodes
       real(kind=sp), dimension(6),     intent(in out) :: Multip
       integer,               optional, intent(in    ) :: Ord
       integer, dimension(:), optional, intent(in    ) :: Ss
       integer,               optional, intent(in    ) :: Ipr

       !---- Local variables ----!
       character (len=1)               :: car
       character (len=1), dimension(6) :: cdd
       character (len=1), dimension(6) ::  strc

       integer                :: i, j, order
       integer, dimension(48) :: ss_ptr
       integer, dimension(6)  :: codd

       real(kind=sp), parameter     :: epss=0.01
       real(kind=sp), dimension(6)  :: cod, mul

       if (all(icodes ==0)) return

       cod=real(icodes)

       do j=1,6
          if (cod(j) < 1.0 .and. abs(multip(j)) > epss)  then
             codini=codini+1
             cod(j) = real(codini)
          end if
       end do

       if (present(ord) .and. present(ss)) then
          order=ord
          ss_ptr(1:order) = ss(1:ord)
       else
          call get_stabilizer(x,Spgr,order,ss_ptr)
       end if

       strc=(/"a","b","c","d","e","f"/)
       cdd = strc

       if (order > 1 ) then
          do j=2,order
             call sym_b_relations(Spgr%SymopSymb(ss_ptr(j)),codd,mul)
             do i=1,6
                if (abs(mul(i)) <= epss) then
                   cdd(i) = "0"
                   multip(i)=0.0
                else
                   cdd(i) = strc(codd(i))
                   multip(i)=mul(i)
                end if
             end do
             strc=cdd
          end do

          car=cdd(1)
          do j=2,6
             if (cdd(j) == "0") then
                cod(j) = 0.0
                betas(j)= 0.0
             end if
             if (cdd(j) == car) then
                cod(j) = cod(1)
                betas(j)= betas(1)*multip(j)
             end if
             if (cdd(j) == "a") then
                cod(j) = cod(1)
                betas(j)= betas(1)*multip(j)
             end if
             if (cdd(j) == "b") then
                cod(j) = cod(2)
                betas(j)= betas(2)*multip(j)
             end if
             if (cdd(j) == "c") then
                cod(j) = cod(3)
                betas(j)= betas(3)*multip(j)
             end if
             if (cdd(j) == "d") then
                cod(j) = cod(4)
                betas(j)= betas(4)*multip(j)
             end if
             if (cdd(j) == "e") then
                cod(j) = cod(5)
                betas(j)= betas(5)*multip(j)
             end if
          end do
       end if

       do j=1,6
          if (multip(j) < epss .or. cdd(j) == "0" ) then
             icodes(j) = 0
          else
             icodes(j) = nint(cod(j))
          end if
       end do

       if (present(ipr)) then
          write(unit=ipr,fmt="(a,6f10.4)")        "     Codes on Betas       : ",real(icodes)
          write(unit=ipr,fmt="(a,6(a,tr1),6f7.3)") "     Codes and multipliers:  ",cdd,multip
       end if

       return
    End Subroutine Get_Atombet_Ctr

    !!----
    !!----  Subroutine Get_Atompos_Ctr(X,Spgr,Codini,ICodes,Multip,Ord,Ss,Ipr)
    !!----     real(kind=sp), dimension(3),     intent(in    ) :: x      !Atom position (fractional coordinates)
    !!----     type(Space_Group_type),          intent(in    ) :: Spgr   !Space Group
    !!----     Integer,                         intent(in out) :: codini !Last attributed parameter
    !!----     integer,       dimension(3),     intent(in out) :: Icodes
    !!----     real(kind=sp), dimension(3),     intent(in out) :: Multip
    !!----     integer,               optional, intent(in    ) :: Ord
    !!----     integer, dimension(:), optional, intent(in    ) :: Ss
    !!----     integer,               optional, intent(in    ) :: Ipr
    !!----
    !!----  Subroutine to get the appropriate constraints in the refinement codes of
    !!----  atoms positions.
    !!----
    !!----  Updated: March - 2005
    !!----
    !!
    Subroutine Get_Atompos_Ctr(X,Spgr,Codini,ICodes,Multip,Ord,Ss,Ipr)
       !---- Arguments ----!
       real(kind=sp), dimension(3),     intent(in    ) :: x
       type(Space_Group_type),          intent(in    ) :: Spgr
       integer,                         intent(in out) :: Codini
       integer,       dimension(3),     intent(in out) :: Icodes
       real(kind=sp), dimension(3),     intent(in out) :: Multip
       integer,               optional, intent(in    ) :: Ord
       integer, dimension(:), optional, intent(in    ) :: Ss
       integer,               optional, intent(in    ) :: Ipr

       !---- Local variables ----!
       real, parameter        :: epss=0.001
       integer                :: j=0, order=0, L=0, L1=0, L2=0, jx=0, ii=0, m=0, ipar=0
       integer, dimension(48) :: ss_ptr
       integer, dimension(3)  :: cdd
       real(kind=sp),    dimension(3)   :: cod
       character (len=40)               :: symbol
       character (len=3),dimension(0:12):: car=(/"  0","  a","  b","  c", &  ! 0     1     2     3
                                                       " -a"," -b"," -c", &  !       4     5     6
                                                       " 2a"," 2b"," 2c", &  !       7     8     9
                                                       "a/2","b/2","c/2"/)   !      10    11    12


       if (present(ord) .and. present(ss)) then
          order=ord
          ss_ptr(1:order) = ss(1:ord)
       else
          call get_stabilizer(x,Spgr,order,ss_ptr)
       end if

       cdd = (/1,2,3/)
       cod=real(icodes)              !Parameter number with sign
       do j=1,3
          if (abs(cod(j)) < 1.0 .and. multip(j) > epss)  then
             codini=codini+1
             cod(j) = real(codini)
          end if
       end do

       if (present(ipr)) then
          write(unit=ipr,fmt="(/,a,3f10.5)") "     Atom Position:",x
          if (order > 1 ) write(unit=98,fmt="( a)")        "     List of symmetry element of the stabilizer:"
       end if

       if (order > 1 ) then
          do j=2,order
             symbol=" "
             call symmetry_symbol(Spgr%SymOp(ss_ptr(j)),symbol)
             ipar=index(symbol,")")
             L =index(symbol(ipar+1:)," ")+ipar
             L1=index(symbol(ipar+1:),",")+ipar
             L2=index(symbol(L1+1:),",")+L1
             if (L1 == 0) L1=1
             if (L2 == 0) L2=1
             if (L  == 0) L=1

             !---- Test fixed parameters on x,y,z-positions ----!
             if (index(symbol(L:L1),"x") == 0 ) then
                cdd(1) = 0
                cod(1) = 0.0
                multip(1) = 0.0
             end if
             if (index(symbol(L1:L2),"x") == 0 .and. index(symbol(L1:L2),"y") == 0 ) then
                cdd(2)=0
                cod(2) = 0.0
                multip(2) = 0.0
             end if
             if (index(symbol(L2:),"x") == 0 .and. index(symbol(L2:),"y") == 0 .and. &
                 index(symbol(L2:),"z") == 0 ) then
                cdd(3)=0
                cod(3) = 0.0
                multip(3) = 0.0
             end if

             !---- Test x on y-position ----!
             if (cdd(2) == 2) then
                ii=index(symbol(L1+1:L2),"x")
                m=ii+L1-1  !absolute position of x - 1
                if (ii /= 0 .and. multip(1) > epss .and. m > 0) then
                   select case (symbol(m:m))
                      case("-")
                         cod(2) = -cod(1)
                         multip(2) =  multip(1)
                         cdd(2) = 4 !-a

                      case(",")
                         jx=index(symbol(L+1:L1),"x")
                         jx=jx+L-1
                         if (jx > 0) then
                            if (symbol(jx:jx) == "2") then
                               cod(2) = cod(1)
                               multip(2) = 0.5*multip(1)
                               cdd(2) = 10  !a/2
                            else
                               cod(2) = cod(1)
                               multip(2) = multip(1)
                               cdd(2) = 1  !a
                            end if
                         end if

                      case ("2")
                         cod(2) = cod(1)
                         multip(2) = 2.0*multip(1)
                         cdd(2) = 7   !2a
                   end select

                else if (ii /= 0) then
                   cod(2) = 0.0
                   multip(2) = 0.0
                   cdd(2) = 0      !0
                end if
             end if  ! cdd(2)==2

             !---- Test x on z-position ----!
             if (cdd(3) == 3) then
                ii=index(symbol(L2+1:),"x")
                m=ii+L2-1
                if (ii /= 0 .and. multip(1) > epss .and. m > 0) then
                   select case (symbol(m:m))
                      case ("-")
                         cod(3) = -cod(1)
                         multip(3) =  multip(1)
                         cdd(3) = 4    !-a

                      case (",")
                         jx=index(symbol(L+1:L1),"x")
                         jx=jx+L-1
                         if (jx > 0) then
                            if (symbol(jx:jx) == "2") then
                               cod(3) = cod(1)
                               multip(3) = 0.5*multip(1)
                               cdd(3) = 10 !a/2
                            else
                               cod(3) = cod(1)
                               multip(3) = multip(1)
                               cdd(3) = 1   !a
                            end if
                         end if

                      case ("2")
                         cod(3) = cod(1)
                         multip(3) = 2.0*multip(1)
                         cdd(3) = 7  !2a
                   end select

                else if (ii /= 0 ) then
                   cod(3) = 0.0
                   multip(3) = 0.0
                   cdd(3) = 0    !0
                end if
             end if ! cdd(3)==0

             !---- Test y on z-position ----!
             if (cdd(3) == 3) then
                ii=index(symbol(L2+1:),"y")
                m=ii+L2-1
                if (ii /= 0 .and. multip(2) > epss .and. m > 0) then
                   select case (symbol(m:m))
                      case ("-")
                         cod(3) = -cod(2)
                         multip(3) =  multip(2)
                         cdd(3) = 5    !-b

                      case (",")
                         jx=index(symbol(L1:L2),"y")
                         jx=jx+L1-1
                         if (j > 0) then
                            if (symbol(jx:jx) == "2") then
                               cod(3) = cod(2)
                               multip(3) = 0.5*multip(2)
                               cdd(3) = 11   !b/2
                            else
                               cod(3) = cod(2)
                               multip(3) = multip(2)
                               cdd(3) = 2   !b
                            end if
                         end if

                      case ("2")
                         cod(3) = cod(2)
                         multip(3) = 2.0*multip(2)
                         cdd(3) = 8  !2b
                   end select

                else if (ii /= 0)  then
                   cod(3) = 0.0
                   multip(3) = 0.0
                   cdd(3) = 0   !0
                end if
             end if ! cdd(3)==0

             if (present(ipr)) then
                write(unit=98,fmt="(a,i2,a,t20,a,t55,a,t90,a,7i4)") "     Operator ",j,": ", &
                     trim(Spgr%SymopSymb(ss_ptr(j))),trim(symbol),"  ipar, L,L1,L2,ii,m,jx:" ,ipar,L,L1,L2,ii,m,jx
             end if
          end do !do j=1,order  over operators of the stabilizer
       end if  !order > 1

       do j=1,3
          if (multip(j) < epss) then
             icodes(j) = 0
          else
             icodes(j)=nint(cod(j))
          end if
       end do

       if (present(ipr)) then
          write(unit=ipr,fmt="(a,3f10.4)") "     Codes positions: ",real(icodes)
          write(unit=ipr,fmt="(5a)")       "     Codes   string : ( ",(car(cdd(j)),j=1,3) ," )"
       end if

       return
    End Subroutine Get_Atompos_Ctr

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_FAtom(Line,FAtom)
    !!--++    character(len=*),         intent(in)     :: Line
    !!--++    integer,                  intent(in)     :: Nat
    !!--++    type(Atom_List_Type),     intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_FAtom(Line,FAtom)
       !---- Arguments ----!
       character(len=*),     intent(in)     :: Line
       type(Atom_List_Type), intent(in out) :: FAtom

       !---- Loval variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: j,ic,n,na,nb,nc,nd,npos
       integer                          :: nl,nl2,iv
       integer, dimension(1)            :: ivet
       real                             :: fac_0,fac_1
       real,dimension(1)                :: vet

       call init_err_refcodes()

       nl=0
       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          err_mess_refcodes="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Father Information----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter

       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          err_mess_refcodes="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if

       na=0
       do j=1,FAtom%Natoms
          if (u_case(FAtom%atom(j)%lab) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do
       if (na == 0) then
          err_refcodes=.true.
          err_mess_refcodes="Atom label not found for "//trim(line)
          return
       end if

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          err_mess_refcodes="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       select case (nb)
          case ( 1:3)
             fac_0=FAtom%atom(na)%mx(nb)
                nl=FAtom%atom(na)%lx(nb)
          case ( 4)
             fac_0=FAtom%atom(na)%mbiso
                nl=FAtom%atom(na)%lbiso
          case ( 5)
             fac_0=FAtom%atom(na)%mocc
                nl=FAtom%atom(na)%locc
          case ( 6:11)
             fac_0=FAtom%atom(na)%mu(nb-5)
                nl=FAtom%atom(na)%lu(nb-5)
          case (12)
             fac_0=FAtom%atom(na)%mu(1)
          case (13:)
             err_refcodes=.true.
             err_mess_refcodes="Incompatible Code-name for parameter name: "//trim(label(1))
             return
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             err_mess_refcodes="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set the rest elements in Contsraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             err_mess_refcodes="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          nc=0
          do j=1,FAtom%Natoms
             if (u_case(FAtom%atom(j)%lab) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do
          if (nc == 0) then
             err_refcodes=.true.
             err_mess_refcodes="Atom label not found for "//trim(label(n))
             return
          end if

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             err_mess_refcodes="Code-name not found for "//trim(label(n))
             return
          end if

          !---- Is there a new multiplier?
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          select case (nd)
             case ( 1:3)
                nl2=FAtom%atom(nc)%lx(nd)
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mx(nd)=fac_1
                FAtom%atom(nc)%lx(nd)=nl
             case ( 4)
                nl2=FAtom%atom(nc)%lbiso
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mbiso=fac_1
                FAtom%atom(nc)%lbiso=nl
             case ( 5)
                nl2=FAtom%atom(nc)%locc
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mocc=fac_1
                FAtom%atom(nc)%locc=nl
             case ( 6:11)
                nl2=FAtom%atom(nc)%lu(nd-5)
                call Delete_refCodes(nl2,FAtom)
                FAtom%atom(nc)%mu(nd-5)=fac_1
                FAtom%atom(nc)%lu(nd-5)=nl
             case (12)
                do j=1,6
                   nl2=FAtom%atom(nc)%lu(j)
                   call Delete_refCodes(nl2,FAtom)
                   FAtom%atom(nc)%mu(j)=fac_1
                   FAtom%atom(nc)%lu(j)=FAtom%atom(na)%lu(j)
                end do
                np_cons=np_cons+5
             case (13:)
                err_refcodes=.true.
                err_mess_refcodes="Incompatible Code-name for parameter name: "//trim(label(1))
                return
          end select ! nb

          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_FAtom

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_Molcrys(Line,Molcrys)
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    integer,                      intent(in)     :: NMol
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_Molcrys(Line,Molcrys)
       !---- Arguments ----!
       character(len=*),             intent(in)     :: Line
       type(molecular_Crystal_type), intent(in out) :: MolCrys

       !---- Loval variables ----!
       character(len=5)                 :: car
       character(len=20), dimension(30) :: label
       integer                          :: i,j,ic,n,na,naa,nb,nc,ncc,nd
       integer                          :: npos, nposm, nmol1,nmol2
       integer                          :: nl,nl2,iv
       integer, dimension(1)            :: ivet
       real                             :: fac_0,fac_1
       real,dimension(1)                :: vet

       call init_err_refcodes()

       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          err_mess_refcodes="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Father Information ----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter
       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          err_mess_refcodes="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if
       nposm=index(u_case(label(1)),"MOL")
       if (nposm /= 0) then
          car=adjustl(label(1)(nposm+3:))
          if (car(1:2) == "  ") then
             nmol1=0
          else
             read(unit=car,fmt="(i2)") nmol1
          end if
       else
          nmol1=-1
       end if

       na=0
       do j=1,molcrys%n_free
          if (u_case(molcrys%atm(j)%lab) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do
       if (na == 0) then
          do i=1,molcrys%n_mol
             do j=1,molcrys%mol(i)%natoms
                if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(1)(npos+1:npos+6))) then
                   if (j > 1) then
                      na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                   else
                      na=molcrys%n_free+j
                   end if
                   exit
                end if
             end do
          end do
       end if

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          err_mess_refcodes="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       !---- Checking ----!
       if (nb <= 12 .and. na==0) then
          err_refcodes=.true.
          err_mess_refcodes="Incompatible option for parameter name: "//trim(label(1))
          return
       end if

       nl=0
       select case (nb)
          case ( 1:3)
             !---- X_, Y_, Z_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mx(nb)
                   nl=molcrys%atm(na)%lx(nb)
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mI_coor(nb,naa)
                      nl=molcrys%mol(i)%lI_coor(nb,naa)
                end do
             end if

          case ( 4)
             !---- Biso_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mbiso
                   nl=molcrys%atm(na)%lbiso
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mbiso(naa)
                      nl=molcrys%mol(i)%lbiso(naa)
                end do
             end if

          case ( 5)
             !---- Occ_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mocc
                   nl=molcrys%atm(na)%locc
             else
                naa=na-molcrys%n_free
                do i=1,molcrys%n_mol
                   if (naa > molcrys%mol(i)%natoms) then
                      naa=naa-molcrys%mol(i)%natoms
                      cycle
                   end if
                   fac_0=molcrys%mol(i)%mocc(naa)
                      nl=molcrys%mol(i)%locc(naa)
                end do
             end if

          case ( 6:11)
             !---- B11_, ..., B23_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mu(nb-5)
                   nl=molcrys%atm(na)%lu(nb-5)
             else
                err_refcodes=.true.
                err_mess_refcodes="Option no valid"
                return
             end if

          case (12)
             !---- Banis_ ----!
             if (na <= molcrys%n_free) then
                fac_0=molcrys%atm(na)%mu(1)
             else
                err_refcodes=.true.
                err_mess_refcodes="Option no valid"
                return
             end if

          case (13:15)
             !---- Xc_, Yc_, Zc_ ----!
             select case (nmol1)
                case (-1)
                   err_refcodes=.true.
                   err_mess_refcodes="Option no valid"
                   return

                case (0)
                   fac_0=molcrys%mol(1)%mxcentre(nb-12)
                   nl=molcrys%mol(1)%lxcentre(nb-12)

                case (1:)
                   fac_0=molcrys%mol(nmol1)%mxcentre(nb-12)
                   nl=molcrys%mol(nmol1)%lxcentre(nb-12)
             end select

          case (16:18)
             !---- Theta_ , Phi_, Chi_ ----!
             select case (nmol1)
                case (-1)
                   err_refcodes=.true.
                   err_mess_refcodes="Option no valid"
                   return

                case (0)
                   fac_0=molcrys%mol(1)%mOrient(nb-15)
                   nl=molcrys%mol(1)%lOrient(nb-15)

                case (1)
                   fac_0=molcrys%mol(nmol1)%mOrient(nb-15)
                   nl=molcrys%mol(nmol1)%lOrient(nb-15)
             end select

          case (19:21)
             !!! Not yet Implemented !!!
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             err_mess_refcodes="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set the rest elements in Constraints ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             err_mess_refcodes="No CrysFML code-name was found for "//trim(label(n))
             return
          end if
          nposm=index(u_case(label(n)),"MOL")
          if (nposm /= 0) then
             car=adjustl(label(n)(nposm+3:))
             if (car(1:2) == "  ") then
                nmol2=0
             else
                read(unit=car,fmt="(i2)") nmol2
             end if
          else
             nmol2=-1
          end if

          nc=0
          do j=1,molcrys%n_free
             if (u_case(molcrys%atm(j)%lab) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do
          if (nc == 0) then
             do i=1,molcrys%n_mol
                do j=1,molcrys%mol(i)%natoms
                   if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(n)(npos+1:npos+6))) then
                      if (j > 1) then
                         nc=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                      else
                         nc=molcrys%n_free+j
                      end if
                      exit
                   end if
                end do
             end do
          end if

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             err_mess_refcodes="Code-name not found for "//trim(label(n))
             return
          end if

          !---- Checking ----!
          if (nd <= 12 .and. nc==0) then
             err_refcodes=.true.
             err_mess_refcodes="Incompatible option for parameter name: "//trim(label(n))
             return
          end if

          !---- Is there a new multiplier? ----!
          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          select case (nd)
             case ( 1:3)
                !---- X_, Y_, Z_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lx(nd)
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mx(nd)=fac_1
                   molcrys%atm(nc)%lx(nd)=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%lI_coor(nd,ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mI_coor(nd,ncc)=fac_1
                      molcrys%mol(i)%lI_coor(nd,ncc)=nl
                   end do
                end if

             case ( 4)
                !---- Biso_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lbiso
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mbiso=fac_1
                   molcrys%atm(nc)%lbiso=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%lbiso(ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mbiso(ncc)=fac_1
                      molcrys%mol(i)%lbiso(ncc)=nl
                   end do
                end if

             case ( 5)
                !---- Occ_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%locc
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mocc=fac_1
                   molcrys%atm(nc)%locc=nl
                else
                   ncc=nc-molcrys%n_free
                   do i=1,molcrys%n_mol
                      if (ncc > molcrys%mol(i)%natoms) then
                         ncc=ncc-molcrys%mol(i)%natoms
                         cycle
                      end if
                      nl2=molcrys%mol(i)%locc(ncc)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(i)%mocc(ncc)=fac_1
                      molcrys%mol(i)%locc(ncc)=nl
                   end do
                end if

             case ( 6:11)
                !---- B11_, ...., B23_ ----!
                if (nc <= molcrys%n_free) then
                   nl2=molcrys%atm(nc)%lu(nd-5)
                   call Delete_RefCodes(nl2,molcrys)
                   molcrys%atm(nc)%mu(nd-5)=fac_1
                   molcrys%atm(nc)%lu(nd-5)=nl
                else
                   err_refcodes=.true.
                   err_mess_refcodes="Option no valid"
                   return
                end if

             case (12)
                !---- Banis_ ----!
                if (nc <= molcrys%n_free) then
                   do j=1,6
                      nl2=molcrys%atm(nc)%lu(j)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%atm(nc)%mu(j)=fac_1
                      molcrys%atm(nc)%lu(j)=molcrys%atm(na)%lu(j)
                   end do
                   np_cons=np_cons+5
                else
                   err_refcodes=.true.
                   err_mess_refcodes="Option no valid"
                   return
                end if

             case (13:15)
                select case (nmol2)
                   case (-1)
                      err_refcodes=.true.
                      err_mess_refcodes="Option no valid"
                      return

                   case (0)
                      do i=1,molcrys%n_mol
                         nl2=molcrys%mol(i)%lxcentre(nc-12)
                         call Delete_RefCodes(nl2,molcrys)
                         molcrys%mol(i)%mxcentre(nc-12)=fac_1
                         molcrys%mol(i)%lxcentre(nc-12)=nl
                      end do
                      np_cons=np_cons+(i-1)

                   case (1:)
                      nl2=molcrys%mol(nmol2)%lxcentre(nc-12)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(nmol2)%mxcentre(nc-12)=fac_1
                      molcrys%mol(nmol2)%lxcentre(nc-12)=nl
                end select

             case (16:18)
                select case (nmol2)
                   case (-1)
                      err_refcodes=.true.
                      err_mess_refcodes="Option no valid"
                      return

                   case (0)
                      do i=1,molcrys%n_mol
                         nl2=molcrys%mol(i)%lorient(nc-15)
                         call Delete_RefCodes(nl2,molcrys)
                         molcrys%mol(i)%morient(nc-15)=fac_1
                         molcrys%mol(i)%lorient(nc-15)=nl
                      end do
                      np_cons=np_cons+(i-1)

                   case (1:)
                      nl2=molcrys%mol(nmol2)%lorient(nc-15)
                      call Delete_RefCodes(nl2,molcrys)
                      molcrys%mol(nmol2)%morient(nc-15)=fac_1
                      molcrys%mol(nmol2)%lorient(nc-15)=nl
                end select

             case (19:21)
                !!! Not yet implemented !!!

          end select ! nb
          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_Molcrys

    !!--++
    !!--++ Subroutine Get_ConCodes_Line_Molec(Line,Molec)
    !!--++    character(len=*),    intent(in)     :: Line
    !!--++    type(molecule_type), intent(in out) :: Molec
    !!--++
    !!--++ Overloaded
    !!--++ Get the Constraints relations
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_ConCodes_Line_Molec(Line,Molec)
       !---- Arguments ----!
       character(len=*),    intent(in)     :: Line
       type(molecule_type), intent(in out) :: Molec

       !---- Loval variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: j,ic,na,nb,nc,nd,npos!,i,naa,ncc
       integer                          :: n,nl,nl2,iv
       integer, dimension(1)            :: ivet
       real                             :: fac_0,fac_1
       real,dimension(1)                :: vet

       call init_err_refcodes()

       call getword(line,label,ic)
       if (ic < 2) then
          err_refcodes=.true.
          err_mess_refcodes="EQUAL keyword needs two labels: "//trim(line)
          return
       end if

       !---- Set Father Information ----!
       !---- Na is the number of atom on List
       !---- Nb is the key (X,Y,Z,Occ,...)
       !---- Fac0 is the multiplier
       !---- Nl is the number of refinement parameter
       npos=index(label(1),"_")
       if (npos ==0) then
          err_refcodes=.true.
          err_mess_refcodes="The name "//trim(label(1))//" does not fit any known code-name of CrysFML "
          return
       end if

       na=0
       do j=1,molec%natoms
          if (u_case(molec%Atname(j)) == u_case(label(1)(npos+1:npos+6))) then
             na=j
             exit
          end if
       end do

       nb=0
       do j=1,ncode
          if (u_case(label(1)(1:npos))==u_case(trim(code_nam(j)))) then
             nb=j
             exit
          end if
       end do
       if (nb == 0) then
          err_refcodes=.true.
          err_mess_refcodes="Code-name not found for parameter name: "//trim(label(1))
          return
       end if

       !---- Checking ----!
       if (nb < 6 .and. na == 0) then
          err_refcodes=.true.
          err_mess_refcodes="Incompatible relation: "//trim(label(1))
          return
       end if

       nl=0
       select case (nb)
          case ( 1:3)
             !---- X_, Y_, Z_ ----!
             fac_0=molec%mI_coor(nb,na)
                nl=molec%lI_coor(nb,na)
          case ( 4)
             !---- Biso_ ----!
             fac_0=molec%mbiso(na)
                nl=molec%lbiso(na)
          case ( 5)
             !---- Occ_ ----!
             fac_0=molec%mocc(na)
                nl=molec%locc(na)
          case ( 6:12)
             !---- Anisotropic Parameters ----!
             err_refcodes=.true.
             err_mess_refcodes="Incompatible Code-name for parameter name: "//trim(label(1))
             return
          case (13:15)
             !---- Xc_, Yc_, Zc_ ----!
             fac_0=molec%mxcentre(nb-12)
                nl=molec%lxcentre(nb-12)
          case (16:18)
             !---- Theta_, Phi_, Chi_ ----!
             fac_0=molec%mxcentre(nb-15)
                nl=molec%lxcentre(nb-15)
          case (19:21)
             !!! Not Yet Implemented !!!
       end select ! nb

       if (nb < ncode) then
          if (nl == 0) then
             err_refcodes=.true.
             err_mess_refcodes="No refinable parameter was selected for "//trim(label(1))
             return
          end if
       end if

       !---- Set Others ----!
       n=1
       do
          n=n+1
          if (n > ic) exit

          npos=index(label(n),"_")
          if (npos ==0) then
             err_refcodes=.true.
             err_mess_refcodes="No CrysFML code-name was found for "//trim(label(n))
             return
          end if

          nc=0
          do j=1,molec%natoms
             if (u_case(molec%Atname(j)) == u_case(label(n)(npos+1:npos+6))) then
                nc=j
                exit
             end if
          end do

          nd=0
          do j=1,ncode
             if (u_case(label(n)(1:npos))==u_case(trim(code_nam(j)))) then
                nd=j
                exit
             end if
          end do
          if (nd == 0) then
             err_refcodes=.true.
             err_mess_refcodes="Code-name not found for "//trim(label(n))
             return
          end if

          n=n+1
          call getnum(label(n),vet,ivet,iv)
          if (iv == 1) then
             fac_1=vet(1)
          else
             fac_1=fac_0
             n=n-1
          end if

          !---- Checking ----!
          if (nd < 6 .and. nc == 0) then
             err_refcodes=.true.
             err_mess_refcodes="Incompatible relation: "//trim(label(n))
             return
          end if

          select case (nd)
             case ( 1:3)
                !---- X_, Y_, Z_ ----!
                nl2=molec%lI_coor(nd,nc)
                call Delete_RefCodes(nl2,molec)
                molec%mI_coor(nd,nc)=fac_1
                molec%lI_coor(nd,nc)=nl
             case ( 4)
                !---- Biso_ ----!
                nl2=molec%lbiso(nc)
                call Delete_RefCodes(nl2,molec)
                molec%mbiso(nc)=fac_1
                molec%lbiso(nc)=nl
             case ( 5)
                !---- Occ_ ----!
                nl2=molec%locc(nc)
                call Delete_RefCodes(nl2,molec)
                molec%mocc(nc)=fac_1
                molec%locc(nc)=nl
             case ( 6:12)
                err_refcodes=.true.
                err_mess_refcodes="Incompatible Code-name for "//trim(label(n))
                return
             case (13:15)
                !---- Xc_, Yc_, Zc_ ----!
                nl2=molec%lxcentre(nc-12)
                call Delete_RefCodes(nl2,molec)
                molec%mxcentre(nc-12)=fac_1
                molec%lxcentre(nc-12)=nl
             case (16:18)
                !---- Theta_, Phi_, Chi_ ----!
                nl2=molec%lorient(nc-15)
                call Delete_RefCodes(nl2,molec)
                molec%morient(nc-15)=fac_1
                molec%lorient(nc-15)=nl
             case (19:21)
                !!! not yet implemented !!!
          end select ! nb
          np_cons=np_cons+1

       end do

       return
    End Subroutine Get_ConCodes_Line_Molec

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_FAtom(Key,Dire,Line,FAtom,Spg)
    !!--++    integer,                 intent(in)     :: Key
    !!--++    character(len=*),        intent(in)     :: Dire
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++    type(space_group_type),  intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_FAtom(Key,Dire,Line,FAtom,Spg)
       !---- Arguments ----!
       integer,                 intent(in)     :: Key
       character(len=*),        intent(in)     :: Dire
       character(len=*),        intent(in)     :: Line
       type(Atom_List_Type),    intent(in out) :: FAtom
       type(space_group_type),  intent(in)     :: Spg

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,n,na,nb,ndir,npos,nlong,ic !,k,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real                             :: x_low,x_up,x_step
       real,dimension(5)                :: vet


       call init_err_refcodes()

       nlong=len_trim(line)

       if (nlong ==0) then
          !---- Default Values ----!
          do i=1,FAtom%natoms
             call Fill_RefCodes(Key,Dire,i,0,0.0,0.0,0.0,0,Fatom,Spg)
          end do

       else
          !---- VARY/FIX Line: [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !--- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   err_mess_refcodes="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             do na=1,FAtom%natoms
                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,fatom,spg)
             end do
             if (err_refcodes) return

          else
             !---- [LABEL, [INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Default values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      err_mess_refcodes="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                do j=1,FAtom%natoms
                   if (u_case(fatom%atom(j)%lab) == u_case(label(n_ini)(npos+1:npos+6))) then
                      na=j
                      exit
                   end if
                end do
                if (na == 0) then
                   err_refcodes=.true.
                   err_mess_refcodes="Atom label not found for "//trim(line)
                   return
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,fatom,spg)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefCodes_Line_FAtom

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_Molcrys(Key,Dire,Line,Molcrys,NMol)
    !!--++   integer,                      intent(in)     :: Key
    !!--++   character(len=*),             intent(in)     :: Dire
    !!--++   character(len=*),             intent(in)     :: Line
    !!--++   type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++   integer,                      intent(in)     :: NMol
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_Molcrys(Key,Dire,Line,Molcrys,NMol)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       character(len=*),             intent(in)     :: Line
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       integer,                      intent(in)     :: NMol

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,k,n,na,nb,nc,ndir,npos,nlong,ic
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real                             :: x_low,x_up,x_step
       real,dimension(5)                :: vet


       call init_err_refcodes()

       nlong=len_trim(line)
       if (nlong ==0) then
          !---- Default values ----!
          select case (NMol)
             case (-1)
                !---- No Molecule Information ----!
                do i=1,molcrys%n_free
                   call Fill_RefCodes(Key,Dire,i,0,0.0,0.0,0.0,0,Molcrys,NMol)
                end do

             case (0)
                do k=1,molcrys%n_mol
                   do i=1,molcrys%mol(k)%natoms
                      if (k > 1) then
                         na=molcrys%n_free+sum(molcrys%mol(1:k-1)%natoms)+i
                      else
                         na=molcrys%n_free+i
                      end if
                      call Fill_RefCodes(Key,Dire,na,0,0.0,0.0,0.0,0,Molcrys,NMol)
                   end do
                end do

             case (1:)
                do i=1,molcrys%mol(nmol)%natoms
                   if (nmol > 1) then
                      na=molcrys%n_free+sum(molcrys%mol(1:nmol-1)%natoms)+i
                   else
                      na=molcrys%n_free+i
                   end if
                   call Fill_RefCodes(Key,Dire,na,0,0.0,0.0,0.0,0,Molcrys,NMol)
                end do
          end select

       else
          !---- VARY/FIX Line: [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !---- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   err_mess_refcodes="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             select case (nmol)
                case (-1)
                   !---- No Molecule Information ----!
                   do na=1,molcrys%n_free
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do

                case ( 0)
                   !---- For all Molecules Defined ----!
                   do n=1,molcrys%n_mol
                      do i=1,molcrys%mol(n)%natoms
                         if (n > 1) then
                            na=molcrys%n_free+sum(molcrys%mol(1:n-1)%natoms)+i
                         else
                            na=molcrys%n_free+i
                         end if
                      end do
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do

                case (1:)
                   !---- Particular molecule ----!
                   do i=1,molcrys%mol(nmol)%natoms
                      if (nmol > 1) then
                         na=molcrys%n_free+sum(molcrys%mol(1:nmol-1)%natoms)+i
                      else
                         na=molcrys%n_free+i
                      end if
                      call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,Molcrys,NMol)
                   end do
             end select
             if (err_refcodes) return

          else
             !---- [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Deafult values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      err_mess_refcodes="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                select case (nmol)
                   case (-1)
                      !---- No molecule Information ----!
                      do j=1,molcrys%n_free
                         if (u_case(molcrys%atm(j)%lab) == u_case(label(n_ini)(npos+1:npos+6))) then
                            na=j
                            exit
                         end if
                      end do

                   case ( 0)
                      !---- For all molecules defined ----!
                      do i=1,molcrys%n_mol
                         do j=1,molcrys%mol(i)%natoms
                            if (u_case(molcrys%mol(i)%Atname(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                               if (j > 1) then
                                  na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                               else
                                  na=molcrys%n_free+j
                               end if
                               exit
                            end if
                         end do
                      end do

                   case (1:)
                      !---- Particular Molecule ----!
                      do j=1,molcrys%mol(nmol)%natoms
                         if (u_case(molcrys%mol(nmol)%Atname(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                            if (j > 1) then
                               na=molcrys%n_free+sum(molcrys%mol(1:j-1)%natoms)+j
                            else
                               na=molcrys%n_free+j
                            end if
                            exit
                         end if
                      end do
                end select
                if (na == 0) then
                   err_refcodes=.true.
                   err_mess_refcodes="Atom label not found for "//trim(line)
                   return
                end if

                !---- Valu List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molcrys,NMol)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if
       end if

       return
    End Subroutine Get_RefCodes_Line_Molcrys

    !!--++
    !!--++ Subroutine Get_RefCodes_Line_Molec(Key,Dire,Line,Molec,Spg)
    !!--++    integer,                      intent(in)     :: Key
    !!--++    character(len=*),             intent(in)     :: Dire
    !!--++    character(len=*),             intent(in)     :: Line
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    type(space_group_type),       intent(in)     :: Spg
    !!--++
    !!--++ Overloaded
    !!--++ Get Refinement Codes for Free atoms type
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RefCodes_Line_Molec(Key,Dire,Line,Molec,Spg)
       !---- Arguments ----!
       integer,                      intent(in)     :: Key
       character(len=*),             intent(in)     :: Dire
       character(len=*),             intent(in)     :: Line
       type(molecule_type),          intent(in out) :: Molec
       type(space_group_type),       intent(in)     :: Spg

       !---- Local Variables ----!
       character(len=20), dimension(30) :: label
       integer                          :: i,j,n,na,nb,ndir,npos,nlong,ic !,k,nc
       integer                          :: icond,iv,n_ini,n_end
       integer, dimension(5)            :: ivet
       integer, dimension(30)           :: ilabel
       real                             :: x_low,x_up,x_step
       real,dimension(5)                :: vet

       call init_err_refcodes()

       nlong=len_trim(line)

       if (nlong ==0) then
          !---- Default values ----!
          do i=1,molec%natoms
             call Fill_RefCodes(key,dire,i,0,0.0,0.0,0.0,0,molec,spg)
          end do

       else
          !---- VARY/FIX Line: [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
          ilabel=0
          call getword(line,label,ic)
          do i=1,ic
             call getnum(label(i),vet,ivet,iv)
             if (iv == 1) ilabel(i)=1
          end do
          ndir=count(ilabel(1:ic) < 1)

          if (ndir <=0) then
             !---- [INF,[SUP,[STEP,[COND]]]] ----!
             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (1)
                   x_low=vet(1)
                    x_up=vet(1)
                  x_step=0.0
                   icond=0

                case (2)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=0.0
                   icond=0

                case (3)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=0

                case (4)
                   x_low=vet(1)
                    x_up=vet(2)
                  x_step=vet(3)
                   icond=ivet(4)

                case default
                   err_refcodes=.true.
                   err_mess_refcodes="Only numbers in "//trim(line)
                   return
             end select

             nb=0
             do na=1,molec%natoms
                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molec,spg)
             end do
             if (err_refcodes) return

          else
             !---- [LABEL,[INF,[SUP,[STEP,[COND]]]]] ----!
             ! If Ilabel(i) == 0 then is a label otherwise is a number
             n_ini=minloc(ilabel,dim=1)
             ilabel(n_ini)=2
             n_end=minloc(ilabel,dim=1)-1

             do n=1,ndir
                na=0
                nb=0

                !---- Deafult values ----!
                x_low =0.0
                x_up  =0.0
                x_step=0.0
                icond =0

                !---- Label ----!
                npos=index(label(n_ini),"_")
                if (npos >0) then
                   do j=1,ncode
                      if (u_case(label(n_ini)(1:npos))==u_case(trim(code_nam(j)))) then
                         nb=j
                         exit
                      end if
                   end do
                   if (nb == 0) then
                      err_refcodes=.true.
                      err_mess_refcodes="Code-name not found for "//trim(label(n_ini))
                      return
                   end if
                end if

                do j=1,molec%natoms
                   if (u_case(molec%AtName(j)) == u_case(label(n_ini)(npos+1:npos+6))) then
                      na=j
                      exit
                   end if
                end do
                if (na == 0) then
                   err_refcodes=.true.
                   err_mess_refcodes="Atom label not found for "//trim(line)
                   return
                end if

                !---- Value List: Inf,Sup,Step,Cond ----!
                i=0
                do j=n_ini+1,n_end
                   i=i+1
                   call getnum(label(j),vet,ivet,iv)
                   select case (i)
                      case (1)
                         x_low=vet(1)
                         x_up =vet(1)

                      case (2)
                         x_up =vet(1)

                      case (3)
                         x_step =vet(1)

                      case (4)
                         icond = ivet(1)
                   end select
                end do

                call fill_refcodes(key,dire,na,nb,x_low,x_up,x_step,icond,molec,spg)
                if (err_refcodes) return

                n_ini=minloc(ilabel,dim=1)
                ilabel(n_ini)=2
                n_end=minloc(ilabel,dim=1)-1

             end do
          end if

       end if

       return
    End Subroutine Get_RefCodes_Line_Molec

    !!--++
    !!--++ Subroutine Get_RestAng_Line_FAtom(Line, FAtom)
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Get Distance Restraints relations for Free atoms type
    !!--++     Line: Angle [sig] At1a At1b At1c At2a At2b At2c....
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RestAng_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter                :: np=30
       character(len=30),dimension(np)   :: dire
       character(len=8), dimension(2,np) :: symtrans
       integer,dimension(3,np)           :: p
       real                              :: ang,sig

       character(len=8), dimension(2)  :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(3)           :: ivet
       real, dimension(3)              :: vet


       if (len_trim(line) == 0) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Angle ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          err_mess_refcodes="Error in AFIX line: "//trim(line)
          return
       end if
       ang=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.2
          n_ini=2
       else
          sig=max(vet(1),0.001)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,3
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             err_mess_refcodes=" The first atom in AFIX command must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car(1)=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if
          npos=index(dire(i+2),"_")
          if (npos /=0) then
             car(2)=dire(i+2)(npos:)
             dire(i+2)=dire(i+2)(1:npos-1)
          end if

          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (trim(u_case(dire(i+2))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(3)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             err_mess_refcodes="  Some atom names in "//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(:,nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Illegal AFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_ang+1
       n_end=np_rest_ang+nr
       ang_rest(n_ini:n_end)%aobs=ang
       ang_rest(n_ini:n_end)%acalc=0.0
       ang_rest(n_ini:n_end)%sigma=sig
       ang_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       ang_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       ang_rest(n_ini:n_end)%p(3) = p(3,1:nr)
       ang_rest(n_ini:n_end)%STCode(1)=symtrans(1,1:nr)
       ang_rest(n_ini:n_end)%STCode(2)=symtrans(2,1:nr)
       np_rest_ang=n_end

       return
    End Subroutine Get_RestAng_Line


    !!--++
    !!--++ Subroutine Get_RestDis_Line_FAtom(Line, FAtom)
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Get Distance Restraints relations for Free atoms type
    !!--++     Line: Dist [sig] At1a At1b At2a At2b......
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RestDis_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter              :: np=20
       character(len=30),dimension(np) :: dire
       character(len=8), dimension(np) :: symtrans
       integer,dimension(2,np)         :: p
       real                            :: dis,sig

       character(len=8)                :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(2)           :: ivet
       real, dimension(2)              :: vet


       if (len_trim(line) == 0) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Dist ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          err_mess_refcodes="Error in DFIX line: "//trim(line)
          return
       end if
       dis=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.02
          n_ini=2
       else
          sig=max(vet(1),0.0001)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,2
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             err_mess_refcodes=" The first atom in DFIX command must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if

          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             err_mess_refcodes="  Some atom names in"//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Illegal DFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_dis+1
       n_end=np_rest_dis+nr
       dis_rest(n_ini:n_end)%dobs=dis
       dis_rest(n_ini:n_end)%dcalc=0.0
       dis_rest(n_ini:n_end)%sigma=sig
       dis_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       dis_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       dis_rest(n_ini:n_end)%STCode=symtrans(1:nr)
       np_rest_dis=n_end

       return
    End Subroutine Get_RestDis_Line

    !!--++
    !!--++ Subroutine Get_RestTor_Line_FAtom(Line, FAtom)
    !!--++    character(len=*),        intent(in)     :: Line
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++
    !!--++ Overloaded
    !!--++ Get Torsion Restraints relations for Free atoms type
    !!--++     Line: Torsion_Angle [sig] At1a At1b At1c At1d ...
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Get_RestTor_Line(Line, FAtom)
       !---- Arguments ----!
       character(len=*),        intent(in) :: Line
       type(Atom_List_Type),    intent(in) :: FAtom

       !---- Local variables ----!
       integer, parameter                :: np=30
       character(len=30),dimension(np)   :: dire
       character(len=8), dimension(3,np) :: symtrans
       integer,dimension(4,np)           :: p
       real                              :: tor,sig

       character(len=8), dimension(3)  :: car
       integer                         :: i,j,iv,nc,nr,n_ini,n_end,npos
       integer, dimension(4)           :: ivet
       real, dimension(4)              :: vet


       if (len_trim(line) == 0) return

       !---- Description for each word ----!
       call getword(line,dire,nc)

       !---- Get Angle ----!
       call getnum(dire(1),vet,ivet,iv)
       if (iv /= 1) then
          err_refcodes=.true.
          err_mess_refcodes="Error in TFIX line: "//trim(line)
          return
       end if
       tor=vet(1)

       !---- Get Sigma ----!
       call getnum(dire(2),vet,ivet,iv)
       if (iv /= 1) then
          sig=0.5
          n_ini=2
       else
          sig=max(vet(1),0.02)
          n_ini=3
       end if

       nr=0
       symtrans=" "
       do i=n_ini,nc,4
          ivet=0
          car=" "
          npos=index(dire(i),"_")
          if (npos /=0) then
             err_refcodes=.true.
             err_mess_refcodes=" The first atom in TFIX must belong to the asymmetric unit: "//trim(Line)
             return
          end if
          npos=index(dire(i+1),"_")
          if (npos /=0) then
             car(1)=dire(i+1)(npos:)
             dire(i+1)=dire(i+1)(1:npos-1)
          end if
          npos=index(dire(i+2),"_")
          if (npos /=0) then
             car(2)=dire(i+2)(npos:)
             dire(i+2)=dire(i+2)(1:npos-1)
          end if
          npos=index(dire(i+3),"_")
          if (npos /=0) then
             car(3)=dire(i+3)(npos:)
             dire(i+3)=dire(i+3)(1:npos-1)
          end if
          do j=1,FAtom%natoms
             if (trim(u_case(dire(i))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(1)=j
             end if
             if (trim(u_case(dire(i+1))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(2)=j
             end if
             if (trim(u_case(dire(i+2))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(3)=j
             end if
             if (trim(u_case(dire(i+3))) == trim(u_case(FAtom%Atom(j)%Lab))) then
                ivet(4)=j
             end if
             if (all(ivet > 0) ) exit
          end do
          if (any(ivet == 0)) then
             err_refcodes=.true.
             err_mess_refcodes="  Some atom names in"//trim(line)//" not found in the asymmetric unit"
             return
          end if

          !---- New Relation ----!
          nr=nr+1
          p(:,nr)=ivet
          symtrans(:,nr)=car
       end do
       if (nr <= 0) then
          err_refcodes=.true.
          err_mess_refcodes="Illegal TFIX command  "//trim(line)
          return
       end if

       !---- Adding relations ----!
       n_ini=np_rest_tor+1
       n_end=np_rest_tor+nr
       tor_rest(n_ini:n_end)%tobs=tor
       tor_rest(n_ini:n_end)%tcalc=0.0
       tor_rest(n_ini:n_end)%sigma=sig
       tor_rest(n_ini:n_end)%p(1) = p(1,1:nr)
       tor_rest(n_ini:n_end)%p(2) = p(2,1:nr)
       tor_rest(n_ini:n_end)%p(3) = p(3,1:nr)
       tor_rest(n_ini:n_end)%p(4) = p(4,1:nr)
       tor_rest(n_ini:n_end)%STCode(1)=symtrans(1,1:nr)
       tor_rest(n_ini:n_end)%STCode(2)=symtrans(2,1:nr)
       tor_rest(n_ini:n_end)%STCode(3)=symtrans(3,1:nr)
       np_rest_tor=n_end

       return
    End Subroutine Get_RestTor_Line

    !!----
    !!---- Subroutine Init_Err_RefCode()
    !!----
    !!----    Initialize the errors flags in RefCodes
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_RefCodes()

       err_refcodes=.false.
       err_mess_refcodes=" "

       return
    End Subroutine Init_Err_RefCodes

    !!--++
    !!--++ Subroutine Init_RefCodes_FAtom(FAtom)
    !!--++    type(Atom_List_Type), intent(in out) :: FAtom  ! Free Atom Object
    !!--++
    !!--++ Overloaded
    !!--++ Initialize all refinement codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Init_RefCodes_FAtom(FAtom)
       !---- Arguments ----!
       type(Atom_List_Type), intent(in out)  :: FAtom

       !---- Local variables ----!
       integer :: i

       NP_Refi=0
       NP_Cons=0

       V_Vec   =0.0
       V_Name=" "
       V_Bounds=0.0
       V_BCon  =0

       do i=1,FAtom%Natoms
          FAtom%atom(i)%mx=0.0
          FAtom%atom(i)%lx=0

          FAtom%atom(i)%mbiso=0.0
          FAtom%atom(i)%lbiso=0

          FAtom%atom(i)%mocc=0.0
          FAtom%atom(i)%locc=0

          FAtom%atom(i)%mu=0.0
          FAtom%atom(i)%lu=0
       end do

       return
    End Subroutine Init_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Init_RefCodes_MolCrys(MolCrys)
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys    ! Molecular Crystal Object
    !!--++
    !!--++ Overloaded
    !!--++ Initialize all refinement codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Init_RefCodes_MolCrys(MolCrys)
       !---- Arguments ----!
       type(molecular_Crystal_type), intent(in out) :: MolCrys    ! Molecular Crystal Object

       !---- Local variables ----!
       integer :: i

       NP_Refi=0
       NP_Cons=0

       V_Vec   =0.0
       V_Name=" "
       V_Bounds=0.0
       V_BCon  =0

       if (MolCrys%N_Free > 0 .and. allocated(MolCrys%Atm)) then
          do i=1,MolCrys%N_Free
             MolCrys%Atm(i)%mx=0.0
             MolCrys%Atm(i)%lu=0

             MolCrys%Atm(i)%mbiso=0.0
             MolCrys%Atm(i)%lbiso=0

             MolCrys%Atm(i)%mocc=0.0
             MolCrys%Atm(i)%locc=0

             MolCrys%Atm(i)%mu=0.0
             MolCrys%Atm(i)%lu=0
          end do
       end if

       if (MolCrys%N_Mol > 0 .and. allocated(MolCrys%Mol)) then
          do i=1,MolCrys%N_Mol
             MolCrys%mol(i)%mxcentre=0.0
             MolCrys%mol(i)%lxcentre=0

             MolCrys%mol(i)%morient=0.0
             MolCrys%mol(i)%lorient=0

             MolCrys%mol(i)%mT_TLS=0.0
             MolCrys%mol(i)%lT_TLS=0

             MolCrys%mol(i)%mL_TLS=0.0
             MolCrys%mol(i)%lL_TLS=0

             MolCrys%mol(i)%mS_TLS=0.0
             MolCrys%mol(i)%lS_TLS=0

             if (MolCrys%Mol(i)%natoms > 0) then
                MolCrys%mol(i)%mI_coor=0.0
                MolCrys%mol(i)%lI_coor=0

                MolCrys%mol(i)%mbiso=0.0
                MolCrys%mol(i)%lbiso=0

                MolCrys%mol(i)%mocc=0.0
                MolCrys%mol(i)%locc=0
             end if
          end do
       end if

       return
    End Subroutine Init_RefCodes_MolCrys

    !!--++
    !!--++ Subroutine Init_RefCodes_Molec(Molec)
    !!--++    type(molecule_type), intent(in out) :: Molec ! Molecule Object
    !!--++
    !!--++ Overloaded
    !!--++ Initialize all refinement codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Init_RefCodes_Molec(Molec)
       !---- Arguments ----!
       type(molecule_type), intent(in out) :: Molec

       !---- Local variables ----!
       integer :: i

       NP_Refi=0
       NP_Cons=0

       V_Vec   =0.0
       V_Name=" "
       V_Bounds=0.0
       V_BCon  =0

       Molec%mxcentre=0.0
       Molec%lxcentre=0

       Molec%morient=0.0
       Molec%lorient=0

       Molec%mT_TLS=0.0
       Molec%lT_TLS=0

       Molec%mL_TLS=0.0
       Molec%lL_TLS=0

       Molec%mS_TLS=0.0
       Molec%lS_TLS=0

       do i=1,molec%natoms
          Molec%mI_coor(:,i)=0.0
          Molec%lI_coor(:,i)=0

          Molec%mbiso(i)=0.0
          Molec%lbiso(i)=0

          Molec%mocc(i)=0.0
          Molec%locc(i)=0
       end do

       return
    End Subroutine Init_RefCodes_Molec

    !!--++
    !!--++ Subroutine Read_RefCodes_File_FAtom(file_dat,n_ini,n_end,FAtom,Spg)
    !!--++    Type(file_list_type),    intent( in)    :: file_dat
    !!--++    integer,                 intent( in)    :: n_ini
    !!--++    integer,                 intent( in)    :: n_end
    !!--++    type(Atom_List_Type),    intent(in out) :: FAtom
    !!--++    type(space_group_type),  intent(in)     :: Spg
    !!--++
    !!--++ Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_FAtom(file_dat,n_ini,n_end,FAtom,Spg)
       !---- Arguments ----!
       Type(file_list_type),     intent( in)    :: file_dat
       integer,                  intent( in)    :: n_ini
       integer,                  intent( in)    :: n_end
       type(Atom_List_Type),     intent(in out) :: fatom
       type(space_group_type),   intent(in)     :: Spg

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle

          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe")
                call cutst(line,nlong)
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                         return       !!!!
                      case ("orie")
                         key=7
                         return       !!!!
                      case ("ther")
                         key=8
                         return       !!!!
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),fatom,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary")
                call cutst(line,nlong)
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                         return     !!!!
                      case ("orie")
                         key=7
                         return     !!!!
                      case ("ther")
                         key=8
                         return     !!!!
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),fatom,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,fatom)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Angles restraints)

             case ("afix") ! AFIX ang sigma    (Angles restraints)
                call cutst(line,nlong)
                call get_restang_line(line,fatom)

             case ("dequ") ! DEQU sigma        (Distance restraints)

             case ("dfix") ! DFIX d sigma      (Distance restraints)
                call cutst(line,nlong)
                call get_restdis_line(line,fatom)

             case ("flat") ! FLAT

             case ("tequ") ! TEQU sigma        (Torsion angle restraints)

             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)
                call cutst(line,nlong)
                call get_resttor_line(line,fatom)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_FAtom

    !!--++
    !!--++ Subroutine Read_RefCodes_File_Molcrys(file_dat,n_ini,n_end,molcrys)
    !!--++    Type(file_list_type),         intent( in)    :: file_dat
    !!--++    integer,                      intent( in)    :: n_ini
    !!--++    integer,                      intent( in)    :: n_end
    !!--++    type(molecular_crystal_type), intent(in out) :: molcrys
    !!--++
    !!--++ Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_Molcrys(file_dat,n_ini,n_end,molcrys)
       !---- Arguments ----!
       Type(file_list_type),         intent( in)    :: file_dat
       integer,                      intent( in)    :: n_ini
       integer,                      intent( in)    :: n_end
       type(molecular_crystal_type), intent(in out) :: molcrys

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,npos,nlong,iv
       integer                         :: nop_in_line,key,nmol
       integer, dimension(1)           :: ivet
       real, dimension(1)              :: vet

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle

          nmol=-1
          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe") ! FIX
                call cutst(line,nlong)

                !---- Molecule Information ----!
                if (u_case(line(1:3)) =="MOL") then
                   npos=index(line," ")
                   call getnum(line(4:npos),vet,ivet,iv)
                   if (iv /= 1) then
                      nmol=0
                   else
                      nmol=ivet(1)
                   end if
                end if
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),molcrys,nmol)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary") ! VARY
                call cutst(line,nlong)

                !---- Molecule Information ----!
                if (u_case(line(1:3)) =="MOL") then
                   npos=index(line," ")
                   call getnum(line(4:npos),vet,ivet,iv)
                   if (iv /= 1) then
                      nmol=0
                   else
                      nmol=ivet(1)
                   end if
                end if
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),molcrys,nmol)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,molcrys)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Angles restraints)
             case ("afix") ! AFIX ang sigma    (Angles restraints)
             case ("dequ") ! DEQU sigma        (Distance restraints)
             case ("dfix") ! DFIX d sigma      (Distance restraints)
             case ("flat") ! FLAT
             case ("tequ") ! TEQU sigma        (Torsion angle restraints)
             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_Molcrys

    !!--++
    !!--++ Subroutine Read_RefCodes_File_Molec(file_dat,n_ini,n_end,molec,spg)
    !!--++    Type(file_list_type),   intent( in)    :: file_dat
    !!--++    integer,                intent( in)    :: n_ini
    !!--++    integer,                intent( in)    :: n_end
    !!--++    type(molecule_type),    intent(in out) :: molec
    !!--++    type(space_group_type), intent(in)     :: Spg
    !!--++
    !!--++ Subroutine for treatment of Codes controls
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Read_RefCodes_File_Molec(file_dat,n_ini,n_end,molec,spg)
       !---- Arguments ----!
       Type(file_list_type),   intent( in)    :: file_dat
       integer,                intent( in)    :: n_ini
       integer,                intent( in)    :: n_end
       type(molecule_type),    intent(in out) :: molec
       type(space_group_type), intent(in)     :: Spg

       !---- Local variables ----!
       character(len=132)              :: line
       character(len=132),dimension(5) :: dire
       integer                         :: i,k,nlong
       integer                         :: nop_in_line,key

       call init_err_refcodes()

       do i=n_ini,n_end
          line=adjustl(file_dat%line(i))
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle

          select case (l_case(line(1:4)))

             !---- Main Directive: FIX ----!
             case ("fix ", "fixe")
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                         return     !!!!
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"fix",dire(k),molec,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: VARY ----!
             case ("vary")
                call cutst(line,nlong)

                !---- Splitting Line ----!
                call split_operations(line,nop_in_line,dire)
                do k=1,nop_in_line
                   select case (l_case(dire(k)(1:4)))
                      case ("xyz ")
                         key=1
                      case ("occ ")
                         key=2
                      case ("biso")
                         key=3
                      case ("bani")
                         key=4
                         return     !!!!
                      case ("all ")
                         key=5
                      case ("cent")
                         key=6
                      case ("orie")
                         key=7
                      case ("ther")
                         key=8
                      case default
                         key=0
                   end select
                   if (key /=0) call cutst(dire(k),nlong)
                   call get_refcodes_line(key,"var",dire(k),molec,spg)
                   if (err_refcodes) return
                end do

             !---- Main Directive: EQUAL ----!
             case ("equa")
                call cutst(line,nlong)
                call get_concodes_line(line,molec)

             !---- Restraints Cases ----!
             case ("aequ") ! AEQU sigma        (Angles restraints)
             case ("afix") ! AFIX ang sigma    (Angles restraints)
             case ("dequ") ! DEQU sigma        (Distance restraints)
             case ("dfix") ! DFIX d sigma      (Distance restraints)
             case ("flat") ! FLAT
             case ("tequ") ! TEQU sigma        (Torsion angle restraints)
             case ("tfix") ! TFIX ang sigma    (Torsion angle restraints)

          end select
       end do

       return
    End Subroutine Read_RefCodes_File_Molec

    !!--++
    !!--++ Subroutine Split_Operations(Line, Ni, S_Lines)
    !!--++    character(len=*),              intent( in) :: line
    !!--++    integer,                       intent(out) :: ni
    !!--++    character(len=*),dimension(:), intent(out) :: s_lines
    !!--++
    !!--++ Private
    !!--++ Split diferent directives according to KEY_CODE Variable
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Split_Operations(Line, Ni, S_Lines)
       !---- Arguments ----!
       character(len=*),              intent( in) :: line
       integer,                       intent(out) :: ni
       character(len=*),dimension(:), intent(out) :: s_lines

       !---- Local variables ----!
       character(len=150)      :: linec
       character(len=10)       :: car
       integer                 :: i,j,npos,nend
       integer,dimension(nkey) :: ip,ipc

       ni=0
       s_lines=" "

       linec=u_case(line)
       nend=len_trim(line)

       !---- Search Subkeys: XYZ,OCC,BIS... ----!
       ip=0
       do i=1,nkey
          npos=index(linec,key_code(i))
          if (npos > 0) then
             car=linec(npos:)
             j=index(car," ")
             if (j > 0) car=car(:j-1)
             j=index(car,"_")
             if (j > 0) cycle
          end if
          ip(i)=npos
       end do

       npos=count(ip > 0)
       if (npos == 0) then
          ni=1
          s_lines(1)=adjustl(line)
       else
          call sort(ip,nkey,ipc)
          do i=1,nkey
             if (ip(ipc(i)) == 0) then
                if (ip(ipc(i+1)) <= 1) cycle
                ni=ni+1
                s_lines(ni)=adjustl(line(1:ip(ipc(i+1))-1))
             else
                ni=ni+1
                if (i < nkey) then
                   s_lines(ni)=adjustl(line(ip(ipc(i)):ip(ipc(i+1))-1))
                else
                   s_lines(ni)=adjustl(line(ip(ipc(i)):))
                end if
             end if
          end do
       end if

       return
    End Subroutine Split_Operations

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_FAtom(FAtom,Mode)
    !!--++    type(Atom_List_Type),       intent(in out) :: FAtom
    !!--++    character(len=*), optional, intent(in)     :: Mode
    !!--++
    !!--++ Overloaded
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine VState_to_AtomsPar_FAtom(FAtom,Mode)
       !---- Arguments ----!
       type(Atom_List_Type),       intent(in out) :: FAtom
       character(len=*), optional, intent(in)     :: Mode

       !---- Local Variables ----!
       integer          :: i,j,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,FAtom%natoms
          !---- XYZ ----!
          do j=1,3
             l=FAtom%atom(i)%lx(j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   FAtom%atom(i)%x(j)=v_vec(l)*FAtom%atom(i)%mx(j)

                case ("s","S") ! Passing Shift
                   FAtom%atom(i)%x(j)=FAtom%atom(i)%x(j)+v_shift(l)*FAtom%atom(i)%mx(j)
             end select
          end do

          !---- BISO ----!
          l=FAtom%atom(i)%lbiso
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
             case ("v","V") ! Passing Value
                FAtom%atom(i)%biso=v_vec(l)*FAtom%atom(i)%mbiso

             case ("s","S") ! Passing Shift
                FAtom%atom(i)%biso=FAtom%atom(i)%biso+v_shift(l)*FAtom%atom(i)%mbiso
          end select

          !---- OCC ----!
          l=FAtom%atom(i)%locc
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
             case ("v","V") ! Passing Value
                FAtom%atom(i)%occ=v_vec(l)*FAtom%atom(i)%mocc

             case ("s","S") ! Passing Shift
                FAtom%atom(i)%occ=FAtom%atom(i)%occ+v_shift(l)*FAtom%atom(i)%mocc
          end select

          !---- BANIS ----!
          do j=1,6
             l=FAtom%atom(i)%lu(i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   FAtom%atom(i)%u(j)=v_vec(l)*FAtom%atom(i)%mu(j)

                case ("s","S") ! Passing Shift
                   FAtom%atom(i)%u(j)=FAtom%atom(i)%u(j)+v_shift(l)*FAtom%atom(i)%mu(j)
             end select
          end do

       end do

       return
    End Subroutine VState_to_AtomsPar_FAtom

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_Molcrys(Molcrys,Mode)
    !!--++    type(molecular_Crystal_type), intent(in out) :: MolCrys
    !!--++    character(len=*), optional,   intent(in)     :: Mode
    !!--++
    !!--++ Overloaded
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine VState_to_AtomsPar_Molcrys(Molcrys,Mode)
       !---- Arguments ----!
       type(molecular_Crystal_type), intent(in out) :: MolCrys
       character(len=*), optional,    intent(in)     :: Mode

       !---- Local variables ----!
       integer          :: i,j,jj,k,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,molcrys%n_free
          !---- XYZ ----!
          do j=1,3
             l=molcrys%atm(i)%lx(j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%x(j)=v_vec(l)*molcrys%atm(i)%mx(j)

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%x(j)=molcrys%atm(i)%x(j)+v_shift(l)*molcrys%atm(i)%mx(j)
             end select
          end do

          !---- BISO ----!
          l=molcrys%atm(i)%lbiso
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
             case ("v","V") ! Passing Value
                molcrys%atm(i)%biso=v_vec(l)*molcrys%atm(i)%mbiso

             case ("s","S") ! Passing Shift
                molcrys%atm(i)%biso=molcrys%atm(i)%biso+v_shift(l)*molcrys%atm(i)%mbiso
          end select

          !---- OCC ----!
          l=molcrys%atm(i)%locc
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
             case ("v","V") ! Passing Value
                molcrys%atm(i)%occ=v_vec(l)*molcrys%atm(i)%mocc

             case ("s","S") ! Passing Shift
                molcrys%atm(i)%occ=molcrys%atm(i)%occ+v_shift(l)*molcrys%atm(i)%mocc
          end select

          !---- BANIS ----!
          do j=1,6
             l=molcrys%atm(i)%lu(i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%atm(i)%u(j)=v_vec(l)*molcrys%atm(i)%mu(j)

                case ("s","S") ! Passing Shift
                   molcrys%atm(i)%u(j)=molcrys%atm(i)%u(j)+v_shift(l)*molcrys%atm(i)%mu(j)
             end select
          end do
       end do

       do k=1,molcrys%n_mol
          do i=1,molcrys%mol(k)%natoms
             !---- Coordinates ----!
             do j=1,3
                l=molcrys%mol(k)%lI_Coor(j,i)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   err_mess_refcodes="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(k)%I_Coor(j,i)=v_vec(l)*molcrys%mol(k)%mI_Coor(j,i)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(k)%I_Coor(j,i)=molcrys%mol(k)%I_Coor(j,i)+v_shift(l)*molcrys%mol(k)%mI_Coor(j,i)
                end select
             end do

             !---- Biso ----!
             l=molcrys%mol(k)%lbiso(i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%mol(k)%biso(i)=v_vec(l)*molcrys%mol(k)%mbiso(i)

                case ("s","S") ! Passing Shift
                   molcrys%mol(k)%biso(i)=molcrys%mol(k)%biso(i)+v_shift(l)*molcrys%mol(k)%mbiso(i)
             end select

             !---- Occ ----!
             l=molcrys%mol(k)%locc(i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molcrys%mol(k)%occ(i)=v_vec(l)*molcrys%mol(k)%mocc(i)

                case ("s","S") ! Passing Shift
                   molcrys%mol(k)%occ(i)=molcrys%mol(k)%occ(i)+v_shift(l)*molcrys%mol(k)%mocc(i)
             end select

             !---- Centre ----!
             do j=1,3
                l=molcrys%mol(i)%lxcentre(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   err_mess_refcodes="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%xcentre(j)=v_vec(l)*molcrys%mol(i)%mxcentre(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%xcentre(j)=molcrys%mol(i)%xcentre(j)+v_shift(l)*molcrys%mol(i)%mxcentre(j)
                end select
             end do

             !---- Orient ----!
             do j=1,3
                l=molcrys%mol(i)%lorient(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   err_mess_refcodes="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%orient(j)=v_vec(l)*molcrys%mol(i)%morient(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%orient(j)=molcrys%mol(i)%orient(j)+v_shift(l)*molcrys%mol(i)%morient(j)
                end select
             end do

             !---- T_TLS ----!
             do j=1,6
                l=molcrys%mol(i)%lT_TLS(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   err_mess_refcodes="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%T_TLS(j)=v_vec(l)*molcrys%mol(i)%mT_TLS(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%T_TLS(j)=molcrys%mol(i)%T_TLS(j)+v_shift(l)*molcrys%mol(i)%mT_TLS(j)
                end select
             end do

             !---- L_TLS ----!
             do j=1,6
                l=molcrys%mol(i)%lL_TLS(j)
                if (l == 0) cycle
                if (l > np_refi) then
                   err_refcodes=.true.
                   err_mess_refcodes="Number of Refinable parameters is wrong"
                   return
                end if
                select case (car)
                   case ("v","V") ! Passing Value
                      molcrys%mol(i)%L_TLS(j)=v_vec(l)*molcrys%mol(i)%mL_TLS(j)

                   case ("s","S") ! Passing Shift
                      molcrys%mol(i)%L_TLS(j)=molcrys%mol(i)%L_TLS(j)+v_shift(l)*molcrys%mol(i)%mL_TLS(j)
                end select
             end do

             !---- S_TLS ----!
             do j=1,3
                do jj=1,3
                   l=molcrys%mol(i)%lS_TLS(j,jj)
                   if (l == 0) cycle
                   if (l > np_refi) then
                      err_refcodes=.true.
                      err_mess_refcodes="Number of Refinable parameters is wrong"
                      return
                   end if
                   select case (car)
                      case ("v","V") ! Passing Value
                         molcrys%mol(i)%S_TLS(j,jj)=v_vec(l)*molcrys%mol(i)%mS_TLS(j,jj)

                      case ("s","S") ! Passing Shift
                         molcrys%mol(i)%S_TLS(j,jj)=molcrys%mol(i)%S_TLS(j,jj)+v_shift(l)*molcrys%mol(i)%mS_TLS(j,jj)
                   end select
                end do
             end do
          end do
       end do

       return
    End Subroutine VState_to_AtomsPar_Molcrys

    !!--++
    !!--++ Subroutine VState_to_AtomsPar_Molec(Molec,Mode)
    !!--++    type(molecule_type),          intent(in out) :: Molec
    !!--++    character(len=*), optional,   intent(in)     :: Mode
    !!--++
    !!--++ Overloaded
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine VState_to_AtomsPar_Molec(Molec,Mode)
       !---- Arguments ----!
       type(molecule_type),          intent(in out) :: Molec
       character(len=*), optional,   intent(in)     :: Mode

       !---- Local variables ----!
       integer          :: i,j,l
       character(len=1) :: car

       call init_err_refcodes()

       car="s"
       if (present(mode)) car=adjustl(mode)

       do i=1,molec%natoms
          !---- Coordinates ----!
          do j=1,3
             l=molec%lI_Coor(j,i)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                case ("v","V") ! Passing Value
                   molec%I_Coor(j,i)=v_vec(l)*molec%mI_Coor(j,i)

                case ("s","S") ! Passing Shift
                   molec%I_Coor(j,i)=molec%I_Coor(j,i)+v_shift(l)*molec%mI_Coor(j,i)
             end select
          end do

          !---- Biso ----!
          l=molec%lbiso(i)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%biso(i)=v_vec(l)*molec%mbiso(i)

              case ("s","S") ! Passing Shift
                 molec%biso(i)=molec%biso(i)+v_shift(l)*molec%mbiso(i)
          end select

          !---- Occ ----!
          l=molec%locc(i)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%occ(i)=v_vec(l)*molec%mocc(i)

              case ("s","S") ! Passing Shift
                 molec%occ(i)=molec%occ(i)+v_shift(l)*molec%mocc(i)
          end select
       end do

       !---- Centre ----!
       do j=1,3
          l=molec%lxcentre(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%xcentre(j)=v_vec(l)*molec%mxcentre(j)

              case ("s","S") ! Passing Shift
                 molec%xcentre(j)=molec%xcentre(j)+v_shift(l)*molec%mxcentre(j)
          end select
       end do

       !---- Orient ----!
       do j=1,3
          l=molec%lorient(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%orient(j)=v_vec(l)*molec%morient(j)

              case ("s","S") ! Passing Shift
                 molec%orient(j)=molec%orient(j)+v_shift(l)*molec%morient(j)
          end select
       end do

       !---- T_TLS ----!
       do j=1,6
          l=molec%lT_TLS(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%T_TLS(j)=v_vec(l)*molec%mT_TLS(j)

              case ("s","S") ! Passing Shift
                 molec%T_TLS(j)=molec%T_TLS(j)+v_shift(l)*molec%mT_TLS(j)
          end select
       end do

       !---- L_TLS ----!
       do j=1,6
          l=molec%lL_TLS(j)
          if (l == 0) cycle
          if (l > np_refi) then
             err_refcodes=.true.
             err_mess_refcodes="Number of Refinable parameters is wrong"
             return
          end if
          select case (car)
              case ("v","V") ! Passing Value
                 molec%L_TLS(j)=v_vec(l)*molec%mL_TLS(j)

              case ("s","S") ! Passing Shift
                 molec%L_TLS(j)=molec%L_TLS(j)+v_shift(l)*molec%mL_TLS(j)
          end select
       end do

       !---- S_TLS ----!
       do i=1,3
          do j=1,3
             l=molec%lS_TLS(i,j)
             if (l == 0) cycle
             if (l > np_refi) then
                err_refcodes=.true.
                err_mess_refcodes="Number of Refinable parameters is wrong"
                return
             end if
             select case (car)
                 case ("v","V") ! Passing Value
                    molec%S_TLS(i,j)=v_vec(l)*molec%mS_TLS(i,j)

                 case ("s","S") ! Passing Shift
                    molec%S_TLS(i,j)=molec%S_TLS(i,j)+v_shift(l)*molec%mS_TLS(i,j)
             end select
          end do
       end do

       return
    End Subroutine VState_to_AtomsPar_Molec

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_FAtom(FAtom, Iunit)
    !!--++    type(Atom_List_Type),  intent(in) :: FAtom
    !!--++    integer, optional,     intent(in) :: Iunit
    !!--++
    !!--++ Overloaded
    !!--++ Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_FAtom(FAtom, Spg, Iunit)
       !---- Arguments ----!
       type(Atom_List_Type),   intent(in) :: FAtom
       type(Space_Group_Type), intent(in) :: Spg
       integer, optional,      intent(in) :: Iunit

       !---- Local variables ----!
       character(len=20)              :: car
       character(len=60)              :: fmt1,fmt2,fmt3,fmt4,fmt5
       Character(len=25),dimension(3) :: symcar
       integer           :: i,j,k,n,na,np,lun,p1,p2,p3,p4
       real              :: mu
       real,dimension(3) :: tr

       !---- Format Zone ----!
       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="
          do i=1,FAtom%natoms
             do j=1,3
                if (FAtom%atom(i)%lx(j) /= 0) then
                   na=FAtom%atom(i)%lx(j)
                   mu=FAtom%atom(i)%mx(j)
                   car=trim(code_nam(j))//trim(FAtom%atom(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (FAtom%atom(i)%locc /=0) then
                na=FAtom%atom(i)%locc
                mu=FAtom%atom(i)%mocc
                car=trim(code_nam(5))//trim(FAtom%atom(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (FAtom%atom(i)%lbiso /=0) then
                na=FAtom%atom(i)%lbiso
                mu=FAtom%atom(i)%mbiso
                car=trim(code_nam(4))//trim(FAtom%atom(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             do j=1,6
                if (FAtom%atom(i)%lu(j) /=0) then
                   na=FAtom%atom(i)%lu(j)
                   mu=FAtom%atom(i)%mu(j)
                   car=trim(code_nam(5+j))//trim(FAtom%atom(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,FAtom%natoms
                n=n+count(FAtom%atom(j)%lx ==i)
                n=n+count(FAtom%atom(j)%lu==i)
                if (FAtom%atom(j)%locc==i) n=n+1
                if (FAtom%atom(j)%lbiso==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,FAtom%natoms
                   do k=1,3
                      if (FAtom%atom(j)%lx(k) == i) then
                        car=trim(code_nam(i))//trim(FAtom%atom(j)%lab)
                        if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car), &
                              trim(V_Name(i)),FAtom%atom(j)%mx(k)
                      end if
                   end do

                   if (FAtom%atom(j)%lbiso == i) then
                      car=trim(code_nam(4))//trim(FAtom%atom(j)%lab)
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car), &
                           trim(V_Name(i)),FAtom%atom(j)%mbiso
                   end if

                   if (FAtom%atom(j)%locc == i) then
                      car=trim(code_nam(5))//trim(FAtom%atom(j)%lab)
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car), &
                           trim(V_Name(i)),FAtom%atom(j)%mocc
                   end if

                   do k=1,6
                      if (FAtom%atom(j)%lu(k) == i) then
                         car=trim(code_nam(5+k))//trim(FAtom%atom(j)%lab)
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car), &
                              trim(V_Name(i)),FAtom%atom(j)%mu(k)
                      end if
                   end do
                end do
             end if

          end do
       end if

       if (NP_Rest_Dis > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Distance Restraints relations: ",np_rest_dis
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "   N.Rest  Distance    Sigma    Atom1     Atom2: Symmetry_Op"
          write(unit=lun, fmt="(a)") " ==============================================================="
          fmt3="(i7,tr3,f8.4,tr4,f6.3,t34,a,t43,a)"
          do i=1,np_rest_dis
             p1=dis_rest(i)%p(1)
             p2=dis_rest(i)%p(2)
             n=0
             tr=0.0
             symcar(1)=" "
             if (len_trim(dis_rest(i)%stcode) > 0) then
                call Read_SymTrans_Code(dis_rest(i)%stcode,n,tr)
                tr=tr+Spg%Symop(n)%Tr
                call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(1))
             end if
             symcar(1)=": "//symcar(1)
             write(unit=lun, fmt=fmt3) i,dis_rest(i)%dobs,dis_rest(i)%sigma,trim(FAtom%Atom(p1)%Lab), &
                  trim(FAtom%Atom(p2)%Lab)//trim(symcar(1))
          end do
       end if

       if (NP_Rest_Ang > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Angle Restraints relations: ",np_rest_ang
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") &
          "   N.Rest   Atom1        Atom2: Symmetry_Op      Atom3: Symmetry_Op              Angle   Sigma"
          write(unit=lun, fmt="(a)") &
          " =============================================================================================="
          fmt4="(i7,tr5,a,t22,a,t50,a,t75,f12.4,tr3,f7.4)"
          do i=1,np_rest_ang
             p1=ang_rest(i)%p(1)
             p2=ang_rest(i)%p(2)
             p3=ang_rest(i)%p(3)
             do j=1,2
                n=0
                tr=0.0
                symcar(j)=" "
                if (len_trim(Ang_rest(i)%stcode(j)) > 0) then
                   call Read_SymTrans_Code(Ang_rest(i)%stcode(j),n,tr)
                   tr=tr+Spg%Symop(n)%Tr
                   call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(j))
                end if
                symcar(j)=": "//symcar(j)
             end do
             write(unit=lun, fmt=fmt4) i,FAtom%Atom(p1)%Lab,trim(FAtom%Atom(p2)%Lab)//trim(symcar(1)), &
                  trim(FAtom%Atom(p3)%Lab)//trim(symcar(2)),ang_rest(i)%aobs,ang_rest(i)%sigma
          end do
       end if

       if (NP_Rest_Tor > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Torsion Angle Restraints relations: ",np_rest_tor
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") &
          "   N.Rest   Atom1        Atom2: Symmetry_Op      Atom3: Symmetry_Op      Atom4: Symmetry_Op              Angle   Sigma"
          write(unit=lun, fmt="(a)") &
          " ======================================================================================================================"
          fmt5="(i7,tr5,a,t26,a,t34,a,t52,a,t75,f12.4,tr3,f7.4)"
          do i=1,np_rest_tor
             p1=tor_rest(i)%p(1)
             p2=tor_rest(i)%p(2)
             p3=tor_rest(i)%p(3)
             p4=tor_rest(i)%p(4)
             do j=1,3
                n=0
                tr=0.0
                symcar(j)=" "
                if (len_trim(Tor_rest(i)%stcode(j)) > 0) then
                   call Read_SymTrans_Code(Tor_rest(i)%stcode(j),n,tr)
                   tr=tr+Spg%Symop(n)%Tr
                   call get_symSymb(Spg%Symop(n)%Rot,Tr,symcar(j))
                end if
                symcar(j)=": "//symcar(j)
             end do
             write(unit=lun, fmt=fmt5) i,trim(FAtom%Atom(p1)%Lab),trim(FAtom%Atom(p2)%Lab)//trim(symcar(1)), &
                  trim(FAtom%Atom(p3)%Lab)//trim(symcar(2)),trim(FAtom%Atom(p4)%Lab)//trim(symcar(3)),       &
                  tor_rest(i)%tobs,tor_rest(i)%sigma
          end do
       end if

       return
    End Subroutine Write_Info_RefCodes_FAtom

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_Molcrys(MolCrys, iunit)
    !!--++    type(molecular_crystal_type), intent(in) :: molcrys
    !!--++    integer, optional,            intent(in) :: iunit
    !!--++
    !!--++ Overloaded
    !!--++ Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_Molcrys(MolCrys,iunit)
       !---- Arguments ----!
       type(molecular_crystal_type), intent(in) :: molcrys
       integer, optional,            intent(in) :: iunit

       !---- Local variables ----!
       character(len=20) :: car
       character(len=60) :: fmt1,fmt2
       integer           :: i,j,k,kk,n,na,np, lun
       real              :: mu

       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " =>  Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="

          do i=1,Molcrys%N_Free
             do j=1,3
                if (molcrys%atm(i)%lx(j) /=0) then
                   na=molcrys%atm(i)%lx(j)
                   mu=molcrys%atm(i)%mx(j)
                   car=trim(code_nam(j))//trim(molcrys%atm(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (molcrys%atm(i)%locc /=0) then
                na=molcrys%atm(i)%locc
                mu=molcrys%atm(i)%mocc
                car=trim(code_nam(5))//trim(molcrys%atm(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (molcrys%atm(i)%lbiso /=0) then
                na=molcrys%atm(i)%lbiso
                mu=molcrys%atm(i)%mbiso
                car=trim(code_nam(4))//trim(molcrys%atm(i)%lab)
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             do j=1,6
                if (molcrys%atm(i)%lu(j) /=0) then
                   na=molcrys%atm(i)%lu(j)
                   mu=molcrys%atm(i)%mu(j)
                   car=trim(code_nam(5+j))//trim(molcrys%atm(i)%lab)
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do

          do k=1,Molcrys%N_Mol

             do i=1,Molcrys%Mol(k)%natoms
                do j=1,3
                   if (Molcrys%Mol(k)%lI_coor(j,i) /=0) then
                      na=Molcrys%Mol(k)%lI_coor(j,i)
                      mu=Molcrys%Mol(k)%mI_coor(j,i)
                      car=trim(code_nam(j))//trim(molcrys%mol(k)%AtName(i))
                      write(unit=lun,fmt=fmt1)  &
                           trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                   end if
                end do

                if (Molcrys%Mol(k)%locc(i) /=0) then
                   na=Molcrys%Mol(k)%locc(i)
                   mu=Molcrys%Mol(k)%mocc(i)
                   car=trim(code_nam(5))//trim(molcrys%mol(k)%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if

                if (Molcrys%Mol(k)%lbiso(i) /=0) then
                   na=Molcrys%Mol(k)%lbiso(i)
                   mu=Molcrys%Mol(k)%mbiso(i)
                   car=trim(code_nam(4))//trim(molcrys%mol(k)%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             write(unit=lun, fmt="(a)") " "

             do j=1,3
                if (Molcrys%Mol(k)%lxcentre(j) /=0) then
                   na=Molcrys%Mol(k)%lxcentre(j)
                   mu=Molcrys%Mol(k)%mxcentre(j)
                   car=trim(code_nam(12+j))//"entre"
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,3
                if (Molcrys%Mol(k)%lOrient(j) /=0) then
                   na=Molcrys%Mol(k)%lOrient(j)
                   mu=Molcrys%Mol(k)%mOrient(j)
                   car=trim(code_nam(15+j))//"Orient"
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lT_TLS(j) /=0) then
                   na=Molcrys%Mol(k)%lT_TLS(j)
                   mu=Molcrys%Mol(k)%mT_TLS(j)
                   car=trim(code_nam(19))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do j=1,6
                if (Molcrys%Mol(k)%lL_TLS(j) /=0) then
                   na=Molcrys%Mol(k)%lL_TLS(j)
                   mu=Molcrys%Mol(k)%mL_TLS(j)
                   car=trim(code_nam(10))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             do i=1,3
                do j=1,3
                   if (Molcrys%Mol(k)%lS_TLS(i,j) /=0) then
                      na=Molcrys%Mol(k)%lS_TLS(i,j)
                      mu=Molcrys%Mol(k)%mS_TLS(i,j)
                      car=trim(code_nam(21))
                      write(unit=lun,fmt=fmt1)  &
                           trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                   end if
                end do
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,Molcrys%N_Free
                n=n+count(molcrys%atm(j)%lx ==i)
                n=n+count(molcrys%atm(j)%lu==i)
                if (molcrys%atm(j)%locc==i) n=n+1
                if (molcrys%atm(j)%lbiso==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,Molcrys%N_Free
                   do k=1,3
                      if (molcrys%atm(j)%lx(k) == i) then
                         car=trim(code_nam(k))//trim(molcrys%atm(j)%lab)
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car),trim(V_Name(i)),molcrys%atm(j)%mx(k)
                      end if
                   end do

                   if (molcrys%atm(j)%lbiso == i) then
                      car=trim(code_nam(4))//trim(molcrys%atm(j)%lab)
                      if (trim(car)==trim(V_name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  &
                           np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mbiso
                   end if

                   if (molcrys%atm(j)%locc == i) then
                      car=trim(code_nam(5))//trim(molcrys%atm(j)%lab)
                      if (trim(car)==trim(V_name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  &
                           np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mocc
                   end if

                   do k=1,6
                      if (molcrys%atm(j)%lu(k) == i) then
                         car=trim(code_nam(5+k))//trim(molcrys%atm(j)%lab)
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%atm(j)%mu(k)
                      end if
                   end do
                end do
             end if
          end do

          do i=1,NP_Refi
             do k=1,Molcrys%N_Mol
                n=count(molcrys%mol(k)%lxcentre  == i)
                n=n+count(molcrys%mol(k)%lorient == i)
                n=n+count(molcrys%mol(k)%lT_TLS  == i)
                n=n+count(molcrys%mol(k)%lL_TLS  == i)
                n=n+count(molcrys%mol(k)%lS_TLS  == i)
                if (n > 1) then
                   do j=1,3
                      if (molcrys%mol(k)%lxcentre(j) == i) then
                         car=trim(code_nam(12+j))//"entre"
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%mxcentre(j)
                      end if
                   end do

                   do j=1,3
                      if (molcrys%mol(k)%lOrient(j) == i) then
                         car=trim(code_nam(15+j))//"Orient"
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%morient(j)
                      end if
                   end do

                   do j=1,6
                      if (molcrys%mol(k)%lT_TLS(j) == i) then
                         car=trim(code_nam(19))
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np,trim(car), trim(V_Name(i)),molcrys%mol(k)%mT_TLS(j)
                      end if
                   end do

                   do j=1,6
                      if (molcrys%mol(k)%lL_TLS(j) == i) then
                         car=trim(code_nam(20))
                         if (trim(car)==trim(V_name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  &
                              np, trim(car), trim(V_Name(i)),molcrys%mol(k)%mL_TLS(j)
                      end if
                   end do

                   do j=1,3
                      do kk=1,3
                         if (molcrys%mol(k)%lS_TLS(j,kk) == i) then
                            car=trim(code_nam(21))
                            if (trim(car)==trim(V_name(i))) cycle
                            np=np+1
                            write(unit=lun,fmt=fmt2)  &
                                 np,trim(car),trim(V_Name(i)),molcrys%mol(k)%mS_TLS(j,kk)
                         end if
                      end do
                   end do

                end if
             end do
          end do
       end if

       return
    End Subroutine Write_Info_RefCodes_Molcrys

    !!--++
    !!--++ Subroutine Write_Info_RefCodes_Molec(Molec, iunit)
    !!--++    type(molecule_type), intent(in) :: molec
    !!--++    integer, optional,   intent(in) :: iunit
    !!--++
    !!--++ Overloaded
    !!--++ Write the Information about Refinement Codes
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Write_Info_RefCodes_Molec(Molec,iunit)
       !---- Arguments ----!
       type(molecule_type), intent(in) :: molec
       integer, optional,   intent(in) :: iunit

       !---- Local variables ----!
       character(len=60) :: fmt1,fmt2
       character(len=20) :: car
       integer           :: i,j,k,n,na,np,lun
       real              :: mu

       fmt1="(t5,a,t16,i3,t27,a,t33,4(tr6,f8.4),tr8,i2,tr6,f8.3,i9)"
       fmt2="(t10,i3,t21,a,t31,a,t39,f8.4)"

       lun=6
       if (present(iunit)) lun=iunit

       if (NP_Refi > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Refinable Parameters: ",np_refi
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,a)")"    Name      N.Param   Code-Name       Value        L.Bound       ",&
                                      "U.Bound         Step    Condition   Multiplier    N.Atom"
          write(unit=lun, fmt="(a,a)")" ==================================================================",&
                                      "========================================================="

          do i=1,Molec%natoms
             do j=1,3
                if (molec%lI_coor(j,i) /=0) then
                   na=molec%lI_coor(j,i)
                   mu=molec%mI_coor(j,i)
                   car=trim(code_nam(j))//trim(molec%AtName(i))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do

             if (molec%locc(i) /=0) then
                na=molec%locc(i)
                mu=molec%mocc(i)
                car=trim(code_nam(4))//trim(molec%AtName(i))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if

             if (molec%lbiso(i) /=0) then
                na=molec%lbiso(i)
                mu=molec%mbiso(i)
                car=trim(code_nam(5))//trim(molec%AtName(i))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          write(unit=lun, fmt="(a)") " "

          do j=1,3
             if (molec%lxcentre(j) /=0) then
                na=molec%lxcentre(j)
                mu=molec%mxcentre(j)
                car=trim(code_nam(12+j))//"entre"
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,3
             if (molec%lOrient(j) /=0) then
                na=molec%lOrient(j)
                mu=molec%mOrient(j)
                car=trim(code_nam(15+j))//"Orient"
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,6
             if (molec%lT_TLS(j) /=0) then
                na=molec%lT_TLS(j)
                mu=molec%mT_TLS(i)
                car=trim(code_nam(19))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do j=1,6
             if (molec%lL_TLS(j) /=0) then
                na=molec%lL_TLS(j)
                mu=molec%mL_TLS(i)
                car=trim(code_nam(20))
                write(unit=lun,fmt=fmt1)  &
                     trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
             end if
          end do

          do i=1,3
             do j=1,3
                if (molec%lS_TLS(i,j) /=0) then
                   na=molec%lS_TLS(i,j)
                   mu=molec%mS_TLS(i,j)
                   car=trim(code_nam(21))
                   write(unit=lun,fmt=fmt1)  &
                        trim(car),na,trim(V_Name(na)),V_Vec(na),V_Bounds(:,na),V_BCon(na),mu,V_List(na)
                end if
             end do
          end do
       end if

       if (NP_Cons > 0) then
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a,i5)") " => Constraints relations: ",np_cons
          write(unit=lun, fmt="(a)") " "
          write(unit=lun, fmt="(a)") "       N.Constr     Name      Father    Factor"
          write(unit=lun, fmt="(a)") "    ============================================="

          np=0
          do i=1,NP_Refi
             n=0
             do j=1,Molec%natoms
                n=n+count(molec%lI_coor(:,j) ==i)
                if (molec%locc(j)==i) n=n+1
                if (molec%lbiso(j)==i) n=n+1
             end do
             if ( n > 1) then
                do j=1,Molec%natoms
                   do k=1,3
                      if (molec%lI_coor(k,j) == i) then
                         car=trim(code_nam(k))//trim(molec%AtName(j))
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mI_coor(k,j)
                      end if
                   end do

                   if (molec%lbiso(j) == i) then
                      car=trim(code_nam(5))//trim(molec%AtName(j))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np, trim(car), trim(V_Name(i)),molec%mbiso(j)
                   end if

                   if (molec%locc(j) == i) then
                      car=trim(code_nam(4))//trim(molec%AtName(j))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mocc(j)
                   end if
                end do
             end if
          end do

          do i=1,NP_Refi
             n=count(molec%lxcentre  == i)
             n=n+count(molec%lorient == i)
             n=n+count(molec%lT_TLS  == i)
             n=n+count(molec%lL_TLS  == i)
             n=n+count(molec%lS_TLS  == i)
             if (n > 1) then
                do j=1,3
                   if (molec%lxcentre(j) == i) then
                      car=trim(code_nam(12+j))//"entre"
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np, trim(car),trim(V_Name(i)),molec%mxcentre(j)
                   end if
                end do

                do j=1,3
                   if (molec%lOrient(j) == i) then
                      car=trim(code_nam(15+j))//"Orient"
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2)  np,trim(car),trim(V_Name(i)),molec%morient(j)
                   end if
                end do

                do j=1,6
                   if (molec%lT_TLS(j) == i) then
                      car=trim(code_nam(19))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mT_TLS(j)
                  end if
                end do

                do j=1,6
                   if (molec%lL_TLS(j) == i) then
                      car=trim(code_nam(20))
                      if (trim(car)==trim(V_Name(i))) cycle
                      np=np+1
                      write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mL_TLS(j)
                   end if
                end do

                do j=1,3
                   do k=1,3
                      if (molec%lS_TLS(j,k) == i) then
                         car=trim(code_nam(21))
                         if (trim(car)==trim(V_Name(i))) cycle
                         np=np+1
                         write(unit=lun,fmt=fmt2) np,trim(car), trim(V_Name(i)),molec%mS_TLS(j,k)
                      end if
                   end do
                end do
             end if
          end do

       end if

       return
    End Subroutine Write_Info_RefCodes_Molec


    !!---
    !!---  Subroutine Write_Restraints_ObsCalc(iunit)
    !!---     integer, optional,   intent(in) :: iunit
    !!---
    !!---
    !!---  Write the current values of the "observed" and calculated
    !!---  restraints, as well as the corresponding cost value.
    !!---
    !!---  Update: April - 2005
    !!
    Subroutine Write_Restraints_ObsCalc(A,iunit)
       !---- Arguments ----!
       type(Atom_List_Type),intent(in) :: A
       integer, optional,   intent(in) :: iunit

       !---- Local variables ----!
       character(len=14) :: car1,car2,car3,car4
       integer           :: i,i1,i2,i3,i4,lun
       real              :: disto,distc,ango,angc,sigm, cost,w, delta

       lun=6
       if (present(iunit)) lun=iunit

       if( NP_Rest_Dis > 0) then
          Write(unit=lun,fmt="(/,a)") " ============================================================"
          Write(unit=lun,fmt="(a)")   "   Distance Restraints: Atoms, Dobs, Dcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ============================================================"
          Write(unit=lun,fmt="(a)") " Rest#    Atom1         Atom2              Dobs        Dcalc       Sigma   (Do-Dc)/Sigma"
          cost=0.0
          do i=1,NP_Rest_Dis
             i1=Dis_rest(i)%p(1)
             i2=Dis_rest(i)%p(2)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//dis_rest(i)%stcode
             disto=Dis_rest(i)%dobs
             distc=Dis_rest(i)%dcalc
             delta=disto-distc
             sigm=Dis_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,2a,4f12.5)") i,car1,car2,disto,distc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Distance Restraints Cost = Sum{[(dobs-dcalc)/Sigma]^2} = ",cost
       end if


       if( NP_Rest_Ang > 0) then
          Write(unit=lun,fmt="(/,a)") " ============================================================="
          Write(unit=lun,fmt="(a)")   "   Angle Restraints: Atoms, Angobs, Angcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ============================================================="
          Write(unit=lun,fmt="(a)") &
          " Rest#    Atom1         Atom2          Atom3            Ang_obs    Ang_calc      Sigma   (Ao-Ac)/Sigma"

          cost=0.0
          do i=1,NP_Rest_Ang
             i1=Ang_rest(i)%p(1)
             i2=Ang_rest(i)%p(2)
             i3=Ang_rest(i)%p(3)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//ang_rest(i)%stcode(1)
             car3=trim(A%Atom(i3)%lab)//ang_rest(i)%stcode(2)
             ango=Ang_rest(i)%Aobs
             angc=Ang_rest(i)%Acalc
             delta=ango-angc
             sigm=Ang_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,3a,4f12.5)") i,car1,car2,car3,ango,angc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Angle Restraints Cost = Sum{[(Ang_obs-Ang_calc)/Sigma]^2} = ",cost
       End If



       if( NP_Rest_tor > 0) then
          Write(unit=lun,fmt="(/,a)") " ====================================================================="
          Write(unit=lun,fmt="(a)")   "   Torsion Angle Restraints: Atoms, Angobs, Angcalc, Sigma, delt/Sigma"
          Write(unit=lun,fmt="(a,/)") " ====================================================================="
          Write(unit=lun,fmt="(a)") " Rest#    Atom1         Atom2          Atom3          Atom4            "//&
               "Ang_obs    Ang_calc      Sigma   (Ao-Ac)/Sigma"

          cost=0.0
          do i=1,NP_Rest_tor
             i1=Tor_rest(i)%p(1)
             i2=Tor_rest(i)%p(2)
             i3=Tor_rest(i)%p(3)
             i4=Tor_rest(i)%p(4)
             car1=trim(A%Atom(i1)%lab)
             car2=trim(A%Atom(i2)%lab)//tor_rest(i)%stcode(1)
             car3=trim(A%Atom(i3)%lab)//tor_rest(i)%stcode(2)
             car4=trim(A%Atom(i4)%lab)//tor_rest(i)%stcode(3)
             ango=Ang_rest(i)%Aobs
             angc=Ang_rest(i)%Acalc
             delta=ango-angc
             sigm=Ang_rest(i)%sigma
             w= 1.0/(sigm*sigm)
             cost= cost+delta*delta*w
             Write(unit=lun,fmt="(i6,tr4,4a,4f12.5)") i,car1,car2,car3,car4,ango,angc,sigm,delta/sigm
          end do

          Write(unit=lun,fmt="(/,a,f12.5)") "   Torsion Angle Restraints Cost = Sum{[(Ang_obs-Ang_calc)/Sigma]^2} = ",cost
       End If

       return
    End Subroutine Write_Restraints_ObsCalc

 End Module Refinement_Codes
