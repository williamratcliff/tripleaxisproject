!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CONFIGURATION_CALCULATIONS
!!----   INFO: Subroutines related to calculations of energy or
!!----         configuration properties depending on the crystal structure: BVS, Energy,....
!!----
!!---- HISTORY
!!----    Update: March - 2005
!!----
!!----
!!---- DEPENDENCIES
!!----
!!----
!!---- VARIABLES
!!----    ATOM_CONF_LIST
!!----    BVS_ANIONS_N
!!----    BVS_ANIONS
!!----    BVS_ANIONS_RION
!!----    BVS_PAR_TYPE
!!----    BVS_SPECIES_N
!!----    BVS_TABLE
!!----    ERR_CONF
!!----    ERR_MESS_CONF
!!----    TABLE_B
!!----    TABLE_D0
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_ATOMS_CONF_LIST
!!----       BOND_VALENCE
!!----       CALC_BVS
!!----       COST_BVS
!!----       COST_BVS_COULOMBREP
!!----       DEALLOCATE_ATOMS_CONF_LIST
!!----       DEALLOCATE_BVS_TABLE
!!----       INIT_ERR_CONF
!!----       SET_BVS_TABLE
!!----       SET_TABLE_D0_B
!!----       SPECIES_ON_LIST
!!----
!!
 Module Configuration_Calculations
    !---- Use Files ----!
    Use Math_Gen,                   only: Sp, Sort_Strings
    use String_Utilities,           only: Getword, U_Case,pack_string, get_logunit
    Use Atom_Module,                only: Atom_type, Init_Atom_type
    Use Scattering_Chemical_Tables, only: Get_Ionic_Radius
    use Geom_Calculations,          only: Coord_Info

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: Allocate_Atoms_Conf_List, Calc_BVS, Deallocate_Atoms_Conf_List, &
              Init_Err_Conf, Set_BVS_Table, Set_Table_d0_b, Species_on_List,  &
              Cost_BVS, Cost_BVS_CoulombRep, Deallocate_BVS_Table

    !---- List of public private ----!
    private :: Bond_Valence

    !---- Definitions ----!

    !!----
    !!---- TYPE :: ATOMS_CONF_LIST_TYPE
    !!--..
    !!---- Type, public :: Atoms_Conf_List_Type
    !!----    integer                                     :: natoms    ! Total number of atoms in the list
    !!----    integer                                     :: N_Spec    ! Number of different species in the list
    !!----    integer                                     :: N_Anions  ! Number of anions in the list
    !!----    integer                                     :: N_Cations ! Number of cations in the list
    !!----    real                                        :: Tol       ! Tolerance(%) for sum of radii conditions
    !!----    real                                        :: totatoms  ! Total number of atoms in the unit cell
    !!----    character(len=4), dimension(:), allocatable :: Species   ! Symbol + valence
    !!----    real,             dimension(:), allocatable :: Radius    !ionic/atomic radius of species
    !!----    type(Atom_Type),  dimension(:), allocatable :: atom
    !!---- End Type Atoms_Conf_List_Type
    !!----
    !!---- Update: March - 2005
    !!
    Type, public :: Atoms_Conf_List_Type
       integer                                     :: natoms    ! Total number of atoms in the list
       integer                                     :: N_Spec    ! Number of different species in the list
       integer                                     :: N_Anions  ! Number of anions in the list
       integer                                     :: N_Cations ! Number of cations in the list
       real                                        :: Tol       ! Tolerance(%) for sum of radii conditions
       real                                        :: totatoms  ! Total number of atoms in the unit cell
       character(len=4), dimension(:), allocatable :: Species
       real,             dimension(:), allocatable :: Radius    !ionic/atomic radius of species
       type(Atom_Type),  dimension(:), allocatable :: Atom
    End type Atoms_Conf_List_Type

    !!----
    !!---- BVS_ANIONS_N
    !!----    integer, parameter, public :: bvs_anions_n
    !!----
    !!----    Number of anions known in BVS Table in O"Keefe, Breese, Brown
    !!----
    !!---- Update: March - 2005
    !!
    integer, parameter, public :: bvs_anions_n=14

    !!----
    !!---- BVS_Anions
    !!----    character(len=4), parameter, dimension(bvs_anions_n) :: bvs_anions
    !!----
    !!----    Anions known from O'Keefe, Bresse, Brown
    !!----
    !!---- Update: March - 2005
    !!
    character(len=*), parameter, dimension(bvs_anions_n), public :: bvs_anions = &
                     (/"O-2 ","F-1 ","CL-1","BR-1","I-1 ","S-2 ","SE-2","TE-2",  &
                       "N-3 ","P-3 ","AS-3","H-1 ","O-1 ","SE-1"/)

    !!----
    !!---- BVS_Anions_Rion
    !!----    real(kind=sp), parameter, dimension(bvs_anions_n) :: bvs_anions_rion
    !!----
    !!----    Radii Ionic for Anions
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=sp), parameter, dimension(bvs_anions_n), public :: bvs_anions_rion = &
                      (/1.40,1.19,1.67,1.95,2.16,1.84,1.98,2.21,1.71,2.12,2.22,2.08,1.35,1.80/)

    !!----
    !!---- TYPE :: BVS_PAR_TYPE
    !!--..
    !!---- Type, public :: Bvs_Par_Type
    !!----    character (len=4)               :: Symb      ! Chemical symbol
    !!----    real,   dimension(bvs_anions_n) :: d0        ! D0 Parameter
    !!----    real,   dimension(bvs_anions_n) :: b_par     ! B Parameter
    !!----    integer,dimension(bvs_anions_n) :: reference ! Integer pointing to the reference paper
    !!---- End Type Bvs_Par_Type
    !!----
    !!----    Definition for BVS Parameters
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Bvs_Par_Type
       character (len=4)                     :: Symb
       real(kind=sp),dimension(bvs_anions_n) :: d0
       real(kind=sp),dimension(bvs_anions_n) :: b_par
       integer      ,dimension(bvs_anions_n) :: refnum
    End Type Bvs_Par_Type

    !!----
    !!---- BVS_SPECIES_N
    !!----    integer, parameter, public :: bvs_species_n
    !!----
    !!----    Number of Maximum species on BVS_Table
    !!----
    !!---- Update: March - 2005
    !!
    integer, parameter, public :: bvs_species_n=247

    !!----
    !!---- BVS_TABLE
    !!----    Type(Bvs_Par_Type), allocatable, dimension(:), public :: BVS_Table
    !!----
    !!----    BVS Parameters for calculations
    !!----
    !!---- Update: March - 2005
    !!
    Type(Bvs_Par_Type), allocatable, dimension(:), public :: BVS_Table

    !!----
    !!---- ERR_CONF
    !!----    logical, public  :: err_conf
    !!----
    !!----    Logical Variable taking the value .true. if an error in the module
    !!----    CONFIGURATIONS_CALCULATIONS occurs.
    !!----
    !!---- Update: March - 2005
    !!
    logical, public  :: err_conf

    !!----
    !!---- ERR_MESS_CONF
    !!----    character(len=150), public:: err_mess_conf
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: March - 2005
    !!
    character(len=150), public :: err_mess_conf

    !!----
    !!---- TABLE_B
    !!----    real(kind=sp),dimension(:,:), allocatable, private :: Table_b
    !!----
    !!----    Matrix N_Species x N_Species of B parameters for BVS
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=sp),dimension(:,:), allocatable, private :: Table_b

    !!----
    !!---- TABLE_D0
    !!----    real(kind=sp),dimension(:,:), allocatable, private :: Table_d0
    !!----
    !!----    Matrix N_Species x N_Species of D0 for BVS
    !!----
    !!---- Update: March - 2005
    !!
    real(kind=sp),dimension(:,:), allocatable, private :: Table_d0

    !!----
    !!---- TABLE_Ref
    !!----    integer,dimension(:,:), allocatable, private :: Table_ref
    !!----
    !!----    Matrix N_Species x N_Species with references for BVS parameters
    !!----
    !!---- Update: March - 2005
    !!
    integer,dimension(:,:), allocatable, private :: Table_ref

    !!----
    !!----  Reference list for BVS parameters
    !!----
    character(len=*),dimension(0:32),parameter, private :: references = (/  &
         "Unknown                                                                         ", &
         "Brown and Altermatt, (1985), Acta Cryst. B41, 244-247 (empirical)               ", &
         "Brese and O'Keeffe, (1991), Acta Cryst. B47, 192-197 (extrapolated)             ", &
         "Adams, 2001, Acta Cryst. B57, 278-287 (includes second neighbours)              ", &
         "Hu et al. (1995) Inorg. Chim. Acta, 232, 161-165.                               ", &
         "I.D.Brown Private communication                                                 ", &
         "Brown et al. (1984) Inorg. Chem. 23, 4506-4508                                  ", &
         "Palenik (1997) Inorg. Chem. 36 4888-4890                                        ", &
         "Kanowitz and Palenik (1998) Inorg. Chem. 37 2086-2088                           ", &
         "Wood and Palenik (1998) Inorg. Chem. 37 4149-4151                               ", &
         "Liu and Thorp (1993) Inorg. Chem. 32 4102-4105                                  ", &
         "Palenik (1997) Inorg. Chem. 36 3394-3397                                        ", &
         "Shields, Raithby, Allen and Motherwell (1999) Acta Cryst.B56, 455-465           ", &
         "Chen, Zhou and Hu (2002) Chinese Sci. Bul. 47, 978-980.                         ", &
         "Kihlbourg (1963) Ark. Kemi 21 471; Schroeder 1975 Acta Cryst. B31, 2294         ", &
         "Allmann (1975) Monatshefte Chem. 106, 779                                       ", &
         "Zachariesen (1978) J.Less Common Metals 62, 1                                   ", &
         "Krivovichev and Brown (2001) Z. Krist. 216, 245                                 ", &
         "Burns, Ewing and Hawthorne (1997) Can. Miner. 35,1551-1570                      ", &
         "Garcia-Rodriguez, et al. (2000) Acta Cryst. B56, 565-569                        ", &
         "Mahapatra et al. (1996) J. Amer.Chem. Soc. 118, 11555                           ", &
         "Wood and Palenik (1999) Inorg. Chem. 38, 1031-1034                              ", &
         "Wood and Palenik (1999) Inorg. Chem. 38, 3926-3930                              ", &
         "Wood, Abboud, Palenik and Palenik (2000) Inorg. Chem. 39, 2065-2068             ", &
         "Tytko, Mehnike and Kurad (1999) Structure and Bonding 93, 1-66                  ", &
         "Gundemann, et al.(1999) J. Phys. Chem. A 103, 4752-4754                         ", &
         "Zocchi (2000) Solid State Sci. 2 383-387                                        ", &
         "Jensen, Palenik and Tiekiak (2001) Polyhedron 20, 2137                          ", &
         "Roulhac and Palenik (2002) Inorg. Chem. 42, 118-121                             ", &
         "Holsa et al.(2002) J.Solid State Chem 165, 48-55                                ", &
         "Trzesowska, Kruszynski & Bartezak (2004) Acta Cryst. B60, 174-178               ", &
         "Locock & Burns (2004) Z.Krist. 219, 267-271                                     ", &
         "J.Rodriguez-Carvajal Private communication                                      " /)

 Contains

    !!----
    !!---- Subroutine Allocate_Atoms_Conf_List(N,A)
    !!----    integer, intent(in)                         :: n    !  In -> Atoms in asymmetric unit
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A    !  In -> Objet to be allocated
    !!----
    !!----    Allocation of objet A of type atom_list_Conf. This subroutine
    !!----    should be called before using an object of type atom_list_Conf.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Allocate_Atoms_Conf_List(n,A)
       !---- Arguments ----!
       integer,                     intent(in)       :: N  !N. atoms in asymmetric unit
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be allocated

       !---- Local variables ----!
       integer :: i

       A%natoms   = n
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0
       A%tol      = 20.0 ! 20% tolerance by default
       A%totatoms = 0.0
       if(.not. allocated(A%atom)) allocate (A%atom(n))
       do i=1,n
          call init_atom_type(A%atom(i))
       end do

       return
    End Subroutine Allocate_Atoms_Conf_List

    !!----
    !!---- Subroutine Bond_Valence(d0,b0,d,sd,bv,sbv)
    !!----    real(kind=sp), intent(in)  :: d0    ! BVS parameter
    !!----    real(kind=sp), intent(in)  :: b0    ! BVS parameter
    !!----    real(kind=sp), intent(in)  :: d     ! Distance
    !!----    real(kind=sp), intent(in)  :: sd    ! Sigma(d)
    !!----    real(kind=sp), intent(out) :: bv    ! BVS
    !!----    real(kind=sp), intent(out) :: sbv   ! Sigma(bv)
    !!----
    !!----    Zachariasen exponential expression of Bond Valence
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Bond_Valence(D0,B0,D,Sd,Bv,Sbv)
       !---- Arguments ----!
       real(kind=sp),  intent(in)  :: d0,b0  !Bond-valence parameters
       real(kind=sp),  intent(in)  :: d,sd   !Bond distance and sigma
       real(kind=sp),  intent(out) :: bv,sbv !Bond-valence and sigma

       bv=EXP((d0-d)/b0)
       sbv=bv*sd/b0

       return
    End Subroutine Bond_Valence

    !!----
    !!---- Subroutine Calc_BVS(A, Ipr, N_BVSm, BVS_M, Filecod)
    !!----    type (Atoms_Conf_List_type),              intent(in)   :: A            !  In -> Object of Atoms_Conf_List_type
    !!----    integer,                        optional, intent(in)   :: Ipr
    !!----    integer,                        optional, intent(in)   :: n_bvsm       !  In -> Number of modifications
    !!----    character(len=*), dimension(:), optional, intent(in)   :: BVS_M        ! In -> Text with BVS parameters
    !!----    character(len=*),               optional, intent(in)   :: Filecod
    !!----
    !!----    Subroutine to calculate Bond-Valence sums.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Calc_Distance_Angles_Sigma" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    Control for error is present.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Calc_BVS(A, Ipr, N_Bvsm, Bvs_M, Filecod)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A      !  In -> Object of Atoms_Conf_List_type
       integer,                      optional, intent(in)  :: Ipr
       integer,                      optional, intent(in)  :: N_bvsm
       character(len=*),dimension(:),optional, intent(in)  :: bvs_m
       character(len=*),             optional, intent(in)  :: Filecod

       !---- Local variables ----!
       integer                                :: i,j,ic,k,n1,n2,icm,l,icn,isoc,kk,is0, &
                                                 isstr,isdav,isigt,ibvs,icat,ian,sig1,sig2
       real(kind=sp)                          :: tol,fact,del2,s0,q2,  &
                                                 dd,sigtot,efcn,sums,dav,sdav,q1,d2,  &
                                                 str2,sstr,ric,r_2,del,perc,spred,disp,  &
                                                 str,rg1,dist,gii_a,gii_b,gii_c,&
                                                 d0_n, b_n
       character(len=4)                       :: rnatom
       character(len=10), dimension(5)        :: dire

       call init_err_conf()

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20


       if (present(ipr)) then
          write(unit=ipr,fmt="(/,a)")   "  ------------------------------------------------"
          write(unit=ipr,fmt="(a)")     "  {--- BOND-VALENCE AND POLYHEDRA DISTORTIONS ---}"
          write(unit=ipr,fmt="(a,/)")   "  ------------------------------------------------"
       end if


       if (present(n_bvsm).and. present(ipr)) then
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (read from external conditions)"
       else
          write(unit=ipr,fmt="(a,/,a,/)")  &
          " Bond-Valence parameters (d0,B0) for Zachariasen formula:  s= exp{(d0-d)/B0}", &
                "   (data read from internal table)"
       end if

       !---- Each line: Cation Anion d0 b
       if (present(n_bvsm) .and. present(bvs_m)) then
          do k=1,N_bvsm
                 dire=" "
                 call getword(bvs_m(k),dire,ic)
                 if (ic <= 0) cycle
                 icat=0
                 ian=0
                 if (ic < 3 ) then
                          err_conf=.true.
                          err_mess_conf="Cation-Anion d0 and B0 parameters must be provided"
                          return
                 end if

                 do i=1,A%N_Spec
                          if (u_case(dire(1)(1:4)) /= A%Species(i)) cycle
                          icat=i
                 end do
                 do i=1,A%N_Spec
                          if (u_case(dire(2)(1:4)) /= A%Species(i)) cycle
                          ian=i
                 end do
                 if (icat == 0 .or. ian == 0) then
                          err_conf=.true.
                    err_mess_conf="The given Cation or Anion cannot be found in the atom list"
                    return
                 end if
                 if (icat > ian) then
                          j=icat
                          icat=ian
                          ian=j
                 end if
                 if (icat > A%N_Cations) then
                          err_conf=.true.
                    err_mess_conf="A given cation is not found in the atom list"
                    return
                 end if

                 read(unit=dire(3),fmt=*) d0_n
                 if (ic > 3) read(unit=dire(4),fmt=*) b_n

                 Table_d0(icat,ian)=d0_n
                 Table_b(icat,ian)=b_n
                 Table_ref(icat,ian)=0

                 Table_d0(ian,icat)=d0_n
                 Table_b(ian,icat)=b_n
                 Table_ref(ian,icat)=0

          end do
       end if

       do n1=1,A%N_Cations
          do j=1,A%N_Anions
             n2=A%N_Cations+j
             if (present(ipr)) then
                write(unit=ipr,fmt="(2(a,i3,a,a4),2(a,f5.3),a)") &
                      "   Type",n1,": ",A%Species(n1)," with type",n2,": ",A%Species(n2),&
                      " d0=",Table_d0(n1,n2),"    B0=",Table_b(n1,n2), "   => Reference: "//trim(references(Table_ref(n1,n2)))
                write(unit=ipr,fmt="(2(a,a,a,f5.3,a),/)") &
                      "   Cation (Eff. radius): ",A%Species(n1),"(",A%Radius(n1),")   ", &
                      "   Anion  (Eff. radius): ",A%Species(n2),"(",A%Radius(n2),")"
             end if
          end do
       end do

       del2=0.0
       call Get_LogUnit(ibvs)

       if (present(filecod)) then
          open(unit=ibvs,file=trim(filecod)//"_sum.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations for file: "//trim(filecod)//".cfl"
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       else
          open(unit=ibvs,file="summary.bvs",status="replace",action="write", position="append")
          write(unit=ibvs,fmt="(a)") "  Subroutine Calc_BVS (JRC-LLB, version: March-2005)"
          write(unit=ibvs,fmt="(a)") "  Title: Summary of Bond-Valence calculations "
          write(unit=ibvs,fmt="(a)") "   Atom      Coord  D_aver Sigm   Distort(x10-4)    Valence    BVSum(Sigma)"
       end if

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_b=0.0
       gii_c=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0,q1)
          icn=0
          efcn=0.0
          sums=0.0
          sigtot=0.0
          dav=0.0
          sdav=0.0
          str2=0.0
          d2=0.0
          isoc=INT(A%atom(i)%VarF(2)*1000.0+0.5)
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,/,a,/,a,a4,a,f5.3,a,i3,a,/,a,/)") &
                  "    ------------------------------------------------------------------",  &
                  " => Bond-valence and coordination of atom: ",A%atom(i)%lab ," occupancy: ",A%atom(i)%VarF(1),"(",isoc,")",  &
                  "    ------------------------------------------------------------------"
          end if
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(i,j))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(i,j))%charge
             sig2= SIGN(1.0,q2)
             if (sig1 == sig2) cycle
             dd=coord_info%dist(i,j)
             if (dd > (A%Radius(l)+A%Radius(k))*(1.0+tol)) cycle
             icn=icn+1
             efcn=efcn+A%Atom(coord_info%n_cooatm(i,j))%VarF(1)/A%Atom(i)%VarF(1)
             kk=k
             d2=d2+dd*dd
             s0=coord_info%s_dist(i,j)
             dav=dav+dd
             sdav=sdav+s0*s0
             rnatom=A%atom(coord_info%n_cooatm(i,j))%lab
             call Bond_valence(Table_d0(l,k),Table_b(l,k),dd,s0,str,sstr)
             str=str*A%Atom(coord_info%n_cooatm(i,j))%VarF(1)/A%Atom(i)%VarF(1)
             sstr=sstr*A%Atom(coord_info%n_cooatm(i,j))%VarF(1)/A%Atom(i)%VarF(1)
             sums=sums+str
             str2=str2+str*str
             sigtot=sigtot+sstr*sstr
             is0=nint(s0*10000)
             isstr=nint(sstr*1000)
             if (present(ipr)) then
                write(unit=ipr,fmt="(a,a4,a,a4,a,f8.4,a,i4,a,f6.3,a,i3,a)")  &
                     " (",A%atom(i)%lab ,")-(",rnatom,") :",dd,"(",is0,")  ",str,"(",isstr,")"
             end if
          end do

          ric=real(icn)
          if (icn == 0) then
            if (present(ipr) ) then
                write(unit=ipr,fmt=*) " => Warning!! atom: ",A%atom(i)%lab ," is non-coordinated"
                write(unit=ipr,fmt=*) "    Increase the tolerance factor for ionic radii"
             end if
             cycle
          end if
          d2=d2/ric
          sigtot=SQRT(sigtot)
          dav=dav/ric
          sdav=SQRT(sdav)/ric

          isdav=INT(sdav*10000+0.5)
          isigt=INT(sigtot*1000+0.5)
          dist=10000.0*(d2/(dav*dav)-1.0)
          r_2=sums/ric
          r_2=SQRT(ABS(str2/ric-r_2*r_2))
          del=sums-ABS(q1)
          del2=del2+del*del
          perc=100.0*ABS(del/q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_b=gii_b+perc*fact
          gii_c=gii_c+del*del*fact

          !efcn=efcn/A%Atom(i)%VarF(1)                    !Division by the site occupancy
          spred=ABS(q1)/efcn                              !Predicted valence of a single bond
          disp=table_d0(l,kk)-table_b(l,kk)*LOG(spred)    !Pred. distance
          if (present(ipr)) then
             write(unit=ipr,fmt="(/,a,i5,a,f5.2,a,a4)") " Coordination number: ",icn, &
                  "      Eff.Coor. number: ",efcn,"  for atom: ",A%atom(i)%lab
             write(unit=ipr,fmt="(a,f8.4,a,i4,a,f8.3,a)")  &
                  " Average distance  :",dav,"(",isdav,")  Distortion:  ",dist," xE-04"
             write(unit=ipr,fmt="(a,f8.4,a,f6.3)") " Predicted distance:",disp, &
                  "        Single bond-valence S=",spred
             write(unit=ipr,fmt="(a,f8.3,/,a,f8.3,a,i3,a)") &
                  "                                      Valence: ",q1,  &
                  "                                         Sums: ",sums, "(",isigt,")"
             write(unit=ipr,fmt="(a,2f8.3,/,a)") &
                  " Deviation from the Valence Sum Rule (r1,%dev):",del,perc,  &
                  " {r1=Sumj(sij)-Vi, %dev=100abs(r1)/Vi} "
             write(unit=ipr,fmt="(a,f8.3,/,a)")  &
                  " Deviation from the Equal Valence Rule    (r2):",r_2,  &
                  " {r2=<sij-<sij>>rms}"
             write(unit=ibvs,fmt="(tr4,a4,tr4,f6.2,f8.4,a,i4,a,f14.3,2f12.3,a,i3,a)") &
                  A%atom(i)%lab,efcn,dav,"(",isdav,")",dist,q1,sums, "(",isigt,")"

          end if
       end do

       rg1=SQRT(del2/real(A%natoms))*100.0
       gii_a=gii_a*100.0/A%totatoms
       gii_b=gii_b/A%totatoms  !*100.0 already multiplied
       gii_c=sqrt(gii_c/A%totatoms)*100.0

       if (present(ipr)) then
          write(unit=ipr,fmt="(/,6(a,/))")  &
             " => Lines concerning predicted average distances and single",  &
             "    bond-valence values, as well as the deviations from the",  &
             "    Equal Valence Rule,  apply only to those central atoms",  &
             "    having N coordination-atoms of the same chemical species. ",  &
             "    (The term 'single bond-valence' refers to the valence value",  &
             "    of a single bond for a regular polyhedron, so S=Valence/N)"
           write(unit=ipr,fmt="(/,4(a,/))")  &
             " => The Old Global Instability Index (GII) is calculated with the atoms of the asymetric unit (Num_Atoms).",&
             "    The normalized GII(a,b,c) below are calculated using the sum over asymmetric unit but multiplying ",&
             "    differences by the multiplicity of the site. N_Atoms_UCell is the total number of atoms in the ", &
             "    conventional unit cell. In all cases the result of the different expressions is multiplied by 100.0"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
               rg1," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
               gii_a," /100"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
               gii_b," %"
          write(unit=ipr,fmt="(a,f7.2,a)") " =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
               gii_c," /100"
       end if

       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Old Global Instability Index (  GII=SQRT{SUM{|BVS-abs(q)|^2}/Num_Atoms} ) =", &
            rg1," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(a)=       SUM {|BVS-abs(q)|  *mult}       /N_Atoms_UCell =", &
            gii_a," /100"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(b)=       SUM {|BVS-abs(q)|  *mult/abs(q)}/N_Atoms_UCell =", &
            gii_b," %"
       write(unit=ibvs,fmt="(a,f7.2,a)")" =>  Normalized   GII(c)= SQRT{ SUM {|BVS-abs(q)|^2*mult}       /N_Atoms_UCell}=", &
            gii_c," /100"

       call flush(ibvs)
       close (unit=ibvs)

       return
    End Subroutine Calc_BVS

    !!----
    !!---- Subroutine Cost_BVS(A, GII, gic)
    !!----    type (Atoms_Conf_List_type),  intent(in)   :: A    !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=sp),                intent(out)  :: GII  !  OUT -> Global instability index
    !!----    character(len=*),   optional, intent(in)   :: gic  ! If present GII_c is put in GII
    !!----
    !!----    Subroutine to calculate the Global Instability Index.
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Calc_TDist_Coordination" in order
    !!----    to update the internal private variables related to distance/angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in A have to
    !!----    be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Cost_BVS(A, GII,gic)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in)  :: A    !  In -> Object of Atoms_Conf_List_type
       real(kind=sp),                intent(out) :: GII  !  GII_a
       character(len=*),   optional, intent(in)  :: gic  !  If present GII_c is put in GII

       !---- Local variables ----!
       integer       :: i,j,k,icm,l,sig1,sig2
       real(kind=sp) :: tol,fact,q2,dd,sums,q1, del, bv,gii_a,gii_c

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20
       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii_a=0.0
       gii_c=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0,q1)
         sums=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(i,j))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(i,j))%charge
             sig2= SIGN(1.0,q2)
             if (sig1 == sig2) cycle
             dd=coord_info%dist(i,j)
             if (dd > (A%radius(l)+A%radius(k))*(1.0+tol)) cycle
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(i,j))%VarF(1) !Occupacy
             sums=sums+bv
          end do

          del=sums-ABS(q1)

          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii_a=gii_a+abs(del)*fact
          gii_c=gii_c+del*del*fact

       end do
       gii_a=gii_a*100.0/A%totatoms
       gii_c=sqrt(gii_c/A%totatoms)*100.0
       GII=gii_a
       if(present(gic)) GII=gii_c
       return
    End Subroutine Cost_BVS

    !!----
    !!---- Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
    !!----    type (Atoms_Conf_List_type),  intent(in)   :: A     !  In  -> Object of Atoms_Conf_List_type
    !!----    real(kind=sp),                intent(out)  :: GII   !  OUT -> Global instability index
    !!----    real(kind=sp),                intent(out) :: ERep   !  Pseudo Repulsion Coulomb "energy"
    !!----
    !!----
    !!----    Subroutine to calculate the Global Instability Index Gii_a and
    !!----    a pseudo Coulomb repulsion energy useful to avoid cation-cation and
    !!----    anion-anion overlap when using this cost function for predicting or
    !!----    solving a ionic crystal structure. It was used in the old program PiXSA,
    !!----    by J. Pannetier, J. Bassas-Alsina, J.Rodriguez-Carvajal and V. Caignaert,
    !!----    in "Prediction of Crystal Structures from Crystal Chemistry Rules by
    !!----    Simulated Annealing", Nature 346, 343-345 (1990).
    !!---
    !!----    Before calling this subroutine it is the responsibility of the calling
    !!----    program to make a previous call to "Calc_TDist_Coordination" in order
    !!----    to update the internal Coord_Info variables related to distance and
    !!----    angle calculations.
    !!----    Needs as input the object A (of type atom_Conf_list_type, that
    !!----    should be allocated in the calling program).
    !!----    All items corresponding to the bond-valence parameters contained in
    !!----    "A" have to be properly set before calling this procedure.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Cost_BVS_CoulombRep(A, GII, ERep)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),  intent(in)  :: A      !  In -> Object of Atoms_Conf_List_type
       real(kind=sp),                intent(out) :: GII    !  GII_a
       real(kind=sp),                intent(out) :: ERep   !  Pseudo Repulsion Coulomb "energy"

       !---- Local variables ----!
       integer        :: i,j,k,icm,l,sig1,sig2
       real(kind=sp)  :: tol,fact,q2,dd,sums,q1, del, bv

       tol=A%tol*0.01
       if (tol <= 0.001) tol=0.20

       !----
       !---- CAUTION: Subroutine Calc_Dist_Angle_Sigma need
       !----          to be called before using this routine
       gii=0.0
       Erep=0.0
       do i=1,A%natoms
          icm=coord_info%coord_num(i)
          l=A%Atom(i)%ind(1)
          q1=A%Atom(i)%charge
          sig1=SIGN(1.0,q1)
         sums=0.0
          do j=1,icm
             k=A%Atom(coord_info%n_cooatm(i,j))%ind(1)
             q2=A%Atom(coord_info%n_cooatm(i,j))%charge
             sig2= SIGN(1.0,q2)
             dd=coord_info%dist(i,j)
             if (sig1 == sig2) then
                Erep=Erep+ q1*q2/dd
                cycle
             end if
             if (dd > (A%radius(l)+A%radius(k))*(1.0+tol)) cycle
             bv=EXP((Table_d0(l,k)-dd)/Table_b(l,k))
             bv=bv*A%Atom(coord_info%n_cooatm(i,j))%VarF(1)
             sums=sums+bv
          end do

          del=sums-ABS(q1)
          fact=A%Atom(i)%VarF(1)*real(A%atom(i)%mult)
          gii=gii+abs(del)*fact

       end do
       gii=gii*100.0/A%totatoms


       return
    End Subroutine Cost_BVS_CoulombRep


    !!----
    !!---- Subroutine Deallocate_Atoms_Conf_List(A)
    !!----    type (Atoms_Conf_List_Type), intent(in out)   :: A  ! In/ Out -> Objet to be deallocated
    !!----
    !!----    De-allocation of objet A of type Atoms_Conf_List. This subroutine
    !!----    should be after using an object of type Atoms_Conf_List that is no
    !!----    more needed.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_Atoms_Conf_List(A)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out)   :: A  !Objet to be deallocated

       if (allocated(A%atom)) deallocate (A%atom)
       A%natoms   = 0
       A%n_spec   = 0
       A%n_anions = 0
       A%n_cations= 0

       return
    End Subroutine Deallocate_Atoms_Conf_List

    !!----
    !!---- Subroutine Deallocate_BVS_Table()
    !!----
    !!----    Deallocating BVS_Table
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Deallocate_BVS_Table()

       if (allocated(BVS_Table)) deallocate(BVS_Table)

       return
    End Subroutine Deallocate_BVS_Table

    !!----
    !!---- Subroutine Init_Err_Conf()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_Conf()

       err_conf=.false.
       err_mess_conf=" "

       return
    End Subroutine Init_Err_Conf

    !!----
    !!----
    !!---- Subroutine Set_BVS_Table()
    !!----
    !!----    Fills the parameters for BVS from O'Keefe, Bresse, Brown
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Set_BVS_Table()

       if (.not. allocated(BVS_Table)) allocate(BVS_Table(bvs_species_n))

     BVS_Table(  1)=BVS_Par_Type("AC+3", &
                   (/ 2.240, 2.130, 2.630, 2.750, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  2)=BVS_Par_Type("AG+1", &
                   (/ 1.842, 1.800, 2.090, 0.000, 0.000, 2.119, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  3)=BVS_Par_Type("AG+2", &
                   (/ 0.000, 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  4)=BVS_Par_Type("AG+3", &
                   (/ 0.000, 1.830, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  5)=BVS_Par_Type("AG+9", &
                   (/ 0.000, 0.000, 0.000, 2.220, 2.380, 0.000, 2.260, 2.510, 1.850, 2.220, 2.300, 1.500, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     0,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(  6)=BVS_Par_Type("AL+3", &
                   (/ 1.620, 1.545, 2.032, 2.200, 2.410, 2.210, 2.270, 2.480, 1.790, 2.240, 2.300, 1.450, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     1,     1,     2,     2,     5,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(  7)=BVS_Par_Type("AM+3", &
                   (/ 2.110, 2.000, 2.480, 2.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  8)=BVS_Par_Type("AM+4", &
                   (/ 2.080, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(  9)=BVS_Par_Type("AM+5", &
                   (/ 2.070, 1.950, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 10)=BVS_Par_Type("AM+6", &
                   (/ 2.050, 1.950, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 11)=BVS_Par_Type("AS+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.240, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 12)=BVS_Par_Type("AS+3", &
                   (/ 1.789, 1.700, 2.160, 2.350, 2.580, 2.272, 2.400, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     5,     5,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 13)=BVS_Par_Type("AS+5", &
                   (/ 1.767, 1.620, 0.000, 0.000, 0.000, 2.280, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 14)=BVS_Par_Type("AU+1", &
                   (/ 0.000, 0.000, 2.020, 0.000, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 15)=BVS_Par_Type("AU+3", &
                   (/ 1.890, 1.890, 2.170, 2.320, 2.540, 2.390, 0.000, 0.000, 1.940, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 16)=BVS_Par_Type("AU+5", &
                   (/ 0.000, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 17)=BVS_Par_Type("AU+9", &
                   (/ 0.000, 0.000, 0.000, 2.120, 2.340, 2.030, 2.180, 2.410, 1.720, 2.140, 2.220, 1.370, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 18)=BVS_Par_Type("B+3 ", &
                   (/ 1.371, 1.281, 1.740, 1.880, 2.100, 1.770, 1.950, 2.200, 1.470, 1.880, 1.970, 1.140, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     5,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 19)=BVS_Par_Type("BA+2", &
                   (/ 2.285, 2.188, 2.690, 2.880, 3.130, 2.769, 2.880, 3.080, 2.470, 2.880, 2.960, 2.220, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 20)=BVS_Par_Type("BE+2", &
                   (/ 1.381, 1.281, 1.760, 1.900, 2.100, 1.830, 1.970, 2.210, 1.500, 1.950, 2.000, 1.110, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 21)=BVS_Par_Type("BI+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.700, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 22)=BVS_Par_Type("BI+3", &
                   (/ 2.094, 1.990, 2.480, 2.590, 2.820, 2.570, 0.000, 0.000, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 23)=BVS_Par_Type("BI+5", &
                   (/ 2.060, 1.970, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 24)=BVS_Par_Type("BI+9", &
                   (/ 0.000, 0.000, 0.000, 2.620, 2.840, 2.550, 2.720, 2.870, 2.240, 2.630, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 25)=BVS_Par_Type("BK+3", &
                   (/ 2.080, 1.960, 2.350, 2.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 26)=BVS_Par_Type("BK+4", &
                   (/ 2.070, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 27)=BVS_Par_Type("BR+3", &
                   (/ 1.900, 1.750, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 28)=BVS_Par_Type("BR+5", &
                   (/ 1.840, 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 29)=BVS_Par_Type("BR+7", &
                   (/ 1.810, 1.720, 2.190, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 30)=BVS_Par_Type("C+2 ", &
                   (/ 1.366, 0.000, 1.410, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 31)=BVS_Par_Type("C+4 ", &
                   (/ 1.390, 1.320, 1.760, 1.910, 0.000, 1.800, 0.000, 0.000, 1.442, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table( 32)=BVS_Par_Type("C+9 ", &
                   (/ 0.000, 0.000, 0.000, 1.900, 2.120, 1.820, 1.970, 2.210, 1.470, 1.890, 1.990, 1.100, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 33)=BVS_Par_Type("CA+2", &
                   (/ 1.967, 1.842, 2.370, 2.507, 2.720, 2.450, 2.560, 2.760, 2.140, 2.550, 2.620, 1.830, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 34)=BVS_Par_Type("CD+2", &
                   (/ 1.904, 1.811, 2.212, 2.350, 2.570, 2.304, 2.400, 2.590, 1.960, 2.340, 2.430, 1.660, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 35)=BVS_Par_Type("CE+3", &
                   (/ 2.151, 2.036, 2.520, 2.650, 2.870, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 36)=BVS_Par_Type("CE+4", &
                   (/ 2.028, 1.995, 0.000, 0.000, 0.000, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 37)=BVS_Par_Type("CE+9", &
                   (/ 0.000, 0.000, 2.410, 2.690, 2.920, 2.620, 2.740, 2.920, 2.340, 2.700, 2.780, 2.040, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 38)=BVS_Par_Type("CF+3", &
                   (/ 2.070, 1.950, 2.450, 2.550, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 39)=BVS_Par_Type("CF+4", &
                   (/ 2.060, 1.920, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 40)=BVS_Par_Type("CL+3", &
                   (/ 1.710, 1.690, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 41)=BVS_Par_Type("CL+5", &
                   (/ 1.670, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 42)=BVS_Par_Type("CL+7", &
                   (/ 1.632, 1.550, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 43)=BVS_Par_Type("CF+3", &
                   (/ 0.000, 0.000, 2.450, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 44)=BVS_Par_Type("CM+3", &
                   (/ 2.230, 2.120, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 45)=BVS_Par_Type("CM+4", &
                   (/ 2.080, 1.940, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 46)=BVS_Par_Type("CO+1", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5,     0,     0 /) )
     BVS_Table( 47)=BVS_Par_Type("CO+2", &
                   (/ 1.692, 1.640, 2.033, 0.000, 0.000, 1.940, 0.000, 0.000, 1.650, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 48)=BVS_Par_Type("CO+3", &
                   (/ 1.637, 1.620, 2.050, 0.000, 0.000, 2.020, 0.000, 0.000, 1.750, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     9,     2,     2,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 49)=BVS_Par_Type("CO+4", &
                   (/ 1.720, 1.550, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 50)=BVS_Par_Type("CO+9", &
                   (/ 1.655, 0.000, 0.000, 2.180, 2.370, 2.060, 2.240, 2.460, 1.840, 2.210, 2.280, 1.440, 0.000, 0.000 /), &
                   (/ 0.420, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 51)=BVS_Par_Type("CR+2", &
                   (/ 1.730, 1.670, 2.090, 2.260, 2.480, 0.000, 0.000, 0.000, 1.830, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 52)=BVS_Par_Type("CR+3", &
                   (/ 1.724, 1.657, 2.080, 2.280, 0.000, 2.162, 0.000, 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 53)=BVS_Par_Type("CR+4", &
                   (/ 1.810, 1.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 54)=BVS_Par_Type("CR+5", &
                   (/ 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    23,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 55)=BVS_Par_Type("CR+6", &
                   (/ 1.794, 1.740, 2.120, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 56)=BVS_Par_Type("CR+9", &
                   (/ 1.790, 0.000, 0.000, 2.260, 2.450, 2.180, 2.290, 2.520, 1.850, 2.270, 2.340, 1.520, 0.000, 0.000 /), &
                   (/ 0.340, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 57)=BVS_Par_Type("CS+1", &
                   (/ 2.417, 2.330, 2.791, 2.950, 3.180, 2.890, 2.980, 3.160, 2.830, 2.930, 3.040, 2.440, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table( 58)=BVS_Par_Type("CU+1", &
                   (/ 1.610, 1.600, 1.858, 2.030, 2.108, 1.898, 1.900, 0.000, 1.520, 1.774, 1.856, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,    12,     5,     1,     1,    12,     0,    12,    12,    12,     0,     0,     0 /) )
     BVS_Table( 59)=BVS_Par_Type("CU+2", &
                   (/ 1.679, 1.594, 2.000, 1.990, 2.160, 2.054, 2.020, 2.270, 1.751, 1.970, 2.080, 1.210, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,    10,     2,     2,     2,     0,     0 /) )
     BVS_Table( 60)=BVS_Par_Type("CU+3", &
                   (/ 1.735, 1.580, 2.078, 0.000, 0.000, 0.000, 0.000, 0.000, 1.768, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    20,     5,    12,     0,     0,     0,     0,     0,    12,     0,     0,     0,     0,     0 /) )
     BVS_Table( 61)=BVS_Par_Type("DY+2", &
                   (/ 1.900, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 62)=BVS_Par_Type("DY+3", &
                   (/ 2.001, 1.922, 2.410, 2.530, 2.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 63)=BVS_Par_Type("DY+9", &
                   (/ 0.000, 0.000, 0.000, 2.560, 2.770, 2.470, 2.610, 2.800, 2.180, 2.570, 2.640, 1.890, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 64)=BVS_Par_Type("ER+2", &
                   (/ 1.880, 0.000, 0.000, 0.000, 0.000, 2.520, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 65)=BVS_Par_Type("ER+3", &
                   (/ 1.988, 1.904, 2.390, 2.510, 2.750, 2.520, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,    16,    16,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 66)=BVS_Par_Type("ER+9", &
                   (/ 0.000, 0.000, 0.000, 2.540, 2.750, 2.460, 2.590, 2.780, 2.160, 2.550, 2.630, 1.860, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 67)=BVS_Par_Type("ES+3", &
                   (/ 2.080, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 68)=BVS_Par_Type("EU+2", &
                   (/ 2.130, 2.040, 2.530, 2.670, 2.900, 2.584, 0.000, 0.000, 2.340, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     5,     5,     1,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table( 69)=BVS_Par_Type("EU+3", &
                   (/ 2.074, 1.961, 2.480, 2.570, 2.790, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     5,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 70)=BVS_Par_Type("EU+9", &
                   (/ 0.000, 0.000, 0.000, 2.610, 2.830, 2.530, 2.660, 2.850, 2.240, 2.620, 2.700, 1.950, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 71)=BVS_Par_Type("FE+2", &
                   (/ 1.734, 1.650, 2.060, 2.210, 2.470, 2.120, 0.000, 0.000, 1.769, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table( 72)=BVS_Par_Type("FE+3", &
                   (/ 1.759, 1.679, 2.090, 0.000, 0.000, 2.149, 0.000, 0.000, 1.815, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     0,     0,     1,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table( 73)=BVS_Par_Type("FE+4", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 74)=BVS_Par_Type("FE+6", &
                   (/ 1.760, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 75)=BVS_Par_Type("FE+9", &
                   (/ 1.740, 0.000, 0.000, 2.260, 2.470, 2.160, 2.280, 2.530, 1.860, 2.270, 2.350, 1.530, 0.000, 0.000 /), &
                   (/ 0.380, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 76)=BVS_Par_Type("GA+1", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.550 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5 /) )
     BVS_Table( 77)=BVS_Par_Type("GA+3", &
                   (/ 1.730, 1.620, 2.070, 2.200, 2.460, 2.163, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 78)=BVS_Par_Type("GA+9", &
                   (/ 0.000, 0.000, 0.000, 2.240, 2.450, 2.170, 2.300, 2.540, 1.840, 2.260, 2.340, 1.510, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 79)=BVS_Par_Type("GD+2", &
                   (/ 2.010, 2.400, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 80)=BVS_Par_Type("GD+3", &
                   (/ 2.065, 1.950, 2.445, 2.560, 2.780, 2.530, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 81)=BVS_Par_Type("GD+9", &
                   (/ 0.000, 0.000, 0.000, 2.600, 2.820, 2.530, 2.650, 2.840, 2.220, 2.610, 2.680, 1.930, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 82)=BVS_Par_Type("GE+4", &
                   (/ 1.748, 1.660, 2.140, 0.000, 0.000, 2.217, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     1,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 83)=BVS_Par_Type("GE+9", &
                   (/ 0.000, 0.000, 0.000, 2.300, 2.500, 2.230, 2.350, 2.560, 1.880, 2.320, 2.430, 1.550, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 84)=BVS_Par_Type("H+1 ", &
                   (/ 0.569, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.940, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 85)=BVS_Par_Type("HF+3", &
                   (/ 0.000, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 86)=BVS_Par_Type("HF+4", &
                   (/ 1.923, 1.850, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 87)=BVS_Par_Type("HF+9", &
                   (/ 0.000, 0.000, 0.000, 2.470, 2.680, 2.390, 2.520, 2.720, 2.090, 2.480, 2.560, 1.780, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 88)=BVS_Par_Type("HG+1", &
                   (/ 1.900, 1.810, 2.280, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 89)=BVS_Par_Type("HG+2", &
                   (/ 1.972, 2.170, 2.280, 2.380, 2.620, 2.308, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     5,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 90)=BVS_Par_Type("HG+9", &
                   (/ 0.000, 0.000, 0.000, 2.400, 2.590, 2.320, 2.470, 2.610, 2.020, 2.420, 2.500, 1.710, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 91)=BVS_Par_Type("HO+3", &
                   (/ 2.025, 1.908, 2.401, 2.520, 2.760, 2.490, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 92)=BVS_Par_Type("HO+9", &
                   (/ 0.000, 0.000, 0.000, 2.550, 2.770, 2.480, 2.610, 2.800, 2.180, 2.560, 2.640, 1.880, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table( 93)=BVS_Par_Type("I+1 ", &
                   (/ 0.000, 2.320, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 94)=BVS_Par_Type("I+3 ", &
                   (/ 2.020, 1.900, 2.390, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 95)=BVS_Par_Type("I+5 ", &
                   (/ 2.003, 1.840, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 96)=BVS_Par_Type("I+7 ", &
                   (/ 1.930, 1.830, 2.310, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 97)=BVS_Par_Type("IN+1", &
                   (/ 0.000, 0.000, 2.560, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 98)=BVS_Par_Type("IN+3", &
                   (/ 1.902, 1.792, 2.280, 2.510, 2.630, 2.370, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table( 99)=BVS_Par_Type("IN+9", &
                   (/ 0.000, 0.000, 0.000, 2.410, 2.630, 2.360, 2.470, 2.690, 2.030, 2.430, 2.510, 1.720, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(100)=BVS_Par_Type("IR+4", &
                   (/ 1.870, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(101)=BVS_Par_Type("IR+5", &
                   (/ 1.916, 1.820, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(102)=BVS_Par_Type("IR+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.380, 2.510, 2.710, 2.060, 2.460, 2.540, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(103)=BVS_Par_Type("K+1 ", &
                   (/ 2.132, 1.992, 2.519, 2.660, 2.880, 2.590, 2.720, 2.930, 2.260, 2.640, 2.830, 2.100, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(104)=BVS_Par_Type("KR+2", &
                   (/ 0.000, 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(105)=BVS_Par_Type("LA+3", &
                   (/ 2.172, 2.020, 2.545, 2.720, 2.930, 2.643, 2.740, 2.940, 2.340, 2.730, 2.800, 2.060, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,    16,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(106)=BVS_Par_Type("LI+1", &
                   (/ 1.466, 1.360, 1.910, 2.020, 2.220, 1.940, 2.090, 2.300, 1.610, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     0,     0,     0,     0,     0 /) )
     BVS_Table(107)=BVS_Par_Type("LU+3", &
                   (/ 1.971, 1.876, 2.361, 2.500, 2.730, 2.430, 2.560, 2.750, 2.110, 2.510, 2.590, 1.820, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(108)=BVS_Par_Type("MG+2", &
                   (/ 1.693, 1.578, 2.080, 2.280, 2.460, 2.180, 2.320, 2.530, 1.850, 2.290, 2.380, 1.530, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(109)=BVS_Par_Type("MN+2", &
                   (/ 1.790, 1.698, 2.133, 2.340, 0.000, 2.220, 0.000, 0.000, 1.849, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     5,     0,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(110)=BVS_Par_Type("MN+3", &
                   (/ 1.760, 1.660, 2.140, 0.000, 0.000, 0.000, 0.000, 0.000, 1.837, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(111)=BVS_Par_Type("MN+4", &
                   (/ 1.753, 1.710, 2.130, 0.000, 0.000, 0.000, 0.000, 0.000, 1.822, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(112)=BVS_Par_Type("MN+6", &
                   (/ 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(113)=BVS_Par_Type("MN+7", &
                   (/ 1.827, 1.720, 2.170, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(114)=BVS_Par_Type("MN+9", &
                   (/ 1.754, 0.000, 0.000, 2.260, 2.490, 2.200, 0.000, 2.550, 1.870, 2.240, 2.360, 1.550, 0.000, 2.320 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     7,     0,     0,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(115)=BVS_Par_Type("MO+3", &
                   (/ 1.834, 1.760, 2.220, 2.340, 0.000, 0.000, 0.000, 0.000, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    13,     5,     5,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(116)=BVS_Par_Type("MO+4", &
                   (/ 1.886, 1.800, 2.170, 0.000, 0.000, 2.235, 0.000, 0.000, 2.043, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    10,     5,     5,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(117)=BVS_Par_Type("MO+5", &
                   (/ 1.907, 0.000, 2.260, 0.000, 0.000, 2.288, 0.000, 0.000, 2.009, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    10,     0,     5,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(118)=BVS_Par_Type("MO+6", &
                   (/ 1.907, 1.810, 2.280, 0.000, 0.000, 2.331, 0.000, 0.000, 2.009, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(119)=BVS_Par_Type("MO+9", &
                   (/ 1.879, 0.000, 0.000, 2.430, 2.640, 2.350, 2.490, 2.690, 2.040, 2.440, 2.520, 1.730, 0.000, 0.000 /), &
                   (/ 0.300, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    26,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(120)=BVS_Par_Type("N+3 ", &
                   (/ 1.361, 1.370, 1.750, 0.000, 0.000, 1.730, 0.000, 0.000, 1.440, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(121)=BVS_Par_Type("N+5 ", &
                   (/ 1.432, 1.360, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(122)=BVS_Par_Type("NA+1", &
                   (/ 1.803, 1.677, 2.150, 2.330, 2.560, 2.300, 2.410, 2.640, 1.930, 2.360, 2.530, 1.680, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(123)=BVS_Par_Type("NB+3", &
                   (/ 1.910, 1.710, 2.200, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(124)=BVS_Par_Type("NB+4", &
                   (/ 1.880, 1.900, 2.260, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(125)=BVS_Par_Type("NB+5", &
                   (/ 1.911, 1.870, 2.270, 0.000, 2.770, 0.000, 0.000, 0.000, 2.010, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(126)=BVS_Par_Type("NB+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.680, 2.370, 2.510, 2.700, 2.060, 2.460, 2.540, 1.750, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(127)=BVS_Par_Type("ND+2", &
                   (/ 1.950, 0.000, 0.000, 0.000, 0.000, 2.600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(128)=BVS_Par_Type("ND+3", &
                   (/ 2.105, 2.008, 2.492, 2.660, 2.870, 2.590, 2.710, 2.890, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0,     0,     0,     0 /) )
     BVS_Table(129)=BVS_Par_Type("NH+1", &
                   (/ 2.226, 2.129, 2.619, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    19,    19,    19,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(130)=BVS_Par_Type("NI+2", &
                   (/ 1.654, 1.596, 2.020, 2.200, 2.400, 1.980, 0.000, 0.000, 1.700, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(131)=BVS_Par_Type("NI+3", &
                   (/ 1.686, 1.580, 0.000, 0.000, 0.000, 2.040, 0.000, 0.000, 1.731, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    32,     5,     0,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(132)=BVS_Par_Type("NI+4", &
                   (/ 1.780, 1.610, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(133)=BVS_Par_Type("NI+9", &
                   (/ 0.000, 0.000, 0.000, 2.160, 2.340, 2.040, 2.140, 2.430, 1.750, 2.170, 2.240, 1.400, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(134)=BVS_Par_Type("NP+3", &
                   (/ 0.000, 2.000, 2.480, 2.620, 2.850, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(135)=BVS_Par_Type("NP+4", &
                   (/ 2.180, 2.020, 2.460, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(136)=BVS_Par_Type("NP+5", &
                   (/ 2.090, 1.970, 2.420, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(137)=BVS_Par_Type("NP+6", &
                   (/ 2.070, 1.970, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(138)=BVS_Par_Type("NP+7", &
                   (/ 2.060, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(139)=BVS_Par_Type("O+2 ", &
                   (/ 1.500, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(140)=BVS_Par_Type("OS+4", &
                   (/ 1.811, 1.720, 2.190, 2.370, 0.000, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(141)=BVS_Par_Type("OS+5", &
                   (/ 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(142)=BVS_Par_Type("OS+6", &
                   (/ 2.030, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(143)=BVS_Par_Type("OS+8", &
                   (/ 1.920, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(144)=BVS_Par_Type("P+3 ", &
                   (/ 1.630, 1.530, 0.000, 0.000, 0.000, 2.120, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(145)=BVS_Par_Type("P+4 ", &
                   (/ 1.640, 1.660, 0.000, 0.000, 0.000, 2.130, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(146)=BVS_Par_Type("P+5 ", &
                   (/ 1.617, 1.540, 2.020, 2.170, 0.000, 2.145, 0.000, 0.000, 1.704, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     5,     5,     0,     1,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(147)=BVS_Par_Type("P+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.150, 2.400, 2.110, 2.260, 2.440, 1.730, 2.190, 2.250, 1.410, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(148)=BVS_Par_Type("PA+4", &
                   (/ 2.150, 2.020, 2.490, 2.660, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(149)=BVS_Par_Type("PA+5", &
                   (/ 2.090, 2.040, 2.450, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(150)=BVS_Par_Type("PB+2", &
                   (/ 1.963, 2.030, 2.530, 2.680, 2.830, 2.541, 2.690, 0.000, 2.180, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.490, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    17,     2,     2,     5,     5,     1,     5,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(151)=BVS_Par_Type("PB+4", &
                   (/ 2.042, 1.940, 2.430, 3.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(152)=BVS_Par_Type("PB+9", &
                   (/ 0.000, 0.000, 0.000, 2.640, 2.780, 2.550, 2.670, 2.840, 2.220, 2.640, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(153)=BVS_Par_Type("PD+2", &
                   (/ 1.792, 1.740, 2.050, 2.200, 2.360, 2.090, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(154)=BVS_Par_Type("PD+4", &
                   (/ 0.000, 1.660, 0.000, 0.000, 0.000, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(155)=BVS_Par_Type("PD+9", &
                   (/ 0.000, 0.000, 0.000, 2.190, 2.380, 2.100, 2.220, 2.480, 1.810, 2.220, 2.300, 1.470, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(156)=BVS_Par_Type("PM+3", &
                   (/ 0.000, 1.960, 2.450, 2.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(157)=BVS_Par_Type("PO+4", &
                   (/ 2.190, 2.380, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(158)=BVS_Par_Type("PR+3", &
                   (/ 2.138, 2.022, 2.500, 2.670, 2.890, 2.600, 0.000, 2.900, 2.300, 2.680, 2.750, 2.020, 0.000, 2.720 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(159)=BVS_Par_Type("PT+2", &
                   (/ 1.768, 1.680, 2.050, 2.200, 0.000, 2.160, 0.000, 0.000, 1.810, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(160)=BVS_Par_Type("PT+3", &
                   (/ 1.870, 0.000, 2.300, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(161)=BVS_Par_Type("PT+4", &
                   (/ 1.879, 1.759, 2.170, 2.600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(162)=BVS_Par_Type("PT+9", &
                   (/ 0.000, 0.000, 0.000, 2.180, 2.370, 2.080, 2.190, 2.450, 1.770, 2.190, 2.260, 1.400, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(163)=BVS_Par_Type("PU+3", &
                   (/ 2.110, 2.000, 2.480, 2.600, 2.840, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(164)=BVS_Par_Type("PU+4", &
                   (/ 2.090, 1.970, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(165)=BVS_Par_Type("PU+5", &
                   (/ 2.110, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(166)=BVS_Par_Type("PU+6", &
                   (/ 2.060, 1.960, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(167)=BVS_Par_Type("PU+7", &
                   (/ 2.050, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    16,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(168)=BVS_Par_Type("RB+1", &
                   (/ 2.263, 2.160, 2.652, 2.780, 3.010, 2.700, 2.810, 3.000, 2.620, 2.760, 2.870, 2.260, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     1,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table(169)=BVS_Par_Type("RE+1", &
                   (/ 0.000, 0.000, 2.620, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(170)=BVS_Par_Type("RE+3", &
                   (/ 1.900, 0.000, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(171)=BVS_Par_Type("RE+4", &
                   (/ 0.000, 1.810, 2.230, 2.350, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(172)=BVS_Par_Type("RE+5", &
                   (/ 1.860, 0.000, 2.240, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(173)=BVS_Par_Type("RE+6", &
                   (/ 0.000, 1.790, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(174)=BVS_Par_Type("RE+7", &
                   (/ 1.970, 1.860, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(175)=BVS_Par_Type("RE+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.610, 2.370, 2.500, 2.700, 2.060, 2.460, 2.540, 1.750, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(176)=BVS_Par_Type("RH+3", &
                   (/ 1.793, 1.710, 2.080, 2.270, 0.000, 0.000, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(177)=BVS_Par_Type("RH+4", &
                   (/ 0.000, 1.590, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(178)=BVS_Par_Type("RH+5", &
                   (/ 0.000, 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(179)=BVS_Par_Type("RH+9", &
                   (/ 0.000, 0.000, 0.000, 2.250, 2.480, 2.150, 0.000, 2.550, 1.880, 2.290, 2.370, 1.550, 0.000, 2.330 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     0,     2,     2,     2,     2,     2,     0,     2 /) )
     BVS_Table(180)=BVS_Par_Type("RU+2", &
                   (/ 0.000, 1.840, 0.000, 0.000, 0.000, 0.000, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(181)=BVS_Par_Type("RU+3", &
                   (/ 1.770, 2.120, 2.250, 0.000, 0.000, 2.200, 0.000, 0.000, 1.820, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     5,     5,     0,     0,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(182)=BVS_Par_Type("RU+4", &
                   (/ 1.834, 1.740, 2.210, 0.000, 0.000, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(183)=BVS_Par_Type("RU+5", &
                   (/ 1.900, 1.820, 2.230, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(184)=BVS_Par_Type("RU+6", &
                   (/ 1.870, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(185)=BVS_Par_Type("RU+7", &
                   (/ 1.990, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(186)=BVS_Par_Type("RU+9", &
                   (/ 0.000, 0.000, 0.000, 2.260, 2.480, 2.160, 2.330, 2.540, 1.880, 2.290, 2.360, 1.610, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(187)=BVS_Par_Type("S+2 ", &
                   (/ 1.740, 0.000, 0.000, 0.000, 0.000, 2.030, 0.000, 0.000, 1.682, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     5,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(188)=BVS_Par_Type("S+4 ", &
                   (/ 1.644, 1.600, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000, 1.762, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     1,     0,     0,     0,     0,     0 /) )
     BVS_Table(189)=BVS_Par_Type("S+6 ", &
                   (/ 1.624, 1.560, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000, 1.720, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(190)=BVS_Par_Type("S+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.170, 2.360, 2.070, 2.210, 2.450, 1.740, 2.150, 2.250, 1.380, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(191)=BVS_Par_Type("SB+3", &
                   (/ 1.973, 1.883, 2.350, 2.510, 2.760, 2.474, 2.600, 0.000, 2.108, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     5,     1,     5,     0,     4,     0,     0,     0,     0,     0 /) )
     BVS_Table(192)=BVS_Par_Type("SB+5", &
                   (/ 1.942, 1.797, 2.300, 2.480, 0.000, 0.000, 0.000, 0.000, 1.990, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     5,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(193)=BVS_Par_Type("SB+9", &
                   (/ 0.000, 0.000, 0.000, 2.500, 2.720, 2.450, 2.570, 2.780, 2.120, 2.520, 2.600, 2.770, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(194)=BVS_Par_Type("SC+3", &
                   (/ 1.849, 1.760, 2.360, 2.380, 2.590, 2.321, 2.440, 2.640, 1.980, 2.400, 2.480, 1.680, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     5,     2,     2,     1,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(195)=BVS_Par_Type("SE+2", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 2.210, 2.330, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(196)=BVS_Par_Type("SE+4", &
                   (/ 1.811, 1.730, 2.220, 2.430, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(197)=BVS_Par_Type("SE+6", &
                   (/ 1.788, 1.690, 2.160, 0.000, 0.000, 0.000, 0.000, 0.000, 1.900, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(198)=BVS_Par_Type("SE+9", &
                   (/ 0.000, 0.000, 0.000, 2.330, 2.540, 2.250, 2.360, 2.550, 0.000, 2.340, 2.420, 1.540, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     0,     2,     2,     2,     0,     0 /) )
     BVS_Table(199)=BVS_Par_Type("SI+4", &
                   (/ 1.624, 1.580, 2.030, 2.200, 2.410, 2.126, 2.260, 2.490, 1.724, 2.230, 2.310, 1.470, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     1,     2,     2,     1,     2,     2,     2,     0,     0 /) )
     BVS_Table(200)=BVS_Par_Type("SM+3", &
                   (/ 2.088, 1.940, 1.977, 2.660, 2.840, 2.550, 2.670, 2.860, 2.240, 2.630, 2.700, 1.960, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,    16,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(201)=BVS_Par_Type("SN+2", &
                   (/ 1.940, 1.925, 2.410, 2.530, 2.810, 2.440, 0.000, 0.000, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,     5,     4,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(202)=BVS_Par_Type("SN+3", &
                   (/ 0.000, 0.000, 2.360, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(203)=BVS_Par_Type("SN+4", &
                   (/ 1.905, 1.843, 2.276, 2.400, 0.000, 2.399, 2.510, 0.000, 2.030, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     1,     5,     0,     1,     5,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(204)=BVS_Par_Type("SN+9", &
                   (/ 0.000, 0.000, 0.000, 2.550, 2.760, 2.390, 2.590, 2.760, 2.060, 2.450, 2.620, 1.850, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,    27,     2,     2,    27,     2,     2,     2,     0,     0 /) )
     BVS_Table(205)=BVS_Par_Type("SR+2", &
                   (/ 2.118, 2.019, 2.510, 2.680, 2.880, 2.590, 2.720, 2.870, 2.230, 2.670, 2.760, 2.010, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(206)=BVS_Par_Type("TA+4", &
                   (/ 2.290, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(207)=BVS_Par_Type("TA+5", &
                   (/ 1.920, 1.880, 2.300, 0.000, 0.000, 2.470, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(208)=BVS_Par_Type("TA+9", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.390, 2.510, 2.700, 2.010, 2.470, 2.550, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(209)=BVS_Par_Type("TB+3", &
                   (/ 2.032, 1.936, 2.427, 2.580, 2.800, 2.510, 2.630, 2.820, 2.200, 2.590, 2.660, 1.910, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(210)=BVS_Par_Type("TC+4", &
                   (/ 0.000, 1.880, 2.210, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(211)=BVS_Par_Type("TC+7", &
                   (/ 1.900, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(212)=BVS_Par_Type("TE+4", &
                   (/ 1.977, 1.870, 2.370, 2.550, 2.787, 2.440, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(213)=BVS_Par_Type("TE+6", &
                   (/ 1.917, 1.820, 2.300, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(214)=BVS_Par_Type("TE+9", &
                   (/ 0.000, 0.000, 0.000, 2.530, 2.760, 2.450, 2.530, 2.760, 2.120, 2.520, 2.600, 1.830, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(215)=BVS_Par_Type("TH+4", &
                   (/ 2.167, 2.068, 2.550, 2.710, 2.930, 2.640, 2.760, 2.940, 2.340, 2.730, 2.800, 2.070, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(216)=BVS_Par_Type("TI+2", &
                   (/ 0.000, 2.150, 2.310, 2.490, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(217)=BVS_Par_Type("TI+3", &
                   (/ 1.791, 1.723, 2.220, 0.000, 2.520, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     5,     0,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(218)=BVS_Par_Type("TI+4", &
                   (/ 1.815, 1.760, 2.190, 2.360, 0.000, 2.290, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(219)=BVS_Par_Type("TI+9", &
                   (/ 1.790, 0.000, 2.184, 2.320, 2.540, 2.240, 2.380, 2.600, 1.930, 2.360, 2.420, 1.610, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    11,     0,    11,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(220)=BVS_Par_Type("TL+1", &
                   (/ 2.124, 2.150, 2.560, 2.690, 2.822, 2.545, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     1,     1,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(221)=BVS_Par_Type("TL+3", &
                   (/ 2.003, 1.880, 2.320, 2.650, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(222)=BVS_Par_Type("TL+9", &
                   (/ 0.000, 0.000, 0.000, 2.700, 2.910, 2.630, 2.700, 2.930, 2.290, 2.710, 2.790, 2.050, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(223)=BVS_Par_Type("TM+3", &
                   (/ 2.000, 1.842, 2.380, 2.530, 2.740, 2.450, 2.580, 2.770, 2.140, 2.530, 2.620, 1.850, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(224)=BVS_Par_Type("U+2 ", &
                   (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.080, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     5,     0 /) )
     BVS_Table(225)=BVS_Par_Type("U+3 ", &
                   (/ 0.000, 2.020, 2.490, 2.640, 2.870, 2.540, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.400, 0.400, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,    16,    16,    16,    16,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(226)=BVS_Par_Type("U+4 ", &
                   (/ 2.112, 2.038, 2.470, 2.600, 2.880, 2.550, 0.000, 0.000, 2.180, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     1,    16,    16,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(227)=BVS_Par_Type("U+5 ", &
                   (/ 2.075, 1.966, 2.460, 2.700, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     2,     2,     2,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(228)=BVS_Par_Type("U+6 ", &
                   (/ 2.051, 1.980, 2.420, 0.000, 0.000, 0.000, 0.000, 0.000, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.519, 0.400, 0.400, 0.370, 0.370, 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    18,    16,    16,     0,     0,     0,     0,     0,     5,     0,     0,     0,     0,     0 /) )
     BVS_Table(229)=BVS_Par_Type("U+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.630, 2.840, 2.560, 2.700, 2.860, 2.240, 2.640, 2.720, 1.970, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(230)=BVS_Par_Type("V+1 ", &
                   (/ 1.880, 0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(231)=BVS_Par_Type("V+2 ", &
                   (/ 1.700, 2.160, 2.440, 0.000, 0.000, 2.110, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(232)=BVS_Par_Type("V+3 ", &
                   (/ 1.743, 1.702, 2.190, 2.330, 0.000, 2.170, 0.000, 0.000, 1.813, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     5,     0,     5,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(233)=BVS_Par_Type("V+4 ", &
                   (/ 1.784, 1.700, 2.160, 0.000, 0.000, 2.226, 0.000, 0.000, 1.875, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,    10,     0,     0,    10,     0,     0,     0,     0,     0 /) )
     BVS_Table(234)=BVS_Par_Type("V+5 ", &
                   (/ 1.803, 1.700, 2.160, 0.000, 0.000, 2.250, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     5,     2,     0,     0,     5,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(235)=BVS_Par_Type("V+9 ", &
                   (/ 1.810, 0.000, 0.000, 2.300, 2.510, 2.230, 2.330, 2.570, 1.860, 2.310, 2.390, 1.580, 0.000, 0.000 /), &
                   (/ 0.340, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/    15,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(236)=BVS_Par_Type("W+5 ", &
                   (/ 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(237)=BVS_Par_Type("W+6 ", &
                   (/ 1.917, 1.830, 2.270, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(238)=BVS_Par_Type("W+9 ", &
                   (/ 0.000, 0.000, 0.000, 2.450, 2.660, 2.390, 2.510, 2.710, 2.060, 2.460, 2.540, 1.760, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     0,     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(239)=BVS_Par_Type("XE+2", &
                   (/ 2.050, 2.020, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.350, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(240)=BVS_Par_Type("XE+4", &
                   (/ 0.000, 1.930, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     0,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(241)=BVS_Par_Type("XE+6", &
                   (/ 2.000, 1.890, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(242)=BVS_Par_Type("XE+8", &
                   (/ 1.940, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(243)=BVS_Par_Type("Y+3 ", &
                   (/ 2.019, 1.904, 2.400, 2.550, 2.770, 2.480, 2.610, 2.800, 2.170, 2.570, 2.640, 1.860, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(244)=BVS_Par_Type("YB+3", &
                   (/ 1.965, 1.875, 2.371, 2.451, 2.720, 2.430, 2.560, 2.760, 2.120, 2.530, 2.590, 1.820, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     BVS_Table(245)=BVS_Par_Type("ZN+2", &
                   (/ 1.704, 1.620, 2.010, 2.150, 2.360, 2.090, 2.220, 2.450, 1.720, 2.150, 2.240, 1.420, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     2,     2,     2,     2,     2,     2,     2,     5,     2,     2,     2,     0,     0 /) )
     BVS_Table(246)=BVS_Par_Type("ZR+2", &
                   (/ 2.340, 2.240, 2.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     5,     5,     5,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 /) )
     BVS_Table(247)=BVS_Par_Type("ZR+4", &
                   (/ 1.928, 1.846, 2.330, 2.480, 2.690, 2.410, 2.530, 2.670, 2.110, 2.520, 2.570, 1.790, 0.000, 0.000 /), &
                   (/ 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370, 0.370 /), &
                   (/     1,     1,     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,     0,     0 /) )
     return
    End Subroutine Set_BVS_Table

    !!----
    !!---- Subroutine Set_Table_d0_b()
    !!----
    !!----
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Set_Table_d0_b(A)
       !---- Arguments ----!
       type (Atoms_Conf_List_type),            intent(in)  :: A

       !---- Local Variables ----!
       integer :: i,j,k,ia,ic

       if (A%N_Spec == 0) then
          err_conf=.true.
          err_mess_conf=" The number of different species was zero on Table_d0"
          return
       end if

       if (allocated(Table_d0))   deallocate(Table_d0)
       if (allocated(Table_b))    deallocate(Table_b)
       if (allocated(Table_ref))  deallocate(Table_ref)

       allocate(Table_d0(A%N_Spec,A%N_Spec))
       allocate(Table_b(A%N_Spec,A%N_Spec))
       allocate(Table_ref(A%N_Spec,A%N_Spec))

       Table_d0=0.0
       Table_b=0.37
       Table_ref = 0

       call Set_BVS_Table()

       do i=1,A%N_Cations
          ic=0
          do j=1,bvs_species_n
             if (A%Species(i) == BVS_Table(j)%Symb) then
                ic=j
                exit
             end if
          end do
          if (ic == 0) then
             err_conf=.true.
             err_mess_conf=" Cation not found on the internal list: "//A%Species(i)
             return
          end if

          do k=1,A%N_Anions
             ia=0
             do j=1,bvs_anions_n
                if (A%Species(A%N_Cations+k) == bvs_anions(j)) then
                   ia=j
                   exit
                end if
             end do
             if (ia == 0) then
                err_conf=.true.
                err_mess_conf=" Anion not found on the internal list: "//A%Species(A%N_Cations+k)
                return
             end if

             Table_d0 (i,A%N_Cations+k)=bvs_table(ic)%d0(ia)
             Table_b  (i,A%N_Cations+k)=bvs_table(ic)%b_par(ia)
             Table_ref(i,A%N_Cations+k)=bvs_table(ic)%refnum(ia)

             Table_d0 (A%N_Cations+k,i)=bvs_table(ic)%d0(ia)
             Table_b  (A%N_Cations+k,i)=bvs_table(ic)%b_par(ia)
             Table_ref(A%N_Cations+k,i)=bvs_table(ic)%refnum(ia)

          end do
       end do

       call Deallocate_BVS_Table()

       return
    End Subroutine Set_Table_d0_b


    !!----
    !!---- Subroutine Species_on_List(A,MulG,tol)
    !!----    type (Atoms_Conf_List_Type), intent(in out) :: A
    !!----    Integer, optional,           intent(in)     :: MulG
    !!----    real, optional,              intent(in)     :: tol
    !!----
    !!----    Determines the different species in the List and,
    !!----    optionally, sets the tolerance factor for ionic radii
    !!----    conditions and provides "corrected" occupation factors
    !!----    (mult/MulG) when the user is using a multiplier. The
    !!----    general multiplicity of the space group MulG must be
    !!----    provided in such a case. This first free variable of the
    !!----    Atom-type A%Atom%VFree(1) is set to the corrected
    !!----    occupation. The first atom in the list must completely
    !!----    occupy its site.
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Species_on_List(A,MulG, tol)
       !---- Arguments ----!
       type (Atoms_Conf_List_Type), intent(in out) :: A
       Integer, optional,           intent(in)     :: MulG
       real, optional,              intent(in)     :: tol

       !---- Local variables ----!
       character(len=4), dimension(50) :: cation,anion,spec
       character(len=2)                :: car,cv
       character(len=4)                :: canio
       integer                         :: i,im,j,v,ns,nc,na
       real                            :: fac1,fact


       if (A%natoms == 0) return

       ns=0
       spec  = " "
       nc=0
       cation= " "
       na=0
       anion = " "

       if(present(tol)) A%tol=tol

       if(present(MulG)) then

         fac1=A%atom(1)%Occ*real(MulG)/real(A%atom(1)%mult)
         fac1=1.0/fac1
         A%totatoms=0.0
         do i=1,a%natoms
            fact=real(MulG)/real(a%atom(i)%mult)
            A%Atom(i)%VarF(1)=A%atom(i)%occ*fact*fac1  !Site Occupancy (=1, full occupation)
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std*fact*fac1  !standard deviation of Site Occupancy
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       else

         A%totatoms=0.0
         do i=1,a%natoms
            A%Atom(i)%VarF(1)=A%atom(i)%occ  !The user has given site occupancy
            A%Atom(i)%VarF(2)=A%atom(i)%occ_std
            A%totatoms=A%totatoms + A%Atom(i)%VarF(1)*real(a%atom(i)%mult) !total number of atoms/conventional unit cell
         end do

       end if

       loop1:do i=1, A%natoms
          car=u_case(a%atom(i)%ChemSymb)
          v=nint(a%atom(i)%charge)
          if (v == 0) then
             err_conf=.true.
             err_mess_conf=" The Atom "//a%atom(i)%lab//"has not charge"
             return
          end if
          write(unit=cv,fmt="(i2)") v
          if (v > 0) cv(1:1)="+"
          canio=car//cv
          canio=pack_string(canio)

          if (v > 0) then
             do j=1,nc
                if (canio == cation(j)) cycle loop1
             end do
             nc=nc+1
             cation(nc)=canio
          else
             do j=1,na
                if (canio == anion(j)) cycle loop1
             end do
             na=na+1
             anion(na)=canio
          end if
          if (na+nc == 50) exit
       end do loop1

       ns=nc+na
       A%N_Spec    = ns
       A%N_Anions  = na
       A%N_Cations = nc

       !---- Order the Species vector ----!
       call sort_strings(cation(1:nc))
       call sort_strings(anion(1:na))
       spec(1:nc)=cation(1:nc)
       spec(nc+1:nc+na)=anion(1:na)

       if (allocated(A%Species)) deallocate(A%Species)
       allocate(A%Species(ns))
       A%Species=spec(1:ns)

       if (allocated(A%Radius)) deallocate(A%Radius)
       allocate(A%Radius(ns))

       do i=1,nc
          im=index(A%Species(i),"+")
          car=A%Species(i)(im+1:im+2)
          read(unit=car,fmt="(i1)") j
          car=A%Species(i)(1:im-1)
          call get_ionic_radius(car,j,A%Radius(i))
          if (A%Radius(i) < 0.01) A%Radius(i)=0.5
       end do

       do i=1,A%N_Anions
          do j=1,bvs_anions_n
             if (A%Species(nc+i) == bvs_anions(j)) then
                  A%Radius(nc+i) = bvs_anions_rion(j)
                exit
             end if
          end do
       end do

       !---- Fix the index on Atom_type point out Species vector ----!
       do i=1, A%natoms
          do j=1,ns
             im=index(A%Species(j),"+")
             if (im == 0) im=index(A%Species(j),"-")
             car=A%Species(j)(1:im-1)
             cv=A%Species(j)(im:im+1)
             if (cv(1:1)=="+") cv(1:1)=" "
                read(unit=cv,fmt="(i2)") v
                if (u_case(A%Atom(i)%ChemSymb) == car .and. nint(A%Atom(i)%charge) == v) then
                A%atom(i)%ind(1)=j
                exit
             end if
          end do
       end do

       return
    End Subroutine Species_on_List

 End Module Configuration_Calculations
