!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: OPTIMIZATION_SAN
!!----   INFO: Module for Global Optimization using Simulated Annealing.
!!----         Currently there is available only a generic Simulated Anneling subroutine
!!----         That must be called with the name of a user-supplied subroutine to calculate
!!----         the cost function as an argument.
!!----         The calling program must define at least two variables of derived types
!!----         SIMANN_CONDITIONS_TYPE and STATE_VECTOR_TYPE respectively.
!!----         The generic simulated annealing procedure can use the constant step algorithm
!!----         or the Corana algorithm depending on the values of the corresponding component
!!----         of the SIMANN_CONDITIONS_TYPE user-defined variable.
!!----
!!---- DEPENDENCIES
!!----     use IO_Messages, mess => write_scroll_text
!!----     use String_Utilities, only: u_case
!!----     use IO_formats,       only: File_List_Type
!!----
!!---- HISTORY
!!----    Update: March - 2005
!!----
!!----            Created by JRC in October -2002
!!----
!!---- VARIABLES
!!----    NP_SAN
!!----    SIMANN_CONDITIONS_TYPE
!!----    STATE_VECTOR_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       SET_SIMANN_COND
!!----       SET_SIMANN_MSTATEV
!!----       SET_SIMANN_STATEV
!!----       SIMANNEAL_GEN
!!----       SIMANNEAL_MULTICONF
!!--++       INIT_RAN        [Private]
!!--++       CHECKM          [Private]
!!--++       CHECK           [Private]
!!----       WRITE_SIMANN_COND
!!----       WRITE_SIMANN_MSTATEV
!!----       WRITE_SIMANN_STATEV
!!
!!
 Module Optimization_SAN
    !---- Use Files ----!
    use IO_Messages, mess => write_scroll_text
    use String_Utilities, only: u_case
    use IO_formats,       only: File_List_Type

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public :: SimAnneal_gen, SimAnneal_MultiConf, Set_SimAnn_Cond, Set_SimAnn_StateV, &
              Write_SimAnn_Cond, Write_SimAnn_StateV, Set_SimAnn_MStateV, Write_SimAnn_MStateV

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private:: checkm, check, init_ran

    !!----
    !!---- NP_SAN
    !!----    integer, parameter, public :: np_SAN=80
    !!----
    !!----    Maximum number of parameters in the model
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public  :: np_SAN=80

    !!----
    !!---- NP_CONF
    !!----    integer, parameter, public :: np_CONF=30
    !!----
    !!----    Maximum number of initial configurations in paralell
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public  :: np_CONF=30


    !---- Definitions ----!

    !!----
    !!---- Type, public       :: SimAnn_Conditions_type
    !!----  real              :: t_ini       ! Initial temperature
    !!----  real              :: anneal      ! Kirpactrick factor for Annealing
    !!----  real              :: accept      ! Minimum percentage of accepted configurations
    !!----  integer           :: initconfig  ! Flag determining if the first configuration is random or read
    !!----  integer           :: nalgor      ! Flag determining if the Corana algorithm is selected (0) or not (/=0)
    !!----  integer           :: nm_cycl     ! Number of Cycles per temp  in SA searchs
    !!----  integer           :: num_temps   ! Maximum number of temperatures in SA
    !!----  integer           :: num_therm   ! Number of thermalization cycles in SA
    !!----  integer           :: num_conf    ! Number of paralell configurations in SA
    !!----  character(len=60) :: Cost_function_name
    !!----  integer           :: seed=0      ! If different from zero, holds the seed
    !!----                                   ! for random number generator
    !!---- End type Opt_Conditions_Type
    !!----
    !!---- Derived type containing the conditions for running the
    !!---- Simulated Annealing Algorithm
    !!----
    !!---- Update: March - 2005
    !!
    Type, public        :: SimAnn_Conditions_type
      real              :: t_ini       ! Initial temperature
      real              :: anneal      ! Kirpactrick factor for Annealing
      real              :: accept      ! Minimum percentage of accepted configurations
      integer           :: initconfig  ! Flag determining if the first configuration is random or read
      integer           :: nalgor      ! Flag determining if the Corana algorithm is selected (0) or not (/=0)
      integer           :: nm_cycl     ! Number of Cycles per temp  in SA searchs
      integer           :: num_temps   ! Maximum number of temperatures in SA
      integer           :: num_therm   ! Number of thermalization cycles in SA
      integer           :: num_conf    ! Number of paralell configurations in SA
      character(len=60) :: Cost_function_name
      integer           :: seed        ! If different from zero, holds the seed for random number generator
    End type SimAnn_Conditions_type


    !!----
    !!----  Type, public  :: MultiState_Vector_Type
    !!----    integer                            :: npar    ! Number of parameters of the model.
    !!----    integer                            :: nconf   ! Number of configurations.
    !!----    integer, dimension(np_SAN,np_CONF) :: code    !=0 fixed parameter, =1 variable parameter
    !!----    integer, dimension(np_SAN)         :: bound   !=0 fixed boundaries,=1 periodic boundaries
    !!----    real,    dimension(np_SAN,np_CONF) :: state   !Vector State with the current configuration
    !!----    real,    dimension(np_SAN,np_CONF) :: stp     !Step vector (one value for each parameter)
    !!----    real,    dimension(np_SAN)         :: low     !Low-limit value of parameters
    !!----    real,    dimension(np_SAN)         :: high    !High-limit value of parameters
    !!----    real,    dimension(np_SAN)         :: cost    !Vector with cost of the different configurations
    !!----    real,    dimension(np_SAN)         :: config  !Vector State with the best configuration
    !!----    character(len=15),dimension(np_SAN):: nampar !name of parameters of the model
    !!----  End Type MultiState_Vector_Type
    !!----
    !!---- Derived type containing the parameters and configurations to be optimized,
    !!---- the limits, steps, names and best configuration to be searched
    !!-----by Simulated Annealing Algorithm
    !!----
    !!---- Update: March - 2005
    !!
    Type, public  :: MultiState_Vector_Type
      integer                                     :: npar    ! Number of parameters of the model. To be supplied in the calling program
      integer                                     :: nconf   ! Number of configurations.
      integer,          dimension(np_SAN)         :: code    !=0 fixed parameter, =1 variable parameter
      integer,          dimension(np_SAN)         :: bound   !=0 fixed boundaries, =1 periodic boundaries
      real,             dimension(np_SAN,np_CONF) :: state   !Vector State characterizing the current configuration
      real,             dimension(np_SAN,np_CONF) :: stp     !Step vector (one value for each parameter)
      real,             dimension(np_SAN)         :: low     !Low-limit value of parameters
      real,             dimension(np_SAN)         :: high    !High-limit value of parameters
      real,             dimension(np_SAN)         :: cost    !Vector with cost of the different configurations
      real,             dimension(np_SAN)         :: config  !Vector State characterizing the current best configuration
      character(len=15),dimension(np_SAN)         :: nampar !name of parameters of the model
    End Type MultiState_Vector_Type

    !!----
    !!----  Type, public  :: State_Vector_Type
    !!----    integer                    :: npar    ! Number of parameters of the model.
    !!----    integer, dimension(np_SAN) :: code    !=0 fixed parameter, =1 variable parameter
    !!----    integer, dimension(np_SAN) :: bound   !=0 fixed boundaries,=1 periodic boundaries
    !!----    real,    dimension(np_SAN) :: state   !Vector State with the current configuration
    !!----    real,    dimension(np_SAN) :: stp     !Step vector (one value for each parameter)
    !!----    real,    dimension(np_SAN) :: low     !Low-limit value of parameters
    !!----    real,    dimension(np_SAN) :: high    !High-limit value of parameters
    !!----    real,    dimension(np_SAN) :: config  !Vector State with the best configuration
    !!----    character(len=15),dimension(np_SAN):: nampar !name of parameters of the model
    !!----  End Type State_Vector_Type
    !!----
    !!---- Derived type containing the parameters to be optimized,
    !!---- the limits, steps, names and best configuration to be searched
    !!-----by Simulated Annealing Algorithm
    !!----
    !!---- Update: March - 2005
    !!
    Type, public  :: State_Vector_Type
      integer                             :: npar    ! Number of parameters of the model. To be supplied in the calling program
      integer,          dimension(np_SAN) :: code    !=0 fixed parameter, =1 variable parameter
      integer,          dimension(np_SAN) :: bound   !=0 fixed boundaries, =1 periodic boundaries
      real,             dimension(np_SAN) :: state   !Vector State characterizing the current configuration
      real,             dimension(np_SAN) :: stp     !Step vector (one value for each parameter)
      real,             dimension(np_SAN) :: low     !Low-limit value of parameters
      real,             dimension(np_SAN) :: high    !High-limit value of parameters
      real,             dimension(np_SAN) :: config  !Vector State characterizing the current best configuration
      character(len=15),dimension(np_SAN) :: nampar !name of parameters of the model
    End Type State_Vector_Type


    !!----
    !!---- ERR_MESS_SAN
    !!----    character (len=150), public  :: err_mess_SAN
    !!----
    !!----    Error message provided if an error condition occurs.
    !!----
    !!---- Update: February - 2005
    !!
    character (len=150), public   :: err_mess_SAN=" "

    !!----
    !!---- ERR_SAN
    !!----    logical , public  :: err_SAN
    !!----
    !!----    Error flag provided if an error condition occurs
    !!----
    !!---- Update: March - 2005
    !!
    logical, public  :: err_SAN =.false.

 Contains

    !---- Functions ----!

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Checkm(c,vs)
    !!--++   type(SimAnn_Conditions_type), intent(in out) :: c
    !!--++   type(State_Vector_Type),      intent(in out) :: vs
    !!--++
    !!--++    (PRIVATE)
    !!--++    Check if all important variables of the SA algorithm have
    !!--++    been supplied. Fix some values
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Checkm(c,vs)
       type(SimAnn_Conditions_type), intent(in out) :: c
       type(MultiState_Vector_Type), intent(in out) :: vs

       !---- Local variables ----!
       integer :: i
       character(len=3) :: num

       err_san=.false.
       if (vs%npar == 0) then
          err_san=.true.
          err_mess_san=" Zero parameters for the vector state: PROVIDE a number > 0 for variable: NPAR"
          return
       end if
       if (c%nm_cycl <= 1) then
          err_san=.true.
          err_mess_san=" Too small value for the number of MCcycles/Temp : PROVIDE a number > 1 for variable: NM_CYCL"
          return
       end if

       !---- Default values ----!
       if (c%num_temps < 1)   c%num_temps=1
       if (c%num_therm < 0)   c%num_therm=0
       if (c%initconfig < 0)  c%initconfig=0
       if (c%nalgor < 0)      c%nalgor=0
       if (c%anneal <= 0.001) c%anneal=0.9
       if (c%t_ini <= 0.0 )   c%t_ini=5.0
       if (c%accept <= 0.0 )  c%accept=0.01
       if (len_trim(c%Cost_function_name) == 0) then
          c%Cost_function_name="Unnamed Cost_function"
       end if

       do i=1,vs%npar
          if (len_trim(vs%nampar(i)) == 0) then
             write(unit=num,fmt="(i3)") i
             num=adjustl(num)
             vs%nampar(i)="Par_"//num
          end if
          if (abs(vs%low(i)-vs%high(i)) < 1.0e-12) vs%code(i)=0
       end do

       return
    End Subroutine Checkm

    !!--++
    !!--++ Subroutine Check(c,vs)
    !!--++   type(SimAnn_Conditions_type), intent(in out) :: c
    !!--++   type(State_Vector_Type),      intent(in out) :: vs
    !!--++
    !!--++    (PRIVATE)
    !!--++    Chek if all important variables of the SA algorithm has
    !!--++    been supplied. Fix some values
    !!--++
    !!--++ Update: March - 2005
    !!
    Subroutine Check(c,vs)
       type(SimAnn_Conditions_type), intent(in out) :: c
       type(State_Vector_Type),      intent(in out) :: vs

       !---- Local variables ----!
       integer :: i
       character(len=3) :: num

       err_san=.false.
       if (vs%npar == 0) then
          err_san=.true.
          err_mess_san=" Zero parameters for the vector state: PROVIDE a number > 0 for variable: NPAR"
          return
       end if
       if (c%nm_cycl <= 1) then
          err_san=.true.
          err_mess_san=" Too small value for the number of MCcycles/Temp : PROVIDE a number > 1 for variable: NM_CYCL"
          return
       end if

       !---- Default values ----!
       if (c%num_temps < 1)   c%num_temps=1
       if (c%num_therm < 0)   c%num_therm=0
       if (c%initconfig < 0)  c%initconfig=0
       if (c%nalgor < 0)      c%nalgor=0
       if (c%anneal <= 0.001) c%anneal=0.9
       if (c%t_ini <= 0.0 )   c%t_ini=5.0
       if (c%accept <= 0.0 )  c%accept=0.01
       if (len_trim(c%Cost_function_name) == 0) then
          c%Cost_function_name="Unnamed Cost_function"
       end if

       do i=1,vs%npar
          if (len_trim(vs%nampar(i)) == 0) then
             write(unit=num,fmt="(i3)") i
             num=adjustl(num)
             vs%nampar(i)="Par_"//num
          end if
          if (abs(vs%low(i)-vs%high(i)) < 1.0e-12) vs%code(i)=0
       end do

       return
    End Subroutine Check

    !!--++
    !!--++ Subroutine Init_Ran(seed)
    !!--++    integer, dimension(1), intent (in) :: seed
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--..    Calling sequence: Call RANDOM_SEED(size=k,put=seed(1:k),get=old(1:k))
    !!--..    where:
    !!--..          size: Integer, intent(out) - default integer size used by the processor to hold the seed
    !!--..          put : Integer, dimension(:), intent(in)- used by the processor to set the seed value
    !!--..          get : Integer, dimension(:), intent(out)- set by the processor to current seed value
    !!--..          size(put) and size(get) >= size
    !!--++
    !!--++ Update: February - 2003
    !!
    Subroutine Init_Ran(Seed)
       !---- Argument ----!
       integer, dimension(1), intent (in) :: seed


       if (seed(1) == 0) then
          call random_seed() !seed selected by the system clock
       else
          call random_seed(put=seed) !seed selected by the user
       end if

       return
    End Subroutine Init_Ran


    !!----
    !!---- Subroutine Set_SimAnn_Cond(file_list,c)
    !!----    type (file_list_type),       intent( in)  :: file_list
    !!----    type(SimAnn_Conditions_type),intent(out)  :: c
    !!----    Type, public        :: SimAnn_Conditions_type
    !!----      real              :: t_ini=0.0   ! Initial temperature
    !!----      real              :: anneal=0.0  ! Kirpactrick factor for Annealing
    !!----      real              :: accept=0.0  ! Minimum percentage of accepted configurations
    !!----      integer           :: initconfig=0! Flag determining if the first configuration is random or read
    !!----      integer           :: nalgor=0    ! Flag determining if the Corana algorithm is selected (0) or not (/=0)
    !!----      integer           :: nm_cycl=0   ! Number of Cycles per temp  in SA searchs
    !!----      integer           :: num_temps=0 ! Maximum number of temperatures in SA
    !!----      integer           :: num_therm=0 ! Number of thermalization cycles in SA
    !!----      character(len=60) :: Cost_function_name=" "
    !!----      integer           :: seed=0      ! If different from zero, holds the seed for random number generator
    !!----    End type SimAnn_Conditions_type
    !!----
    !!---- Subroutine for reading and set up the SimAnn_Conditions_type
    !!---- variable "c"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Set_SimAnn_Cond(file_list,c)
       !---- Arguments ----!
       type(file_list_type),        intent( in)  :: file_list
       type(SimAnn_Conditions_type),intent(out)  :: c

       !--- Local Variables ---!
       integer :: i,j,ier
       logical :: notset,tempar,algor,noinst
       character(len=80) :: line

        notset=.true.
        tempar=.false.
        algor =.false.

        c%t_ini=5.0   ! Initial temperature
        c%anneal=0.9  ! Kirpactrick factor for Annealing
        c%accept=0.01 ! Minimum percentage of accepted configurations
        c%initconfig=0! Flag determining if the first configuration is random or read
        c%nalgor=0    ! Flag determining if the Corana algorithm is selected (0) or not (/=0)
        c%nm_cycl=0   ! Number of Cycles per temp  in SA searchs
        c%num_temps=1 ! Maximum number of temperatures in SA
        c%num_therm=0 ! Number of thermalization cycles in SA
        c%num_conf=1  ! Number of paralell configurations in SA
        c%Cost_function_name=" Unnamed Cost Function"
        c%seed=0      ! If different from zero, holds the seed for random number generator

        do_read: do i=1,file_list%nlines
          if(file_list%line(i)(1:7) == "SIM_ANN") then

            do j=i+1,file_list%nlines
               line=u_case(file_list%line(j))

               Select Case (line(1:7))

                  Case("COSTNAM")
                     c%Cost_function_name=adjustl(file_list%line(j)(8:))

                  Case("TEMPARM")
                     read(unit=line(8:),fmt=*,iostat=ier) c%T_ini,c%anneal,c%num_temps
                     if(ier /= 0) exit do_read
                     tempar=.true.

                  Case("ALGOR_T")
                     read(unit=line(8:),fmt=*,iostat=ier) c%nalgor,c%num_conf,c%nm_cycl, c%num_therm, c%accept
                     if(ier /= 0) exit do_read
                     algor=.true.

                  Case("SEEDVAL")
                     read(unit=line(8:),fmt=*,iostat=ier) c%seed
                     if(ier /= 0) c%seed=0

                  Case("INITCON")
                     line = adjustl(line(8:))
                     if(line(1:3) /= "RAN") c%initconfig = 1

               End Select
               if(tempar .and. algor) notset=.false.
            end do
            exit do_read
          end if
          noinst=.true.
        end do do_read
        if(notset) then
          err_SAN =.true.
          if(noinst) then
            err_mess_SAN=" => No Simulated Annealing conditions in input file "
          else if(.not. tempar) then
            err_mess_SAN=&
            " => Unable to set Simulated Annealing conditions (Error in line: TemParM T_ini anneal num_temps num_therm)"
          else if(.not. algor) then
            err_mess_SAN= &
            " => Unable to set Simulated Annealing conditions (Error in line: Algor_T  nalgor  num_conf  nm_cycl   num_therm)"
          else
            err_mess_SAN=" Unable to set Simulated Annealing conditions (Error in CFL file) "
          end if
        end if
       return
    End Subroutine Set_SimAnn_Cond

    !!----
    !!---- Subroutine Set_SimAnn_MStateV(n,nsol,Con,Bounds,VNam,Vec,vs,cod)
    !!----    integer,                      intent(in) :: n,nsol !number of parameters & configurations
    !!----    integer,         dimension(:),intent(in) :: con    !Boundary conditions
    !!----    real,          dimension(:,:),intent(in) :: Bounds ! (1,:)-> Low, (2,:) -> High, (3,:) -> Step
    !!----    character(len=*),dimension(:),intent(in) :: VNam   !Names of parameters
    !!----    real,            dimension(:),intent(in) :: Vec    !Initial value of parameters
    !!----    type(MultiState_Vector_Type), intent(out):: vs     !Initial State vector
    !!----    integer,optional,dimension(:),intent(in) :: cod    !If present, cod(i)=0 fix the "i" parameter
    !!----
    !!----
    !!---- Subroutine for setting up the State_Vector_type
    !!---- variable "vs"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Set_SimAnn_MStateV(n,nsol,Con,Bounds,VNam,Vec,vs,cod)
       !---- Arguments ----!
       integer,                      intent(in) :: n,nsol
       integer,         dimension(:),intent(in) :: con
       real,          dimension(:,:),intent(in) :: Bounds
       character(len=*),dimension(:),intent(in) :: VNam
       real,            dimension(:),intent(in) :: Vec
       type(MultiState_Vector_Type), intent(out):: vs
       integer,optional,dimension(:),intent(in) :: cod

       integer :: j
       if(n > np_SAN) then
         err_SAN =.true.
         write(unit=err_mess_SAN,fmt="(a,i4,a)") " => Too many parameters! Simulated Anneling module limited to ",&
              np_SAN," parameters"
         return
       end if
       vs%nconf= 0
       vs%npar = 0
       vs%code(:) = 1
       vs%bound(:) = 0
       vs%state(:,:) = 0.0
       vs%low(:) = 0.0
       vs%high(:) = 0.0
       vs%stp(:,:) = 0.0
       vs%config(:) = 0.0
       vs%cost(:) =0.0
       vs%Nampar(:) = " "

         ! Copying arguments in state vector
           vs%npar        = n
           vs%nconf       = nsol
           vs%code(1:n)   = 1
          vs%bound(1:n)   = Con(1:n)
            vs%low(1:n)   = Bounds(1,1:n)
           vs%high(1:n)   = Bounds(2,1:n)
           do j=1,nsol
              vs%stp(1:n,j) = Bounds(3,1:n)
            vs%state(1:n,j) = vec(1:n)
           end do
         vs%config(1:n)   = vec(1:n)
         vs%Nampar(1:n)   = vNam(1:n)
         vs%Nampar(1:n)   = adjustr(vs%Nampar(1:n))
         if(present(cod)) then
           vs%code(1:n) = Cod(1:n)
         end if

       return
    End Subroutine Set_SimAnn_MStateV

    !!----
    !!---- Subroutine Set_SimAnn_StateV(n,Con,Bounds,VNam,Vec,vs,cod)
    !!----    integer,                      intent(in) :: n      !number of parameters
    !!----    integer,         dimension(:),intent(in) :: con    !Boundary conditions
    !!----    real,          dimension(:,:),intent(in) :: Bounds ! (1,:)-> Low, (2,:) -> High, (3,:) -> Step
    !!----    character(len=*),dimension(:),intent(in) :: VNam   !Names of parameters
    !!----    real,            dimension(:),intent(in) :: Vec    !Initial value of parameters
    !!----    type(State_Vector_Type),      intent(out):: vs     !Initial State vector
    !!----    integer,optional,dimension(:),intent(in) :: cod    !If present, cod(i)=0 fix the "i" parameter
    !!----
    !!----
    !!---- Subroutine for setting up the State_Vector_type
    !!---- variable "vs"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Set_SimAnn_StateV(n,Con,Bounds,VNam,Vec,vs,cod)
       !---- Arguments ----!
       integer,                      intent(in) :: n
       integer,         dimension(:),intent(in) :: con
       real,          dimension(:,:),intent(in) :: Bounds
       character(len=*),dimension(:),intent(in) :: VNam
       real,            dimension(:),intent(in) :: Vec
       type(State_Vector_Type),      intent(out):: vs
       integer,optional,dimension(:),intent(in) :: cod

       if(n > np_SAN) then
         err_SAN =.true.
         write(unit=err_mess_SAN,fmt="(a,i4,a)") " => Too many parameters! Simulated Anneling module limited to ",&
              np_SAN," parameters"
         return
       end if

       vs%npar = 0
       vs%code(:) = 1
       vs%bound(:) = 0
       vs%state(:) = 0.0
       vs%low(:) = 0.0
       vs%high(:) = 0.0
       vs%stp(:) = 0.0
       vs%config(:) = 0.0
       vs%Nampar(:) = " "

         ! Copying arguments in state vector
           vs%npar      = n
           vs%code(1:n) = 1
          vs%bound(1:n) = Con(1:n)
          vs%state(1:n) = vec(1:n)
            vs%low(1:n) = Bounds(1,1:n)
           vs%high(1:n) = Bounds(2,1:n)
            vs%stp(1:n) = Bounds(3,1:n)
         vs%config(1:n) = vec(1:n)
         vs%Nampar(1:n) = vNam(1:n)
         vs%Nampar(1:n) = adjustr(vs%Nampar(1:n))
         if(present(cod)) then
           vs%code(1:n) = Cod(1:n)
         end if

       return
    End Subroutine Set_SimAnn_StateV


    !!----
    !!---- Subroutine Simanneal_Gen(Model_Funct,c,vs,Ipr,fileSav)
    !!----    type(SimAnn_Conditions_type),intent(in out)  :: c
    !!----    type(State_Vector_Type),     intent(in out)  :: vs
    !!----    integer,                     intent(in)      :: Ipr
    !!----    character(len=*), optional,  intent(in)      :: filesav
    !!----
    !!----    Interface
    !!----       Subroutine Model_Funct(v,cost)
    !!----          real,dimension(:),    intent( in):: v
    !!----          real,                 intent(out):: cost
    !!----       End Subroutine Model_Funct
    !!----    End Interface
    !!----
    !!---- Update: March - 2005
    !!

    Subroutine Simanneal_Gen(Model_Funct,c,vs,Ipr,fileSav)
       !---- Arguments ----!
       type(SimAnn_Conditions_type),intent(in out)  :: c
       type(State_Vector_Type),     intent(in out)  :: vs
       integer,                     intent(in)      :: Ipr
       character(len=*), optional,  intent(in)      :: filesav

       Interface
          Subroutine Model_Funct(v,cost)
             real,dimension(:),    intent( in):: v
             real,                 intent(out):: cost
          End Subroutine Model_Funct
       End Interface

       !--- Local Variables ---!
       character (len=132) :: messag, strings
       integer             :: i, j, neval, ncf, jk, naj, ntp, naver, last, jj
       real                :: temp, cost, cost1, cost2, costop, ep, ener, cp, dsen, dsen2, &
                              energ, paj, prob, rav, rati, plage, shift, stepav,random
       integer, parameter      :: i_conf=99
       integer, dimension(1)   :: seed
       real,    dimension(np_SAN) :: stateo     !Vector State characterizing the old configuration
       real,    dimension(np_SAN) :: cv         !constant vector used in the Corana algorithm
       integer, dimension(np_SAN) :: nacp       !number of accepted moves for parameter i
       real,    dimension(np_SAN) :: raver      !Vector State characterizing the average configuration
       real,    dimension(np_SAN) :: sigp       !Standard deviations of the average configuration

       call check(c,vs)
       if (err_san) then
          call mess(err_mess_san)
          return
       end if
       seed(1)=c%seed
       call init_ran(seed)
       jj=min(vs%npar,26)

       if (present(filesav)) then
          open(unit=i_conf,file=trim(filesav)//".anl",status="replace",action="write")
          write(unit=i_conf,fmt="(a,a)")"! Simulated Anneling Run with cost function: ", trim(c%Cost_function_name)
          write(unit=i_conf,fmt="(a,26(a,a))") &
                   "! NT    Temp  <Cost-val>   %Accpt     <Step>     Cp   ",("  ",vs%nampar(i),i=1,jj)
          call flush(i_conf)
          close(unit=i_conf)
       end if

       if (c%nalgor == 0) then
          vs%stp(1:vs%npar)=vs%high(1:vs%npar)-vs%low(1:vs%npar)
       end if
       cv(:)=2.0

       !---- Get the initial configuration in config(1:msz) ----!
       if (c%initconfig == 1) then
          do i=1,vs%npar
             vs%state(i) = vs%config(i)
          end do
       else
          do i=1,vs%npar
             if (vs%code(i)==1) then
                call random_number(random)
                vs%state(i) = vs%low(i) + random*(vs%high(i)-vs%low(i))
            else
                vs%state(i)=vs%config(i)
             end if
          end do
       end if
       stateo(:)=vs%state(:)

       !---- Determine the initial value of the cost function -> cost1 ----!
       call Model_Funct(stateo,Cost1)

       costop=cost1
       vs%config(:)=stateo(:)

       messag=" ---- Simulated Annealing to minimize General Cost Functions ----"
       write(unit=ipr,fmt="(/,a,/,a,/)") messag, "     Cost-function name:  "//trim(c%Cost_function_name)
       call mess(" ")
       call mess(messag)
       call mess(" ")

       strings=" "
       write(unit=strings,fmt="(a,f16.3)") " => Initial configuration cost: ",cost1
       call mess(strings)
       write(unit=ipr,fmt="(a)") strings

       strings=" "
       write(unit=strings,fmt="(a)") " => Initial configuration state vector: "
       call mess(strings)
       write(unit=ipr,fmt="(a)") strings

       do i=1,vs%npar
          write(unit=strings,fmt="(i6,a,F16.5)")  i,vs%nampar(i), vs%state(i)
          write(unit=ipr,fmt="(a)") strings
          call mess(strings)
       end do

       !---- Loop over temperatures ----!
       naver=c%nm_cycl-c%num_therm
       rav= real(naver*vs%npar)
       temp=c%t_ini/c%anneal
       neval=0
       nacp=0
       do ntp=1,c%num_temps     ! Global DO for changing temperature
          naj=0
          cost=0.0
          energ=0.0
          sigp(:)=0.0
          raver(:)=0.0
          dsen=0.0
          dsen2=0.0

          temp=c%anneal*Temp  ! Current temperature

          do ncf=1,c%nm_cycl
             last=vs%npar

             cyc_par: do i=1,vs%npar        !loop on the components of the vector state
                if (vs%code(i) == 0) cycle cyc_par
                !---- Set new configuration and store the previous one in config(1:msz) ----!
                call random_number(random)
                vs%state(i)=stateo(i)+vs%stp(i)*(2.0*random-1.0)
                plage=vs%high(i)-vs%low(i)
                if (vs%bound(i) == 0) then   !boundary conditions
                   if (vs%state(i) < vs%low(i) .or. vs%state(i) > vs%high(i))  &  !Fixed boundaries
                      call random_number(random)
                   vs%state(i)=vs%low(i)+ random*plage
                else
                   if (vs%state(i) < vs%low(i) ) vs%state(i)=vs%state(i)+plage  !Periodic boundary conditions
                   if (vs%state(i) > vs%high(i) ) vs%state(i)=vs%state(i)-plage
                end if
                shift=vs%state(i)-stateo(i)

                !---- Calculate the cost function ----!
                call Model_Funct(vs%state,Cost2)
                neval=neval+1
                ener= cost2-cost1

                !---- Metropolis test ----!
                if (ener > 0.0 ) then
                   ep=ener/temp
                   if (ep > 88.7228) then
                      prob=0.0
                   else
                      prob=exp(-ep)
                   end if
                   call random_number(random)
                   if (prob <= random) then
                      !---- Restore the old configuration ----!
                      vs%state(i)=stateo(i)
                      cycle cyc_par   !cycle and try another configuration
                   end if
                end if

                !---- Accepted configuration ----!
                cost1=cost2               !update the cost function
                stateo(i)=vs%state(i)     !update configuration
                if (cost1 < costop) then  !the best current configuration is found
                   costop=cost1
                   vs%config(:)=stateo(:)
                end if

                if (ncf <= c%Num_therm) cycle cyc_par
                nacp(i)=nacp(i)+1
                naj=naj+1   !number of accepted jumps
                energ=energ+abs(ener)
                cost=cost+cost2
                dsen=dsen+ener
                dsen2=dsen2+ener*ener

             end do cyc_par  !end loop over components of state vector

             do j=1,vs%npar
                sigp(j)=sigp(j)+stateo(j)*stateo(j)
                raver(j)=raver(j)+stateo(j)
             end do

          end do    !End loop over Montecarlo moves at fixed temperature

          !---- Statistic  and average values for previous temperature ----!
          if (naj == 0) then
             naj=1
             cost=cost1
             energ=abs(ener)
             dsen=ener
             dsen2=ener*ener
          end if
          paj=100.0*real(naj)/rav

          energ=energ/naj
          cost=cost/naj
          dsen=dsen/naj
          dsen2=dsen2/naj
          cp=(dsen2-dsen*dsen)/(temp*temp)
          do j=1,vs%npar
             raver(j)=raver(j)/naver
             sigp(j)=sqrt(abs(sigp(j)/naver-raver(j)*raver(j)))
          end do
          stepav=sum(vs%stp(:))/real(vs%npar)

          !---- Writing partial results ----!
          strings=" "
          write(unit=strings,fmt="(a,i4,2(a,f8.3),a,f10.5,a,f12.4)")  &
               " => NT:",ntp," Temp:",temp," (%Acc):",paj,"  <Step>:",stepav,"  <"//trim(c%Cost_function_name)//">:",cost
          call mess(strings)
          write(unit=ipr,fmt="(/,a)") strings
          strings=" "
          write(unit=strings,fmt="(a,f12.5,a,f16.5,a,i10 )")  &
               "    Cost-val:",energ, "     Cp:", cp, "     Num-Cost-Evaluations: ",neval
          write(unit=ipr,fmt="(a)") strings

          write(unit=ipr,fmt="(a)") " => Average value of parameters, sigmas and steps:"
          do i=1,vs%npar
            write(unit=strings,fmt="(i6,a,3F16.5)")  i,vs%nampar(i), raver(i),sigp(i),vs%stp(i)
            write(unit=ipr,fmt="(a)") strings
          end do


          jk= min(vs%npar,26)
          if (present(filesav)) then
             open(unit=i_conf,file=trim(filesav)//".anl",status="old",action="write", position="append")
             write(unit=i_conf,fmt="(i4,31f10.4)")  ntp,temp,cost,paj,stepav,cp,(raver(jj),jj=1,jk)
             call flush(i_conf)
             close(unit=i_conf)
          end if
          if (paj-c%accept <=0 ) exit   !Convergence criterium

          !---- Adjust STP so that approximately half of all evaluations are accepted ----!
          !---- (Corana's algorithm)

          if (c%nalgor == 0 .or. c%nalgor == 1) then
             do i=1,vs%npar
                jj=0
                if (vs%stp(i) > abs(0.01*c%accept*raver(i))) exit
                jj=1
             end do
             if (jj == 1) exit
             do i = 1, vs%npar
                plage=vs%high(i)-vs%low(i)
                rati = real(nacp(i)) /real(naver)
                if (rati > 0.6) then
                   vs%stp(i) = vs%stp(i)*(1.0 + cv(i)*(rati - 0.6)/0.4)
                else if (rati < 0.4) then
                   vs%stp(i) = vs%stp(i)/(1.0 + cv(i)*((0.4 - rati)/0.4))
                end if
                if (vs%stp(i) > plage) then
                   vs%stp(i) = plage
                end if
             end do
          end if
          nacp(:) = 0

       end do   !ntp=1,c%num_temps

       !---- Re-calculate the cost function for the best configuration ----!
       call Model_Funct(vs%config,Cost2)

       messag="  "
       call mess(messag)
       messag=" => Best configuration found by Simulated Annealing (Sigma of the last Montecarlo Cycles):"
       call mess(messag)
       write(unit=ipr,fmt="(/,a,a,/)") "     ",trim(messag)
       strings=" "
       write(unit=strings,fmt="(a,f10.4,a)")" -> Best Solution Cost =",cost2," :: "
       call mess(strings)
       write(unit=ipr,fmt="(/,a)") strings
       messag=" -> Configuration parameters :"
       call mess(messag)
       write(unit=ipr,fmt="(a,/)") trim(messag)
       do i=1,vs%npar
         write(unit=strings,fmt="(i6,a,2F16.5)")  i,vs%nampar(i), vs%config(i),sigp(i)
         write(unit=ipr,fmt="(a)") strings
         call mess(strings)
       end do

       return
   End Subroutine Simanneal_Gen

    !!----
    !!---- Subroutine SimAnneal_MultiConf(Model_Funct,Nsol,c,vs,Ipr,fileSav)
    !!----    integer,                       intent(in out)  :: Nsol
    !!----    type(SimAnn_Conditions_type),  intent(in out)  :: c
    !!----    type(MultiState_Vector_Type),  intent(in out)  :: vs
    !!----    integer,                       intent(in)      :: Ipr
    !!----    character(len=*), optional,    intent(in)      :: filesav
    !!----
    !!----    Interface
    !!----       Subroutine Model_Funct(v,cost)
    !!----          real,dimension(:),    intent( in):: v
    !!----          real,                 intent(out):: cost
    !!----       End Subroutine Model_Funct
    !!----    End Interface
    !!----
    !!---- Update: March - 2005
    !!

    Subroutine SimAnneal_MultiConf(Model_Funct,c,vs,Ipr,fileSav)
       !---- Arguments ----!
       type(SimAnn_Conditions_type),  intent(in out)  :: c
       type(MultiState_Vector_Type),  intent(in out)  :: vs
       integer,                       intent(in)      :: Ipr
       character(len=*), optional,    intent(in)      :: filesav

       Interface
          Subroutine Model_Funct(v,cost)
             real,dimension(:),    intent( in):: v
             real,                 intent(out):: cost
          End Subroutine Model_Funct
       End Interface

       !--- Local Variables ---!
       character (len=256) :: messag, strings
       integer             :: i, j, k, neval, ncf, ntp, naver, last, jj, &
                              survive,jopt
       real                :: temp, ep, ener, costop, half_init_avstp, sumdel,sumsig, &
                              prob, rav, rati, plage, shift, stepav,random
       integer, parameter                 :: i_conf=99
       integer, dimension(1)              :: seed
       logical, dimension(np_CONF)        :: dead
       integer, dimension(np_CONF)        :: naj
       real,    dimension(np_CONF)        :: cost, cost1, cost2, paj, costav
       real,    dimension(np_SAN,np_CONF) :: stateo     !Vector State characterizing the old configuration
       real,    dimension(np_SAN)         :: cv         !constant vector used in the Corana algorithm
       integer, dimension(np_SAN,np_CONF) :: nacp       !number of accepted moves for parameter i
       real,    dimension(np_SAN,np_CONF) :: raver      !Vector State characterizing the average configuration
       real,    dimension(np_SAN,np_CONF) :: sigp       !Standard deviations of the average configuration

       call checkm(c,vs)
       if (err_san) then
          call mess(err_mess_san)
          return
       end if

       seed(1)=c%seed
       call init_ran(seed)
       jj=min(vs%npar,26)

       if (present(filesav)) then
          open(unit=i_conf,file=trim(filesav)//".anl",status="replace",action="write")
          write(unit=i_conf,fmt="(a,a)")"! Simulated Anneling Run with cost function: ", trim(c%Cost_function_name)
          write(unit=i_conf,fmt="(a,26(a,a))") &
                   "! NT    Temp  <Cost-val>   %Accpt     <Step>     Cp   ",("  ",vs%nampar(i),i=1,jj)
          call flush(i_conf)
          close(unit=i_conf)
       end if

       cv(1:vs%npar)=vs%high(1:vs%npar)-vs%low(1:vs%npar)
       half_init_avstp=0.5*sum(cv(1:vs%npar)/real(vs%npar))

       if (c%nalgor == 0) then
         do j=1,vs%nconf
          vs%stp(1:vs%npar,j)=cv(1:vs%npar)
         end do
       else
         half_init_avstp=4.0*half_init_avstp
       end if

       cv(:)=2.0
       stateo=0.0
       !---- Get the initial configuration in config(1:msz) ----!
       dead(:)=.false.
       if (c%initconfig == 1) then
          do j=1,vs%nconf
            do i=1,vs%npar
               vs%state(i,j) = vs%config(i)
            end do
          end do
       else
          do j=1,vs%nconf
            do i=1,vs%npar
               if (vs%code(i)==1) then
                  call random_number(random)
                  vs%state(i,j) = vs%low(i) + random*(vs%high(i)-vs%low(i))
               else
                  vs%state(i,j)=vs%config(i)
               end if
               stateo(i,j)=vs%state(i,j)
            end do
          end do
       end if

       !---- Determine the initial values of the cost function -> cost1 ----!
       do j=1,vs%nconf
         call Model_Funct(stateo(:,j),Cost1(j))
       end do
         j=minloc(cost1(1:vs%nconf),dim=1)
         vs%config(:)=stateo(:,j)  !Best configuration for the moment
         costop=cost1(j)
         jopt=j

       messag=" ---- MultiConfiguration Simulated Annealing to minimize General Cost Functions ----"
       write(unit=ipr,fmt="(/,a,/,a,/)") messag, "     Cost-function name:  "//trim(c%Cost_function_name)
       call mess(" ")
       call mess(messag)
       call mess(" ")

       do j=1,vs%nconf
         strings=" "
         write(unit=strings,fmt="(a,i2,a,f17.4)") " => Initial configuration cost(",j,"): ",cost1(j)
         call mess(strings)
         write(unit=ipr,fmt="(a)") strings
       end do
       strings=" "
       write(unit=strings,fmt="(a)") " => Initial Best configuration state vector: "
       call mess(strings)
       write(unit=ipr,fmt="(a)") strings

       strings=" "
       do i=1,vs%npar
         write(unit=strings,fmt="(i6,a,f16.5)")  i,vs%nampar(i), vs%config(i)
         write(unit=ipr,fmt="(a)") strings
         call mess(strings)
       end do

       !---- Loop over temperatures ----!
       naver=c%nm_cycl-c%num_therm
       rav= real(naver*vs%npar)
       temp=c%t_ini/c%anneal
       neval=0
       nacp=0
       costav(:)=0.0
       do ntp=1,c%num_temps     ! Global DO for changing temperature
            naj(:)   =0
           cost(:)   =0.0
           sigp(:,:) =0.0
          raver(:,:) =0.0

          temp=c%anneal*Temp  ! Current temperature

          strings=" "
          write(unit=strings,fmt="(a,f9.5,a,i5,a,i8)") "  => New Temp:",temp,"  NT: ",ntp, &
               "       Number of function evaluations:",neval
          call mess(strings)
          write(unit=ipr,fmt="(/,a)") strings

          do j=1,vs%nconf   ! loop over different configuration state vectors

             if(dead(j)) cycle  !skip is configuration is dead

             do ncf=1,c%nm_cycl
                last=vs%npar

                cyc_par: do i=1,vs%npar        !loop on the components of the vector state
                   if (vs%code(i) == 0) cycle cyc_par
                   !---- Set new configuration and store the previous one in config(1:msz) ----!
                   call random_number(random)
                   vs%state(i,j)=stateo(i,j)+vs%stp(i,j)*(2.0*random-1.0)
                   plage=vs%high(i)-vs%low(i)
                   if (vs%bound(i) == 0) then   !boundary conditions
                      if (vs%state(i,j) < vs%low(i) .or. vs%state(i,j) > vs%high(i))  &  !Fixed boundaries
                         call random_number(random)
                      vs%state(i,j)=vs%low(i)+ random*plage
                   else
                      if (vs%state(i,j) < vs%low(i) )  vs%state(i,j)=vs%state(i,j)+plage  !Periodic boundary conditions
                      if (vs%state(i,j) > vs%high(i) ) vs%state(i,j)=vs%state(i,j)-plage
                   end if
                   shift=vs%state(i,j)-stateo(i,j)

                   !---- Calculate the cost function ----!
                   call Model_Funct(vs%state(:,j),Cost2(j))
                   neval=neval+1
                   ener= cost2(j)-cost1(j)

                   !---- Metropolis test ----!
                   if (ener > 0.0 ) then
                      ep=ener/temp
                      if (ep > 88.7228) then
                         prob=0.0
                      else
                         prob=exp(-ep)
                      end if
                      call random_number(random)
                      if (prob <= random) then
                         !---- Restore the old configuration ----!
                         vs%state(i,j)=stateo(i,j)
                         cycle cyc_par   !cycle and try another configuration
                      end if
                   end if

                   !---- Accepted configuration ----!
                   cost1(j)=cost2(j)               !update the cost function
                   stateo(i,j)=vs%state(i,j)       !update configuration
                   if (cost1(j) < costop) then     !the best current configuration is found
                      costop=cost1(j)
                      vs%config(:)=stateo(:,j)
                      jopt=j
                   end if

                   if (ncf <= c%Num_therm) cycle cyc_par
                   nacp(i,j)=nacp(i,j)+1
                   naj(j)=naj(j)+1   !number of accepted jumps for state vector j
                   cost(j)=cost(j)+cost2(j)

                end do cyc_par  !end loop over components of state vector


                do k=1,vs%npar
                   sigp(k,j)= sigp(k,j)+stateo(k,j)*stateo(k,j)
                  raver(k,j)=raver(k,j)+stateo(k,j)
                end do

             end do    !End loop over Montecarlo moves at fixed temperature

          end do   !end loop over configurations

          vs%cost(1:vs%nconf)=cost1(1:vs%nconf)
          if(modulo(ntp,4) == 0) then
            costav=costav/3.0
          end if

          !---- Statistic  and average values for previous temperature ----!
          do j=1,vs%nconf
             if(dead(j)) cycle
             if (naj(j) == 0) then
                naj(j)=1
                cost(j)=cost1(j)
             end if
             paj(j)=100.0*real(naj(j))/rav

             cost(j)=cost(j)/naj(j)

             do i=1,vs%npar
                raver(i,j)=raver(i,j)/naver
                sigp(i,j)=sqrt(abs(sigp(i,j)/naver-raver(i,j)*raver(i,j)))
             end do
             stepav=sum(vs%stp(:,j))/real(vs%npar)

             !---- Writing partial results ----!
             strings=" "
             write(unit=strings,fmt="(a,i4,a,f8.3,a,i2,a,f8.3,2(a,f12.3))")  &
             "     Conf:",j,"  (%Acc):",paj(j),"  <Step(",j,")>:",stepav, &
             "  <"//trim(c%Cost_function_name)//">:",cost(j),"  -> Current Cost:", cost1(j)
             call mess(strings)
             write(unit=ipr,fmt="(a)") strings
             write(unit=ipr,fmt="(a,i10 )")  "     Num-Cost-Evaluations: ",neval

            !Apply test of convergence and suppress bad configurations

            if (paj(j) <= c%accept ) then
                 dead(j)=.true. !Convergence criterium
                 write(unit=strings,fmt="(a,i3,a)") " => Configuration #",j," converged => dead in the algorithm!"
                 call mess(strings)
                 write(unit=ipr,fmt="(a)") strings
            end if

            if(abs(costav(j)-cost(j)) < 0.20 .and. modulo(ntp,4) == 0 .and. stepav < half_init_avstp ) then
                 dead(j)=.true. !Convergence criterium
                 write(unit=strings,fmt="(a,i3,a)") " => Configuration #",j," do not change anymore => dead in the algorithm!"
                 call mess(strings)
                 write(unit=ipr,fmt="(a)") strings
            end if

            do k=j+1,vs%nconf
              if(dead(k)) cycle
            !  sumdel = sum(abs(vs%state(:,j)-vs%state(:,k)))
              sumdel = sum(abs(raver(:,j)- raver(:,k)))
              sumsig= min(0.2, 0.5*sum(sigp(:,j)+sigp(:,k)))
              if( sumdel < sumsig )  then
                 dead(k) = .true.
                 write(unit=strings,fmt="(2(a,i3))") " => Configuration #",k," dead, because it is equal to Configuration #",j
                 call mess(strings)
                 write(unit=ipr,fmt="(a)") strings
              end if
            end do

          end do !j=1,vs%nconf

          if (present(filesav)) then
             open(unit=i_conf,file=trim(filesav)//".anl",status="old",action="write", position="append")
             write(unit=i_conf,fmt="(i4,31f10.4)")  ntp,temp,costop,vs%config(1:vs%npar)
             call flush(i_conf)
             close(unit=i_conf)
          end if


          survive=0
          do j=1,vs%nconf
           if(dead(j)) cycle
            survive=survive+1
          end do
          if(survive == 0) then
            strings = " => Convergence reached, look the list of configurations"
            call mess(strings)
            write(unit=ipr,fmt="(a)") strings
            exit
          end if

          !---- Adjust STP so that approximately half of all evaluations are accepted ----!
          !---- (Corana's algorithm)

          if (c%nalgor == 0 .or. c%nalgor == 1) then
             do j=1,vs%nconf
                if(dead(j)) cycle
                do i=1,vs%npar
                   jj=0
                   if (vs%stp(i,j) > abs(0.01*c%accept*raver(i,j))) exit
                   jj=1
                end do
                if (jj == 1) exit
                do i = 1, vs%npar
                   plage=vs%high(i)-vs%low(i)
                   rati = real(nacp(i,j)) /real(naver)
                   if (rati > 0.6) then
                      vs%stp(i,j) = vs%stp(i,j)*(1.0 + cv(i)*(rati - 0.6)/0.4)
                   else if (rati < 0.4) then
                      vs%stp(i,j) = vs%stp(i,j)/(1.0 + cv(i)*((0.4 - rati)/0.4))
                   end if
                   if (vs%stp(i,j) > plage) then
                      vs%stp(i,j) = plage
                   end if
                end do
             end do !j=1,vs%nconf
          end if
          nacp(:,:) = 0
          if(modulo(ntp,4) == 0) then
             costav(:)=0.0
          else
             costav(:)=costav(:)+ cost(:)
          end if

       end do   !ntp=1,c%num_temps

       !---- Re-calculate the cost function for the best configuration ----!
       call Model_Funct(vs%config,Costop)

       messag=" "
       call mess(messag)
       messag=" => Best configuration found by Simulated Annealing (Sigma of the last Montecarlo Cycles):"
       call mess(messag)
       write(unit=ipr,fmt="(/,a,a,/)") "     ",trim(messag)
       strings=" "
       write(unit=strings,fmt="(a,f10.4,a)")" -> Best Solution Cost = ",costop," :: "
       call mess(strings)
       write(unit=ipr,fmt="(/,a)") strings
       messag=" -> Configuration parameters :"
       call mess(messag)
       write(unit=ipr,fmt="(a,/)") trim(messag)
       do i=1,vs%npar
         write(unit=strings,fmt="(i6,a,2F16.5)")  i,vs%nampar(i), vs%config(i),sigp(i,jopt)
         write(unit=ipr,fmt="(a)") strings
         call mess(strings)
       end do

       return
   End Subroutine SimAnneal_MultiConf

    !!----
    !!---- Subroutine Write_SimAnn_Cond(ipr,c)
    !!----    integer,                     intent(in)  :: ipr !Logical unit for writing
    !!----    type(SimAnn_Conditions_type),intent(in)  :: c   !SAN Conditions
    !!----
    !!---- Subroutine for Writing in unit=ipr the SimAnn_Conditions_type
    !!---- variable "c"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Write_SimAnn_Cond(ipr,c)
       !---- Arguments ----!
       integer,                     intent(in)  :: ipr
       type(SimAnn_Conditions_type),intent(in)  :: c

       !--- Local Variables ---!


       write(unit=ipr,fmt="(/,a)") "     ===================================="
       write(unit=ipr,fmt="(a)")   "     =  SIMULATED ANNEALING CONDITIONS  ="
       write(unit=ipr,fmt="(a,/)") "     ===================================="

       write(unit=ipr,fmt="(a)")        " =>               Cost Function name: "//trim(c%Cost_function_name)
       write(unit=ipr,fmt="(a,f8.3)")   " =>              Initial Temperature: ",c%T_ini
       write(unit=ipr,fmt="(a,f8.3)")   " =>                 Annealing Factor: ",c%anneal
       write(unit=ipr,fmt="(a,i4  )")   " =>   Maximum number of Temperatures: ",c%num_temps
       if(c%nalgor == 0) then
          write(unit=ipr,fmt="(a)")     " => Corana's variable step algorithm  "
       else if(c%nalgor == 1) then
          write(unit=ipr,fmt="(a)")     " => Corana's variable step algorithm  "
       else
          write(unit=ipr,fmt="(a)")     " => Const. step conventional algorithm "
       end if
       write(unit=ipr,fmt="(a,i4  )")   " => Number of Montecarlo cycles/Temp: ",c%nm_cycl
       write(unit=ipr,fmt="(a,i4  )")   " =>  Number of thermalization cycles: ",c%num_therm
       write(unit=ipr,fmt="(a,i8)")     " =>Number of paralell configurations: ",c%num_conf
       if(c%seed == 0) then
          write(unit=ipr,fmt="(a)")     " => Random Seed selected from system clock "
       else
          write(unit=ipr,fmt="(a,i8)")  " =>    Initial Seed selected by user: ",c%seed
       end if
       if(c%initconfig == 0) then
          write(unit=ipr,fmt="(a)")     " =>    Initial random configuration "
       else
          write(unit=ipr,fmt="(a)")     " =>    Initial configuration as given  "
       end if
       return
    End Subroutine Write_SimAnn_Cond


    !!----
    !!---- Subroutine Write_SimAnn_MStateV(ipr,vs,text)
    !!----    integer,                     intent(in) :: ipr   !Logical unit for writing
    !!----    type(MultiState_Vector_Type),intent(in) :: vs    !State vector
    !!----    character(len=*),            intent(in) :: text
    !!----
    !!---- Subroutine for Writing in unit=ipr the SimAnn_Conditions_type
    !!---- variable "c"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Write_SimAnn_MStateV(ipr,vs,text)
       !---- Arguments ----!
       integer,                     intent(in) :: ipr
       type(MultiState_Vector_Type),intent(in) :: vs
       character(len=*),            intent(in) :: text

       !--- Local Variables ---!
       integer :: i,j
       character(len=30) :: forma1,forma2
       character(len=15) :: namep

       forma1="(a,   (a,i2))"
       write(unit=forma1(4:6),fmt="(i3)") vs%nconf
       forma2="(i6,a15,2f12.5,2i6,   f12.5))"
       write(unit=forma2(20:22),fmt="(i3)") 2*vs%nconf+1

       write(unit=ipr,fmt="(/,a)") "     ========================================================="
       write(unit=ipr,fmt="(  a)") "     => SIMULATED ANNEALING MULTI-STATE VECTOR: "//trim(text)
       write(unit=ipr,fmt="(a,/)") "     ========================================================="

       write(unit=ipr,fmt=forma1) "   Num           Name     Low_lim    High_Lim  BCond Code    BestConf   ",&
                                  ("Value&Step Conf#",j, j=1,vs%nconf)
       do i=1,vs%npar
         namep=adjustr(vs%nampar(i))
         write(unit=ipr,fmt=forma2) i,namep,vs%low(i),vs%high(i),vs%bound(i),vs%code(i),vs%config(i),&
                                   (vs%state(i,j),vs%stp(i,j),j=1,vs%nconf)
       end do

       forma1="(   (a,i2,a,f10.4))"
       write(unit=forma1(2:4),fmt="(i3)") vs%nconf
       write(unit=ipr,fmt="(/,a,/)") "  ==> Cost Function values of the different configurations"
       write(unit=ipr,fmt=forma1)   ("      Cost#",j,":",vs%Cost(j),j=1,vs%nconf)

       return
    End Subroutine Write_SimAnn_MStateV

    !!----
    !!---- Subroutine Write_SimAnn_StateV(ipr,vs)
    !!----    integer,                intent(in)  :: ipr   !Logical unit for writing
    !!----    type(State_Vector_Type),intent(in)  :: vs    !State vector
    !!----
    !!---- Subroutine for Writing in unit=ipr the SimAnn_Conditions_type
    !!---- variable "c"
    !!----
    !!---- Update: April - 2005
    !!

    Subroutine Write_SimAnn_StateV(ipr,vs,text)
       !---- Arguments ----!
       integer,                intent(in)  :: ipr
       type(State_Vector_Type),intent(in)  :: vs
       character(len=*),       intent(in)  :: text

       !--- Local Variables ---!
       integer :: i
       character(len=30) :: forma
       character(len=15)  :: namep
       forma="(i6,a8,3f12.5,2(i7,f12.5))"

       write(unit=ipr,fmt="(/,a)") "     ================================================="
       write(unit=ipr,fmt="(  a)") "     => SIMULATED ANNEALING STATE VECTOR: "//trim(text)
       write(unit=ipr,fmt="(a,/)") "     ================================================="
       write(unit=ipr,fmt="(a)") "   Num    Name       Value     Low_lim    High_Lim    BCond      Step    Code   BestConf"
       do i=1,vs%npar
         namep=adjustr(vs%nampar(i))
         write(unit=ipr,fmt=forma)  i,namep,vs%state(i),vs%low(i),vs%high(i),vs%bound(i),vs%stp(i),vs%code(i),vs%config(i)
       end do
       return
    End Subroutine Write_SimAnn_StateV


 End Module Optimization_SAN
