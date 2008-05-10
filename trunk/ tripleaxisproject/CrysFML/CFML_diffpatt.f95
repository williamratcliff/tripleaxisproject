!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: DIFFRACTION_PATTERNS_MOD
!!----   INFO: Diffraction Patterns Information
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----            January - 2004  Created by JRC
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    DIFFRACTION_PATTERN_TYPE
!!----    ERR_DIFFPATT
!!----    ERR_MESS_DIFFPATT
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ALLOCATE_DIFFRACTION_PATTERN
!!----       INIT_ERR_DIFFPATT
!!----       PURGE_DIFFRACTION_PATTERN
!!----       READ_BACKGROUND_FILE
!!----       READ_PATTERN
!!--++       READ_PATTERN_D1A_D2B           [Private]
!!--++       READ_PATTERN_D1A_D2B_OLD       [Private]
!!--++       READ_PATTERN_D1B_D20           [Private]
!!--++       READ_PATTERN_DMC               [Private]
!!--++       READ_PATTERN_FREE              [Private]
!!--++       READ_PATTERN_G41               [Private]
!!--++       READ_PATTERN_GSAS              [Private]
!!--++       READ_PATTERN_ISIS_M            [Private]
!!--++       READ_PATTERN_MULT              [Overloaded]
!!--++       READ_PATTERN_NLS               [Private]
!!--++       READ_PATTERN_ONE               [Overloaded]
!!--++       READ_PATTERN_PANALYTICAL_CSV   [Private]
!!--++       READ_PATTERN_PANALYTICAL_JCP   [Private]
!!--++       READ_PATTERN_PANALYTICAL_UDF   [Private]
!!--++       READ_PATTERN_PANALYTICAL_XRDML [Private]
!!--++       READ_PATTERN_SOCABIM           [Private]
!!--++       READ_PATTERN_TIME_VARIABLE     [Private]
!!--++       READ_PATTERN_XYSIGMA           [Private]
!!--++       SET_BACKGROUND_INTER           [Private]
!!--++       SET_BACKGROUND_POLY            [Private]
!!----
!!
 Module Diffraction_patterns_mod
    !---- Use Modules ----!
    Use Math_gen,         only : sp, spline, splint, locate
    use String_Utilities, only : FindFmt,  Init_FindFmt , ierr_fmt, &
                                 get_logunit, u_case, getword !, mess_findfmt

    implicit none

    private

    !---- List of public subroutines ----!
    public ::  Init_Err_DiffPatt, Read_Background_File, Read_Pattern, &
               Purge_Diffraction_Pattern, Allocate_Diffraction_Pattern

    !---- List of private subroutines ----!
    private :: Read_Pattern_D1A_D2B, Read_Pattern_D1A_D2B_Old, Read_Pattern_D1B_D20,       &
               Read_Pattern_Dmc, Read_Pattern_Free, Read_Pattern_G41, Read_Pattern_Gsas,   &
               Read_Pattern_Isis_M, Read_Pattern_Mult, Read_Pattern_Nls, Read_Pattern_One, &
               Read_Pattern_Panalytical_Csv, Read_Pattern_Panalytical_Jcp,                 &
               Read_Pattern_Panalytical_Udf, Read_Pattern_Panalytical_Xrdml,               &
               Read_Pattern_Socabim, Read_Pattern_Time_Variable, Read_Pattern_Xysigma,     &
               Set_Background_Inter, Set_Background_Poly

    !---- Definitions ----!

    !!----
    !!---- TYPE :: DIFFRACTION_PATTERN_TYPE
    !!--..
    !!---- Type, public :: Diffraction_Pattern_Type
    !!----    character(len=80)                           :: Title         !Identification of the pattern
    !!----    character(len=20)                           :: diff_kind     !type of radiation
    !!----    character(len=20)                           :: scat_var      !x-space: 2theta, TOF, Q, S, d, etc
    !!----    character(len=20)                           :: instr         !file type
    !!----    real(kind=sp)                               :: xmin
    !!----    real(kind=sp)                               :: xmax
    !!----    real(kind=sp)                               :: ymin
    !!----    real(kind=sp)                               :: ymax
    !!----    real(kind=sp)                               :: scal
    !!----    real(kind=sp)                               :: step
    !!----    integer                                     :: npts          ! Number of points
    !!----    logical                                     :: ct_step       ! Constant step
    !!----    logical                                     :: gy,gycalc,&
    !!----                                                   gbgr,gsigma   !logicals for graphics
    !!----
    !!----    logical                                     :: al_x,al_y,&
    !!----                                                   al_ycalc, &   !logicals for allocation
    !!----                                                   al_bgr,   &
    !!----                                                   al_sigma, &
    !!----                                                   al_istat
    !!----
    !!----    real(kind=sp), dimension (3)                :: conv          ! Wavelengths or Dtt1, Dtt2 for converting to Q,d, etc
    !!----    real(kind=sp), dimension (:), allocatable   :: x             ! Scattering variable (2theta...)
    !!----    real(kind=sp), dimension (:), allocatable   :: y             ! Experimental intensity
    !!----    real(kind=sp), dimension (:), allocatable   :: sigma         ! observations VARIANCE (it is the square of sigma!)
    !!----    integer,       dimension (:), allocatable   :: istat         ! Information about the point "i"
    !!----    real(kind=sp), dimension (:), allocatable   :: ycalc         ! Calculated intensity
    !!----    real(kind=sp), dimension (:), allocatable   :: bgr           ! Background
    !!----
    !!----
    !!----
    !!----
    !!----
    !!----
    !!---- End Type Diffraction_Pattern_Type
    !!----
    !!----    Definition for Diffraction Pattern
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Diffraction_Pattern_Type
       character(len=80)                           :: Title         !Identification of the pattern
       character(len=20)                           :: diff_kind     !type of radiation
       character(len=20)                           :: scat_var      !x-space: 2theta, TOF, Q, S, d, etc
       character(len=20)                           :: instr         !file type
       real(kind=sp)                               :: xmin
       real(kind=sp)                               :: xmax
       real(kind=sp)                               :: ymin
       real(kind=sp)                               :: ymax
       real(kind=sp)                               :: scal
       real(kind=sp)                               :: step
       integer                                     :: npts          ! Number of points
       logical                                     :: ct_step       ! Constant step
       logical                                     :: gy,gycalc,&
                                                      gbgr,gsigma   !logicals for graphics

       logical                                     :: al_x,al_y,&
                                                      al_ycalc, &   !logicals for allocation
                                                      al_bgr,   &
                                                      al_sigma, &
                                                      al_istat

       real(kind=sp), dimension (3)                :: conv          ! Wavelengths or Dtt1, Dtt2 for converting to Q,d, etc
       real(kind=sp), dimension (:), allocatable   :: x             ! Scattering variable (2theta...)
       real(kind=sp), dimension (:), allocatable   :: y             ! Experimental intensity
       real(kind=sp), dimension (:), allocatable   :: sigma         ! observations VARIANCE (it is the square of sigma!)
       integer,       dimension (:), allocatable   :: istat         ! Information about the point "i"
       real(kind=sp), dimension (:), allocatable   :: ycalc         ! Calculated intensity
       real(kind=sp), dimension (:), allocatable   :: bgr           ! Background
    End Type Diffraction_Pattern_Type

    !!----
    !!---- ERR_DIFFPATT
    !!----    logical, public :: Err_Diffpatt
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_diffpatt=.false.

    !!----
    !!---- ERR_MESS_DIFFPATT
    !!----    character(len=150), public :: Err_Mess_Diffpatt
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_diffpatt=" "


    !---- Interfaces - Overlapp ----!
    Interface Read_Pattern
       Module procedure Read_Pattern_Mult
       Module procedure Read_Pattern_One
    End Interface

 Contains

    !---------------------!
    !---- Subroutines ----!
    !---------------------!
    !!----
    !!---- Subroutine Allocate_Diffraction_Pattern(pat,n)
    !!----    type(Diffraction_Pattern_Type), intent (in out) :: pat
    !!----    Integer,                        intent (in)     :: n
    !!----
    !!----    Allocate the object pat of type Diffraction_Pattern_Type
    !!----
    !!---- Update: December - 2005
    !!
    Subroutine Allocate_Diffraction_Pattern(pat,npts)
       type(Diffraction_Pattern_Type), intent (in out) :: pat
       Integer, optional,              intent (in)     :: npts
       !----- Local variables ----!
       integer :: n

       if(present(npts)) then
         pat%npts=npts
         n=npts
       else
         n=pat%npts
       end if

       if(n <= 0) then
         err_diffpatt=.true.
         err_mess_diffpatt=" Attempt to allocate Diffraction_Pattern with 0-dimension "
         return
       end if

       if(allocated(pat%y) ) deallocate(pat%y)
       allocate(pat%y(n))
       pat%y=0.0
       pat%gy=.true.
       pat%al_y=.true.

       if(allocated(pat%ycalc) ) deallocate(pat%ycalc)
       allocate(pat%ycalc(n))
       pat%ycalc=0.0
       pat%gycalc=.true.
       pat%al_ycalc=.true.

       if(allocated(pat%bgr) ) deallocate(pat%bgr)
       allocate(pat%bgr(n))
       pat%bgr=0.0
       pat%gbgr=.true.
       pat%al_bgr=.true.

       if(allocated(pat%x) ) deallocate(pat%x)
       allocate(pat%x(n))
       pat%x=0.0
       pat%al_x=.true.

       if(allocated(pat%sigma) ) deallocate(pat%sigma)
       allocate(pat%sigma(n))
       pat%sigma=0.0
       pat%gsigma=.true.
       pat%al_sigma=.true.

       if(allocated(pat%istat) ) deallocate(pat%istat)
       allocate(pat%istat(n))
       pat%istat=1
       pat%al_istat=.true.

       return
    End Subroutine Allocate_Diffraction_Pattern

    !!----
    !!---- Subroutine Init_Err_DiffPatt()
    !!----
    !!----    Initialize the errors flags in DiffPatt
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_DiffPatt()

       err_diffpatt=.false.
       err_mess_diffpatt=" "

       return
    End Subroutine Init_Err_Diffpatt

    !!----
    !!---- Subroutine Purge_Diffraction_Pattern(pat,MODE)
    !!----    type(Diffraction_Pattern_Type), intent (in out) :: pat
    !!----    Character(len=*),               intent (in)     :: MODE
    !!----
    !!----    De-Allocate components of the object "pat", of type Diffraction_Pattern_Type
    !!----    depending on the value of the MODE string. At present the following MODE
    !!----    values are available:
    !!----      "DATA " -> x,y remain allocated                  (purge sigma,ycalc,bgr,istat)
    !!----      "DATAS" -> x,y,sigma remain allocated            (purge ycalc,bgr,istat)
    !!----      "RIETV" -> x,y,sigma,ycalc,bgr remain allocated  (purge istat)
    !!----      "GRAPH" -> x,y,sigma,istat remain allocated      (purge ycalc, bgr)
    !!----      "PRF  " -> x,y,ycalc,bgr,istat, remain allocated (purge sigma)
    !!----
    !!----
    !!---- Update: December - 2005
    !!
    Subroutine Purge_Diffraction_Pattern(pat,MODE)
       type(Diffraction_Pattern_Type), intent (in out) :: pat
       Character(len=*),               intent (in)     :: MODE

         !Modes: "DATA " -> only x,y are allocated
         !       "DATAS" -> x,y,sigma are allocated
         !       "RIETV" -> x,y,sigma,ycalc,bgr are allocated
         !       "GRAPH" -> x,y,sigma,istat are allocated
         !       "PRF  " -> x,y,ycalc,bgr,istat, are allocated


       Select Case (u_case(MODE))

         Case("DATA ")    !Mode: "DATA " -> only x,y remain allocated

            if(allocated(pat%ycalc)) deallocate(pat%ycalc)
            pat%gycalc=.false.
            pat%al_ycalc=.false.

            if(allocated(pat%bgr)) deallocate(pat%bgr)
            pat%gbgr=.false.
            pat%al_bgr=.false.

            if(allocated(pat%sigma)) deallocate(pat%sigma)
            pat%gsigma=.false.
            pat%al_sigma=.false.

            if(allocated(pat%istat)) deallocate(pat%istat)

         Case("DATAS")    !Mode: "DATAS" -> only x,y, sigma remain allocated

            if(allocated(pat%ycalc)) deallocate(pat%ycalc)
            pat%gycalc=.false.
            pat%al_ycalc=.false.

            if(allocated(pat%bgr)) deallocate(pat%bgr)
            pat%gbgr=.false.
            pat%al_bgr=.false.

            if(allocated(pat%istat)) deallocate(pat%istat)

         Case("RIETV")   !Mode: "RIETV" -> x,y,sigma,ycalc,bgr remain allocated

            if(allocated(pat%istat)) deallocate(pat%istat)

         Case("GRAPH")   !Mode: "GRAPH" -> x,y,sigma,istat remain allocated

            if(allocated(pat%ycalc)) deallocate(pat%ycalc)
            pat%gycalc=.false.
            pat%al_ycalc=.false.

            if(allocated(pat%bgr)) deallocate(pat%bgr)
            pat%gbgr=.false.
            pat%al_bgr=.false.

         Case("PRF  ")

            if(allocated(pat%sigma)) deallocate(pat%sigma)
            pat%gsigma=.false.
            pat%al_sigma=.false.

       End Select
       return
    End Subroutine Purge_Diffraction_Pattern

    !!----
    !!---- Subroutine Read_Backgound_File(bck_file, bck_mode, dif_pat)
    !!----    character (len=*),               intent(in   )    :: bck_file
    !!----    character (len=*),               intent(in   )    :: bck_mode
    !!----    type (diffraction_pattern_type), intent(out)      :: dif_Pat
    !!----
    !!----    Read background from a file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Background_File( Bck_File, Bck_Mode, Dif_Pat)
       !---- Arguments ----!
       character (len=*),               intent(in   )    :: bck_file
       character (len=*),               intent(in   )    :: bck_mode
       type (diffraction_pattern_type), intent(in out)   :: dif_pat

       !---- local variables ----!
       logical                               :: esta
       character (len=132)                   :: line
       integer                               :: bck_points
       integer                               :: i,j,i_bck
       integer                               :: ier, alloc_error
       real , dimension (:), allocatable     :: bck_v
       real , dimension (:), allocatable     :: bck_p

       call init_err_diffpatt()

       inquire(file=bck_file, exist =esta)
       if (.not. esta) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" The file "//trim(bck_file)//" doesn't exist"
          return
       else
          call get_logunit(i_bck)
          open(unit=i_bck,file=trim(bck_file),status="old",action="read",position="rewind",iostat=ier)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error opening the file: "//trim(bck_file)
             return
          end if
       end if

       i=0
       do
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit
          if (len_trim(line) == 0) cycle
          if (index(line,"!") /= 0) cycle
          i=i+1
       end do
       bck_points=i
       rewind(unit = i_bck)

       if (allocated(bck_v)) deallocate(bck_v)
       allocate(bck_v(bck_points+1),stat= alloc_error)
       if (alloc_error /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Allocation error reading background points"
          return
       end if

       if (allocated(bck_p)) deallocate(bck_p)
       allocate(bck_p(bck_points+1), stat= alloc_error)
       if (alloc_error /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Allocation error reading background points"
          return
       end if

       read(unit=i_bck,fmt="(a)",iostat=ier) line
       read(unit=i_bck,fmt="(a)",iostat=ier) line

       do j=1, bck_points
          read(unit=i_bck,fmt="(a)",iostat=ier) line
          if (ier /= 0) exit
          if (len_trim(line) == 0) cycle
          read(unit=line, fmt=*, iostat=ier)  bck_p(j), bck_v(j)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading background file!"
             return
          end if
       end do

       select case (u_case(bck_mode(1:3)))
          case ("POL") ! Polynomial
             call set_background_poly (dif_pat,50.0, bck_p,bck_points )

          case ("INT") ! Interpolation
             call  set_background_inter (dif_pat, bck_v,bck_p, bck_points )

          case default
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Not a valid mode"
             return
       end select

       close(unit=i_bck,iostat=ier)
       if (ier/=0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Problems closing data file"
          return
       end if

       return
    End Subroutine Read_Background_File

    !!----
    !!---- Subroutine Read_Pattern(Filename, Dif_Pat, Mode)
    !!----                   or   (Filename, Dif_Pat, NumPat, Mode)
    !!--<<    character(len=*),                              intent (in) :: Filename
    !!----    type (diffraction_pattern_type),               intent (out):: Dif_Pat
    !!----    character(len=*), optional,                    intent (in) :: mode
    !!----
    !!----    character(len=*),                              intent (in) :: Filename
    !!----    type (diffraction_pattern_type), dimension(:), intent (out):: Dif_Pat
    !!----    integer,                                       intent (out):: numpat
    !!----    character(len=*), optional,                    intent (in) :: mode
    !!-->>
    !!----    Read one pattern from a Filename
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Read_Pattern_D1A_D2B(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    (PRIVATE)
    !!--++    Read a pattern for D1A, D2B
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_D1A_D2B(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=132)                           :: txt1
       integer                                      :: i, nlines, j, no, ier
       integer, dimension(:), allocatable           :: iww
       real                                         :: rmoni, rmoniold, cnorm

       call init_err_diffpatt()

       read(unit=i_dat,fmt="(a)",iostat=ier) txt1
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if
       pat%title=txt1

       read(unit=i_dat,fmt="(tr16,F8.3)",iostat=ier) pat%step
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(F8.3)",iostat=ier)pat%xmin
       if (ier /= 0)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(2f8.0)",iostat=ier) rmoni,rmoniold
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       if (rmoniold < 1.0) then
          cnorm=1.00
          rmoniold=rmoni
       else
          cnorm=rmoni/rmoniold
       end if

       nlines = nint(18.0/pat%step)
       pat%npts  = 10*nlines
       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww) ) deallocate(iww)
       allocate(iww(pat%npts))

       j=0
       do i=1,nlines
          read(unit=i_dat,fmt="(10(i2,f6.0))",iostat=ier)(iww(j+no),pat%y(j+no),no=1,10)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
             return
          end if
          if(abs(pat%y(j+1)+1000.0) < 1.0e-03) exit
          j = j+10
       end do
       j=j-10
       pat%npts=j
       pat%xmax = pat%xmin+(pat%npts-1)*pat%step
       do i=1,pat%npts
          if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) == 0) iww(i) = 1
          pat%sigma(i) = cnorm*pat%y(i)/real(iww(i))
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_D1A_D2B

    !!--++
    !!--++ Subroutine Read_Pattern_D1A_D2B_OLD(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for D1A, D2B (Old Format)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_D1A_D2B_OLD(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       integer                                      :: ier,i
       integer, dimension(:), allocatable           :: iww

       call init_err_diffpatt()

       read(unit=i_dat,fmt=*,iostat=ier)pat%xmin,pat%step,pat%xmax
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if
       pat%title=" No title: data format -> old D1A"
       pat%npts = (pat%xmax-pat%xmin)/pat%step+1.5
       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww) ) deallocate(iww)
       allocate(iww(pat%npts))

       read(unit=i_dat,fmt="(10(i2,f6.0))",iostat=ier)(iww(i),pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       do i=1,pat%npts
          if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) == 0) iww(i) = 1
          pat%sigma(i) = pat%y(i)/real(iww(i))
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End Subroutine Read_Pattern_D1A_D2B_Old

    !!--++
    !!--++ Subroutine Read_Pattern_D1B_D20(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for D1B or D20
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_D1B_D20(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=132)                           :: line
       integer                                      :: i,j,npunt,no,ier
       integer, dimension(:), allocatable           :: iww

       call init_err_diffpatt()

       do i=1,3
          read(unit=i_dat,fmt=*, iostat=ier)line
          if (ier /= 0 )then
             Err_diffpatt=.true.
             Err_mess_diffpatt=" Error in  Intensity file, check your instr parameter!"
             return
          end if
          if( i == 1) pat%title=line
       end do

       read(unit=i_dat,fmt="(tr23,f8.3,tr45, f8.3)  ",iostat=ier) pat%xmin,pat%step
       if (ier /= 0 ) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(i4)",iostat=ier) pat%npts
       if (ier /= 0 )then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww)) deallocate(iww)
       allocate(iww(pat%npts))

       j = 0
       npunt = nint(REAL(pat%npts)/10.0)
       pat%xmax = pat%xmin+(pat%npts-1)*pat%step
       do i=1,npunt
          read(unit=i_dat,fmt="(10(i2,f8.0))",iostat=ier) (iww(j+no),pat%y(j+no),no=1,10)
          if (ier /= 0 )then
             Err_diffpatt=.true.
             Err_mess_diffpatt=" Error in  Intensity file, check your instr parameter!"
             return
          end if
          j = j+10
       end do

       read(unit=i_dat,fmt=*,iostat=ier)line
       if (ier /= 0 )then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       do i=1,pat%npts
          if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
          if (iww(i) <= 0) iww(i) = 1
          pat%sigma(i) = pat%y(i)/REAL(iww(i))
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_D1B_D20

    !!--++
    !!--++ Subroutine Read_Pattern_Dmc(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for DMC
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Dmc(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=132)                           :: txt1
       integer                                      :: ier, i

       call init_err_diffpatt()

       read(unit=i_dat,fmt="(A)",iostat=ier)txt1
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if
       pat%title=txt1

       read(unit=i_dat,fmt="(A)",iostat=ier)txt1
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt=*,iostat=ier) pat%xmin,pat%step,pat%xmax
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if
       pat%npts = (pat%xmax - pat%xmin)/pat%step + 1.005
       if (pat%npts < 20)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       read(unit=i_dat,fmt="(10f8.0)",iostat=ier)(pat%y(i),i=1,pat%npts)
       if (ier > 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(10f8.0)",iostat=ier)(pat%sigma(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       do i=1,pat%npts
         pat%sigma(i) = pat%sigma(i)*pat%sigma(i)
         pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_Dmc

    !!--++
    !!--++ Subroutine Read_Pattern_Free(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Free
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Free(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       integer                                      :: i,no,ier,inum
       character(len=132)                           :: aline
       logical                                      :: title_given

       call init_err_diffpatt()
       title_given=.false.
       no=0
       do
          read(unit=i_dat,fmt="(a)",iostat=ier) aline
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" End of file *.dat"
             return
          else
             if(.not. title_given) then
               Pat%title=aline(2:)
               title_given=.true.
             end if
             if (aline(1:1) == "!" .or. aline(1:1) == "#") cycle
             if (aline(1:4) == "BANK") then
                read(unit=aline(5:41),fmt=*) inum,pat%npts
                read(unit=aline(47:90),fmt=*) pat%xmin,pat%step
                pat%xmax=pat%xmin+(pat%npts-1)*pat%step
             else
                read(unit=aline,fmt=*,iostat=ier)pat%xmin,pat%step,pat%xmax
                if (ier /= 0) then
                   no=no+1
                   if (no > 7)then
                      Err_diffpatt=.true.
                      Err_Mess_DiffPatt=" Error on Intensity file, check your instr parameter "
                      return
                   else
                      cycle
                   end if
                end if
                pat%npts = (pat%xmax-pat%xmin)/pat%step+1.5
             end if
          end if
          exit
       end do

       if (pat%npts <= 0)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       read(unit=i_dat,fmt=*,iostat=ier)(pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter"
          return
       end if
       do i=1,pat%npts
          pat%sigma(i) = pat%y(i)
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_Free

    !!--++
    !!--++ Subroutine Read_Pattern_G41(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for G41
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_G41(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=132)                           :: txt1, txt2, txt3
       integer                                      :: i, ier, ivari
       real                                         :: tsamp, cnorm
       real(kind=sp)                                :: treg, rmon1, rmon2


       call init_err_diffpatt()

       read(unit=i_dat,fmt="(A)",iostat=ier)txt1                  !1
       pat%title=txt1
       if (ier /= 0 ) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier)txt2                  !2
       if (ier /= 0 ) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier)txt3                  !3
       if (ier /= 0 ) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       read(unit=i_dat,fmt="(I6,tr1,2F10.3,i5,2f10.1)",iostat=ier)  pat%npts,tsamp,treg,ivari,rmon1,rmon2
       if (ier /= 0 )then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)


       read(unit=i_dat,fmt="(3F10.0)",iostat=ier)pat%xmin,pat%step,pat%xmax              !5
       if (ier /= 0 ) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       read(unit=i_dat,fmt=*,iostat=ier)(pat%y(i),i=1, pat%npts)
       if (ier /= 0 )then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
          return
       end if

       if (ivari /= 0) then          !IVARI
          read(unit=i_dat,fmt=*,iostat=ier)(pat%sigma(i),i=1, pat%npts)
          if (ier /= 0 ) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter! "
             return
          end if
          cnorm=0.0
          do i=1,pat%npts
             IF (pat%y(i) < 0.0001) pat%y(i)=0.0001
             pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
             IF (pat%sigma(i) < 0.000001) pat%sigma(i)=1.0
             pat%x(i)= pat%xmin+(i-1)*pat%step
             cnorm=cnorm+pat%sigma(i)/pat%y(i)
          end do
          cnorm=cnorm/REAL(pat%npts)
       else                         !ivari
          if (rmon1 > 1.0 .and. rmon2 > 1.0) then
             cnorm=rmon1/rmon2
          else
             cnorm=1.0
          end if
          do i=1,pat%npts
             pat%sigma(i)=pat%y(i)*cnorm
             pat%x(i)= pat%xmin+(i-1)*pat%step
          end do
       end if                        !IVARI
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_G41

    !!--++
    !!--++ Subroutine Read_Pattern_Gsas(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for GSAS
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Gsas(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       logical                                      :: previous, bank_missed
       logical, save                                :: keep_open=.false.
       character (len=80)                           :: line
       character (len=8 )                           :: bintyp,datyp
       integer                                      :: items,i
       integer                                      :: ibank,nchan,nrec, ier !, jobtyp
       integer,          dimension(:), allocatable  :: iww
       integer,          dimension(40)              :: pointi, pointf
       real(kind=sp),    dimension(4)               :: bcoef
       real(kind=sp)                                :: divi
       real                                         :: cnorm

       call init_err_diffpatt()

       divi=100.0
       if ( .not. keep_open) then
          bank_missed=.true.

          do i=1,7
             read(unit=i_dat,fmt="(a)") line
             if(i == 1) pat%title=line
             if (line(1:4) == "BANK") then
                bank_missed=.false.
                exit
             end if
          end do

          if (bank_missed) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" => Error in the input GSAS-file: BANK not found!"
             return
          end if
       else
          read(unit=i_dat,fmt="(a)") line
       end if

       items=0
       previous=.false.
       do i=5,80
          if (line(i:i) /= " ") then
             if (.not. previous) then
                items=items+1
                pointi(items)=i
                previous=.true.
             end if
          else
             if (items > 0 .and. previous) pointf(items)=i-1
             previous=.false.
          end if
       end do
       IF (items > 0) read(unit=line(pointi(1):pointf(1)),fmt=*) ibank
       IF (items > 1) read(unit=line(pointi(2):pointf(2)),fmt=*) nchan
       IF (items > 2) read(unit=line(pointi(3):pointf(3)),fmt=*) nrec
       IF (items > 3) read(unit=line(pointi(4):pointf(4)),fmt="(a)") bintyp
       IF (items > 4) read(unit=line(pointi(5):pointf(5)),fmt=*) bcoef(1)
       IF (items > 5) read(unit=line(pointi(6):pointf(6)),fmt=*) bcoef(2)
       IF (items > 6) read(unit=line(pointi(7):pointf(7)),fmt=*) bcoef(3)
       IF (items > 7) read(unit=line(pointi(8):pointf(8)),fmt=*) bcoef(4)
       datyp="STD"
       IF (items > 8) read(unit=line(pointi(9):pointf(9)),fmt="(a)") datyp
       divi=100.0
       pat%npts=nchan
       if (pat%npts <= 0)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(iww) ) deallocate(iww)
       allocate(iww(pat%npts))

       if (datyp == "STD") then
          pat%ct_step  = .true.
          if (bintyp == "CONST") then
             pat%xmin=bcoef(1)/divi !divide by 100 for CW
             pat%step=bcoef(2)/divi  !divide by 100 for CW
             pat%xmax=pat%xmin+(pat%npts-1)*pat%step
          else
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" => Only BINTYP=CONST is allowed for ESD data"
             return
          end if

          read(unit=i_dat,fmt="(10(i2,f6.0))", iostat=ier) (iww(i),pat%y(i),i=1,pat%npts)
          if (ier /= 0) then
             backspace (unit=i_dat)
          end if
          do i=1,pat%npts
             if (pat%y(i) <= 0.00001) pat%y(i) = 1.0
             if (iww(i) == 0) iww(i) = 1
             pat%sigma(i) = pat%y(i)/real(iww(i))
             pat%x(i)=pat%xmin+(i-1)*pat%step
          end do
          cnorm=1.0
          return

       else if(datyp == "ESD") then
          if (bintyp == "CONST") then
             pat%ct_step  = .true.
             pat%xmin=bcoef(1)/divi !divide by 100 for CW
             pat%step=bcoef(2)/divi  !divide by 100 for CW
             pat%xmax=pat%xmin+(pat%npts-1)*pat%step
             read(unit=i_dat,fmt="(10f8.0)",iostat=ier) (pat%y(i),pat%sigma(i),i=1,pat%npts)
             if (ier /= 0) then
                backspace (unit=i_dat)
             end if
             cnorm=0.0
             do i=1,pat%npts
                pat%x(i)=pat%xmin+(i-1)*pat%step
                pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
                cnorm=cnorm+pat%sigma(i)/max(pat%y(i),0.001)
             end do
             cnorm=cnorm/real(pat%npts)
          else
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" => Only BINTYP=CONST is allowed for ESD data"
             return
          end if

       else if(datyp == "ALT") then
          if (bintyp == "RALF") then
             pat%ct_step  = .false.
             read(unit=i_dat,fmt="(4(f8.0,f7.4,f5.4))",iostat=ier)(pat%x(i),pat%y(i),pat%sigma(i),i=1,pat%npts)
             if (ier /= 0) then
                backspace (unit=i_dat)
             end if
             pat%x=pat%x/32.0
             cnorm=0.0
             do i=1,pat%npts-1
                divi=pat%x(i+1)-pat%x(i)
                pat%y(i)=1000.0*pat%y(i)/divi
                pat%sigma(i)=1000.0*pat%sigma(i)/divi
                pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
                cnorm=cnorm+pat%sigma(i)/max(pat%y(i),0.001)
             end do
             cnorm=cnorm/real(pat%npts)
             pat%npts=pat%npts-1
             pat%xmin=bcoef(1)/32.0
             pat%step=bcoef(2)/32.0
             pat%xmax=pat%x(pat%npts)

          else if(bintyp == "CONST") then
             pat%ct_step  = .true.
             read(unit=i_dat,fmt="(4(f8.0,f7.4,f5.4))", iostat=ier)(pat%x(i),pat%y(i),pat%sigma(i),i=1,pat%npts)
             if (ier /= 0) then
                backspace (unit=i_dat)
             end if
             pat%x=pat%x/32.0
             cnorm=0.0
             do i=1,pat%npts
                pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
                cnorm=cnorm+pat%sigma(i)/max(pat%y(i),0.001)
             end do
             cnorm=cnorm/real(pat%npts)
             pat%xmin=bcoef(1)
             pat%step=bcoef(2)
             pat%xmax=pat%x(pat%npts)
          else
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" =>  Only BINTYP=RALF or CONST is allowed for ALT data"
          end if
       end if

       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       !---- Checking range and re-select the usable diffraction pattern

       return
    End Subroutine Read_Pattern_Gsas

    !!--++
    !!--++ Subroutine Read_Pattern_Isis_M(i_dat,Pat,NPat)
    !!--++    integer,                         intent(in    )  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(   out) :: Pat
    !!--++    integer,                         intent(in out) :: Npat
    !!--++
    !!--++    Read a pattern for ISIS
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Isis_m(i_dat,Pat,NPat)
       !---- Arguments ----!
       integer,                                                    intent(in    ) :: i_dat
       type (diffraction_pattern_type),  dimension(:), allocatable,intent(   out) :: pat
       integer,                                                    intent(in out) :: npat

       !---- Local Variables ----!
       real                                            :: fac_y
       real                                            :: cnorm
       real                                            :: sumavar
       integer                                         :: ntt, i, j, ier
       integer                                         :: n_pat      !index of current pattern
       integer, dimension(npat)                        :: npp        !number of points per pattern
       character(len=120)                              :: txt1
       character(len=132)                              :: aline
       real(kind=sp)                                   :: divi
       real(kind=sp), parameter                        :: eps1=1.0e-1
       logical                                         :: bankfound
       logical, save                                   :: ralf_type, title_given

       call init_err_diffpatt()

       fac_y=1000.0
       npp(:)=0
       n_pat=0
       bankfound=.false.
       title_given=.false.

       do
          read(unit=i_dat,fmt="(a)", iostat = ier) txt1
          if (ier /= 0) exit
          txt1=adjustl(txt1)
          if (txt1(1:4) == "BANK") then
             n_pat=n_pat+1
             npp(n_pat)=0
             bankfound=.true.
             cycle
          end if
          if (bankfound) npp(n_pat)=npp(n_pat)+1
       end do
       rewind(unit=i_dat)

       pat%ct_step = .false.

       if (npat <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       if (allocated(pat)) deallocate (pat)
       allocate (pat(npat))

       do n_pat=1,npat
          call Allocate_Diffraction_Pattern(pat(n_pat))
       end do

       do
          read(unit=i_dat,fmt="(a)") txt1
          if(.not. title_given) then
            Pat(1)%title=txt1
            title_given=.true.
          end if
          if (txt1(1:4) == "BANK") then
             IF (index(txt1,"RALF") /= 0) ralf_type =.true.
             exit
          end if
          i=index(txt1,"fac_y")
          if (i /= 0) then
             read(unit=txt1(i+5:),fmt=*) fac_y
          end if
       end do

       do n_pat=1,npat
          i=0
          ntt=0
          sumavar=0.0
          cnorm=0.0
          Pat(n_pat)%title=Pat(1)%title

          if (ralf_type) then
             do j=1,npp(n_pat)+1
                read(unit=i_dat,fmt="(a)",iostat=ier) aline
                if (ier /= 0)  exit
                if (aline(1:1) == "!" .or. aline(1:1) == "#") cycle
                if (aline(1:4) == "BANK") exit
                if (len_trim(aline)==0)exit
                i=i+1
                read(unit=aline,fmt=*,iostat=ier) pat(n_pat)%x(i),pat(n_pat)%y(i),pat(n_pat)%sigma(i)
                if (ier /= 0) then
                   Err_diffpatt=.true.
                   Err_Mess_DiffPatt=" Error reading an ISIS profile DATA file"
                   return
                end if

                if (abs(pat(n_pat)%x(i)) < eps1 .and. pat(n_pat)%y(i) < eps1 .and. pat(n_pat)%sigma(i) < eps1) exit
                pat(n_pat)%y(i)=pat(n_pat)%y(i)*fac_y
                pat(n_pat)%sigma(i)=pat(n_pat)%sigma(i)*fac_y
                pat(n_pat)%sigma(i)=pat(n_pat)%sigma(i)*pat(n_pat)%sigma(i)
                sumavar=sumavar+pat(n_pat)%sigma(i)
                if (pat(n_pat)%sigma(i) < eps1) pat(n_pat)%sigma(i) =fac_y
                if (pat(n_pat)%y(i) < eps1) then
                   pat(n_pat)%y(i)   = eps1
                   pat(n_pat)%sigma(i) = fac_y
                end if
                cnorm=cnorm+pat(n_pat)%sigma(i)/max(pat(n_pat)%y(i),0.001)
                if (i > 1) then
                   pat(n_pat)%step=pat(n_pat)%step+pat(n_pat)%x(i)-pat(n_pat)%x(i-1)
                   ntt=ntt+1
                end if
             end do
             do i=1,ntt
                divi=pat(n_pat)%x(i+1)-pat(n_pat)%x(i)
                pat(n_pat)%y(i)=pat(n_pat)%y(i)/divi
                pat(n_pat)%sigma(i)=pat(n_pat)%sigma(i)/divi/divi
             end do
             ntt=ntt-1

          else
             do j=1,npp(n_pat)
                read(unit=i_dat,fmt="(a)",iostat=ier) aline
                if (ier /= 0) exit
                if (aline(1:1) == "!" .or. aline(1:1) == "#") cycle
                if (aline(1:4) == "BANK") exit
                i=i+1
                read(unit=aline,fmt=*,iostat=ier) pat(n_pat)%x(i),pat(n_pat)%y(i),pat(n_pat)%sigma(i)
                if (ier /= 0) then
                   Err_diffpatt=.true.
                   Err_Mess_DiffPatt=" Error reading an ISIS profile DATA file"
                   return
                end if
                if(abs(pat(n_pat)%x(i)) < eps1 .and. pat(n_pat)%y(i) < eps1 .and. pat(n_pat)%sigma(i) < eps1) exit
                pat(n_pat)%x(i)=pat(n_pat)%x(i)
                pat(n_pat)%y(i)=pat(n_pat)%y(i)*fac_y
                pat(n_pat)%sigma(i)=pat(n_pat)%sigma(i)*fac_y
                pat(n_pat)%sigma(i)=pat(n_pat)%sigma(i)*pat(n_pat)%sigma(i)
                sumavar=sumavar+pat(n_pat)%sigma(i)
                if (pat(n_pat)%sigma(i) < eps1) pat(n_pat)%sigma(i) =fac_y
                if (pat(n_pat)%y(i) < eps1) then
                   pat(n_pat)%y(i)   = eps1
                   pat(n_pat)%sigma(i) = fac_y
                end if
                cnorm=cnorm+pat(n_pat)%sigma(i)/max(pat(n_pat)%y(i),0.001)
                if (i > 1) then
                   pat(n_pat)%step=pat(n_pat)%step+pat(n_pat)%x(i)-pat(n_pat)%x(i-1)
                   ntt=ntt+1
                end if
             end do
          end if  !RALF question

          pat(n_pat)%npts=ntt+1
          pat(n_pat)%xmin=pat(n_pat)%x(1)
          pat(n_pat)%xmax=pat(n_pat)%x(pat(n_pat)%npts)
          cnorm=cnorm/real(pat(n_pat)%npts)
          pat(n_pat)%step=pat(n_pat)%step/real(ntt)
          if (sumavar < eps1) then
             do i=1,pat(n_pat)%npts
                pat(n_pat)%sigma(i)=pat(n_pat)%y(i)
             end do
             cnorm=1.0
          end if

       end do !n_pat
       pat(n_pat)%ymin=minval(pat(n_pat)%y(1:pat(n_pat)%npts))
       pat(n_pat)%ymax=maxval(pat(n_pat)%y(1:pat(n_pat)%npts))
       return
    End Subroutine Read_Pattern_Isis_M

    !!--++
    !!--++ Subroutine Read_Pattern_Mult(Filename,Dif_Pat, NumPat, Mode)
    !!--++    character(len=*),                              intent (in) :: filename
    !!--++    type (diffraction_pattern_type), dimension(:), intent (out):: Dif_Pat
    !!--++    integer,                                       intent (out):: numpat
    !!--++    character(len=*), optional,                    intent (in) :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read one pattern from a Filename
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Mult(filename, dif_pat, numpat, mode)
       !---- Arguments ----!
       character(len=*),                                          intent (in)   :: filename
       type (diffraction_pattern_type), dimension(:),allocatable, intent (out)  :: dif_pat
       integer,                                                   intent (out)  :: numpat
       character(len=*), optional,                                intent (in)   :: mode

       !---- Local variables ----!
       logical :: esta
       integer :: i_dat, ier

       call init_err_diffpatt()

       inquire(file=filename,exist=esta)
       if ( .not. esta) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" The file "//trim(filename)//" doesn't exist"
          return
       else
          call get_logunit(i_dat)
          open(unit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error opening the file "//trim(filename)
             return
          end if
       end if

       if (present(mode)) then
          select case (u_case(mode))
              case ("XYSIGMA")
                 !   call  Read_Pattern_xysigma_m(dif_pat,npat)

              case ("ISIS")
                 call Read_Pattern_isis_m(i_dat,dif_pat,numpat)
                 dif_pat%diff_kind = "neutrontof"
                 dif_pat%scat_var =  "microsec"
                 dif_pat%instr  = mode

              case ("GSAS")
                 !   call Read_Pattern_gsas_m(dif_pat,npat)      ! GSAS file

              case default
                 Err_diffpatt=.true.
                 Err_Mess_DiffPatt="Invalid Mode"
                 return
          end select
          return
       end if
       close(unit=i_dat,iostat=ier)

       if (ier/=0) then
           Err_diffpatt=.true.
           Err_Mess_DiffPatt=" Problems closing data file"
       end if

       return
    End Subroutine Read_Pattern_Mult

    !!--++
    !!--++ Subroutine Read_Pattern_Nls(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for NLS
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Nls(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) ::  pat

       !---- Local Variables ----!
       character(len=132)                           :: aline
       integer                                      :: nlines,j,i,ier, no
       logical                                      :: title_given

       call init_err_diffpatt()
       title_given=.false.
       do
          read(unit=i_dat,fmt="(a)") aline
          aline=adjustl(aline)
          if(.not. title_given) then
            Pat%title=aline
            title_given=.true.
          end if
          if (aline(1:1) == "!") cycle
          read(unit=aline,fmt=*,iostat=ier) pat%xmin,pat%step,pat%xmax
          if (ier /= 0 ) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
             return
          end if
          exit
       end do

       pat%npts = (pat%xmax-pat%xmin)/pat%step+1.5
       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if

       nlines = pat%npts/10+1

       call Allocate_Diffraction_Pattern(pat)

       j = 0
       do i=1,nlines
          read(unit=i_dat,fmt="(10F8.0)",iostat=ier)(pat%y(j+no),no=1,10)
          if (ier /= 0 ) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in (NLS) Intensity file, check your instr parameter!1"
             return
          end if
          read(unit=i_dat,fmt="(10F8.0)",iostat=ier)(pat%sigma(j+no),no=1,10)
          if (ier /= 0 ) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in (NLS) Intensity file, check your instr parameter!2"
             return
          end if
          j = j+10
       end do

       pat%sigma(1) = pat%sigma(1)**2
       pat%x(1)=pat%xmin

       do i=2,pat%npts
          if (  pat%y(i) < 0.00001) pat%y(i) = pat%y(i-1)
          if (pat%sigma(i) < 0.00001) pat%sigma(i) = 1.0
          pat%sigma(i) = pat%sigma(i)**2
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))
       return
    End Subroutine Read_Pattern_Nls

    !!--++
    !!--++ Subroutine Read_Pattern_One(Filename,Dif_Pat, Mode)
    !!--++    character(len=*),                intent (in) :: filename
    !!--++    type (diffraction_pattern_type), intent(out) :: Dif_Pat
    !!--++    character(len=*), optional,      intent (in) :: mode
    !!--++
    !!--++    Read one pattern from a Filename
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_One(Filename,Dif_Pat, Mode)
       !---- Arguments ----!
       character(len=*),                intent (in)   :: filename
       type (diffraction_pattern_type), intent (out)  :: dif_pat
       character(len=*), optional,      intent (in)   :: mode

       !---- Local Variables ----!
       character(len=6)                               :: extdat !extension of panalytical file
       character(len=12)                              :: modem !extension of panalytical file
       logical                                        :: esta
       integer                                        :: i, i_dat,ier

       call init_err_diffpatt()

       inquire(file=filename,exist=esta)
       if (.not. esta) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" The file "//trim(filename)//" doesn't exist"
          return
       else
          call get_logunit(i_dat)
          open(unit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error opening the file: "//trim(filename)
             return
          end if
       end if

       if (present(mode)) then
          modem=u_case(mode)
       else
          modem="DEFAULT"
       end if

       select case (modem)
          case ("D1B" , "D20")
             call Read_Pattern_d1b_d20(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("NLS")                   ! Data from N.L.S (Brookhaven) Synchrotron Radiation  ,data from synchrotron source and correct data for dead time
             call Read_Pattern_nls(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("G41")                   ! Data from general format of two axis instruments with fixed step in twotheta
             call Read_Pattern_g41(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("D1A","D2B","3T2","G42")
             call Read_Pattern_d1a_d2b(i_dat,dif_pat)     ! Data from D1A,D2B  (Files *.sum, renamed *.dat, as prepared by D1ASUM or D2BSUM programs)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("D1AOLD", "D1BOLD","OLDD1A", "OLDD2B")
             call Read_Pattern_d1a_d2b_old(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("DMC","HRPT")                   ! Data from DMC,HRPT
             call Read_Pattern_dmc(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("SOCABIM")
             call  Read_Pattern_socabim(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("XYSIGMA")            !T.O.F. corrected data file from ISIS
             call  Read_Pattern_xysigma(i_dat, dif_pat)
             dif_pat%diff_kind = "neutrons_tof"
             dif_pat%scat_var =  "microsec"
             dif_pat%instr  = mode
             dif_pat%ct_step = .true.

          case ("GSAS")
             call Read_Pattern_gsas(i_dat,dif_pat)         ! GSAS file
             dif_pat%diff_kind = "neutrons_tof"
             dif_pat%scat_var =  "microsec"
             dif_pat%instr  = mode

          case ("PANALYTICAL")
             i=index(filename,".",back=.true.)
             extdat=u_case(filename(i:))
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode

             select case (extdat)
                case(".CSV")
                   CALL Read_Pattern_PANalytical_CSV(i_dat,dif_pat)

                case(".UDF")
                   CALL Read_Pattern_PANalytical_UDF(i_dat,dif_pat)

                case(".JCP")
                   CALL Read_Pattern_PANalytical_JCP(i_dat,dif_pat)

                case(".XRDML")
                   CALL Read_Pattern_PANalytical_XRDML(i_dat,dif_pat)
             end select

          case ("TIMEVARIABLE")
             call Read_Pattern_time_variable(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%instr  = mode

          case default
             call Read_Pattern_free(i_dat,dif_pat)
             dif_pat%instr  = "free format"

       end select

       close(unit=i_dat,iostat=ier)
       if (ier/=0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Problems closing the data file: "//trim(filename)
       end if

       return
    End Subroutine Read_Pattern_One

    !!--++
    !!--++ Subroutine Read_Pattern_Panalytical_CSV(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Panalitical Format CSV
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Panalytical_Csv(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: pat

       !---- Local Variables ----!
       character (len=132)                          :: line
       integer                                      :: i, j, long, ier
       real(kind=sp)                                :: alpha1, alpha2, ratio_I


       call init_err_diffpatt()

       !---- lecture fichier Philips X"celerator
       pat%ct_step = .false.
       do
          read(unit=i_dat,fmt="(a)",IOSTAT=ier) line
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile CSV-DATA file: end of file"
             return
          end if
          long=LEN_TRIM(line)

          if (line(1:7) == "Title 1") then
             pat%title=line(8:)
          else if(line(1:19) =="K-Alpha1 wavelength") then
             read(unit=line(21:long),fmt=*, IOSTAT=ier) alpha1
               pat%conv(1) = alpha1
          else if(line(1:19) =="K-Alpha2 wavelength") then
             read(unit=line(21:long),fmt=*, IOSTAT=ier) alpha2
               pat%conv(2) = alpha2
          else if(line(1:23) =="Ratio K-Alpha2/K-Alpha1") then
             read(unit=line(25:long),fmt=*, IOSTAT=ier) ratio_I
               pat%conv(3) = ratio_I
          else if(line(1:16) =="Data angle range") then
             read(unit=line(18:long),fmt=*)  pat%xmin  , pat%xmax

          else if(line(1:14) =="Scan step size") then
             read(unit=line(16:long),fmt=*, IOSTAT=ier) pat%step

          else if(line(1:13) =="No. of points") then
             read(unit=line(15:long),fmt=*, IOSTAT=ier) pat%npts

          else if(line(1:13) =="[Scan points]") then
             read(unit=i_dat,fmt="(a)",IOSTAT=ier) line     ! lecture de la ligne Angle,Intensity
             if (ier/=0) return
             exit
          end if
       end do

       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in (Csv)Intensity file, Number of poits negative or zero!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       i=0
       do
          i=i+1
          if (i > pat%npts) then
             i=i-1
             exit
          end if
          read(unit=i_dat,fmt="(a)",IOSTAT=ier) line
          j=index(line,",")
          if (j == 0) then
             read(unit=line,fmt=*,IOSTAT=ier) pat%x(i),pat%y(i)
          else
             read(unit=line(1:j-1),fmt=*,IOSTAT=ier) pat%x(i)
             read(unit=line(j+1:),fmt=*,IOSTAT=ier) pat%y(i)
          end if
          if (ier /=0) exit
          pat%sigma(i) = pat%y(i)
          if (pat%sigma(i)   <= 0.00001) pat%sigma(i) = 1.0
       end do

       pat%npts = i
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End Subroutine Read_Pattern_Panalytical_Csv

    !!--++
    !!--++ Subroutine Read_Pattern_Panalytical_JCP(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Panalitical Format JCP
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Panalytical_Jcp(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: pat

       !---- Local Variables ----!
       character (len=132)                          :: line
       integer                                      :: i, j, long , k, ier
       real(kind=sp)                                :: alpha1, alpha2, ratio_I

       call init_err_diffpatt()

       !---- lecture fichier JCP
       pat%ct_step = .false.
       k=0
       do
          k=k+1
          read(unit=i_dat,fmt="(a)",IOSTAT=ier) line
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile JCP-DATA file: end of file"
             return
          end if
          if( k == 1) pat%title=line
          long=LEN_TRIM(line)

          if (line(1:7) == "## END=") exit

          if (line(1:21) =="##$WAVELENGTH ALPHA1=") then
             read(unit=line(22:long),fmt=*, IOSTAT=ier) alpha1
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile  file: end of file"
                return
             end if
             pat%conv(1)= alpha1

          else if(line(1:21) =="##$WAVELENGTH ALPHA2=") then
             read(unit=line(22:long),fmt=*, IOSTAT=ier) alpha2
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if
              pat%conv(2)= alpha2

          else if(line(1:33) =="##$INTENSITY RATIO ALPHA2/ALPHA1=") then
             read(unit=line(34:long),fmt=*, IOSTAT=ier) ratio_I
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if
             pat%conv(3)= ratio_I

          else if(line(1:10) =="## FIRSTX=") then
             read(unit=line(11:long),fmt=*, IOSTAT=ier) pat%xmin
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if

          else if(line(1:10) =="## LASTX=") then
             read(unit=line(11:long),fmt=*, IOSTAT=ier) pat%xmax
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if

          else if(line(1:10) =="## DELTAX=") then
             read(unit=line(11:long),fmt=*, IOSTAT=ier) pat%step
             if (ier /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if

          else if(line(1:11) =="## NPOINTS=") then
             read(unit=line(12:long),fmt=*, IOSTAT=ier) pat%npts
             if (ier /= 0 .or. pat%npts <=0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                return
             end if

          else if(line(1:20) =="## XYDATA= X++<Y..Y>") then

             call Allocate_Diffraction_Pattern(pat)

             i=1
             do
                read(unit=i_dat,fmt="(f9.3,tr1,5f11.3)",iostat=ier) pat%x(i),(pat%y(i+j),j=0,4)
                if (ier /= 0) then
                   Err_diffpatt=.true.
                   Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
                   return
                end if

                do j=1,4
                   pat%x(i+j) = pat%x(i) + real(j)*pat%step
                end do
                if (i+5 > pat%npts ) exit
                i=i+5
             end do
             pat%npts = i

          end if
       end do ! File

       do i=1,pat%npts
          pat%sigma(i) = pat%y(i)
          if (pat%sigma(i)   <= 0.00001) pat%sigma(i) = 1.0
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End Subroutine Read_Pattern_Panalytical_Jcp

    !!--++
    !!--++ Subroutine Read_Pattern_Panalytical_UDF(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Panalitical Format UDF
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Panalytical_Udf(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: pat

       !---- Local Variables ----!
       character (len=132)                            :: line, newline
       integer                                        :: i, j, long, ier, n, nb_lignes, np
       real(kind=sp)                                  :: alpha1, alpha2, ratio !, ratio_I
       logical                                        :: title_given

       call init_err_diffpatt()

       !---- lecture fichier UDF
       pat%ct_step = .true.
       title_given = .false.

       do
          read(unit=i_dat,fmt="(a)",IOSTAT=ier) line
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile UDF-DATA file: end of file"
             return
          end if
          if(.not. title_given) then
            pat%title=line
            title_given=.true.
          end if
          long=LEN_TRIM(line)

          if (line(1:12) =="LabdaAlpha1,") then
             read(unit=line(23:long-2),fmt=*, IOSTAT=ier) alpha1
             pat%conv(1)=  alpha1
          else if(line(1:12) =="LabdaAlpha2,") then
             read(unit=line(23:long-2),fmt=*, IOSTAT=ier) alpha2
             pat%conv(2)=  alpha2

          else if(line(1:13) =="RatioAlpha21,") then
             read(unit=line(14:long-2),fmt=*, IOSTAT=ier) ratio
             pat%conv(3)= ratio

          else if(line(1:15) =="DataAngleRange,") then
             write(unit=newline,fmt="(a)")  line(16:long-2)
             i = INDEX(NewLine,",")
             long=LEN_TRIM(NewLine)
             read(unit=NewLine(1:i-1),fmt=*, IOSTAT=ier)   pat%xmin
             read(unit=NewLine(i+1:long),fmt=*,IOSTAT=ier) pat%xmax

          else if(line(1:13) =="ScanStepSize,") then
             read(unit=line(14:long-2),fmt=*, IOSTAT=ier) pat%step
             pat%npts=(pat%xmax-pat%xmin)/pat%step+1.2
             if (pat%npts <= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading a profile UDF-DATA file: end of file"
                return
             end if

          else if(line(1:7) =="RawScan") then

             call Allocate_Diffraction_Pattern(pat)

             nb_lignes = int(pat%npts/8)
             n = 0
             do j=1, nb_lignes
                read(unit=i_dat,fmt= "(7(f8.0,tr1),F8.0)", IOSTAT=ier) (pat%y(i+n),i=1,7), pat%y(n+8)
                if (ier /= 0) then
                   Err_diffpatt=.true.
                   Err_Mess_DiffPatt=" Error reading a profile UDF-DATA file: end of file"
                   return
                end if
                n = n + 8
             end do
             np = pat%npts - n

             if (np /= 0) then
                read(unit=i_dat, fmt = "(7(f8.0,tr1),F8.0)") (pat%y(i), i=n+1, pat%npts)
             endif
             exit

          end if
       end do ! file

       do i=1,pat%npts
          pat%x(i)=pat%xmin+real(i-1)*pat%step
          pat%sigma(i) = pat%y(i)
          if (pat%sigma(i)   <= 0.00001) pat%sigma(i) = 1.0
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End Subroutine Read_Pattern_Panalytical_Udf

    !!--++
    !!--++ Subroutine Read_Pattern_Panalytical_XRDML(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Panalitical Format XRDML
    !!--++
    !!--++ Update: January - 2005
    !!
    Subroutine Read_Pattern_Panalytical_Xrdml(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: pat

       !---- Local Variables ----!
       character (len=256000), allocatable, dimension(:)  :: XRDML_line, XRDML_intensities_line
       integer                                            :: i, i1, i2, nl, ier

       call init_err_diffpatt()

       if (allocated(XRDML_line))             deallocate(XRDML_line)
       if (allocated(XRDML_intensities_line)) deallocate(XRDML_intensities_line)

       allocate(XRDML_line(1))
       allocate(XRDML_intensities_line(1))

       !---- recherche de "<positions axis="2Theta" unit="deg">"
       do
          read(unit=i_dat, fmt="(a)", iostat=ier) XRDML_line(1)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
             return
          end if
          i1= index(XRDML_line(1), "<positions axis=""2Theta"" unit=""deg"">")
          if (i1/=0) exit
       end do

       !---- recherche de 2theta_min
       do
          read(unit=i_dat, fmt="(a)", iostat=ier) XRDML_line(1)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
             return
          end if
          i1= index(XRDML_line(1), "<startPosition>")
          if (i1==0) cycle
          i2= index(XRDML_line(1), "</startPosition>")
          read(unit=XRDML_line(1)(i1+15:i2-1), fmt=*) pat%xmin
          exit
       end do

       !---- recherche de 2theta_max
       do
          read(unit=i_dat, fmt="(a)", iostat=ier) XRDML_line(1)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
             return
          end if
          i1= index(XRDML_line(1), "<endPosition>")
          if (i1==0) cycle
          i2= index(XRDML_line(1), "</endPosition>")
          read(unit=XRDML_line(1)(i1+13:i2-1), fmt=*) pat%xmax
          exit
       end do

       !---- recherche des intensites
       do
          read(unit=i_dat, fmt="(a)", iostat=ier) XRDML_line(1)
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
             return
          end if
          i1= index(XRDML_line(1), "<intensities unit=""counts"">")
          if (i1==0) cycle
          i2= index(XRDML_line(1), "</intensities>")
          XRDML_intensities_line(1) = XRDML_line(1)(i1+27:i2-1)
          exit
       end do
       XRDML_intensities_line(1)=adjustl(XRDML_intensities_line(1))
       nl=LEN_TRIM(XRDML_intensities_line(1))
       i1=1
       do i=2,nl
          if (XRDML_intensities_line(1)(i:i) /= " ") then
             if (XRDML_intensities_line(1)(i-1:i-1) == " ") then
                i1=i1+1
             end if
          end if
       end do
       pat%npts=i1
       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       pat%step = (pat%xmax - pat%xmin) / real(pat%npts-1)
       read(unit=XRDML_intensities_line(1), fmt=*, iostat=ier) (pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error reading a profile XRDML-DATA file: end of file"
          return
       end if
       do i=1,pat%npts
          pat%x(i)=pat%xmin+real(i-1)*pat%step
          pat%sigma(i) = pat%y(i)
          IF (pat%sigma(i)   <= 0.00001) pat%sigma(i) = 1.0
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       if (allocated(XRDML_line))             deallocate(XRDML_line)
       if (allocated(XRDML_intensities_line)) deallocate(XRDML_intensities_line)

       return
    End Subroutine Read_Pattern_Panalytical_Xrdml

    !!--++
    !!--++ Subroutine Read_Pattern_Socabim(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Socabim
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Socabim(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) ::  pat

       !---- Local Variables ----!
       logical                                      :: string_counts, string_2thetacounts, string_2thetacps ,free_format
       character (len=132)                          :: line
       character(len=20),dimension(30)              :: dire
       character(len=1)                             :: separateur
       integer                                      :: i, j, i1, long, nb_sep, nb_col, n, ier
       real                                         :: step_time


       call init_err_diffpatt()

       string_COUNTS       = .false.
       string_2THETACOUNTS = .false.
       string_2THETACPS    = .false.
       free_format         = .false.

       nb_sep = 0    ! nombre de separateurs
       nb_col = 0    ! nombre de colonnes

       !---- recherche du type de donnees et de divers parametres (step, 2theta_min ...) ----!

        DO
           read(unit=i_dat,fmt="(a)",IOSTAT=ier) line
           if (ier/=0) then
              Err_diffpatt=.true.
              Err_Mess_DiffPatt=" Error on Socabim UXD Intensity file, check your mode parameter!"
              return
           end if
           IF (line(1:7) == "_COUNTS") THEN
              string_COUNTS    = .true.
              EXIT
           ELSE IF (line(1:13) =="_2THETACOUNTS") then
              string_2THETACOUNTS = .true.
              exit
           ELSE IF (line(1:10) == "_2THETACPS") THEN
              string_2THETACPS = .true.
              EXIT
           ELSE IF (line(1:8) == "_2THETA=") THEN
              i = INDEX(line,"=")
              j = LEN_TRIM(line)
              if (LEN_TRIM(line(i+1:j)) /=0) then
                 read(unit=line(i+1:j),fmt=*)  pat%xmin
              end if
           ELSE IF (line(1:10) == "_STEPSIZE=") THEN
              i=INDEX(line,"=")
              j = LEN_TRIM(line)
              if (LEN_TRIM(line(i+1:j)) /=0) then
                 read(unit=line(i+1:j),fmt=*)  pat%step
              end if
           ELSE IF (line(1:9) == "_STEPTIME") then
              i=INDEX(line,"=")
              j = LEN_TRIM(line)
              if (LEN_TRIM(line(i+1:j)) /=0) then
                 read(unit=line(i+1:j),fmt=*)  step_time
              end if
           ELSE IF (line(1:11) == "_STEPCOUNT=") THEN
              i=INDEX(line,"=")
              j = LEN_TRIM(line)
              if (LEN_TRIM(line(i+1:j)) /=0) then
                 read(unit=line(i+1:j),fmt=*)  pat%npts
              end if
           END IF
        END DO

        if (pat%npts <= 0) then
           Err_diffpatt=.true.
           Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
           return
        end if

        call Allocate_Diffraction_Pattern(pat)

        !---- lecture de la premiere ligne de donnees pour determiner le
        !---- format: format libre, type de separateur
        read(unit=i_dat,fmt= "(a)", IOSTAT=ier) line
        if (ier/=0) then
           Err_diffpatt=.true.
           Err_Mess_DiffPatt=" Error on Socabim UXD Intensity file, check your instr parameter!"
           return
        end if
        i1 = INDEX(line, CHAR(9))      ! "TAB" ?
        if (i1/=0) then
           separateur=CHAR(9)
        else
           i1 = INDEX(line, ";")         ! ";" ?
           if (i1/=0) then
              separateur=";"
           else
              i1 = INDEX(line,",")         ! ","
              if (i1/=0) separateur = ","
           end if
        end if

        if (i1==0) then   ! format libre  (separateur= caractere blanc)
           call getword(line,dire,nb_col)
           if (nb_col ==0) then
              Err_diffpatt=.true.
              Err_Mess_DiffPatt=" Error on Socabim UXD Intensity file, check your instr parameter!"
              return
           end if
           free_format = .true.
        else
           !---- determination du nombre de tabulations
           do
              nb_sep = nb_sep + 1
              long=LEN_TRIM(line)
              line = line(i1+1:long)
              i1=INDEX(line,separateur)
              if (i1==0) exit
           end do
           nb_col = nb_sep + 1
        end if

        if (string_2THETACOUNTS  .or. string_2THETACPS) nb_col = nb_col -1
        BACKSPACE(unit=i_dat)   ! on remonte d"une ligne

        !---- lecture des donnees
        j = 0       ! indice de la ligne
        n = 0       ! indice du comptage

        do
           j = j+1
           read(unit=i_dat,fmt= "(a)", IOSTAT=ier) line
           if (ier /= 0) exit
           IF (free_format) then
              call getword(line,dire,nb_col)
              if (nb_col==0) then
                 Err_diffpatt=.true.
                 Err_Mess_DiffPatt=" Error on Socabim UXD Intensity file, check your instr parameter!"
                 return
              end if
              if (string_2THETACOUNTS  .or. string_2THETACPS)then
                 nb_col = nb_col - 1           !  << corrected 14.03.02
                 read(unit=line,fmt=*,IOSTAT=ier) pat%x(n+1), (pat%Y(n+i),i=1,nb_col)
                 if (ier/=0) then
                    n=n-1
                    exit
                 end if
              else
                 read(unit=line,fmt=*, IOSTAT=ier) (pat%Y(n + i),i=1,nb_col)
                 if (ier/=0) then
                    n=n-1
                    exit
                 end if
              end if
              n = n + nb_col

           else
              if (string_2THETACOUNTS  .or. string_2THETACPS)then
                 i1=INDEX(line,separateur)
                 long=LEN_TRIM(line)
                 read(unit=line(1:i1-1),fmt=*, IOSTAT=ier)pat%x(1+nb_col*(j-1))
                 if (ier/=0)  exit
                 line = line(i1+1:long)
              end if

              !---- lecture des comptages d'une ligne
              if (nb_sep > 1) then
                 do i =1, nb_sep
                    n=n+1
                    i1=INDEX(line, separateur)
                    long=LEN_TRIM(line)
                    if (i1==0) then
                       n=n-1
                       exit
                    end if
                    read(unit=line(1:i1-1), fmt=*, IOSTAT=ier) pat%Y(n)

                    if (ier/=0) then
                       n=n-1
                       exit
                    end if
                    j=j+1
                    line= line(i1+1:long)
                 end do
              end if

              !---- lecture du dernier point de la ligne
              n =n + 1
              read(unit=line, fmt=*, IOSTAT=ier) pat%Y(n)
              if (ier/=0) exit
           end if
        end do

        pat%npts = n
        pat%xmax = pat%xmin + pat%step*(pat%npts-1)  !! TR 28.11.02

        !---- creation des abcisses
        !---- modif. des comptages si necessaire et calculs sigmas_comptages

        if (string_COUNTS .or. string_2THETACOUNTS ) then
           pat%sigma(1:n ) = pat%Y(1:n )
        else  ! data in CPS
           pat%sigma(1:n ) = pat%Y(1:n )/ step_time
        endif

        where (pat%sigma(:) <= 0.00001) pat%sigma(:) = 1.0

        do i=1,pat%npts
           pat%x(i)= pat%xmin+(i-1)*pat%step
        end do
        pat%ymin=minval(pat%y(1:pat%npts))
        pat%ymax=maxval(pat%y(1:pat%npts))

        return
     End subroutine Read_Pattern_Socabim

    !!--++
    !!--++ Subroutine Read_Pattern_Time_Variable(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for Time Variable
    !!--++
    !!--++ Update: January - 2005
    !!
    Subroutine Read_Pattern_Time_Variable(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) ::  pat

       !---- Local Variables ----!
       character(len=132)                           :: txt1
       character(len=132)                           :: txt2
       character(len=132)                           :: txt3
       real, dimension(:), allocatable              :: bk
       real(kind=sp)                                :: cnorma
       integer                                      :: i,ier

       call init_err_diffpatt()

       read(unit=i_dat,fmt="(A)",iostat=ier)txt1   !1
       pat%title=txt1
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier)txt2   !2
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier)txt3   !3
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt="(A)",iostat=ier)txt3                  !4
       if (ier /= 0)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       read(unit=i_dat,fmt=*,iostat=ier)pat%xmin,pat%step,pat%xmax
       if (ier /= 0)then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       pat%npts = (pat%xmax-pat%xmin)/pat%step+1.5
       if (pat%npts <=0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(bk) ) deallocate(bk)
       allocate(bk(pat%npts))

       read(unit=i_dat,fmt="(5(F6.0,F10.0))",iostat=ier)(bk(i),pat%y(i),i=1,pat%npts)
       if (ier /= 0) then
          Err_diffpatt=.true.
          Err_Mess_DiffPatt=" Error in  Intensity file, check your instr parameter!"
          return
       end if

       !---- Normalize data to constant time
       cnorma=0.0
       DO i=1,pat%npts
          IF (bk(i) < 1.0E-06) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Zero time in *.DAT "
             return
          end if
          cnorma=cnorma+bk(i)
          pat%x(i)= pat%xmin+(i-1)*pat%step
       end do
       cnorma=cnorma/real(pat%npts)
       do i=1,pat%npts
          pat%y(i)=pat%y(i)*cnorma/bk(i)
          pat%sigma(i)=pat%y(i)
          bk(i)=0.0
       end do
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End subroutine Read_Pattern_Time_Variable

    !!--++
    !!--++ Subroutine Read_Pattern_XYSigma(i_dat,Pat)
    !!--++    integer,                         intent(in)  :: i_dat
    !!--++    type (diffraction_pattern_type), intent(out) :: Pat
    !!--++
    !!--++    Read a pattern for X,Y,Sigma
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_XYSigma(i_dat,Pat)
       !---- Arguments ----!
       integer,                         intent(in)  :: i_dat
       type (diffraction_pattern_type), intent(out) :: Pat

       !---- Local Variables ----!
       character(len=132)                           :: txt1, aline, fmtfields, fmtformat
       character (len=5)                            :: date1
       integer                                      :: line_da, ntt, interpol, i, j,ier,npp
       real                                         :: fac_x, fac_y,  yp1, sumavar, cnorm
       real(kind=sp)                                :: ycor, xt, stepin, ypn
       real(kind=sp), parameter                     :: eps1=1.0E-1
       real, dimension(:), allocatable              :: yc, bk

       call init_err_diffpatt()

       !---- Or X,Y sigma data ----!
       fac_x=1.0
       fac_y=1.0
       interpol=0
       line_da=1
       npp=0
       ntt=0

       do
          read(unit=i_dat,fmt="(a)", iostat=ier) txt1
          if (ier /= 0)   exit
          npp=npp+1
       end do

       pat%npts  =   npp -6
       if (pat%npts <= 0) then
          Err_diffpatt=.true.
          Err_mess_diffpatt=" Error in Intensity file, Number of poits negative or zero!"
          return
       end if
       rewind(unit=i_dat)

       read(unit=i_dat,fmt="(a)") txt1
       pat%title=txt1
       do i=1,5
          read(unit=i_dat,fmt="(a)", iostat=ier) txt1
          if (ier /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile DATA file"
             return
          end if

          line_da=line_da+1
          if (txt1(1:5) == "INTER") then !Interpolation possible!
             backspace (unit=i_dat)
             line_da=line_da-2
             call init_findfmt(line_da)
             fmtfields = "5ffif"
             call findfmt(i_dat,aline,fmtfields,fmtformat)
             if (ierr_fmt /= 0) then
                Err_diffpatt=.true.
                Err_Mess_DiffPatt=" Error reading"
                return
             end if

             read(unit=aline,fmt=fmtformat) date1,fac_x,fac_y,interpol,stepin
             if (fac_x <= 0.0) fac_x=1.0
             if (fac_y <= 0.0) fac_y=1.0
          end if

          if (txt1(1:4) == "TEMP") then
             read(unit=txt1(5:80),fmt=*) aline
          end if
       end do

       if (interpol == 0) then
          pat%ct_step = .false.
       else if(interpol == 1) then
          pat%ct_step = .true.
       else if(interpol == 2) then
          pat%ct_step = .true.
       end if

       call Allocate_Diffraction_Pattern(pat)

       if(allocated(bk) ) deallocate(bk)
       allocate(bk(pat%npts))
       if(allocated(yc) ) deallocate(yc)
       allocate(yc(pat%npts))


       fmtfields = "fff"
       sumavar=0.0
       cnorm=0.0
       i=0
       do j=1,pat%npts
          call findfmt(i_dat,aline,fmtfields,fmtformat)
          if (ierr_fmt == -1) exit
          if (ierr_fmt /= 0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error reading a profile DATA file"
             return
          end if
          i=i+1
          read(unit=aline,fmt = fmtformat, iostat=ier )pat%x(i),pat%y(i),pat%sigma(i)
          if (ier /=0) then
             Err_diffpatt=.true.
             Err_Mess_DiffPatt=" Error in Intensity file, check your instr parameter!"
             return
          end if
          IF (ABS(pat%x(i)) < eps1 .AND. pat%y(i) < eps1 .AND.  pat%sigma(i) < eps1) exit
          pat%x(i)=pat%x(i)*fac_x
          pat%y(i)=pat%y(i)*fac_y
          pat%sigma(i)=pat%sigma(i)*fac_y
          pat%sigma(i)=pat%sigma(i)*pat%sigma(i)
          sumavar=sumavar+pat%sigma(i)
          if(pat%sigma(i) < eps1) pat%sigma(i) =1.0
          if(pat%y(i) < eps1) then
             pat%y(i)   = eps1
             pat%sigma(i) =1.0
          end if
          cnorm=cnorm+pat%sigma(i)/MAX(pat%y(i),0.001)
          if(i > 1) then
            pat%step=pat%step+pat%x(i)-pat%x(i-1)
          end if
       end do

       pat%xmin=pat%x(1)
       pat%xmax=pat%x(pat%npts)
       cnorm=cnorm/REAL(pat%npts)
       if (sumavar < eps1) then
          do i=1,pat%npts
             pat%sigma(i)=pat%y(i)
          end do
          cnorm=1.0
       end if

       if (interpol == 0) then      !if interpol
          pat%step=pat%step/real(ntt)
       else                        !else interpol
          ntt=ntt+1
          pat%step=stepin
          pat%npts=(pat%x(ntt)-pat%x(1))/pat%step+1.05

          !---- First interpolate the raw intensities ----!
          yp1=9.9E+32
          ypn=9.9E+32
          call spline(pat%x(:),pat%y(:),ntt,yp1,ypn,bk(:))
          do i=1,pat%npts
             xt=pat%x(1)+(i-1)*pat%step
             call splint(pat%x(:),pat%y(:),bk(:),ntt,xt,ycor)
             yc(i)=MAX(1.0,ycor)
          end do
          do i=1,pat%npts
             pat%y(i)=yc(i)
             yc(i)=0.0
             bk(i)=0.0
          end do

          !---- Second interpolate the sigmas ----!
          call spline(pat%x(:),pat%sigma(:),ntt,yp1,ypn,bk(:))
          do i=1,pat%npts
             xt=pat%x(1)+(i-1)*pat%step
             call splint(pat%x(:),pat%sigma(:),bk(:),ntt,xt,ycor)
             yc(i)=max(1.0,ycor)
          end do
          do i=1,pat%npts
             pat%sigma(i)=yc(i)
             yc(i)=0.0
             bk(i)=0.0
          end do
          pat%xmax=pat%xmin+pat%step*(pat%npts-1)
       end if                       !End If interpol
       pat%ymin=minval(pat%y(1:pat%npts))
       pat%ymax=maxval(pat%y(1:pat%npts))

       return
    End Subroutine Read_Pattern_XYSigma

    !!--++
    !!--++ Subroutine Set_Background_Inter(Difpat,Bcky,Bckx,N)
    !!--++    type (diffraction_pattern_type), intent(in out)  :: difPat
    !!--++    real (kind=sp), dimension(:),    intent(in out ) :: bcky
    !!--++    real (kind=sp), dimension(:),    intent(in out ) :: bckx
    !!--++    integer,                         intent(in    )  :: n
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a Background
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Set_Background_Inter(Difpat,Bcky,Bckx,N)
       !---- Arguments ----!
       type (diffraction_pattern_type),intent(in out) :: difpat
       real (kind=sp), dimension(:),   intent(in out) :: bcky
       real (kind=sp), dimension(:),   intent(in out) :: bckx
       integer,                        intent(in    ) :: n

       !---- Local variables ----!
       integer                                        :: nbx, nbac1 , i , j  , nxx
       real                                           :: difl, difr , thx , delt, slope, bstep,p


       nbx=1
       nbac1=n

       difl=bckx(1)-difpat%xmin
       difr=bckx(n)-difpat%xmax

       if (difl >= 0) then
          if (difpat%ct_step) then
             nbx=difl/difpat%step+1.5
          else
             nbx=locate(difpat%x(:),difpat%npts,bckx(1))
             if (nbx <= 0) nbx=1
          end if
          do i=1,nbx
             difpat%bgr(i)=bcky(1)
          end do
       end if

       if (difr <= 0) then
          nbac1=n+1
          bckx(nbac1)=difpat%xmax
          bcky(nbac1)=bcky(n)
       end if

       nxx=2
       do_i: do i=nbx,difpat%npts
          thx=difpat%x(i)
          do j=nxx,nbac1
             delt=bckx(j)-thx
             if (delt > 0.0) then
                p=bckx(j)-bckx(j-1)
                if (abs(p) > 0.0001) then
                   slope=(bcky(j)-bcky(j-1))/p
                else
                   slope=0.0
                end if
                bstep=(thx-bckx(j-1))*slope
                difpat%bgr(i)=bcky(j-1)+bstep
                nxx=j-1
                cycle do_i
             end if
          end do
       end do  do_i

       return
    End Subroutine Set_Background_Inter

    !!--++
    !!--++ Subroutine Set_Background_Poly( Difpat,Bkpos,Bckx,N)
    !!--++    type (diffraction_pattern_type), intent(in out) :: difPat
    !!--++    real (kind=sp),                  intent(in    ) :: bkpos
    !!--++    real (kind=sp), dimension(:),    intent(in    ) :: bckx
    !!--++    integer,                         intent(in    ) :: n
    !!--++
    !!--++    (PRIVATE)
    !!--++    Define a Background
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Set_Background_Poly( Difpat,Bkpos,Bckx,N)
       !---- Arguments ----!
       type (diffraction_pattern_type), intent(in out) :: difpat
       real (kind=sp),                  intent(in    ) :: bkpos
       real (kind=sp), dimension(:),    intent(in    ) :: bckx
       integer,                         intent(in    ) :: n

       !---- Local Variables ----!
       integer                                         :: i,j

       if (allocated(difpat%bgr) ) deallocate(difpat%bgr)
       allocate(difpat%bgr(difpat%npts))

       do i=1, difpat%npts
          difpat%bgr(i)=0
          do j=1,n
             difpat%bgr(i)= difpat%bgr(i)+ bckx(j)*((difpat%x(i)/bkpos-1.0)**(j-1))
          end do
       end do

       return
    End Subroutine Set_Background_Poly

 End Module Diffraction_Patterns_Mod


