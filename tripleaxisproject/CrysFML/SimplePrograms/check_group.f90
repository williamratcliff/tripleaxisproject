  Program Check_Group
    use crystallographic_symmetry, only: Space_Group_Type, set_spacegroup
    use reflections_utilities, only: Hkl_Absent
    use string_utilities, only: u_case
    use Symmetry_Tables, only: spgr_info, Set_Spgr_Info
    implicit none

    type(Space_Group_Type)             :: Spacegroup
    character(len=80)                  :: title
    character(len=256)                 :: filehkl
    character(len=12)                  :: hms,system_n
    character(len=15)                  :: hall
    character(len=1)                   :: csystem,cent
    logical                            :: check_cent=.true., absent
    integer,         dimension(200)    :: num_group
    integer(kind=2), dimension(200)    :: merit
    integer(kind=4), dimension(3,8000) :: hkl
    integer(kind=2), dimension(  8000) :: good = 0
    real,    dimension(  5000)         :: Intensity, sigma, bragg, fwhm
    integer, parameter                 :: Num_Spgr_Info = 612, Monoc=  15, Ortor = 163, Hexag  = 527, &
                                          Cubic = 554, Tetra = 410, Trigo = 495
    real                               ::  angmax,maxint,dfw,threshold,dd
    integer                            ::  nhkl,nmax,ngood
    integer                            ::  i,j,m,ier,i1,i2,numg

    character(len=10),dimension(6):: argument = " "
    integer :: narg,iargc

    !---- Arguments on the command line ----!
    !
    !  It is supposed that the program is invoked as follows
    !  > check_group  filehkl Csystem Cent AngMax DFW Threshold
    !
    !  The three first arguments are compulsory (characters separated by spaces)
    !  Thre three numerical values may be omitted. Default values are available.

    narg=iargc()

    if(narg /= 0) then
      if(narg > 6) narg=6
      do i=1,narg
        call getarg(i,argument(i))
      end do
      filehkl=argument(1)          !Code of hkl file
      csystem = argument(2)        !Crystal system  (M,O,T,H,R,C)
      cent=argument(3)             !Check centred lattices? (y/n)
      if(narg > 3) then
        read(unit=argument(4),fmt=*) angmax     !Maximum 2theta
        if(narg > 4) then
          read(unit=argument(5),fmt=*) dfw      !Minimum FWHM distance between adjacent reflections
          if(narg > 5) then
             read(unit=argument(6),fmt=*) threshold  !% of the maximum intensity to consider observed reflections
          else
             threshold = 1.0
          end if
        else
          dfw=3.0
          threshold = 1.0
        end if
      else
        angmax=180.0
        dfw=3.0
        threshold = 1.0
      end if

    else
       write(unit=*,fmt=*) " => Please enter the code of the hkl-file: "
       read(unit=*,fmt=*) filehkl
       write(unit=*,fmt=*) " => Please enter the crystallographic system (M,O,T,R,H,C): "
       read(unit=*,fmt="(a)") csystem
       write(unit=*,fmt=*) " => Check centred lattices? (y=<cr>/n): "
       read(unit=*,fmt="(a)") cent
       write(unit=*,fmt=*) " => Please enter the maximum value of the Scattering angle (2theta): "
       read(unit=*,fmt=*) angmax
       write(unit=*,fmt=*) " => Minimum FWHM distance between adjacent reflections (e.g.:1,2...): "
       read(unit=*,fmt=*) dfw
       write(unit=*,fmt=*) " => % of the maximum intensity to consider observed reflections (e.g.: 1): "
       read(unit=*,fmt=*) threshold
    end if

    if(len_trim(csystem) == 0) stop
    csystem= u_case(csystem)
    if(index("MOTRHC",csystem) == 0) stop
    if(len_trim(cent) == 0) cent="Y"
    cent=u_case(cent)
    if(cent /= "Y" .and. cent /= "N") cent="Y"
    if(cent == "N") check_cent=.false.

    open(unit=1,file=trim(filehkl)//".hkl",status="old",action="read", position="rewind")
    open(unit=2,file=trim(filehkl)//".spg",status="replace",action="write")
    read(unit=1,fmt="(a)") title
    read(unit=1,fmt=*) nhkl
    maxint=0.0
    do i=1,nhkl
      read(unit=1,fmt=*,iostat=ier) hkl(:,i),m,intensity(i),sigma(i),bragg(i),fwhm(i)
      if(maxint < intensity(i))  maxint = intensity(i)
      if(ier /= 0 .or. (bragg(i) > angmax+dfw*fwhm(i)) ) then
        nhkl=i-1
        exit
      end if
    end do
    if(nhkl < 3) then
      write(unit=*,fmt=*) " => Too few reflections! "
      stop
    end if
    Select Case(csystem)
      case("M")
        i1=monoc
        i2=ortor-1
        system_n="Monoclinic"
      case("O")
        i1=ortor
        i2=tetra-1
        system_n="Orthorhombic"
      case("T")
        i1=tetra
        i2=trigo-1
        system_n="Tetragonal"
      case("R")
        i1=trigo
        i2=hexag-1
        system_n="Rhombohedral"
      case("H")
        i1=hexag
        i2=cubic-1
        system_n="Hexagonal"
      case("C")
        i1=cubic
        i2=Num_Spgr_Info
        system_n="Cubic"
    End Select


    !Pointer to useful reflections
     good(1)=1
     ngood=1
     write(unit=2,fmt="(a)") "   ----------------------------------------------------------------------- "
     write(unit=2,fmt="(a)") "     PROGRAM CHECK_GROUP: attempt to select the possible space groups from "
     write(unit=2,fmt="(a)") "                          an experimental Powder Diffraction Pattern       "
     write(unit=2,fmt="(a)") "   ----------------------------------------------------------------------- "
     write(unit=2,fmt="(a)") "     Author: J.Rodriguez-Carvajal (version 0.01, based on CrysFML)         "
     write(unit=2,fmt="(a)") "   ----------------------------------------------------------------------- "
     write(unit=2,fmt="(a)") "   "
     write(unit=2,fmt="(a)")     "   Conditions: "
     write(unit=2,fmt="(a,a)")   "                Input hkl-file      : ", trim(filehkl)//".hkl"
     write(unit=2,fmt="(a,a)")   "                Crystal System      : ", system_n
     write(unit=2,fmt="(a,a)")   "                Check centred cells?: ", cent
     write(unit=2,fmt="(a,f8.4)")"                Maximum angle       : ", angmax
     write(unit=2,fmt="(a,f8.4)")"                Number of FWHMs     : ", dfw
     write(unit=2,fmt="(a,f8.4)")"                Threshold in %      : ", threshold
     write(unit=2,fmt="(a)") "   "
     write(unit=2,fmt="(a,/)") " => List of read reflections:"
     write(unit=2,fmt="(a)")   "       h   k   l   Intensity       Sigma      2theta       FWHM  Good?"
     write(unit=2,fmt="(a,3i4,4f12.4,i4)") "    ", hkl(:,1), Intensity(1), Sigma(1), Bragg(1), fwhm(1),good(1)
     do i=2,nhkl-1
        dd=dfw*(fwhm(i-1)+fwhm(i)+fwhm(i+1))/3.0
        if(bragg(i-1) < bragg(i)-dd  .and.  bragg(i) < bragg(i+1)-dd) then
           good(i)=1
           ngood=ngood+1
        end if
        write(unit=2,fmt="(a,3i4,4f12.4,i4)") "    ", hkl(:,i), Intensity(i), Sigma(i), Bragg(i), fwhm(i),good(i)
     end do

    if(ngood < 5) then
      write(unit=*,fmt=*) " => Too few GOOD reflections!  "
      write(unit=*,fmt=*) " => ... diminish de value of DFW, now DFW= ",dfw
      stop
    end if

    call Set_Spgr_Info()

    threshold=threshold*maxint*0.01
        m=0
        do_group: do i=i1,i2
          hms=adjustl(spgr_info(i)%HM)
          hall=spgr_info(i)%hall
          if( hms(1:1) /= "P" .and. .not. check_cent ) cycle do_group ! Skip centred groups
          call set_spacegroup(hall,Spacegroup,Force_Hall="y")
          do j=1,nhkl
             if(good(j) == 0) cycle  !Skip reflections that are not good (overlap) for checking
             absent=Hkl_Absent(hkl(:,j), Spacegroup)
             if(absent .and. intensity(j) > threshold) cycle do_group !Group not allowed
          end do
          ! Passing here means that all reflections are allowed in the group -> Possible group!
          m=m+1
          num_group(m)=i
        end do  do_group

        write(unit=2,fmt="(/a,i3)")      " => Number of good reflections  : ",ngood
        write(unit=2,fmt="( a,f12.4)")   "    Maximum intensity           : ",maxint
        write(unit=2,fmt="( a,f12.4)")   "    Minimum (for observed)      : ",threshold
        write(unit=2,fmt="( a,i3   )")   "    Number of Space Group tested: ",i2-i1+1
        write(unit=2,fmt="(/a,i3,a)") " => LIST OF POSSIBLE SPACE GROUPS, a total of ",m," groups are possible"
        write(unit=2,fmt="(a)") "     ------------------------------------------------------"
        write(unit=2,fmt="(a)") "     Number(IT)      Hermann-Mauguin Symbol     Hall Symbol"
        write(unit=2,fmt="(a)") "     ------------------------------------------------------"

        write(unit=*,fmt=*) " => LIST OF POSSIBLE SPACE GROUPS, a total of ",m," groups are possible"
        write(unit=*,fmt=*) "     ------------------------------------------------------"
        write(unit=*,fmt=*) "     Number(IT)      Hermann-Mauguin Symbol     Hall Symbol"
        write(unit=*,fmt=*) "     ------------------------------------------------------"
        do i=1,m
          j=num_group(i)
          hms=adjustl(spgr_info(j)%HM)
          hall=spgr_info(j)%hall
          numg=spgr_info(j)%N
          write(unit=*,fmt="(i10,4a)")  numg,"                   ",hms,"          ",hall
          write(unit=2,fmt="(i10,4a)")  numg,"                   ",hms,"          ",hall
        end do
     stop
  End Program Check_Group
