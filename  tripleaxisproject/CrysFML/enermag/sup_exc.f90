Module Super_Exchange

  Use String_Utilities, only: Frac_Trans_1Dig, Pack_String

  implicit none
  private
  public ::  Get_Expo, rkky, equiv_jotas, init_exchange_interaction, &
             write_exchange_interaction, Get_vect

  integer, parameter, public :: num_se=50, num_sse=100

  interface  Get_Expo
    module procedure Get_Expo_c
    module procedure Get_Expo_r
  end interface

  Type, public :: Exchange_interaction
     character(len=3)     :: J       ! Name of J
     character(len=4)     :: nam1    ! Nature of atom1
     character(len=4)     :: nam2    ! Nature of atom2
     integer              :: ns      ! Number of superexchange paths
     integer              :: nss     ! Number of super-superexchange paths
     real                 :: valj    ! value of the exchange integral
     real                 :: dist    ! distance betwen atoms atom1-atom2
     character(len=40)               :: de_nam  ! atom1-atom2(transl)
     character(len=40),dimension(6)  :: s_nam   ! atom1-atom2-atom3(transl)
     real,             dimension(3,6):: se_geo  ! super-exchange geometry
                                                ! Up to 6 s-e paths:(d1,d2,theta)1,(d1,d2,theta)2, ...(d1,d2,theta)6
     character(len=50),dimension(6)  :: ss_nam  ! atom1-atom2-atom3-atom4(transl)
     real,             dimension(6,6):: sse_geo ! super-exchange geometry (up to six ss-e) (d1,d2,d3,ang1,ang2,dihed,dtot)
  End Type  Exchange_interaction

  Type, public :: Direct_Ex_path
     character(len=10)    :: nam1  ! atom1
     character(len=10)    :: nam2  ! atom2
     character(len=40)    :: nam   ! atom1-atom2(transl)
     real                 :: dist  ! Distance atom1-atom2
     real, dimension(3,2) :: coord ! Fractional coordinates of atom k:  coord(:,k) (k=1,2)
     real, dimension(3,2) :: carte ! Cartesian  coordinates of atom k:  carte(:,k) (k=1,2)
  End Type  Direct_Ex_path

  Type, public :: S_Ex_path
     character(len=10)    :: nam1  ! atom1
     character(len=10)    :: nam2  ! atom2
     character(len=10)    :: nam3  ! atom3
     character(len=40)    :: nam   ! atom1-atom2-atom3(transl)
     real, dimension(4)   :: geom  !(d12, d23, angle(123), dtot )
     real, dimension(3,3) :: coord !Fractional coordinates of atom k:  coord(:,k) (k=1,2,3)
     real, dimension(3,3) :: carte !Cartesian  coordinates of atom k:  carte(:,k) (k=1,2,3)
  End Type  S_Ex_path

  Type, public :: SS_Ex_path
     character(len=10)    :: nam1  ! atom1
     character(len=10)    :: nam2  ! atom2
     character(len=10)    :: nam3  ! atom3
     character(len=10)    :: nam4  ! atom4
     character(len=50)    :: nam   ! atom1-atom2-atom3-atom4(transl)
     real, dimension(7)   :: geom  !(d12, d23, d24,angle(123),angle(234), Dihedral(1234), dtot)
     real, dimension(3,4) :: coord !Fractional coordinates of atom k:  coord(:,k) (k=1,2,3,4)
     real, dimension(3,4) :: carte !Cartesian  coordinates of atom k:  carte(:,k) (k=1,2,3,4)
  End Type  SS_Ex_path

  Type, public :: SE_Connection
    integer :: nd              !Integer to indicate if direct-Exchange between two atoms is allowed
    integer :: ns              !Number of Super-Exchange paths between two atoms
    integer :: nss             !Number of Super-Super-Exchange paths between two atoms
    Type (Direct_Ex_path)                 :: DE
    Type (S_Ex_path), dimension(num_se)   :: SE
    Type (SS_Ex_path), dimension(num_sse) :: SSE
  End Type SE_Connection


  contains

   Subroutine Get_vect(trans,vect)
    character (len=*),  intent(in)  :: trans
    real, dimension(3), intent(out) :: vect
    character (len=5), dimension(3)  :: XX
    integer:: i1,i2,i
    real   :: num, denom
      i1=index(trans,",")
      i2=index(trans,",",back=.true.)
      XX(1)=trans(2:i1-1)
      XX(2)=trans(i1+1:i2-1)
      XX(3)=trans(i2+1:len_trim(trans)-1)
      do i=1,3
        i1=index(XX(i),"/")
        if(i1 /= 0) then
          read(unit=xx(i)(1:i1-1),fmt=*) num
          read(unit=xx(i)(i1+1:),fmt=*) denom
          vect(i)= real(num)/real(denom)
        else
          read(unit=xx(i),fmt=*) vect(i)
        end if
      end do
      return
   End Subroutine Get_vect

   Subroutine Get_Expo_c(trans, expo)
    character (len=*),  intent(in)  :: trans
    character (len=30), intent(out) :: expo
                                 ! 1234567890123456
                                 ! (-X/2+ Y/4-8Z/9) -> -X/2+Y/4-8Z/9
                                 ! (-1/2, 1/4,-8/9) -> -X/2+Y/4-8Z/9
    character (len=16) :: ex
    character (len=5)  :: XX,YY,ZZ
    integer:: i1,i2
      i1=index(trans,",")
      i2=index(trans,",",back=.true.)
      XX=trans(2:i1-1)
      YY=trans(i1+1:i2-1)
      ZZ=trans(i2+1:len_trim(trans)-1)
   !!  XX
      i1=index(XX,"/")
      if (i1 /= 0) then
        XX=XX(1:i1-1)//"X"//XX(i1:)
      else
        XX=trim(XX)//"X"
      end if
      i1=index(XX,"X")
      if(XX(i1-1:i1-1) == "1") XX(i1-1:i1-1)=" "
      if(XX(1:1) == "0")  XX="     "

   !!  YY
      i1=index(YY,"/")
      if (i1 /= 0) then
          YY=YY(1:i1-1)//"Y"//YY(i1:)
      else
          YY=trim(YY)//"Y"
      end if
      if(YY(1:1) /= "-") YY="+"//YY
      i1=index(YY,"Y")
      if(YY(i1-1:i1-1) == "1") YY(i1-1:i1-1)=" "
      if(YY(i1-1:i1-1) == "0") YY="     "

   !!  ZZ
      i1=index(ZZ,"/")
      if (i1 /= 0) then
          ZZ=ZZ(1:i1-1)//"Z"//ZZ(i1:)
      else
          ZZ=trim(ZZ)//"Z"
      end if
      if(ZZ(1:1) /= "-") ZZ="+"//ZZ
      i1=index(ZZ,"Z")
      if(ZZ(i1-1:i1-1) == "1") ZZ(i1-1:i1-1)=" "
      if(ZZ(i1-1:i1-1) == "0")  ZZ="     "

      ex=XX//YY//ZZ
     if(len_trim(ex) == 0) then
      expo(1:30)=" "
     else
      ex= Pack_String(ex)
      if(ex(1:1) == "+") ex(1:1) = " "
      expo="exp{-2pi("//trim(ex)//")}"
     end if
   End Subroutine Get_Expo_c

   Subroutine Get_Expo_r(v, expo)
    real, dimension(3), intent(in)  :: v
    character (len=30), intent(out) :: expo
    character (len=16) :: transla, ex
                                   ! 1234567890123456
                                   ! (-X/2+ Y/4-8Z/9) -> -X/2+Y/4-8Z/9
     call Frac_Trans_1Dig(v,transla)        ! (-1/2, 1/4,-8/9) -> -X/2+Y/4-8Z/9
     ex=transla
        ex(16:16)=" "    !remove parentheses
        ex(1:1)=" "
   !!  XXXXXX
     if(ex(3:3) /= "1") then
        ex(1:1)=ex(2:2)
        ex(2:2)=ex(3:3)
     end if
        ex(3:3)="X"
        if(ex(2:2) == "0") ex(1:5)=" "
   !!  YYYYYYY
     if(ex(8:8) /= "1") then
        ex(6:6)=ex(7:7)
        ex(7:7)=ex(8:8)
        if(ex(6:6) == " ") ex(6:6)="+"
     else
        ex(6:6)=" "
        if(ex(7:7) == " ") ex(7:7)="+"
     end if
        ex(8:8)="Y"
        if(ex(7:7) == "0") ex(6:10)=" "
   !!  ZZZZZZZ
     if(ex(13:13) /= "1") then
        ex(11:11)=ex(12:12)
        ex(12:12)=ex(13:13)
        if(ex(11:11) == " ") ex(11:11)="+"
     else
        ex(11:11)=" "
        if(ex(12:12) == " ") ex(12:12)="+"
     end if
        ex(13:13)="Z"
        if(ex(12:12) == "0") ex(11:15)=" "
     if(len_trim(ex) == 0) then
      expo(1:30)=" "
     else
      transla= Pack_String(ex)
      if(transla(1:1) == "+") transla(1:1) = " "
      expo="exp{-2pi("//trim(transla)//")}"
     end if
   End Subroutine Get_Expo_r

   Pure Function rkky(dist,kf) result(jval)
     real, intent( in) :: dist
     real, intent( in) :: kf
     real :: jval

     real :: x
     x=kf*dist
     jval= 100.0 * (sin(x)-x*cos(x))/(x*x*x*x)
     return

   end Function rkky

   Pure Function Equiv_jotas(j1,j2) result(equiv)
     Type (Exchange_interaction), intent( in) :: j1
     Type (Exchange_interaction), intent( in) :: j2
     Logical :: equiv
     !Local variables
     real, parameter :: eps=0.002, epa=0.02
     integer, dimension(6) :: neq
     real :: d1,d2,d3,ang1,ang2,ang3
     integer :: i,j
     equiv=.false.
     !Compare two exchange interactions

     !Compare first the nature of the extreme atoms
     if(.not. ( (j1%nam1 == j2%nam1 .and. j1%nam2 == j2%nam2) .or. &
                (j1%nam1 == j2%nam2 .and. j1%nam2 == j2%nam1))        ) return

     !Compare the number of super-exchange and super-super-exchange paths
     if(.not. (j1%ns == j2%ns .and. j1%nss == j2%nss)) return

     !Compare the super-exchange paths
     if( j1%ns /= 0) then
       neq=0
       do i=1,j1%ns
          d1=j1%se_geo(1,i)
          d2=j1%se_geo(2,i)
         ang1=j1%se_geo(3,i)
         do j=1,j2%ns
          if(neq(j) > 0) cycle
          if(  (abs(d1-j2%se_geo(1,j)) < eps .and. abs(d2-j2%se_geo(2,j)) < eps .and. abs(ang1-j2%se_geo(3,j)) < epa) .or. &
               (abs(d1-j2%se_geo(2,j)) < eps .and. abs(d2-j2%se_geo(1,j)) < eps .and. abs(ang1-j2%se_geo(3,j)) < epa) ) then
           neq(j)=neq(j) + 1
           exit
          end if
         end do
       end do
       do j=1,j1%ns
         if(neq(j) == 0) return   !A difference in a path has been found
       end do
     end if

     if( j1%nss /= 0) then
       !Compare the super-super-exchange paths
          neq=0
       do i=1,j1%nss
          d1=j1%sse_geo(1,i)
          d2=j1%sse_geo(2,i)
          d3=j1%sse_geo(3,i)
          ang1=j1%sse_geo(4,i)
          ang2=j1%sse_geo(5,i)
          ang3=j1%sse_geo(6,i)
         do j=1,j2%nss
          if(neq(j) > 0) cycle
          if(abs(d2-j2%sse_geo(2,j)) < eps ) then      !Comparing A-A' distance
            if(  (abs(d1-j2%sse_geo(1,j)) < eps .and. abs(d3-j2%sse_geo(3,j)) < eps) .or. &    !comparing M-A and A'-M'
                 (abs(d1-j2%sse_geo(3,j)) < eps .and. abs(d3-j2%sse_geo(1,j)) < eps)   ) then
                 if(abs(abs(ang3)-abs(j2%sse_geo(6,j))) < epa) then
                    if ((abs(ang1-j2%sse_geo(4,j)) < epa .and. abs(ang2-j2%sse_geo(5,j)) < epa) .or. &
                        (abs(ang2-j2%sse_geo(4,j)) < epa .and. abs(ang1-j2%sse_geo(5,j)) < epa)         ) then
                      neq(j)=neq(j) + 1
                      exit
                    end if
                 end if
            end if
          end if
         end do
       end do
       do j=1,j1%nss
         if(neq(j) == 0) return   !A difference in a path has been found
       end do
     end if
     equiv=.true.
     return

   end Function Equiv_jotas

   subroutine init_exchange_interaction(J)
    type (Exchange_interaction), intent (in out) :: J
    !Initializa the object J
    J%J=" "
    J%nam1=" "
    J%nam2=" "
    J%ns=0
    J%nss=0
    J%valj=0.0
    J%dist=0.0
    J%de_nam=" "
    J%s_nam=" "
    J%ss_nam=" "
    J%se_geo=0.0
    J%sse_geo=0.0
    return



   end subroutine init_exchange_interaction

   subroutine write_exchange_interaction(lun,J)
      integer,                     intent (in) :: lun
      type (Exchange_interaction), intent (in) :: J
      !Write the object J
      integer :: i
      write(unit=lun,fmt="(/,a,a,f8.4)") " => Exchange interaction: "//j%J//" between "//j%nam1//" and ", &
                                            j%nam2//" -> Distance: ",j%dist
      do i=1,j%ns
        write(unit=lun,fmt="(a,a,2f8.4,f8.2,a)")   "    S-Exchange (d1,d2,theta): ", &
                                                   trim(j%s_nam(i))//"  (",j%se_geo(:,i)," )"
      end do
      do i=1,j%nss
        write(unit=lun,fmt="(a,a,3f8.4,3f8.2,a)")  "    SS-Exchange (d1,d2,d3,ang1,ang2,dihed): ",&
                                                   trim(j%ss_nam(i))//"  (",j%sse_geo(:,i)," )"
      end do
      return
   end subroutine write_exchange_interaction

End Module Super_Exchange

