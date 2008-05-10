!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: STRING_UTILITIES
!!----   INFO: Manipulation of strings with alfanumeric characters
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----            October - 1999: Reorder the subroutines and functions
!!----                            All routines have general I/O parameters
!!----
!!---- DEPENDENCIES
!!--++    Use Math_Gen, only: Sp, Negligible Zbelong
!!----
!!---- VARIABLES
!!--++    CTAB                    [Private]
!!--++    DIGIT                   [Private]
!!----    ERR_MESS_STRING
!!----    ERR_STRING
!!----    ERR_TEXT_TYPE
!!--++    IENDFMT                 [Private]
!!----    IERR_FMT
!!--++    IERRCHARBEGG            [Private]
!!--++    IERREFRMT               [Private]
!!--++    IERREOF                 [Private]
!!--++    IERREMPTYFIELD          [Private]
!!--++    IERRFIELDTYPE           [Private]
!!--++    IERRFIELDS              [Private]
!!--++    IERRINVALC              [Private]
!!--++    IERRINVALCHAR           [Private]
!!--++    IERRINVALFIELD          [Private]
!!--++    IERRIO                  [Private]
!!--++    IERRNONE                [Private]
!!--++    IERRNUMBER              [Private]
!!--++    IERRSEPMISS             [Private]
!!--++    IERRSTRLENGTH           [Private]
!!--++    IINTE                   [Private]
!!--++    IREAL                   [Private]
!!--++    I_NINE                  [Private]
!!--++    I_ONE                   [Private]
!!--++    I_ZERO                  [Private]
!!--++    LINE_NB                 [Private]
!!----    MESS_FINDFMT
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       EQUAL_SETS_TEXT
!!----       L_CASE
!!----       PACK_STRING
!!----       U_CASE
!!----
!!----    Subroutines:
!!--++       BUILDFMT             [Private]
!!----       CUTST
!!----       FINDFMT
!!--++       FINDFMT_ERR          [Private]
!!----       FRAC_TRANS_1DIG
!!----       FRAC_TRANS_2DIG
!!----       GET_FRACTION_1DIG
!!----       GET_FRACTION_2DIG
!!----       GET_LOGUNIT
!!----       GETNUM
!!----       GETNUM_STD
!!----       GETWORD
!!----       INIT_ERR_STRING
!!----       INIT_FINDFMT
!!----       LCASE
!!----       NUMBER_LINES
!!----       READ_KEY_STR
!!----       READ_KEY_STRVAL
!!----       READ_KEY_VALUE
!!----       READ_KEY_VALUEST
!!----       READING_LINES
!!----       SETNUM_STD
!!--++       SGETFTMFIELD         [Private]
!!--++       TREATMCHARFIELD      [Private]
!!--++       TREATNUMERFIELD      [Private]
!!----       UCASE
!!----
!!
 Module String_Utilities
    !---- Use Modules ----!
    use Math_Gen, only: Sp, Negligible, Zbelong

    implicit none

    private

    !---- List of public functions ----!
    public :: Equal_Sets_Text,  L_Case, Pack_String, U_Case

    !---- List of public subroutines ----!
    public :: Cutst, Get_Fraction_1Dig, Get_Fraction_2Dig, Getnum, Getnum_std, Getword, &
              Init_err_String, lcase, Number_lines, Read_Key_str, Read_Key_strVal,      &
              Read_Key_Value, Read_Key_ValueST, Reading_Lines, Setnum_std, Ucase,       &
              FindFmt,  Init_FindFmt, Frac_Trans_1Dig, Frac_Trans_2Dig, get_logunit

    !---- List of private subroutines ----!
    private :: BuildFmt, TreatNumerField, TreatMCharField, SgetFtmField, FindFmt_Err


    !---- Definitions ----!

    !!--++
    !!--++ CTAB
    !!--++    character (len=*), private, parameter :: cTab=Char(9)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Character parameter for TAB
    !!--++
    !!--++ Update: February - 2005
    !!
    character (len=*), private, parameter :: cTab=Char(9)

    !!--++
    !!--++ DIGIT
    !!--++    character (len=*), private, parameter :: digit="0123456789.-"
    !!--++
    !!--++    (PRIVATE)
    !!--++    Character parameter for numbers
    !!--++
    !!--++ Update: February - 2005
    !!
    character (len=*), private, parameter :: digit="0123456789.-"

    !!----
    !!---- ERR_MESS_STRING
    !!----    character(len=150) :: err_mess_string
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_string

    !!----
    !!---- ERR_STRING
    !!----    logical :: err_string
    !!----
    !!----    Logical Variable indicating an error in STRING_UTILITIES module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_string

    !!----
    !!---- TYPE :: ERR_TEXT_TYPE
    !!--..
    !!---- Type :: Err_Text_Type
    !!----    integer :: nlines
    !!----    character (len=132), dimension(5) :: txt
    !!---- End Type Err_Text_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, Public :: Err_Text_Type
       integer :: nlines
       character (len=132), dimension(5) :: txt
    End Type Err_Text_Type

    !!--++
    !!--++ IENDFMT
    !!--++    integer, paramater, private :: iEndFMT=0
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for EndFMT
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iEndFMT=0

    !!----
    !!---- IERR_FMT
    !!----    integer :: ierr_fmt
    !!----
    !!----    Integer signaling if an error has occurred (/=0) in using the procedure findFMT
    !!----
    !!---- Update: February - 2005
    !!
    integer, public :: iErr_fmt  ! Error code value (should be normally = 0)

    !!--++
    !!--++ IERRCHARBEGG
    !!--++    integer, paramater, private :: iErrCharBegg=4
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrCharBegg=4

    !!--++
    !!--++ IERREFRMT
    !!--++    integer, paramater, private :: iErrEfrmt=11
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter  :: iErrEfrmt=11

    !!--++
    !!--++ IERREOF
    !!--++    integer, paramater, private :: iErrEof=-1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrEof=-1

    !!--++
    !!--++ IERREMPTYFIELD
    !!--++    integer, paramater, private :: iErrEmptyField=8
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrEmptyField=8

    !!--++
    !!--++ IERRFIELDTYPE
    !!--++    integer, paramater, private :: iErrFieldType=3
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter  :: iErrFieldType=3

    !!--++
    !!--++ IERRFIELDS
    !!--++    integer, paramater, private :: iErrFields=1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter  :: iErrFields=1

    !!--++
    !!--++ IERRINVALC
    !!--++    integer, paramater, private :: iErrInvalC=5
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrInvalC=5

    !!--++
    !!--++ IERRINVALCHAR
    !!--++    integer, paramater, private :: iErrInvalChar=7
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrInvalChar=7

    !!--++
    !!--++ IERRINVALFIELD
    !!--++    integer, paramater, private :: iErrInvalField=6
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter  :: iErrInvalField=6

    !!--++
    !!--++ IERRIO
    !!--++    integer, paramater, private :: iErrIO=2
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrIO=2

    !!--++
    !!--++ IERRNONE
    !!--++    integer, paramater, private :: iErrNone=0
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrNone=0

    !!--++
    !!--++ IERRNUMBER
    !!--++    integer, paramater, private :: iErrNumber=12
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrNumber=12

    !!--++
    !!--++ IERRSEPMISS
    !!--++    integer, paramater, private :: iErrSepMiss=10
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrSepMiss=10

    !!--++
    !!--++ iErrStrLength
    !!--++    integer, paramater, private :: iErrStrLength=9
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for Error code
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iErrStrLength=9

    !!--++
    !!--++ IINTE
    !!--++    integer, paramater, private :: iInte=-1
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for iInte
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iInte=-1

    !!--++
    !!--++ IREAL
    !!--++    integer, paramater, private :: iReal=-2
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for iReal
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: iReal=-2

    !!--++
    !!--++ I_NINE
    !!--++    integer, paramater, private :: i_Nine=57
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for ASCII code for Nine
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: i_Nine=57

    !!--++
    !!--++ I_ONE
    !!--++    integer, paramater, private :: i_One=49
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for ASCII code for One
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: i_One=49

    !!--++
    !!--++ I_ZERO
    !!--++    integer, paramater, private :: i_Zero=48
    !!--++
    !!--++    (PRIVATE)
    !!--++    Integer parameter for ASCII code for Zero
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer , private, parameter :: i_Zero=48

    !!--++
    !!--++ LINE_NB
    !!--++    integer :: Line_Nb
    !!--++
    !!--++    (PRIVATE)
    !!--++    Line number updated each time the procedure findFMT is called
    !!--++    To initialize LINE_NB, the subroutine Init_FindFMT should be called.
    !!--++
    !!--++ Update: February - 2005
    !!
    Integer, private :: Line_Nb   ! Line number

    !!----
    !!---- MESS_FINDFMT
    !!----    Type (Err_Text_Type) :: Mess_FindFMT
    !!----
    !!----    Text composed of a maximum of 5 lines to inform about position or error
    !!----    in free format reading (used by procedure findFMT)
    !!----
    !!---- Update: February - 2005
    !!
    Type (Err_Text_Type), public :: Mess_FindFMT = Err_Text_Type(0,(/" "," "," "," "," "/))

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Logical Function Equal_Sets_Text(Text1,N1,Text2,N2) Result(Equal_sets_texto)
    !!----    character(len=*), dimension(:), intent(in) :: Text1   ! In -> String array
    !!----    integer,                        intent(in) :: N1      ! In -> Lines on Text1 variable
    !!----    character(len=*), dimension(:), intent(in) :: Text2   ! In -> String array
    !!----    integer,                        intent(in) :: N2      ! In -> Lines on Text2 variable
    !!----    logical                                    :: Equal_sets_texto
    !!----
    !!----    Determine if two sets of text lines are equal irrespective of the
    !!----    order of the lines. The function is true if the two sets of text
    !!----    have the same lines in whatever order.  Two lines are equal only
    !!----    if they have the same length and all their component characters
    !!----    are equal and in the same order.
    !!----
    !!---- Update: February - 2005
    !!
    Function Equal_Sets_Text(text1,n1,text2,n2) result(Equal_sets_texto)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: text1,text2
       integer,                        intent(in) :: n1,n2
       logical                                    :: Equal_sets_texto

       !---- Local variables ----!
       integer :: i,j
       logical :: info

       Equal_sets_texto=.false.

       if (n1 /= n2) return
       if (len(text1) /= len(text2)) return

       do i=1,n1
          info=.false.
          do j=1,n2
             if (text1(i) == text2(j)) then
                info=.true.
                exit
             end if
          end do
          if (.not. info) return
       end do

       Equal_sets_texto=.true.

       return
    End Function Equal_Sets_Text

    !!----
    !!---- Character Function L_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String: "InPUT Line"
    !!----    character (len=len(text))     :: mtex   ! Out -> String: "input line"
    !!----
    !!----    Conversion to lower case, text is not modified
    !!----
    !!---- Update: February - 2005
    !!
    Function L_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "A" .and. mtext(pos:pos) <= "Z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) - inc )
       end do

       return
    End Function L_Case


    !!----
    !!---- Character Function Pack_String(String) Result (Strp)
    !!----    character (len=*), intent(in) :: String
    !!----    character (len=*)             :: Strp
    !!----
    !!----    Pack a string: the function provides a string without empty spaces
    !!----
    !!---- Update: February - 2005
    !!
    Function Pack_String(String) Result (Strp)
       !---- Argument ----!
       character (len=*), intent(in)    :: string
       character (len=len_trim(string)) :: strp

       !---- Local variables ----!
       integer ::  i,n

       n=0
       strp=" "
       do i=1,len(string)
          if (string(i:i) /= " ") then
             n=n+1
             strp(n:n)=string(i:i)
          end if
       end do

       return
    End Function Pack_String

    !!----
    !!---- Character Function U_Case(Text) Result (Mtext)
    !!----    character (len=*), intent(in) :: text   !  In -> String:"Input Line"
    !!----    character (len=len(text))     :: mtext  ! Out -> String:"INPUT LINE"
    !!----
    !!----    Conversion to upper case, text is not modified
    !!----
    !!---- Update: February - 2005
    !!
    Function U_Case(Text) Result (Mtext)
       !---- Argument ----!
       character (len=*), intent(in) :: text
       character (len=len(text))     :: mtext

       !---- Local variables ----!
       integer, parameter :: inc = ICHAR("A") - ICHAR("a")
       integer            :: leng, pos

       mtext=text
       leng=len_trim(mtext)
       do pos=1,leng
          if (mtext(pos:pos) >= "a" .and. mtext(pos:pos) <= "z")           &
              mtext(pos:pos) = CHAR ( ICHAR(mtext(pos:pos)) + inc )
       end do

       return
    End Function U_Case

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!--++
    !!--++ Subroutine BuildFMT(iFld,nCar,nStr,FMTstring)
    !!--++    Integer,           intent(in    ) ::   iFld       -> Format type
    !!--++    Integer,           intent(in out) ::   nCar       -> integer/real field: number of characters in field
    !!--++                                                      -> character field: number of characters to skip before A field
    !!--++    Integer,           intent(in out) ::   nStr      <-> current character number in FMTstring
    !!--++    Character (len=*) ,intent(in out) ::   FMTstring <-> FORTRAN format string
    !!--++
    !!--++    (PRIVATE)
    !!--++    Add a new field to the FMT string
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine BuildFMT(iFld,nCar,nStr,FMTstring)
       !---- Arguments ----!
       Integer,           intent(in    ) ::   iFld
       Integer,           intent(in out) ::   nCar
       Integer,           intent(in out) ::   nStr
       Character (len=*) ,intent(in out) ::   FMTstring

       !---- Local variables ----!
       Integer ::  N

       !---- heading symbol "F"
       nStr = nStr + 1
       if (nStr > Len(FMTstring)) then
          iErr_fmt = iErrStrLength          ! format string length exceeded
          return
       end if

       if (iFld == iInte) then
          FMTstring(nStr:nStr)  = "i"   !descriptor are in lower case to be F-compatible
       else if (iFld == iReal) then
          FMTstring(nStr:nStr)  = "f"
       else if (iFld > 0) then
          if (nCar == 0) then
             FMTstring(nStr:nStr)  = "a"
          else
             if (nCar < 10) then
                write(unit=FMTstring(nStr:),fmt="(a,i1,a)") "tr",nCar,",a"
             else
                write(unit=FMTstring(nStr:),fmt="(a,i2,a)") "tr",nCar,",a"
             end if
             nStr=len_trim(FMTstring)
          end if
       end if

       !---- numeric part of Integer and real fields
       if (iFld < 0) then
          !---- hundredth ----!
          if (nCar >= 100) then
             N = Int(nCar/100)
             nStr = nStr + 1
             if (nStr > Len(FMTstring)) then
                iErr_fmt = iErrStrLength          ! format string length exceeded
                return
             end if
             FMTstring(nStr:nStr) = Char(N+48)
             nCar = nCar - N*100
          end if

          !---- tenth ----!
          if (nCar >= 10) then
             N = Int(nCar/10)
             nStr = nStr + 1
             if (nStr > Len(FMTstring)) then
                iErr_fmt = iErrStrLength          ! format string length exceeded
                return
             end if
             FMTstring(nStr:nStr) = Char(N+48)
             nCar = nCar - N*10
          end if

          !---- units ----!
          nStr = nStr + 1
          if (nStr > Len(FMTstring)) then
             iErr_fmt = iErrStrLength          ! format string length exceeded
             return
          end if
          FMTstring(nStr:nStr) = Char(nCar+48)

          !---- Add ".0" to the end of real fields ----!
          if (iFld == iReal) then
             nStr = nStr + 2
             if (nStr > Len(FMTstring)) then
                iErr_fmt = iErrStrLength          ! format string length exceeded
                return
             end if
             FMTstring(nStr-1:nStr) = ".0"
          end if

       else if (iFld > 0) then
          !---- numeric part of "A" fields ----!
          nStr = nStr + 1
          if (nStr > Len(FMTstring)) then
             iErr_fmt = iErrStrLength          ! format string length exceeded
             return
          end if
          FMTstring(nStr:nStr)   = Char(iFld)
       end if

       !---- Add a separator "," after each new FORTRAN field ----!
       nStr = nStr + 1
       if (nStr > Len(FMTstring)) then
          iErr_fmt = iErrStrLength          ! format string length exceeded
          return
       end if
       FMTstring(nStr:nStr) = ","

       return
    End Subroutine BuildFMT


    !!----
    !!---- Subroutine Cutst(Line1, Nlong1, Line2, Nlong2)
    !!----    character(len=*),           intent(in out) :: Line1   !  In -> Input string
    !!----                                                          ! Out -> Input string without the first word
    !!----    integer,          optional, intent(   out) :: Nlong1  ! Out -> Give the length of Line1 on Output
    !!----    character(len=*), optional, intent(   out) :: Line2   ! Out -> The first word of String on Input
    !!----    integer,          optional, intent(   out) :: Nlong2  ! Out -> Give the length of Line2 on Output
    !!----
    !!----    Removes the first word of the input String.
    !!----    Provides (optionally) a string with the first word.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Cutst(line1,nlong1,line2,nlong2)
       !---- Argument ----!
       character (len=*),           intent(in out) :: line1
       character (len=*), optional, intent(   out) :: line2
       integer,           optional, intent(   out) :: nlong1
       integer,           optional, intent(   out) :: nlong2

       !---- Local variables ----!
       integer  :: k,iniz1

       !---- Initializing variables ----!
       if (present(nlong1)) nlong1=0
       if (present(nlong2)) nlong2=0

       !---- Initializing to blank the directive ----!
       if (present(line2)) line2=" "

       !---- Elimination of possible blanks on the left ----!
       line1=adjustl(line1)
       if (len_trim(line1) <= 0) return

       k=len(line1)
       iniz1=index(line1," ")

       if (k ==1) then
          if (present(line2)) line2=line1
          if (present(nlong2)) nlong2=1
          line1=" "
       else
          if (iniz1 > 0) then
             if (present(line2))  line2=line1(1:iniz1-1)
             if (present(nlong2)) nlong2=len_trim(line1(1:iniz1-1))
             line1=line1(iniz1:)
          else
             if (present(line2))  line2=line1
             if (present(nlong2)) nlong2=len_trim(line1)
             line1=" "
          end if
       end if

       line1=adjustl(line1)
       if(present(nlong1)) nlong1=len_trim(line1)

       return
    End Subroutine Cutst

    !!----
    !!---- Subroutine FindFmt(Lun,aLine,FMTfields,FMTstring)
    !!----    Integer ,           intent(in    ) ::  Lun         !  -> Logical unit number
    !!----    Character (len=*) , intent(in out) ::  aLine       ! <-> character string to be decoded
    !!----    Character (len=*) , intent(in    ) ::  FMTfields   ! <-> description of the format fields (e.g. IIFIF)
    !!----    Character (len=*) , intent(   out) ::  FMTstring   ! <-  format of the line (e.g. (I5,I1,F8.0,I4,F7.0,) )
    !!--<<
    !!----    The routine "FindFmt" emulates the free format data input
    !!----    Read(unit=String1,fmt="(a,i,2f,..)") aString,i1,R1,R2,...
    !!----    but with additional error checking. Thus, given a description
    !!----    of the expected fields "FindFmt" returns the format of the line
    !!----    to be decoded. Valid field descriptors are:
    !!----    I:integer; R:real; A:free A format; 1 to 5:A1 to A5
    !!----
    !!----    This routine have an associated FindFMT error code (iErr_fmt)
    !!----      -2 : FORTRAN read error
    !!----      -1 : End of file
    !!----       0 : No Error
    !!----       1 : empty format descriptor (0 field)
    !!----       2 : data string read error
    !!----       3 : integer field found real !
    !!----       4 : begged dot, sign or "e" character !
    !!----       5 : invalid character in an integer field !
    !!----       6 : invalid field in format descriptor !
    !!----       7 : invalid character in a numeric field !
    !!----       8 : 0 character in current field !
    !!----       9 : format string length exceeded !
    !!----      10 : separator missing !
    !!----      11 : incomplete E or D format !
    !!----      12 : incomplete number !
    !!----
    !!----   An error message is generated ans written to the public variable "Mess_FindFMT"
    !!----   Consult the structure of Mess_FindFMT that is of type: Err_Text_Type.
    !!-->>
    !!--..   Example of use:
    !!--..       Character aLine*(*),FMTfields*(*),FMTstring*(*),String*5
    !!--..       Parameter (iLun=30)       ! input logical unit number
    !!--..
    !!--..    !-- Usual fixed format input (e.g.)
    !!--..    Read(unit=iLun,fmt="(4x,a5,i3,1x,2f8.2,i5)") String,i1,R1,R2,i2
    !!--..
    !!--..    !-- Free format input (Read performed by FindFMT)
    !!--..       FMTfields = "5iffi"
    !!--..       Call FindFmt(Lun,aLine,FMTfields,FMTstring)
    !!--..       if (iErr_fmt == -1) GoTo 998  ! End of Line| Block treating
    !!--..       if (iErr_fmt /= 0)  GoTo 999  ! input error|   errors
    !!--..       Read(unit=aLine,fmt=FMTstring) String,i1,R1,R2,i2
    !!--..
    !!--..    !-- Free format input (Read performed by calling routine)
    !!--..       Read(unit=iLun,fmt="(a)") aLine
    !!--..       FMTfields = "5iffi"
    !!--..       Call FindFmt(0,aLine,FMTfields,FMTstring,nC_L)
    !!--..       if (iErr_fmt == -1) GoTo 998 ! End of Line | Block treating
    !!--..       if (iErr_fmt /= 0)  GoTo 999 ! input error |   errors
    !!--..       Read(unit=aLine,fmt=FMTstring) String,i1,R1,R2,i2
    !!--..       ......
    !!--..   998 Continue ! End of file encountered
    !!--..       ......
    !!--..    !-- Output error message if any
    !!--..   999 Continue
    !!--..        if(ierr_fmt /= 0 .and. Mess_FindFMT%nlines > 0) then
    !!--..          do i=1,Mess_FindFMT%nlines
    !!--..           Write(unit=lun,fmt="(a)") Mess_FindFMT%txt(i)
    !!--..          end do
    !!--..        end if
    !!--..        ........
    !!--..
    !!---- Update: February - 2005
    !!
    Subroutine FindFmt(Lun,aLine,FMTfields,FMTstring)
       !---- Arguments ----!
       Character (len=*) , intent(in out) ::  aLine
       Character (len=*) , intent(in    ) ::  FMTfields
       Character (len=*) , intent(   out) ::  FMTstring
       Integer ,           intent(in    ) ::  Lun      ! Logical unit number

       !---- Local variables ----!
       Character (len=len(FMTfields)) ::  UFMTfields
       Integer  :: nC_L     ! counts characters in Line
       Integer  :: ioS      ! Fortran status code
       Integer  :: L_Fields ! true length of format descriptor
       Integer  :: L_Line   ! true length of data line
       Integer  :: nCar     ! counts characters in current format field
       Integer  :: nFld     ! counts format fields in FMTfields
       Integer  :: nStr     ! counts characters in FMTstring
       Integer  :: iFld     ! field type -1:integer;-2:real;>0:A1 to A9
       Integer  :: GetFTMfield     ! old function now argument of a subroutine
       Logical  :: ifSearchEnd

       !---- Initialize ----!
       nC_L = 0
       nFld = 0
       FMTstring = "()"     ! will receive FORTRAN format
       nStr = 1             ! at least a right parentheses in FMTstring
       iErr_fmt = iErrNone
       L_Fields  = Len_trim(FMTfields)
       line_nb = line_nb + 1  ! Update the line number

       !---- Format descriptor in upper case ----!
       if (FMTfields == " ") then
          iErr_fmt = iErrFields           ! empty FMT format descriptor
          Call FindFMT_Err(aLine,nC_L)
          Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
          Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
               " => Please check your input file at line: ",Line_Nb," !"
               return
       end if
       UFMTfields=FMTfields
       Call UCase(UFMTfields)

       !---- (Get and) verify data line ----!
       if (Lun > 0) then
          do
             Read(unit=Lun,fmt="(a)",ioStat=ioS) aLine
             if (ioS == -1) then
                iErr_fmt = iErrEof            ! End Of File
                Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
                Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4)") " => Non FATAL End of file !,  logical unit: ",Lun
                return                    !leave reading routine to handle end of file

             else if (ioS > 0) then
                iErr_fmt = -ioS-100           ! FORTRAN read error
                Call FindFMT_Err(aLine,nC_L)
                Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
                Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                     " => Please check your input file at line: ",Line_Nb," !"
                return
             end if

             l_line = len_trim(aLine)    ! true length without trailing spaces
             if (aLine(1:1) == "!" .or. aLine(1:1) == "#" .or. L_line == 0) then
                Line_Nb=Line_Nb+1
             else
                exit
             end if
          end do
       end if

       !---- Start decoding line character by character ----!
       ifSearchEnd = .false.

       do
          if (ifSearchEnd) exit

          !---- Get a new format field type ----!
          nCar = 0                    ! new format field
          call SGetFTMfield(GetFTMfield,UFMTfields, nFld, L_fields)
          iFld = GetFTMfield
          if (iErr_fmt /= iErrNone) then ! Error in field definition
             Call FindFMT_Err(aLine,nC_L)
             Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
             Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                  " => Please check your input file at line: ",Line_Nb," !"
             return
          end if
          if (iFld == iEndFMT) then   ! format exhausted
             if (nFld == 0) then
                iErr_fmt = iErrInvalField   ! invalid field in FMTfields
                Call FindFMT_Err(aLine,nC_L)
                Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
                Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                     " => Please check your input file at line: ",Line_Nb," !"
                return
             else
                exit                    ! scan end
             end if
          end if

          !---- Decode current field (character or numeric ?) ----!
          if (iFld > iEndFMT) then
             Call TreatMCharField(iFld,aLine,L_Line,nC_L,nCar)
          else if (iFld == iEndFMT) then    ! format exhausted
             exit
          else if (iFld < iEndFMT) then
             Call TreatNumerField(iFld,aLine,L_Line,nC_L,nCar)
          end if
          if (iErr_fmt /= iErrNone) then
             Call FindFMT_Err(aLine,nC_L)
             Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
             Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                  " => Please check your input file at line: ",Line_Nb," !"
             return
          end if
          if ((iFld < iEndFMT.and.nCar == 0) .or. iFld == 0) then
             iErr_fmt = iErrEmptyField           ! no characters in field
             return
          end if

          !---- Build current FMT element ----!
          Call BuildFMT(iFld,nCar,nStr,FMTstring)
          if (iErr_fmt /= iErrNone) then   ! format string length exceeded
             Call FindFMT_Err(aLine,nC_L)
             Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
             Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                  " => Please check your input file at line: ",Line_Nb," !"
             return
          end if

          !---- End of data Line ? ----!
          if (nC_L >= L_Line) ifSearchEnd = .true.
       end do

       !---- Terminates and close the format field ----!

       !---- If FMT not exhausted we append the remaining fields to ----!
       !---- the format string                                      ----!
       if (iErr_fmt == iErrNone .and. nFld < L_Fields) then
          !do while (iFld /= iEndFMT)
          do
             if (iFld == iEndFMT) exit
             call SGetFTMfield(GetFTMfield,UFMTfields, nFld, L_fields)
             iFld = GetFTMfield
             if (iErr_fmt /= iErrNone) then   ! Error in field definition
                Call FindFMT_Err(aLine,nC_L)
                Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
                Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                     " => Please check your input file at line: ",Line_Nb," !"
                return
             end if
             if (iFld /= iEndFMT) then
                nCar=1     !Put ==1 because BuildFMT required INOUT arg.
                Call BuildFMT(iFld,nCar,nStr,FMTstring)
                if (iErr_fmt /= iErrNone) then ! format string length exceeded
                   Call FindFMT_Err(aLine,nC_L)
                   Mess_FindFMT%nlines=Mess_FindFMT%nlines+1
                   Write(unit=Mess_FindFMT%txt(Mess_FindFMT%nlines),fmt="(a,i4,a)")    &
                        " => Please check your input file at line: ",Line_Nb," !"
                   return
                end if
             end if
          end do
       end if

       !---- Close format string ----!
       FMTstring(nStr:nStr) = ")"

       return
    End Subroutine FindFmt

    !!--++
    !!--++ Subroutine FindFMT_Err(aLine,nC_L)
    !!--++    character(len=*), intent(in) :: aLine   !  In -> Current data line
    !!--++    integer,          intent(in) :: nC_L    !  In -> location of last character treated
    !!--++
    !!--++    (PRIVATE)
    !!--++    Output the error messages from FindFMT
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine FindFMT_Err(aLine,nC_L)
       !---- Arguments ----!
       Character(len=*), intent(in) ::   aLine
       Integer,         intent (in) ::   nC_L

       !---- Local variables ----!
       Integer, parameter                             :: MssgBeg=-2   ! lower message number
       Integer, parameter                             :: MssgEnd=12   ! upper message number
       Character (len=48), dimension(MssgBeg:MssgEnd) :: Message=(/ &
                                                         "FindFMT: data line FORTRAN read error nber:     ",          &
                                                         "FindFMT: End of file !                          ",          &
                                                         "FindFMT: no error                               ",          &
                                                         "FindFMT: empty format descriptor (0 field) !    ",          &
                                                         "FindFMT: data string, read error !              ",          &
                                                         "FindFMT: integer field found real !             ",          &
                                                         "FindFMT: begged dot, sign or 'e' character !    ",          &
                                                         "FindFMT: invalid character in an integer field !",          &
                                                         "FindFMT: invalid field in format descriptor !   ",          &
                                                         "FindFMT: invalid character in a numeric field ! ",          &
                                                         "FindFMT: 0 character in current field !         ",          &
                                                         "FindFMT: format string length exceeded !        ",          &
                                                         "FindFMT: separator missing !                    ",          &
                                                         "FindFMT: incomplete E or D format !             ",          &
                                                         "FindFMT: incomplete number !                    "/)

       Integer                                         :: Ln, i
       Character (len=40)                              :: LaMarque

       !---- Error message ----!
       if (iErr_fmt == iErrNone .or. iErr_fmt == iErrEof) then
          Return
       else if (iErr_fmt < iErrEof) then
          Mess_FindFMT%nlines=1
          Write(unit=Mess_FindFMT%txt(1),fmt="(a,i4)") " "//Message(-2)(1:Len_trim(Message(-2))), -(iErr_fmt+100)
       else if (iErr_fmt < MssgBeg .or. iErr_fmt > MssgEnd) then
          Mess_FindFMT%nlines=1
          Write(unit=Mess_FindFMT%txt(1),fmt="(a,i2)") " FMT decode error number:",iErr_fmt
       else
          Mess_FindFMT%nlines=1
          Write(unit=Mess_FindFMT%txt(1),fmt="(a)") " "//Message(iErr_fmt)(1:Len_trim(Message(iErr_fmt)))
       end if

       !---- Output data line and print a mark at error location ----!
       Ln = max(Len_trim(aLine),1)
       if (Ln <= 129) then
          Mess_FindFMT%nlines=2
          Write(unit=Mess_FindFMT%txt(2),fmt="(tr1,a)") "'"//aLine(1:Ln)//"'"
          if (nC_L == 1) then
             Mess_FindFMT%nlines=3
             Write(unit=Mess_FindFMT%txt(3),fmt="(tr1,a)")  "  ^----"
          else if (nC_L > 1) then
             Write(unit=LaMarque,fmt="(a,i3,a)")  "(a,", nC_L, "a,a)"
             Mess_FindFMT%nlines=3
             write(unit=Mess_FindFMT%txt(3),fmt=LaMarque)  " ",("-",i=1,nC_L),"^"
          end if
       else
          Mess_FindFMT%nlines=2
          Write(unit=Mess_FindFMT%txt(2),fmt="(a)") " "//aLine(1:Ln)
          Write(unit=LaMarque,fmt="(a,i3,a)")  "(a,", nC_L-1, "a,a)"
          Mess_FindFMT%nlines=3
          Write(unit=Mess_FindFMT%txt(3),fmt=LaMarque) " ",("-",i=1,nC_L-1),"^"
       end if

       return
    End Subroutine FindFMT_Err

    !!----
    !!---- Subroutine Frac_Trans_1Dig(v,CharF)
    !!----    real(kind=sp), dimension(3), intent( in)   :: V     !In -> Vector: v(1)=0.25, v(2)=-0.4, v(3)=0.33333
    !!----    character (len=* ),          intent(out)   :: CharF ! Out -> String: "(1/4,-2/5,1/3)"
    !!----
    !!----    Subroutine returning a string of length=16 describing a
    !!----    3D translation vector written in fractional form as quotient
    !!----    of 1-digit integers with sign.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Frac_Trans_1Dig(v,CharF)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in)   :: v
       character (len=* ),          intent(out)   :: CharF

       !---- Local Variables ----!
       character (len=4), dimension(3)   :: Frac
       integer                           :: i,j

       CharF="(    ,    ,    )"
       do i=1,3
          call Get_Fraction_1Dig(v(i),Frac(i))
          j=index(Frac(i),"+")
          if (j /= 0) Frac(i)(j:j) = " "
       end do
       CharF(2:5)  =Frac(1)
       CharF(7:10) =Frac(2)
       CharF(12:15)=Frac(3)

       return
    End Subroutine Frac_Trans_1Dig

    !!----
    !!---- Subroutine Frac_Trans_2Dig(v,CharF)
    !!----    real(kind=sp), dimension(3), intent( in) :: V       !  In -> Vector: v(1)=0.3, v(2)=-0.4, v(3)=-5.5
    !!----    character (len=* ),          intent(out) :: CharF   ! Out -> String: "(3/10,-2/5,-11/2)"
    !!----
    !!----    Function returning a string of length=22 describing a
    !!----    3D translation vector written in fractional form as quotient
    !!----    of 2-digit integers with sign.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Frac_Trans_2Dig(v,CharF)
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: v
       character (len=* ),          intent(out) :: CharF

       !---- Local Variables ----!
       character (len=6), dimension(3) :: Frac
       character (len=22)              :: str
       integer                         :: i,j

       str="(      ,      ,      )"
       do i=1,3
          call Get_Fraction_2Dig(v(i),Frac(i))
          j=index(Frac(i),"+")
          if (j /= 0) Frac(i)(j:j) = " "
       end do
       str( 2: 7) =Frac(1)
       str( 9:14) =Frac(2)
       str(16:21) =Frac(3)
       CharF=Pack_String(str)

       return
    End Subroutine Frac_Trans_2Dig

    !!----
   !!---- Subroutine Get_Fraction_1Dig(V,Fracc)
    !!----    real(kind=sp),      intent( in) :: V       !  In -> Input real(kind=sp) number
    !!----    character (len=*),  intent(out) :: Fracc   ! Out -> Fracction in character form
    !!----
    !!----    Get a string with the most simple fraction that uses single digits
    !!----    in numerator and denominator. Used, for instance, to get a character
    !!----    representation of symmetry operators.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Fraction_1Dig(V,Fracc)
       !---- Argument ----!
       real(kind=sp),    intent( in) :: v
       character(len=*), intent(out) :: fracc

       !---- Local variables ----!
       integer          ::  numerator, denominator
       real(kind=sp)    ::  num, denom, frac

       fracc="**/*"
       if (Zbelong(v)) then
          fracc="    "
          if (v > 0.0) then
             write(unit=fracc, fmt="(a,i1)") "+", nint(v)
          else
             write(unit=fracc, fmt="(i2)") nint(v)
          end if
       else
          do numerator=1,9
             num=numerator
             do denominator=2,9
                denom=denominator
                frac=num/denom
                if (Negligible(frac-abs(v))) then
                   fracc="    "
                   if (v > 0.0) then
                      write(unit=fracc, fmt="(2(a,i1))") "+",numerator,"/",denominator
                   else
                      write(unit=fracc, fmt="(2(a,i1))") "-",numerator,"/",denominator
                   end if
                   return
                end if
             end do
          end do
       end if

       return
    End Subroutine Get_Fraction_1Dig

    !!----
    !!---- Subroutine Get_Fraction_2Dig(V,Fracc)
    !!----    real(kind=sp),      intent( in) :: V       !  In -> Input real(kind=sp) number
    !!----    character (len=*),  intent(out) :: Fracc   ! Out -> Fracction in character form
    !!----
    !!----    Get a string with the most simple fraction that uses up to two
    !!----    digits in numerator and denominator. Used, for instance, to get a
    !!----    character representation of symmetry operators.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Fraction_2Dig(v,fracc)
       !---- Argument ----!
       real(kind=sp),    intent( in) :: v
       character(len=*), intent(out) :: fracc

       !---- Local variables ----!
       character (len=16) :: formm
       real(kind=sp)      :: num, denom, frac
       integer            :: numerator, denominator

       fracc="***/**"
       if (Zbelong(v)) then
          fracc="      "
          if (v > 0.0_sp) then
             formm="(a,i1)"
             if(v >=10.0_sp) formm="(a,i2)"
             write(unit=fracc,fmt=formm) "+", nint(v)
          else
             formm="(i2)"
             if(v >=10.0_sp) formm="(i3)"
             write(unit=fracc,fmt=formm) nint(v)
          end if
       else
          do numerator=1,24
             num=numerator
             do denominator=2,24
                denom=denominator
                frac=num/denom
                if (Negligible(frac-abs(v))) then
                   fracc=" "
                   formm="(a1,i1,a1,i1)"
                   if(numerator >=10 .and. denominator <=  9) formm="(a1,i2,a1,i1)"
                   if(numerator >=10 .and. denominator >= 10) formm="(a1,i2,a1,i2)"
                   if(numerator <= 9 .and. denominator >= 10) formm="(a1,i1,a1,i2)"
                   if (v > 0.0_sp) then
                      write(unit=fracc,fmt=formm) "+",numerator,"/",denominator
                   else
                      write(unit=fracc,fmt=formm) "-",numerator,"/",denominator
                   end if
                   return
                end if
             end do
          end do
       end if

       return
    End Subroutine Get_Fraction_2Dig

    !!----
    !!---- Subroutine Get_LogUnit(lun)
    !!----   integer,     intent(out) :: lun !First logical unit available
    !!----
    !!----   Provides the number of the first logical unit that is not opened.
    !!----   Useful for getting a logical unit to a file that should be opened
    !!----   of the flight.
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Get_LogUnit(lun)
       !---- Arguments ----!
       integer,  intent(out) :: lun

       !---- Local variables ----!
       logical :: op
       integer, parameter :: max_iunits=500

       lun=1
       do
          inquire(unit=lun,opened=op)
          if (.not. op) exit
          lun=lun+1
          if (lun == max_iunits) then
             lun=-1
             exit
          end if
       end do

       return
    End Subroutine Get_LogUnit

    !!----
    !!---- Subroutine Getnum(Line, Vet, Ivet, Iv)
    !!----    character(len=*),              intent( in) :: Line    !  In -> Input String to convert
    !!----    real(kind=sp), dimension(:),   intent(out) :: Vet     ! Out -> Vector of real(kind=sp) numbers
    !!----    integer,dimension(:),          intent(out) :: Ivet    ! Out -> Vector of integer numbers
    !!----    integer,                       intent(out) :: Iv      ! Out -> Number of numbers in Vet/Ivet
    !!----
    !!----    Converts a string to numbers and write on VET/IVET if real/integer. Control
    !!----    of errors is possible by inquiring the global variables ERR_STRING and
    !!----    ERR_MESS_STRING
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Getnum(line,vet,ivet,iv)
       !---- Argument ----!
       character (len=*),          intent ( in) :: line
       real(kind=sp), dimension(:),intent (out) :: vet
       integer, dimension(:),      intent (out) :: ivet
       integer,                    intent (out) :: iv

       !---- Local variables ----!
       logical                   :: numero
       character (len=len(line)) :: resto,cifre
       integer                   :: i,isum,ncharl,nchard,isegno,iniz,ipoi,idec,idig
       integer                   :: nchart, npos,nchard1,isum_exp,ioper
       real(kind=sp)             :: suma,segno,dec
       real(kind=sp)             :: sum_m

       !---- Initializing variables ----!
       call init_err_string()
       iv=0
       ivet=0
       vet=0.0

       resto=u_case(line)

       do
          ioper=0
          isum_exp=0
          nchard1=0
          sum_m=0.0
          suma=0.0
          isum=0
          call cutst(resto,ncharl,cifre,nchard)
          if (nchard <= 0) exit

          !---- Is a number ----!
          numero=.true.
          do i=1,nchard
             if (cifre(i:i) =='E') cycle
             npos=index(digit,cifre(i:i))
             if (npos /= 0) cycle
             numero=.false.
          end do
          if (.not. numero) then
             err_string=.true.
             err_mess_string="The variable cannot be computed as a number in GETNUM "
             return
          end if

          !---- Positive or Negative number ----!
          segno=1.0
          isegno=1
          iniz=1
          if (cifre(1:1) == digit(12:12)) then
             segno=-1.0
             isegno=-1
             iniz=2
          end if

          !---- Decimal Number ----!
          ipoi=index(cifre(1:nchard),digit(11:11))

          !---- Exponential Number ----!
          nchard1=index(cifre(1:nchard),"E")
          if (nchard1 /= 0) then
             nchart=nchard
             nchard=nchard1-1
          end if

          if (ipoi == 0) ipoi=nchard+1
          dec=real(ipoi-1-iniz)
          idec=ipoi-1-iniz
          do i=iniz,nchard
             idig=index(digit,cifre(i:i))
             if (idig >= 1 .and. idig <= 11)  then
                if (idig <= 10)  then
                   suma=suma+real(idig-1)*10.0**dec
                   if (idec >= 0) isum=isum*10+(idig-1)
                   dec=dec-1.0
                   idec=idec-1
                end if
             else
                err_string=.true.
                err_mess_string="Limits of digit variable exceeded in GETNUM"
                return
             end if
          end do

          if (nchard1 /= 0) then
             nchard1=nchard1+1
             select case (cifre(nchard1:nchard1))
                case ("-")
                   ioper=1
                   nchard1=nchard1+1

                case ("+")
                   nchard1=nchard1+1
             end select

             do i=nchard1,nchart
                idig=index(digit,cifre(i:i))
                if (idig >= 1 .and. idig <= 10)  then
                   isum_exp=isum_exp*10+(idig-1)
                else
                   err_string=.true.
                   err_mess_string="Limits of digit variable exceeded in GETNUM"
                   return
                end if
             end do
          end if

          iv=iv+1
          vet(iv)=suma*segno
          ivet(iv)=isum*isegno

          if (nchard1 /= 0) then
             select case (ioper)
                case (0)
                   sum_m=10.0**isum_exp

                case (1)
                   sum_m=10.0**isum_exp
                   sum_m=1.0/sum_m
             end select
             vet(iv)=vet(iv)*sum_m
          end if

          if (ncharl <= 0) then
             exit
          end if
       end do

       return
    End Subroutine Getnum

    !!----
    !!---- Subroutine Getnum_Std(Line, Value, Std, Ic)
    !!----    character(len=*),            intent( in) :: Line    !  In -> Input String
    !!----    real(kind=sp), dimension(:), intent(out) :: Value   ! Out -> Vector of values with real(kind=sp) numbers
    !!----    real(kind=sp), dimension(:), intent(out) :: Std     ! Out -> Vector of standard deviation values
    !!----    integer,                     intent(out) :: Ic      ! Out -> Number of components of vector Value
    !!----
    !!----    Converts a string to a numbers with standard deviation with format: x.fffff(s)
    !!----    Control of errors is possible by inquiring the global variables ERR_STRING
    !!----    and ERR_MESS_STRING.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine GetNum_Std(line, value, std, ic)
       !----Arguments ----!
       character(len=*),             intent( in) :: line
       real(kind=sp), dimension(:),  intent(out) :: value
       real(kind=sp), dimension(:),  intent(out) :: std
       integer,                      intent(out) :: ic

       !---- Local Variables ----!
       character(len=len(line))               :: resto,dire,numm
       integer                                :: iv,nlong
       integer                                :: np, np1, np2
       integer, dimension(size(value))        :: ivet
       real(kind=sp), dimension(size(value))  :: vet

       value=0.0
       std  =0.0
       ic   =0
       call init_err_string()

       !---- Initial Checks ----!
       if (len_trim(line) == 0) then
          err_string=.true.
          err_mess_string="Blank line"
          return
       end if
       resto=adjustl(line)

       do
          if (len_trim(resto) == 0) exit
          call cutst(resto,nlong,dire)
          np1=index(dire,"(")
          np2=index(dire,")")

          if ( (np2 < np1) .or.               &  ! ")" before than "("
               (np1==0 .and. np2 >0) .or.     &  ! "(" don"t exists
               (np2==0 .and. np1 >0) ) then      ! ")" don"t exists
             err_string=.true.
             err_mess_string="Wrong format using Standard values"
             return
          end if

          if (np1 == 0 .and. np2 ==0) then
             call getnum(dire,vet,ivet,iv)
             if (iv /= 1 .or. err_string) then
                err_string=.true.
                err_mess_string="Bad format"
                return
             end if
             ic=ic+1
             value(ic)=vet(1)
          else
             numm=dire(1:np1-1)
             np=index(numm,".")
             if (np == 0) then
                call getnum(numm,vet,ivet,iv)
                if (iv /= 1 .or. err_string) then
                   err_string=.true.
                   err_mess_string="Bad format"
                   return
                end if
                ic=ic+1
                value(ic)=vet(1)
                numm=dire(np1+1:np2-1)
                call getnum(numm,vet,ivet,iv)
                if (iv /= 1) then
                   err_string=.true.
                   err_mess_string="Bad format"
                   return
                end if
                std(ic)=vet(1)
             else
                np=np1-np-1
                call getnum(numm,vet,ivet,iv)
                if (iv /= 1 .or. err_string) then
                   err_string=.true.
                   err_mess_string="Bad format"
                   return
                end if
                ic=ic+1
                value(ic)=vet(1)
                numm=dire(np1+1:np2-1)
                call getnum(numm,vet,ivet,iv)
                if (iv /= 1 .or. err_string) then
                   err_string=.true.
                   err_mess_string="Bad format"
                   return
                end if
                std(ic)=vet(1)/(10.0**np)
             end if
          end if
       end do

       return
    End Subroutine GetNum_Std

    !!----
    !!---- Subroutine Getword(Line, Dire, Ic)
    !!----    character(len=*),              intent( in) :: Line   !  In -> Input String
    !!----    character(len=*),dimension(:), intent(out) :: Dire   ! Out -> Vector of Words
    !!----    integer,                       intent(out) :: Ic     ! Out -> Number of words
    !!----
    !!----    Determines the number of words (Ic) in the string "Line" and generates a
    !!----    character vector "Dire" with separated words.
    !!----    Control of errors is possible by inquiring the global variables ERR_STRING
    !!----    and ERR_MESS_STRING.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Getword(line,dire,ic)
       !---- Argument ----!
       character (len=*),                 intent ( in) :: line
       character (len=*), dimension(:),   intent (out) :: dire
       integer,                           intent (out) :: ic

       !---- Local variables ----!
       character (len=len(line)) :: line1,line2
       integer                   :: nlong2
       integer                   :: ndim

       call init_err_string()
       ic=0
       ndim=size(dire)
       line1=line

       do
          line1=adjustl(line1)
          call cutst(line1,line2=line2,nlong2=nlong2)
          if (nlong2 == 0) exit
          ic=ic+1
          if (ic > ndim) then
             err_string=.true.
             err_mess_string="Dimension of DIRE exceeded"
             exit
          end if
          dire(ic)=line2(:nlong2)
       end do

       return
    End Subroutine Getword

    !!----
    !!---- Subroutine Init_Err_String()
    !!----
    !!----    Initializes general error variables for this module as:
    !!----    ERR_STRING=.false. ;  ERR_MESS_STRING=" "
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_String()

       err_string=.false.
       err_mess_string=" "

       return
    End Subroutine Init_Err_String

    !!----
    !!---- Subroutine Init_FindFMT(nline)
    !!----   integer, optional, intent(in) :: nline
    !!----
    !!----    Initializes the subroutine FindFMT.
    !!----    Mess_FindFMT (of type Err_Text_Type) is initialized to zero lines.
    !!----    Line_nb is initialized to zero (current line in the file),
    !!----    or Line_nb=line if the optional argument "line" is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_FindFMT(nline)
       !---- Arguments ----!
       integer, optional, intent(in) :: nline

       line_nb=0
       if(present(nline)) line_nb=nline
       Mess_FindFMT = Err_Text_Type(0,(/" "," "," "," "," "/))

       return
    End Subroutine Init_FindFMT

    !!----
    !!---- Subroutine Lcase(Line)
    !!----    character(len=*) :: Line
    !!----
    !!----    Conversion to lower case. Line is modified
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Lcase(line)
       !---- Argument ----!
       character (len=*), intent(in out) :: line

       line=l_case(line)

       return
    End Subroutine Lcase

    !!----
    !!---- Subroutine Number_Lines(Filename,n)
    !!----    character(len=*), intent(in) :: Filename     !  In -> Name of the file
    !!----    integer        , intent(out) :: N            ! Out -> Number of lines in the file
    !!----
    !!----    Return the number of lines contained in a file. If the file is
    !!----    open, a rewind procedure is made.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Number_Lines(filename,n)
       !---- Arguments ----!
       character(len=*), intent(in)  :: filename
       integer,          intent(out) :: n

       !---- Local Variables ----!
       logical            :: info
       integer            :: lun,cond

       !---- Init ----!
       info=.false.
       call get_logunit(lun)
       n=0
       cond=0

       !---- Exist filename ? ----!
       inquire (file=filename,exist=info)
       if (.not. info) return

       !---- Is it Open? ----!
       inquire (file=filename,opened=info)
       if (.not. info) then
          open(unit=lun,file=filename, status="old",action="read", position="rewind")
       else
          inquire(file=filename,number=lun)
          rewind(unit=lun)
       end if

       !---- Counting lines ----!
       do
          read(unit=lun,fmt="(a)",iostat=cond)
          if (cond /= 0) exit
          n=n+1
       end do

       !---- Was Opened? ----!
       if (.not. info) then
          close(unit=lun)
       else
          rewind(unit=lun)
       end if

       return
    End Subroutine Number_Lines

    !!----
    !!---- Subroutine Read_Key_Str(Filevar,Nline_Ini,Nline_End,Keyword,String)
    !!----    character(len=*),dimension(:), intent(in)      :: Filevar      !  In -> Input vector of String
    !!----    integer,                       intent(in out)  :: Nline_Ini    !  In -> Pointer to initial position to search
    !!----                                                                   ! Out -> Pointer to final position in search
    !!----    integer,                       intent(in)      :: Nline_End    !  In -> Pointer to final position to search
    !!----    character(len=*),              intent(in)      :: Keyword      !  In -> Word to search
    !!----    character(len=*),              intent(out)     :: String       ! Out -> Rest of the input string
    !!----
    !!----    Read a string on "filevar" starting with a particular "keyword" between lines "nline_ini" and
    !!----    "nline_end".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Key_Str(filevar,nline_ini,nline_end,keyword,string)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)      :: filevar
       integer,                        intent(in out)  :: nline_ini
       integer,                        intent(in)      :: nline_end
       character(len=*),               intent(in)      :: keyword
       character(len=*),               intent(out)     :: string

       !---- Local Variable ----!
       character(len=len(filevar(1))) :: line,linec
       character(len=len(keyword))    :: key
       integer                        :: i,np,nt

       !---- Initial value ----!
       nt=min(size(filevar),nline_end)
       string=" "
       key =adjustl(keyword)
       call lcase(key)

       do i=nline_ini,nt
          line=adjustl(filevar(i))
          if (len_trim(line) == 0 .or. line(1:1) == "!") cycle
          linec=line
          call lcase(line)
          np=index(line,key)
          if (np == 0) cycle
          linec=linec(np:)
          call cutst(linec)
          string=linec
          nline_ini=i
          exit
       end do

       return
    End Subroutine Read_Key_Str

    !!----
    !!---- Subroutine Read_Key_Strval(Filevar,Nline_Ini,Nline_End,Keyword,String,Vet,Ivet,Iv)
    !!----    character(len=*),dimension(:),          intent(in)      :: Filevar      !  In -> Input vector of String
    !!----    integer,                                intent(in out)  :: Nline_Ini    !  In -> Pointer to initial position to search
    !!----                                                                            ! Out -> Pointer to final position in search
    !!----    integer,                                intent(in)      :: Nline_End    !  In -> Pointer to final position to search
    !!----    character(len=*),                       intent(in)      :: Keyword      !  In -> Word to search
    !!----    character(len=*),                       intent(out)     :: String       ! Out -> Rest of the input string
    !!----    real(kind=sp),dimension(:),   optional, intent(out)     :: Vet          ! Out -> Vector for real(kind=sp) numbers
    !!----    integer,dimension(:),         optional  intent(out)     :: Ivet         ! Out -> Vector for integer numbers
    !!----    integer,                      optional, intent(out)     :: Iv           ! Out -> Number of numbers
    !!----
    !!----    Read a string on "filevar" starting with a particular "keyword" between lines "nline_ini" and
    !!----    "nline_end". If the string contains numbers they are read and put into "vet/ivet". The variable
    !!----    "string" contains the input string without the "keyword".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Key_StrVal(filevar,nline_ini,nline_end,keyword,string,vet,ivet,iv)
       !---- Arguments ----!
       character(len=*), dimension(:),           intent(in)      :: filevar
       integer,                                  intent(in out)  :: nline_ini
       integer,                                  intent(in)      :: nline_end
       character(len=*),                         intent(in)      :: keyword
       character(len=*),                         intent(out)     :: string
       real(kind=sp),dimension(:),     optional, intent(out)     :: vet
       integer,dimension(:),           optional, intent(out)     :: ivet
       integer,                        optional, intent(out)     :: iv

       !---- Local Variable ----!
       logical                        :: sval
       character(len=len(filevar(1))) :: line,linec
       character(len=len(keyword))    :: key
       integer                        :: i,np,nt

       !---- Initial value ----!
       nt=min(size(filevar),nline_end)
       string=" "
       key =adjustl(keyword)
       call lcase(key)
       sval=.false.
       if (present(vet) .and. present(ivet) .and. present(iv)) sval=.true.
       if (sval) then
          vet=0.0
         ivet=0
           iv=0
       end if

       do i=nline_ini,nt
          line=adjustl(filevar(i))
          if (len_trim(line) == 0) cycle
          linec=line
          call lcase(line)
          np=index(line,key)
          if (np == 0) cycle
          linec=linec(np:)
          call cutst(linec)
          string=linec
          nline_ini=i
          exit
       end do

       if (sval .and. (len_trim(string) > 0) ) then
          line=string

          !---- String Value ----!
          call cutst(line,np,string)

          !---- Values ----!
          call getnum(line,vet,ivet,iv)
          if (iv <=0) then
              vet=0.0
             ivet=0
          end if
       end if

       return
    End Subroutine Read_Key_StrVal

    !!----
    !!---- Subroutine Read_Key_Value(Filevar,Nline_Ini,Nline_End,Keyword,Vet,Ivet,Iv)
    !!----    character(len=*),dimension(:), intent(in)      :: Filevar     !  In -> Input vector of String
    !!----    integer,                       intent(in out)  :: Nline_Ini   !  In -> Pointer to initial position to search
    !!----                                                                  ! Out -> Pointer to final position in search
    !!----    integer,                       intent(in)      :: Nline_End   !  In -> Pointer to final position to search
    !!----    character(len=*),              intent(in)      :: Keyword     !  In -> Word to search
    !!----    real(kind=sp),dimension(:),    intent(out)     :: Vet         ! Out -> Vector for real(kind=sp) numbers
    !!----    integer,dimension(:),          intent(out)     :: Ivet        ! Out -> Vector for integer numbers
    !!----    integer,                       intent(out)     :: Iv          ! Out -> Number of components
    !!----
    !!----    Read a string on "filevar" starting with a particular "keyword" between lines "nline_ini" and
    !!----    "nline_end". If the string contains numbers they are read and put into "vet/ivet".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Key_Value(filevar,nline_ini,nline_end,keyword,vet,ivet,iv)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)     :: filevar
       integer,                        intent(in out) :: nline_ini
       integer,                        intent(in)     :: nline_end
       character(len=*),               intent(in)     :: keyword
       real(kind=sp),dimension(:),     intent(out)    :: vet
       integer,dimension(:),           intent(out)    :: ivet
       integer,                        intent(out)    :: iv

       !---- Local Variable ----!
       character(len=len(filevar(1))) :: line
       character(len=len(keyword))    :: key
       integer                        :: i,np,nt

       !---- Initial value ----!
       nt=min(size(filevar),nline_end)
       iv  = 0
       vet = 0.0
       ivet= 0
       key =adjustl(keyword)
       call lcase(key)

       do i=nline_ini,nt
          np=0
          line=adjustl(filevar(i))
          if (len_trim(line) == 0 .or. line(1:1) == "!") cycle
          call lcase(line)
          np=index(line,key)
          if (np == 0) cycle
          line=line(np:)
          call cutst(line)
          call getnum(line,vet,ivet,iv)
          if (err_string) exit
          nline_ini=i
          exit
       end do

       return
    End Subroutine Read_Key_Value

    !!----
    !!---- Subroutine Read_Key_Valuest(Filevar,Nline_Ini,Nline_End,Keyword,Vet1,Vet2,Iv)
    !!----    character(len=*),dimension(:),  intent(in)     :: Filevar      !  In -> Input vector of String
    !!----    integer,                        intent(in out) :: Nline_Ini    !  In -> Pointer to initial position to search
    !!----                                                                   ! Out -> Pointer to final position in search
    !!----    integer,                        intent(in)     :: Nline_End    !  In -> Pointer to final position to search
    !!----    character(len=*),               intent(in)     :: Keyword      !  In -> Word to search
    !!----    real(kind=sp),dimension(:),     intent(out)    :: Vet1         ! Out -> Vector of real(kind=sp) numbers
    !!----    real(kind=sp),dimension(:),     intent(out)    :: Vet2         ! Out -> Vector of standard deviations
    !!----    integer,                        intent(out)    :: Iv           ! Out -> Number of components
    !!----
    !!----    Read parameters and standard deviation on the line of "filevar" starting with a particular "keyword".
    !!----    The search is done between lines "nline_ini" and "nline_end".
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Key_ValueST(filevar,nline_ini,nline_end,keyword,vet1,vet2,iv)
       !---- Arguments ----!
       character(len=*), dimension(:),  intent(in)     :: filevar
       integer,                         intent(in out) :: nline_ini
       integer,                         intent(in)     :: nline_end
       character(len=*),                intent(in)     :: keyword
       real(kind=sp),dimension(:),      intent(out)    :: vet1
       real(kind=sp),dimension(:),      intent(out)    :: vet2
       integer,                         intent(out)    :: iv

       !---- Local Variable ----!
       character(len=len(filevar(1))) :: line
       character(len=len(keyword))    :: key
       integer                        :: i,np,nt

       !---- Initial value ----!
       nt=min(size(filevar),nline_end)
       iv  = 0
       vet1 = 0.0
       vet2 = 0.0
       key =adjustl(keyword)
       call lcase(key)

       do i=nline_ini,nt
          line=adjustl(filevar(i))
          if (len_trim(line) == 0 .or. line(1:1) == "!") cycle
          call lcase(line)
          np=index(line,key)
          if (np == 0) cycle
          line=line(np:)
          call cutst(line)
          call getnum_std(line,vet1,vet2,iv)
          if (err_string) exit
          nline_ini=i
          exit
       end do

       return
    End Subroutine Read_Key_ValueST

    !!----
    !!---- Subroutine Reading_Lines(Filename,Nlines,Filevar)
    !!----    character(len= *), intent(in)                :: Filename   !  In -> Filename
    !!----    integer,           intent(in)                :: Nlines     !  In -> Number of lines to read
    !!----    character(len= *), dimension(:), intent(out) :: Filevar    ! Out -> String vector
    !!----
    !!----    Read nlines of the file and put the information on Filevar. If the file
    !!----    is open, the a rewind procedure is made.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Reading_Lines(filename,nlines,filevar)
       !---- Arguments ----!
       character(len=*),               intent( in) :: filename
       integer,                        intent( in) :: nlines
       character(len=*), dimension(:), intent(out) :: filevar

       !---- Local Variables ----!
       logical :: info
       integer :: lun,i

       !---- Init ----!
       call init_err_string()
       info=.false.
       call get_logunit(lun)

       !---- Exist filename ? ----!
       inquire (file=filename,exist=info)
       if (.not. info) then
          err_string=.true.
          err_mess_string="Not exist the file"
          return
       end if

       !---- Is it Open? ----!
       inquire (file=filename,opened=info)
       if (.not. info) then
          open(unit=lun,file=filename, status="old",action="read", position="rewind")
       else
          inquire(file=filename,number=lun)
          rewind(unit=lun)
       end if

       !---- Reading... ----!
       do i=1,nlines
          read(unit=lun,fmt="(a)") filevar(i)
       end do

       !---- Was Opened? ----!
       if (.not. info) then
          close(unit=lun)
       else
          rewind(unit=lun)
       end if

       return
    End Subroutine Reading_Lines

    !!----
    !!----
    !!---- Subroutine SetNum_Std(Value,Std,Line)
    !!----    real(kind=sp),            intent(in)  :: Value
    !!----    real(kind=sp),            intent(in)  :: Std
    !!----    character(len=*),intent (out):: Line
    !!----
    !!----    String with real(kind=sp) value and standar deviation
    !!----    between parenthesis
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine SetNum_Std(Value, Std, Line)
       !---- Argument ----!
       real(kind=sp),   intent(in)  :: Value
       real(kind=sp),   intent(in)  :: Std
       character(len=*),intent (out):: Line

       !---- Local Variables ----!
       character(len=10) :: fmtcar
       integer           :: n,np,iy
       real(kind=sp)     :: y

       if (abs(std) < 0.0000001) then
          if (abs(value) > 999999.0) then
             write(unit=line,fmt=*) value
          else
             write(unit=line,fmt="(f14.5)") value
          end if
          line=adjustl(line)
          if (line(1:1) /= "-") line=" "//trim(line)
          return
       end if

       np=0
       y=std
       do
          if (y >= 2.0) exit
          np=np+1
          y=y*10.0
       end do
       iy=nint(y)

       write(unit=line,fmt=*) value
       n=len_trim(line)
       fmtcar="f"
       if (n < 10) then
          write(unit=fmtcar(2:2),fmt="(i1)") n
       else
          write(unit=fmtcar(2:3),fmt="(i2)") n
       end if

       fmtcar=trim(fmtcar)//"."
       n=len_trim(fmtcar)
       if (np < 10) then
          write(unit=fmtcar(n+1:),fmt="(i1)") np
       else
          write(unit=fmtcar(n+1:),fmt="(i2)") np
       end if
       fmtcar="("//trim(fmtcar)//")"

       write(unit=line,fmt=fmtcar) value
       line=trim(line)
       n=len_trim(line)
       if (line(n:n) == ".") then
          line(n:n)=" "
       end if
       line=trim(line)//"("
       n=len_trim(line)
       write(unit=line(n+1:),fmt=*) iy
       line=pack_string(line)
       line=trim(line)//")"
       if(line(1:1) /= "-") line=" "//trim(line)

       return
    End Subroutine SetNum_Std

    !!--++
    !!--++ Subroutine SGetFTMfield(GetFTMfield,FMTfields,nFld,nFldMax)
    !!--++    Integer ,          intent(out)    ::  GetFTMfield
    !!--++    Character (len=*) ,intent( in)    ::  FMTfields     !  -> format descriptor
    !!--++    Integer ,          intent(in out) ::  nFld          ! <-> current field in format descriptor
    !!--++    Integer ,          intent( in)    ::  nFldMax       !  -> max. number of fields in format descriptor
    !!--++
    !!--++    (PRIVATE)
    !!--++    Get current field type
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine SGetFTMfield(GetFTMfield,FMTfields,nFld,nFldMax)
       !---- Arguments ----!
       Character (len=*) ,intent( in)    ::  FMTfields
       Integer ,          intent(in out) ::  nFld
       Integer ,          intent( in)    ::  nFldMax
       Integer ,          intent(out)    ::  GetFTMfield

       !---- Local variables ----!
       character (len=1) ::  Car

       nFld = nFld + 1
       if (nFld > nFldMax) then
          GetFTMfield = iEndFMT
       else
          Car = FMTfields(nFld:nFld)
          if (Car == "I") then
             GetFTMfield = iInte
          else if (Car == "F") then
             GetFTMfield = iReal
          else if (iChar(Car) >= i_One .and. iChar(Car) <= i_Nine) then
             GetFTMfield = iChar(Car)
          else
             GetFTMfield = iEndFMT
             iErr_fmt = iErrInvalField         ! Error in field definition
          end if
       end if

       return
    End Subroutine SGetFTMfield

    !!--++
    !!--++ Subroutine TreatMCharField(iFld,aLine,L_Line,nC_L,nC_X)
    !!--++    Integer,          intent(in out)  :: iFld   ! <-> "A" format size (1 to 9)
    !!--++    Character(len=*), intent(in)      :: aLine  !  -> data line to be analysed
    !!--++    Integer,          intent(in)      :: L_Line !  -> true length of data Line
    !!--++    Integer,          intent(in out)  :: nC_L   ! <-> current character in data line
    !!--++    Integer,          intent(out)     :: nC_X   ! <-  number of characters in X format field (now nx -> trn)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Fixed length "A1 to A9" field : A<iFld-48>
    !!--++    Leading spaces are ignored; separators : space and Tab
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine TreatMCharField(iFld,aLine,L_Line,nC_L,nC_X)
       !---- Arguments ----!
       Integer,           intent(in out)  :: iFld
       Character (len=*), intent(in)      :: aLine
       Integer,           intent(in)      :: L_Line
       Integer,           intent(in out)  :: nC_L
       Integer,           intent(out)     :: nC_X

       !---- Local variables ----!
       Character (len=1) ::   Car
       Integer           ::   nCar
       Logical           ::   ifEnd

       nC_X = 0
       iErr_fmt = 0

       !---- End of ligne ----!
       if (nC_L >= L_Line) return

       !---- if not 1rst field, 1rst character must be a separator ----!
       if (nC_L > 1) Then
          nC_L = nC_L+1
          Car  = aLine(nC_L:nC_L)
          if (Car /= " " .and. Car /= cTab) then
             iErr_fmt = iErrSepMiss              ! separator missing
             return
          end if
          nC_X = nC_X+1
       end if

       !---- Remove leading spaces ----!
       ifEnd = .false.
       do
          if (ifEnd) exit
          if (nC_L >= L_Line) return        ! end of line
          nC_L = nC_L+1
          Car  = aLine(nC_L:nC_L)
          if (Car == " ") then
             nC_X = nC_X+1                   ! count leading spaces
          else
             ifEnd = .true.                  ! 1rst valid character
             nC_L = nC_L-1
          end if
       end do

       !---- Count characters until next Tab or end of line ----!
       nCar = 0
       ifEnd = .false.
       do
          if (ifEnd) exit
          if (nC_L < L_Line .and. nCar < (iFld-48)) then
             nC_L = nC_L+1
             nCar = nCar+1
             Car = aLine(nC_L:nC_L)
             if (Car == " " .or. Car == cTab) then
                ifEnd = .true.                ! separator found
                nCar  = nCar - 1
                nC_L  = nC_L - 1
             end if
          else
             ifEnd = .true.                  ! end of line
          end if
       end do

       !---- Load size of the A format field ----!
       if (nCar == 0) then
          iErr_fmt = iErrEmptyField             ! no charac. in field
       else
          iFld = nCar+48                    ! true size of the A field
       end if

       return
    End Subroutine TreatMCharField

    !!--++
    !!--++ Subroutine TreatNumerField(iFld,aLine,L_Line,nC_L,nCar)
    !!--++    Integer ,          intent( in)    ::  iFld   !  -> field type
    !!--++    Character (len=*), intent(in out) ::  aLine  ! <-> data line
    !!--++    Integer ,          intent( in)    ::  L_Line !  -> true length of the data line
    !!--++    Integer ,          intent(in out) ::  nC_L   ! <-> counts characters in data line
    !!--++    Integer ,          intent(in out) ::  nCar   ! <-> counts characters in format field
    !!--++
    !!--++    (PRIVATE)
    !!--++    Free "I" and "F" formats
    !!--++    Look for a separator (space or Tab) after any valid character
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine TreatNumerField(iFld,aLine,L_Line,nC_L,nCar)
       !---- Arguments ----!
       Integer ,          intent( in)    ::  iFld   ! field type
       Character (len=*), intent(in out) ::  aLine
       Integer ,          intent( in)    ::  L_Line ! true length of the data line
       Integer ,          intent(in out) ::  nC_L   ! counts characters in data line
       Integer ,          intent(in out) ::  nCar   ! counts characters in format field

       !---- Local variables ----!
       Character (len=1)   ::  Car,Car_
       Integer             ::  nCar1                ! 1st usefull character in field
       Integer             ::  nPosi                ! number of 1st character in field
       Logical             ::  ifEnd,ifDot,ifSign

       iErr_fmt   = 0
       nCar   = 0
       ifDot  = .false.
       ifSign = .false.
       nPosi  = nC_L

       !---- Skip previous separator (space, Tab or sign) and leading spaces ----!
       ifEnd = .false.
       do
          if (ifEnd) exit
          nC_L = nC_L+1
          if (nC_L <= L_Line) then
             nCar = nCar+1
             Car = aLine(nC_L:nC_L)

             !---- Tab character ----!
             if (Car == cTab) Then
                if (nCar == 1 .and. nC_L > 1) then
                   aLine(nC_L:nC_L) = " "      ! previous separator
                else
                   if (ifSign) then
                      iErr_fmt = iErrNumber         ! incomplete number
                      return
                   end if
                   nC_L = nC_L-1               ! new separator
                   nCar = nCar-1
                   return
                end if

             else if (Car == "+" .or. Car == "-") then
                !---- a sign ----!
                ifSign = .true.

             else if (Car == " ") then
                !---- a space ----!
                if (ifSign) then
                   iErr_fmt = iErrNumber           ! incomplete number
                   return
                end if

             else
                !---- any other character ----!
                ifEnd = .true.
             end if
          else
             return                          ! end of line
          end if
       end do

       !---- No valid previous separator found (Except for 1st field) ----!
       if (nPosi > 1 .and. nCar == 1) then
          iErr_fmt = iErrSepMiss                ! separator missing
          return
       end if

       !---- Check first character and initialize search ----!

       !---- Decimal point -> valid in real fields only ----!
       if (Car == ".") then
          ifDot = .true.
          if (iFld /= iReal)  then
             iErr_fmt = iErrFieldType            ! not an integer field
             Return
          end if

       else if(Car == "E".or.Car == "e".or.Car == "d".or.Car == "D") then
          !---- e,E,d,D -> always invalid at this position ----!
          if (iFld == iReal) then
             iErr_fmt = iErrEfrmt                ! incomplete E or D format
          else
             iErr_fmt = iErrInvalC               ! invalid char in int. field
          end if
          return

       else if (iChar(Car) < i_Zero .or. iChar(Car) > i_Nine) then
          !---- invalid if not a sign or a digit ----!
          iErr_fmt = iErrInvalChar        ! invalid character
          return
       end if

       !---- save position of first character ----!
       nCar1 = nCar

       !---- Count characters in number ----!
       ifEnd = .false.

       do
          if (ifEnd) exit
          Car_ = Car      ! save previous character
          nC_L = nC_L+1
          if (nC_L <= L_Line) then
             nCar = nCar+1
             Car = aLine(nC_L:nC_L)

             !---- Current character is a decimal point ----!
             if (Car == ".") then
                if (ifDot) then
                   iErr_fmt = iErrCharBegg         ! begged character (dot)
                   Return
                else if (iFld /= iReal) then
                   iErr_fmt = iErrFieldType        ! not an integer field
                   Return
                else
                   ifDot = .true.
                end if

             else if (Car == " " .or. Car == cTab) then
                !---- Current character is a space or Tab (separator) ----!
                if (Car_ == "+" .or. Car_ == "-") then
                   iErr_fmt = iErrNumber           ! incomplete number
                   return
                end if
                ifEnd = .true.
                nCar  = nCar - 1
                nC_L  = nC_L - 1

             else if (Car == "+" .or. Car == "-") then
                !---- Current character is a sign ----!
                if (Car_ == "+" .or. Car_ == "-") then
                   iErr_fmt = iErrCharBegg         ! begged character
                   return
                else if (nCar > nCar1) then
                   if (iFld == iReal) then
                      if (Car_ /= "E" .and. Car_ /= "e" .and. Car_ /= "D" .and. Car_ /= "d") then
                         ifEnd = .true.          ! Sign is a valid separator
                         nCar  = nCar - 1
                         nC_L  = nC_L - 1
                         Return
                      end if
                   else                        ! Sign is a valid separator
                      ifEnd = .true.
                      nCar  = nCar - 1
                      nC_L  = nC_L - 1
                      Return
                   end if
                end if

             else if (Car == "E" .or. Car == "e" .or. Car == "d" .or. Car == "D") then
                !---- Current character is a "e E d D" ----!
                if (nCar == nCar1 .or. Car_ == "-" .or. Car_ == "+") then
                   iErr_fmt = iErrEfrmt            ! incomplete E or D format
                   return
                else if (Car_ == Car) then
                   iErr_fmt = iErrCharBegg         ! begged character
                   return
                end if

             else if (iChar(Car) < i_Zero .or. iChar(Car) > i_Nine) then
                !---- Ccurrent character is not a valid one ? ----!
                iErr_fmt = iErrInvalChar          ! invalid character
                Return
             end if
          else
             ifEnd = .true.                  ! end of line
          end if
       end do

       return
    End Subroutine TreatNumerField

    !!----
    !!---- Subroutine Ucase(Line)
    !!----    character(len=*) :: Line
    !!----
    !!----    Conversion to upper case. Line is modified
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Ucase(line)
       !---- Argument ----!
       character (len=*), intent(in out) :: line

       line=u_case(line)

       return
    End Subroutine Ucase

 End Module String_Utilities
