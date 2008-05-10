!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: IO_MESSAGES
!!----   INFO: Input / Output General Messages. It is convenient to use these intermediate procedures instead of
!!----         Fortran Write(*,*) or Print*, because it is much more simple to modify a program for making a GUI.
!!----         Usually GUI tools and libraries need special calls to dialog boxes for screen messages. These
!!----         calls may be implemented within this module using the same name procedures. The subroutines
!!----         ERROR_MESSAGE and INFO_MESSAGE are just wrappers for the actual calls.
!!--..
!!--..         NON-GRAPHICS ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: February - 2005
!!----            June    - 1999   Updated by JGP
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       PRINT_MESS
!!----       WAIT
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module IO_Messages
    !---- Use Modules ----!

    !---- Definitions ----!
    implicit none

    !---- List of public subroutines ----!
    public :: Info_Message, Error_Message, Print_Mess, Wait, write_scroll_text


 Contains

    !!----
    !!---- Subroutine Error_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an error message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(line,iunit)
       !---- Arguments ----!
       character(len=*), intent(in)           :: line
       integer,          intent(in), optional :: iunit

       !---- Local Variables ----!
       integer :: lun

       lun=6
       if (present(iunit)) lun=iunit

       write(unit=lun,fmt="(a)") " ****"
       write(unit=lun,fmt="(a)") " **** ERROR: "//line
       write(unit=lun,fmt="(a)") " ****"
       write(unit=lun,fmt="(a)") "  "

       return
    End Subroutine Error_Message

    !!----
    !!---- Subroutine Info_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Info_Message(line, iunit, scroll_window)
       !---- Arguments ----!
       character(len=*), intent(in)           :: line
       integer,          intent(in), optional :: iunit
       integer,          intent(in), optional :: scroll_window

       !---- Local Variables ----!
       integer :: lun

       lun=6
       if (present(iunit)) lun=iunit
       if (present(scroll_window)) lun=6
       write(unit=lun,fmt="(a)") "  "//line

       return
    End Subroutine Info_Message

    !!----
    !!---- Subroutine Print_Mess(Warning)
    !!----    character(len=*), intent(in)  :: Warning    !  In -> Print information
    !!----
    !!----    Print an message on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Print_Mess(Warning)
       !---- Arguments ----!
       character(len=*),intent(in) ::  Warning

       !---- Local Variables ----!
       integer :: lon

       lon=len_trim(Warning)
       if (lon == 0) then
          write(unit=*,fmt="(a)") "  "
       else
          if (warning(1:1) == "=" .or. warning(2:2) == "=") then
             write(unit=*,fmt="(a)") warning(1:lon)
          else
             write(unit=*,fmt="(a,a)")" =>", warning(1:lon)
          end if
       end if

       return
    End Subroutine Print_Mess

    !!----
    !!---- Subroutine Wait(Message)
    !!----    character(len=*), optional, intent(in) :: Message
    !!----
    !!----    Similar to Pause for Console version
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Wait(Message)
       !---- Argument ----!
       character(len=*), optional, intent(in) :: Message

       !---- Local variable ----!
       character(len=1) :: car

       write(unit=*,fmt="(a)") " "
       if (present(message)) write(unit=*,fmt="(a)", advance="no") message
       read(unit=*,fmt="(a)") car

       return
    End Subroutine Wait

    !!----
    !!---- SUBROUTINE WRITE_SCROLL_TEXT(Line)
    !!----    character(len=*), intent(in)           :: Line
    !!----
    !!----    Print the string in a the scroll window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Write_Scroll_Text(Line)
       !---- Argument ----!
       character(len=*), intent(in) :: line

       write(unit=*, fmt="(a)") trim(line)

       return
    End Subroutine Write_Scroll_Text

 End Module IO_Messages
