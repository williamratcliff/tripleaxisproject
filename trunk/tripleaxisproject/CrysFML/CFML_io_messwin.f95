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
!!--.. WINTERACTER ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: February - 2005
!!----            June    - 1999   Updated by JGP
!!----
!!---- DEPENDENCIES
!!--++    Winteracter or X/Winteracter Library
!!----
!!---- VARIABLES
!!--++    ICWINDOW                     [Private]
!!--++    WSCROLL                      [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       CLOSE_SCROLL_WINDOW
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       QUESTION_MESSAGE
!!----       WARNING_MESSAGE
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module IO_Messages
    !---- Use Modules ----!
    use Winteracter, only: YesNo,OKOnly,CommonYes, CommonOK, Modeless, ViewOnly,         &
                           WordWrap, NoMenu, NoToolbar, SystemFixed, EditTextLength,     &
                           ScreenHeight, InformationIcon, ExclamationIcon, QuestionIcon, &
                           WMessageBox, WindowCloseChild, WindowOpenChild, WEditFile,    &
                           WEditPutTextPart, WindowSelect, WInfoEditor, WInfoScreen

    !---- Definitions ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: close_scroll_window, error_message, info_message, question_message, warning_message, &
              write_scroll_text

    !---- Definitions ----!
    !!--++
    !!--++ ICWINDOW
    !!--++    integer, private :: icwindow
    !!--++
    !!--++    Code number for Scroll Window
    !!--++
    !!--++ Update: March - 2005
    !!
    integer, private :: icwindow= -1

    !!--++
    !!--++ WSCROLL
    !!--++    logical, private :: wscroll
    !!--++
    !!--++    Logical variable to indicate if the Scroll Window is
    !!--++    active or not.
    !!--++
    !!--++ Update: March - 2005
    !!
    logical, private :: wscroll = .false.

 Contains
    !!----
    !!---- SUBROUTINE CLOSE_SCROLL_WINDOW()
    !!----
    !!----    Close the Scroll Window
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Close_Scroll_Window()

       if (wscroll) call WindowCloseChild(icwindow)
       wscroll=.false.

       return
    End Subroutine Close_Scroll_Window

    !!----
    !!---- Subroutine Error_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an error message on the screen and in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(line,iunit)
       !---- Arguments ----!
       character(len=*), intent(in) :: line
       integer, optional,intent(in) :: iunit

       call WMessageBox(OKOnly, ExclamationIcon, CommonOK, line,"Error Message")

       if (present(iunit)) then
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") "**** ERROR: "//line
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") " "
       end if

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

       if (present(scroll_window) ) then
         !write on a file to be edited by Winteracter
       else
         call WMessageBox(OKOnly, InformationIcon, CommonOK, line,"Information Message")
       end if
       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") " "
          write(unit=iunit,fmt="(1x,a)") " "//line
          write(unit=iunit,fmt="(1x,a)") " "
       end if

       return
    End Subroutine Info_Message

    !!----
    !!---- SUBROUTINE Question_Message(line)
    !!----    character(len=*)  :: line
    !!----
    !!----    Print an question on the screen
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Question_Message(line)
       !---- Argument ----!
       character (len=*), intent(in) :: line

       call WMessageBox(YesNo,QuestionIcon,CommonYes,line,"Question")

       return
    End Subroutine Question_Message

    !!----
    !!---- SUBROUTINE WARNING_MESSAGE(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Info information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an message on the screen or in "Iunit" if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Warning_Message(line, iunit)
       !---- Arguments ----!
       character(len=*), intent(in) :: line
       integer, optional,intent(in) :: iunit

       call WMessageBox(OKOnly,ExclamationIcon,CommonOK, line,"Warning Message")

       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") "****"
          write(unit=iunit,fmt="(1x,a)") "**** WARNING: "//line
          write(unit=iunit,fmt="(1x,a)") "****"
       end if

       return
    End Subroutine Warning_Message

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

       !---- Local variables ----!
       character(len=2), parameter :: newline = char(13)//char(10)
       integer                     :: iendbuf

       !---- Open the Scroll Window if necessary ----!
       if (.not. wscroll) then
          call WindowOpenChild(icwindow,height=nint(WInfoScreen(ScreenHeight)*0.5), &
                               title='Info Window')
          call WEditFile(" ",IFlags=ViewOnly+NoMenu+NoToolbar,IFont=SystemFixed)
          wscroll=.true.
       end if

       iendbuf=WInfoEditor(icwindow,EditTextLength)+1
       call WindowSelect(icwindow)
       call WEditPutTextPart(trim(line)//newline,iendbuf)
       call WindowSelect(0)

       return
    End Subroutine Write_Scroll_Text

 End Module IO_Messages
