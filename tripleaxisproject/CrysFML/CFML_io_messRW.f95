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
!!----
!!--..
!!--..         REALWIN ZONE
!!--..
!!---- HISTORY
!!----
!!----    Update: February - 2005
!!----            June    - 1999   Updated by JGP
!!----
!!---- DEPENDENCIES
!!--++    RealWin Library
!!----
!!---- VARIABLES
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       ERROR_MESSAGE
!!----       INFO_MESSAGE
!!----       WRITE_SCROLL_TEXT
!!----
!!
 Module IO_Messages
    !---- Use Modules ----!
    use Realwin, only: message_box, write_scroll_line,SCROLL_TEXT, create_window,select_font

    !---- Definitions ----!
    implicit none

    !---- List private variables ----!
    integer, private :: messval

    !---- List of public subroutines ----!
    public :: error_message, info_message, write_scroll_text

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
    !!---- Subroutine Error_Message(Line, Iunit)
    !!----    character(len=*), intent(in)           :: Line    !  In -> Error information
    !!----    integer,          intent(in), optional :: Iunit   !  In -> Write information on Iunit unit
    !!----
    !!----    Print an error message on the screen and in 'Iunit' if present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Error_Message(line, iunit)
       !---- Arguments ----!
       character(len=*), intent(in)  :: line
       integer, optional, intent(in) :: iunit

       messval=message_box(text="Error: "//trim(line),title="Stop/Warning Box")  !rw

       if (present(iunit)) then
          write(unit=iunit,fmt="(1x,a)") "***"
          write(unit=iunit,fmt="(1x,a)") "*** ERROR: "//line
          write(unit=iunit,fmt="(1x,a)") "***"
          write(unit=iunit,fmt="(1x,a)") " "
       end if

       return
    End Subroutine Error_Message

    !!----
    !!---- Subroutine Info_Message(line, iunit, scroll_window)
    !!----    character(len=*),  intent(in) :: line           !  In -> Info information
    !!----    integer, optional, intent(in) :: iunit          !  In -> Write information on Iunit unit
    !!----    integer, optional, intent(in) :: scroll_window  !  In -> Write information on scroll windows
    !!----
    !!----    Print an message on the screen and in 'Iunit' if present
    !!----
    !!---- Update: February - 2005
    !!

    Subroutine Info_Message(line, iunit, scroll_window)
       character(len=*), intent(in)           :: line
       integer,          intent(in), optional :: iunit
       integer,          intent(in), optional :: scroll_window

       if(present(scroll_window)) then
         call write_scroll_line(scroll_window,text=trim(line))
       else
         messval=message_box(text=trim(line),title="Info Message")  !rw
       end if
       if (present(iunit) ) then
          write(unit=iunit,fmt="(1x,a)") " "
          write(unit=iunit,fmt="(1x,a)") " "//trim(line)
          write(unit=iunit,fmt="(1x,a)") " "
       end if
       return
    End Subroutine Info_Message

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


       !---- Open the Scroll Window if necessary ----!
       if (.not. wscroll) then
          icwindow = create_window(window_name="Scroll Text Window",x=0.15,y=0.07, &
                            width= 0.7, height=0.42, &
                            paint_code=SCROLL_TEXT,&
                            text_font=select_font(typeface='courier',point=8))
         wscroll=.true.
       end if
       call write_scroll_line(icwindow,text=trim(line))

       return
    End Subroutine Write_Scroll_Text

 End Module IO_Messages
