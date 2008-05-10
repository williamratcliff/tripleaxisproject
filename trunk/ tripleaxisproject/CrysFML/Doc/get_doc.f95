!!----
!!---- Copyright(C) 1999-2003,              Version: 2.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- PROGRAM : GET_DOC
!!----           Program to get a Document file from fortran 90 sources
!!----           and generate a HTML Document.
!!----
!!
 Module HTML_File
    !---- Variables ----!
    implicit none

    logical                                     :: private_version
    logical                                     :: modify_html
    logical                                     :: new_html

    character(len=5),  parameter                :: version=" 2.00"
    character(len=15), parameter                :: Date   ="March-2005     "

    character(len=60),dimension(2000)           :: modules_name
    integer                                     :: nmod    ! Number of Modules

    character(len=150), dimension(80000)        :: doc_fortran ! Fortran lines
    integer                                     :: nt_fortran  ! Number of lines of fortran code

    character(len=150), dimension(80000)        :: doc_html    ! Html Lines
    integer                                     :: nt_html     ! Number of lines of Html code

    character(len=80), dimension(2000)          :: var_name    ! Variables Name
    integer                                     :: nvar        ! Number of Variables
    character(len=80), dimension(2000)          :: fun_name    ! Functions Name
    integer                                     :: nfun        ! Number of Functions
    character(len=80), dimension(2000)          :: sub_name    ! Subroutines Name
    integer                                     :: nsub        ! Number of Subroutines
    integer, dimension(2,2000)                  :: doc_bloc    ! Marcador de bloques
    integer                                     :: nt_bloc     ! numeros de bloques

 Contains

    !!----
    !!---- Subroutine: Actualize_HTML ----!!
    !!----
    Subroutine Actualize_Html()
       !---- Local Variables ----!
       character(len=80) :: mname,fname
       character(len=4)  :: car
       integer           :: i,j,n,ipos,imod
       integer           :: n_ini


       !---- Updating CrysFML_Cont File ----!
       do n_ini=1,nt_fortran
          ipos=index(doc_fortran(n_ini),'MODULE:')
          if (ipos ==0) cycle
          mname=adjustl(doc_fortran(n_ini)(ipos+7:))
          if (nmod ==0) then
             nmod=1
             modules_name(nmod)='CrysFML_'//trim(mname)//'.html#mod1'
          else
             n=0
             do i=1,nmod
                ipos=index(modules_name(i),trim(mname))
                if (ipos==0) cycle
                n=i
                exit
             end do
             if (n==0) then
                nmod=nmod+1
                write(car,fmt=*) nmod
                modules_name(nmod)='CrysFML_'//trim(mname)//'.html#mod'//adjustl(car)
             end if
          end if
       end do

       call Create_Html_Cont()

       !---- Updating HTML Files of Modules ----!
       nvar=0
       nfun=0
       nsub=0
       doc_html=' '
       nt_html=0

       n_ini=0
       do
          n_ini=n_ini+1
          if (n_ini > nt_fortran) exit
          ipos=index(doc_fortran(n_ini),'MODULE:')
          if (ipos ==0) cycle
          mname=adjustl(doc_fortran(n_ini)(ipos+7:))
          fname='CrysFML_'//trim(mname)//'.html'
          imod=0
          do j=1,nmod
             ipos=index(modules_name(j),trim(mname))
             if (ipos==0) cycle
             imod=j
             exit
          end do

          !---- Loading Information on doc_html ----!
          write(car,fmt=*) imod
          car=adjustl(car)
          nt_html=1
          doc_html(nt_html)='<h2> <a name="mod'//trim(car)//'" > '//trim(mname)//' </a> </h2>'

          !---- Glossary_List ----!
          call Block_Glossary_List(n_ini)

          !---- Variables List ----!
          call Block_Variables(imod,n_ini)

          !---- Functions List ----!
          call Block_Functions(imod,n_ini)

          !---- Subroutines List ----!
          call Block_Subroutines(imod,n_ini)

          call Create_Html_MFiles(fname)
       end do

       return
    End Subroutine Actualize_Html

    !!----
    !!---- Subroutine: Block_Functions ----!!
    !!----
    Subroutine Block_Functions(imod,n_ini)
       !---- Arguments ----!
       integer, intent(in    ) :: imod
       integer, intent(in out) :: n_ini

       !---- Local variables ----!
       logical           :: lbloc
       character(len=150):: line,line2
       character(len=10) :: car
       integer           :: i,j,jj,k,p_ini,p_fin,p_con,np,nb

       do k=1,nfun
          p_ini=0
          p_fin=0
          do i=n_ini,nt_fortran
             line=adjustl(doc_fortran(i))
             if (len_trim(line) == 0) cycle
             p_ini=i
             do j=i,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (line2(1:7) == "Update:") then
                   p_fin=j
                   exit
                end if
             end do
             if (p_ini ==0 .or. p_fin ==0) then
                write(*,'(a)') ' Function not found!'
                write(*,'(/a)') ' '
                stop
             end if

             !---- Function Name ----!
             select case (k)
                case (:9)
                   car='"fun "  > '
                   write(car(5:5),'(i1)') k

                case (10:99)
                   car='"fun  " > '
                   write(car(5:6),'(i2)') k

                case (100:)
                   car='"fun   "> '
                   write(car(5:7),'(i3)') k
             end select

             !---- Header ----!
             nt_html=nt_html+1
             line=adjustl(doc_fortran(p_ini))
             doc_html(nt_html)="<h4> <a name="//car//trim(line)//" </a> </h4>"
             p_ini=p_ini+1

             !---- Arguments of the Function ----!
             p_con=0
             do j=p_ini,p_fin-1
                if (len_trim(doc_fortran(j)) ==0) then
                   p_con=j
                   exit
                end if
             end do

             if (p_con > p_ini) then
                nt_html=nt_html+1
                doc_html(nt_html)="<p> <pre>"
                do j=p_ini,p_con
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(doc_fortran(j))
                end do
                p_ini=p_con+1
                nt_html=nt_html+1
                doc_html(nt_html)="</pre> </p> "
             end if

             !---- General Information ----!
             lbloc=.false.
             nt_html=nt_html+1
             doc_html(nt_html)="<p> "
             p_for:do j=p_ini,p_fin-1

                if (nt_bloc == 0) then
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                else
                   do jj=1,nt_bloc
                      if (j >= doc_bloc(1,jj) .and. j <= doc_bloc(2,jj)) then

                         if (j == doc_bloc(1,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)="<p> <pre>"
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            lbloc=.true.
                            cycle p_for
                         end if

                         if (j == doc_bloc(2,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            nt_html=nt_html+1
                            doc_html(nt_html)="</pre> </p> "
                            lbloc=.false.
                            cycle p_for

                         end if

                         nt_html=nt_html+1
                         if (lbloc) then
                            doc_html(nt_html)=trim(doc_fortran(j))
                         else
                            doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                         end if
                         cycle p_for
                      end if
                   end do
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                end if

             end do p_for

             nt_html=nt_html+1
             doc_html(nt_html)="</p> "

             write(car,fmt=*) imod
             car=adjustl(car)
             nt_html=nt_html+1
             doc_html(nt_html)=' <a href="#Mod'//trim(car)//'"> [Top Document] </a> '
             nt_html=nt_html+1
             doc_html(nt_html)='<hr size="1" noshade> '
             exit
          end do
          n_ini=p_fin+1
       end do

       return
    End Subroutine Block_Functions

    !!----
    !!---- Subroutine: Block_Glossary_List ----!!
    !!----
    Subroutine Block_Glossary_List(n_ini)
       !---- Arguments ----!
       integer, intent (in out) :: n_ini

       !---- Local variables ----!
       character(len=150) :: line,line2
       character(len=11)  :: car
       integer            :: i, j

       !---- Info of Module ----!
       j=0
       do i=n_ini,nt_fortran
          line=adjustl(doc_fortran(i))
          if (line(1:5) == "INFO:") then
             nt_html=nt_html+1
             doc_html(nt_html)="<dd> "//line(6:)
             do j=i+1,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (len_trim(line2) /= 0) then
                   nt_html=nt_html+1
                   doc_html(nt_html)=line2
                else
                   doc_html(nt_html)=trim(doc_html(nt_html))//" </dd> "
                   exit
                end if
             end do
             exit
          end if
       end do
       if (j > 0) then
          n_ini=j
       else
          n_ini=i
       end if

       !---- Public variables ----!
       nt_html=nt_html+1
       doc_html(nt_html)="<dd> <p> <b> <i> Variables </i> </b> </p> </dd>"

       nt_html=nt_html+1
       doc_html(nt_html)="<ul> "
       do i=n_ini,nt_fortran
          line=adjustl(doc_fortran(i))
          if (line(1:6) == "VARIAB") then
             do j=i+1,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (len_trim(line2) /= 0) then
                   nt_html=nt_html+1
                   nvar=nvar+1
                   var_name(nvar)=trim(line2)
                   select case (nvar)
                      case (:9)
                         car='"#var "  > '
                         write(car(6:6),'(i1)') nvar

                      case (10:99)
                         car='"#var  " > '
                         write(car(6:7),'(i2)') nvar

                      case (100:)
                         car='"#var   "> '
                         write(car(6:8),'(i3)') nvar
                   end select
                   doc_html(nt_html)="  <li> <a href="//car//trim(var_name(nvar))//" </a> </li> "
                else
                   n_ini=j
                   exit
                end if
             end do
             exit
          end if
       end do
       nt_html=nt_html+1
       doc_html(nt_html)="</ul> "

       !---- Public functions ----!
       nt_html=nt_html+1
       doc_html(nt_html)="<dd> <p> <b> <i> Functions </i> </b> </p> </dd>"

       nt_html=nt_html+1
       doc_html(nt_html)="<ul> "
       do i=n_ini,nt_fortran
          line=adjustl(doc_fortran(i))
          if (line(1:10) == "Functions:") then
             do j=i+1,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (len_trim(line2) /= 0) then
                   nt_html=nt_html+1
                   nfun=nfun+1
                   fun_name(nfun)=trim(line2)
                   select case (nfun)
                      case (:9)
                         car='"#fun "  > '
                         write(car(6:6),'(i1)') nfun

                      case (10:99)
                         car='"#fun  " > '
                         write(car(6:7),'(i2)') nfun

                      case (100:)
                         car='"#fun   "> '
                         write(car(6:8),'(i3)') nfun
                   end select
                   doc_html(nt_html)="  <li> <a href="//car//trim(fun_name(nfun))//" </a> </li> "
                else
                   n_ini=j
                   exit
                end if
             end do
             exit
          end if
       end do
       nt_html=nt_html+1
       doc_html(nt_html)="</ul> "

       !---- Public subroutines ----!
       nt_html=nt_html+1
       doc_html(nt_html)="<dd> <p> <b> <i> Subroutines </i> </b> </p> </dd>"

       nt_html=nt_html+1
       doc_html(nt_html)="<ul> "
       do i=n_ini,nt_fortran
          line=adjustl(doc_fortran(i))
          if (line(1:12) == "Subroutines:") then
             do j=i+1,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (len_trim(line2) /= 0) then
                   nt_html=nt_html+1
                   nsub=nsub+1
                   sub_name(nsub)=trim(line2)
                   select case (nsub)
                      case (:9)
                         car='"#sub "  > '
                         write(car(6:6),'(i1)') nsub

                      case (10:99)
                         car='"#sub  " > '
                         write(car(6:7),'(i2)') nsub

                      case (100:)
                         car='"#sub   "> '
                         write(car(6:8),'(i3)') nsub
                   end select
                   doc_html(nt_html)="  <li> <a href="//car//trim(sub_name(nsub))//" </a> </li> "
                else
                   n_ini=j
                   exit
                end if
             end do
             exit
          end if
       end do
       nt_html=nt_html+1
       doc_html(nt_html)="</ul> "

       nt_html=nt_html+1
       doc_html(nt_html)="</dl> "

       nt_html=nt_html+1
       doc_html(nt_html)='<hr size="1" noshade> '

       return
    End Subroutine Block_Glossary_List

    !!----
    !!---- Subroutine: Block_Subroutines ----!!
    !!----
    Subroutine Block_Subroutines(imod,n_ini)
       !---- Arguments ----!
       integer, intent (in    ) :: imod
       integer, intent (in out) :: n_ini

       !---- Local variables ----!
       logical           :: lbloc
       character(len=150):: line,line2
       character(len=10) :: car
       integer           :: i,j,jj,k,p_ini,p_fin,p_con,np

       do k=1,nsub
          p_ini=0
          p_fin=0
          do i=n_ini,nt_fortran
             line=adjustl(doc_fortran(i))
             if (len_trim(line) == 0) cycle
             p_ini=i
             do j=i,nt_fortran
                line2=adjustl(doc_fortran(j))
                if (line2(1:7) == "Update:") then
                   p_fin=j
                   exit
                end if
             end do
             if (p_ini ==0 .or. p_fin ==0) then
                write(*,'(a)') ' Subroutine not found!'
                write(*,'(/a)') ' '
                stop
             end if

             !---- Subroutine Name ----!
             select case (k)
                case (:9)
                   car='"sub "  > '
                   write(car(5:5),'(i1)') k

                case (10:99)
                   car='"sub  " > '
                   write(car(5:6),'(i2)') k

                case (100:)
                   car='"sub   "> '
                   write(car(5:7),'(i3)') k
             end select

             !---- Header ----!
             nt_html=nt_html+1
             line2=adjustl(doc_fortran(p_ini))
             doc_html(nt_html)="<h4> <a name="//car//trim(line2)//" </a> </h4>"
             p_ini=p_ini+1

             !---- Arguments of the subroutine ----!
             p_con=0
             do j=p_ini,p_fin-1
                if (len_trim(doc_fortran(j)) ==0) then
                   p_con=j
                   exit
                end if
             end do

             if (p_con > p_ini) then
                nt_html=nt_html+1
                doc_html(nt_html)="<p> <pre>"
                do j=p_ini,p_con
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(doc_fortran(j))
                end do
                p_ini=p_con+1
                nt_html=nt_html+1
                doc_html(nt_html)="</pre> </p> "
             end if

             !---- General Information ----!
             lbloc=.false.
             nt_html=nt_html+1
             doc_html(nt_html)="<p> "
             p_for:do j=p_ini,p_fin-1

                if (nt_bloc == 0) then
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                else
                   do jj=1,nt_bloc
                      if (j >= doc_bloc(1,jj) .and. j <= doc_bloc(2,jj)) then

                         if (j == doc_bloc(1,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)="<p> <pre>"
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            lbloc=.true.
                            cycle p_for
                         end if

                         if (j == doc_bloc(2,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            nt_html=nt_html+1
                            doc_html(nt_html)="</pre> </p> "
                            lbloc=.false.
                            cycle p_for

                         end if

                         nt_html=nt_html+1
                         if (lbloc) then
                            doc_html(nt_html)=trim(doc_fortran(j))
                         else
                            doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                         end if
                         cycle p_for
                      end if
                   end do
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                end if

             end do p_for

             nt_html=nt_html+1
             doc_html(nt_html)="</p> "

             write(car,fmt=*) imod
             car=adjustl(car)
             nt_html=nt_html+1
             doc_html(nt_html)=' <a href="#Mod'//trim(car)//'"> [Top Document] </a> '
             nt_html=nt_html+1
             doc_html(nt_html)='<hr size="1" noshade> '

             exit
          end do
          n_ini=p_fin+1
       end do

       return
    End Subroutine Block_Subroutines

    !!----
    !!---- Subroutine: Block_Variables ----!!
    !!----
    Subroutine Block_Variables(imod,n_ini)
       !---- Arguments ----!
       integer, intent (in    ) :: imod
       integer, intent (in out) :: n_ini

       !---- Local variables ----!
       logical            :: lbloc
       character(len=150) :: line,line2
       character(len=10)  :: car
       integer            :: i,j,jj,k,p_ini,p_fin,p_con,np,nb

       do k=1,nvar
          p_ini=0
          p_fin=0
          do i=n_ini,nt_fortran
             if (len_trim(doc_fortran(i)) == 0) cycle
             p_ini=i
             do j=i,nt_fortran
                line=adjustl(doc_fortran(j))
                if (line(1:7) == "Update:") then
                   p_fin=j
                   exit
                end if
             end do
             if (p_ini ==0 .or. p_fin ==0) then
                write(unit=*,fmt='(a)') ' Variable not found!'
                write(unit=*,fmt='(/a)') ' '
                stop
             end if

             !---- Variable Name ----!
             select case (k)
                case (:9)
                   car='"var "  > '
                   write(car(5:5),'(i1)') k

                case (10:99)
                   car='"var  " > '
                   write(car(5:6),'(i2)') k

                case (100:)
                   car='"var   "> '
                   write(car(5:7),'(i3)') k
             end select

             !---- Header ----!
             nt_html=nt_html+1
             line=adjustl(doc_fortran(p_ini))
             doc_html(nt_html)="<h4> <a name="//car//trim(line)//" </a> </h4>"
             p_ini=p_ini+1

             !---- Definition of Variables ----!
             p_con=0
             do j=p_ini,p_fin-1
                if (len_trim(doc_fortran(j)) ==0) then
                   p_con=j
                   exit
                end if
             end do

             if (p_con > p_ini) then
                nt_html=nt_html+1
                doc_html(nt_html)="<p> <pre>"
                do j=p_ini,p_con
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(doc_fortran(j))
                end do
                p_ini=p_con+1
                nt_html=nt_html+1
                doc_html(nt_html)="</pre> </p> "
             end if

             !---- General Information ----!
             lbloc=.false.
             nt_html=nt_html+1
             doc_html(nt_html)="<p> "
             p_for:do j=p_ini,p_fin-1

                if (nt_bloc == 0) then
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                else
                   do jj=1,nt_bloc
                      if (j >= doc_bloc(1,jj) .and. j <= doc_bloc(2,jj)) then

                         if (j == doc_bloc(1,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)="<p> <pre>"
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            lbloc=.true.
                            cycle p_for
                         end if

                         if (j == doc_bloc(2,jj)) then
                            nt_html=nt_html+1
                            doc_html(nt_html)=trim(doc_fortran(j))
                            nt_html=nt_html+1
                            doc_html(nt_html)="</pre> </p> "
                            lbloc=.false.
                            cycle p_for

                         end if

                         nt_html=nt_html+1
                         if (lbloc) then
                            doc_html(nt_html)=trim(doc_fortran(j))
                         else
                            doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                         end if
                         cycle p_for
                      end if
                   end do
                   nt_html=nt_html+1
                   doc_html(nt_html)=trim(adjustl(doc_fortran(j)))
                end if

             end do p_for

             nt_html=nt_html+1
             doc_html(nt_html)="</p> "

             write(car,fmt=*) imod
             car=adjustl(car)
             nt_html=nt_html+1
             doc_html(nt_html)=' <a href="#Mod'//trim(car)//'"> [Top Document] </a> '
             nt_html=nt_html+1
             doc_html(nt_html)='<hr size="1" noshade> '
             exit
          end do
          n_ini=p_fin+1
       end do

       return
    End Subroutine Block_Variables

    !!----
    !!---- Subroutine: Create_HTML_Cont
    !!----
    Subroutine Create_Html_Cont()
       !---- Local Variables ----!
       character(len=80)  :: fname_html="CrysFML_Cont.html"
       character(len=80)  :: mname
       integer, parameter :: iunit=10
       integer            :: i,ipos, ipos1,ipos2

       !---- Init ----!
       open(unit=iunit,file=fname_html,status='replace')

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "<html>"

       !---- HEAD Part ----!
       write(unit=iunit,fmt='(a)') "  <head>"
       write(unit=iunit,fmt='(a)') "    <title> CrysFML </title>"
       write(unit=iunit,fmt='(a)') '    <base target="principal">'
       write(unit=iunit,fmt='(a)') "  </head>"

       !---- BODY Part ----!
       write(unit=iunit,fmt='(a)') '  <body bgcolor="#FFFFCC">'
       write(unit=iunit,fmt='(a)') '    <hr size="1" noshade>'
       write(unit=iunit,fmt='(a)') '    <p> <h2> Modules Information </h2> </p>'
       write(unit=iunit,fmt='(a)') '    <ul>'
       do i=1,nmod
          ipos=index(modules_name(i),'CrysFML_')
          ipos1=ipos+8
          ipos2=index(modules_name(i),'.html')
          mname=modules_name(i)(ipos1:ipos2-1)
          write(unit=iunit,fmt='(a)') '     <li> <a target="principal" href="'//trim(modules_name(i))// &
                                      '" > '//trim(mname)//' </a> </li>'
       end do
       write(unit=iunit,fmt='(a)') '    </ul>'
       write(unit=iunit,fmt='(a)') '    <hr size="1" noshade>'
       write(unit=iunit,fmt='(a)') '    <p><b>Authors: </b> <br>'
       write(unit=iunit,fmt='(a)') '    <i>Juan Rodriguez-Carvajal</i> <a href="mailto:juan@llb.saclay.cea.fr"> juan@llb.saclay.cea.fr </a><br>'
       write(unit=iunit,fmt='(a)') '    <i>Javier Gonzalez-Platas</i> <a href="mailto:jplatas@ull.es"> jplatas@ull.es </a> <br>'
       write(unit=iunit,fmt='(a)') "  </body>"

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "</html>"

       close(unit=iunit)

       return
    End Subroutine Create_Html_Cont

    !!----
    !!---- Subroutine: Create_HTML_Main ----!!
    !!----
    Subroutine Create_Html_Main()

       !---- Local Variables ----!
       character(len=80)  :: fname_html="CrysFML_Doc.html"
       integer, parameter :: iunit=10

       !---- Init ----!
       open(unit=iunit,file=fname_html,status='replace')

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "<html>"

       !---- HEAD Part ----!
       write(unit=iunit,fmt='(a)') "  <head>"
       write(unit=iunit,fmt='(a)') "    <title>  CrysFML</title>"
       write(unit=iunit,fmt='(a)') "  </head>"

       !---- FRAME Part ----!
       write(unit=iunit,fmt='(a)') '  <frameset rows="116,*">'
       write(unit=iunit,fmt='(a)') '    <frame name="titular" scrolling="no" noresize target="contenido" src="CrysFML_Title.html">'
       write(unit=iunit,fmt='(a)') '    <frameset cols="270,*">'
       write(unit=iunit,fmt='(a)') '      <frame name="contenido" target="principal" src="CrysFML_Cont.html" scrolling="auto">'
       write(unit=iunit,fmt='(a)') '      <frame name="principal" src="CrysFML_Presen.html" scrolling="auto">'
       write(unit=iunit,fmt='(a)') '    </frameset>'
       write(unit=iunit,fmt='(a)') '  </frameset>'

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "</html>"

       close(unit=iunit)

       return
    End Subroutine Create_Html_Main

    !!----
    !!---- Subroutine: Create_HTML_MFiles ----!!
    !!----
    Subroutine Create_Html_MFiles(fname)
       !---- Arguments ----!
       character(len=*), intent(in) :: fname

       !---- Local Variables ----!
       integer, parameter :: iunit=10
       integer            :: i,ipos

       !---- Init ----!
       open(unit=iunit,file=fname,status='replace')

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "<html>"

       !---- HEAD Part ----!
       write(unit=iunit,fmt='(a)') "  <head>"
       write(unit=iunit,fmt='(a)') "    <title> CrysFML </title>"
       write(unit=iunit,fmt='(a)') "  </head>"

       !---- BODY Part ----!
       write(unit=iunit,fmt='(a)') '  <body bgcolor="#FFFFCC">'
       do i=1,nt_html
          write(unit=iunit,fmt='(a)') '    '//doc_html(i)
       end do
       write(unit=iunit,fmt='(a)') "  </body>"

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "</html>"

       close(unit=iunit)

       return
    End Subroutine Create_Html_MFiles

    !!----
    !!---- Subroutine: Create_HTML_Presentation ----!!
    !!----
    Subroutine Create_Html_Presentation()

       !---- Local Variables ----!
       character(len=80)  :: fname_html="CrysFML_Presen.html"
       integer, parameter :: iunit=10

       !---- Init ----!
       open(unit=iunit,file=fname_html,status='replace')

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "<html>"

       !---- HEAD Part ----!
       write(unit=iunit,fmt='(a)') "  <head>"
       write(unit=iunit,fmt='(a)') "    <title> CrysFML</title>"
       write(unit=iunit,fmt='(a)') "  </head>"

       !---- BODY Part ----!
       write(unit=iunit,fmt='(a)') '  <body bgcolor="#FFFFCC">'
       write(unit=iunit,fmt='(a)') "    <p></p>"
       write(unit=iunit,fmt='(a)') "    <p><h1><b><center>Crystallographic Fortran Modules Library (CrysFML)</center></b></h1></p>"
       write(unit=iunit,fmt='(a)') "    <p><h2><center> A simple toolbox for crystallographic computing programs</center></h2></p><br>"
       write(unit=iunit,fmt='(a)') "    <p></p>"
       write(unit=iunit,fmt='(a)') "    <p><h4> This is a set of modules written in a subset of Fortran95 (F-language)"
       write(unit=iunit,fmt='(a)') "            to facilitate the design of crystallographic computing programs</h4></p>"
       write(unit=iunit,fmt='(a)') "  </body>"

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "</html>"

       close(unit=iunit)

       return
    End Subroutine Create_Html_Presentation

    !!----
    !!---- Subroutine: Create_HTML_Title ----!!
    !!----
    Subroutine Create_Html_Title()

       !---- Local Variables ----!
       character(len=80)  :: fname_html="CrysFML_Title.html"
       integer, parameter :: iunit=10

       !---- Init ----!
       open(unit=iunit,file=fname_html,status='replace')

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "<html>"

       !---- HEAD Part ----!
       write(unit=iunit,fmt='(a)') "  <head>"
       write(unit=iunit,fmt='(a)') "    <title>  Crystallographic Fortran Modules Library (CrysFML) </title>"
       write(unit=iunit,fmt='(a)') "  </head>"

       !---- BODY Part ----!
       write(unit=iunit,fmt='(a)') '  <body bgcolor="#FFFFCC">'
       write(unit=iunit,fmt='(a)') "    <p><h1><b>Crystallographic Fortran Modules Library (CrysFML)</b></h1></p>"
       write(unit=iunit,fmt='(a)') "    <p><b>JRC - JGP / "//trim(date)//" / Revision:"//trim(version)//" </b></p>"
       write(unit=iunit,fmt='(a)') "  </body>"

       !---- HTML Tag ----!
       write(unit=iunit,fmt='(a)') "</html>"

       close(unit=iunit)

       return
    End Subroutine Create_Html_Title

    !!----
    !!---- Subroutine: Read_HTML_Cont
    !!----
    Subroutine Read_Html_Cont()
       !---- Local variables ----!
       character(len=150) :: line
       integer, parameter :: iunit=10
       integer            :: i,ierror
       integer            :: ipos, ipos1,ipos2

       !---- Read Html File ----!
       open(unit=iunit,file='CrysFML_Cont.html',status='old')
       nt_html=0
       doc_html=' '
       do
          read(unit=iunit,fmt='(a)',iostat=ierror) line
          if (ierror /=0) exit

          nt_html=nt_html+1
          if (nt_html <= 50000) then
             doc_html(nt_html)=line
          else
             write(unit=*,fmt='(a)') ' Exceed lines from HTML Document! '
             stop
          end if
       end do
       close(unit=iunit)

       !---- Load Module List ----!
       nmod=0
       modules_name=' '
       do i=1,nt_html
          line=adjustl(doc_html(i))
          if (line(1:4)=='<li>') then
             nmod=nmod+1
             ipos=index(line,'href=')
             ipos1=ipos+6
             ipos2=index(line,'" >')
             modules_name(nmod)=adjustl(line(ipos1:ipos2-1))
          end if
       end do

       return
    End Subroutine Read_Html_Cont

    !!----
    !!---- Subroutine: Write_HTML ----!!
    !!----
    Subroutine Write_Html()
       !---- Local variables ----!
       integer            :: i

       if (new_html) then
          call Create_Html_Main()
          call Create_Html_Title()
          call Create_Html_Presentation()
       end if

       if (modify_html) then
          call Read_Html_Cont()
       end if

       call Actualize_Html()

       return
    End Subroutine Write_Html

 End Module HTML_File

 !!--------------------------------------------------------------!!
 !!--------------------------------------------------------------!!
 !!----                                                      ----!!
 !!----                 Program Get_Doc                      ----!!
 !!----                                                      ----!!
 !!----   Aim: Produce HTML files from special comment       ----!!
 !!----        lines in Fortran 95/2003 source code files    ----!!
 !!----                                                      ----!!
 !!----                                                      ----!!
 !!----   Use of the program:                                ----!!
 !!----                                                      ----!!
 !!----   My_Prompt> get_doc arg1 arg2   <cr>                ----!!
 !!----                                                      ----!!
 !!----   where                                              ----!!
 !!----   arg1: is "private" or whatever other string        ----!!
 !!----         In such a case all the information lines     ----!!
 !!----         making reference to private or overloaded    ----!!
 !!----         procedures are not output into the final     ----!!
 !!----         HTML file                                    ----!!
 !!----                                                      ----!!
 !!----   arg2: is the name of a file containing a list      ----!!
 !!----         of Fortran files from which we want to get   ----!!
 !!----         the documentation lines in HTML.             ----!!
 !!----                                                      ----!!
 !!----   <cr>  denotes press the enter key                  ----!!
 !!----                                                      ----!!
 !!----   If one invokes de program just with its name, the  ----!!
 !!----   program prompts for a fortran file in the current  ----!!
 !!----   directory and outputs an HTML documentation file.  ----!!
 !!----                                                      ----!!
 !!--------------------------------------------------------------!!
 !!--------------------------------------------------------------!!
 Program Get_Doc
    !---- Use File ----!
    use HTML_File

    !---- Local Variables ----!
    implicit none

    logical            :: info
    logical            :: add_lines
    character(len=10)  :: dire
    character(len=80)  :: filef95, list_file
    character(len=150) :: line
    character(len=1)   :: car

    integer, parameter :: iunit=10
    integer            :: narg, iargc
    integer            :: i,ipos,ierror,num_lines,num_info, tot_lines

    !---- Init Variables ----!
    private_version=.false.
    modify_html    =.false.
    new_html       =.true.

    !---- Presentation ----!
    write(unit=*,fmt='(/a)') ' *******************************************************************'
    write(unit=*,fmt='(a)')  ' **** GET_DOC Program (2.0)              Copyright(C) 1998-2005 ****'
    write(unit=*,fmt='(a)')  ' **** Juan Rodriguez-Carvajal      &     Javier Gonzalez-Platas ****'
    write(unit=*,fmt='(a)')  ' *******************************************************************'
    write(unit=*,fmt='(a)')  ' '

    !---- Arguments ----!
    narg=iargc()
    if (narg > 0) call getarg(1,dire)
    if (dire(1:7) == 'private') private_version=.true.

    !---- HTML Document ----!
    inquire(file="CrysFML_Doc.html",exist=info)
    if (info) then
       !---- Updating ----!
       modify_html=.true.
       new_html   =.false.
    else
       !---- New ----!
       new_html   =.true.
       modify_html=.false.

       nmod=0
       nvar=0
       nfun=0
       nsub=0
       modules_name=' '
    end if

    !---- Reading Fortran Files ----!

    if (narg >= 1) then   !a list of fortran files is provided in "list_file"

       call getarg(2,list_file)
       inquire(file=list_file,exist=info)
       if(info) then
         open(unit=1,file=list_file, status="old",action="read",position="rewind")
       end if

       num_lines=0
       num_info=0
       tot_lines=0

       do
          read(unit=1,fmt="(a)",iostat=ierror) filef95
          if(ierror /= 0 ) exit

          if(len_trim(filef95) == 0) cycle

          ipos=index(filef95,'.',back=.true.)
          if (ipos == 0) then
             filef95=trim(filef95)//'.f95'
          end if

          !---- Exist ? ----!
          inquire(file=filef95,exist=info)
          if (info) then
             open(unit=iunit,file=filef95,status='old',action="read",position="rewind")
          else
             write(unit=*,fmt='(a)') ' File '//trim(filef95)//' not found!'
             write(unit=*,fmt='(/a)') ' '
             cycle    !continue with the next line
          end if

          write(unit=*,fmt='(a)') ' => Treating File: '//trim(filef95)

          !---- Read Fortran File ----!
          nt_fortran=0
          nt_bloc=0
          add_lines=.false.
          do
             read(unit=iunit,fmt='(a)',iostat=ierror) line
             if (ierror /=0) exit
             tot_lines=tot_lines+1
             line=adjustl(line)
             if(len_trim(line) == 0 .or. line(1:2) == "! ") cycle
             num_lines=num_lines+1
             if (line(1:4) /= '!!--') cycle

             select case (line(5:6))

                case ('--')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)

                case ('++')
                   if (private_version) then
                      nt_fortran=nt_fortran+1
                      doc_fortran(nt_fortran)=line(7:)
                   end if

                case ('<<')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)
                   nt_bloc=nt_bloc+1
                   doc_bloc(1,nt_bloc)=nt_fortran
                   add_lines=.true.

                case ('>>')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)
                   doc_bloc(2,nt_bloc)=nt_fortran
                   add_lines=.false.

                case default

                   if (add_lines) then
                      nt_fortran=nt_fortran+1
                      doc_fortran(nt_fortran)=line(7:)
                   end if

             end select

          end do
          close(unit=iunit)
          num_info=num_info + nt_fortran

          call Write_Html()

       end do
       close(unit=1)



    else

       num_lines=0
       num_info=0
       tot_lines=0

       do
          !---- Name of File ----!
          write(unit=*,fmt='(a)') " Name of the Fortran 90 input file: "
          read(unit=*,fmt='(a)') filef95
          if (len_trim(filef95) == 0) exit
          ipos=index(filef95,'.',back=.true.)
          if (ipos == 0) then
             filef95=trim(filef95)//'.f95'
          end if

          !---- Exist ? ----!
          inquire(file=filef95,exist=info)
          if (info) then
             open(unit=iunit,file=filef95,status='old')
          else
             write(unit=*,fmt='(a)') ' File '//trim(filef95)//' not found!'
             write(unit=*,fmt='(/a)') ' '
             stop
          end if

          !---- Read Fortran File ----!

          add_lines=.false.
          nt_fortran=0
          nt_bloc=0

          do
             read(unit=iunit,fmt='(a)',iostat=ierror) line
             if (ierror /=0) exit
             tot_lines=tot_lines+1
             line=adjustl(line)
             if(len_trim(line) == 0 .or. line(1:2) == "! ") cycle
             num_lines=num_lines+1

             if (line(1:4) /= '!!--') cycle
             select case (line(5:6))

                case ('--')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)

                case ('++')
                   if (private_version) then
                      nt_fortran=nt_fortran+1
                      doc_fortran(nt_fortran)=line(7:)
                   end if

                case ('<<')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)
                   nt_bloc=nt_bloc+1
                   doc_bloc(1,nt_bloc)=nt_fortran
                   add_lines=.true.

                case ('>>')
                   nt_fortran=nt_fortran+1
                   doc_fortran(nt_fortran)=line(7:)
                   doc_bloc(2,nt_bloc)=nt_fortran
                   add_lines=.false.

                case default
                   if (add_lines) then
                      nt_fortran=nt_fortran+1
                      doc_fortran(nt_fortran)=line(7:)
                   end if

             end select

          end do
          close(unit=iunit)
          num_info=num_info + nt_fortran
             !---- Write HTML Document ----!
          call Write_Html()

       end do
    end if

    if(num_info == 0) then
      write(unit=*,fmt='(a)') " => ERROR: no information lines, in fortran files, have been found!"
    else
      write(unit=*,fmt='(/,a,i6,a)') " => All Fortran files have ",num_info," information lines"
      num_lines=num_lines-num_info
      write(unit=*,fmt='(a,i7)')     " => The total number of lines in the code   is ",num_lines
      write(unit=*,fmt='(a,i7)')     " => The total number of lines read in files is ",tot_lines
    end if


 End Program Get_Doc
