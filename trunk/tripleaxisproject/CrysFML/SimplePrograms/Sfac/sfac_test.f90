Program Calc_structure_factors

 use crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
 use Atom_Module,              only: Atom_List_Type, Write_Atom_List
 use crystal_types,            only: Crystal_Cell_Type, Write_Crystal_Cell
 use Reflections_Utilities,    only: Reflection_List_Type, Hkl_Uni, get_maxnumref
 use IO_Formats,               only: Readn_set_Xtal_Structure,err_mess_form,err_form,file_list_type
 use Structure_Factor_Module,  only: Structure_Factors, Write_Structure_Factors, &
                                     Init_Structure_Factors,Calc_StrFactor

 implicit none

 type (file_list_type)       :: fich_cfl
 type (space_group_type)     :: SpG
 type (Atom_list_Type)       :: A
 type (Crystal_Cell_Type)    :: Cell
 type (Reflection_List_Type) :: hkl

 character(len=256)          :: filcod     !Name of the input file
 character(len=15)           :: sinthlamb  !String with stlmax (2nd cmdline argument)
 real                        :: stlmax     !Maximum Sin(Theta)/Lambda
 real                        :: sn,sf2
 integer                     :: MaxNumRef, Num, lun=1, ier,i, ln, lr
 complex                     :: fc

 integer                     :: narg,iargc
 Logical                     :: esta, arggiven=.false.,sthlgiven=.false.
      !---- Arguments on the command line ----!
      lr=0
      ln=0
      narg=iargc()

      if(narg > 0) then
              call getarg(1,filcod)
              arggiven=.true.
      end if
      if(narg > 1) then
              call getarg(2,sinthlamb)
              read(unit=sinthlamb,fmt=*,iostat=ier) stlmax
              if(ier == 0) sthlgiven=.true.
      end if

      write(unit=*,fmt="(/,/,6(a,/))")                                                  &
           "            ------ PROGRAM STRUCTURE FACTORS ------"                      , &
           "               ---- Version 0.2 August-2004----"                          , &
           "    *******************************************************************"  , &
           "    * Calculates structure factors reading a *.CFL or a *.CIF file    *"  , &
           "    *******************************************************************"  , &
           "                      (JRC- December 2003 )"
    write(unit=*,fmt=*) " "

     if(.not. arggiven) then
       write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cif(cfl) (give xx): "
       read(unit=*,fmt="(a)") filcod
       if(len_trim(filcod) == 0) stop
     end if
     if(.not. sthlgiven) then
       write(unit=*,fmt="(a)",advance="no") " => Maximum sinTheta/Lambda: "
       read(unit=*,fmt=*) stlmax
     end if

     open(unit=lun,file=trim(filcod)//".sfa", status="replace",action="write")
      write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
           "            ------ PROGRAM STRUCTURE FACTORS ------"                      , &
           "               ---- Version 0.2 August-2004----"                          , &
           "    *******************************************************************"  , &
           "    * Calculates structure factors reading a *.CFL or a *.CIF file    *"  , &
           "    *******************************************************************"  , &
           "                      (JRC- December 2003 )"

     inquire(file=trim(filcod)//".cif",exist=esta)
     if(esta) then
       call Readn_set_Xtal_Structure(trim(filcod)//".cif",Cell,SpG,A,Mode="CIF")
     else
       inquire(file=trim(filcod)//".cfl",exist=esta)
       if( .not. esta) then
         write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) does'nt exist!"
         stop
       end if
       call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)
     end if

     If(err_form) then

       write(unit=*,fmt="(a)") trim(err_mess_form)

     else

       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun)
       call Write_Atom_List(A,lun=lun)

       MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=SpG%Multip)
       call Hkl_Uni(Cell,Spg,.true.,0.0,stlmax,"s",MaxNumRef,hkl)

     !Calculation for neutron scattering
       call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
       call Structure_Factors(A,SpG,hkl,mode="NUC")
       call Write_Structure_Factors(lun,hkl,mode="NUC")

     !Test of the new structure factor subroutine
       write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fc|       Phase        F-Real      F-Imag      Num"
       do i=1, hkl%nref
         sn=hkl%ref(i)%s * hkl%ref(i)%s
         call Calc_StrFactor("P","N",i,sn,A,Spg,sf2,fc=fc)
         write(unit=lun,fmt="(3i4,i5,5f12.5,i8,f12.5)") hkl%ref(i)%h, hkl%ref(i)%mult, &
                                 hkl%ref(i)%S, hkl%ref(i)%Fc, hkl%ref(i)%Phase,   &
                                 real(fc), aimag(fc), i, sqrt(sf2)
       end do

     !Calculation for X-rays
       call Init_Structure_Factors(hkl,A,Spg,lun=lun)
       call Structure_Factors(A,SpG,hkl)
       call Write_Structure_Factors(lun,hkl)

     !Test of the new structure factor subroutine
       write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fc|       Phase        F-Real      F-Imag      Num"
       do i=1, hkl%nref
         sn=hkl%ref(i)%s * hkl%ref(i)%s
         call Calc_StrFactor("P","X",i,sn,A,Spg,sf2,fc=fc)
         write(unit=lun,fmt="(3i4,i5,5f12.5,i8,f12.5)") hkl%ref(i)%h, hkl%ref(i)%mult, &
                                 hkl%ref(i)%S, hkl%ref(i)%Fc, hkl%ref(i)%Phase,   &
                                 real(fc), aimag(fc), i, sqrt(sf2)
       end do

       write(unit=*,fmt="(a)") " Normal End of: PROGRAM STRUCTURE FACTORS "
       write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".sfa"
     end if

     close(unit=lun)
     stop

End Program Calc_structure_factors

