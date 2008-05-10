Program Optimizing_structure
   !---- Use Modules ----!
   use crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
   use string_utilities,         only: u_case
   use Atom_Module,              only: Atom_List_Type, Write_Atom_List
   use crystal_types,            only: Crystal_Cell_Type, Write_Crystal_Cell
   use Reflections_Utilities,    only: Reflection_List_Type
   use IO_Formats,               only: Readn_set_Xtal_Structure,err_mess_form,err_form,file_list_type
   use Structure_Factor_Module,  only: Structure_Factors, Write_Structure_Factors, &
                                       Init_Structure_Factors,Calc_StrFactor, err_sfac,err_mess_sfac
   use observed_reflections,     only: Read_observations,Observation_Type,Observation_List_Type, &
                                       Write_ObsCalc_SFactors,err_observ,err_mess_observ, SumGobs
   use Refinement_Codes,         only: NP_Max,NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,V_Shift,&
                                       V_list,Allocate_VParam,Init_RefCodes, Read_RefCodes_File, &
                                       Write_Info_RefCodes, Err_RefCodes, Err_Mess_RefCodes, allocate_restparam
   use cost_functions,           only: Cell,A,SpG,hkl,Oh, Cost_function_FobsFcal,Cost_function_MFobsFcal, &
                                       Cost_function_MulFobsFcal
   use Optimization_SAN,         only: np_SAN, SimAnn_Conditions_type, state_Vector_Type, Multistate_Vector_Type, &
                                       err_mess_SAN,err_SAN, Simanneal_Gen,Set_SimAnn_Cond, &
                                       Set_SimAnn_StateV,Write_SimAnn_Cond, Write_SimAnn_StateV, &
                                       Write_SimAnn_MStateV, SimAnneal_MultiConf,Set_SimAnn_MStateV

   !---- Local Variables ----!
   implicit none

   type (file_list_type)              :: fich_cfl
   type (SimAnn_Conditions_type),save :: c
   type (state_Vector_Type)           :: vs
   type (Multistate_Vector_Type)      :: mvs

   character(len=256)                 :: filcod     !Name of the input file
   character(len=256)                 :: filhkl     !Name of the hkl-file
   character(len=256)                 :: line       !Text line
   integer                            :: MaxNumRef, Num, lun=1, ier,i, i_hkl=2, n
   real                               :: start,fin, mindspc, maxsintl
   integer                            :: narg,iargc
   Logical                            :: esta, arggiven=.false.,sthlgiven=.false.

    !---- Arguments on the command line ----!
    narg=iargc()

    if(narg > 0) then
            call getarg(1,filcod)
            arggiven=.true.
    end if

    write (unit=*,fmt="(/,/,6(a,/))")                                                  &
            "            ------ PROGRAM FOR OPTIMIZING X-TAL STRUCTURES ------"            , &
            "                    ---- Version 0.1 April-2005----"                          , &
            "    ***********************************************************************"  , &
            "    * Optimizes a X-tal structure reading F^2 and a *.CFL or a *.CIF file *"  , &
            "    ***********************************************************************"  , &
            "                        (JRC- April 2005 )"
   write (unit=*,fmt=*) " "

   if(.not. arggiven) then
     write(unit=*,fmt="(a)") " => Code of the file xx.cfl (give xx): "
     read(unit=*,fmt="(a)") filcod
     if(len_trim(filcod) == 0) stop
   end if

   open(unit=lun,file=trim(filcod)//".out", status="replace",action="write")
    write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
         "            ------ PROGRAM FOR OPTIMIZING X-TAL STRUCTURES ------"            , &
         "                    ---- Version 0.1 April-2005----"                          , &
         "    ***********************************************************************"  , &
         "    * Optimizes a X-tal structure reading F^2 and a *.CFL or a *.CIF file *"  , &
         "    ***********************************************************************"  , &
         "                        (JRC- April 2005 )"

   inquire(file=trim(filcod)//".cfl",exist=esta)
   if( .not. esta) then
     write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cif (or .cfl) does'nt exist!"
     stop
   end if
   call Readn_set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

   If(err_form) then
     write(unit=*,fmt="(a)") trim(err_mess_form)
   else
     call Write_Crystal_Cell(Cell,lun)
     call Write_SpaceGroup(SpG,lun)
     call Write_Atom_List(A,lun=lun)

     np_max= A%natoms*11  ! x,y,z,biso,occ,betas
     call allocate_vparam(np_max)
     call allocate_restparam(fich_cfl)
     call Init_RefCodes(A)
     call Read_RefCodes_File(fich_cfl,1,fich_cfl%nlines,A,Spg)
     if(Err_RefCodes) then
       write(unit=*,fmt="(a)") trim(Err_Mess_RefCodes)
     end if
     call Write_Info_RefCodes(A,Spg,lun)

    ! Reading observed structure factors squared and construct hkl and
    ! First detect the name of the hkl file
      esta=.false.
      do i=1,fich_cfl%nlines
        line=u_case(fich_cfl%line(i))
        if(line(1:7) == "HKL-OBS") then
          filhkl= adjustl(fich_cfl%line(i)(9:))
          inquire(file=trim(filhkl),exist=esta)
          exit
        end if
      end do

      if(.not. esta) then
        write(unit=*,fmt="(a)") " => No hkl-file (or wrong name!) has been provided! "
        stop
      end if

      mindspc=0.0       !Read minimun d-spacing to be considered in the list
      do i=1,fich_cfl%nlines
        line=u_case(fich_cfl%line(i))
        if(line(1:12) == "MIN-DSPACING") then
          read(unit=line(14:),fmt=*,iostat=ier) mindspc
          if(ier /= 0) mindspc=0.0
          exit
        end if
      end do

     !call read_observations(trim(filhkl),Cell,SpG,.true.,Oh,hkl)
      call read_observations(trim(filhkl),Cell,SpG,.true.,hkl)
      if(err_observ) then
         write(unit=*,fmt="(a)") " => Error reading observations"
         write(unit=*,fmt="(a)") trim(err_mess_observ)
        stop
      end if

    ! Change the number of reflections to be considered according to the
    ! resolution (mindspc) read above. Change also SumGobs
     if( mindspc > 0.001) then
        n=0
        maxsintl=0.5/mindspc
        maxsintl=maxsintl*1.01
        do i=1,hkl%Nref
          n=n+1
          if(hkl%Ref(i)%s > maxsintl) then
            hkl%Nref=n
            exit
          end if
        end do
        SumGobs=0.0
        do i=1,n
          SumGobs=SumGobs+abs(hkl%Ref(i)%Fo)
        end do
        write(unit= * ,fmt="(a,i5,/,a,f8.4)") " => Number of reflections to consider: ",n,&
                                              "    Resolution: ",mindspc
        write(unit=lun,fmt="(a,i5,/,a,f8.4)") " => Number of reflections to consider: ",n,&
                                              "    Resolution: ",mindspc
     end if

    !Set up Structure factor tables for neutron scattering and initial values
     call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
     call Structure_Factors(A,SpG,hkl,mode="NUC")
      if(err_sfac) then
         write(unit=*,fmt="(a)") " => Error in Structure factor calculations"
         write(unit=*,fmt="(a)") trim(err_mess_sfac)
        stop
      end if

    !Set up the simulated annealing conditions
     call Set_SimAnn_Cond(fich_cfl,c)
      if(err_SAN) then
         write(unit=*,fmt="(a)") " => Error setting Simulated Annealing conditions"
         write(unit=*,fmt="(a)") trim(err_mess_SAN)
        stop
      end if
     call Write_SimAnn_Cond(lun,c)

    !Set up the Simulated Annealing State Vector
    call Set_SimAnn_StateV(NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,vs)
      if(err_SAN) then
         write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
         write(unit=*,fmt="(a)") trim(err_mess_SAN)
        stop
      end if
     call Write_SimAnn_StateV(lun,vs,"INITIAL STATE")
   !call Set_SimAnn_MStateV(NP_Refi,5,V_BCon,V_Bounds,V_Name,V_Vec,mvs)
   !  if(err_SAN) then
   !     write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
   !     write(unit=*,fmt="(a)") trim(err_mess_SAN)
   !    stop
   !  end if
   ! call Write_SimAnn_MStateV(lun,mvs,"INITIAL STATE")

     call cpu_time(start)
      call Simanneal_Gen(Cost_function_FobsFcal,c,vs,lun)
    ! call Simanneal_Gen(Cost_function_MFobsFcal,c,vs,lun)
    ! call SimAnneal_MultiConf(Cost_function_MulFobsFcal,c,mvs,lun)
     call cpu_time(fin)
     write(unit=*,fmt="(a,f10.2,a)")  "  CPU-Time: ", fin-start," seconds"

     call Write_SimAnn_StateV(lun,vs,"FINAL STATE")
    ! call Write_SimAnn_MStateV(lun,mvs,"FINAL STATE")
     call Write_Atom_List(A,lun=lun)
     call Write_ObsCalc_SFactors(lun,hkl,mode="NUC")

     write(unit=*,fmt="(a)") " Normal End of: PROGRAM FOR OPTIMIZING X-TAL STRUCTURES "
     write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".out"
   end if

   close(unit=lun)
   stop
End Program Optimizing_structure

