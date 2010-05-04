Program Optimizing_structure
   !---- Use Modules ----!
   use crystallographic_symmetry,only: space_group_type, Write_SpaceGroup
   use string_utilities,         only: u_case
   use Atom_Module,              only: Atom_List_Type, Write_Atom_List
   use crystal_types,            only: Crystal_Cell_Type, Write_Crystal_Cell
   use IO_Formats,               only: Readn_set_Xtal_Structure,err_mess_form,err_form,file_list_type
   use Refinement_Codes,         only: NP_Max,NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,V_Shift,&
                                       V_list,Allocate_VParam,Init_RefCodes, Read_RefCodes_File, &
                                       Write_Info_RefCodes, Err_RefCodes, Err_Mess_RefCodes,       &
                                       Read_RefCodes_File, Read_RefCodes_File, Allocate_restparam, &
                                        NP_REst_dis,NP_REst_ang,NP_REst_tor, Write_restraints_ObsCalc

   use Configuration_Calculations,only : calc_BVS
   use cost_functions,           only: Cell,A,SpG,Ac, dmax, max_coor, General_Cost_function, &
                                       Icost, Wcost, P_cost, Readn_Set_CostFunctPars, err_cost,&
                                       err_mess_cost, Write_CostFunctPars,Write_FinalCost
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

     Call Readn_Set_CostFunctPars(fich_cfl)
     if(Err_cost) then
       write(unit=*,fmt="(a)") trim(Err_Mess_cost)
       stop
     end if
     Call Write_CostFunctPars(lun)

   !  if(icost(2) == 1 .or. icost(3) == 1 .or. icost(4) == 1) then
       call allocate_restparam(fich_cfl)
   !  end if

     call Init_RefCodes(A)
     call Read_RefCodes_File(fich_cfl,1,fich_cfl%nlines,A,Spg)
     if(Err_RefCodes) then
       write(unit=*,fmt="(a)") trim(Err_Mess_RefCodes)
     end if

     call Write_Info_RefCodes(A,SpG, lun)

    !Set up the simulated annealing conditions
     call Set_SimAnn_Cond(fich_cfl,c)
      if(err_SAN) then
         write(unit=*,fmt="(a)") " => Error setting Simulated Annealing conditions"
         write(unit=*,fmt="(a)") trim(err_mess_SAN)
        stop
      end if
     call Write_SimAnn_Cond(lun,c)


     call cpu_time(start)
     if(c%num_conf == 1 ) then    !Conventional single configuration  algorithm

        !Set up the Simulated Annealing State Vector
         call Set_SimAnn_StateV(NP_Refi,V_BCon,V_Bounds,V_Name,V_Vec,vs)
         if(err_SAN) then
            write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
            write(unit=*,fmt="(a)") trim(err_mess_SAN)
           stop
         end if
         call Write_SimAnn_StateV(lun,vs,"INITIAL STATE")
        !call Simanneal_Gen(Cost_function_FobsFcal,c,vs,lun)
        !call Simanneal_Gen(Cost_function_MFobsFcal,c,vs,lun)
         call Simanneal_Gen(General_Cost_function,c,vs,lun)
         call Write_SimAnn_StateV(lun,vs,"FINAL STATE")

     else     !Multiconfigurational algorithm

         call Set_SimAnn_MStateV(NP_Refi,c%Num_Conf,V_BCon,V_Bounds,V_Name,V_Vec,mvs)
         if(err_SAN) then
            write(unit=*,fmt="(a)") " => Error setting Simulated Annealing State Vector"
            write(unit=*,fmt="(a)") trim(err_mess_SAN)
           stop
         end if
         call Write_SimAnn_MStateV(lun,mvs,"INITIAL STATE")
         !call SimAnneal_MultiConf(Cost_function_MulFobsFcal,c,mvs,lun)
         !call SimAnneal_MultiConf(Cost_function_Restraints,c,mvs,lun)
         call SimAnneal_MultiConf(General_Cost_function,c,mvs,lun)
         call Write_SimAnn_MStateV(lun,mvs,"FINAL STATE")
     end if

     call cpu_time(fin)
     write(unit=*,fmt="(a,f10.2,a)")  "   CPU-Time: ", fin-start," seconds"

     call Write_Atom_List(A,lun=lun)
     if(icost(2) == 1 .or. icost(3) == 1 .or. icost(4) == 1) then
        call Write_restraints_ObsCalc(A,lun)
     end if
     if(icost(5) == 1 .or. icost(6) == 1) then
       call Calc_BVS(Ac,lun)
     end if

     call Write_FinalCost(lun)
     call Write_FinalCost(6)

     write(unit=*,fmt="(a)") " Normal End of: PROGRAM FOR OPTIMIZING X-TAL STRUCTURES "
     write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".out"
   end if

   close(unit=lun)
   stop
End Program Optimizing_structure

