  Module cost_functions
      Use crystallographic_symmetry  ,only: space_group_type, ApplySO, Read_SymTrans_Code
      Use Geom_Calculations,          only: distance,angle_uv,Angle_dihedral, Set_tdist_coordination,&
                                            allocate_coordination_type
      Use Configuration_Calculations, only: Cost_BVS, Atoms_Conf_list_Type, species_on_list, &
                                            Allocate_Atoms_Conf_list, err_conf, err_mess_conf, &
                                            Cost_BVS_CoulombRep, set_Table_d0_B

      use IO_Formats,                 only: file_list_type
      use string_utilities,           only: l_case
      Use Atom_Module,                only: Atom_List_Type
      Use crystal_types,              only: Crystal_Cell_Type
      Use Reflections_Utilities,      only: Reflection_List_Type
      Use Structure_Factor_Module,    only: Structure_Factors,Init_Structure_Factors, &
                                            Calc_StrFactor, Modify_SF
      Use observed_reflections,       only: Observation_Type,Observation_List_Type, SumGobs,ScaleFact
      Use Molecular_Crystals,         only: Molecular_Crystal_Type

      Use Refinement_Codes,           only: NP_Max,NP_Refi, v_Vec,v_Shift,v_Bounds,v_BCon,v_Name, v_List, &
                                            VState_to_AtomsPar, &
                                            Distance_Restraint_Type, Dis_Rest, NP_Rest_Dis, &
                                            Angle_Restraint_Type, Ang_Rest, NP_Rest_Ang,    &
                                            Torsion_Restraint_Type, Tor_Rest, NP_Rest_Tor
      implicit none
      private

      public::  Cost_function_FobsFcal, Cost_function_MFobsFcal, Cost_function_MulFobsFcal,Cost_function_Restraints, &
                General_cost_function, Readn_Set_CostFunctPars, Write_CostFunctPars, Write_FinalCost

      logical,                      public :: err_cost=.false.
      character(len=132),           public :: err_mess_cost=" "
      type (space_group_type),      public :: SpG
      type (Atom_list_Type),        public :: A
      type (Atoms_Conf_list_Type),  public :: Ac
      type (Crystal_Cell_Type),     public :: Cell
      type (Reflection_List_Type),  public :: hkl
      type (Observation_List_Type), public :: Oh

      Integer, parameter,           public :: N_costf=7
      Integer,dimension(N_costf),   public :: Icost
      real,   dimension(N_costf),   public :: Wcost
      real,   dimension(N_costf),   public :: P_cost !Partial cost

      integer,  public :: Max_Coor
      real,     public :: Dmax
      character(len=3) :: diff_mode="NUC"   ! XRA for x-rays

    Contains

    Subroutine Readn_Set_CostFunctPars(file_dat)
       !---- Arguments ----!
       Type(file_list_type),   intent( in)  :: file_dat
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i,ier,j


       Icost=0
       wcost=0.0
       dmax=3.2

       do j=1,file_dat%nlines
          line=adjustl(file_dat%line(j))
          line=l_case(line)
          if (line(1:1) ==" ") cycle
          if (line(1:1) =="!") cycle

          select case (line(1:4))

              case ("opti")    !Optimization

                 i=index(line,"Fobs-Fcal")
                 if(i == 0) then
                   Icost(1)=0; wcost(1)=0.0
                 else
                   Icost(1)=1
                   read(unit=line(10:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(1)=1.0
                   else
                      wcost(1)= w
                   end if
                 end if

                 i=index(line,"dis-restr")
                 if(i == 0) then
                   Icost(2)=0; wcost(2)=0.0
                 else
                   Icost(2)=1
                   read(unit=line(10:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(2)=1.0
                   else
                      wcost(2)= w
                   end if
                 end if

                 i=index(line,"ang-restr")
                 if(i == 0) then
                   Icost(3)=0; wcost(3)=0.0
                 else
                   Icost(3)=1
                   read(unit=line(10:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(3)=1.0
                   else
                      wcost(3)= w
                   end if
                 end if

                 i=index(line,"tor-restr")
                 if(i == 0) then
                   Icost(4)=0; wcost(4)=0.0
                 else
                   Icost(4)=1
                   read(unit=line(10:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(4)=1.0
                   else
                      wcost(4)= w
                   end if
                 end if

                 i=index(line,"bond-valence")
                 if(i == 0) then
                   Icost(5)=0; wcost(5)=0.0
                 else
                   Icost(5)=1
                   read(unit=line(13:),fmt=*,iostat=ier) w
                   if(ier /= 0) then
                      wcost(5)=1.0
                   else
                      wcost(5)= w
                   end if
                 end if

                 i=index(line,"bvs_coulomb")
                 if(i == 0) then
                   Icost(6)=0; wcost(6)=0.0
                 else
                   Icost(6)=1
                   read(unit=line(13:),fmt=*,iostat=ier) w,tol
                   if(ier /= 0) then
                      wcost(5)=1.0 ; wcost(6)=1.0
                   else
                      wcost(5)= w ; wcost(6)= tol
                   end if
                 end if

              case ("dmax")

                 read(unit=line(5:),fmt=*,iostat=ier) w
                 if( ier == 0) dmax=w

              case ("tol ")

                 read(unit=line(5:),fmt=*,iostat=ier) w
                 if( ier == 0) tol=w

              case ("radi")    !Radiation
                 i=index(line,"xra")
                 if(i /= 0) diff_mode="XRA"

          end select

       end do

       !Normalize weight vector
        w=sum(Wcost)
        Wcost=Wcost/w

       !Allocate the necesary types and arrays
        bond: if(Icost(5) == 1 .or. Icost(6) == 1) then !Bond-Valence parameters

          if(A%natoms == 0) then

             err_cost=.true.
             err_mess_cost=" Fatal Error: 0 atoms! => Atom-list must be allocated before calling Readn_Set_CostFunctPars"
             return

          else
             !Construction of the Atom_Conf_List variable "Ac"
             call Allocate_Atoms_Conf_list(A%natoms,Ac)
             Ac%atom=A%atom
             Call Species_on_List(Ac,Spg%Multip,tol)
             if(err_conf) then
               err_cost=.true.
               err_mess_cost=err_mess_conf
               return
             end if

             call set_Table_d0_B(Ac)
             if(err_conf) then
               err_cost=.true.
               err_mess_cost=err_mess_conf
               return
             end if
             !Allocation of the global variable "Coord_Info" in Geom_Calculations
             !to be used in distance calculations and Bond-valence
             call Allocate_Coordination_type(A%natoms,Spg%Multip,dmax,Max_coor)

          end if

        end if  bond

       return
    End Subroutine Readn_Set_CostFunctPars

    Subroutine Write_CostFunctPars(lun)
       !---- Arguments ----!
       integer,   intent( in)    :: lun
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i,ier


       Write(unit=lun,fmt="(/,a)")    "=================================="
       Write(unit=lun,fmt="(a)")      "= Cost functions to be minimized ="
       Write(unit=lun,fmt="(a,/)")    "=================================="

       do i=1,n_costf

          if(icost(i) == 0) cycle
          select case (i)

              case (1)    !Optimization "Fobs-Fcal"

                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C1=Sum|Fobs-Fcal|/Sum|Fobs|, with weight: ",wcost(i)

              case (2)     !Optimization "dis-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C2=Sum{w(dobs-dcal)^2}, w=1/var(d),with weight: ", wcost(i)

              case (3)     !Optimization "Ang-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C3=Sum{w(Ang_obs-Ang_cal)^2}, w=1/var(Ang),with weight: ", wcost(i)

              case (4)     !Optimization "Tor-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C4=Sum{w(Tor_obs-Tor_cal)^2}, w=1/var(Tor),with weight: ", wcost(i)

              case (5)     !Optimization "bond-valence"
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C5=Sum{|q-BVS|/tot_Atoms}, with weight: ", wcost(i)

              case (6)     !Optimization "bvs_coulomb"
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C5=Sum{|q-BVS|/tot_Atoms}}, with weight: ", wcost(5)
                 Write(unit=lun,fmt="(a,i2,a,f8.4)") &
                 "  => Cost #",i,": Optimization of C6=Sum{qi qj/dij},          with weight: ", wcost(6)

          end select

       end do


       return
    End Subroutine Write_CostFunctPars

    Subroutine Write_FinalCost(lun)
       !---- Arguments ----!
       integer,   intent( in)    :: lun
       !---- Local variables ----!
       character(len=132)   :: line
       real                 :: w,tol
       integer              :: i,ier


       Write(unit=lun,fmt="(/,a)")    "====================================="
       Write(unit=lun,fmt="(a)")      "= Final Cost of minimized functions ="
       Write(unit=lun,fmt="(a,/)")    "====================================="

       do i=1,n_costf

          if(icost(i) == 0) cycle
          select case (i)

              case (1)    !Optimization "Fobs-Fcal"

                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C1=Sum|Fobs-Fcal|/Sum|Fobs|, with weight: ",wcost(i),&
                 "  Final Cost: ",P_cost(1)

              case (2)     !Optimization "dis-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C2=Sum{w(dobs-dcal)^2}, w=1/var(d),with weight: ", wcost(i),&
                 "  Final Cost: ",P_cost(2)


              case (3)     !Optimization "Ang-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C3=Sum{w(Ang_obs-Ang_cal)^2}, w=1/var(Ang),with weight: ", wcost(i),&
                 "  Final Cost: ",P_cost(3)

              case (4)     !Optimization "Tor-restr"
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C4=Sum{w(Tor_obs-Tor_cal)^2}, w=1/var(Tor),with weight: ", wcost(i),&
                 "  Final Cost: ",P_cost(4)

              case (5)     !Optimization "bond-valence"
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C5=Sum{|q-BVS|/tot_Atoms}, with weight: ", wcost(i),&
                 "  Final Cost: ",P_cost(5)

              case (6)     !Optimization "bvs_coulomb"
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C5=Sum{|q-BVS|/tot_Atoms}}, with weight: ", wcost(5),&
                 "  Final Cost: ",P_cost(5)
                 Write(unit=lun,fmt="(a,i2,a,f8.4,a,f12.4)") &
                 "  => Cost #",i,": Optimization of C6=Sum{qi qj/dij},          with weight: ", wcost(6),&
                 "  Final Cost: ",P_cost(6)

          end select

       end do


       return
    End Subroutine Write_FinalCost

    Subroutine General_Cost_function(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer :: i,ic, nop, nlist=1, numv
      integer, dimension(1) :: List




      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      numv=count(abs(v_shift) > 0.00001, dim=1)

      if(numv == 1) then
         i=maxloc(abs(v_shift),dim=1)
         List(1)=v_list(i)
         call VState_to_AtomsPar(A) !Update Atomic parameters with the proper constraints

         cost=0.0

         do ic=1,N_costf

            if(Icost(ic) == 0) cycle

            Select Case(ic)

               case(1)      !Experimental Gobs-Gcalc diffraction pattern
                     call Modify_SF(hkl,A,SpG,List,Nlist,partyp="CO",mode=diff_mode)
                     call Cost_FobsFcal(P_cost(1))
                     cost=cost+ P_cost(1)* WCost(1)

               case(2)      !Distance restraints
                     call Cost_Dis_Rest_Partial(List(1),P_cost(2))
                     cost=cost+ P_cost(2)* WCost(2)

               case(3)      !Angle restraints
                     call Cost_Ang_Rest_Partial(List(1),P_cost(3))
                     cost=cost+ P_cost(3)* WCost(3)

               case(4)      !Torsion angle restraints
                     call Cost_Tor_Rest_Partial(List(1),P_cost(4))
                     cost=cost+ P_cost(4)* WCost(4)

               case(5)      !Bond-Valence
                     call Set_TDist_Coordination(max_coor, Dmax, Cell, SpG, A)
                     call Cost_BVS(Ac, P_cost(5))
                     cost=cost+ P_cost(5)* WCost(5)

               case(6)      !Bond-Valence+ Coulomb repulsion
                     call Set_TDist_Coordination(max_coor, Dmax, Cell, SpG, A)
                     call Cost_BVS_CoulombRep(Ac, P_cost(5),P_cost(6))
                     cost=cost+ P_cost(5)* WCost(5)+P_cost(6)* WCost(6)

               case(7)      !Anti-Bumb functions
               case(8)      !Potential

            End Select
         end do

      else            !New configuration

         call VState_to_AtomsPar(A,mode="V")   !Update Atomic parameters with the proper constraints
         cost=0.0

         do ic=1,N_costf

            if(Icost(ic) == 0) cycle

            Select Case(ic)

               case(1)      !Experimental Gobs-Gcalc diffraction pattern
                     call Structure_Factors(A,SpG,hkl,mode=diff_mode)
                     call Cost_FobsFcal(P_cost(1))
                     cost=cost+ P_cost(1)* WCost(1)

               case(2)      !Distance restraints
                     call Cost_Dis_Rest(P_cost(2))
                     cost=cost+ P_cost(2)* WCost(2)

               case(3)      !Angle restraints
                     call Cost_Ang_Rest(P_cost(3))
                     cost=cost+ P_cost(3)* WCost(3)

               case(4)      !Torsion angle restraints
                     call Cost_Tor_Rest(P_cost(4))
                     cost=cost+ P_cost(4)* WCost(4)

               case(5)      !Bond-Valence
                     call Set_TDist_Coordination(max_coor, Dmax, Cell, SpG, A)
                     call Cost_BVS(Ac, P_cost(5))
                     cost=cost+ P_cost(5)* WCost(5)

               case(6)      !Bond-Valence+ Coulomb repulsion
                     call Set_TDist_Coordination(max_coor, Dmax, Cell, SpG, A)
                     call Cost_BVS_CoulombRep(Ac, P_cost(5),P_cost(6))
                     cost=cost+ P_cost(5)* WCost(5)+P_cost(6)* WCost(6)

               case(7)      !Anti-Bumb functions
               case(8)      !Potential

            End Select
         end do
      end if

      return
    End Subroutine General_Cost_function

    Subroutine Cost_FobsFcal(cost)
       real,                 intent(out):: cost
       !---- Local variables ----!
       integer              :: i,n
       real                 :: delta,sumcal

       n=hkl%Nref
       sumcal=sum(abs(hkl%ref(1:n)%Fc))
       ScaleFact=1.0
       if(sumcal > 0.0000001) ScaleFact=SumGobs/sumcal
       cost=0.0
       do i=1,hkl%Nref
         delta=hkl%ref(i)%Fo-ScaleFact*hkl%ref(i)%Fc
         cost=cost+abs(delta)
       end do
       cost=100.0*cost/SumGobs
       return
    End Subroutine Cost_FobsFcal

    Subroutine Cost_Dis_Rest_Partial(List,cost)
       integer,   intent(in) :: List
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2
       real    :: w, delta
       real, dimension(3) :: x1,x2,tr

       cost=0.0
       do i=1,NP_Rest_Dis
          i1=Dis_rest(i)%p(1)
          i2=Dis_rest(i)%p(2)
          if(list == i1 .or. list == i2) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(dis_rest(i)%stcode,nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr
              Dis_rest(i)%dcalc=distance(x1,x2,cell)
          end if
          delta=Dis_rest(i)%dobs-Dis_rest(i)%dcalc
          w= 1.0/(Dis_rest(i)%sigma*Dis_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do

       return
    End Subroutine Cost_Dis_Rest_Partial

    Subroutine Cost_Dis_Rest(cost)
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2
       real    :: w, delta
       real, dimension(3) :: x1,x2,tr

       cost=0.0
       do i=1,NP_Rest_Dis
          i1=Dis_rest(i)%p(1)
          i2=Dis_rest(i)%p(2)
          x1=A%Atom(i1)%x
          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(dis_rest(i)%stcode,nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr
          Dis_rest(i)%dcalc=distance(x1,x2,cell)
          delta=Dis_rest(i)%dobs-Dis_rest(i)%dcalc
          w= 1.0/(Dis_rest(i)%sigma*Dis_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do
       return
    End Subroutine Cost_Dis_Rest

    Subroutine Cost_Ang_Rest_Partial(List,cost)
       integer,   intent(in) :: List
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2,i3
       real    :: w, delta
       real, dimension(3) :: x1,x2,x3,tr

       cost=0.0

       do i=1,NP_Rest_Ang
          i1=Ang_rest(i)%p(1)
          i2=Ang_rest(i)%p(2)
          i3=Ang_rest(i)%p(3)
          if(list == i1 .or. list == i2 .or. list == i3 ) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(ang_rest(i)%stcode(1),nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr
              x3=A%Atom(i3)%x
              call Read_SymTrans_Code(ang_rest(i)%stcode(2),nop,tr)
              x3=ApplySO(SpG%SymOP(nop),x3)+tr
              Ang_rest(i)%Acalc=Angle_UV(x1-x2,x3-x2,cell%GD)
          end if
          delta=Ang_rest(i)%Aobs-Ang_rest(i)%Acalc
          w= 1.0/(Ang_rest(i)%sigma*Ang_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do

       return
    End Subroutine Cost_Ang_Rest_Partial

    Subroutine Cost_Ang_Rest(cost)
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2,i3
       real    :: w, delta
       real, dimension(3) :: x1,x2,x3,tr

       cost=0.0

       do i=1,NP_Rest_Ang
          i1=Ang_rest(i)%p(1)
          i2=Ang_rest(i)%p(2)
          i3=Ang_rest(i)%p(3)
          x1=A%Atom(i1)%x
          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(ang_rest(i)%stcode(1),nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr
          x3=A%Atom(i3)%x
          call Read_SymTrans_Code(ang_rest(i)%stcode(2),nop,tr)
          x3=ApplySO(SpG%SymOP(nop),x3)+tr
          Ang_rest(i)%Acalc=Angle_UV(x1-x2,x3-x2,cell%GD)
          delta=Ang_rest(i)%Aobs-Ang_rest(i)%Acalc
          w= 1.0/(Ang_rest(i)%sigma*Ang_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do

       return
    End Subroutine Cost_Ang_Rest

    Subroutine Cost_Tor_Rest_Partial(List,cost)
       integer,   intent(in) :: List
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2,i3,i4
       real    :: w, delta
       real, dimension(3) :: x1,x2,x3,x4,tr

       cost=0.0

       do i=1,NP_Rest_tor
          i1=Tor_rest(i)%p(1)
          i2=Tor_rest(i)%p(2)
          i3=Tor_rest(i)%p(3)
          i4=Tor_rest(i)%p(4)
          if(list == i1 .or. list == i2 .or. list == i3 .or. list == i4 ) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(1),nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr

              x3=A%Atom(i3)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(2),nop,tr)
              x3=ApplySO(SpG%SymOP(nop),x3)+tr

              x4=A%Atom(i3)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(3),nop,tr)
              x4=ApplySO(SpG%SymOP(nop),x4)+tr

              x1=matmul(Cell%Cr_Orth_cel,x1)
              x2=matmul(Cell%Cr_Orth_cel,x2)
              x3=matmul(Cell%Cr_Orth_cel,x3)
              x4=matmul(Cell%Cr_Orth_cel,x4)

              tor_rest(i)%Tcalc=Angle_Dihedral(x1,x2,x3,x4)
          end if
          delta=tor_rest(i)%Tobs-tor_rest(i)%Tcalc
          w= 1.0/(Tor_rest(i)%sigma*Tor_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do

       return
    End Subroutine Cost_Tor_Rest_Partial

    Subroutine Cost_Tor_Rest(cost)
       real,      intent(out):: cost
       !---- Local variables ----!
       integer :: i, nop, i1,i2,i3,i4
       real    :: w, delta
       real, dimension(3) :: x1,x2,x3,x4,tr

       cost=0.0

       do i=1,NP_Rest_tor
          i1=Tor_rest(i)%p(1)
          i2=Tor_rest(i)%p(2)
          i3=Tor_rest(i)%p(3)
          i4=Tor_rest(i)%p(4)
          x1=A%Atom(i1)%x
          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(1),nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr

          x3=A%Atom(i3)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(2),nop,tr)
          x3=ApplySO(SpG%SymOP(nop),x3)+tr

          x4=A%Atom(i3)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(3),nop,tr)
          x4=ApplySO(SpG%SymOP(nop),x4)+tr

          x1=matmul(Cell%Cr_Orth_cel,x1)
          x2=matmul(Cell%Cr_Orth_cel,x2)
          x3=matmul(Cell%Cr_Orth_cel,x3)
          x4=matmul(Cell%Cr_Orth_cel,x4)

          tor_rest(i)%Tcalc=Angle_Dihedral(x1,x2,x3,x4)
          delta=tor_rest(i)%Tobs-tor_rest(i)%Tcalc
          w= 1.0/(Tor_rest(i)%sigma*Tor_rest(i)%sigma)
          cost= cost+delta*delta*w
       end do

       return
    End Subroutine Cost_Tor_Rest



!    ==========================================

    Subroutine Cost_function_Restraints(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer :: i, nop, i1,i2,i3,i4, list, numv
      real    :: w, delta
      real, dimension(3) :: x1,x2,x3,x4,tr


      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      numv=count(abs(v_shift) > 0.00001, dim=1)
      if(numv == 1) then
         i=maxloc(abs(v_shift),dim=1)
         List=v_list(i)
         call VState_to_AtomsPar(A)                        !Update Atomic parameters with the proper constraints
      else
         call VState_to_AtomsPar(A,mode="V")               !Update Atomic parameters with the proper constraints
      end if

      cost=0.0

      if(numv == 1) then

         do i=1,NP_Rest_Dis
          i1=Dis_rest(i)%p(1)
          i2=Dis_rest(i)%p(2)
          if(list == i1 .or. list == i2) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(dis_rest(i)%stcode,nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr
              Dis_rest(i)%dcalc=distance(x1,x2,cell)
          end if
          delta=Dis_rest(i)%dobs-Dis_rest(i)%dcalc
          w= 1.0/(Dis_rest(i)%sigma*Dis_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do

         do i=1,NP_Rest_Ang
          i1=Ang_rest(i)%p(1)
          i2=Ang_rest(i)%p(2)
          i3=Ang_rest(i)%p(3)
          if(list == i1 .or. list == i2 .or. list == i3 ) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(ang_rest(i)%stcode(1),nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr
              x3=A%Atom(i3)%x
              call Read_SymTrans_Code(ang_rest(i)%stcode(2),nop,tr)
              x3=ApplySO(SpG%SymOP(nop),x3)+tr
              Ang_rest(i)%Acalc=Angle_UV(x1-x2,x3-x2,cell%GD)
          end if
          delta=Ang_rest(i)%Aobs-Ang_rest(i)%Acalc
          w= 1.0/(Ang_rest(i)%sigma*Ang_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do


         do i=1,NP_Rest_tor
          i1=Tor_rest(i)%p(1)
          i2=Tor_rest(i)%p(2)
          i3=Tor_rest(i)%p(3)
          i4=Tor_rest(i)%p(4)
          if(list == i1 .or. list == i2 .or. list == i3 .or. list == i4 ) then
              x1=A%Atom(i1)%x
              x2=A%Atom(i2)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(1),nop,tr)
              x2=ApplySO(SpG%SymOP(nop),x2)+tr

              x3=A%Atom(i3)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(2),nop,tr)
              x3=ApplySO(SpG%SymOP(nop),x3)+tr

              x4=A%Atom(i3)%x
              call Read_SymTrans_Code(tor_rest(i)%stcode(3),nop,tr)
              x4=ApplySO(SpG%SymOP(nop),x4)+tr

              x1=matmul(Cell%Cr_Orth_cel,x1)
              x2=matmul(Cell%Cr_Orth_cel,x2)
              x3=matmul(Cell%Cr_Orth_cel,x3)
              x4=matmul(Cell%Cr_Orth_cel,x4)

              tor_rest(i)%Tcalc=Angle_Dihedral(x1,x2,x3,x4)
          end if
          delta=tor_rest(i)%Tobs-tor_rest(i)%Tcalc
          w= 1.0/(Tor_rest(i)%sigma*Tor_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do
      else

         do i=1,NP_Rest_Dis
          i1=Dis_rest(i)%p(1)
          i2=Dis_rest(i)%p(2)
          x1=A%Atom(i1)%x
          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(dis_rest(i)%stcode,nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr

          Dis_rest(i)%dcalc=distance(x1,x2,cell)
          delta=Dis_rest(i)%dobs-Dis_rest(i)%dcalc
          w= 1.0/(Dis_rest(i)%sigma*Dis_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do

         do i=1,NP_Rest_Ang
          i1=Ang_rest(i)%p(1)
          i2=Ang_rest(i)%p(2)
          i3=Ang_rest(i)%p(3)
          x1=A%Atom(i1)%x
          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(ang_rest(i)%stcode(1),nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr
          x3=A%Atom(i3)%x
          call Read_SymTrans_Code(ang_rest(i)%stcode(2),nop,tr)
          x3=ApplySO(SpG%SymOP(nop),x3)+tr
          Ang_rest(i)%Acalc=Angle_UV(x1-x2,x3-x2,cell%GD)
          delta=Ang_rest(i)%Aobs-Ang_rest(i)%Acalc
          w= 1.0/(Ang_rest(i)%sigma*Ang_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do


         do i=1,NP_Rest_tor
          i1=Tor_rest(i)%p(1)
          i2=Tor_rest(i)%p(2)
          i3=Tor_rest(i)%p(3)
          i4=Tor_rest(i)%p(4)
          x1=A%Atom(i1)%x

          x2=A%Atom(i2)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(1),nop,tr)
          x2=ApplySO(SpG%SymOP(nop),x2)+tr

          x3=A%Atom(i3)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(2),nop,tr)
          x3=ApplySO(SpG%SymOP(nop),x3)+tr

          x4=A%Atom(i3)%x
          call Read_SymTrans_Code(tor_rest(i)%stcode(3),nop,tr)
          x4=ApplySO(SpG%SymOP(nop),x4)+tr

          x1=matmul(Cell%Cr_Orth_cel,x1)
          x2=matmul(Cell%Cr_Orth_cel,x2)
          x3=matmul(Cell%Cr_Orth_cel,x3)
          x4=matmul(Cell%Cr_Orth_cel,x4)

          tor_rest(i)%Tcalc=Angle_Dihedral(x1,x2,x3,x4)
          delta=tor_rest(i)%Tobs-tor_rest(i)%Tcalc
          w= 1.0/(Tor_rest(i)%sigma*Tor_rest(i)%sigma)
          cost= cost+delta*delta*w
         end do
      end if
      return
    End Subroutine Cost_function_Restraints

    Subroutine Cost_function_FobsFcal(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer :: i, n
      real    :: delta,sumcal


      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      call VState_to_AtomsPar(A)                        !Update Atomic parameters with the proper constraints

      call Structure_Factors(A,SpG,hkl,mode="NUC")

      cost=0.0
      n=hkl%Nref
      sumcal=sum(abs(hkl%ref(1:n)%Fc))
      ScaleFact=1.0
      if(sumcal > 0.0000001) ScaleFact=SumGobs/sumcal
      do i=1,hkl%Nref
       delta=hkl%ref(i)%Fo-ScaleFact*hkl%ref(i)%Fc
       cost=cost+abs(delta)
      end do
      cost=100.0*cost/SumGobs

    End Subroutine Cost_function_FobsFcal

    Subroutine Cost_function_MFobsFcal(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer              :: i,nlist=1,n
      real                 :: delta,sumcal
      integer,dimension(1) :: List


      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      i=maxloc(abs(v_shift),dim=1)
      List(1)=v_list(i)

      call VState_to_AtomsPar(A)                        !Update Atomic parameters with the proper constraints
      call Modify_SF(hkl,A,SpG,List,Nlist,partyp="CO",mode="NUC")

      n=hkl%Nref
      sumcal=sum(abs(hkl%ref(1:n)%Fc))
      ScaleFact=1.0
      if(sumcal > 0.0000001) ScaleFact=SumGobs/sumcal
      cost=0.0
      do i=1,hkl%Nref
       delta=hkl%ref(i)%Fo-ScaleFact*hkl%ref(i)%Fc
       cost=cost+abs(delta)
      end do
      cost=100.0*cost/SumGobs

    End Subroutine Cost_function_MFobsFcal


    Subroutine Cost_function_MulFobsFcal(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer              :: i,nlist=1,n,numv
      real                 :: delta,sumcal
      integer,dimension(1) :: List


      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration

      v_vec(1:NP_Refi)=v(1:NP_Refi)       !Update the State vector

      numv=count(abs(v_shift) > 0.00001, dim=1)
      if(numv == 1) then
         i=maxloc(abs(v_shift),dim=1)
         List(1)=v_list(i)
         call VState_to_AtomsPar(A)                        !Update Atomic parameters with the proper constraints
         call Modify_SF(hkl,A,SpG,List,Nlist,partyp="CO",mode="NUC")
      else
         call VState_to_AtomsPar(A,mode="V")                !Update Atomic parameters with the proper constraints
         call Structure_Factors(A,SpG,hkl,mode="NUC")
      end if

      n=hkl%Nref
      sumcal=sum(abs(hkl%ref(1:n)%Fc))
      ScaleFact=1.0
      if(sumcal > 0.0000001) ScaleFact=SumGobs/sumcal
      cost=0.0
      do i=1,hkl%Nref
       delta=hkl%ref(i)%Fo-ScaleFact*hkl%ref(i)%Fc
       cost=cost+abs(delta)
      end do
      cost=100.0*cost/SumGobs

    End Subroutine Cost_function_MulFobsFcal


    Subroutine LSQ_FobsFcal(v,cost)
      real,dimension(:),    intent( in):: v
      real,                 intent(out):: cost
      !---- Local variables ----!
      integer :: i,n
      real    :: delta,sumcal


      v_shift(1:NP_Refi)=v(1:NP_Refi)-v_vec(1:NP_Refi)  !Calculate the shifts w.r.t. old configuration
      v_vec(1:NP_Refi)=v(1:NP_Refi)

      call VState_to_AtomsPar(A)                        !Update Atomic parameters with the proper constraints

      call Structure_Factors(A,SpG,hkl,mode="NUC")

      n=hkl%Nref
      sumcal=sum(abs(hkl%ref(1:n)%Fc))
      ScaleFact=1.0
      if(sumcal > 0.0000001) ScaleFact=SumGobs/sumcal
      cost=0.0
      do i=1,hkl%Nref
       delta=hkl%ref(i)%Fo-ScaleFact*hkl%ref(i)%Fc
       cost=cost+abs(delta)
      end do
      cost=100.0*cost/SumGobs

    End Subroutine LSQ_FobsFcal

  End Module cost_functions
