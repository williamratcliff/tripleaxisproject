!-----------------------------------------
!  Example of simple program using CrysFML
!-----------------------------------------

  Program Generate_T_SubGroups
    use crystallographic_symmetry, only:  get_T_SubGroups,space_group_type, &
                set_spacegroup, write_spacegroup, Lattice_trans,similar_transf_SG
    use math_3d, only : determ_a

    character(len=20)      :: spg_symb
    real, dimension (3,3)  :: trans
    real, dimension (3  )  :: orig
    integer                :: nsg, i, j, ng, l
    real                   :: det
    type(space_group_type) :: SpG,SpGn
    type(space_group_type), dimension(48) :: Subgroup

    do
      write(unit=*,fmt="(a)")          "    "
      write(unit=*,fmt="(a)")          "   ----------------------------- "
      write(unit=*,fmt="(a)")          "   Program: Generate_T_SubGroups "
      write(unit=*,fmt="(a)")          "   ----------------------------- "
      write(unit=*,fmt="(a)")          "    "
      write(unit=*,fmt="(a)",advance="no")  " => Please enter a space group (H-M/Hall/number): "
      read(unit=*,fmt="(a)") spg_symb
      if(len_trim(spg_symb) == 0) exit

      call set_spacegroup(spg_symb,SPG)       !Constructing the space group SPG
      write(unit=*,fmt="(a)",advance="no")  " => Please enter a transformation matrix: "
      read(unit=*,fmt=*) (trans(i,:),i=1,3)
      det=determ_a(trans)
      det=abs(det)
      write(unit=*,fmt="(a)",advance="no")  " => Please enter the new origin: "
      read(unit=*,fmt=*) orig

      call similar_transf_SG(trans,orig,SpG,SpGn)  !Construct the subgroup of SPG that is compatible
      call write_spacegroup(SPGn,full=.true.)      !with the transformation matrix and change of origin
                                                   !give above

     !Determine all subgroups of the new space group
      call get_T_SubGroups(SPGn,SubGroup,nsg)

      write(unit=*,fmt="(/,a,/)") " => LIST of Translationengleiche Subgroups: "
      do i=1,nsg
        j=SPGn%Multip/SubGroup(i)%multip
        ng=SubGroup(i)%numops
        write(unit=*,fmt="(4a,i2,30a)") " => ", SubGroup(i)%Spg_Symb, SubGroup(i)%hall,&
          " Index: [",j,"]   ->  { ", (trim(SubGroup(i)%SymopSymb(l))//" : ",l=1,ng-1),&
          trim(SubGroup(i)%SymopSymb(ng))," }    ", trim(SubGroup(i)%centre)
      end do
    end do
    stop
  End Program Generate_T_SubGroups


