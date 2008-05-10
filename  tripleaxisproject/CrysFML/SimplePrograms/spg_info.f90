!-----------------------------------------
!  Example of simple program using CrysFML
!-----------------------------------------
  Program Get_SPG_info
    use crystallographic_symmetry, only:  &
        space_group_type,set_spacegroup, write_spacegroup
    character(len=20)      :: spg_symb
    type(space_group_type) :: SPG

    do
      write(unit=*,fmt="(a)",advance="no") &
      " => Please enter a space group (H-M/Hall/number): "
      read(unit=*,fmt="(a)") spg_symb
      if(len_trim(spg_symb) == 0) exit
      call set_spacegroup(spg_symb,SPG)
      call write_spacegroup(SPG,full=.true.)
    end do
    stop
  End Program Get_SPG_info
