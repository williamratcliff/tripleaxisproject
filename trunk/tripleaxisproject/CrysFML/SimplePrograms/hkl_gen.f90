!---------------------------------------
!  Example of simple program using CFML
!---------------------------------------
  Program Test_HKL_GEN
   use Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                        Write_SpaceGroup
   use Crystal_types, only: Crystal_Cell_Type,set_crystal_cell,Write_Crystal_Cell
   use Reflections_utilities, only: Reflect_Type, get_maxnumref, HKL_uni !HKL_gen

   implicit none

        character(len=1)           :: car
        character(len=12)          :: name_file
        character(len=20)          :: spgr
        character(len=35)          :: texto
        type (Space_Group_Type)    :: grp_espacial
        type (Crystal_Cell_Type)   :: cell
        integer                    :: i,num , ier, MaxNumRef
        real,dimension(3)          :: celda, angulo
        real                       :: val1,val2
        logical                    :: info

        type (Reflect_Type),allocatable, dimension(:) :: reflections
        info=.true.
        num=0
        do
           call system("cls")
           write(unit=*,fmt="(a)") "     List of unique reflections     "
           write(unit=*,fmt="(a)") " ==================================="
           write(unit=*,fmt="(a)") " "

           if (info) then
              write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
              read(unit=*,fmt="(a)") spgr
              if (len_trim(spgr)==0) exit
              call set_spacegroup(spgr,grp_espacial)
              write(unit=*,fmt="(a)",advance="no") " => Cell (6 reals): "
              read(unit=*,fmt=*,iostat=ier) celda(1),celda(2),celda(3),angulo(1),angulo(2),angulo(3)
              if( ier /= 0 ) cycle
              call set_crystal_cell(celda,angulo,cell)
              info=.false.
           end if
           write(unit=*,fmt=*) " "
           write(unit=*,fmt="(a)",advance="no") " Interval in Sin_Theta/Lambda (2 reals, val_1/2=0 => stops): "
           read(unit=*,fmt=*,iostat=ier) val1,val2
           if(ier /= 0) cycle
           if (val1*val2 < 0.0001) exit
           texto = "  1/Angtroms (Sin_Theta/Lambda)"
           car="s"

           MaxNumRef = get_maxnumref(val2,Cell%CellVol,mult=grp_espacial%Multip)
           write(unit=*,fmt="(a,i10)") " => Maximum number of reflections: ", MaxNumRef
           if (allocated(reflections)) then
             deallocate(reflections)
             allocate (reflections(MaxNumRef))
           else
             allocate (reflections(MaxNumRef))
           end if

        !   call HKL_GEN(cell,grp_espacial,.true.,val1,val2,num,reflections) !Not ordered
           call HKL_UNI(cell,grp_espacial,.true.,val1,val2,car,num,reflections) !Ordered

           write(unit=*,fmt="(a)",advance="no") " Name of the output file: "
           read(unit=*,fmt="(a)") name_file

           open(unit=1,file=trim(name_file),status="REPLACE",action="WRITE")
              write(unit=1,fmt="(a)")   "    LIST OF UNIQUE REFLECTIONS"
              write(unit=1,fmt="(a,/)") "    =========================="
              call Write_SpaceGroup(grp_espacial,1)
              call Write_Crystal_Cell(Cell,1)
              write(unit=1,fmt="(/,a,2f8.4,a,/)") " => List of reflections within: ",val1,val2,texto
           do i=1,num
              write(unit=1,fmt="(3i4,i5,f10.5,i8)") reflections(i)%h, reflections(i)%mult, &
                                                    reflections(i)%S, i
           end do
           close(unit=1)

        end do
    stop
  End Program Test_HKL_GEN
