  Module observed_reflections
    Use Math_gen,                 only: sp,sort
    use Reflections_Utilities,    only: Reflection_List_Type, hkl_mult, hkl_equiv,hkl_s
    use String_Utilities,         only: number_lines, u_Case, get_logunit
    use crystallographic_symmetry,only: Space_Group_Type
    use crystal_types,            only: Crystal_Cell_Type

    implicit none

    private

    public :: Read_observations, Write_ObsCalc_SFactors


    logical,             public ::    err_observ=.false.
    character (len=256), public ::    err_mess_observ="  "
    Real,                public ::    SumGobs, ScaleFact

    !!----
    !!---- TYPE :: OBSERVATION_TYPE
    !!--..
    !!---- Type, public :: Observation_Type
    !!----    integer                 :: ncont ! number of contributing reflections
    !!-----   integer,dimension(3,12) :: Hd    ! Indices of contributing reflections
    !!-----   integer,dimension(  12) :: icod  ! Pointer to the number of domain
    !!-----   integer,dimension(  12) :: p     ! Pointer to the list of independent reflections
    !!----    real(kind=sp)           :: Gobs  ! Observed Structure Factor squared
    !!----    real(kind=sp)           :: SGobs ! Sigma of  Gobs
    !!---- End Type Observation_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Observation_Type
       integer                            :: ncont ! number of contributing reflections
       integer,allocatable,dimension(:,:) :: Hd    ! (3,ncont)Indices of contributing reflections
       integer,allocatable,dimension(:)   :: icod  ! (ncont)  Pointer to the number of domain
       integer,allocatable,dimension(:)   :: p     ! (ncont)  Pointer to the list of independent reflections
       real(kind=sp)                      :: Gobs  ! Observed Structure Factor squared
       real(kind=sp)                      :: SGobs ! Sigma of  Gobs
    End Type Observation_Type

    !!----
    !!---- TYPE :: OBSERVATION_LIST_TYPE
    !!--..
    !!---- Type, public :: Observation_List_Type
    !!----    integer                                         :: Nobs  ! Number of Observations
    !!----    type(observation_type),allocatable,dimension(:) :: Ob    ! Observation (F2, etc...)
    !!---- End Type Observation_List_Type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: Observation_List_Type
       integer                                         :: Nobs  ! Number of Observations
       type(observation_type),allocatable,dimension(:) :: Ob    ! Observation (F2, etc...)
    End Type Observation_List_Type

    Interface Read_observations
       Module Procedure Read_observations_clusters
       Module Procedure Read_observations_reflections
    End Interface Read_observations



    Contains

    !!----
    !!---- Subroutine Read_observations_reflections(file_hkl,Cell,Spg,Friedel,Rf)
    !!----   character(len=*),            intent (in)  :: file_hkl  !Name of the hkl-file
    !!----   type(Crystal_Cell_Type),     intent (in)  :: Cell
    !!----   type(Space_Group_Type),      intent (in)  :: Spg
    !!----   logical,                     intent (in)  :: Friedel
    !!----   type(Reflection_List_Type),  intent (out) :: Rf
    !!----
    !!----   Subroutine for reading a list of independent reflections
    !!----   The unit cell and space group need to be set before calling this subroutine.
    !!----   Construct the Reflection_List_Type variable "Rf", containing the experimental
    !!----   data and the value of hkl, s, etc. The reflections are reordered by increasing
    !!----   sinTheta/Lambda before return.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Read_observations_reflections(file_hkl,Cell,Spg,Friedel,Rf)
      character(len=*),            intent (in)  :: file_hkl
      type(Crystal_Cell_Type),     intent (in)  :: Cell
      type(Space_Group_Type),      intent (in)  :: Spg
      logical,                     intent (in)  :: Friedel
      type(Reflection_List_Type),  intent (out) :: Rf
      !--- Local Variables ---!
      integer :: ier, i,j,nlines,i_hkl, itemp
      integer, allocatable,dimension(:,:) :: hkl
      integer, allocatable,dimension(  :) :: ic
      real,    allocatable,dimension(  :) :: sv,fobs,sigma



      call get_logunit(i_hkl)
      open(unit=i_hkl, file=trim(file_hkl), status="old", action="read",position="rewind",iostat=ier)
      if(ier /= 0) then
        err_observ=.true.
        err_mess_observ="  Error opening the file: "//trim(file_hkl)
        return
      end if

      call number_lines(trim(file_hkl),nlines)

      !Allocate local types to the maximum possible value

      if(allocated(sv)) deallocate(sv)
      allocate(sv(nlines))
      if(allocated(fobs)) deallocate(fobs)
      allocate(fobs(nlines))
      if(allocated(sigma)) deallocate(sigma)
      allocate(sigma(nlines))
      if(allocated(ic)) deallocate(ic)
      allocate(ic(nlines))
      if(allocated(hkl)) deallocate(hkl)
      allocate(hkl(3,nlines))

      do i=1,nlines
          read(unit=i_hkl,fmt=*, iostat=ier) hkl(:,i),fobs(i),sigma(i)
          if(ier /= 0) exit
      end do
      close(unit=i_hkl)

      !Now ordering of all reflections
      do i=1,nlines
       sv(i)=hkl_s(hkl(:,i),Cell)
      end do
      call sort(sv,nlines,ic) !use ic for pointer ordering

      call get_logunit(itemp)
      open(unit=itemp,status="scratch",form="unformatted",action="readwrite")
      do i=1,nlines
        j=ic(i)
        write(itemp) hkl(:,j),sv(j),fobs(j),sigma(j)
      end do
      rewind(unit=itemp)
      do i=1,nlines
        read(itemp) hkl(:,i),sv(i),fobs(i),sigma(i)
      end do

      if(allocated(Rf%Ref)) deallocate(Rf%Ref)
      allocate(Rf%Ref(nlines))
      Rf%Nref=nlines
      SumGobs=0.0
      do i=1,nlines
        SumGobs=SumGobs+abs(fobs(i))
        Rf%Ref(i)%h    = hkl(:,i)
        Rf%Ref(i)%Mult = hkl_mult(hkl(:,i),Spg,Friedel)
        Rf%Ref(i)%Fc   = 0.0
        Rf%Ref(i)%Fo   = fobs(i)
        Rf%Ref(i)%SFo  = sigma(i)
        Rf%Ref(i)%S    = sv(i)
        Rf%Ref(i)%W    = 0.0
        Rf%Ref(i)%Phase= 0.0
        Rf%Ref(i)%A    = 0.0
        Rf%Ref(i)%B    = 0.0
        Rf%Ref(i)%AA   = 0.0
        Rf%Ref(i)%BB   = 0.0
      end do
      return
    End Subroutine Read_observations_reflections


    !!----
    !!---- Subroutine Read_observations_clusters(file_hkl,Cell,Spg,Friedel,Obs,Rf)
    !!----   character(len=*),            intent (in)  :: file_hkl  !Name of the hkl-file
    !!----   type(Crystal_Cell_Type),     intent (in)  :: Cell
    !!----   type(Space_Group_Type),      intent (in)  :: Spg
    !!----   logical,                     intent (in)  :: Friedel
    !!----   type(Observation_List_Type), intent (out) :: Obs
    !!----   type(Reflection_List_Type),  intent (out) :: Rf
    !!----
    !!----   Subroutine for reading a list of observations (cluster of reflections)
    !!----   The unit cell and space group need to be set before calling this subroutine.
    !!----   Construct the Observation_List_Type variable "Obs", containing the experimental
    !!----   data and pointers to the strictly independent reflections that are gathered in
    !!----   the intent out Reflection_List_Type variable "Rf". This last variable is suited
    !!----   for using the calculation of structure factors.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Read_observations_clusters(file_hkl,Cell,Spg,Friedel,Obs,Rf)
      character(len=*),            intent (in)  :: file_hkl
      type(Crystal_Cell_Type),     intent (in)  :: Cell
      type(Space_Group_Type),      intent (in)  :: Spg
      logical,                     intent (in)  :: Friedel
      type(Observation_List_Type), intent (out) :: Obs
      type(Reflection_List_Type),  intent (out) :: Rf
      !--- Local Variables ---!
      integer :: ier, i,j,k,nlines,i_hkl,icod, nover, itemp, nri, iobs,nr,n
      integer, dimension(3)               :: h
      integer, allocatable,dimension(  :) :: ic
      integer, allocatable,dimension(:,:) :: hov
      integer, allocatable,dimension(:,:) :: hkl
      integer, allocatable,dimension(:,:) :: point
      real,    allocatable,dimension(  :) :: sv,sq
      real                                :: fobs,sigma
      type(Observation_List_Type):: O



      call get_logunit(i_hkl)
      open(unit=i_hkl, file=trim(file_hkl), status="old", action="read",position="rewind",iostat=ier)
      if(ier /= 0) then
        err_observ=.true.
        err_mess_observ="  Error opening the file: "//trim(file_hkl)
        return
      end if

      call number_lines(trim(file_hkl),nlines)

      !Allocate local types to the maximum possible value

      if(allocated(O%Ob)) deallocate(O%Ob)
      allocate(O%Ob(nlines))
      if(allocated(ic)) deallocate(ic)
      allocate(ic(nlines))
      if(allocated(sv)) deallocate(sv)
      allocate(sv(nlines))
      if(allocated(sq)) deallocate(sq)
      allocate(sq(nlines))
      if(allocated(hov)) deallocate(hov)
      allocate(hov(3,nlines))
      if(allocated(hkl)) deallocate(hkl)
      allocate(hkl(3,nlines))
      if(allocated(point)) deallocate(point)
      allocate(point(2,nlines))


      iobs=0
      nr=0
      do_lines:do i=1,nlines
        nover=0
        do
          read(unit=i_hkl,fmt=*, iostat=ier) h,fobs,sigma,icod
          if(ier /= 0) exit do_lines
            nr          = nr+1
            nover       = nover+1
            point(1,nr) = iobs
            point(2,nr) = nover
          if(fobs < 0 ) then
            hov(:,nover)= h
            hkl(:,nr)   = h
            ic(nover)   = icod
          else
            hov(:,nover)= h
            hkl(:,nr)   = h
            iobs=iobs+1
            O%Ob(iobs)%Gobs  = fobs
            O%Ob(iobs)%SGobs = sigma
            O%Ob(iobs)%ncont = nover
            allocate(O%Ob(iobs)%hd(3,nover))

            do j=1,nover
              O%Ob(iobs)%hd(:,j) = hov(:,j)
              O%Ob(iobs)%icod(j) = ic(j)
            end do
            exit
          end if
        end do
      end do do_lines
      close(unit=i_hkl)
      Obs%Nobs=iobs
      if(allocated(Obs%Ob)) deallocate(Obs%Ob)
      allocate(Obs%Ob(iobs))

      do i=1,iobs
        Obs%Ob(i)%Gobs =O%Ob(i)%Gobs
        Obs%Ob(i)%SGobs=O%Ob(i)%SGobs
        Obs%Ob(i)%ncont=O%Ob(i)%ncont
        n=O%Ob(i)%ncont
        allocate(Obs%Ob(i)%hd(3,n))
        allocate(Obs%Ob(i)%icod(n))
        allocate(Obs%Ob(i)%p(n))
        do j=1,n
         Obs%Ob(i)%hd(:,j)=O%Ob(i)%hd(:,j)
         Obs%Ob(i)%icod(j)=O%Ob(i)%icod(j)
        end do
      end do

      !Now ordering of all reflections read and put the  equivalents
      do i=1,nr
       sv(i)=hkl_s(hkl(:,i),Cell)
      end do
      call sort(sv,nr,ic) !use ic for pointer ordering

      call get_logunit(itemp)
      open(unit=itemp,status="scratch",form="unformatted",action="readwrite")
      do i=1,nr
        j=ic(i)
        write(itemp) hkl(:,j),sv(j),point(:,j)
      end do
      rewind(unit=itemp)
      do i=1,nr
        read(itemp) hkl(:,i),sv(i),point(:,i)
      end do
      ic=0 !nullify these vector for another use
      hov=0
      i=0
      nri=0
      do
        i=i+1
        if(i >= nr) exit
        iobs = point(1,i)
        nover= point(2,i)
        if(ic(i) == 0) then
          ic(i) = 1
          Obs%Ob(iobs)%p(nover)=i
          nri=nri+1
          hov(:,nri) = hkl(:,i)
           sq(  nri) = sv(i)
        end if

        do j=i+1,nr
          iobs = point(1,j)
          nover= point(2,j)
          if( hkl_equiv(hkl(:,j),hkl(:,i),Spg,Friedel) ) then
            Obs%Ob(iobs)%p(nover)=i
            ic(j)=1
          else
            Obs%Ob(iobs)%p(nover)=j
            ic(j)=1
            i=j
            nri=nri+1
            hov(:,nri) = hkl(:,j)
            sq(  nri) = sv(j)
            exit
          end if
        end do

      end do

      if(allocated(Rf%Ref)) deallocate(Rf%Ref)
      allocate(Rf%Ref(nri))
      Rf%Nref=nri
      do i=1,nri
        Rf%Ref(i)%h    = hov(:,i)
        Rf%Ref(i)%Mult = hkl_mult(hov(:,i),Spg,Friedel)
        Rf%Ref(i)%Fo   = 0.0
        Rf%Ref(i)%Fc   = 0.0
        Rf%Ref(i)%SFo  = 0.0
        Rf%Ref(i)%S    = sq(i)
        Rf%Ref(i)%W    = 0.0
        Rf%Ref(i)%Phase= 0.0
        Rf%Ref(i)%A    = 0.0
        Rf%Ref(i)%B    = 0.0
        Rf%Ref(i)%AA   = 0.0
        Rf%Ref(i)%BB   = 0.0
      end do
      return
    End Subroutine Read_observations_clusters
    !!----
    !!---- Subroutine Write_ObsCalc_SFactors(lun,Reflex,Mode)
    !!----    integer,                            intent(in) :: lun
    !!----    type(reflection_list_type),         intent(in) :: Reflex
    !!----    Character(len=*), optional,         intent(in) :: Mode
    !!----
    !!----    Writes in logical unit=lun the list of observed versus
    !!----    calculated structure factors contained in the array hkl
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Write_ObsCalc_SFactors(lun,Reflex,Mode)
       !---- Argument ----!
       integer,                            intent(in) :: lun
       type(reflection_list_type),         intent(in) :: Reflex
       Character(len=*), optional,         intent(in) :: Mode
       !---- Local Variables ----!
       integer :: i
       real    :: delta, deltaSig,R

       If(present(mode)) then
         Select Case (mode(1:3))
           Case("NUC","nuc")
             write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(NEUTRONS)"
             write(unit=lun,fmt="(a)")     "    ==================================================="
           Case default
             write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
             write(unit=lun,fmt="(a)")     "    ================================================="
         End Select
       else
         write(unit=lun,fmt="(a)")   "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
         write(unit=lun,fmt="(a)")   "    ================================================="
       end if

       write(unit=lun,fmt="(/,a,/)") "   H   K   L   Mult  SinTh/Lda    |Fobs|      SFobs        |Fc|       Delta     Delta/Sigma"
       R=0.0
       do i=1,reflex%Nref
        delta=reflex%ref(i)%Fo-ScaleFact*reflex%ref(i)%Fc
        R=R+abs(delta)
        write(unit=lun,fmt="(3i4,i5,6f12.5)") reflex%ref(i)%h, reflex%ref(i)%mult,     &
              reflex%ref(i)%S, reflex%ref(i)%Fo/ScaleFact,reflex%ref(i)%SFo, reflex%ref(i)%Fc,   &
              delta, delta/max(0.0001,reflex%ref(i)%SFo)
       end do
       R=100.0*R/sum(reflex%ref(1:reflex%Nref)%Fo)
       write(unit=lun,fmt="(a,f12.5)") "  =>  R-Factor(%) (Sum {|Fo-Fc|}/Sum{|Fo|} = ",R
       return
    End Subroutine Write_ObsCalc_SFactors

  End Module observed_reflections
