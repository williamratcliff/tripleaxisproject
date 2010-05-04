!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: FFT_GEN
!!----   INFO: FFT Calculations Routines
!!----
!!---- HISTORY
!!--..    Update: January - 2005
!!--..
!!--..    23 - November - 2000  Updated by JGP
!!--..
!!--..    Multivariate Fast Fourier Transform
!!--..    Fortran 90 (ELF90) Implementation of Singleton's mixed-radix
!!--..    algorithm, RC Singleton, Stanford Research Institute, Sept. 1968.
!!--..    Adapted from fftn.c, translated from Fortran 66 to C by Mark Olesen
!!--..    and John Beale.
!!--..    Fourier transforms can be computed either in place, using assumed
!!--..    size arguments, or by generic function, using assumed shape arguments.
!!--..
!!--..    FFT(array, dim, inv)                 generic transform function
!!--..    COMPLEX(fftkind), DIMENSION(:,...,:), INTENT(IN)           :: array
!!--..    INTEGER,          DIMENSION(:),       INTENT(IN),  OPTIONAL:: dim
!!--..    LOGICAL,                              INTENT(IN),  OPTIONAL:: inv
!!--..
!!--..    Formal Parameters:
!!--..
!!--..    array    The complex array to be transformed. array can be of arbitrary
!!--..             rank (i.e. up to seven).
!!--..
!!--..    shape    With subroutine fftn, the shape of the array to be transformed
!!--..             has to be passed separately, since fftradix - the internal trans-
!!--..             formation routine - will treat array always as one dimensional.
!!--..             The product of elements in shape must be the number of
!!--..             elements in array.
!!--..             Although passing array with assumed shape would have been nicer,
!!--..             I prefered assumed size in order to prevent the compiler from
!!--..             using a copy-in-copy-out mechanism. That would generally be
!!--..             necessary with fftn passing array to fftradix and with fftn
!!--..             being prepared for accepting non consecutive array sections.
!!--..             Using assumed size, it's up to the user to pass an array argu-
!!--..             ment, that can be addressed as continous one dimensional array
!!--..             without copying. Otherwise, transformation will not really be
!!--..             performed in place.
!!--..             On the other hand, since the rank of array and the size of
!!--..             shape needn't match, fftn is appropriate for handling more than
!!--..             seven dimensions.
!!--..             As far as function fft is concerned all this doesn't matter,
!!--..             because the argument will be copied anyway. Thus no extra
!!--..             shape argument is needed for fft.
!!--..
!!--..    Optional Parameters:
!!--..
!!--..    dim      One dimensional integer array, containing the dimensions to be
!!--..             transformed. Default is (/1,...,N/) with N being the rank of
!!--..             array, i.e. complete transform. dim can restrict transformation
!!--..             to a subset of available dimensions. Its size must not exceed the
!!--..             rank of array or the size of shape respectivly.
!!--..
!!--..    inv      If .true., inverse transformation will be performed. Default is
!!--..             .false., i.e. forward transformation.
!!--..
!!--..    stat     If present, a system dependent nonzero status value will be
!!--..             returned in stat, if allocation of temporary storage failed.
!!--..             For functions, the integer variable status is used.
!!--..
!!--..    Scaling:
!!--..             Transformation results will always be scaled by the square
!!--..             root of the product of sizes of each dimension in dim.
!!--..             (See examples below)
!!--..
!!--..    Examples:
!!--..             Let A be a L*M*N three dimensional complex array. Then
!!--..             result = fft(A) will produce a three dimensional transform,
!!--..             scaled by sqrt(L*M*N), while call fftn(A, SHAPE(A)) will do
!!--..             the same in place.
!!--..
!!--..             result = fft(A, dim=(/1,3/)) will transform with respect to
!!--..             the first and the third dimension, scaled by sqrt(L*N).
!!--..
!!--..             result = fft(fft(A), inv=.true.) should (approximately)
!!--..             reproduce A. With B having the same shape as A
!!--..
!!--..             result = fft(fft(A) * CONJG(fft(B)), inv=.true.) will
!!--..             correlate A and B.
!!--..
!!--..     Remarks:
!!--..
!!--..             Following changes have been introduced with respect to fftn.c:
!!--..             - complex arguments and results are of type complex, rather
!!--..               than real an imaginary part separately.
!!--..             - increment parameter (magnitude of isign) has been dropped,
!!--..               inc is always one, direction of transform is given by inv.
!!--..             - maxf and maxp have been dropped. The amount of temporary
!!--..               storage needed is determined by the fftradix routine.
!!--..               Both fftn and fft can handle any size of array. (Maybe they
!!--..               take a lot of time and memory, but they will do it)
!!--..
!!--..             Redesigning fftradix in a way, that it handles assumed shape
!!--..             arrays would have been desirable. However, I found it rather
!!--..             hard to do this in an efficient way. Problems were:
!!--..             - to prevent stride multiplications when indexing arrays. At
!!--..               least our compiler was not clever enough to discover that
!!--..               in fact additions would do the job as well. On the other
!!--..               hand, I haven't been clever enough to find an implementation
!!--..               using array operations.
!!--..             - fftradix is rather large and different versions would be
!!--..               necessaray for each possible rank of array.
!!--..             Consequently, in place transformation still needs the
!!--..             argument stored in a consecutive bunch of memory and can't be
!!--..             performed on array sections like A(100:199:-3, 50:1020).
!!--..             Calling fftn with such sections will most probably imply
!!--..             copy-in-copy-out. However, the function fft works with
!!--..             everything it gets and should be convenient to use.
!!--..
!!--..             To enable this module to be used with ELF90 it appears to be
!!--..             necessary to allocate a 1-D work array into which the
!!--..             multi-dimensional array is copied, and then to copy the
!!--..             results back from the 1-D array to the multi-dimensional
!!--..             array ft.
!!--..
!!--..             Unfortunately, ELF90 will not allow a function to return more
!!--..             than one output variable.   The variable `stat' has been
!!--..             dropped from the function arguments. Users should examine the
!!--..             value of the variable `StatusF' instead. This is a PUBLIC
!!--..             variable declared in this module.
!!--..
!!--..             Michael Steffens, 09.12.96,
!!--..             <Michael.Steffens@mbox.muk.uni-hannover.de>
!!--..
!!--..             ELF90-compatible version by Alan Miller, 29 April 1997
!!--..             Alan.Miller @ mel.dms.csiro.au
!!--..
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!--++    FFTKIND                  [Private]
!!--++    STATUSF                  [Private]
!!--++    COS72                    [Private]
!!--++    SIN72                    [Private]
!!--++    SIN60                    [Private]
!!--++    PI                       [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       FFT
!!--++       FFT1D                 [Overloaded]
!!--++       FFT2D                 [Overloaded]
!!--++       FFT3D                 [Overloaded]
!!--++       FFT4D                 [Overloaded]
!!--++       FFT5D                 [Overloaded]
!!--++       FFT6D                 [Overloaded]
!!--++       FFT7D                 [Overloaded]
!!----
!!----    Subroutines:
!!--++       FFTN                  [Private]
!!--++       FFTRADIX
!!--..       FACTORIZE
!!--..       PERMUTE
!!--..       TRANSFORM
!!----
!!
 Module Fft_Gen

    !---- Use Modules ----!

    !---- Local Variables ----!
    implicit None

    private

    !---- List of public variables ----!

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!
    public :: fft

    !---- List of public subroutines ----!

    !---- List of public overloaded procedures: subroutines ----!

    !---- List of private functions ----!
    private :: Fft1D, Fft2D, Fft3D, Fft4D, Fft5D, Fft6D, FFt7D

    !---- List of private subroutines ----!
    private :: fftn

    !---- Definitions ----!

    !!--++
    !!--++ FFTKIND
    !!--++    integer, private, parameter:: Fftkind = Kind(0.0)
    !!--++
    !!--++    (PRIVATE)
    !!--++    Default Precicion for FFT Variables
    !!--++
    !!--++ Update: February - 2005
    !!
    integer, private, parameter:: Fftkind = Kind(0.0) !--- Adjust Here For Other Precisions

    !!--++
    !!--++ STATUSF
    !!--++    integer, private, Save     :: StatusF
    !!--++
    !!--++    Information on FFT Routines
    !!--++
    !!--++ Update: February - 2005
    !!
    integer, private, Save     :: StatusF    !--- Shifted To Here As Elf90 Does Not Allow
                                             !    Arguments To Be Intent(Out)
    !!--++
    !!--++ COS72
    !!--++    Real(Fftkind), Parameter:: Cos72
    !!--++
    !!--++    (PRIVATE)
    !!--++    Cos72 = 0.30901699437494742_Fftkind
    !!--++
    !!--++ Update: February - 2005
    !!
    Real(Fftkind), Parameter:: Cos72 = 0.30901699437494742_Fftkind

    !!--++
    !!--++ SIN72
    !!--++    Real(Fftkind), Parameter:: Sin72
    !!--++
    !!--++    (PRIVATE)
    !!--++    Sin72 = 0.95105651629515357_Fftkind
    !!--++
    !!--++ Update: February - 2005
    !!
    Real(Fftkind), Parameter:: Sin72 = 0.95105651629515357_Fftkind

    !!--++
    !!--++ SIN60
    !!--++    Real(Fftkind), Parameter:: Sin60
    !!--++
    !!--++    (PRIVATE)
    !!--++    Sin60 = 0.86602540378443865_Fftkind
    !!--++
    !!--++ Update: February - 2005
    !!
    Real(Fftkind), Parameter:: Sin60 = 0.86602540378443865_Fftkind

    !!--++
    !!--++ PI
    !!--++    Real(Fftkind), Parameter:: Pi
    !!--++
    !!--++    (PRIVATE)
    !!--++    Pi    = 3.14159265358979323_Fftkind
    !!--++
    !!--++ Update: February - 2005
    !!
    Real(Fftkind), Parameter:: Pi    = 3.14159265358979323_Fftkind

    !---- Interfaces - Overlapp ----!
    Interface Fft
       Module Procedure Fft1D
       Module Procedure Fft2D
       Module Procedure Fft3D
       Module Procedure Fft4D
       Module Procedure Fft5D
       Module Procedure Fft6D
       Module Procedure Fft7D
    End Interface

 Contains

    !---- Functions ----!
    !!----
    !!---- Function Fft(Array, Dim, Inv) Result(Ft)
    !!----    complex(fftkind), dimension(:), intent(in)            :: array  !  In -> Complex array
    !!----    integer,          dimension(:), intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!----                                                                             to be transformed
    !!----    logical,                        intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!----                                                                             Default is .false., i.e. forward transformation.
    !!----    Calculation of FFT from 1 to up 7 dimensions
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Fft1D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:), intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                             to be transformed
    !!--++    logical,                        intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                             Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT one dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft1D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:), intent(in)           :: array
       integer,          dimension(:), intent(in),  optional:: dim
       logical,                        intent(in),  optional:: inv

       !--- function result
       complex(fftkind), dimension(size(array, 1)):: ft
       ft = array
       call fftn(ft, shape(array), dim, inv = inv, stat = StatusF)

       return
    End Function Fft1D

    !!--++
    !!--++ Function Fft2D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),   intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                               to be transformed
    !!--++    logical,                          intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                               Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT two dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft2D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:), intent(in)           :: array
       integer,          dimension(:),   intent(in),  optional:: dim
       logical,                          intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension(size(array, 1), size(array, 2)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft2D

    !!--++
    !!--++ Function Fft3D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),     intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                 to be transformed
    !!--++    logical,                            intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                 Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT three dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft3D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:), intent(in)           :: array
       integer,          dimension(:),     intent(in),  optional:: dim
       logical,                            intent(in),  optional:: inv
       !--- function result
       complex(fftkind), &
            dimension(size(array, 1), size(array, 2), size(array, 3)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft3D

    !!--++
    !!--++ Function Fft4D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:), intent(in)            :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT four dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft4D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:), intent(in)           :: array
       integer,          dimension(:),       intent(in),  optional:: dim
       logical,                              intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4) /))
       if (allocated(work)) deallocate(work)

       return
    End Function Fft4D

    !!--++
    !!--++ Function Fft5D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:), intent(in)          :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT five dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft5D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),         intent(in),  optional:: dim
       logical,                                intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft5D

    !!--++
    !!--++ Function Fft6D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:,:), intent(in)        :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (Overloaded)
    !!--++    FFT six dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft6D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),           intent(in),  optional:: dim
       logical,                                  intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5), size(array, 6) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft6D

    !!--++
    !!--++ Function Fft7D(Array, Dim, Inv) Result(Ft)
    !!--++    complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)      :: array  !  In -> Complex array
    !!--++    integer,          dimension(:),       intent(in),  optional :: dim    !  In -> array containing the dimensions
    !!--++                                                                                   to be transformed
    !!--++    logical,                              intent(in),  optional :: inv    !  In -> If .true., inverse transformation will be performed.
    !!--++                                                                                   Default is .false., i.e. forward transformation.
    !!--++    (OVERLOADED)
    !!--++    FFT seven dimension
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fft7D(Array, Dim, Inv) Result(Ft)
       !--- formal parameters
       complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)           :: array
       integer,          dimension(:),             intent(in),  optional:: dim
       logical,                                    intent(in),  optional:: inv
       !--- function result
       complex(fftkind), dimension( &
            size(array, 1), size(array, 2), size(array, 3), size(array, 4), &
            size(array, 5), size(array, 6), size(array, 7)):: ft

       !--- allocate 1-d array to be used by routine fftn
       complex(fftkind), dimension(:), allocatable :: work

       allocate( work(size(array)) )
       work = reshape(array, (/ size(array) /))
       call fftn(work, shape(array), dim, inv, stat = StatusF)
       ft = reshape(work, (/ size(array, 1), size(array, 2), size(array, 3), &
                             size(array, 4), size(array, 5), size(array, 6), &
                             size(array, 7) /))

       if (allocated(work)) deallocate(work)

       return
    End Function Fft7D

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Fftn(Array, Shape, Dim, Inv, Stat)
    !!--++    complex(fftkind), dimension(:), intent(in out)       :: array
    !!--++    integer,          dimension(:), intent(in)           :: shape
    !!--++    integer,          dimension(:), intent(in), optional :: dim
    !!--++    logical,                        intent(in), optional :: inv
    !!--++    integer,                        intent(in), optional :: stat
    !!--++
    !!--++    (PRIVATE)
    !!--++    General routine for FFT calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Fftn(Array, Shape, Dim, Inv, Stat)
       !--- formal parameters
       complex(fftkind), dimension(:), intent(in out)       :: array
       integer,          dimension(:), intent(in)           :: shape
       integer,          dimension(:), intent(in),  optional:: dim
       logical,                        intent(in),  optional:: inv
       integer,                        intent(out), optional:: stat

       !--- local arrays
       integer, dimension(size(shape)):: d
       !--- local scalars
       logical      :: inverse
       integer      :: i, ndim, ntotal
       real(fftkind):: scal

       !--- optional parameter settings
       if (present(inv)) then
          inverse = inv
       else
          inverse = .false.
       end if
       if (present(dim)) then
          ndim = min(size(dim), size(d))
          d(1:ndim) = dim(1:ndim)
       else
          ndim = size(d)
          d = (/(i, i = 1, size(d))/)
       end if
       ntotal = product(shape)
       scal = sqrt(1.0_fftkind / product(shape(d(1:ndim))))
       array(1:ntotal) = array(1:ntotal) * scal
       do i = 1, ndim
          call fftradix(array, ntotal, shape(d(i)), product(shape(1:d(i))), &
               inverse, stat)
          if (present(stat)) then
             if (stat /=0) return
          end if
       end do

       return
    End Subroutine Fftn

    !!--++
    !!--++ Subroutine Fftradix(Array, Ntotal, Npass, Nspan, Inv, Stat)
    !!--++    complex(fftkind), dimension(:), intent(in out)       :: array
    !!--++    integer,                        intent(in)           :: ntotal
    !!--++    integer,                        intent(in)           :: npass
    !!--++    integer,                        intent(in)           :: nspan
    !!--++    logical,                        intent(in)           :: inv
    !!--++    integer,                        intent(in), optional :: stat
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Fftradix(Array, Ntotal, Npass, Nspan, Inv, Stat)
       !--- formal parameters
       INTEGER,                        INTENT(IN)           :: ntotal, npass, nspan
       COMPLEX(fftkind), DIMENSION(:), INTENT(IN OUT)       :: array
       LOGICAL,                        INTENT(IN)           :: inv
       INTEGER,                        INTENT(OUT), OPTIONAL:: stat

       !--- local arrays
       INTEGER,          DIMENSION(BIT_SIZE(0))     :: factor
       COMPLEX(fftkind), DIMENSION(:), ALLOCATABLE  :: ctmp
       REAL(fftkind),    DIMENSION(:), ALLOCATABLE  :: sine, cosine
       INTEGER,          DIMENSION(:), ALLOCATABLE  :: perm

       !--- local scalars
       INTEGER         :: ii, kspan, ispan
       INTEGER         :: j, jc, jf, jj, k, k1, k2, k3, k4, kk, kt, nn, ns, nt
       INTEGER         :: maxfactor, nfactor, nperm
       REAL(fftkind)   :: s60, c72, s72, pi2
       REAL(fftkind)   :: radf
       REAL(fftkind)   :: c1, c2, c3, cd, ak
       REAL(fftkind)   :: s1, s2, s3, sd
       COMPLEX(fftkind):: cc, cj, ck, cjp, cjm, ckp, ckm

       IF (npass <= 1) RETURN
       c72 = cos72
       IF (inv) THEN
          s72 = sin72
          s60 = sin60
          pi2 = pi
       ELSE
          s72 = -sin72
          s60 = -sin60
          pi2 = -pi
       END IF
       nt = ntotal
       ns = nspan
       kspan = ns
       nn = nt - 1
       jc = ns / npass
       radf = pi2 * jc
       pi2 = pi2 * 2.0_fftkind !-- use 2 PI from here on
       CALL factorize()
       maxfactor = MAXVAL(factor (:nfactor))
       IF (nfactor - ISHFT(kt, 1) > 0) THEN
          nperm = MAX(nfactor + 1, PRODUCT(factor(kt+1: nfactor-kt)) - 1)
       ELSE
          nperm = nfactor + 1
       END IF
       IF (PRESENT(stat)) THEN
          ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor), STAT=stat)
          IF (stat /= 0) RETURN
          CALL transform()
          DEALLOCATE(sine, cosine, STAT=stat)
          IF (stat /= 0) RETURN
          ALLOCATE(perm(nperm), STAT=stat)
          IF (stat /= 0) RETURN
          CALL permute()
          DEALLOCATE(perm, ctmp, STAT=stat)
          IF (stat /= 0) RETURN
       ELSE
          ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor))
          CALL transform()
          DEALLOCATE(sine, cosine)
          ALLOCATE(perm(nperm))
          CALL permute()
          DEALLOCATE(perm, ctmp)
       END IF
       RETURN

    Contains

       !!--++
       !!--++ Subroutine factorize()
       !!--++
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE factorize()

          nfactor = 0
          k = npass
          DO WHILE (MOD(k, 16) == 0)
             nfactor = nfactor + 1
             factor (nfactor) = 4
             k = k / 16
          END DO

          j = 3
          jj = 9
          DO
             DO WHILE (MOD(k, jj) == 0)
                nfactor=nfactor + 1
                factor (nfactor) = j
                k = k / jj
             END DO
             j = j + 2
             jj = j * j
             IF (jj > k) EXIT
          END DO

          IF (k <= 4) THEN
             kt = nfactor
             factor (nfactor + 1) = k
             IF (k /= 1) nfactor = nfactor + 1
          ELSE
             IF (k - ISHFT(k / 4, 2) == 0) THEN
                nfactor = nfactor + 1
                factor (nfactor) = 2
                k = k / 4
             END IF
             kt = nfactor
             j = 2
             DO
                IF (MOD(k, j) == 0) THEN
                   nfactor = nfactor + 1
                   factor (nfactor) = j
                   k = k / j
                END IF
                j = ISHFT((j + 1)/2, 1) + 1
                IF (j > k) EXIT
             END DO
          END IF

          IF (kt > 0) THEN
             j = kt
             DO
                nfactor = nfactor + 1
                factor (nfactor) = factor (j)
                j = j - 1
                IF (j==0) EXIT
             END DO
          END IF

          RETURN
       END SUBROUTINE factorize

       !!--++
       !!--++ Subroutine transform()
       !!--++
       !!--++    compute fourier transform
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE transform()

          ii = 0
          jf = 0
          DO
             sd = radf / kspan
             cd = SIN(sd)
             cd = 2.0_fftkind * cd * cd
             sd = SIN(sd + sd)
             kk = 1
             ii = ii + 1
             SELECT CASE (factor (ii))
                CASE (2)
                   !-- transform for factor of 2 (including rotation factor)
                   kspan = kspan / 2
                   k1 = kspan + 2
                   DO
                      DO
                         k2 = kk + kspan
                         ck = array(k2)
                         array(k2) = array(kk) - ck
                         array(kk) = array(kk) + ck
                         kk = k2 + kspan
                         IF (kk > nn) EXIT
                      END DO
                      kk = kk - nn
                      IF (kk > jc) EXIT
                   END DO

                   IF (kk > kspan) RETURN
                   DO
                      c1 = 1.0_fftkind - cd
                      s1 = sd
                      DO
                         DO
                            DO
                               k2 = kk + kspan
                               ck = array(kk) - array(k2)
                               array(kk) = array(kk) + array(k2)
                               array(k2) = ck * CMPLX(c1, s1, kind=fftkind)
                               kk = k2 + kspan
                               IF (kk >= nt) EXIT
                            END DO
                            k2 = kk - nt
                            c1 = -c1
                            kk = k1 - k2
                            IF (kk <= k2) EXIT
                         END DO
                         ak = c1 - (cd * c1 + sd * s1)
                         s1 = sd * c1 - cd * s1 + s1
                         c1 = 2.0_fftkind - (ak * ak + s1 * s1)
                         s1 = s1 * c1
                         c1 = c1 * ak
                         kk = kk + jc
                         IF (kk >= k2) EXIT
                      END DO

                      k1 = k1 + 1 + 1
                      kk = (k1 - kspan) / 2 + jc
                      IF (kk > jc + jc) EXIT
                   END DO

                CASE (4) !-- transform for factor of 4
                   ispan = kspan
                   kspan = kspan / 4
                   DO
                      c1 = 1.0_fftkind
                      s1 = 0.0_fftkind
                      DO
                         DO
                            k1 = kk + kspan
                            k2 = k1 + kspan
                            k3 = k2 + kspan
                            ckp = array(kk) + array(k2)
                            ckm = array(kk) - array(k2)
                            cjp = array(k1) + array(k3)
                            cjm = array(k1) - array(k3)
                            array(kk) = ckp + cjp
                            cjp = ckp - cjp
                            IF (inv) THEN
                               ckp = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), kind=fftkind)
                               ckm = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), kind=fftkind)
                            ELSE
                               ckp = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), kind=fftkind)
                               ckm = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), kind=fftkind)
                            END IF
                            !-- avoid useless multiplies
                            IF (s1 == 0.0_fftkind) THEN
                               array(k1) = ckp
                               array(k2) = cjp
                               array(k3) = ckm
                            ELSE
                               array(k1) = ckp * CMPLX(c1, s1, kind=fftkind)
                               array(k2) = cjp * CMPLX(c2, s2, kind=fftkind)
                               array(k3) = ckm * CMPLX(c3, s3, kind=fftkind)
                            END IF
                            kk = k3 + kspan
                            IF (kk > nt) EXIT
                         END DO
                         c2 = c1 - (cd * c1 + sd * s1)
                         s1 = sd * c1 - cd * s1 + s1
                         c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                         s1 = s1 * c1
                         c1 = c1 * c2
                         !-- values of c2, c3, s2, s3 that will get used next time
                         c2 = c1 * c1 - s1 * s1
                         s2 = 2.0_fftkind * c1 * s1
                         c3 = c2 * c1 - s2 * s1
                         s3 = c2 * s1 + s2 * c1
                         kk = kk - nt + jc
                         IF (kk > kspan) EXIT
                      END DO
                      kk = kk - kspan + 1
                      IF (kk > jc) EXIT
                   END DO
                   IF (kspan == jc) RETURN

                CASE default
                   !-- transform for odd factors
                   k = factor (ii)
                   ispan = kspan
                   kspan = kspan / k
                   SELECT CASE (k)
                      CASE (3) !-- transform for factor of 3 (optional code)
                         DO
                            DO
                               k1 = kk + kspan
                               k2 = k1 + kspan
                               ck = array(kk)
                               cj = array(k1) + array(k2)
                               array(kk) = ck + cj
                               ck = ck - 0.5_fftkind * cj
                               cj = (array(k1) - array(k2)) * s60
                               array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               kk = k2 + kspan
                               IF (kk >= nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO

                      CASE (5) !-- transform for factor of 5 (optional code)
                         c2 = c72 * c72 - s72 * s72
                         s2 = 2.0_fftkind * c72 * s72
                         DO
                            DO
                               k1 = kk + kspan
                               k2 = k1 + kspan
                               k3 = k2 + kspan
                               k4 = k3 + kspan
                               ckp = array(k1) + array(k4)
                               ckm = array(k1) - array(k4)
                               cjp = array(k2) + array(k3)
                               cjm = array(k2) - array(k3)
                               cc = array(kk)
                               array(kk) = cc + ckp + cjp
                               ck = ckp * c72 + cjp * c2 + cc
                               cj = ckm * s72 + cjm * s2
                               array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k4) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               ck = ckp * c2 + cjp * c72 + cc
                               cj = ckm * s2 - cjm * s72
                               array(k2) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                               array(k3) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                               kk = k4 + kspan
                               IF (kk >= nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO

                      CASE default
                         IF (k /= jf) THEN
                            jf = k
                            s1 = pi2 / k
                            c1 = COS(s1)
                            s1 = SIN(s1)
                            cosine (jf) = 1.0_fftkind
                            sine (jf) = 0.0_fftkind
                            j = 1
                            DO
                               cosine (j) = cosine (k) * c1 + sine (k) * s1
                               sine (j) = cosine (k) * s1 - sine (k) * c1
                               k = k-1
                               cosine (k) = cosine (j)
                               sine (k) = -sine (j)
                               j = j + 1
                               IF (j >= k) EXIT
                            END DO
                         END IF

                         DO
                            DO
                               k1 = kk
                               k2 = kk + ispan
                               cc = array(kk)
                               ck = cc
                               j = 1
                               k1 = k1 + kspan
                               DO
                                  k2 = k2 - kspan
                                  j = j + 1
                                  ctmp(j) = array(k1) + array(k2)
                                  ck = ck + ctmp(j)
                                  j = j + 1
                                  ctmp(j) = array(k1) - array(k2)
                                  k1 = k1 + kspan
                                  IF (k1 >= k2) EXIT
                               END DO
                               array(kk) = ck
                               k1 = kk
                               k2 = kk + ispan
                               j = 1
                               DO
                                  k1 = k1 + kspan
                                  k2 = k2 - kspan
                                  jj = j
                                  ck = cc
                                  cj = (0.0_fftkind, 0.0_fftkind)
                                  k = 1
                                  DO
                                     k = k + 1
                                     ck = ck + ctmp(k) * cosine (jj)
                                     k = k + 1
                                     cj = cj + ctmp(k) * sine (jj)
                                     jj = jj + j
                                     IF (jj > jf) jj = jj - jf
                                     IF (k >= jf) EXIT
                                  END DO
                                  k = jf - j
                                  array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), kind=fftkind)
                                  array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), kind=fftkind)
                                  j = j + 1
                                  IF (j >= k) EXIT
                               END DO
                               kk = kk + ispan
                               IF (kk > nn) EXIT
                            END DO
                            kk = kk - nn
                            IF (kk > kspan) EXIT
                         END DO
                   END SELECT !k

                   !--  multiply by rotation factor (except for factors of 2 and 4)
                   IF (ii == nfactor) RETURN
                   kk = jc + 1
                   DO
                      c2 = 1.0_fftkind - cd
                      s1 = sd
                      DO
                         c1 = c2
                         s2 = s1
                         kk = kk + kspan
                         DO
                            DO
                               array(kk) = CMPLX(c2, s2, kind=fftkind) * array(kk)
                               kk = kk + ispan
                               IF (kk > nt) EXIT
                            END DO
                            ak = s1 * s2
                            s2 = s1 * c2 + c1 * s2
                            c2 = c1 * c2 - ak
                            kk = kk - nt + kspan
                            IF (kk > ispan) EXIT
                         END DO
                         c2 = c1 - (cd * c1 + sd * s1)
                         s1 = s1 + sd * c1 - cd * s1
                         c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                         s1 = s1 * c1
                         c2 = c2 * c1
                         kk = kk - ispan + jc
                         IF (kk > kspan) EXIT
                      END DO
                      kk = kk - kspan + jc + 1
                      IF (kk > jc + jc) EXIT
                   END DO
             END SELECT ! Factor
          END DO

          RETURN
       END SUBROUTINE transform

       !!--++
       !!--++ Subroutine permute()
       !!--++
       !!--++    permute the results to normal order---done in two stages
       !!--++    permutation for square factors of n
       !!--++
       !!--++    (PRIVATE)
       !!--++
       !!--++ Update: February - 2005
       !!
       SUBROUTINE permute()

          perm (1) = ns
          IF (kt > 0) THEN
             k = kt + kt + 1
             IF (nfactor < k) k = k - 1
             j = 1
             perm (k + 1) = jc
             DO
                perm (j + 1) = perm (j) / factor (j)
                perm (k) = perm (k + 1) * factor (j)
                j = j + 1
                k = k - 1
                IF (j >= k) EXIT
             END DO
             k3 = perm (k + 1)
             kspan = perm (2)
             kk = jc + 1
             k2 = kspan + 1
             j = 1
             IF (npass /= ntotal) THEN
                permute_multi: DO
                   DO
                      DO
                         k = kk + jc
                         DO
                            !-- swap array(kk) <> array(k2)
                            ck = array(kk)
                            array(kk) = array(k2)
                            array(k2) = ck
                            kk = kk + 1
                            k2 = k2 + 1
                            IF (kk >= k) EXIT
                         END DO
                         kk = kk + ns - jc
                         k2 = k2 + ns - jc
                         IF (kk >= nt) EXIT
                      END DO
                      kk = kk - nt + jc
                      k2 = k2 - nt + kspan
                      IF (k2 >= ns) EXIT
                   END DO

                   DO
                      DO
                         k2 = k2 - perm (j)
                         j = j + 1
                         k2 = perm (j + 1) + k2
                         IF (k2 <= perm (j)) EXIT
                      END DO
                      j = 1

                      DO
                         IF (kk < k2) CYCLE permute_multi
                         kk = kk + jc
                         k2 = k2 + kspan
                         IF (k2 >= ns) EXIT
                      END DO
                      IF (kk >= ns) EXIT
                   END DO

                   EXIT
                END DO permute_multi
             ELSE
                permute_single: DO
                   DO
                      !-- swap array(kk) <> array(k2)
                      ck = array(kk)
                      array(kk) = array(k2)
                      array(k2) = ck
                      kk = kk + 1
                      k2 = k2 + kspan
                      IF (k2 >= ns) EXIT
                   END DO

                   DO
                      DO
                         k2 = k2 - perm (j)
                         j = j + 1
                         k2 = perm (j + 1) + k2
                         IF (k2 <= perm (j)) EXIT
                      END DO
                      j = 1
                      DO
                         IF (kk < k2) CYCLE permute_single
                         kk = kk + 1
                         k2 = k2 + kspan
                         IF (k2 >= ns) EXIT
                      END DO
                      IF (kk >= ns) EXIT
                   END DO
                   EXIT
                END DO permute_single
             END IF
             jc = k3
          END IF

          IF (ISHFT(kt, 1) + 1 >= nfactor) RETURN
          ispan = perm (kt + 1)

          !-- permutation for square-free factors of n
          j = nfactor - kt
          factor (j + 1) = 1
          DO
             factor(j) = factor(j) * factor(j+1)
             j = j - 1
             IF (j == kt) EXIT
          END DO
          kt = kt + 1
          nn = factor(kt) - 1
          j = 0
          jj = 0
          DO
             k = kt + 1
             k2 = factor(kt)
             kk = factor(k)
             j = j + 1
             IF (j > nn) EXIT !-- exit infinite loop
             jj = jj + kk
             DO WHILE (jj >= k2)
                jj = jj - k2
                k2 = kk
                k = k + 1
                kk = factor(k)
                jj = jj + kk
             END DO
             perm (j) = jj
          END DO

          !--  determine the permutation cycles of length greater than 1
          j = 0
          DO
             DO
                j = j + 1
                kk = perm(j)
                IF (kk >= 0) EXIT
             END DO
             IF (kk /= j) THEN
                DO
                   k = kk
                   kk = perm (k)
                   perm (k) = -kk
                   IF (kk == j) EXIT
                END DO
                k3 = kk
             ELSE
                perm (j) = -j
                IF (j == nn) EXIT !-- exit infinite loop
             END IF
          END DO

          !--  reorder a and b, following the permutation cycles
          DO
             j = k3 + 1
             nt = nt - ispan
             ii = nt - 1 + 1
             IF (nt < 0) EXIT !-- exit infinite loop
             DO
                DO
                   j = j-1
                   IF (perm(j) >= 0) EXIT
                END DO
                jj = jc

                DO
                   kspan = jj
                   IF (jj > maxfactor) kspan = maxfactor
                   jj = jj - kspan
                   k = perm(j)
                   kk = jc * k + ii + jj
                   k1 = kk + kspan
                   k2 = 0
                   DO
                      k2 = k2 + 1
                      ctmp(k2) = array(k1)
                      k1 = k1 - 1
                      IF (k1 == kk) EXIT
                   END DO
                   DO
                      k1 = kk + kspan
                      k2 = k1 - jc * (k + perm(k))
                      k = -perm(k)
                      DO
                         array(k1) = array(k2)
                         k1 = k1 - 1
                         k2 = k2 - 1
                         IF (k1 == kk) EXIT
                      END DO
                      kk = k2
                      IF (k == j) EXIT
                   END DO
                   k1 = kk + kspan
                   k2 = 0
                   DO
                      k2 = k2 + 1
                      array(k1) = ctmp(k2)
                      k1 = k1 - 1
                      IF (k1 == kk) EXIT
                   END DO
                   IF (jj == 0) EXIT
                END DO
                IF (j == 1) EXIT
             END DO
          END DO

          RETURN
       END SUBROUTINE permute

    End Subroutine Fftradix

 End Module Fft_Gen
