  Module J_k_exchange
    Implicit none

    public :: j_k, spin_conf, genj, genk  !, read_exchange
    private :: sin2p, cos2p, spin_confr, spin_confc

    integer, public, parameter :: nat=25, id=50, njo=180000, ijo=15, nv=125000
    integer, public                           :: natcel !number of atoms/primitive cell
    real,    public, dimension        (3,nat) :: xyz    !Coordinates of atoms in the cell
    real,    public, dimension        (  nat) :: s      !Spin of atoms
    real,    public, dimension        (  id ) :: valj   !Exchange interactions of bond id
    integer, public, dimension      (nat,nat) :: nterm  !number of terms of J(k)(i,j)
    integer, public, dimension   (nat,nat,id) :: nvalj  !Pointer to the value of Jid
    real   , public, dimension (3,nat,nat,id) :: trans  !Translation associated to the bond "id"
                                                        !in the position term (i,j) of J(k)
    character (len=6), public, dimension(nat) :: rnam   !Name of the atoms


    Interface spin_conf
      Module procedure spin_confc
      Module procedure spin_confr
    End Interface

    CONTAINS

!----------------------------------------------------------------------------------------------
!      subroutine read_exchange()

!      end subroutine read_exchange
!----------------------------------------------------------------------------------------------
      SUBROUTINE j_k(rk,exchk,iex)
      real,    INTENT(IN)     , dimension(3)       :: rk    !propagation vector
      COMPLEX, INTENT(OUT)    , dimension(nat,nat) :: exchk !J(k)
      integer, INTENT(OUT)                         :: iex   !=0 if real J(k)
      integer :: i,j,k,n
      real :: arg

      iex=0             !Imaginary part equal to zero
      DO i=1,natcel
        DO j=1,natcel
          exchk(i,j)=cmplx(0.0,0.0)
          DO n=1,nterm(i,j)
            arg=0.0
            DO k=1,3
              arg=arg+rk(k)*trans(k,i,j,n)
            END DO
            exchk(i,j)=exchk(i,j)-valj(nvalj(i,j,n))*cmplx( cos2p(-arg),sin2p(-arg)) ! -J(k) to minimize!!!
          END DO
          exchk(i,j)=exchk(i,j)*s(i)*s(j)       !Spin values absorbed in the Energy
          IF(ABS(aimag(exchk(i,j))) > 1.e-04) iex=1
        END DO
      END DO
      RETURN
      END SUBROUTINE j_k
!-----------------------------------------------------------------------
        FUNCTION cos2p(x) result(f_val)
         real, intent (in) :: x
         real :: f_val
         f_val=COS(6.283185307*x)
         RETURN
        END FUNCTION cos2p
        FUNCTION sin2p(x)  result(f_val)
         real, intent (in) :: x
         real :: f_val
         f_val=SIN(6.283185307*x)
         RETURN
        END FUNCTION sin2p
!-----------------------------------------------------------------------
      SUBROUTINE spin_confc(lun,rk,energy,d)
      integer, INTENT(IN)               :: lun
      real, INTENT(IN),dimension(3)     :: rk
      real, INTENT(IN)                  :: energy
      complex, INTENT(IN),dimension(nat):: d
      real, dimension(nat)              :: rmo, phas
      integer :: i,j
      write(unit=lun,fmt="(/,a,f12.4,/,a,3f8.4,/,a,/)")  &
          " => Spin configurations of Energy: ",energy,  &
          "    And propagation vector       : ",(rk(i),i=1,3),  &
          "    can be obtained from the following expressions:"

        DO i=1,natcel
          rmo(i)=abs(d(i))
          IF(rmo(i) > 1.e-15) THEN
            phas(i)=ATAN2(aimag(d(i)),real(d(i)))/6.2831853
          ELSE
            phas(i)=0.0
          END IF
        END DO
      DO i=1,natcel
        write(unit=lun,fmt="(i4,a,a6,4f9.5)")i,"  ",rnam(i),(xyz(j,i),j=1,3),s(i)
        IF(rmo(i) > 1.e-15) THEN
          if(phas(i) < 0.0 ) then
           write(unit=lun,fmt="(a,f7.4,a)")"               Mx =  Cx . cos{2pi(k*Rn",phas(i),")}"
           write(unit=lun,fmt="(a,f7.4,a)")"               My =  Cy . sin{2pi(k*Rn",phas(i),")}"
          else
           write(unit=lun,fmt="(a,f7.4,a)")"               Mx =  Cx . cos{2pi(k*Rn +",phas(i),")}"
           write(unit=lun,fmt="(a,f7.4,a)")"               My =  Cy . sin{2pi(k*Rn +",phas(i),")}"
          end if
           write(unit=lun,fmt="(a)")       "               Mz =  0"
        ELSE
           write(unit=lun,fmt="(a)")       "               M =  0    <--- Atom with zero moment"
        END IF
      END DO
      RETURN
      END SUBROUTINE spin_confc

      SUBROUTINE spin_confr(lun,rk,energy,d)
      integer, INTENT(IN)               :: lun
      real, INTENT(IN),dimension(3)     :: rk
      real, INTENT(IN)                  :: energy
      real, INTENT(IN),dimension(nat):: d
      real, dimension(nat)              :: rmo, phas
      integer :: i,j
      write(unit=lun,fmt="(/,a,f12.4,/,a,3f8.4,/,a,/)")  &
          " => Spin configurations of Energy: ",energy,  &
          "    And propagation vector       : ",(rk(i),i=1,3),  &
          "    can be obtained from the following expressions:"

        DO i=1,natcel
          IF(real(d(i)) >= 0.0) THEN
            phas(i)=0.0
            rmo(i)=d(i)
          ELSE
            phas(i)=0.5
            rmo(i)=-d(i)
          END IF
        END DO
      DO i=1,natcel
        write(unit=lun,fmt="(i4,a,a6,4f9.5)")i,"  ",rnam(i),(xyz(j,i),j=1,3),s(i)
        IF(rmo(i) > 1.e-15) THEN
          if(phas(i) < 0.0 ) then
           write(unit=lun,fmt="(a,f7.4,a)")"               Mx =  Cx . cos{2pi(k*Rn",phas(i),")}"
           write(unit=lun,fmt="(a,f7.4,a)")"               My =  Cy . sin{2pi(k*Rn",phas(i),")}"
          else
           write(unit=lun,fmt="(a,f7.4,a)")"               Mx =  Cx . cos{2pi(k*Rn +",phas(i),")}"
           write(unit=lun,fmt="(a,f7.4,a)")"               My =  Cy . sin{2pi(k*Rn +",phas(i),")}"
          end if
           write(unit=lun,fmt="(a)")       "               Mz =  0"
        ELSE
           write(unit=lun,fmt="(a)")       "               M =  0    <--- Atom with zero moment"
        END IF
      END DO
      RETURN
      END SUBROUTINE spin_confr
!-----------------------------------------------------------------------

      SUBROUTINE genj(nojvar,njotas,rang1,rang2,npoi,valjj)

      integer, INTENT(IN)                     :: nojvar
      integer, INTENT(OUT)                    :: njotas
      real,    INTENT(IN) , dimension(ijo)    :: rang1
      real,    INTENT(IN) , dimension(ijo)    :: rang2
      integer, INTENT(IN) , dimension(ijo)    :: npoi
      real,    INTENT(OUT), dimension(ijo,njo):: valjj

      integer, dimension(ijo) ::  nj
      real,    dimension(ijo) ::  step

      integer :: i,i1,i2,i3,i4,i5, jj

      DO i=1,nojvar
        if(npoi(i) == 1) then
          step(i) = 0.0
        else
          step(i)=(rang2(i)-rang1(i))/real(npoi(i)-1)
        end if
      END DO
      valjj(1:nojvar,1:njo)=0.0
      jj=0
      DO i1=1,npoi(1)
        nj(1)=i1
        DO i2=1,npoi(2)
          nj(2)=i2
          DO i3=1,npoi(3)
            nj(3)=i3
            DO i4=1,npoi(4)
              nj(4)=i4
              DO i5=1,npoi(5)
                nj(5)=i5          !Increase the number of loops if ijo>5
                jj=jj+1
                DO i=1,nojvar
                  valjj(i,jj)=rang1(i)+real(nj(i)-1)*step(i)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      njotas=jj
      RETURN
      END SUBROUTINE genj
!-----------------------------------------------------------------------

      SUBROUTINE genk(lun,ln,infil,nvk,vk,iprint)

      integer, intent(in)                      :: lun, ln
      character (len=256), intent(in out)      :: infil
      integer, intent(out)                     :: nvk
      real, intent(out), dimension(3,nv)       :: vk
      integer, intent(out)                     :: iprint
      integer :: ier,i,j,k
      character (len=1)    :: ans

      logical :: double_k
      real, dimension(3)   :: xran
      real, dimension(3,84), parameter :: sk = reshape ( (/         &
          0.00000,0.00000,0.00000,  0.50000,0.00000,0.00000,  0.00000,0.50000,0.00000,  &
          0.00000,0.00000,0.50000,  0.50000,0.50000,0.00000,  0.50000,0.00000,0.50000,  &
          0.00000,0.50000,0.50000,  0.50000,0.50000,0.50000,  0.25000,0.00000,0.00000,  &
          0.00000,0.25000,0.00000,  0.00000,0.00000,0.25000,  0.25000,0.25000,0.00000,  &
          0.25000,0.00000,0.25000,  0.00000,0.25000,0.25000,  0.25000,0.25000,0.25000,  &
          0.50000,0.25000,0.25000,  0.25000,0.50000,0.25000,  0.25000,0.25000,0.50000,  &
          0.50000,0.50000,0.25000,  0.25000,0.50000,0.50000,  0.50000,0.25000,0.50000,  &
          0.50000,0.00000,0.25000,  0.00000,0.50000,0.25000,  0.50000,0.25000,0.00000,  &
          0.00000,0.25000,0.50000,  0.25000,0.00000,0.50000,  0.25000,0.50000,0.00000,  &
          0.33333,0.00000,0.00000,  0.00000,0.33333,0.00000,  0.00000,0.00000,0.33333,  &
          0.33333,0.33333,0.00000,  0.33333,0.00000,0.33333,  0.00000,0.33333,0.33333,  &
          0.33333,0.33333,0.33333,  0.33333,0.66667,0.00000,  0.66667,0.33333,0.00000,  &
          0.66667,0.66667,0.00000,  0.66667,0.00000,0.66667,  0.00000,0.66667,0.66667,  &
          0.66667,0.33333,0.66667,  0.66667,0.66667,0.33333,  0.33333,0.66667,0.66667,  &
          0.12500,0.00000,0.00000,  0.00000,0.12500,0.00000,  0.00000,0.00000,0.12500,  &
          0.12500,0.12500,0.00000,  0.12500,0.00000,0.12500,  0.00000,0.12500,0.12500,  &
          0.12500,0.12500,0.12500,  0.12500,0.25000,0.00000,  0.12500,0.00000,0.25000,  &
          0.00000,0.12500,0.25000,  0.25000,0.12500,0.00000,  0.25000,0.00000,0.12500,  &
          0.00000,0.25000,0.12500,  0.25000,0.25000,0.12500,  0.25000,0.12500,0.25000,  &
          0.12500,0.25000,0.25000,  0.12500,0.25000,0.12500,  0.12500,0.12500,0.25000,  &
          0.00000,0.50000,0.12500,  0.50000,0.00000,0.12500,  0.50000,0.12500,0.00000,  &
          0.00000,0.12500,0.50000,  0.12500,0.00000,0.50000,  0.12500,0.50000,0.00000,  &
          0.12500,0.50000,0.50000,  0.50000,0.12500,0.50000,  0.50000,0.50000,0.12500,  &
          0.50000,0.12500,0.12500,  0.12500,0.50000,0.12500,  0.12500,0.12500,0.50000,  &
          0.33333,0.50000,0.00000,  0.33333,0.00000,0.50000,  0.00000,0.33333,0.50000,  &
          0.50000,0.33333,0.00000,  0.50000,0.00000,0.33333,  0.00000,0.50000,0.33333,  &
          0.50000,0.33333,0.50000,  0.50000,0.50000,0.33333,  0.33333,0.50000,0.50000,  &
          0.33333,0.50000,0.33333,  0.33333,0.33333,0.50000,  0.50000,0.33333,0.33333   &
                     /),(/3,84/) )


      real,    dimension(3) ::  x1,x2,x3,fs,dk
      integer, dimension(3) ::  ngrid
      integer :: iop, n, l, i1,i2
      iprint=0
      double_k=.false.
      call random_seed()
      do
       write(unit=*,fmt="(/,a,/)")"    Vectors k in BZ, components in [-1,1]"
       write(unit=*,fmt="(a)")" => Options:"
       write(unit=*,fmt="(a)")"             0: Special k-vectors"
       write(unit=*,fmt="(a)")"             1: Give individual k-vectors"
       write(unit=*,fmt="(a)")"             2: k-vectors along a line"
       write(unit=*,fmt="(a)")"             3: k-vectors in a plane"
       write(unit=*,fmt="(a)")"             4: General + Special k-vectors"
       write(unit=*,fmt="(a)")"             5: Grid inside the Brillouin Zone"
       write(unit=*,fmt="(a)")"   "
       write(unit=*,fmt="(a)",advance="no")" => Give the option: "
       read(unit=*,fmt=*,iostat=ier) iop
        if(ier /=0) cycle
        if(iop < 0 .or. iop > 5) cycle
        exit
      end do
      Select Case (iop)

      Case (0)
        iprint=1
        nvk=84
        DO i=1,84
          DO j=1,3
            vk(j,i)=sk(j,i)
          END DO
        END DO
        write(unit=lun,fmt="(/,a)") " => Only special k-vectors stored in the program will be tested"
        write(unit=lun,fmt="( a)") "    These vectors are:"
         do i=1,84,4
           write(unit=lun,fmt="(4(a,i2,a,3f8.4))") ("  #",i+l,":",vk(:,i+l),l=0,3)
         end do

      Case (1)
        iprint=1
        write(unit=*,fmt="(a,a,a)")  &
            " => Interactive (i) or read from ",infil(1:ln),".kve file (k): "
        read(unit=*,fmt="(a)") ans
        IF(ans == "k" .OR. ans == "K") THEN
          open(unit=20,file=infil(1:ln)//".kve",status="replace",action="write")
          read(unit=20,fmt=*) nvk
          DO j=1,nvk
            read(unit=20,fmt=*) (vk(i,j),i=1,3)
          END DO
          close(unit=20)
        ELSE
         do
          write(unit=*,fmt="(a,i6,a)",advance="no")" => Number of k-vectors (<",nv,"): "
          read(unit=*,fmt=*,iostat=ier) nvk
          if(ier /=0) cycle
          exit
         end do
          DO j=1,nvk
           do
            write(unit=*,fmt="(a,i2,a)",advance="no")" => K-vector No ",j,": "
            read(unit=*,fmt=*,iostat=ier) (vk(i,j),i=1,3)
            if(ier /=0) cycle
            exit
           end do
          END DO
        END IF

        write(unit=lun,fmt="(/,a)",advance="no") " => Selected set of k-vectors given interactively/in a file "
        write(unit=lun,fmt="(  a)") "    These vectors are:"
         do i=1,nvk,4
           write(unit=lun,fmt="(4(a,i2,a,3f8.4))") ("  #",i+l,":",vk(:,i+l),l=0,3)
         end do

      Case (2)
        iprint=1
         do
          write(unit=*,fmt="(a,i6,a)",advance="no")" => Number of k-vectors (<",nv,"): "
          read(unit=*,fmt=*,iostat=ier) nvk
          if(ier /=0) cycle
          exit
         end do
        write(unit=*,fmt="(a)")" => Origin and extreme of the line (in r.l.u.): "
        read(unit=*,fmt=*) (x1(j),j=1,3),(x2(k),k=1,3)
        DO i=1,3
          x3(i)=(x2(i)-x1(i))/real(nvk-1)
        END DO
        DO j=1,nvk
          DO i=1,3
            vk(i,j)=x1(i)+real(j-1)*x3(i)
          END DO
        END DO

        write(unit=lun,fmt="(/,a)")      " => Test of k-vectors along a line "
        write(unit=lun,fmt="( a,3f8.4)") "    These vectors are between:",(x1(j),j=1,3)
        write(unit=lun,fmt="( a,3f8.4)") "                         and :",(x2(j),j=1,3)
        write(unit=lun,fmt="( a,i6   )") "    Total number of k-vectors:",nvk

      Case (3)
        iprint=1
         do
          write(unit=*,fmt="(a,i6,a)",advance="no")" => Number of k-vectors (<",nv,"): "
          read(unit=*,fmt=*,iostat=ier) nvk
          if(ier /=0) cycle
          exit
         end do
        n=SQRT(real(nvk))
        write(unit=*,fmt="(a)")  " => Plane defined by K=m1*U+m2*V  "
        write(unit=*,fmt="(a)")  "    Give the two vectors U,V (in r.l.u.) defining the plane: "
        read(unit=*,fmt=*) (x1(j),j=1,3),(x2(k),k=1,3)
        DO i=1,3
          x1(i)=x1(i)/real(n-1)
          x2(i)=x2(i)/real(n-1)
        END DO
        nvk=0
        DO i1=1,n
          DO i2=1,n
            DO i=1,3
              x3(i)=real(i1-1)*x1(i)+real(i2-1)*x2(i)
              IF(ABS(x3(i)) > 1.0) x3(i)=SIGN(1.0,x3(i))
            END DO
            nvk=nvk+1
            DO i=1,3
              vk(i,nvk)=x3(i)
            END DO
          END DO
        END DO

        write(unit=lun,fmt="(/,a)")      " => Test of k-vectors within a plane defined by K=m1*U+m2*V  "
        write(unit=lun,fmt="( a,3f8.4)") "     The vectors U and V are :",(x1(j),j=1,3)
        write(unit=lun,fmt="( a,3f8.4)") "                         and :",(x2(j),j=1,3)
        write(unit=lun,fmt="( a,i6   )") "    Total number of k-vectors:",nvk
        write(unit=*,fmt="(/,a)")      " => Test of k-vectors within a plane defined by K=m1*U+m2*V  "
        write(unit=*,fmt="( a,3f8.4)") "     The vectors U and V are :",(x1(j),j=1,3)
        write(unit=*,fmt="( a,3f8.4)") "                         and :",(x2(j),j=1,3)
        write(unit=*,fmt="( a,i6   )") "    Total number of k-vectors:",nvk

      Case (4)
         do
          write(unit=*,fmt="(a,i6,a)",advance="no")" => Number of k-vectors (<",nv,"): "
          read(unit=*,fmt=*,iostat=ier) nvk
          if(ier /=0) cycle
          exit
         end do
        write(unit=*,fmt="(a,/,a)",advance="no")" => Give region of k-vectors (in r.l.u.)",  &
                        "              (kx1,kx2,ky1,ky2,kz1,kz2): "
        read(unit=*,fmt=*) (x1(i),x2(i),i=1,3)
        do i=1,3
          if(x1(i) < -0.000001) then
             double_k=.true.
             i1=i
          end if
        end do
        DO i=1,84
          DO j=1,3
            vk(j,i)=sk(j,i)
          END DO
        END DO
        i2=84
        if(double_k) then
          DO i=2,84
            i2=i2+1
            vk(:,i2)=sk(:,i)
            vk(i1,i2)=-vk(i1,i2)
          END DO
        end if

        DO j=i2+1,nvk
            call random_number(xran)
            vk(:,j)=x1(:)+(x2(:)-x1(:))*xran(:)
        END DO

        write(unit=lun,fmt="(/,a)")      " => k-vectors generated at random in the region:   "
        write(unit=lun,fmt="( a,2f8.4)") "    Kx :",x1(1),x2(1)
        write(unit=lun,fmt="( a,2f8.4)") "    Ky :",x1(2),x2(2)
        write(unit=lun,fmt="( a,2f8.4)") "    Kz :",x1(3),x2(3)
        if(double_k) then
          write(unit=lun,fmt="(a)") "    Plus the special 167= 84+83  vectors stored in the program"
        else
          write(unit=lun,fmt="(a)") "    Plus the special 84 vectors stored in the program"
        end if
        write(unit=lun,fmt="( a,i6   )") "    Total number of k-vectors:",nvk

      Case (5)
        do
         write(unit=*,fmt="(a)")" => Region for k-vectors (min,max,num.of divisions):"
         write(unit=*,fmt="(a)",advance="no")"     Give kx_min,kx_max, xdiv:"
         read(unit=*,fmt=*)  x1(1),x2(1),ngrid(1)
         write(unit=*,fmt="(a)",advance="no")"     Give ky_min,ky_max, ydiv:"
         read(unit=*,fmt=*)  x1(2),x2(2),ngrid(2)
         write(unit=*,fmt="(a)",advance="no")"     Give kz_min,kz_max, zdiv:"
         read(unit=*,fmt=*)  x1(3),x2(3),ngrid(3)
         if ((ngrid(1)+1)*(ngrid(2)+1)*(ngrid(3)+1)> nv) then
          write(unit=*,fmt="(a)") " => Too many k-vectors ... reduce the grid!"
          cycle
         else
          exit
         end if
        end do

         dk(:)=(x2(:)-x1(:))/ngrid(:)  !steps in k-space
         nvk=0
         DO i=0,ngrid(1)
           fs(1)=real(i)
           DO j=0,ngrid(2)
             fs(2)=real(j)
            DO k=0,ngrid(3)
              fs(3)=real(k)
               nvk=nvk+1
               vk(:,nvk)=x1(:)+fs(:)*dk(:)       !k-vector components
            END DO
           END DO
         END DO

        write(unit=*,fmt="( a,i6)") " =>  Total number of k_vectors :", nvk
        write(unit=lun,fmt="(/,a)")       " => k-vectors generated in the region:   "
        write(unit=lun,fmt="( a,2f8.4,a,i5,a)") "    Kx :",x1(1),x2(1), " with",ngrid(1)+1," points"
        write(unit=lun,fmt="( a,2f8.4,a,i5,a)") "    Ky :",x1(2),x2(2), " with",ngrid(2)+1," points"
        write(unit=lun,fmt="( a,2f8.4,a,i5,a)") "    Kz :",x1(3),x2(3), " with",ngrid(3)+1," points"
        write(unit=lun,fmt="( a,i6   )") "    Total number of k-vectors:",nvk

      Case default
        nvk=1
        DO i=1,3
          vk(i,1)=0.0
        END DO
      End Select
     ! if(iprint == 1) then
     !   open(unit=11,file=infil(1:ln)//".val",status="replace",action="write")
     ! end if
      RETURN
      END SUBROUTINE genk
!------------------------------------------------------------------------
  end Module J_k_exchange

!----------------------------------------------------------------------------


!-----------------------------------------------------------------------
! Program to calculate the classical magnetic energy for a set of
! k-vectors and exchange parameters. The input file is created by
! the program SIMBO.
! Test version beta-0. Isotropic exchange.
! Created in December 1994 (Juan Rodriguez-Carvajal, LLB)
! New things added in May 1997
! New things added in June 1999
! Transformed to F90 in July 1999
! Some modifications in August 2001
!-----------------------------------------------------------------------
   PROGRAM enermag
      !use Mod_fun
      use Math_gen, only: diagonalize_sh
      use j_k_exchange
      use Super_Exchange
      Implicit none

      real, dimension(3,nv)  :: vk
      integer, dimension(nv) :: nvect
      real, dimension(3) :: ktar, rktar, rk
      CHARACTER(LEN=30), dimension(id ) ::  Jotas
      CHARACTER(LEN=80)   ::  TITLE, cmdline
      CHARACTER(LEN=1)    ::  ANS
      CHARACTER(LEN=20)   ::  spg
      CHARACTER(LEN=132)  ::  forma, line
      CHARACTER (LEN=256) ::  resfil,infil,outfil
      CHARACTER (LEN=4)   ::  jota
      CHARACTER (LEN=30)  ::  expo,name_Jota
      real, dimension(6)  :: ad
      LOGICAL :: ansf, compl, compl_tar
      integer, PARAMETER :: lun1=1, lun2=2, lun7=7
      COMPLEX, dimension(nat,nat) :: exchk,aop,eigen_cvectors
      COMPLEX, dimension(nat)  :: dtarc
      real, dimension(nat,nat) :: aux,eigen_rvectors
      real, dimension(nat)     :: eigen_val,dtarr
      real, dimension(ijo)     :: rang1,rang2
      real, dimension(ijo,njo) :: valjj
      integer, dimension(ijo)  :: ivar,npoi
      integer :: ier,i,j,nex,ln,nt,jnn,inn,lr,njjj,njotas,m,jj,ik,iopt, &
                 iprint,inum,iphase,nojvar,ishifp,nvk,iex,n,nop,noptar, &
                 iwrt, k, narg, iargc
      real    :: dmax,exch,emin,energy,rvalj,energytar,cpu_fin, cpu_ini,cpu_par, &
                 cpu, remaining
!-----------------------------------------------------------------------
      nex=0
      ansf=.false.
      ln=0
      lr=0
      narg=iargc()               !it does'nt work on F for Windows
      if(narg > 0) then
              call getarg(1,infil)
              ln=len_trim(infil)
              outfil=infil
              lr=ln
      end if
      if(narg > 1) then
              call getarg(2,outfil)
              lr=len_trim(outfil)
      end if

      write(unit=*,fmt="(a)")"   *****************************************************"
      write(unit=*,fmt="(a)")"   ****               PROGRAM ENERMAG               ****"
      write(unit=*,fmt="(a)")"   *****************************************************"
      write(unit=*,fmt="(a)")"  "
      write(unit=*,fmt="(a)")"                 *** Version 1.5 *** "
      write(unit=*,fmt="(a)")"   *****************************************************"
      write(unit=*,fmt="(a)")"   * Calculates the Magnetic Energy for k-vectors in BZ*"
      write(unit=*,fmt="(a)")"   *      Uses a Classical Heisenberg Hamiltonian      *"
      write(unit=*,fmt="(a)")"   *        and isotropic exchange interactions        *"
      write(unit=*,fmt="(a)")"   *****************************************************"
      write(unit=*,fmt="(a)") "                   (JRC August-2001, LLB)"
      write(unit=*,fmt="(a)") " "
      write(unit=*,fmt="(a)") " "
      write(unit=*,fmt="(a)") " "
      write(unit=*,fmt="(a)") " "

      if(lr == 0) then
        write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.exc (give xx): "
        read(unit=*,fmt="(a)") infil
        write(unit=*,fmt="(a,a,a)",advance="no") " => Code of the .out file ( <cr>= ",trim(infil),") :"
        read(unit=*,fmt="(a)") outfil
        lr=len_trim(outfil)
        IF(lr == 0) outfil=infil
      end if
!
      open(unit=1,file=trim(infil)//".exc",status="old",action="read",position="rewind",iostat=ier)
      if(ier /= 0) then
        write(unit=*,fmt="(a)") " => File "//trim(infil)//".exc not found !!!"
        STOP
      end if
      open(unit=2,file=trim(outfil)//".out",status="replace",action="write")
      open(unit=7,file=trim(outfil)//".mom",status="replace",action="write")
      write(unit=*,fmt="(a)",advance="no") " => Do you want to save data in *.res file (y/n)?: "
      read(unit=*,fmt="(a)") ans
      IF(ans == "y" .OR. ans == "Y") THEN
        write(unit=*,fmt="(a,a,a)",advance="no") " => Code of the .res file ( <cr>= ",trim(outfil),") :"
        read(unit=*,fmt="(a)") resfil
        lr=len_trim(resfil)
        IF(lr == 0) resfil=outfil
        ansf=.true.
        open(unit=4,file=trim(resfil)//".res",status="replace",action="write")
        write(unit=*,fmt="(a)",advance="no")" => Give a comment for the .res file: "
        read(unit=*,fmt="(a)") title
        write(unit=4,fmt="(a)") title
      END IF

! Read input file infil.exc


!      Call read_exchange(lun1)

      write(unit=lun2,fmt="(a)")"   *****************************************************"
      write(unit=lun2,fmt="(a)")"   ****               PROGRAM ENERMAG               ****"
      write(unit=lun2,fmt="(a)")"   *****************************************************"
      write(unit=lun2,fmt="(a)")"  "
      write(unit=lun2,fmt="(a)")"                 *** Version 1.5 *** "
      write(unit=lun2,fmt="(a)")"   *****************************************************"
      write(unit=lun2,fmt="(a)")"   * Calculates the Magnetic Energy for k-vectors in BZ*"
      write(unit=lun2,fmt="(a)")"   *      Uses a Classical Heisenberg Hamiltonian      *"
      write(unit=lun2,fmt="(a)")"   *        and isotropic exchange interactions        *"
      write(unit=lun2,fmt="(a)")"   *****************************************************"
      write(unit=lun2,fmt="(a)") "                   (JRC August-2001, LLB)"
      write(unit=lun2,fmt="(a)") " "
      write(unit=lun2,fmt="(a)") " "
      write(unit=lun2,fmt="(a)") " "

      read(unit=1,fmt="(a)") title
      write(unit=*,fmt="(a)") title
      IF(ansf) write(unit=4,fmt="(a)") title
      write(unit=lun2,fmt="(/,a,/)") title
      write(unit=lun7,fmt="(a)")    "      OUTPUT OF THE SPIN CONFIGURATIONS"
      write(unit=lun7,fmt="(/,a,/)")  " => Warning!, this output is provisional!"
      write(unit=lun7,fmt="(a,a,/)") " => Title:",title
      read(unit=1,fmt="(i4,f8.4)")natcel,dmax
      IF(ansf) write(unit=4,fmt="(i5,a)") natcel, "   <--- Number of atoms"
      write(unit=lun2,fmt="(a,i3)") " => Number of atoms: ",natcel
      write(unit=*,fmt="(a,i3)") " => Number of atoms: ",natcel
      write(unit=lun2,fmt="(a,f8.3)") " => Maximum interatomic distance: ",dmax
      write(unit=*,fmt="(a,f8.3)") " => Maximum interatomic distance: ",dmax
      read(unit=1,fmt="(a)") spg
      write(unit=lun2,fmt="(a,a)") " => Space group: ",spg
      write(unit=*,fmt="(a,a)") " => Space group: ",spg
      read(unit=1,fmt="(3f8.4,3f8.3)") (ad(j),j=1,6)
      write(unit=lun2,fmt="(a,/,a,3f8.4,3f8.3)") "     ",   &
          " => Unit cell parameters (Angstroms & degrees):",(ad(j),j=1,6)
      write(unit=*,fmt="(a,/,a,3f8.4,3f8.3)") "     ", &
          " => Unit cell parameters (Angstroms & degrees):",(ad(j),j=1,6)
      write(unit=lun2,fmt="(/,a,/)") " => Atoms in the unit cell"
      write(unit=*,fmt="(/,a,/)") " => Atoms in the unit cell"
      write(unit=lun2,fmt="(a,/)") "  No. Name      x        y        z       Spin"
      write(unit=*,fmt="(a,/)") "  No. Name      x        y        z       Spin"
      DO i=1,natcel
        read(unit=1,fmt="(a6,4f9.5)")rnam(i),(xyz(j,i),j=1,3),s(i)
        IF(s(i) == 0.0) s(i)=1.0
        write(unit=lun2,fmt="(i4,a,a6,4f9.5)")i,"  ",rnam(i),(xyz(j,i),j=1,3),s(i)
        write(unit=*,fmt="(i4,a,a6,4f9.5)")i,"  ",rnam(i),(xyz(j,i),j=1,3),s(i)
      END DO
      write(unit=lun2,fmt="(/,a)") "--------------------------------------------"
      write(unit=lun2,fmt="(  a)") " => Input elements for constructing J_k(i,j)"
      write(unit=*,   fmt="(  a)") " => Input elements for constructing J_k(i,j)"
      write(unit=lun2,fmt="(a,/)") "--------------------------------------------"
      DO i=1,natcel
        DO j=1,natcel
          read(unit=1,fmt="(2i3,i4)") inn,jnn,nterm(i,j)
          IF(inn /= i .AND. jnn /= j)  &
              write(unit=*,fmt="(a)") " => Warning!, check the input exchange matrix!"
          write(unit=lun2,fmt="(a,3(i2,a))") " => J(",i,",",j,")[K]   (",nterm(i,j), " terms)"
          write(unit=*,fmt="(a,3(i2,a))") " => J(",i,",",j,")[K]   (",nterm(i,j), " terms)"
          IF(nterm(i,j) /= 0) THEN
            DO nt=1,nterm(i,j)
              read(unit=1,fmt="(a)")  line
              k=index(line,"J")
              if(k == 0) then
               write(unit=*,fmt="(a)") " => Error reading J in exchange file "
               stop
              end if
              read(unit=line(1:k-1),fmt=*) (trans(m,i,j,nt),m=1,3),exch
              jota=line(k:k+3)
              read(unit=line(k+1:k+3),fmt=*) nvalj(i,j,nt)
              name_jota=line(k+4:)
              !read(unit=1,fmt="(t5,3f9.5,f10.3,t45,2a)")  &
              !    (trans(m,i,j,nt),m=1,3),exch,jota,name_jota
              CALL get_expo(trans(:,i,j,nt),expo)
              write(unit=lun2,fmt="(a,a,3f6.2,a,f8.2,a,a,a)") "  ", &
                  " Rn=(", (trans(m,i,j,nt),m=1,3)," ) J= ",exch, " --> ",jota,expo
              write(unit=*,fmt="(a,a,3f6.2,a,f8.2,a,a,a)") "  ", &
                  " Rn=(", (trans(m,i,j,nt),m=1,3)," ) J= ",exch, " --> ",jota,expo
              IF(nvalj(i,j,nt) > nex) nex=nvalj(i,j,nt)
              valj(nvalj(i,j,nt))=exch
              jotas(nvalj(i,j,nt))=name_jota
            END DO
          END IF
        END DO
      END DO
!  End reading input file

      write(unit=*,fmt="(a)") " => Input file read!!!"
!  For small number of exchange parameters select a strategy
!  of changing them for constructing a phase diagram
    DO
      write(unit=lun2,fmt="(/,a)") "------------------------------------"
      write(unit=lun2,fmt="(  a)") " => Input set of exchange parameters"
      write(unit=lun2,fmt="(a,/)") "------------------------------------"
      write(unit=*,fmt="(/,a)") "------------------------------------"
      write(unit=*,fmt="(  a)") " => Input set of exchange parameters"
      write(unit=*,fmt="(a,/)") "------------------------------------"
      IF(ansf) THEN
        write(unit=4,fmt="(/,a)") "------------------------------------"
        write(unit=4,fmt="(  a)") " => Input set of exchange parameters"
        write(unit=4,fmt="(a,/)") "------------------------------------"
      END IF
      IF(nex < 10) THEN
        DO j=1,nex
          write(unit=*,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          write(unit=lun2,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          IF(ansf) write(unit=4,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
        END DO
        iphase=0
        write(unit=*,fmt="(a)",advance="no") " => Do you want to calculate a phase diagram (y/n)?: "
        read(unit=*,fmt="(a)") ans
        IF(ans == "y".OR.ans == "Y") iphase=1
      ELSE
        DO j=1,9
          write(unit=*,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          write(unit=lun2,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          IF(ansf) write(unit=4,fmt="(a,i1,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
        END DO
        DO j=10,nex
          write(unit=*,fmt="(a,i2,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          write(unit=lun2,fmt="(a,i2,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
          IF(ansf) write(unit=4,fmt="(a,i2,a,f8.3,a,a)") "  J",j,"=",valj(j),"  ",jotas(j)
        END DO
        iphase=0
      END IF

! Ask for changing J"s if no phase diagram is generated

      IF(iphase == 0) THEN
        write(unit=*,fmt="(a)",advance="no") " => Change exchange paramenters (y/n)?: "
        read(unit=*,fmt="(a)") ans
        IF(ans == "y".OR.ans == "Y") THEN
          DO j=1,id
            write(unit=*,fmt="(a)",advance="no") " -> Enter number of J and value (0 0: stop): "
            read(unit=*,fmt=*) inum, rvalj
            IF(inum == 0) EXIT
            valj(inum)=rvalj
          END DO
          cycle
        END IF
      END IF
      IF(iphase == 1) THEN
        write(unit=*,fmt="(/,a,/)")" => Phase diagram will be calculated "
        write(unit=lun2,fmt="(/,a,/)")" => Phase diagram will be calculated "
        write(unit=*,fmt="(a)",advance="no")" => Enter a target k-vector (3 reals): "
        read(unit=*,fmt=*) (ktar(i),i=1,3)
        DO i=1,ijo
          rang1(i)=0
          rang2(i)=0
          npoi(i)=1
        END DO
        write(unit=*,fmt="(a,i2,a)",advance="no") " => Number of J's TO BE VARIED (<=",ijo,"): "
        read(unit=*,fmt=*) nojvar
        do
          njjj=1
          DO i=1,nojvar
            write(unit=*,fmt="(a,i1,a)",advance="no")  &
                " -> Number & Range of varied J-parameter ",i," & No of points: "
            read(unit=*,fmt=*) ivar(i),rang1(i),rang2(i),npoi(i)
            njjj=njjj*npoi(i)
          END DO
          IF(njjj > njo) THEN
            write(unit=*,fmt="(a,i5)")  &
                " => Too many J-parameters! Reduce the number of points!"
            cycle
          END IF
          exit
        end do
        CALL genj(nojvar,njotas,rang1,rang2,npoi,valjj)
      ELSE
        njotas=1
        nojvar=0
      END IF
      IF(ansf) THEN
        write(unit=4,fmt="(i7,a)")  njotas,"   <--- Number of sets of J "
        DO i=1,nojvar
          write(unit=4,fmt="(a,i1,a, i2,2f10.4,i6)")  &
              " -> Number & Range of varied J-parameter ",i," & No of points: ",  &
              ivar(i),rang1(i),rang2(i),npoi(i)
        END DO
      END IF
      ishifp=njotas/24
      if(ishifp == 0) ishifp = 1
      IF(iphase == 1) then
          write(unit=*,fmt="(a,i7,a)") " => Phase diagram for: ",  &
          njotas," sets of exchange integrals"
          write(unit=lun2,fmt="(a,i5,a)") " => Phase diagram for: ",  &
          njotas," sets of exchange integrals"
        DO i=1,nojvar
          write(unit=lun2,fmt="(a,i1,a, i2,2f10.4,i6)")  &
              " -> Number & Range of varied J-parameter ",i," & No of points: ",  &
              ivar(i),rang1(i),rang2(i),npoi(i)
        END DO
      END IF
      iwrt=0
!------------------------------------------------------------------
      DO jj=1,njotas            !Start phase diagram
!------------------------------------------------------------------

        DO j=1,nojvar         !this loop is executed only if iphase>0
          valj(ivar(j))=valjj(j,jj)
        END DO
        IF(iphase == 0) THEN
          write(unit=lun7,fmt="(5(a,i1))") ("        J",j,j=1,nex)
          write(unit=lun7,fmt="(5f10.2)") (valj(j),j=1,nex)
        END IF
        emin=9999999.0E+20
!-------------------------------------------
!  Select strategy of searching of k-vectors
!-------------------------------------------
        IF(jj == 1) THEN
          CALL genk(lun2,ln,infil,nvk,vk,iprint)
          DO ik=1,nvk
            nvect(ik)=0
          END DO
         IF(iphase == 1) THEN
          forma="(/,a, (a,i1),a,/)"
          write(unit=forma(6:6),fmt="(i1)") nojvar
          if(ansf) then
           write(unit=4,fmt=forma)"       Energy   Vector  Kx    Ky    Kz "  &
              , ("      J",ivar(j),j=1,nojvar),"                    Eigenvectors "
           write(unit=lun2,fmt=forma)"       Energy   Targetk Kx    Ky    Kz "  &
              , ("      J",ivar(j),j=1,nojvar),"                    Eigenvectors "
          END IF
         END IF
         call cpu_time(cpu_ini)
         cpu=cpu_ini
        END IF
        IF(iphase == 0) THEN
          write(unit=lun2,fmt="(/,a)")"-------------------------------------------------"
          write(unit=lun2,fmt="(  a)")" => Calculation of energy for different k-vectors"
          write(unit=lun2,fmt="(a,/)")"-------------------------------------------------"
          IF(ansf) write(unit=4,fmt="(i6,a)") nvk, "   <--- Number of k-vectors "
          iprint=1
        END IF
!--------------------
        DO ik=1,nvk    !Loop over k-vectors
!--------------------
!----------------------------
!  Construct the matrix J(k)
!----------------------------
          rk(:)=vk(:,ik)
          CALL j_k(rk,exchk,iex)
!-----------------------------------------------------------
!  Diagonalize exc, first check if imaginary component exist
!-----------------------------------------------------------
          n=natcel
          IF(iex == 0) THEN
            aux(1:natcel,1:natcel)=real(exchk(1:natcel,1:natcel))
            if(iprint == 0 ) then
              Call diagonalize_sh(aux,natcel,eigen_val)
            else
              Call diagonalize_sh(aux,natcel,eigen_val,eigen_rvectors)
            end if
          ELSE
            if(iprint == 0) then
              Call diagonalize_sh(exchk,natcel,eigen_val)
            else
              Call diagonalize_sh(exchk,natcel,eigen_val,eigen_cvectors)
            end if
          END IF

          IF(ansf .AND. iphase == 0) THEN
            forma="(3f8.4, a,  f14.2)"
            write(unit=forma(11:12),fmt="(i2)") natcel
            write(unit=4,fmt=forma)(rk(i),i=1,3)," ",(eigen_val(i),i=1,natcel)
          END IF

          IF(eigen_val(natcel) < emin ) THEN
            IF(iphase == 0) THEN
              write(unit=*,fmt="(a,3f8.4,a)") " => Eigenvalues for k = (",(rk(i),i=1,3)," )"
              write(unit=*,fmt="(99f12.3)") eigen_val(1:natcel)
              write(unit=lun2,fmt="(a,3f8.4,a)") " => Eigenvalues & Eigenvectors for k = (",(rk(i),i=1,3)," )"
            !  write(unit=11,fmt="(a,3f8.4,a)") " => Eigenvalues & Eigenvectors for k = (",(rk(i),i=1,3)," )"
              forma="(f14.4,a,  f12.4,a)"
              j=natcel
              if(iex /= 0) j=2*j
              if(j < 10) then
                write(unit=forma(11:11),fmt="(i1)") j
              else
                write(unit=forma(10:11),fmt="(i2)") j
              end if
              if(iex /= 0) then
               do i=1,natcel
                 write(unit=lun2,fmt=forma ) eigen_val(i)," : (",eigen_cvectors(1:natcel,i),")"
               !  write(unit=11,fmt=forma ) eigen_val(i)," : (",eigen_cvectors(1:natcel,i),")"
               end do
              else
               do i=1,natcel
                 write(unit=lun2,fmt=forma ) eigen_val(i)," : (",eigen_rvectors(1:natcel,i),")"
               !  write(unit=11,fmt=forma ) eigen_val(i)," : (",eigen_rvectors(1:natcel,i),")"
               end do
              end if
            END IF
            IF(eigen_val(natcel) < emin) THEN     !Store the eigenvalue and eigenvectors
              iopt=ik                             !corresponding to the minimum energy
              emin=eigen_val(natcel)              !Update emin
              nop=natcel
              compl=.false.
              if(iex /= 0) then
                nop=nop*2
                compl=.true.
              end if
              aop(1:natcel,1:natcel)=exchk(1:natcel,1:natcel)
            END IF
          END IF
!--------------------
        END DO    !End Loop over k-vectors ik=1,nvk
!--------------------
!-------------------------------------------------------------
! Calculation of the spin configuration for the optimum energy
!-------------------------------------------------------------
        IF(compl) THEN
          Call diagonalize_sh(aop,natcel,eigen_val,eigen_cvectors)
        ELSE
          aux(1:natcel,1:natcel)=real(aop(1:natcel,1:natcel))
          Call diagonalize_sh(aux,natcel,eigen_val,eigen_rvectors)
        END IF

        IF(iphase == 0) THEN   !write only if no phase diagram is calculated

          write(unit=lun2,fmt="(/,a)") "--------------------------------------------------------"
          write(unit=lun2,fmt="(  a)") " => Exchange interaction matrix for the optimum k-vector"
          write(unit=lun2,fmt="(a,/)") "--------------------------------------------------------"
          if(compl) then
            forma="( a,  (i15,tr5))"
          else
            forma="( a,  i10)"
          end if
          write(unit=forma(5:6),fmt="(i2)") natcel
          write(unit=lun2,fmt=forma) "     ",(i,i=1,natcel)
          forma="(i3, a,  f10.2)"
          write(unit=forma(8:9),fmt="(i2)") nop
          if(compl) then
            DO i=1,natcel
              write(unit=lun2,fmt=forma)i,"  ",(aop(i,j),j=1,natcel)
            END DO
          else
            DO i=1,natcel
              write(unit=lun2,fmt=forma)i,"  ",(real(aop(i,j)),j=1,natcel)
            END DO
          end if
          write(unit=lun2,fmt="(/,a)") "-------------------------------------------------------"
          write(unit=lun2,fmt="(  a)") " => Eigenvalues and eigenvector of the optimum k-vector"
          write(unit=lun2,fmt="(a,/)") "-------------------------------------------------------"
          write(unit=lun2,fmt="(a,3f8.4,a)") " => Optimum k-vector -> (",(vk(i,iopt),i=1,3)," )"
          forma="(a,i2,a,  f10.4,a)"
          write(unit=forma(9:10),fmt="(i2)") nop
          if(compl) then
             DO i=1,natcel
               write(unit=lun2,fmt="(a,i2,a,f14.4)")  "    Eigenvalue  E(k,",i,") = ",eigen_val(i)
               write(unit=lun2,fmt=forma) "    Eigenvector V(k,",i,") = (",(eigen_cvectors(j,i),j=1,natcel),")"
             END DO
          else
             DO i=1,natcel
               write(unit=lun2,fmt="(a,i2,a,f14.4)")  "    Eigenvalue  E(k,",i,") = ",eigen_val(i)
               write(unit=lun2,fmt=forma) "    Eigenvector V(k,",i,") = (",(eigen_rvectors(j,i),j=1,natcel),")"
             END DO
          end if
        ELSE
!-------------- Below is for a phase diagram to be calculated
          IF(.NOT. ansf) THEN
            write(unit=lun2,fmt="(5(a,i1))") ("        J",j,j=1,nex)
            write(unit=lun2,fmt="(5f10.2)")     (valj(j),j=1,nex)
            write(unit=lun2,fmt="(a,3f8.4,a)")  &
                " => Optimum k-vector -> (",(vk(i,iopt),i=1,3)," )"
            write(unit=lun2,fmt="(a,f14.4)")  "    Emin = ",eigen_val(nop)
            forma="(a,i2,a,  f10.4,a)"
            write(unit=forma(9:10),fmt="(i2)") nop
            if(compl) then
                write(unit=lun2,fmt=forma)  &
                "    Eigenvector V(k,",natcel,") = (",(eigen_cvectors(j,natcel),j=1,natcel),")"
            else
                write(unit=lun2,fmt=forma)  &
                "    Eigenvector V(k,",natcel,") = (",(eigen_rvectors(j,natcel),j=1,natcel),")"
            end if
!--------------
          ELSE
            forma="(f14.2,   i6, a,3f6.2, a,  f8.2,   a,  f10.4,a)"
            write(unit=forma(26:27),fmt="(i2)") nojvar
            write(unit=forma(38:39),fmt="(i2)") nop
            if(compl) then
                write(unit=4,fmt=forma) eigen_val(nop),iopt," ",(vk(i,iopt),i=1,3),"  ",(valj(ivar(j)),j=1,nojvar),  &
                "  (",(sign(1.0,real(eigen_cvectors(1,natcel)))*eigen_cvectors(j,natcel),j=1,natcel),")"
            else
                write(unit=4,fmt=forma) eigen_val(nop),iopt," ",(vk(i,iopt),i=1,3),"  ",(valj(ivar(j)),j=1,nojvar),  &
                "  (",(sign(1.0,eigen_rvectors(1,natcel))*eigen_rvectors(j,natcel),j=1,natcel),")"
            end if
            nvect(iopt)=nvect(iopt)+1
            IF(MOD(jj,ishifp) == 0) THEN
              iwrt=iwrt+1
              call cpu_time(cpu_par)
              cpu= (cpu_par-cpu)/60.0
              remaining=(24-iwrt)*cpu/60.0
              cmdline(1:80)=" "
              write(unit=cmdline(1:10),fmt="(i3,a)") int(remaining)," hours,"
              remaining=(remaining-int(remaining))*60.0
              write(unit=cmdline(12:21),fmt="(i2,a)") int(remaining)," minutes"
              forma="(a,i6, a,i4, a,3f6.2,  f8.2)"
              write(unit=forma(22:23),fmt="(i2)") nojvar
              write(unit=*,fmt=forma) " -> Point #:",jj,"  ",iopt," ",(vk(i,iopt),i=1,3),  &
                                  (valj(ivar(j)),j=1,nojvar)
              write(unit=*,fmt="(a,f8.2,a,a,/)")"    Partial time:",cpu, &
              " minutes  ->  Approx. remaining time: ",trim(cmdline)
              cpu=cpu_par
            END IF
            ! Test for target vector
            IF(sum(ABS(ktar(:)-vk(:,iopt))) < 0.02) THEN
!                  12345678901234567890123456789012345678901234567890
             forma="(f14.4, a,i4, a,3f6.2, a,  f8.2,   a,  f10.4,a)"
             write(unit=forma(26:27),fmt="(i2)") nojvar
             write(unit=forma(38:39),fmt="(i2)") nop
             if(compl) then     !complex

              if (real(eigen_cvectors(1,nop))  < 0.0 ) then
               write(unit=lun2,fmt=forma) eigen_val(natcel),"  ",iopt," ",(vk(i,iopt),i=1,3),"  ",  &
                 (valj(ivar(j)),j=1,nojvar),"  (",(-eigen_cvectors(j,natcel),j=1,natcel),")"
             !  forma="(f14.4,a,a,"//trim(forma(38:))
             ! write(*,*) " forma: ", forma
             ! do k=natcel-1,1,-1
             !   write(unit=lun2,fmt=forma) eigen_val(k),"                                                   ",  &
             !   "  (",(-eigen_cvectors(j,k),j=1,natcel),")"
             ! end do
              else
               write(unit=lun2,fmt=forma) eigen_val(natcel),"  ",iopt," ",(vk(i,iopt),i=1,3),"  ",  &
                 (valj(ivar(j)),j=1,nojvar),"  (",(eigen_cvectors(j,natcel),j=1,natcel),")"
             ! forma="(f14.4,a,a,"//trim(forma(38:))
             ! write(*,*) " forma: ", forma
             ! do k=natcel-1,1,-1
             !   write(unit=lun2,fmt=forma) eigen_val(k),"                                                   ",  &
             !   "  (",(eigen_cvectors(j,k),j=1,natcel),")"
             ! end do
              end if

             else ! no complex

              if (eigen_rvectors(1,nop)  < 0.0 ) then
               write(unit=lun2,fmt=forma) eigen_val(natcel),"  ",iopt," ",(vk(i,iopt),i=1,3),"  ",  &
                 (valj(ivar(j)),j=1,nojvar),"  (",(-eigen_rvectors(j,natcel),j=1,natcel),")"
             ! forma="(f14.4,a,a,"//trim(forma(38:))
             ! write(*,*) " forma: ", forma
             ! do k=natcel-1,1,-1
             !   write(unit=lun2,fmt=forma) eigen_val(k),"                                                   ",  &
             !   "  (",(-eigen_cvectors(j,k),j=1,natcel),")"
             ! end do
              else
               write(unit=lun2,fmt=forma) eigen_val(natcel),"  ",iopt," ",(vk(i,iopt),i=1,3),"  ",  &
                 (valj(ivar(j)),j=1,nojvar),"  (",(eigen_rvectors(j,natcel),j=1,natcel),")"
             ! forma="(f14.4,a,a,"//trim(forma(38:))
             ! write(*,*) " forma: ", forma
             ! do k=natcel-1,1,-1
             !   write(unit=lun2,fmt=forma) eigen_val(k),"                                                   ",  &
             !   "  (",(eigen_cvectors(j,k),j=1,natcel),")"
             ! end do
              end if

             end if  !complex
             rktar(:)=vk(:,iopt)
             compl_tar=compl
             if(compl_tar) then
               dtarc(1:natcel)=eigen_cvectors(1:natcel,natcel)
             else
               dtarr(1:natcel)=eigen_rvectors(1:natcel,natcel)
             end if
             energytar= eigen_val(natcel)
             noptar=nop
            END IF
          END IF
        END IF         !End of write only if no phase diagr...
        IF(iphase == 0) THEN
          rk(:)=vk(:,iopt)
          energy=eigen_val(natcel)
          if(compl) then
            CALL spin_conf(lun7,rk,energy,eigen_cvectors(1:natcel,natcel))
          else
            CALL spin_conf(lun7,rk,energy,eigen_rvectors(1:natcel,natcel))
          end if
        END IF
!----------------------------------
      END DO              !njotas
!----------------------------------
      IF(iphase == 1 .and. ansf) THEN

          if(compl) then
            CALL spin_conf(lun7,rktar,energytar,dtarc)
          else
            CALL spin_conf(lun7,rktar,energytar,dtarr)
          end if
        write(unit=lun2,fmt="(a)") " => Frequency of optimum k-vectors"
        DO ik=1,nv
          IF(nvect(ik) > 1) THEN
            write(unit=lun2,fmt="(a,3f6.2,a,i6,a,i7,a)")  &
                " k = (",(vk(j,ik),j=1,3),")   #:",ik," -> ",nvect(ik), " times"
          END IF
        END DO
      END IF
      call cpu_time(cpu_fin)
      !---- CPU times ----!
      write(unit=*,fmt="(/,/,t8,a,f8.2,a)") "Cpu Time: ", (cpu_fin - cpu_ini)/60.0, " min."
      write(unit=lun2,fmt="(/,/,t8,a,f8.2,a)") "Cpu Time: ", (cpu_fin - cpu_ini)/60.0, " min."
      write(unit=*,fmt="(/,a)",advance="no")" => Do you want to continue (y/n)?: "
      read(unit=*,fmt="(a)") ans
      IF(ans == "y".OR.ans == "Y") cycle
      exit
    END DO
      STOP
   END PROGRAM enermag
!-----------------------------------------------------------------------

