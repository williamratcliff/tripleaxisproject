!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: STRUCTURE_FACTOR_MODULE
!!----   INFO: Main module for Structure Factors Calculations
!!----
!!---- HISTORY
!!----    Update: January - 2004
!!----
!!----
!!---- DEPENDENCIES
!!----
!!----     Use Scattering_Chemical_Tables
!!----     Use Crystallographic_Symmetry,   only: Space_Group_Type
!!----     Use Reflections_Utilities,       only: Reflection_List_Type, HKL_R
!!----     Use Atom_Module,                 only: atom_list_type
!!----     Use Math_Gen,                    only: sp, tpi, atan2d
!!----     Use String_Utilities,            only: L_Case,U_Case
!!----
!!---- VARIABLES
!!--++    AF0                          [Private]
!!--++    AFP                          [Private]
!!--++    AFPP                         [Private]
!!--++    AJH                          [Private]
!!--++    BJH                          [Private]
!!----    ERR_MESS_SFAC
!!----    ERR_SFAC
!!--++    HR_TYPE                      [Private]
!!--++    HR                           [Private]
!!--++    HT                           [Private]
!!--++    SF_INITIALIZED               [Private]
!!--++    TH                           [Private]
!!----
!!---- PUBLIC PROCEDURES
!!----    Functions:
!!--++       FJ                        [Private]
!!----
!!----    Subroutines:
!!--++       CALC_TABLE_AB             [Private]
!!--++       CALC_TABLE_TH             [Private]
!!----       CALC_STRFACTOR
!!--++       CREATE_TABLE_AF0_XRAY     [Private]
!!--++       CREATE_TABLE_AFP_NEUTNUC  [Private]
!!--++       CREATE_TABLE_HR_HT        [Private]
!!----       INIT_STRUCTURE_FACTORS
!!----       MODIFY_SF
!!--++       SET_FIXED_TABLES          [Private]
!!----       STRUCTURE_FACTORS
!!--++       SUM_AB                    [Private]
!!--++       SUM_AB_NEUTNUC            [Private]
!!----       WRITE_STRUCTURE_FACTORS
!!----
!!
 Module Structure_Factor_Module

    !---- Use Modules ----!
    Use Math_Gen,                    only: sp, tpi, atan2d
    Use String_Utilities,            only: L_Case,U_Case
    Use Scattering_Chemical_Tables
    Use Crystallographic_Symmetry,   only: Space_Group_Type
    Use Reflections_Utilities,       only: Reflection_List_Type, HKL_R
    Use Atom_Module,                 only: atom_list_type

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Init_Structure_Factors, Structure_Factors, Modify_SF, &
              Write_Structure_Factors,Calc_StrFactor

    !---- List of private functions ----!
    private :: Fj

    !---- List of private subroutines ----!
    private :: Calc_Table_AB, Create_Table_AF0_Xray, Create_Table_AFP_NeutNuc, &
               Create_Table_HR_HT, Set_Fixed_Tables, Calc_Table_TH, Sum_AB,    &
               Sum_AB_NeutNuc

    !---- Definitions ----!

    !!--++
    !!--++ AF0
    !!--++     real(kind=sp), dimension(:,:), allocatable, private :: AF0
    !!--++
    !!--++     Array for Atomic Factor. The dimensions are
    !!--++           AF0(Natoms,NRef)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=sp), dimension(:,:), allocatable, private :: AF0

    !!--++
    !!--++ AFP
    !!--++     real(kind=sp), dimension(:), allocatable, private :: AFP
    !!--++
    !!--++     Array for real part of anomalous scattering form factor.
    !!--++     The dimension is: AFP(Natoms)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=sp), dimension(:), allocatable, private :: AFP

    !!--++
    !!--++ AFPP
    !!--++     real(kind=sp), dimension(:), allocatable, private :: AFPP
    !!--++
    !!--++     Array for imaginary part of anomalous scattering form factor.
    !!--++     The dimension is: AFPP(Natoms)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=sp), dimension(:), allocatable, private :: AFPP

    !!--++
    !!--++ AJH
    !!--++     real(kind=sp), dimension(:,:), allocatable, private :: Ajh
    !!--++
    !!--++     Array for Aj(h). The dimensions are
    !!--++           Ajh(Natoms,Nref)
    !!--++     where
    !!--++           F(h)=Sum_j[Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=sp), dimension(:,:), allocatable, private :: AJH

    !!--++
    !!--++ BJH
    !!--++     real(kind=sp), dimension(:,:), allocatable, private :: Bjh
    !!--++
    !!--++     Array for Bj(h). The dimensions are
    !!--++           Bjh(Natoms,Nref)
    !!--++     where
    !!--++           F(h)=Sum_j[Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=sp), dimension(:,:), allocatable, private :: BJH

    !!----
    !!---- ERR_MESS_SFAC
    !!----    character(len=150), public :: err_mess_sfac
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_sfac

    !!----
    !!---- ERR_SFAC
    !!----    logical, public :: err_sfac
    !!----
    !!----    Logical Variable indicating an error in SF_Calculation module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_sfac

    !!--++
    !!--++    Type :: HR_Type
    !!--++       integer,dimension(3) :: H
    !!--++    End Type HR_Type
    !!--++
    !!--++    (Private)
    !!--++    Define a H vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Type, Private :: HR_Type
       integer, dimension(3) :: H
    End Type HR_Type

    !!--++
    !!--++ HR
    !!--++     type(HR_Type), dimension(:,:), allocatable, private :: Hr
    !!--++
    !!--++     Array for HR Calculations. The dimension are
    !!--++           HR(Natoms,NRef)
    !!--++
    !!--++ Update: February - 2005
    !!
    type(HR_Type), dimension(:,:), allocatable, private :: HR

    !!--++
    !!--++ HT
    !!--++    real(kind=sp), dimension(:,:), allocatable, private :: Ht
    !!--++
    !!--++    Array for HT Calculations. The dimension are
    !!--++          HT(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(:,:), allocatable, private :: HT

    !!----
    !!---- SF_Initialized
    !!----    logical, private :: SF_Initialized
    !!----
    !!----  Logical Variable indicating if the module has been initialized.
    !!----
    !!---- Update: February - 2005
    !!
    logical, private :: SF_Initialized=.false.

    !!--++
    !!--++ TH
    !!--++    real(kind=sp), dimension(:,:), allocatable, private :: Th
    !!--++
    !!--++    Array for TH Calculations. The dimension are
    !!--++          TH(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=sp), dimension(:,:), allocatable, private :: TH

 Contains

    !---- Functions ----!

    !!--++
    !!--++ Function Fj(s,a,b,c)
    !!--++    real(kind=sp),             intent(in) :: s
    !!--++    real(kind=sp),dimension(4),intent(in) :: a
    !!--++    real(kind=sp),dimension(4),intent(in) :: b
    !!--++    real(kind=sp),             intent(in) :: c
    !!--++
    !!--++    (Private)
    !!--++    Atomic scattering factor calculation according to:
    !!--++       Fj(s)=Sum_i[Ai*exp(-Bi*s*s)] + C (i=1..4)
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Fj(s,a,b,c) Result(res)
       !---- Arguments ----!
       real(kind=sp),             intent(in) :: s
       real(kind=sp),dimension(4),intent(in) :: a
       real(kind=sp),dimension(4),intent(in) :: b
       real(kind=sp),             intent(in) :: c
       real(kind=sp)                         :: res

       !---- Local variables ----!
       integer :: i

       res=0.0
       do i=1,4
          res=res + a(i)*exp(-b(i)*s*s)
       end do
       res=res+c

       return
    End Function Fj

    !---- Subroutines ----!

    !!--++
    !!--++ Subroutine Calc_Table_AB(Nref,Atm,Grp)
    !!--++    integer,                            intent(in) :: Nref
    !!--++    type(atom_list_type),              intent(in) :: Atm
    !!--++    type(space_group_type),             intent(in) :: Grp
    !!--++
    !!--++    (Private)
    !!--++    Calculate Table with Aj(h) and Bj(h) values
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Calc_Table_AB(Nref,Atm,Grp)
       !---- Arguments ----!
       integer,                            intent(in) :: Nref
       type(atom_list_type),              intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp

       !---- Local Variables ----!
       integer                       :: i,j,k
       real(kind=sp)                 :: arg,anis
       real(kind=sp),dimension(3)    :: h
       real(kind=sp),dimension(6)    :: beta

       Ajh=0.0
       Bjh=0.0
       if(Grp%Centred == 2) then
         do j=1,Nref
            do i=1,Atm%natoms
               arg=0.0
               do k=1,grp%NumOps
                  h=hr(k,j)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  Ajh(i,j)=Ajh(i,j)+cos(arg)*anis
               end do ! symmetry
            end do ! Atoms
         end do ! Reflections
       else
         do j=1,Nref
            do i=1,Atm%natoms
               arg=0.0
               do k=1,grp%NumOps
                  h=hr(k,j)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  Ajh(i,j)=Ajh(i,j)+cos(arg)*anis
                  Bjh(i,j)=Bjh(i,j)+sin(arg)*anis
               end do ! symmetry
            end do ! Atoms
         end do ! Reflections
       end if

       return
    End Subroutine Calc_Table_AB

    !!----
    !!---- Subroutine Calc_StrFactor(mode,rad,nn,sn,Atm,Grp,sf2,deriv,fc)
    !!----    character(len=*),                   intent(in) :: mode !S-XTAL (S) or Powder (P)
    !!----    character(len=*),                   intent(in) :: rad  !Radiation: X-rays, Neutrons
    !!----    integer,                            intent(in) :: nn
    !!----    real,                               intent(in) :: sn !(sinTheta/Lambda)**2
    !!----    type(atom_list_type),               intent(in) :: Atm
    !!----    type(space_group_type),             intent(in) :: Grp
    !!----    real,                               intent(out):: sf2
    !!----    real,dimension(:),optional,         intent(out):: deriv
    !!----    complex, optional,                  intent(out):: fc
    !!----
    !!----    (Private)
    !!----    Calculate Structure Factor for reflection "nn" in the list
    !!----    and derivatives with respect to refined parameters
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_StrFactor(mode,rad,nn,sn,Atm,Grp,sf2,deriv,fc)
       !---- Arguments ----!
       character(len=*),                   intent(in) :: mode
       character(len=*),                   intent(in) :: rad
       integer,                            intent(in) :: nn
       real,                               intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       real,                               intent(out):: sf2
       real,dimension(:),optional,         intent(out):: deriv
       complex, optional,                  intent(out):: fc

       !---- Local Variables ----!
       character(len=1)                      :: modi
       integer                               :: i,j,k,m
       real(kind=sp)                         :: arg,anis,cosr,sinr,scosr,ssinr,fr,fi,der
       real(kind=sp)                         :: a1,a2,a3,a4,b1,b2,b3,b4,av,bv,f
      real(kind=sp),dimension(3)            :: h
       real(kind=sp),dimension(6)            :: beta
       real(kind=sp),dimension(Atm%natoms)   :: frc,frs,otr,oti,afpxn
       real(kind=sp),dimension(9,Atm%natoms) :: drs,drc

       !--- Initialising local variables
       a1=0.0
       a2=0.0
       a3=0.0
       a4=0.0
       b1=0.0
       b2=0.0
       b3=0.0
       b4=0.0
       av=0.0
       bv=0.0
       fr=1.0
       fi=0.0
       frc=0.0
       frs=0.0
       otr=0.0
       oti=0.0
       modi=u_case(mode(1:1))
       if(rad(1:1) == "N") then
         afpxn(:)=afp(:)
       else
       	 afpxn(:)=af0(:,nn)
       end if

       if(Grp%Centred == 2) then
            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=hr(k,nn)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,nn))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     sinr=SIN(arg)*anis*fr   !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr
                  end if

               end do ! symmetry

               frc(i) = scosr
               otr(i) = afpxn(i)*th(i,nn)
               oti(i) =  afpp(i)*th(i,nn)
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       else

            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=hr(k,nn)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,nn))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  sinr=SIN(arg)*anis*fr     !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
                  ssinr=ssinr+sinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drs(1:3,i)=drs(1:3,i)+h(1:3)*cosr      ! +

                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr

                     drs(4,i)=drs(4,i)+h(1)*h(1)*sinr
                     drs(5,i)=drs(5,i)+h(2)*h(2)*sinr
                     drs(6,i)=drs(6,i)+h(3)*h(3)*sinr
                     drs(7,i)=drs(7,i)+h(1)*h(2)*sinr
                     drs(8,i)=drs(8,i)+h(1)*h(3)*sinr
                     drs(9,i)=drs(9,i)+h(2)*h(3)*sinr
                  end if

               end do ! symmetry

               frc(i) = scosr
               frs(i) = ssinr
               otr(i) = afpxn(i)*th(i,nn)
               oti(i) =  afpp(i)*th(i,nn)
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)
               a3 = a3 + oti(i)*frs(i)
               b3 = b3 + otr(i)*frs(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       end if

       If(modi == "P") then
          sf2 = a1*a1 + a2*a2 + a3*a3 + a4*a4 + b1*b1 + b2*b2 + b3*b3 + b4*b4
          sf2 = sf2 + 2.0*(b1*b4 -  a1*a4 + a2*a3 - b2*b3)
       else
          sf2= av*av+bv*bv
       End if

       if(present(fc)) then
         fc=cmplx(av,bv)
       end if

       if(present(deriv)) then

         if(modi == "P") then

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  POWDER
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der= otr(i)*(-a1*drc(m,i)+b3*drs(m,i))+oti(i)*(-b1*drc(m,i)+a3*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der= otr(i)*(a1*frc(i) +b3*frs(i))+oti(i)*(b1*frc(i) +a3*frs(i))
                   der=-2.0*der*sn
                   deriv(k) = sign(1.0,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor   POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der= otr(i)*(a1*frc(i)+b3*frs(i))+oti(i)*(b1*frc(i)+a3*frs(i))
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0,f)*der+deriv(k)
                 end if

                 do m=4,9      !Derivatives w.r.t. anisotropic temperature factors   POWDER
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=  otr(i)*(a1*drc(i,j)+b3*drs(m,i))+oti(i)*(b1*drc(m,i)+a3*drs(m,i))
                      der=-2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0,f)*der+deriv(k)
                    end if
                 end do

             end do

         else

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  S-XTAL
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der=   -av*(otr(i)*drc(m,i) + oti(i)*drs(m,i))
                     der=der-bv*(oti(i)*drc(m,i) - otr(i)*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der=   -av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der-bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der*sn
                   deriv(k) = sign(1.0,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der=    av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der+bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0,f)*der+deriv(k)
                 end if

                 do m=4,9        !Derivatives w.r.t. anisotropic temperature factors S-XTAL
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=   -av*(otr(i)*drc(m,i) - oti(i)*drs(m,i))
                      der=der-bv*(oti(i)*drc(m,i) + otr(i)*drs(m,i))
                      der=2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0,f)*der+deriv(k)
                    end if
                 end do

             end do
         end if !modi
       end if  !derivatives

       return
    End Subroutine Calc_StrFactor

    !!--++
    !!--++ Subroutine Calc_Table_TH(Reflex,Atm)
    !!--++    type(reflection_List_type),   intent(in) :: Reflex
    !!--++    type(atom_list_type),        intent(in) :: Atm
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Isotropinc Thermal contribution and occupation
    !!--..         TH(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Calc_Table_TH(Reflex,Atm)
       !---- Argument ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(atom_list_type),      intent(in) :: Atm

       !---- Local variables ----!
       integer          :: i,j
       real(kind=sp)    :: b,s

       !---- Isotropic model ----!
       do j=1,reflex%nref
          s=reflex%ref(j)%s
          do i=1,atm%natoms
             b=atm%atom(i)%biso
             th(i,j)= atm%atom(i)%occ * exp(-b*s*s)
          end do
       end do

       return
    End Subroutine Calc_Table_TH

    !!--++
    !!--++ Subroutine Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
    !!--++    type(reflection_List_type), intent(in) :: Reflex
    !!--++    type(atom_list_type),      intent(in) :: Atm
    !!--++    real(kind=sp), optiona      intent(in) :: lambda
    !!--++    integer, optional,          intent(in) :: lun
    !!--++
    !!--++    Calculate a Table of Atomic Factors for X-Ray
    !!--..      AF0(Natoms,Nref), AFP(Natoms), AFPP(Natoms)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
       !---- Arguments ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(atom_list_type),      intent(in) :: Atm
       real(kind=sp), optional,    intent(in) :: lambda
       integer, optional,          intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)               :: symbcar
       integer                        :: i,j, k,n,L
       integer, dimension(atm%natoms) :: ix,jx,ia
       real(kind=sp)                  :: dmin,d

       !---- Init ----!
       err_sfac=.false.

       !---- Load form factor values for XRay ----!
       call Set_Xray_Form()

       !---- Found Species on Xray_Form ----!
       ix=0
       jx=0
       n=0
       do i=1,atm%natoms
          symbcar=l_case(atm%atom(i)%chemsymb)
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             ix(i)=j
             if(any(jx == j) ) exit
             n=n+1
             jx(n)=j
             ia(n)=i
             exit
          end do
       end do

       if (present(lun)) then
          write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)") "  ==================================================="
       End if
       if (present(lambda)) then
          !---- Load anomalous scattering form factor values for XRays ----!
          call Set_Delta_Fp_Fpp()

          !---- Select wavelength (by default is CuKalpha1: k=5 in the list) ----!
          dmin=1000.0
          do i=1,5
             d=abs(lambda-Xray_Wavelengths(i)%Kalfa(1))
             if (d < dmin) then
                dmin=d
                k=i        !Selection of the index for fp and fpp lists
             end if
          end do

          !---- Found Species on Anomalous_ScFac ----!
          do i=1,atm%natoms
             symbcar=l_case(atm%atom(i)%chemsymb)
             do j=1,Num_Delta_Fp
                if (symbcar /= Anomalous_ScFac(j)%Symb) cycle
                afp(i)=Anomalous_ScFac(j)%fp(k)
                afpp(i)=Anomalous_ScFac(j)%fpp(k)
                exit
             end do
          end do
          call Remove_Delta_Fp_Fpp()
       else
           if (present(lun)) then
             write(unit=lun,fmt="(a)")    "  Missed lambda, anomalous dipersion corrections not applied   "
             write(unit=lun,fmt="(a)")    "  The default wavelength is that of Cu-Kalpha1 spectral line  "
           end if
       end if

       if (any(ix==0)) then
          err_sfac=.true.
          err_mess_sfac="The Species "//symbcar//" was not found"
       else
          !---- Fill AF Table ----!
          do j=1,reflex%nref
             do i=1,atm%natoms
                af0(i,j)=fj(reflex%ref(j)%s,xray_form(ix(i))%a,xray_form(ix(i))%b,xray_form(ix(i))%c)+afp(i)
             end do
          end do
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  Dfp  Dfpp "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               "   Atom     a1       b1       a2       b2       a3       b3       a4       b4        c      Dfp     Dfpp"
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,11F9.5)")    &
                           "     "//atm%atom(i)%chemsymb, &
                           (xray_form(j)%a(L),xray_form(j)%b(L), L=1,4), xray_form(j)%c, &
                           afp(i), afpp(i)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_Xray_Form()

       return
    End Subroutine Create_Table_AF0_Xray

    !!--++
    !!--++ Subroutine Create_Table_AFP_NeutNuc(Atm,lun)
    !!--++    type(atom_list_type),              intent(in) :: Atm
    !!--++    integer, optional,                  intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Setting a Table of Fermi Lengths for Neutron Nuclear Scattering
    !!--..      AFP(Natoms)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_AFP_NeutNuc(Atm,lun)
       !---- Arguments ----!
       type(atom_list_type),              intent(in) :: Atm
       integer, optional,                  intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)                        :: symbcar
       integer                                 :: i,k,n
       character(len=4), dimension(atm%natoms) :: symb
       real(kind=sp),    dimension(atm%natoms) :: bs
       real(kind=sp)                           :: b

       !---- Init ----!
       err_sfac=.false.

       !---- Load chemical information ----!
       call set_chem_info()

       !---- Getting Fermi Lengths of atoms ----!
       symb="    "
       bs=0.0
       n=0
       do i=1,atm%natoms
          symbcar=u_case(atm%atom(i)%chemsymb)
          call Get_Fermi_Length(symbcar,b)
          if (abs(b) < 0.0001) then
             err_sfac=.true.
             err_mess_sfac="The Fermi Length of Species "//symbcar//" was not found"
             return
          else
             afp(i) = b
             if(any(symb == symbcar)) cycle
             n=n+1
             symb(n)=symbcar
             bs(n) = b
          end if
       end do

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED NEUTRON SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)")  "  ==================================================="
          write(unit=lun,fmt="(a)")    "  FERMI LENGTHS "
          write(unit=lun,fmt="(a,i3)") "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)")  "   Atom     Fermi Length [10^(-12) cm]"
          do k=1,n
             write(unit=lun,fmt="(a,F15.6)")  "     "//symb(k), bs(k)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_chem_info()

       return
    End Subroutine Create_Table_AFP_NeutNuc

    !!--++
    !!--++ Subroutine Create_Table_HR_HT(Reflex,Grp)
    !!--++    type(reflection_list_type), intent(in) :: Hkl
    !!--++    type(space_group_type),     intent(in) :: Grp
    !!--++
    !!--++    Calculate a Table with HR and HT values
    !!--..       Hr(Grp%Numops,Reflex%Nref)
    !!--..       HT(Grp%Numops,Reflex%Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_HR_HT(Reflex,Grp)
       !---- Arguments ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(space_group_type),     intent(in) :: Grp

       !---- Local Variables ----!
       integer :: i,j

       do j=1,reflex%nref
          do i=1,grp%NumOps
             hr(i,j)%h=Hkl_R(reflex%ref(j)%h,grp%symop(i))
             ht(i,j)=dot_product(real(reflex%ref(j)%h),Grp%SymOp(i)%Tr)
          end do
       end do

       return
    End Subroutine Create_Table_HR_HT

    !!----
    !!---- Subroutine Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda,lun)
    !!----    type(reflection_list_type),          intent(in) :: Reflex
    !!----    type(atom_list_type),               intent(in) :: Atm
    !!----    type(space_group_type),              intent(in) :: Grp
    !!----    character(len=*),          optional, intent(in) :: Mode
    !!----    real(kind=sp),             optional, intent(in) :: lambda
    !!----    integer,                   optional, intent(in) :: lun  !Logical unit for writing scatt-factors
    !!----
    !!----    Allocates and initializes arrays for Structure Factors calculations.
    !!----    A calculation of fixed tables is also performed.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda,lun)
       !---Arguments ---!
       type(reflection_list_type),          intent(in) :: Reflex
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),              intent(in) :: Grp
       character(len=*),          optional, intent(in) :: Mode
       real(kind=sp),             optional, intent(in) :: lambda
       integer,                   optional, intent(in) :: lun

       !--- Local variables ---!

       integer :: Natm, Multr
       integer :: ierr

       err_sfac=.false.
       Natm = Atm%natoms
       Multr= Grp%Numops

       !---- Scattering factor tables ----!
       if (allocated(AF0)) deallocate(AF0)
       allocate(AF0(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for AF0"
          return
       end if
       AF0=0.0

       !---- Anomalous Scattering factor tables ----!
       if (allocated(AFP)) deallocate(AFP)
       allocate(AFP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for AFP"
          return
       end if
       AFP=0.0

       if (allocated(AFPP)) deallocate(AFPP)
       allocate(AFPP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for AFPP"
          return
       end if
       AFPP=0.0

       !---- HR Table ----!
       if (allocated(HR)) deallocate(HR)
       allocate(HR(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for HR"
          return
       end if
       HR=HR_Type(0)

       !---- HT Table ----!
       if (allocated(HT)) deallocate(HT)
       allocate(HT(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for HTR"
          return
       end if
       HT=0.0

       if (allocated(TH)) deallocate(TH)
       allocate(TH(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error on memory for HTR"
          return
       end if
       TH=0.0

       if (allocated(Ajh)) deallocate(Ajh)
       allocate(Ajh(Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error in Memory for Aj(h)"
          return
       end if
       Ajh=0.0

       if (allocated(Bjh)) deallocate(Bjh)
       allocate(Bjh(Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          err_mess_sfac="Error in Memory for Bj(h)"
          return
       end if
       Bjh=0.0

       if (present(mode)) then
          if (present(lambda)) then
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda,lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda)
             end if
          else
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode)
             end if
          end if
       else
          if (present(lambda)) then
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,lambda=lambda,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,lambda=lambda)
             end if
          else
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp)
             end if
          end if
       end if

       if (.not. err_sfac) SF_Initialized=.true.

       return
    End Subroutine Init_Structure_Factors

    !!----
    !!---- Subroutine Modify_SF(Reflex,Atm,Grp,List,Nlist,Mode)
    !!----    type(reflection_list_type),         intent(in out) :: Reflex
    !!----    type(atom_list_type),              intent(in)     :: Atm
    !!----    type(space_group_type),             intent(in)     :: Grp
    !!----    integer,dimension(:),               intent(in)     :: List
    !!----    integer,                            intent(in)     :: Nlist
    !!----    character(len=*),optional,          intent(in)     :: Mode
    !!----
    !!----    Recalculation of Structure Factors because a list of Atoms
    !!----    parameters were modified. List variable
    !!----    contains the number of atoms to be changed.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Modify_SF(Reflex,Atm,Grp,List,Nlist,partyp,Mode)
       !---- Arguments ----!
       type(reflection_list_type),   intent(in out) :: Reflex
       type(atom_list_type),         intent(in)     :: Atm
       type(space_group_type),       intent(in)     :: Grp
       integer,dimension(:),         intent(in)     :: List
       integer,                      intent(in)     :: NList
       character(len=*),optional,    intent(in)     :: partyp
       character(len=*),optional,    intent(in)     :: Mode

       !---- Local variables ----!
       character(len=2) :: typ
       integer          :: i,j,k,ii
       real(kind=sp)    :: arg,b,s

       typ="CO"
       if (present(partyp)) typ=adjustl(partyp)
       typ=U_Case(typ)

       select case (typ)

          case ("CO") ! by coordinates

            if(Grp%Centred == 2) then

               do j=1,Reflex%Nref
                  do ii=1,Nlist
                     i=list(ii)
                     Ajh(i,j)=0.0
                     arg=0.0
                     do k=1,grp%NumOps
                        arg=tpi*(dot_product(hr(k,j)%h,Atm%atom(i)%x)+ht(k,j))
                        Ajh(i,j)=Ajh(i,j)+cos(arg)
                     end do ! symmetry
                  end do ! NList
               end do ! Reflections

            else

               do j=1,Reflex%Nref
                  do ii=1,Nlist
                     i=list(ii)
                     arg=0.0
                     Ajh(i,j)=0.0
                     Bjh(i,j)=0.0
                     do k=1,grp%NumOps
                        arg=tpi*(dot_product(hr(k,j)%h,Atm%atom(i)%x)+ht(k,j))
                        Ajh(i,j)=Ajh(i,j)+cos(arg)
                        Bjh(i,j)=Bjh(i,j)+sin(arg)
                     end do ! symmetry
                  end do ! NList
               end do ! Reflections

            end if

          case ("TH") ! by thermal parameter or occupation number

             do j=1,Reflex%Nref
                s=reflex%ref(j)%s
                do ii=1,Nlist
                   i=list(ii)
                   b=atm%atom(i)%biso
                   th(i,j)=atm%atom(i)%occ*exp(-b*s*s)
                end do ! NList
             end do ! Reflections

       end select

       !---- Recalculation of SF ----!
       if(present(mode)) then
         if(mode == "XRA") then
            call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "NUC") then
            call Sum_AB_NeutNuc(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "MAG") then
         end if
       else
         call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
       end if


       return
    End Subroutine Modify_SF

    !!--++
    !!--++ Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,mode,lambda,lun)
    !!--++    type(reflection_list_type),         intent(in) :: Reflex
    !!--++    type(atom_list_type),              intent(in) :: Atm
    !!--++    type(space_group_type),             intent(in) :: Grp
    !!--++    character(len=*), optional,         intent(in) :: Mode
    !!--++    real(kind=sp), optional,            intent(in) :: lambda
    !!--++    integer, optional,                  intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculates arrays that are fixed during all further
    !!--++    calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda,lun)
       !---- Arguments ----!
       type(reflection_list_type),         intent(in) :: Reflex
       type(atom_list_type),              intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       character(len=*), optional,         intent(in) :: Mode
       real(kind=sp), optional,            intent(in) :: lambda
       integer, optional,                  intent(in) :: lun

       !---- Local variables ----!
       character(len=3) :: tipo

       tipo="XRA"
       if (present(mode)) tipo=adjustl(mode)
       tipo=U_Case(tipo)

       !---- Table HR - HT ----!
       call Create_Table_HR_HT(Reflex,Grp)

       !---- Table AF0 ----!
       select case (tipo)
          case ("XRA")
             if (present(lambda)) then
                if (present(lun)) then
                   call Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
                else
                   call Create_Table_AF0_Xray(Reflex,Atm,lambda)
                end if
             else
                if (present(lun)) then
                   call Create_Table_AF0_Xray(Reflex,Atm,lun=lun)
                else
                   call Create_Table_AF0_Xray(Reflex,Atm)
                end if
             end if

             !---- Modify the scattering factor tables to include the
             !---- multipliers factors concerning centre of symmetry and
             !---- centred translations
             if (Grp%Centred == 2) then
                af0=2.0*af0
                afpp=2.0*afpp
             end if

             if (Grp%NumLat  > 1) then
                af0=Grp%NumLat*af0
                afpp=Grp%NumLat*afpp
             end if

          case ("NUC")
             if (present(lun)) then
                call Create_Table_AFP_NeutNuc(Atm,lun=lun)
             else
                call Create_Table_AFP_NeutNuc(Atm)
             end if
             if (Grp%Centred == 2) afp=2.0*afp
             if (Grp%NumLat  > 1) afp=Grp%NumLat*afp

          case ("MAG")
       end select

       return
    End Subroutine Set_Fixed_Tables

    !!----
    !!---- Subroutine Structure_Factors(Atm,Grp,Reflex,Mode,lambda)
    !!----    type(atom_list_type),               intent(in)     :: Atm
    !!----    type(space_group_type),             intent(in)     :: Grp
    !!----    type(reflection_list_type),         intent(in out) :: Reflex
    !!----    character(len=*), optional,         intent(in)     :: Mode
    !!----    real(kind=sp), optional,            intent(in)     :: lambda
    !!----
    !!----    Calculate the Structure Factors from a list of Atoms
    !!----    and a set of reflections. A call to Init_Structure_Factors
    !!----    is a pre-requisite for using this subroutine. In any case
    !!----    the subroutine calls Init_Structure_Factors if SF_initialized=.false.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Structure_Factors(Atm,Grp,Reflex,Mode,lambda)
       !---- Arguments ----!
       type(atom_list_type),               intent(in)     :: Atm
       type(space_group_type),             intent(in)     :: Grp
       type(reflection_list_type),         intent(in out) :: Reflex
       character(len=*), optional,         intent(in)     :: Mode
       real(kind=sp), optional,            intent(in)     :: lambda

       !Provisional items
       ! integer::i,j
       !---------------
       if(present(Mode)) then
          if(present(lambda)) then
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda)
          else
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Mode)
          end if
       else
          if(present(lambda)) then
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Lambda=lambda)
          else
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp)
          end if
       end if

       !---- Table TH ----!
       Call Calc_Table_TH(Reflex,Atm)

       !---- Table AB ----!
       call Calc_Table_AB(Reflex%Nref,Atm,Grp)

       !Provisional items
       !open(unit=111,file="stfac.inf",status="replace",action="write")
       !do j=1,Nref
       !  write(111,"(a,3i4)") " Reflection:  ",hkl(j)%h
       !
       !  write(111,"(a)") " Atom              F0         occ*W         Ajh         Bjh"
       !  do i=1,Atm%natoms
       !       write(111,"(a,4f12.4)") "  "//atm%atom(i)%lab, af0(i,j),th(i,j),Ajh(i,j),Bjh(i,j)
       !  end do ! Atoms
       !end do ! Reflections
       !close(unit=111)
       !End Provisional items

       !---- Final Calculation ----!
       if(present(mode)) then
         if(mode == "XRA") then
            call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "NUC") then
            call Sum_AB_NeutNuc(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "MAG") then
         end if
       else
         call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
       end if
       return
    End Subroutine Structure_Factors

    !!--++
    !!--++ Subroutine Sum_AB(Reflex,Natm,icent)
    !!--++    type(reflection_list_type), intent(in out) :: Reflex
    !!--++    integer,                    intent(in)     :: Natm
    !!--++    integer,                    intent(in)     :: icent
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Final Sum for Structure Factors calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sum_AB(Reflex,Natm,icent)
       !---- Arguments ----!
       type(reflection_list_type), intent(in out)  :: Reflex
       integer,                    intent(in)      :: Natm
       integer,                    intent(in)      :: icent

       !---- Local Variables ----!
       integer                                     :: i,j
       real(kind=sp)                               :: a,b, ph
       real(kind=sp), dimension(natm,reflex%nref)  :: aa,bb,cc,dd


       ! A(h)=SIG(i){(f0+Deltaf')*OCC*Tiso*Ag}    asfa=a-d
       ! C(h)=SIG(i){    Deltaf" *OCC*Tiso*Ag}    bsfa=b+c

       ! B(h)=SIG(i){(f0+Deltaf')*OCC*Tiso*Bg}
       ! D(h)=SIG(i){    Deltaf" *OCC*Tiso*Bg}

       !---- Fj(h)*Aj(h) ----!

       aa=af0*th*ajh

       if (icent == 2) then    !Calculation for centrosymmetric structures
          do j=1,reflex%nref
             cc(:,j)= afpp(:)*th(:,j)*ajh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             b=sum(cc(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do

       else       !Calculation for non-centrosymmetric structures
          !---- Fj(h)*Bj(h) ----!
          bb=af0*th*bjh

          do j=1,reflex%nref
             cc(:,j)= afpp(:)*th(:,j)*ajh(:,j)
             dd(:,j)= afpp(:)*th(:,j)*bjh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i)-dd(:,i))
             b=sum(bb(:,i)+cc(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do
       end if

       return
    End Subroutine Sum_AB

    !!--++
    !!--++ Subroutine Sum_AB_NeutNuc(Reflex,Natm,icent)
    !!--++    type(reflection_list_type),         intent(in out) :: Reflex
    !!--++    integer,                            intent(in)     :: Natm
    !!--++    integer,                            intent(in)     :: icent
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Final Sum for Structure Factors calculations
    !!--++    Adapted for Neutron Nuclear Scattering (real scattering lengths)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sum_AB_NeutNuc(Reflex,Natm,icent)
       !---- Arguments ----!
       type(reflection_list_type),   intent(in out) :: Reflex
       integer,                      intent(in)     :: Natm
       integer,                      intent(in)     :: icent

       !---- Local Variables ----!
       integer                                     :: i,j
       real(kind=sp)                               :: a,b, ph
       real(kind=sp), dimension(natm,reflex%nref)  :: aa,bb

       if (icent == 2) then    !Calculation for centrosymmetric structures

          !---- Fj(h)*Aj(h) ----!
          do j=1,reflex%nref
             aa(:,j)= afp(:)*th(:,j)*ajh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             reflex%ref(i)%Fc=abs(a)
             reflex%ref(i)%Phase = 90.0 - 90.0 * sign(1.0,a)
             reflex%ref(i)%A=a
             reflex%ref(i)%B=0.0
          end do

       else       !Calculation for non-centrosymmetric structures
          !---- Fj(h)*Bj(h) ----!
          !---- Fj(h)*Aj(h) ----!
          do j=1,reflex%nref
             aa(:,j)= afp(:)*th(:,j)*ajh(:,j)
             bb(:,j)= afp(:)*th(:,j)*bjh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             b=sum(bb(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do
       end if

       return
    End Subroutine Sum_AB_NeutNuc

    !!----
    !!---- Subroutine Write_Structure_Factors(lun,Reflex,Mode)
    !!----    integer,                            intent(in) :: lun
    !!----    type(reflection_list_type),         intent(in) :: Reflex
    !!----    Character(len=*), optional,         intent(in) :: Mode
    !!----
    !!----    Writes in logical unit=lun the list of structure factors
    !!----    contained in the array hkl
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Structure_Factors(lun,Reflex,Mode)
       !---- Argument ----!
       integer,                            intent(in) :: lun
       type(reflection_list_type),         intent(in) :: Reflex
       Character(len=*), optional,         intent(in) :: Mode
       !---- Local Variables ----!
       integer :: i

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

       write(unit=lun,fmt="(/,a,/)") &
            "   H   K   L   Mult  SinTh/Lda     dspc        |Fc|       Phase        F-Real      F-Imag       |Fc|^2    Num"
       do i=1,reflex%Nref
             write(unit=lun,fmt="(3i4,i5,7f12.5,i8)") reflex%ref(i)%h, reflex%ref(i)%mult, &
                                 reflex%ref(i)%S,0.5/reflex%ref(i)%S, reflex%ref(i)%Fc, reflex%ref(i)%Phase,   &
                                 reflex%ref(i)%a, reflex%ref(i)%b, reflex%ref(i)%Fc*reflex%ref(i)%Fc,i
       end do
       return
    End Subroutine Write_Structure_Factors

 End Module Structure_Factor_Module
