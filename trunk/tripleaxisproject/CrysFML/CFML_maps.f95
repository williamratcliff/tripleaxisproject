!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: MAPS_CALCULATIONS
!!----   INFO: Subroutines related to operations on the array's
!!----
!!---- HISTORY
!!----    Update: January - 2004
!!----
!!----
!!---- DEPENDENCIES
!!----
!!--++    Use Math_Gen,                  only: sp
!!--++    Use Crystallographic_Symmetry, only: Space_Group_Type
!!--++    Use Crystal_Types,             only: Crystal_Cell_Type
!!--++    Use Geom_Calculations,         only: Distance
!!----
!!---- VARIABLES
!!----    ERR_MAPS
!!----    ERR_MESS_MAPS
!!----
!!---- PROCEDURES
!!----    Functions:
!!--++       PEAK_POSITION       [Private]
!!----
!!----    Subroutines:
!!----       INIT_ERR_MAPS
!!----       LOAD_EXTENDEDMAP
!!----       LOAD_SECTION
!!--++       PEAK_LIST           [Private]
!!----       SEARCH_PEAKS
!!----       STATISTIC_MAP
!!----
!!
 Module Maps_Calculations

    !---- Use Modules ----!
    !Use Mod_fun    !To be commented for non-F compilers
    use Math_Gen,                  only: sp
    use Crystallographic_Symmetry, only: Space_Group_Type, ApplySO
    use Crystal_Types,             only: Crystal_Cell_Type
    use Geom_Calculations,         only: Distance

    implicit none

    private

    !---- List of public subroutines ----!
    public :: Init_Err_Maps, Load_ExtendedMap, Load_Section, Search_Peaks, Statistic_Map

    !---- List of private functions ----!
    private :: Peak_Position

    !---- List of private subroutines ----!
    private :: Peak_List

    !---- Definitions ----!

    !!----
    !!---- ERR_MAPS
    !!----    logical, public  :: err_maps
    !!----
    !!----    Logical Variable indicating an error in MAPS_CALCULATIONS module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: err_maps

    !!----
    !!---- ERR_MESS_MAPS
    !!----    character(len=150), public :: err_mess_maps
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: err_mess_maps

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!--++
    !!--++ Function Peak_Position(Rho,i,j,k) Result(Pto)
    !!--++    integer,dimension(:,:,:), intent(in) :: Rho     ! Density map scaled as integer values
    !!--++    integer,                  intent(in) :: i
    !!--++    integer,                  intent(in) :: j       ! (i,j,k) is the central point
    !!--++    integer,                  intent(in) :: k
    !!--++
    !!--++    (Private)
    !!--++    Return the position of the peak
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Peak_Position(nr3d,i,j,k) Result(Pto)
       !---- Arguments ----!
       integer, dimension(:,:,:),intent(in) :: nr3d
       integer,                  intent(in) :: i
       integer,                  intent(in) :: j
       integer,                  intent(in) :: k
       real(kind=sp), dimension(3)          :: pto

       !---- Local variables ----!
       integer       :: ntx,nty,ntz
       integer       :: i1,i2,i3,j1,j2,j3,k1,k2,k3
       real(kind=sp) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19
       real(kind=sp) :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19
       real(kind=sp) :: b,c,d,h,kk,l,e,f,g,det,deltax,deltay,deltaz,x,y,z,dx,dy,dz

       !---- Calculation of the peak position ----!
       ntx=size(nr3d,1)
       nty=size(nr3d,2)
       ntz=size(nr3d,3)
       dx=1.0/real(ntx)
       dy=1.0/real(nty)
       dz=1.0/real(ntz)
       pto=0.0

       i2=i
       i1=i-1
       i3=i+1
       if (i1 <= 0)  i1=ntx+i1
       if (i3 > ntx) i3=i3-ntx

       j2=j
       j1=j-1
       j3=j+1
       if (j1 <= 0) j1=nty+j1
       if (j3 > nty) j3=j3-nty

       k2=k
       k1=k-1
       k3=k+1
       if (k1 <= 0) k1=ntz+k1
       if (k3 > ntz) k3=k3-ntz

        a1=real( MAX( nr3d(i2,j2,k2),1 ) )
        a2=real( MAX( nr3d(i1,j2,k2),1 ) )
        a3=real( MAX( nr3d(i3,j2,k2),1 ) )
        a4=real( MAX( nr3d(i3,j1,k2),1 ) )
        a5=real( MAX( nr3d(i2,j1,k2),1 ) )
        a6=real( MAX( nr3d(i1,j1,k2),1 ) )
        a7=real( MAX( nr3d(i1,j3,k2),1 ) )
        a8=real( MAX( nr3d(i2,j3,k2),1 ) )
        a9=real( MAX( nr3d(i3,j3,k2),1 ) )
       a10=real( MAX( nr3d(i2,j2,k1),1 ) )
       a11=real( MAX( nr3d(i1,j2,k1),1 ) )
       a12=real( MAX( nr3d(i3,j2,k1),1 ) )
       a13=real( MAX( nr3d(i2,j1,k1),1 ) )
       a14=real( MAX( nr3d(i2,j3,k1),1 ) )
       a15=real( MAX( nr3d(i2,j2,k3),1 ) )
       a16=real( MAX( nr3d(i1,j2,k3),1 ) )
       a17=real( MAX( nr3d(i3,j2,k3),1 ) )
       a18=real( MAX( nr3d(i2,j1,k3),1 ) )
       a19=real( MAX( nr3d(i2,j3,k3),1 ) )

        b1=LOG(a1)
        b2=LOG(a2)
        b3=LOG(a3)
        b4=LOG(a4)
        b5=LOG(a5)
        b6=LOG(a6)
        b7=LOG(a7)
        b8=LOG(a8)
        b9=LOG(a9)
       b10=LOG(a10)
       b11=LOG(a11)
       b12=LOG(a12)
       b13=LOG(a13)
       b14=LOG(a14)
       b15=LOG(a15)
       b16=LOG(a16)
       b17=LOG(a17)
       b18=LOG(a18)
       b19=LOG(a19)

       b=(b3+b4+b9+b12+b17-b2-b7-b6-b11-b16)/10.0
       c=(b7+b8+b9+b14+b19-b4-b5-b6-b13-b18)/10.0
       d=(b15+b16+b17+b18+b19-b10-b11-b13-b12-b14)/10.0
       h=(b13+b19-b18-b14)/4.0
       kk=(b11+b17-b16-b12)/4.0
       l=(b6+b9-b4-b7)/4.0
       e=(-10*b1+b2+b3+5*(b4+b6+b7+b9+b11+b12+b16+b17) &
                                   -6*(b5+b8+b10+b15)-2*(b13+b14+b18+b19))/21.0
       f=(-10*b1+b5+b8-6*(b2+b3+b10+b15)+5*(b4+b6+b7+b9+b13+b14+b18+b19) &
                                                     -2*(b11+b12+b16+b17))/21.0
       g=(-10*b1+b10+b15-6*(b2+b3+b5+b8)+5*(b11+b12+b13+b14+b16+b17+b18+b19) &
                                                         -2*(b4+b6+b7+b9))/21.0

       det=e*f*g+2*h*kk*l-kk*kk*f-h*h*e-l*l*g

       deltax=(-b*g*f-h*l*d-h*c*kk+kk*f*d+h*h*b+c*l*g)/det
       if (abs(deltax)-1.0 <= 0.0) then
          deltay=(-e*c*g-b*h*kk-l*d*kk+kk*kk*c+d*h*e+b*l*g)/det
          if (abs(deltay)-1.0 <= 0.0) then
             deltaz=(-e*f*d-l*c*kk-b*l*h+b*f*kk+l*l*d+h*c*e)/det
             if (abs(deltaz)-1.0 <= 0.0) then
                deltax=deltax*dx
                deltay=deltay*dy
                deltaz=deltaz*dz
             else
                deltax=0.0
                deltay=0.0
                deltaz=0.0
             end if
          else
            deltax=0.0
             deltay=0.0
             deltaz=0.0
          end if
       else
          deltax=0.0
          deltay=0.0
          deltaz=0.0
       end if

       x=(i-1)*dx
       y=(j-1)*dy
       z=(k-1)*dz

       x=x+deltax
       y=y+deltay
       z=z+deltaz

       x=mod(x+10.0,1.0)
       y=mod(y+10.0,1.0)
       z=mod(z+10.0,1.0)

       if (abs(x) <= 0.001) x=0.0
       if (abs(y) <= 0.001) y=0.0
       if (abs(z) <= 0.001) z=0.0
       if (abs(1.0-x) <= 0.001 ) x=0.0
       if (abs(1.0-y) <= 0.001 ) y=0.0
       if (abs(1.0-z) <= 0.001 ) z=0.0

       pto(1)=x
       pto(2)=y
       pto(3)=z

       return
    End Function Peak_Position


    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Maps()
    !!----
    !!----    Initialize the errors flags in Maps_Calculations
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Maps()

       err_maps=.false.
       err_mess_maps=" "

       return
    End Subroutine Init_Err_Maps

    !!----
    !!---- Subroutine Load_ExtentedMap(Rho,Ngrid,Limits,Rhonew)
    !!----    real, dimension(:,:,:), intent(in) :: rho
    !!----    integer, dimension(3),  intent(in) :: ngrid
    !!----    real,dimension(2,3),    intent(in) :: limits
    !!----    real, dimension(:,:,:), intent(out):: rhonew
    !!----
    !!----    Rhonew has one dimension more in each dimension that Rho
    !!----    This routine is useful for 2D representation.
    !!----        Rho(nx,ny,nz) -> Rhonew(nx+1,ny+1,nz+1)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Load_ExtendedMap(Rho,Ngrid,Limits,Rhonew)
       !---- Arguments ----!
       real, dimension(:,:,:), intent(in) :: rho
       integer, dimension(3),  intent(in) :: ngrid
       real,dimension(2,3),    intent(in) :: limits
       real, dimension(:,:,:), intent(out):: rhonew

       !---- Local Variables ----!
       integer            :: nx,ny,nz,nx1,ny1,nz1
       integer            :: i,j,k !,ii,jj,kk
       integer            :: ii1,ii2,jj1,jj2,kk1,kk2
       real               :: dx, dy, dz !, xinc, yinc, zinc
       real               :: xval,yval,zval
       real               :: v1,v2,v3,v4,v5,v6,v7,v8,r,s,t,valor
       real               :: x1,y1,z1,xx1,xx2,yy1,yy2,zz1,zz2
       real, dimension(3) :: rinc
       real, parameter    :: eps=0.0001

       nx=ngrid(1)
       ny=ngrid(2)
       nz=ngrid(3)
       nx1=nx+1
       ny1=ny+1
       nz1=nz+1

       dx=1.0/real(nx)
       dy=1.0/real(ny)
       dz=1.0/real(nz)

       rinc(1)=abs(limits(2,1)-limits(1,1))*dx
       rinc(2)=abs(limits(2,2)-limits(1,2))*dy
       rinc(3)=abs(limits(2,3)-limits(1,3))*dz

       !---- Loading RhoNew from Rho ----!
       if (abs(limits(1,1)-0.0) <= eps .and. abs(limits(2,1)-1.0) <= eps .and. &
           abs(limits(1,2)-0.0) <= eps .and. abs(limits(2,2)-1.0) <= eps .and. &
           abs(limits(1,3)-0.0) <= eps .and. abs(limits(2,3)-1.0) <= eps ) then

           do i=1,nx1
               do j=1,ny1
                  do k=1,nz1

                     ii1=mod(i,nx1)
                     if (ii1 == 0) ii1=1
                     jj1=mod(j,ny1)
                     if (jj1 == 0) jj1=1
                     kk1=mod(k,nz1)
                     if (kk1 == 0) kk1=1

                     rhonew(i,j,k)=rho(ii1,jj1,kk1)

                  end do
               end do
            end do

       else

          do i=1,nx1
             do j=1,ny1
                do k=1,nz1
                   !---- Pto equivalente en (0,1) ----!
                   xval=limits(1,1) + (i-1)*rinc(1)
                   yval=limits(1,2) + (j-1)*rinc(2)
                   zval=limits(1,3) + (k-1)*rinc(3)

                   if (abs(xval) <= eps) then
                      xval=0.0
                   elseif (abs(xval-1.0) <= eps) then
                      xval=1.0
                   else
                      if (xval < 0.0) xval=xval+1.0
                      if (xval > 1.0) xval=xval-1.0
                   end if

                   if (abs(yval) <= eps) then
                      yval=0.0
                   elseif (abs(yval-1.0) <= eps) then
                      yval=1.0
                   else
                      if (yval < 0.0) yval=yval+1.0
                      if (yval > 1.0) yval=yval-1.0
                   end if

                   if (abs(zval) <= eps) then
                      zval=0.0
                   elseif (abs(zval-1.0) <= eps) then
                      zval=1.0
                   else
                      if (zval < 0.0) zval=zval+1.0
                      if (zval > 1.0) zval=zval-1.0
                   end if

                   !---- Entre que planos de X esta el Pto ----!
                   do ii1=1,nx
                      ii2=ii1+1
                      xx1=(ii1-1)*dx
                      xx2=ii1*dx
                      if (abs(xx1-xval) <= eps) then
                         exit
                      elseif (xval > xx1 .and. xval <= xx2) then
                         exit
                      end if
                   end do

                   !---- Entre que planos de Y esta el Pto ----!
                   do jj1=1,ny
                      jj2=jj1+1
                      yy1=(jj1-1)*dy
                      yy2=jj1*dy
                      if (abs(yy1-yval) <= eps) then
                         exit
                      elseif (yval > yy1 .and. yval <= yy2) then
                         exit
                      end if
                   end do

                   !---- Entre que planos de Z esta el Pto ----!
                   do kk1=1,nz
                      kk2=kk1+1
                      zz1=(kk1-1)*dz
                      zz2=kk1*dz
                      if (abs(zz1-zval) <= eps) then
                         exit
                      elseif (zval > zz1 .and. zval <= zz2) then
                         exit
                      end if
                   end do

                   ii1=mod(ii1,nx+1)
                   if (ii1 == 0) ii1=1
                   ii2=mod(ii2,nx+1)
                   if (ii2 == 0) ii2=1

                   jj1=mod(jj1,ny+1)
                   if (jj1 == 0) jj1=1
                   jj2=mod(jj2,ny+1)
                   if (jj2 == 0) jj2=1

                   kk1=mod(kk1,nz+1)
                   if (kk1 == 0) kk1=1
                   kk2=mod(kk2,nz+1)
                   if (kk2 == 0) kk2=1

                   v1=rho(ii1,jj1,kk1)
                   v2=rho(ii2,jj1,kk1)
                   v3=rho(ii2,jj2,kk1)
                   v4=rho(ii1,jj2,kk1)
                   v5=rho(ii1,jj1,kk2)
                   v6=rho(ii2,jj1,kk2)
                   v7=rho(ii2,jj2,kk2)
                   v8=rho(ii1,jj2,kk2)

                   !---- Defino puntos del cubo ----!
                   x1=(ii1-1)*dx
                   y1=(jj1-1)*dy
                   z1=(kk1-1)*dz

                   r=(xval-x1)/dx
                   s=(yval-y1)/dy
                   t=(zval-z1)/dz

                   if (abs(r) <= eps) r=0.0
                   if (abs(s) <= eps) s=0.0
                   if (abs(t) <= eps) t=0.0

                   !---- Interpolacion tridimensional a 8 puntos ----!
                   valor=(1.0-r)*(1.0-s)*(1.0-t)* v1 + &
                              r *(1.0-s)*(1.0-t)* v2 + &
                              r *     s *(1.0-t)* v3 + &
                         (1.0-r)*     s *(1.0-t)* v4 + &
                         (1.0-r)*(1.0-s)*     t * v5 + &
                              r *(1.0-s)*     t * v6 + &
                              r *     s *     t * v7 + &
                         (1.0-r)*     s *     t * v8

                   rhonew(i,j,k)=valor

                end do
             end do
          end do

       end if

       return
    End Subroutine Load_ExtendedMap

    !!----
    !!---- SUBROUTINE Load_Section(Rho,ngrid,imap,section,limits,ngrid2,dmap)
    !!----    real, dimension(:,:,:), intent(in) :: rho
    !!----    integer, dimension(3),  intent(in) :: ngrid
    !!----    integer,                intent(in) :: imap
    !!----    integer,                intent(in) :: section
    !!----    real, dimension(2,2),   intent(in) :: limits
    !!----    integer, dimension(2),  intent(in) :: ngrid2
    !!----    real, dimension(:,:),   intent(out):: dmap
    !!----
    !!----    This routine only works with fractional coordinates
    !!----
    !!---- Updated: February - 2005
    !!
    Subroutine Load_Section(Rho,ngrid,imap,section,limits,ngrid2,dmap)
       !---- Arguments ----!
       real, dimension(:,:,:), intent(in) :: rho
       integer, dimension(3),  intent(in) :: ngrid
       integer,                intent(in) :: imap
       integer,                intent(in) :: section
       real, dimension(2,2),   intent(in) :: limits
       integer, dimension(2),  intent(in) :: ngrid2
       real, dimension(:,:),   intent(out):: dmap

       !---- Local Variables ----!
       integer :: ndimx, ndimy
       integer :: i,j,ii1,ii2,jj1,jj2,kk

       real            :: uinc,vinc
       real            :: uu,vv,x1,y1
       real            :: dx,dy,f1,f2,f3,z1,z2,z3,z4
       real            :: xx1,xx2,yy1,yy2
       real, parameter :: eps=0.0001

       !---- Contorno ----!
       dmap = 0.0

       select case (imap)
          case (1)
             ndimx=ngrid(1)+1
             ndimy=ngrid(2)+1
             dx=1.0/real(ngrid(1))
             dy=1.0/real(ngrid(2))

             kk=mod(section,ngrid(3)+1)
             if (kk == 0) kk=1

          case (2)
             ndimx=ngrid(2)+1
             ndimy=ngrid(3)+1
             dx=1.0/real(ngrid(2))
             dy=1.0/real(ngrid(3))

             kk=mod(section,ngrid(1)+1)
             if (kk == 0) kk=1

          case (3)
             ndimx=ngrid(3)+1
             ndimy=ngrid(1)+1
             dx=1.0/real(ngrid(3))
             dy=1.0/real(ngrid(1))

             kk=mod(section,ngrid(2)+1)
             if (kk == 0) kk=1
       end select

       uinc=(limits(1,2) - limits(1,1))/real(ngrid2(1)-1)
       vinc=(limits(2,2) - limits(2,1))/real(ngrid2(2)-1)

       do i=1,ngrid2(1)
          uu=limits(1,1)+uinc*(i-1)

          do j=1,ngrid2(2)
             vv=limits(2,1)+vinc*(j-1)

             x1=mod(uu+10.0,1.0)
             y1=mod(vv+10.0,1.0)

             !---- Entre que nodos X ----!
             do ii1=1,ndimx-1
                xx1=(ii1-1)*dx
                xx2=xx1+dx
                if (abs(xx1-x1) <= 0.0001) then
                   exit
                else if (x1 > xx1 .and. x1 <= xx2) then
                     exit
                end if
             end do
             ii2=ii1+1
             if (ii2 == ndimx) ii2=1

             !---- Entre que nodos Y ----!
             do jj1=1,ndimy-1
                yy1=(jj1-1)*dy
                yy2=yy1+dy
                if (abs(yy1-y1) <= 0.0001) then
                   exit
                else if (y1 > yy1 .and. y1 <= yy2) then
                   exit
                end if
             end do
             jj2=jj1+1
             if (jj2 == ndimy) jj2=1

             select case (imap)
                case (1)
                   z1=rho(ii1,jj1,kk)
                   z2=rho(ii2,jj1,kk)
                   z3=rho(ii1,jj2,kk)
                   z4=rho(ii2,jj2,kk)
                case (2)
                   z1=rho(kk,ii1,jj1)
                   z2=rho(kk,ii2,jj1)
                   z3=rho(kk,ii1,jj2)
                   z4=rho(kk,ii2,jj2)
                case (3)
                   z1=rho(jj1,kk,ii1)
                   z2=rho(jj1,kk,ii2)
                   z3=rho(jj2,kk,ii1)
                   z4=rho(jj2,kk,ii2)
             end select

             f1=( (x1-xx1)/(xx2-xx1) )*(z2-z1) + z1
             f2=( (y1-yy1)/(yy2-yy1) )*(z3-z1)
             f3=( (x1-xx1)/(xx2-xx1) )*( (y1-yy1)/(yy2-yy1) )*(z1-z2-z3+z4)

             dmap(i,j)=dmap(i,j) + f1+f2+f3

          end do
       end do

       return
    End Subroutine Load_Section

    !!--++
    !!--++ Subroutine Peak_List(Pto,Grp,Cell,Npeaks,Peaks)
    !!--++    real(kind=sp), dimension(4),   intent(in)     :: Pto      ! New position to add on the List
    !!--++    type(space_group_type),        intent(in)     :: Grp      ! Space Group
    !!--++    type(crystal_cell_type),       intent(in)     :: Cell     ! Cell
    !!--++    integer,                       intent(in out) :: NPeaks   ! Number of peaks on the list
    !!--++    real(kind=sp), dimension(:,:), intent(in out) :: Peaks    ! List of Peaks
    !!--++
    !!--++    (Private)
    !!--++    Add a new peak position on the list if there aren't
    !!--++    a closed peak (< 0.25).
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Peak_List(Pto,Grp,Cell,Npks,Pks)
       !---- Arguments ----!
       real(kind=sp), dimension(4),   intent(in)     :: Pto
       type(space_group_type),        intent(in)     :: Grp
       type(crystal_cell_type),       intent(in)     :: Cell
       integer,                       intent(in out) :: NPks
       real(kind=sp), dimension(:,:), intent(in out) :: Pks

       !---- Local variables ----!
       integer                      :: i,j
       real(kind=sp), parameter     :: eps=0.002
       real(kind=sp), dimension (3) :: pto1,pto2

       !---- Pto in asymmetric unit ----!
       do i=1,grp%multip
          pto1=ApplySO(grp%Symop(i),pto(1:3))
          pto1=mod(pto1+10.0,1.0)
          if (pto1(1) < grp%R_Asym_Unit(1,1)     .or. &
              pto1(1) > grp%R_Asym_Unit(1,2)+eps .or. &
              pto1(2) < grp%R_Asym_Unit(2,1)     .or. &
              pto1(2) > grp%R_Asym_Unit(2,2)+eps .or. &
              pto1(3) < grp%R_Asym_Unit(3,1)     .or. &
              pto1(3) > grp%R_Asym_Unit(3,2)+eps ) cycle
          exit
       end do

       if (npks == 0) then
          !---- First Peak ---!
          npks=1
          pks(1:3,1)=pto1
          pks(4,1)=pto(4)
       else
          !---- Searching if the peak in in the list ----!
          do j=1,npks
             do i=1,grp%multip
                pto2=ApplySO(grp%Symop(i),pks(1:3,j))
                pto2=mod(pto2+10.0,1.0)
                if (distance(pto1,pto2,cell) <= 0.25) return
             end do
          end do
          npks=npks+1
          pks(1:3,npks)=pto1
          pks(4,npks)=pto(4)
       end if

       return
    End Subroutine Peak_List

    !!----
    !!---- Subroutine Search_Peaks(Rho,Grp,Cell,Npeaks_to_Found,Peaks,Abs_Code)
    !!----    real(kind=sp), dimension(:,:,:),    intent(in)      :: Rho         ! The Map
    !!----    type(space_group_type),             intent(in)      :: Grp         ! Space Group
    !!----    type(crystal_cell_type),            intent(in)      :: Celda       ! Cell
    !!----    integer,                            intent(in out)  :: NPFound     ! Number of peaks to found
    !!----    real(kind=sp), dimension(4,NPfound),intent(out)     :: Peaks       ! Peak List
    !!----    logical, optional,                  intent(in)      :: Abs_Code    ! logical to use absolute value on Rho
    !!----
    !!----    General procedure to search peaks on Rho
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Search_Peaks(Rho,Grp,Cell,NPFound,Pks,Abs_Code)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:,:),    intent(in)      :: Rho
       type(space_group_type),             intent(in)      :: Grp
       type(crystal_cell_type),            intent(in)      :: Cell
       integer,                            intent(in out)  :: NPFound
       real(kind=sp), dimension(:,:),      intent(out)     :: Pks
       logical, optional,                  intent(in)      :: Abs_Code

       !---- Local Variables ----!
       logical                                    :: mode_abs
       integer                                    :: nscan,npks,ierror
       integer                                    :: nu,nv,nw
       integer                                    :: ii_ini,jj_ini,kk_ini,ii_fin,jj_fin,kk_fin
       integer                                    :: ia,ja,ka
       integer                                    :: i,j,k,ji,ji11,ji21,i11,i21
       integer                                    :: izt1t,izt2t,izt3t
       integer, allocatable,dimension (:,:,:)     :: nr3d
       real(kind=sp), parameter                   :: fac_max=1.0e3
       real(kind=sp)                              :: denabsmax,fac_scale
       real(kind=sp), dimension(4)                :: pto

       !---- Init ----!
       call init_err_maps()
       mode_abs=.false.
       if (present(abs_code)) mode_abs=abs_code
       npks=0
       pks =0.0

       !---- Memory for NR3D ----!
       nu=size(rho,1)
       nv=size(rho,2)
       nw=size(rho,3)
       if (allocated(nr3d)) deallocate(nr3d)
       allocate (nr3d(nu,nv,nw),stat=ierror)
       if (ierror /= 0) then
          err_maps=.true.
          err_mess_maps="Problems allocating memory for Nr3D"
          return
       end if
       nr3d=0

       !---- Loading Rho on Nr3D ----!
       denabsmax=max(maxval(rho),abs(minval(rho)))
       fac_scale=fac_max/denabsmax
       nr3d=nint(rho*fac_scale)
       if (mode_abs) nr3d=abs(nr3d)

       !---- Searching Zone ----!
       if (grp%R_Asym_Unit(2,1) < 0.0 ) then
          ii_ini=nint(grp%R_Asym_Unit(2,1)*nv)
       else
          ii_ini=nint(grp%R_Asym_Unit(2,1)*nv) + 1
       end if
       if (grp%R_Asym_Unit(1,1) < 0.0 ) then
          jj_ini=nint(grp%R_Asym_Unit(1,1)*nu)
       else
          jj_ini=nint(grp%R_Asym_Unit(1,1)*nu) + 1
       end if
       if (grp%R_Asym_Unit(3,1) < 0.0 ) then
          kk_ini=nint(grp%R_Asym_Unit(3,1)*nw)
       else
          kk_ini=nint(grp%R_Asym_Unit(3,1)*nw) + 1
       end if
       ii_fin=nint(grp%R_Asym_Unit(2,2)*nv) + 1
       jj_fin=nint(grp%R_Asym_Unit(1,2)*nu) + 1
       kk_fin=nint(grp%R_Asym_Unit(3,2)*nw) + 1
       ii_fin=min(ii_fin,nv)
       jj_fin=min(jj_fin,nu)
       kk_fin=min(kk_fin,nw)

       !---- Searching Procedure ----!
       search:do nscan=nint(fac_max),0,-10
          do ka=kk_ini,kk_fin
             if (ka <= 0) then
                k=nw+ka
             elseif (ka > nw) then
                k=ka-nw
             else
                k=ka
             end if
             do ia=ii_ini,ii_fin
                if (ia <= 0) then
                   i=nv+ia
                elseif (ia > nv) then
                   i=ia-nv
                else
                   i=ia
                end if
                do ja=jj_ini,jj_fin
                   if (ja <= 0) then
                      j=nu+ja
                   elseif(ja > nu) then
                      j=ja-nu
                   else
                      j=ja
                   end if

                   if (nr3d(j,i,k)-nscan <= 0) cycle
                   ji=j
                   ji11=j-1
                   ji21=j+1
                   i11=i-1
                   i21=i+1
                   izt2t=k
                   izt1t=k-1
                   izt3t=k+1

                   if (ji11  <= 0) ji11=nu+ji11
                   if (i11   <= 0) i11=nv+i11
                   if (izt1t <= 0) izt1t=nw+izt1t

                   if (ji21  > nu) ji21=ji21-nu
                   if (i21   > nv) i21=i21-nv
                   if (izt3t > nw) izt3t=izt3t-nw

                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt2t) <= 0) cycle   !Punto 122 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt2t) <  0) cycle   !      322 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,i11,izt2t) <= 0) cycle   !      312 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt2t) <= 0) cycle   !      212 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,i11,izt2t) <= 0) cycle   !      112 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,i21,izt2t) <  0) cycle   !      132 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt2t) <  0) cycle   !      232 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,i21,izt2t) <  0) cycle   !      332 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,  i,izt1t) <= 0) cycle   !      221 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt1t) <= 0) cycle   !      121 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt1t) <= 0) cycle   !      321 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt1t) <= 0) cycle   !      211 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt1t) <= 0) cycle   !      231 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,  i,izt3t) <  0) cycle   !      223 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt3t) <  0) cycle   !      123 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt3t) <  0) cycle   !      323 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt3t) <  0) cycle   !      213 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt3t) <  0) cycle   !      233 >=

                   !---- Position of the peak ----!
                   pto(1:3)=peak_position(nr3d,j,i,k)
                   pto(4)=rho(j,i,k)

                   !---- Good peak? ----!
                   call peak_list(pto,Grp,Cell,npks,pks)
                   if (npks == npfound) exit search
                end do
             end do
          end do
       end do search! nscan

       if (allocated(nr3d)) deallocate(nr3d)
       npfound=npks

       return
    End Subroutine Search_Peaks

    !!----
    !!---- Subroutine Statistic_Map(Rho,MaxV,MinV,AveV,SigmaV)
    !!----    real(kind=sp), dimension(:,:,:), intent(in) :: Rho
    !!----    real(kind=sp),                   intent(out):: MaxV      ! Maximum value of Rho
    !!----    real(kind=sp),                   intent(out):: MinV      ! Minimum value of Rho
    !!----    real(kind=sp),                   intent(out):: AveV      ! Average value of Rho
    !!----    real(kind=sp),                   intent(out):: SigmaV    ! Sigma value of Rho
    !!----
    !!----    Some statistic parameters of the map
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Statistic_Map(Rho,MaxV,MinV,AveV,SigmaV)
       !---- Arguments ----!
       real(kind=sp), dimension(:,:,:), intent(in) :: Rho
       real(kind=sp),                   intent(out):: MaxV
       real(kind=sp),                   intent(out):: MinV
       real(kind=sp),                   intent(out):: AveV
       real(kind=sp),                   intent(out):: SigmaV

       !---- Local Variables ----!
       integer :: i,j,k
       integer :: nu,nv,nw

       call init_err_maps()

       nu=size(rho,1)
       nv=size(rho,2)
       nw=size(rho,3)
       if (nu*nv*nw == 0) then
          err_maps=.true.
          err_mess_maps="Any dimension on Rho is zero"
          return
       end if

       MaxV  =maxval(rho)
       MinV  =minval(rho)
       AveV  = 0.0
       SigmaV= 0.0

       do i=1,nu
         do j=1,nv
             do k=1,nw
                avev=avev + rho(i,j,k)
                sigmav= sigmav + rho(i,j,k)*rho(i,j,k)
             end do
          end do
       end do
       avev  = avev/real(nu*nv*nw)
       sigmav= sigmav/real(nu*nv*nw) - avev*avev
       if (sigmav > 0.0001) then
          sigmav=sqrt(sigmav)
       else
          sigmav=0.0
       end if

       return
    End Subroutine Statistic_Map


 End Module Maps_Calculations