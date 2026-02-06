module variables
implicit none
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14,300)
!------------------------------------
! Options

#ifdef fdx
 real(dp), PARAMETER  :: dx = fdx
#else
 real(dp), PARAMETER  :: dx = 10.0
#endif

#ifdef fnVerticalLayers
 INTEGER, PARAMETER  :: nz = fnVerticalLayers
#else
 INTEGER, PARAMETER  :: nz = 80
#endif

#ifdef fnHorisontalColumns
 INTEGER, PARAMETER  :: ny = fnHorisontalColumns
#else
 INTEGER, PARAMETER  :: ny = 99
#endif

#ifdef fzIndexPlume
 INTEGER, PARAMETER  :: point_z = fzIndexPlume
#else
 INTEGER, PARAMETER  :: point_z = 12
#endif

#ifdef fyIndexPlume
 INTEGER, PARAMETER  :: point_y = fyIndexPlume
#else
 INTEGER, PARAMETER  :: point_y = 49
#endif

! duration of one mixing process = dt_chem
#ifdef ftotalRuntime
  INTEGER, PARAMETER  :: mixstep = ftotalRuntime
#else
  INTEGER, PARAMETER  :: mixstep = 30
#endif

#ifdef fmixingTimestep
real(dp), PARAMETER  :: dt_mete = fmixingTimestep
#else
real(dp), PARAMETER  :: dt_mete = 1
#endif

#ifdef fnSpec
 INTEGER, PARAMETER  :: NSPEC = fnSpec
#else
 INTEGER, PARAMETER  :: NSPEC = 2409
#endif

#ifdef fPath
 character(len=*), PARAMETER :: path = fPath
#else
 character(len=*), PARAMETER :: path = 'cells/'
#endif

#ifdef fuse_namelist
real(dp) :: KyKz = 2.0
INTEGER  :: randomize_k = 1
real(dp) :: friction_vel = 0.15
real(dp) :: pblh = 600.0
NAMELIST /NML_meteo/ KyKz,randomize_k,friction_vel,pblh
#else
  real(dp), PARAMETER  :: KyKz = fuse_KyKz
  INTEGER,  PARAMETER  :: randomize_k = fuse_randomize_k
  real(dp), PARAMETER  :: friction_vel = fuse_friction_vel
  real(dp), PARAMETER  :: pblh = fuse_pblh
#endif

!
! #ifdef f
! INTEGER, PARAMETER  ::  = f
! #else
! INTEGER, PARAMETER  ::  =
! #endif
!

! real(dp),PARAMETER  :: KyKz = 2.0 !, sigmaK = 1.0 ! horisontal:vertical diffusivity
! logical, PARAMETER  :: randomize_k = .true.
logical, PARAMETER  :: initW_Previous = .true.
! logical, PARAMETER  :: plume = .false.
real(dp),PARAMETER  :: point_multi = 1d4
! real(dp),PARAMETER  :: friction_vel = 0.15
! real(dp),PARAMETER  :: pblh = 600.0

!------------------------------------

character(len=128)  :: nmlfile ='empty'
character(len=2)    :: numb(ny)
character(len=1)    :: f_init
INTEGER             :: i_init=0,i_file=0, advance=0
INTEGER             :: iz,iy,k, ic,im,ioi
real(dp)            :: C(nz,ny)=1d0
real(dp)            :: Z(nz)=1d0, dz(2:nz)=1d0, Y(ny)=1d0, dy(2:ny)=1d0
real(dp)            :: daz(nz)=1d0, dbz(nz)=1d0, dcz(nz)=1d0, dfz(nz)=1d0
real(dp)            :: day(ny)=1d0, dby(ny)=1d0, dcy(ny)=1d0, dfy(ny)=1d0
real(dp)            :: Khtz(nz)=1d0, Khbz(nz)=1d0 !,kt(nz)=1d0
real(dp)            :: Khty(ny)=1d0, Khby(ny)=1d0 !,kt(nz)
real(dp)            :: termaz(nz)=1d0,termay(ny,nz)=1d0
real(dp)            :: termcz(nz)=1d0,termcy(ny,nz)=1d0
real(dp)            :: termbz(nz)=1d0,termby(ny,nz)=1d0
real(dp)            :: sigmaz0(nz)=1d0,sigmaz(nz)=1d0
real(dp)            :: sigmahor(nz,ny)=1d0
real(dp)            :: ch_gas(NSPEC)=1d0
real(dp)            :: Cmajor(NSPEC,nz,ny),Cmajor1(NSPEC,nz,ny)
real(dp)            :: random_z(nz)=0.5d0,random_y(ny)=0.5d0

contains
subroutine read_initfile
implicit none

#ifdef fuse_namelist
OPEN(UNIT=888, FILE=TRIM(nmlfile), STATUS='OLD', ACTION='READ', iostat=ioi)
if (ioi == 0) then
  READ(888,NML = NML_meteo, IOSTAT=ioi)
  close(888)
end if
#endif
end subroutine read_initfile


end module variables
