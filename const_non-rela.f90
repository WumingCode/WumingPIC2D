module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx     = 16001     ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny     = 32      ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs   = 2         ! START POINT IN X
  integer, parameter :: nxge   = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs   = 2         ! START POINT IN Y
  integer, parameter :: nyge   = nygs+ny-1 ! END POINT
  integer, parameter :: ndim   = 5         ! DIMENSION OF PHASE SPACE--5:2D-3V, 6:2D-3V+ID FOR TRACKING
  integer, parameter :: np     = 100*nx    ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp    = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc  = 4       ! NUMBER OF PROCESSES
  integer, parameter :: nroot  = 0         ! ROOT PROC. NUMBER

!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs    = nxgs      ! START POINT IN X
  integer            :: nxe    = nxgs+(nx-1)*0.2 ! INITIAL SIZE

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 50   !NUMBER OF ITERATION
  integer            :: it0    = 0	!0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 50    !INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl2 = 30        !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3 = 160        !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  integer, parameter :: intvl4 = 10        !INTERVAL FOR RECORDING MOMENT DATA
  integer, parameter :: intvl5 = 1000        !INTERVAL FOR RECORDING ORBIT DATA
  character(len=128) :: dir    = './dat/'      !DIRECTORY FOR OUTPUT
  character(len=128) :: dir2   = './dat/'      !DIRECTORY FOR OUTPUT
  character(len=128) :: file9  = 'init_param.dat' !FILENAME OF INIT CONDITIONS
!  real(8), parameter :: etlim  = 23.5*60*60 !MAX. ELAPSE TIME IN SEC.
  real(8), parameter :: etlim  = 1.*60*60 !MAX. ELAPSE TIME IN SEC.

!! OTHER CONSTANTS
  real(8)            :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl   = 1.0D0/sqrt(2.0D0) !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!      mr : ION-TO-ELECTRON MASS RATIO
!!   alpha : wpe/wge = c/vth_e * sqrt(beta_e)
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!!      ma : ALFVEN MACH NUMBER ~ (SIGMA*GAMMA)^(-1/2) FOR V0~C
!!   theta : SHOCK ANGLE
!!     phi : INCLINATION ANGLE FROM X-Y PLANE FOR BY & BZ
!!    ldmp : FOR INITIAL VELOCITY PROFILE POSITION NEAR X=0 IN UNIT OF C/WPI
  integer, parameter :: n0    = 20
  real(8), parameter :: mr    = 400.0D0
  real(8), parameter :: alpha = 10.0D0, beta = 0.25D0, rtemp=1.0D0
  real(8), parameter :: ma    = 5.0D0*2.0D0/3.0D0 !MA IN THE SIM. FRAME
  real(8), parameter :: theta = 85.0D0/180.D0*pi 
  real(8), parameter :: phi   = 0.0D0/180.D0*pi
  real(8)            :: ldmp  = ma
!!TRACKING PARTICLES INITIALLY XRS <= X <= XRE ONLY WHEN NDIM=6
  real(8), parameter :: xrs   = 530.0D0 
  real(8), parameter :: xre   = 540.0D0

end module

