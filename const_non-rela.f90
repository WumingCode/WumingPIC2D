module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx     = 30001     ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny     = 2040      ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs   = 2         ! START POINT IN X
  integer, parameter :: nxge   = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs   = 2         ! START POINT IN Y
  integer, parameter :: nyge   = nygs+ny-1 ! END POINT
  integer, parameter :: np     = 200*nx    ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp    = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc  = 85        ! NUMBER OF PROCESSES
  integer, parameter :: nroot  = 0         ! ROOT PROC. NUMBER

!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs    = nxgs      ! START POINT IN X
  integer            :: nxe    = nxgs+(nx-1)*0.2 ! INITIAL SIZE

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 630000    !NUMBER OF ITERATION
  integer            :: it0    = 9999999   !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 90000     !INTERVAL FOR PARTICLES & FIELDS STORAGE          
  integer, parameter :: intvl2 = 1         !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3 = 20        !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  character(len=128) :: dir    = './pic2d/shock/run5@xc/'!DIRECTORY FOR OUTPUT
  character(len=128) :: file9  = 'init_param.dat' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: etlim  = 7.75*60*60 !MAX. ELAPSE TIME IN SEC.

!! OTHER CONSTANTS
  real(8)            :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 0.5D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl   = 1.0       !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*datan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!      mr : ION-TO-ELECTRON MASS RATIO
!!   alpha : wpe/wge = c/vth_e * sqrt(beta_e)
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!!      ma : ALFVEN MACH NUMBER ~ (SIGMA*GAMMA)^(-1/2) FOR V0~C
!!   theta : SHOCK ANGLE
!!   phi   : INCLINATION ANGLE FROM X-Y PLANE FOR BY & BZ
  integer, parameter :: n0     = 20
  real(8), parameter :: mr     = 225.0D0
  real(8), parameter :: alpha  = 10.0D0, beta = 0.5D0, rtemp=1.0D0
  real(8), parameter :: ma     = 30.0D0
  real(8), parameter :: theta  = 90.0D0/180.D0*pi 
  real(8), parameter :: phi    = 90.0D0/180.D0*pi

end module

