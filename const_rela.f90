module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx     = 3001      ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny     = 24        ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs   = 2         ! START POINT IN X
  integer, parameter :: nxge   = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs   = 2         ! START POINT IN Y
  integer, parameter :: nyge   = nygs+ny-1 ! END POINT
  integer, parameter :: np     = 40*nx     ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp    = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc  = 2         ! NUMBER OF PROCESSES
  integer, parameter :: nroot  = 0         ! ROOT PROC. NUMBER

!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs    = nxgs      ! START POINT IN X
  integer            :: nxe    = nxgs+(nx-1)*0.2 ! INITIAL SIZE

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 3000      !NUMBER OF ITERATION
  integer            :: it0    = 0         !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 3000      !INTERVAL FOR PARTICLES & FIELDS STORAGE          
  integer, parameter :: intvl2 = 1         !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3 = 1         !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  character(len=128) :: dir    = './pic2d/shock/test/'!DIRECTORY FOR OUTPUT
  character(len=128) :: file9  = 'init_param.dat' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: etlim  = 7.75*60*60 !MAX. ELAPSE TIME IN SEC.

!! OTHER CONSTANTS
  real(8)            :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl   = 1.0       !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*datan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!      mr : ION-TO-ELECTRON MASS RATIO
!!    gam0 : UPSTREAM BULK LORENTZ FACTOR
!!    sig0 : UPSTREAM SIGMA PARAMETER USING ION MASS
!!   alpha : wpe/wge = c/vth_e * sqrt(beta_e)
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!!      ma : ALFVEN MACH NUMBER
!!   theta : SHOCK ANGLE
!!   phi   : INCLINATION ANGLE FROM X-Y PLANE FOR BY & BZ
  integer, parameter :: n0     = 10
  real(8), parameter :: mr     = 100.0D0
  real(8), parameter :: gam0   = 100.0D0
  real(8), parameter :: sig0   = 1.0D-6 
  real(8), parameter :: alpha  = dsqrt(1.0D0/(gam0*sig0*mr)), beta = 0.5D0, rtemp=1.0D0
  real(8), parameter :: ma     = alpha*dsqrt(ma)*dsqrt(gam0*gam0-1.D0)/gam0
  real(8), parameter :: theta  = 90.0D0/180.D0*pi 
  real(8), parameter :: phi    = 0.0D0/180.D0*pi

end module

