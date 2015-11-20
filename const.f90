module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx    = 3001     ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny    = 1224     ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs  = 2         ! START POINT IN X
  integer, parameter :: nxge  = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs  = 2         ! START POINT IN Y
  integer, parameter :: nyge  = nygs+ny-1 ! END POINT
  integer, parameter :: np    = 200*nx    ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp   = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc = 51        ! NUMBER OF PROCESSES
!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs   = nxgs
  integer            :: nxe   = nxgs+(nx-1)*0.2
!! itmax  : NUMBER OF ITERATION
!! it0    : BASE COUNT
!! intvl1 : INTERVAL FOR PARTICLES & FIELDS STORAGE 
!! intvl3 : INTERVAL FOR INJECTING PARTICLES
!! intvl4 : INTERVAL FOR UPDATING PHYSICAL REGION IN X
  integer, parameter :: itmax  = 63000 
  integer            :: it0    = 0                !0:INITAL, 9999999: RESTART DATA
  integer, parameter :: intvl1 = 9000
  integer, parameter :: intvl3 = 1
  integer, parameter :: intvl4 = 20
  character(len=128) :: dir    = './pic2d/shock/test/'!PARENT DIRECTORY FOR OUTPUT
  character(len=128) :: file9  = 'init_param.dat' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: pi     = 4.0d0*atan(1.0)
  real(8), parameter :: gfac   = 0.501d0          !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 0.5d0            !CFL CONDITION FOR LIGHT SPEED
  real(8)            :: delx   = 1.0d0            !CELL WIDTH
  real(8), parameter :: rdbl   = 1.0              !DEBYE LENGTH / CELL WIDTH

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!       c : SPEED OF LIGHT
!!      ma : ION-TO-ELECTRON MASS RATIO
!!   alpha : wpe/wge           
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!!      ma : ALFVEN MACH NUMBER
!!   theta : SHOCK ANGLE
  integer, parameter :: n0    = 20
  real(8)            :: c     = 1.0d0
  real(8), parameter :: mr    = 225.0d0 
  real(8), parameter :: alpha = 10.0d0, beta = 0.5d0, rtemp=1.0d0
  real(8), parameter :: ma    = 30.0d0
  real(8), parameter :: theta = 90.D0 /360.D0*2.*pi 

end module

