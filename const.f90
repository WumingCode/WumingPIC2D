module const

  implicit none
  integer, parameter :: nx    = 40001     ! number of grid points in x
  integer, parameter :: ny    = 8192      ! number of grid points in y
  integer, parameter :: nxgs  = 2         ! start point in x
  integer, parameter :: nxge  = nxgs+nx-1 ! end point
  integer, parameter :: nygs  = 2         ! start point in y
  integer, parameter :: nyge  = nygs+ny-1 ! end point
  integer, parameter :: np    = 200*nx    ! number of particles in each cell
  integer, parameter :: nsp   = 2         ! number of particle species
  integer, parameter :: nproc = 1024       ! number of processors

end module

