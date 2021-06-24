module app
  use iso_fortran_env, only: int64
  use mpi
  use wuming2d
  implicit none
  private

  ! main simulation loop
  public:: app__main


  ! configuration file and initial parameter files
  character(len=*), parameter   :: config_default = 'config.json'
  character(len=:), allocatable :: config_string
  character(len=:), allocatable :: config
  character(len=*), parameter   :: param = 'init_param'

  ! read from "config" section
  logical                       :: restart
  character(len=128)            :: restart_file
  character(len=:), allocatable :: datadir
  real(8)                       :: max_elapsed
  integer                       :: max_it
  integer                       :: intvl_ptcl
  integer                       :: intvl_mom
  integer                       :: intvl_orb
  integer                       :: intvl_expand

  ! read from "parameter" section
  integer :: num_process
  integer :: n_ppc
  integer :: n_x
  integer :: n_x_ini
  integer :: n_y
  real(8) :: u_inject
  real(8) :: mass_ratio
  real(8) :: sigma_e
  real(8) :: omega_pe
  real(8) :: v_the
  real(8) :: v_thi
  real(8) :: theta_bn
  real(8) :: phi_bn
  real(8) :: l_damp_ini

  integer :: nproc
  integer :: np
  integer :: n0
  integer :: nx, nxgs, nxge, nxs, nxe
  integer :: ny, nygs, nyge !, nys, nye
  integer :: mpierr

  !************************ NUMERICAL CONSTANTS ***********************************!
  !integer, parameter :: nx     = 1000      ! NUMBER OF GRID POINTS IN X
  !integer, parameter :: ny     = 100       ! NUMBER OF GRID POINTS IN Y
  !integer, parameter :: nxgs   = 2         ! START POINT IN X
  !integer, parameter :: nxge   = nxgs+nx-1 ! END POINT
  !integer, parameter :: nygs   = 2         ! START POINT IN Y
  !integer, parameter :: nyge   = nygs+ny-1 ! END POINT
  integer, parameter :: ndim   = 6         ! DIMENSION OF PHASE SPACE--5:2D-3V, 6:2D-3V+ID FOR TRACKING
  !integer, parameter :: np     = 6*nx      ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp    = 2         ! NUMBER OF PARTICLE SPECIES
  !integer, parameter :: nproc  = 4         ! NUMBER OF PROCESSES
  integer, parameter :: nroot  = 0         ! ROOT PROC. NUMBER

  ! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  !integer            :: nxs    = nxgs      ! START POINT IN X
  !integer            :: nxe    = nxgs+500 ! INITIAL SIZE

  ! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  !integer, parameter :: itmax  = 1000      !NUMBER OF ITERATION
  integer            :: it0    = 0         !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 100       !INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl2 = 1         !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3 = 1         !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  integer, parameter :: intvl4 = 10        !INTERVAL FOR RECORDING MOMENT DATA
  integer, parameter :: intvl5 = 5         !INTERVAL FOR RECORDING ORBIT DATA
  !character(len=128) :: dir    = './'      !DIRECTORY FOR OUTPUT
  !character(len=128) :: dir2   = './'      !DIRECTORY FOR OUTPUT
  !character(len=128) :: file9  = 'init_param' !FILENAME OF INIT CONDITIONS
  !real(8), parameter :: etlim  = 23.5*60*60 !MAX. ELAPSE TIME IN SEC.

  ! OTHER CONSTANTS
  real(8), parameter :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

  !************************ PHYSICAL CONSTANTS ***********************************!
  !      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
  !      mr : ION-TO-ELECTRON MASS RATIO
  !    gam0 : UPSTREAM BULK LORENTZ FACTOR
  !    sig0 : UPSTREAM SIGMA PARAMETER USING ION MASS
  !   rtemp : Te/Ti
  !     vt? : (ION/ELECTRON) THERMAL SPEED
  !   theta : SHOCK ANGLE
  !     phi : INCLINATION ANGLE FROM X-Y PLANE FOR BY & BZ
  !    ldmp : FOR INITIAL VELOCITY PROFILE POSITION NEAR X=0 IN UNIT OF C/WPI
  !integer, parameter :: n0     = 2
  !real(8), parameter :: mr     = 1.0D0
  !real(8), parameter :: gam0   = 40.0D0
  !real(8), parameter :: sig0   = 1.D-1
  !real(8), parameter :: rtemp  = 1.0D0
  !real(8), parameter :: vte    = 0.0D0
  !real(8), parameter :: vti    = sqrt(1.0D0/(rtemp*mr))*vte
  !real(8), parameter :: theta  = 90.0D0/180.D0*pi
  !real(8), parameter :: phi    = 90.0D0/180.D0*pi
  !real(8)            :: ldmp   = 10.0D0

  ! TRACKING PARTICLES INITIALLY XRS <= X <= XRE ONLY WHEN NDIM=6
  real(8), parameter :: xrs   = 450.0D0
  real(8), parameter :: xre   = 500.0D0

  !
  ! main variables
  !
  integer, allocatable, public :: np2(:,:), cumcnt(:,:,:)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  real(8), allocatable, public :: den(:,:,:),vel(:,:,:,:),temp(:,:,:,:)
  real(8), save                :: u0, v0, b0, gam0, delt

contains
  !
  ! main simulation loop
  !
  subroutine app__main()
    implicit none
    integer :: it, it_last
    real(8) :: etime, etime0

    ! initialization
    call load_config()
    call init()

    ! current clock
    etime0 = get_etime()

    ! main loop
    do it = it0, max_it
       ! update
       call particle__solv(gp, up, uf, cumcnt, nxs, nxe)
       call boundary__particle_injection(gp, np2, nxs, nxe, u0)
       call field__fdtd_i(uf, up, gp, cumcnt, nxs, nxe, &
            & boundary__dfield, &
            & boundary__curre, &
            & boundary__phi)
       call boundary__particle_y(gp, np2)
       call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)

       ! injection
       if(mod(it,intvl2) == 0) call inject(it)
       if(mod(it,intvl3) == 0) call relocate(it)

       ! output entire particles
       if ( mod(it, intvl_ptcl) == 0 ) then
          call io__ptcl(up, uf, np2, it)
       end if

       ! ouput tracer particles
       if ( mod(it, intvl_orb) == 0 ) then
          call io__orb(up, uf, np2, it)
       end if

       ! output moments and electromagnetic fields
       if ( mod(it, intvl_mom) == 0 ) then
          call mom_calc__accl(gp, up,uf, cumcnt, nxs, nxe)
          call mom_calc__nvt(den, vel, temp, gp, np2)
          call boundary__mom(den, vel, temp)
          call io__mom(den, vel, temp, uf, it)
       endif

       ! check elapsed time
       etime = get_etime() - etime0
       if ( etime >= max_elapsed ) then
          if ( nrank == nroot ) then
             write(0,*) '*** elapse time over ***', it, etime
          end if

          ! save state for restart
          write(restart_file, '(i7.7, "_restart")') it
          call save_restart(up, uf, np2, nxs, nxe, it, restart_file)
          call finalize()
          stop
       endif
    enddo

    ! save final state
    it = max_it + 1
    write(restart_file, '(i7.7, "_restart")') it
    call save_restart(up, uf, np2, nxs, nxe, it, restart_file)
    call finalize()

  end subroutine app__main

  !
  ! parse command line argument and load configuration file
  !
  subroutine load_config()
    implicit none

    logical :: status, found
    integer :: arg_count, npp
    character(len=:), allocatable :: filename

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    ! check ndim
    if( ndim /= 6 ) then
       write(0,*) 'Error: ndim must be 6'
       stop
    end if

    !
    ! process command line to find a configuration file
    !
    arg_count = command_argument_count()

    if( arg_count == 0 ) then
       config = config_default
    else
       ! only the first argument is relevant
       call get_command_argument(1, config)
    end if

    ! check init file
    inquire(file=trim(config), exist=status)

    if( .not. status ) then
       write(0, '("Error: ", a, " does not exists")') trim(config)
       stop
    end if

    ! try loading init file
    call json%initialize()
    call file%initialize()
    call file%load(trim(config))

    if( file%failed() ) then
       write(0, '("Error: failed to load ", a)') trim(config)
       stop
    end if

    ! read "config" section
    call file%get(root)
    call json%get(root, 'config', p)

    call json%get(p, 'datadir', datadir)
    call json%get(p, 'max_elapsed', max_elapsed)
    call json%get(p, 'max_it', max_it)
    call json%get(p, 'intvl_ptcl', intvl_ptcl)
    call json%get(p, 'intvl_mom', intvl_mom)
    call json%get(p, 'intvl_orb', intvl_orb)
    call json%get(p, 'intvl_expand', intvl_expand)

    ! make sure this is a directory
    datadir = trim(datadir) // '/'

    ! restart file
    call json%get(p, 'restart_file', filename, found)

    if ( found ) then
       restart = .true.
       restart_file = filename
    endif

    ! read "parameter" section and initialize
    call json%get(root, 'parameter', p)
    call json%get(p, 'num_process', num_process)
    call json%get(p, 'n_ppc', n_ppc)
    call json%get(p, 'n_x', n_x)
    call json%get(p, 'n_y', n_y)
    call json%get(p, 'n_x_ini', n_x_ini)
    call json%get(p, 'u_inject', u_inject)
    call json%get(p, 'mass_ratio', mass_ratio)
    call json%get(p, 'sigma_e', sigma_e)
    call json%get(p, 'omega_pe', omega_pe)
    call json%get(p, 'v_the', v_the)
    call json%get(p, 'v_thi', v_thi)
    call json%get(p, 'theta_bn', theta_bn)
    call json%get(p, 'phi_bn', phi_bn)
    call json%get(p, 'l_damp_ini', l_damp_ini)

    nproc = num_process
    nx    = n_x
    ny    = n_y
    np    = n_ppc * nx * 5
    n0    = n_ppc
    nxgs  = 2
    nxge  = nxgs + nx - 1
    nygs  = 2
    nyge  = nygs + ny - 1
    nxs   = nxgs
    nxe   = nxgs + n_x_ini

    ! degree to radian
    theta_bn = theta_bn * pi / 180
    phi_bn   = phi_bn   * pi / 180

    call json%serialize(root, config_string)
    call json%destroy()
    call file%destroy()

  end subroutine load_config

  !
  ! initialize simulation
  !
  subroutine init()
    implicit none
    integer :: n, isp, i, j, ndim_in
    real(8) :: r(nsp), q(nsp)
    real(8) :: wpe, wpi, wge, wgi, vte, vti

    ! MPI
    call mpi_set__init(nygs, nyge, nproc)

    ! random number
    call init_random_seed()

    ! allocate memory and initialize everything by zero
    allocate(np2(nys:nye,nsp))
    allocate(cumcnt(nxgs:nxge,nys:nye,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2))
    allocate(up(ndim,np,nys:nye,nsp))
    allocate(gp(ndim,np,nys:nye,nsp))
    allocate(den(nxgs-1:nxge+1,nys-1:nye+1,1:nsp))
    allocate(vel(nxgs-1:nxge+1,nys-1:nye+1,1:3,1:nsp))
    allocate(temp(nxgs-1:nxge+1,nys-1:nye+1,1:3,1:nsp))
    np2    = 0
    cumcnt = 0
    uf     = 0
    up     = 0
    gp     = 0
    den    = 0
    vel    = 0
    temp   = 0

    ! set physical parameters
    delt = cfl*delx/c
    u0   =-abs(u_inject)
    gam0 = sqrt(1 + (u0*u0)/(c*c))
    v0   = u0/gam0
    wpe  = omega_pe * sqrt(gam0)
    wge  = omega_pe * sqrt(sigma_e)
    wpi  = wpe / sqrt(mass_ratio)
    wgi  = wge / mass_ratio
    vte  = v_the
    vti  = v_thi
    r(1) = mass_ratio
    r(2) = 1.0d0
    q(1) =+sqrt(r(1) / (4*pi*n0)) * wpi
    q(2) =-sqrt(r(2) / (4*pi*n0)) * wpe
    b0   = gam0*r(1)*c / q(1) * wgi

    ! number of particles
    np2(nys:nye,1:nsp) = n0*(nxe-nxs)
    if ( nrank == nroot ) then
       if ( n0*(nxge-nxgs) > np ) then
          write(0,*) 'Error: Too large number of particles'
          stop
       endif
    endif

    ! preparation of sort
    do isp = 1, nsp
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j = nys, nye
          cumcnt(nxs,j,isp) = 0
          do i = nxs+1, nxe
             cumcnt(i,j,isp) = cumcnt(i-1,j,isp) + n0
          enddo
          if ( cumcnt(nxe,j,isp) /= np2(j,isp) ) then
             write(0,*) 'Error: invalid values encounterd for cumcnt'
             stop
          endif
       enddo
       !$OMP END PARALLEL DO
    enddo

    ! initialize modules
    call boundary__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & nup, ndown, mnpi, mnpr, ncomw, nerr, nstat, delx, delt, c)
    call particle__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & delx, delt, c, q, r)
    call field__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & mnpr, ncomw, opsum, nerr, delx, delt, c, q, r, gfac)
    call sort__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye)
    call io__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & nproc, nrank, delx, delt, c, q, r, datadir)
    call mom_calc__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & delx, delt, c, q, r)

    if ( restart ) then
       ! restart
       call io__input(up, uf, np2, nxs, nxe, it0, restart_file)
       call sort__bucket(up, up, cumcnt, np2, nxs, nxe)
    else
       ! output parameters and set initial condition
       call save_param(n0, wpe, wpi, wge, wgi, vti, vte, param)
       call set_initial_condition()
    endif

    ! copy
    gp = up

  end subroutine init

  !
  ! finalize simulation
  !
  subroutine finalize()
    implicit none

    call io__finalize()
    call MPI_Finalize(mpierr)

  end subroutine finalize

  !
  ! set initial condition for field and particles
  !
  subroutine set_initial_condition()
    implicit none
    integer :: i, j, ii, isp
    real(8) :: v1, gam1, gamp, sd(nsp)

    !
    ! electromagnetic field
    !
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys-2,nye+2
       do i=nxgs-2,nxge+2
          uf(1,i,j) = b0*cos(theta_bn)
          uf(2,i,j) = b0*sin(theta_bn)*cos(phi_bn)
          uf(3,i,j) = b0*sin(theta_bn)*sin(phi_bn)
          uf(4,i,j) = 0.0D0
          uf(5,i,j) =+vprofile(i*delx)*uf(3,i,j)/c
          uf(6,i,j) =-vprofile(i*delx)*uf(2,i,j)/c
       enddo
    enddo
    !$OMP END PARALLEL DO

    !
    ! particle position
    !
    isp = 1
    !$OMP PARALLEL DO PRIVATE(ii,j)
    do j=nys,nye
       do ii=1,np2(j,isp)
          up(1,ii,j,1) = nxs*delx+(nxe-nxs)*delx*ii/(np2(j,isp)+1)
          up(2,ii,j,1) = (j + uniform_rand()) * delx
          up(1,ii,j,2) = up(1,ii,j,1)
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !
    ! particle velocity
    !
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp
       !$OMP PARALLEL DO PRIVATE(ii,j,v1,gam1,gamp)
       do j=nys,nye
          do ii=1,np2(j,isp)
             ! Maxwellian in fluid rest frame
             up(3,ii,j,isp) = sd(isp) * normal_rand()
             up(4,ii,j,isp) = sd(isp) * normal_rand()
             up(5,ii,j,isp) = sd(isp) * normal_rand()

             ! Lorentz transform to lab frame
             v1   = vprofile(up(1,ii,j,isp))
             gam1 = 1/sqrt(1-(v1/c)**2)
             gamp = sqrt(1 + ( &
                  & up(3,ii,j,isp)**2 + &
                  & up(4,ii,j,isp)**2 + &
                  & up(5,ii,j,isp)**2 )/c**2)
             up(3,ii,j,isp) = gam1*(up(3,ii,j,isp) + v1*gamp)
          enddo
       enddo
       !$OMP END PARALLEL DO
    enddo

    ! set particle IDs
    call set_particle_ids()

  end subroutine set_initial_condition

  !
  ! set initial particle IDs
  !
  subroutine set_particle_ids()
    implicit none
    integer :: isp, i, j

    integer(8) :: gcumsum(nproc+1,nsp), lcumsum(nys:nye+1,nsp), pid

    if( ndim /= 6 ) then
       return
    end if

    ! calculate the first particle IDs
    call get_global_cumsum(np2, gcumsum)

    do isp = 1, nsp
       lcumsum(nys,isp) = gcumsum(nrank+1,isp)
       do j = nys, nye
          lcumsum(j+1,isp) = lcumsum(j,isp) + np2(j,isp)
       end do
    end do

    ! unique ID as 64bit integer (negative by default)
    do isp = 1, nsp
       !$OMP PARALLEL DO PRIVATE(i,j,pid)
       do j = nys, nye
          do i = 1, np2(j,isp)
             pid = lcumsum(j,isp) + i
             up(6,i,j,isp) = transfer(-pid, 1.0_8)
          end do
       end do
       !$OMP END PARALLEL DO
    end do

    ! make particle ID positive for output
    do isp = 1, nsp
       !$OMP PARALLEL DO PRIVATE(i,j,pid)
       do j = nys, nye
          do i = 1, np2(j,isp)
             if ( up(1,i,j,isp) >= xrs .and. up(1,i,j,isp) <= xre ) then
                pid = transfer(up(6,i,j,isp), 1_8)
                up(6,i,j,isp) = transfer(sign(pid, +1_8), 1.0_8)
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO
    enddo

  end subroutine set_particle_ids

  !
  ! save everything for restart
  !
  subroutine save_restart(up, uf, np2, nxs, nxe, it, restart_file)
    implicit none
    integer, intent(in)          :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in)          :: it
    real(8), intent(in)          :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    character(len=*), intent(in) :: restart_file

    logical :: found
    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    call json%initialize()
    call file%initialize()
    call file%deserialize(config_string)
    call file%get(root)
    call json%get(root, 'config', p)
    call json%update(p, 'restart_file', trim(restart_file), found)

    ! write data to the disk
    call io__output(up, uf, np2, nxs, nxe, it, trim(restart_file))

    if ( nrank == nroot ) then
       call json%print(root, config)
    end if

    call json%destroy()
    call file%destroy()

  end subroutine save_restart

  !
  ! save parameters
  !
  subroutine save_param(n0, wpe, wpi, wge, wgi, vti, vte, filename)
    implicit none
    integer, intent(in)          :: n0
    real(8), intent(in)          :: wpe, wpi, wge, wgi, vti, vte
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp
    integer :: fh

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    ! save deafult parameters
    call io__param(n0, wpe, wpi, wge, wgi, vti, vte, filename)

    ! save additional parameters
    datafile = trim(datadir) // trim(filename) // '.raw'
    jsonfile = trim(datadir) // trim(filename) // '.json'

    ! open json file
    call file%initialize()
    call json%initialize()
    call file%load(jsonfile)
    call file%get(root)

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'a')

    ! put attributes
    call json%get(root, 'attribute', p)

    call jsonio_put_attribute(json, p, u0, 'u0', disp, '')
    call mpiio_write_atomic(fh, disp, u0)

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine save_param


  subroutine relocate(it)
    implicit none
    integer, intent(in) :: it
    integer :: isp, j, ii, ii2 ,ii3, dn
    !real(8) :: aa, bb, cc, gamp, dx
    real(8) :: dx, v1, gam1, gamp, sd(nsp)

    integer(8) :: gcumsum(nproc+1,nsp), nptotal(nsp), pid

    if(nxe==nxge) return

    dx  = v0*delt*mod(it,intvl2)/delx
    dn  = int(n0*abs(dx)+0.5)
    nxe  = nxe+1

    ! get particle number for ID
    call get_global_cumsum(np2, gcumsum)
    nptotal(1:nsp) = gcumsum(nproc+1,1:nsp)

    !PARTICLE POSITION
    !$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          up(1,ii2,j,1) = (nxe-1.D0+dx)*delx-dx*delx*ii/(dn+1.D0)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          up(2,ii2,j,1) = (j + uniform_rand()) * delx
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
       do ii=1,n0
          ii2 = np2(j,1)+dn+ii
          ii3 = np2(j,2)+dn+ii

          up(1,ii2,j,1) = (nxe-1.)*delx+delx*ii/(n0+1.)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          up(2,ii2,j,1) = (j + uniform_rand()) * delx
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp

       !$OMP PARALLEL DO PRIVATE(ii,j,v1,gam1,gamp,pid)
       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn+n0
             ! Maxwellian in fluid rest frame
             up(3,ii,j,isp) = sd(isp) * normal_rand()
             up(4,ii,j,isp) = sd(isp) * normal_rand()
             up(5,ii,j,isp) = sd(isp) * normal_rand()

             ! Lorentz transform to lab frame
             v1   = vprofile(up(1,ii,j,isp))
             gam1 = 1/sqrt(1-(v1/c)**2)
             gamp = sqrt(1 + ( &
                  & up(3,ii,j,isp)**2 + &
                  & up(4,ii,j,isp)**2 + &
                  & up(5,ii,j,isp)**2 )/c**2)
             up(3,ii,j,isp) = gam1*(up(3,ii,j,isp) + v1*gamp)

             ! particle ID
             pid = ii - np2(j,isp) + (j-nygs)*(dn+n0) + nptotal(isp)
             up(6,ii,j,isp) = transfer(-pid, 1.0_8)
          enddo
          np2(j,isp) = np2(j,isp)+dn+n0
          cumcnt(nxe-1,j,isp) = cumcnt(nxe-1,j,isp)+dn
          cumcnt(nxe,j,isp) = cumcnt(nxe-1,j,isp)+n0
       enddo
       !$OMP END PARALLEL DO
    enddo

    !$OMP PARALLEL DO PRIVATE(j)
    do j=nys-2,nye+2
       uf(2,nxe-1,j) = b0*sin(theta_bn)*cos(phi_bn)
       uf(3,nxe-1,j) = b0*sin(theta_bn)*sin(phi_bn)
       uf(5,nxe-1,j) = v0*uf(3,nxe-1,j)/c
       uf(6,nxe-1,j) = -v0*uf(2,nxe-1,j)/c

       uf(2,nxe,j) = b0*sin(theta_bn)*cos(phi_bn)
       uf(3,nxe,j) = b0*sin(theta_bn)*sin(phi_bn)
    enddo
    !$OMP END PARALLEL DO

  end subroutine relocate


  subroutine inject(it)
    implicit none
    integer :: it, isp, ii, ii2, ii3, j, dn
    !real(8) :: aa, bb, cc, gamp, dx
    real(8) :: dx, v1, gam1, gamp, sd(nsp)

    integer(8) :: gcumsum(nproc+1,nsp), nptotal(nsp), pid

    !INJECT PARTICLES IN x=nxe-v0*dt~nxe*dt
    dx  = v0*delt*(it-max(int(it/intvl3)*intvl3,it-intvl2))/delx
    if(dx == 0.0D0) dx = v0*delt*min(intvl2,intvl3)/delx
    if(nxe == nxge) dx = v0*delt*intvl2/delx
    dn  = int(abs(n0*dx)+0.5)

    ! get particle number for ID
    call get_global_cumsum(np2, gcumsum)
    nptotal(1:nsp) = gcumsum(nproc+1,1:nsp)

    !$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          up(1,ii2,j,1) = nxe*delx+dx*(dn-ii+1)/(dn+1)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          up(2,ii2,j,1) = (j + uniform_rand()) * delx
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp
       !$OMP PARALLEL DO PRIVATE(ii,j,v1,gam1,gamp,pid)
       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn
             ! Maxwellian in fluid rest frame
             up(3,ii,j,isp) = sd(isp) * normal_rand()
             up(4,ii,j,isp) = sd(isp) * normal_rand()
             up(5,ii,j,isp) = sd(isp) * normal_rand()

             ! Lorentz transform to lab frame
             v1   = vprofile(up(1,ii,j,isp))
             gam1 = 1/sqrt(1-(v1/c)**2)
             gamp = sqrt(1 + ( &
                  & up(3,ii,j,isp)**2 + &
                  & up(4,ii,j,isp)**2 + &
                  & up(5,ii,j,isp)**2 )/c**2)
             up(3,ii,j,isp) = gam1*(up(3,ii,j,isp) + v1*gamp)

             ! particle ID
             pid = ii - np2(j,isp) + (j-nygs)*dn + nptotal(isp)
             up(6,ii,j,isp) = transfer(-pid, 1.0_8)
          enddo
       enddo
       !$OMP END PARALLEL DO
    enddo

    do isp=1,nsp
       !$OMP WORKSHARE
       np2(nys:nye,isp) = np2(nys:nye,isp)+dn
       cumcnt(nxe,nys:nye,isp) = cumcnt(nxe,nys:nye,isp)+dn
       !$OMP END WORKSHARE
    enddo

    !$OMP PARALLEL DO PRIVATE(j)
    do j=nys-2,nye+2
       uf(2,nxe-1,j) = b0*sin(theta_bn)*cos(phi_bn)
       uf(3,nxe-1,j) = b0*sin(theta_bn)*sin(phi_bn)
       uf(5,nxe-1,j) = v0*uf(3,nxe-1,j)/c
       uf(6,nxe-1,j) = -v0*uf(2,nxe-1,j)/c

       uf(2,nxe,j) = b0*sin(theta_bn)*cos(phi_bn)
       uf(3,nxe,j) = b0*sin(theta_bn)*sin(phi_bn)
    enddo
    !$OMP END PARALLEL DO

  end subroutine inject

  !
  ! get global cumulative sum of particle numbers
  !
  subroutine get_global_cumsum(np2, cumsum)
    implicit none
    integer, intent(in)       :: np2(nys:nye,nsp)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)

    integer :: i, isp, mpierr
    integer(8) :: lcount(nsp), gcount(nsp, nproc)

    ! get number of particles for each proces
    lcount(1:nsp) = sum(np2(nys:nye,1:nsp), dim=1)
    call MPI_Allgather(lcount, nsp, MPI_INTEGER8, gcount, nsp, MPI_INTEGER8, &
         & MPI_COMM_WORLD, mpierr)

    ! calculate cumulative sum
    do isp = 1, nsp
       cumsum(1,isp) = 0
       do i = 1, nproc
          cumsum(i+1,isp) = cumsum(i,isp) + gcount(isp,i)
       end do
    end do

  end subroutine get_global_cumsum

  !
  ! initial velocity profile
  !
  function vprofile(x) result(y)
    implicit none
    real(8), intent(in) :: x
    real(8) :: y
    real(8) :: x0, xs

    x0 = nxgs*delx + l_damp_ini
    xs = l_damp_ini * 0.25
    y  = 0.5d0 * v0 * (1 + tanh((x-x0)/xs))

  end function vprofile

  !
  ! initialize raondom seed
  !
  subroutine init_random_seed()
    implicit none
    integer :: n
    integer, allocatable :: seed(:)

    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)*(nrank+1)
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_random_seed

  !
  ! uniform random number
  !
  function uniform_rand() result(y)
    implicit none
    real(8) :: y

    call random_number(y)

  end function uniform_rand

  !
  ! normal random number via Box-Muller method
  !
  function normal_rand() result(y)
    implicit none
    real(8) :: y
    real(8) :: xx(2), rr
    real(8), save :: yy
    logical, save :: has_saved

    if ( has_saved ) then
       y = yy
       has_saved = .false.
    else
       call random_number(xx)
       rr = sqrt(-2*log(xx(1)))
       yy = rr * cos(2*pi*xx(2))
       y  = rr * sin(2*pi*xx(2))
       has_saved = .true.
    end if

  end function normal_rand

  !
  ! get elapsed time
  !
  function get_etime() result(y)
    implicit none
    real(8) :: y

    if ( nrank == nroot ) then
       y = MPI_Wtime()
    end if

    call MPI_Bcast(y, 1, MPI_REAL8, nroot, MPI_COMM_WORLD, mpierr)

  end function get_etime

end module app
