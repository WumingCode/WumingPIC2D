module project
  use wuming2d
  implicit none
  private

  ! main simulation loop
  public:: project__main

  !************************ NUMERICAL CONSTANTS ***********************************!
  integer, parameter :: nx     = 1000      ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny     = 100       ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs   = 2         ! START POINT IN X
  integer, parameter :: nxge   = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs   = 2         ! START POINT IN Y
  integer, parameter :: nyge   = nygs+ny-1 ! END POINT
  integer, parameter :: ndim   = 6         ! DIMENSION OF PHASE SPACE--5:2D-3V, 6:2D-3V+ID FOR TRACKING
  integer, parameter :: np     = 6*nx      ! MAX. NUMBER OF PARTICLES IN CONUMN AT Y
  integer, parameter :: nsp    = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc  = 4         ! NUMBER OF PROCESSES
  integer, parameter :: nroot  = 0         ! ROOT PROC. NUMBER

  ! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs    = nxgs      ! START POINT IN X
  integer            :: nxe    = nxgs+500 ! INITIAL SIZE

  ! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 1000      !NUMBER OF ITERATION
  integer            :: it0    = 0         !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 10        !INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl2 = 1         !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3 = 1         !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  integer, parameter :: intvl4 = 10        !INTERVAL FOR RECORDING MOMENT DATA
  integer, parameter :: intvl5 = 5         !INTERVAL FOR RECORDING ORBIT DATA
  character(len=128) :: dir    = './'      !DIRECTORY FOR OUTPUT
  character(len=128) :: dir2   = './'      !DIRECTORY FOR OUTPUT
  character(len=128) :: file9  = 'init_param' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: etlim  = 23.5*60*60 !MAX. ELAPSE TIME IN SEC.

  ! OTHER CONSTANTS
  real(8)            :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl   = 10.0      !ELECTRON SKIN DEPTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*datan(1.0D0)

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
  integer, parameter :: n0     = 2
  real(8), parameter :: mr     = 1.0D0
  real(8), parameter :: gam0   = 40.0D0
  real(8), parameter :: sig0   = 1.D-1
  real(8), parameter :: rtemp  = 1.0D0
  real(8), parameter :: vte    = 0.0D0
  real(8), parameter :: vti    = sqrt(1.0D0/(rtemp*mr))*vte
  real(8), parameter :: theta  = 90.0D0/180.D0*pi
  real(8), parameter :: phi    = 90.0D0/180.D0*pi
  real(8)            :: ldmp   = 10.0D0

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
  real(8), save                :: u0, v0, b0, delt

contains
  !
  ! main simulation loop
  !
  subroutine project__main()
    implicit none
    integer :: it=0
    real(8) :: etime, etime0, omp_get_wtime

    etime0 = omp_get_wtime()
    call set_param()
    call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)

    do it=1,itmax-it0
       if(nrank == nroot) etime = omp_get_wtime()
       call MPI_BCAST(etime,1,mnpr,nroot,ncomw,nerr)
       if(etime-etime0 >= etlim) then
          call h5io__output(up,uf,np2,nxs,nxe,it-1+it0,.true.)
          if(nrank == nroot) write(*,*) '*** elapse time over ***',it-1+it0,etime-etime0
          exit
       endif

       call particle__solv(gp,up,uf,cumcnt,nxs,nxe)
       call boundary__particle_injection(gp,np2,nxs,nxe, &
            mod(it+it0-1,intvl2),mod(it+it0-1,intvl3))
       call field__fdtd_i(uf,up,gp,cumcnt,nxs,nxe, &
            & boundary__dfield, &
            & boundary__curre, &
            & boundary__phi)
       call boundary__particle_y(gp,np2)
       call sort__bucket(up,gp,cumcnt,np2,nxs,nxe)

       if(mod(it+it0,intvl2) == 0) call inject(it+it0)
       if(mod(it+it0,intvl3) == 0) call relocate(it+it0)

       if(mod(it+it0,intvl1) == 0) call h5io__output(up,uf,np2,nxs,nxe,it+it0,.false.)
       if(ndim == 6 .and. mod(it+it0,intvl5) == 0) call h5io__orb(up,uf,np2,it+it0)

       if(mod(it+it0,intvl4) == 0)then
          call mom_calc__accl(gp,up,uf,cumcnt,nxs,nxe)
          call mom_calc__nvt(den,vel,temp,gp,np2)
          call boundary__mom(den,vel,temp)
          call h5io__mom(den,vel,temp,uf,it+it0)
       endif

    enddo

    call h5util_finalize()
    call MPI_FINALIZE(nerr)

  end subroutine project__main

  subroutine set_param()
    implicit none
    integer              :: n, isp, i, j, ndim_in
    integer, allocatable :: seed(:)
    real(8)              :: fpe
    real(8)              :: r(nsp), q(nsp)
    character(len=128)   :: file11

    !************** MPI SETTINGS *******************!
    call mpi_set__init(nygs,nyge,nproc)
    !*********** End of MPI settings  ***************!

    !*********** MEMORY ALLOCATIONS ****************!
    allocate(np2(nys:nye,nsp))
    allocate(cumcnt(nxgs:nxge,nys:nye,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2))
    allocate(up(ndim,np,nys:nye,nsp))
    allocate(gp(ndim,np,nys:nye,nsp))
    allocate(den(nxgs-1:nxge+1,nys-1:nye+1,1:nsp))
    allocate(vel(nxgs-1:nxge+1,nys-1:nye+1,1:3,1:nsp))
    allocate(temp(nxgs-1:nxge+1,nys-1:nye+1,1:3,1:nsp))
    !***************** END OF  **********************!

    !*********** RANDOM SEED *************!
    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)*(nrank+1)
    call random_seed(put=seed)
    deallocate(seed)
    !***********   END OF    *************!

    !**** SETTING OTHER NUMERICAL & PHYSICAL CONSTANTS ****!
    r(2) = 1.0D0    ! ELECTRON MASS
    r(1) = r(2)*mr  ! ION MASS
    delt = cfl*delx/c   ! TIME STEP SIZE
    fpe  = c/delx/rdbl*sqrt(gam0)   ! ELE. PLASMA FREQ. IN SIM. FRAME
    ldmp = c/fpe*sqrt(gam0)*ldmp    ! INITIAL VELOCITY PROFILE NEAR X=0

    !CHARGE
    q(2) = -fpe*sqrt(r(2)/(4.0D0*pi*n0/delx**2))
    q(1) = -q(2)

    !MAGNETIC FIELD STRENGTH
    b0 = sqrt(sig0*4.0D0*pi*r(1)*n0*gam0*c**2)

    !INJECTOR IS ON THE RIGHT-HAND-SIDE; MINUS SIGN IS NECESSARY
    v0   = -c*sqrt(1.0D0-1.0D0/gam0**2)
    u0   = v0*gam0

    !NUMBER OF PARTICLES IN CELL COLUMN IN X AT Y
    np2(nys:nye,1:nsp) = n0*(nxe-nxs)
    if(nrank == nroot)then
       if(n0*(nxge-nxgs) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    !PREPAREATION FOR SORT
    do isp=1,nsp
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j=nys,nye
          cumcnt(nxs,j,isp) = 0
          do i=nxs+1,nxe
             cumcnt(i,j,isp) = cumcnt(i-1,j,isp)+n0
          enddo
          if(cumcnt(nxe,j,isp) /= np2(j,isp))then
             write(*,*)'error in cumcnt'
             stop
          endif
       enddo
       !$OMP END PARALLEL DO
    enddo

    !INITIALIZATION OF SUBROUTINES
    call boundary__init(ndim,np,nsp,                          &
         nxgs,nxge,nygs,nyge,nys,nye,          &
         nup,ndown,mnpi,mnpr,ncomw,nerr,nstat, &
         delx,delt,u0,c)
    call particle__init(ndim,np,nsp,                 &
         nxgs,nxge,nygs,nyge,nys,nye, &
         delx,delt,c,q,r)
    call field__init(ndim,np,nsp,                 &
         nxgs,nxge,nygs,nyge,nys,nye, &
         mnpr,ncomw,opsum,nerr,       &
         delx,delt,c,q,r,gfac)
    call sort__init(ndim,np,nsp,                &
         nxgs,nxge,nygs,nyge,nys,nye)
    call h5io__init(ndim,np,nsp,                 &
         nxgs,nxge,nygs,nyge,nys,nye, &
         nproc,nrank,                 &
         delx,delt,c,q,r,dir)
    call mom_calc__init(ndim,np,nsp,nxgs,nxge,nygs,nyge,nys,nye, &
         delx,delt,c,q,r)

    if(it0 /= 0)then
       !RESTART FROM THE PAST CALCULATION
       write(file11,'(i7.7,a)')it0,'.h5'
       call h5io__input(gp,uf,np2,ndim_in,nxs,nxe,it0,file11)
       call sort__bucket(up,gp,cumcnt,np2,nxs,nxe)
       if(ndim_in == 5 .and. ndim == 6) call indexpos()
       return
    endif

    call loading()
    if(ndim == 6) call indexpos()

    call h5io__param(n0,np2,                                     &
         0.5*r(1)*vti**2,rtemp,fpe,q(1)*b0/(r(2)*c), &
         rdbl*delx, file9,                           &
         nroot)

  end subroutine set_param


  subroutine loading()
    implicit none
    integer :: i, j, ii, isp
    real(8) :: sd, aa, bb, cc, gamp, v1, gam1
    real(8) :: vfunc,x0,eps=1d-40

    !INITIAL VELOCITY PROFILE
    vfunc(x0) = 0.5*v0*(1.d0+tanh( (x0-ldmp)/(0.1*ldmp) ))

    !--- SETTING OF INITIAL FIELDS ---!
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys-2,nye+2
       do i=nxgs-2,nxge+2
          uf(1,i,j) = b0*cos(theta)
          uf(2,i,j) = b0*sin(theta)*cos(phi)
          uf(3,i,j) = b0*sin(theta)*sin(phi)
          uf(4,i,j) = 0.0D0
          uf(5,i,j) = vfunc(i*delx)*uf(3,i,j)/c
          uf(6,i,j) = -vfunc(i*delx)*uf(2,i,j)/c
       enddo
    enddo
    !$OMP END PARALLEL DO

    !particle position
    isp = 1
    !$OMP PARALLEL DO PRIVATE(ii,j,aa)
    do j=nys,nye
       do ii=1,np2(j,isp)
          up(1,ii,j,1) = nxs*delx+(nxe-nxs)*delx*ii/(np2(j,isp)+1)
          up(1,ii,j,2) = up(1,ii,j,1)

          call random_number(aa)
          up(2,ii,j,1) = dble(j)*delx+delx*aa
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then
          sd = vti/sqrt(2.)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
       endif

       !$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc,gamp,v1,gam1)
       do j=nys,nye
          do ii=1,np2(j,isp)

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)
             up(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = dsqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/c**2)

             v1 = vfunc(up(1,ii,j,isp))
             gam1 = 1.0D0/dsqrt(1.0D0-v1**2/c**2)

             call random_number(cc)

             if(up(3,ii,j,isp)*v1 >= 0.)then
                up(3,ii,j,isp) = (+up(3,ii,j,isp)+v1*gamp)*gam1
             else
                if(cc < (-v1*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = (-up(3,ii,j,isp)+v1*gamp)*gam1
                else
                   up(3,ii,j,isp) = (+up(3,ii,j,isp)+v1*gamp)*gam1
                endif
             endif

          enddo
       enddo
       !$OMP END PARALLEL DO

    enddo

  end subroutine loading


  subroutine relocate(it)
    implicit none
    integer, intent(in) :: it
    integer :: isp, j, ii, ii2 ,ii3, dn
    real(8) :: aa, bb, cc, sd, gamp, dx

    if(nxe==nxge) return

    dx  = v0*delt*mod(it,intvl2)/delx
    dn  = int(n0*abs(dx)+0.5)
    nxe  = nxe+1

    !PARTICLE POSITION
    !$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          up(1,ii2,j,1) = (nxe-1.D0+dx)*delx-dx*delx*ii/(dn+1.D0)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          call random_number(aa)
          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
       do ii=1,n0
          ii2 = np2(j,1)+dn+ii
          ii3 = np2(j,2)+dn+ii

          up(1,ii2,j,1) = (nxe-1.)*delx+delx*ii/(n0+1.)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          call random_number(aa)
          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then
          sd = vti/sqrt(2.)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
       endif

       !$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc,gamp)
       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn+n0

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)
             up(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = dsqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/c**2)

             call random_number(cc)

             if(up(3,ii,j,isp)*v0 >= 0.)then
                up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = (-up(3,ii,j,isp)+v0*gamp)*gam0
                else
                   up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
                endif
             endif
          enddo
          np2(j,isp) = np2(j,isp)+dn+n0
          cumcnt(nxe-1,j,isp) = cumcnt(nxe-1,j,isp)+dn
          cumcnt(nxe,j,isp) = cumcnt(nxe-1,j,isp)+n0
       enddo
       !$OMP END PARALLEL DO

    enddo

    !$OMP PARALLEL DO PRIVATE(j)
    do j=nys-2,nye+2
       uf(2,nxe-1,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe-1,j) = b0*sin(theta)*sin(phi)
       uf(5,nxe-1,j) = v0*uf(3,nxe-1,j)/c
       uf(6,nxe-1,j) = -v0*uf(2,nxe-1,j)/c

       uf(2,nxe,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe,j) = b0*sin(theta)*sin(phi)
    enddo
    !$OMP END PARALLEL DO

  end subroutine relocate


  subroutine inject(it)
    implicit none
    integer :: it, isp, ii, ii2, ii3, j, dn
    real(8) :: sd, aa, bb, cc, gamp, dx

    !INJECT PARTICLES IN x=nxe-v0*dt~nxe*dt
    dx  = v0*delt*(it-max(int(it/intvl3)*intvl3,it-intvl2))/delx
    if(dx == 0.0D0) dx = v0*delt*min(intvl2,intvl3)/delx
    if(nxe == nxge) dx = v0*delt*intvl2/delx
    dn  = int(abs(n0*dx)+0.5)

    !$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          up(1,ii2,j,1) = nxe*delx+dx*(dn-ii+1)/(dn+1)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          call random_number(aa)
          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)
       endif

       !$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc,gamp)
       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)
             up(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = dsqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/c**2)

             call random_number(cc)

             if(up(3,ii,j,isp)*v0 >= 0.)then
                up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = (-up(3,ii,j,isp)+v0*gamp)*gam0
                else
                   up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
                endif
             endif
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
       uf(2,nxe-1,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe-1,j) = b0*sin(theta)*sin(phi)
       uf(5,nxe-1,j) = v0*uf(3,nxe-1,j)/c
       uf(6,nxe-1,j) = -v0*uf(2,nxe-1,j)/c

       uf(2,nxe,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe,j) = b0*sin(theta)*sin(phi)
    enddo
    !$OMP END PARALLEL DO

  end subroutine inject


  subroutine indexpos()

    integer :: isp, j, ii
    real(8) :: aa

    !index particle position
    do isp=1,nsp

       !$OMP PARALLEL DO PRIVATE(ii,j,aa)
       do j=nys,nye
          do ii=1,np2(j,isp)
             call random_number(aa)
             if(up(1,ii,j,isp) >= xrs .and. up(1,ii,j,isp) <= xre)then
                up(6,ii,j,isp) = aa
             else
                up(6,ii,j,isp) = aa-1.D0
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO

    enddo

  end subroutine indexpos

end module project
