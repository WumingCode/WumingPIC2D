module init

  use const
  use mpi_set

  implicit none

  private

  public :: init__set_param, init__inject

  integer, public, parameter   :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public              :: itmax, it0, intvl1, intvl2
  real(8), public              :: delx, delt, gfac
  real(8), public              :: c, q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=128), public   :: dir
  character(len=128), public   :: file12
  real(8), save                :: pi, n0, u0, v0, b0, vti, vte


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param

    integer              :: n
    integer, allocatable :: seed(:)
    real(8)              :: fgi, fpi, alpha, beta, va, fpe, fge, rgi, rge, ldb, rtemp
    character(len=128)   :: file9 
    character(len=128)   :: file11

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)

    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs1:nxe1,nys1:nye1))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
!*********** End of MPI settings  ***************!

!*********** Random seed *************!
    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)+nrank
    call random_seed(put=seed)
    deallocate(seed)
!***********   End of    *************!

!*********************************************************************
!   time0   : start time (if time0 < 0, initial data from input.f)
!   itmax   : number of iteration
!   it0     : base count
!   intvl1  : storage interval for particles & fields
!   intvl2  : printing interval for energy variation
!   dir     : directory name for data output
!   file??  : output file name for unit number ??
!           :  9 - initial parameters
!           : 10 - for saving all data
!           : 11 - for starting from saved data
!           : 12 - for saving energy history
!   gfac    : implicit factor
!             gfac < 0.5 : unstable
!             gfac = 0.5 : no implicit
!             gfac = 1.0 : full implicit
!*********************************************************************
    pi     = 4.0*atan(1.0)
    itmax  = 100000
    intvl1 = 5000
    intvl2 = 5000
!!$    dir    = '../../dat/shock/test/'          !for pc
!!$    dir    = './pic/shock/test/'              !for hx600
    dir    = '/large/m/m082/pic/shock/run1/'   !for fx1@jaxa
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.505
    it0    = 1

!*********************************************************************
!   r(1)  : ion mass             r(2)  : electron mass
!   q(1)  : ion charge           q(2)  : electron charge
!   c     : speed of light       ldb   : debye length
!
!   rgi   : ion Larmor radius    rge   : electron Larmor radius
!   fgi   : ion gyro-frequency   fge   : electron gyro-frequency
!   vti   : ion thermal speed    vte   : electron thermal speed
!   b0    : magnetic field       
!  
!   alpha : wpe/wge
!   beta  : ion plasma beta
!   rtemp : Te/Ti
!*********************************************************************
    delx = 1.0
    c    = 1.0
    delt = 1.0
    ldb  = delx

    r(1) = 100.0
    r(2) = 1.0

    alpha = 10.0
    beta  = 0.5
    rtemp = 1.0

    fpe = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge = fpe/alpha

    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = alpha*ldb*dsqrt(2.D0)
    rgi = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)
    vte = rge*fge
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)
    v0  = 10.0*va
    u0  = v0/dsqrt(1.-(v0/c)**2)

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    !average number density at x=nxgs (magnetosheath)
    n0 = 40.

    if(nrank == nroot)then
       if(n0*(nxge+bc-nxgs+1) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    !number of particles in each cell in y
    np2(nys:nye,1:nsp) = n0*(nxge-nxgs)*delx

    !charge
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(a,i3.3,a)')'9999999_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,                             &
                       np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,bc,nproc,nrank, &
                       dir,file11)
       return
    endif

    call init__loading
    call fio__param(np,nsp,np2,                             &
                    nxgs,nxge,nygs,nyge,nys,nye,            &
                    c,q,r,n0,0.5*r(1)*vti**2,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,                &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__field

    integer :: i, j, ii, isp
    real(8) :: sd, aa, bb, cc

    !*** setting of fields ***!
    !magnetic field
!$OMP PARALLEL

!$OMP DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = 0.0D0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = 0.0D0
       uf(3,i,j) = b0
    enddo
    enddo
!$OMP END DO NOWAIT

    !electric field
!$OMP DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe
       uf(4,i,j) = 0.0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(5,i,j) = v0*b0/c
       uf(6,i,j) = 0.0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call boundary__field(uf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)
    !*** end of ***!

    !particle position
    isp = 1
!$OMP PARALLEL DO PRIVATE(ii,j,aa)
    do j=nys,nye
       do ii=1,np2(j,isp)
          call random_number(aa)
          up(1,ii,j,1) = nxs*delx+aa*delx*(nxe+bc-nxs+1.)
          up(1,ii,j,2) = up(1,ii,j,1)

          call random_number(aa)
          up(2,ii,j,1) = dble(j)*delx+delx*aa
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo
!$OMP END PARALLEL DO

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(ii,j,sd,aa,bb,cc)
       do j=nys,nye
          do ii=1,np2(j,isp)
             if(isp .eq. 1) then 
                sd = vti/sqrt(2.)
             endif
             if(isp .eq. 2) then
                sd = vte/sqrt(2.)
             endif

             aa = 0.0D0
             do while(aa == 0.0D0)
                call random_number(aa)
             enddo
             sd = sd*dsqrt(-2.*dlog(aa))
             call random_number(bb)
             call random_number(cc)

             up(3,ii,j,isp) = sd*(2.*bb-1)+u0
             up(4,ii,j,isp) = sd*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = sd*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
          enddo
       enddo
!$OMP END PARALLEL DO
    enddo

  end subroutine init__loading


  subroutine init__inject

    use boundary, only : boundary__field

    integer :: isp, ii, ii2, ii3, j, dn
    real(8) :: sd, aa, bb, cc, dx

    !Inject particles in x=nxs~nxs+v0*dt

    dx  = v0*delt/delx
    dn  = n0*dx

!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii
          call random_number(aa)
          up(1,ii2,j,1) = nxs*delx+aa*dx
          up(1,ii3,j,2) = up(1,ii2,j,1)

          call random_number(aa)
          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
!$OMP END PARALLEL DO

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc)
       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn
             aa = 0.0D0
             do while(aa == 0.0D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)+u0
             up(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
          enddo
       enddo
!$OMP END PARALLEL DO
    enddo

    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(j)
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+dn
       enddo
!$OMP END PARALLEL DO
    enddo

    !set Ex and Bz
!$OMP PARALLEL DO PRIVATE(j)
    do j=nys,nye
       uf(3,nxs,j)  = b0
       uf(5,nxs,j)  = v0*b0/c
    enddo
!$OMP END PARALLEL DO

    call boundary__field(uf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

  end subroutine init__inject


end module init
