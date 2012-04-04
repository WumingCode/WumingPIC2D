module init

  use const
  use mpi_set

  implicit none

  private

  public :: init__set_param, init__inject, init__relocate

  integer, public, parameter   :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public              :: itmax, it0, intvl1, intvl2, intvl3, intvl4
  integer, public              :: nxs
  integer, public              :: nxe
  integer, public              :: nxs1
  integer, public              :: nxe1
  real(8), public              :: delx, delt, gfac
  real(8), public              :: c, q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=128), public   :: dir
  character(len=128), public   :: file12
  real(8), save                :: pi, n0, u0, v0, b0, vti, vte, gam0


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param

    integer              :: n
    integer, allocatable :: seed(:)
    real(8)              :: fgi, fpi, alpha, beta, va, fpe, fge, rgi, rge, ldb, rtemp
    character(len=128)   :: file9 
    character(len=128)   :: file11

!************** MPI settings  *******************!
    call mpi_set__init(nygs,nyge,nproc)
!*********** End of MPI settings  ***************!

!************* Physical region ******************!
    nxs  = nxgs
    nxs1 = nxs-1
    nxe  = nxge
    nxe1 = nxe+1
!!$    nxs  = nxgs
!!$    nxs1 = nxs-1
!!$    nxe  = nxs+nx*0.2-1
!!$    nxe1 = nxe+1
!****************   End of  * *******************!

!*********** Memory Allocations  ****************!
    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxgs-1:nxge+1,nys1:nye1))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
!***************** ENd of  **********************!

!!$!*********** Random seed *************!
    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)*(nrank+1)
    call random_seed(put=seed)
    deallocate(seed)
!!$!***********   End of    *************!

!*********************************************************************
!   time0   : start time (if time0 < 0, initial data from input.f)
!   itmax   : number of iteration
!   it0     : base count
!   intvl1  : storage interval for particles & fields
!   intvl2  : printing interval for energy variation
!   intvl3  : interval for injecting particles
!   intvl4  : interval for updating physical region in x
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
    itmax  = 720000
    intvl1 = 2000
    intvl2 = 90000
    intvl3 = 3
    intvl4 = 20
!!$    dir    = '../../dat/shock/test/'          !for pc
!!$    dir    = './pic/shock/test/'              !for hx600@nagoya, xt@nao
    dir    = '/large/m/m082/pic/shock/run5/'   !for fx1@jaxa
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.505
    it0    = 0

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
    delt = 0.25
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
    v0  = -20.0*va
    u0  = v0/dsqrt(1.-(v0/c)**2)
    gam0 = dsqrt(1.+u0**2/c**2)

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    !average number density at x=nxgs (magnetosheath)
    n0 = 40

    if(nrank == nroot)then
       if(n0*(nxge-nxgs) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    !number of particles in each cell in y
    np2(nys:nye,1:nsp) = n0*(nxe-nxs)*delx

    !charge
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(a,i3.3,a)')'9999999_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,nxs,nxe,c,q,r,delt,delx,it0,          &
                       np,nxgs,nxge,nygs,nyge,nys,nye,nsp,nproc,nrank, &
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
    real(8) :: sd, aa, bb, cc, gamp

    !*** setting of fields ***!
    !magnetic field
!$OMP PARALLEL

!$OMP DO PRIVATE(i,j)
    do j=nys-1,nye+1
    do i=nxgs,nxge-1
       uf(1,i,j) = 0.0D0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j)
    do j=nys-1,nye+1
    do i=nxgs,nxge
       uf(2,i,j) = 0.0D0
       uf(3,i,j) = b0
    enddo
    enddo
!$OMP END DO NOWAIT

    !electric field
!$OMP DO PRIVATE(i,j)
    do j=nys-1,nye+1
    do i=nxgs,nxge
       uf(4,i,j) = 0.0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j)
    do j=nys-1,nye+1
    do i=nxgs,nxge-1
       uf(5,i,j) = v0*b0/c
       uf(6,i,j) = 0.0
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
    !*** end of ***!

    !particle position
    isp = 1
!$OMP PARALLEL DO PRIVATE(ii,j,aa)
    do j=nys,nye
       do ii=1,np2(j,isp)
          up(1,ii,j,1) = nxs*delx+(nxe-nxs)*delx*ii/(np2(j,isp)+1)
          up(1,ii,j,2) = up(1,ii,j,1)

          aa = 0.0D0
          do while(aa==1.D0)
             call random_number(aa)
          enddo

          up(2,ii,j,1) = dble(j)*delx+delx*aa
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo
!$OMP END PARALLEL DO

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/sqrt(2.)
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc)
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

  end subroutine init__loading


  subroutine init__relocate

    use boundary, only : boundary__field
!$  use omp_lib

    integer :: dn, isp, j, ii, ii2 ,ii3
    real(8) :: aa, bb, cc, sd, gamp

!!$    if(nxs==nxgs) return
    if(nxe==nxge) return

!!$    nxs  = nxs-1
!!$    nxs1 = nxs-1
    nxe  = nxe+1
    nxe1 = nxe1+1

    dn = n0

    !particle position
!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

!!$          up(1,ii2,j,1) = nxs*delx+delx*ii/(dn+1)
          up(1,ii2,j,1) = (nxe-1)*delx+delx*ii/(dn+1)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          aa = 0.0D0
          do while(aa==1.D0)
             call random_number(aa)
          enddo

          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
!$OMP END PARALLEL DO

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/sqrt(2.)
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc)
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

!$OMP PARALLEL WORKSHARE
       np2(nys:nye,isp) = np2(nys:nye,isp)+dn
!$OMP END PARALLEL WORKSHARE

    enddo

!$OMP PARALLEL DO PRIVATE(j)
    do j=nys,nye
       uf(1,nxe-1,j) = 0.0D0
       uf(2,nxe-1,j) = 0.0D0
       uf(3,nxe-1,j) = b0
       uf(4,nxe-1,j) = 0.0D0
       uf(5,nxe-1,j) = v0*b0/c
       uf(6,nxe-1,j) = 0.0D0
       uf(2,nxe,j)   = 0.0D0
       uf(3,nxe,j)   = b0
       uf(4,nxe,j)   = 0.0D0
    enddo
!$OMP END PARALLEL DO

    call boundary__field(uf,                        &
                         nxgs,nxge,nxs,nxe,nys,nye, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

  end subroutine init__relocate


  subroutine init__inject

    use boundary, only : boundary__field

    integer :: isp, ii, ii2, ii3, j, dn
    real(8) :: sd, aa, bb, cc, dx, fac, gamp

    !Inject particles in x=nxe-v0*dt~nxe*dt

    dx  = v0*delt*intvl3/delx
    dn  = abs(n0*dx)
    fac = (dn)/(0.30+dn)

!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

!!$          up(1,ii2,j,1) = nxs*delx+dx*ii/(dn+1)
          up(1,ii2,j,1) = nxe*delx+dx*(dn-ii+1)/(dn+1)
          up(1,ii3,j,2) = up(1,ii2,j,1)

          aa = 0.0D0
          do while(aa==1.D0)
             call random_number(aa)
          enddo

          up(2,ii2,j,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo
!$OMP END PARALLEL DO

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/dsqrt(2.0D0)/fac
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)/fac
          gamp = dsqrt(1.D0+sd*sd/(c*c))
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,aa,bb)
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
!$OMP PARALLEL DO PRIVATE(j)
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+dn
       enddo
!$OMP END PARALLEL DO
    enddo

    !set Ex and Bz
!$OMP PARALLEL DO PRIVATE(j)
    do j=nys,nye
!!$       uf(3,nxs,j)   = b0
!!$       uf(5,nxs,j)   = v0*b0/c
!!$       uf(3,nxs+1,j) = b0
       uf(3,nxe-1,j) = b0
       uf(5,nxe-1,j) = v0*b0/c
       uf(3,nxe,j)   = b0
    enddo
!$OMP END PARALLEL DO

    call boundary__field(uf,                        &
                         nxgs,nxge,nxs,nxe,nys,nye, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

  end subroutine init__inject


end module init
