module init

  use const
  use mpi_set
  use sort, only : sort__bucket

  implicit none

  private

  public :: init__set_param, init__inject, init__relocate

  integer, allocatable, public :: np2(:,:), cumcnt(:,:,:)
  real(8), public              :: delt
  real(8), public              :: r(nsp), q(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  real(8), save                :: u0, v0, b0, vti, vte


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param

    integer              :: n, isp, i, j
    integer, allocatable :: seed(:)
    real(8)              :: ldb, fgi, fpi, va, fpe, fge, rgi, rge
    character(len=128)   :: file11

!************** MPI SETTINGS *******************!
    call mpi_set__init(nygs,nyge,nproc)
!*********** End of MPI settings  ***************!

!*********** MEMORY ALLOCATIONS ****************!
    allocate(np2(nys:nye,nsp))
    allocate(cumcnt(nxgs:nxge,nys:nye,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
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
    r(2) = 1.0D0      ! ELECTRON MASS
    r(1) = r(2)*mr    ! ION MASS
    delt = cfl*delx/c ! TIME STEP SIZE
    ldb  = delx*rdbl
    fpe  = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge  = fpe/alpha
    fgi  = fge*r(2)/r(1)
    fpi  = fpe*dsqrt(r(2)/r(1))
    va   = fge/fpe*c*dsqrt(r(2)/r(1))
    rge  = alpha*ldb*dsqrt(2.D0)
    rgi  = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)
    vte  = rge*fge
    vti  = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)

    !CHARGE
    q(1) = fpi*dsqrt(r(1)/(4.0D0*pi*n0))
    q(2) = -q(1)

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)

    !INJECTOR IS ON THE RIGHT-HAND-SIDE; MINUS SIGN IS NECESSARY
    v0   = -ma*va
    u0   = v0/dsqrt(1.-(v0/c)**2)

    !NUMBER OF PARTICLES IN CELL COLUMN IN X AT Y
    np2(nys:nye,1:nsp) = n0*(nxe-nxs)*delx
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

    if(it0 /= 0)then
       !RESTART FROM THE PAST CALCULATION
       write(file11,'(i7.7,a,i3.3,a)')it0,'_rank=',nrank,'.dat'
       call fio__input(gp,uf,np2,nxs,nxe,c,q,r,delt,delx,it0,          &
                       np,nxgs,nxge,nygs,nyge,nys,nye,nsp,nproc,nrank, &
                       dir,file11)
       call sort__bucket(up,gp,cumcnt,np,nsp,np2,nxgs,nxge,nxs,nxe,nys,nye)
       return
    endif

    call init__loading
    call fio__param(np,n0,nsp,np2,                       &
                    nxgs,nxge,nygs,nyge,nys,nye,         &
                    c,q,r,0.5*r(1)*vti**2,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,             &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading
!$  use omp_lib

    integer :: i, j, ii, isp
    real(8) :: sd, aa, bb, cc, gamp

!--- SETTING OF INITIAL FIELDS ---!
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys-2,nye+2
    do i=nxgs-2,nxge+2
       uf(1,i,j) = b0*cos(theta)
       uf(2,i,j) = b0*sin(theta)*cos(phi)
       uf(3,i,j) = b0*sin(theta)*sin(phi)
       uf(4,i,j) = 0.0D0
       uf(5,i,j) = v0*uf(3,i,j)/c
       uf(6,i,j) = -v0*uf(2,i,j)/c
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

!$OMP PARALLEL DO PRIVATE(ii,j,aa,bb,cc,gamp)
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

!$  use omp_lib

    integer :: dn, isp, j, ii, ii2 ,ii3
    real(8) :: aa, bb, cc, sd, gamp

    if(nxe==nxge) return
    nxe  = nxe+1

    dn = n0

    !PARTICLE POSITION
!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          gp(1,ii2,j,1) = (nxe-1)*delx+delx*ii/(dn+1)
          gp(1,ii3,j,2) = gp(1,ii2,j,1)

          call random_number(aa)
          gp(2,ii2,j,1) = dble(j)*delx+delx*aa
          gp(2,ii3,j,2) = gp(2,ii2,j,1)
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
          do ii=np2(j,isp)+1,np2(j,isp)+dn

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)
             
             gp(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)
             gp(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             gp(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = dsqrt(1.D0+(gp(3,ii,j,isp)**2+gp(4,ii,j,isp)**2+gp(5,ii,j,isp)**2)/c**2)

             call random_number(cc)

             if(gp(3,ii,j,isp)*v0 >= 0.)then
                gp(3,ii,j,isp) = (+gp(3,ii,j,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*gp(3,ii,j,isp)/gamp))then
                   gp(3,ii,j,isp) = (-gp(3,ii,j,isp)+v0*gamp)*gam0
                else
                   gp(3,ii,j,isp) = (+gp(3,ii,j,isp)+v0*gamp)*gam0
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
    do j=nys-2,nye+2
       uf(2,nxe-1,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe-1,j) = b0*sin(theta)*sin(phi)
       uf(5,nxe-1,j) = v0*uf(3,nxe-1,j)/c
       uf(6,nxe-1,j) = -v0*uf(2,nxe-1,j)/c

       uf(2,nxe,j) = b0*sin(theta)*cos(phi)
       uf(3,nxe,j) = b0*sin(theta)*sin(phi)
    enddo
!$OMP END PARALLEL DO

  end subroutine init__relocate


  subroutine init__inject

    integer :: isp, ii, ii2, ii3, j, dn
    real(8) :: sd, aa, bb, cc, dx, gamp

    !INJECT PARTICLES IN x=nxe-v0*dt~nxe*dt
    dx  = v0*delt*intvl2/delx
    dn  = abs(n0*dx)+0.5

!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,aa)
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii

          gp(1,ii2,j,1) = nxe*delx+dx*(dn-ii+1)/(dn+1)
          gp(1,ii3,j,2) = gp(1,ii2,j,1)

          call random_number(aa)
          gp(2,ii2,j,1) = dble(j)*delx+delx*aa
          gp(2,ii3,j,2) = gp(2,ii2,j,1)
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

             gp(3,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*(2.*bb-1)
             gp(4,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             gp(5,ii,j,isp) = sd*dsqrt(-2.*dlog(aa))*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = dsqrt(1.D0+(gp(3,ii,j,isp)**2+gp(4,ii,j,isp)**2+gp(5,ii,j,isp)**2)/c**2)

             call random_number(cc)

             if(gp(3,ii,j,isp)*v0 >= 0.)then
                gp(3,ii,j,isp) = (+gp(3,ii,j,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*gp(3,ii,j,isp)/gamp))then
                   gp(3,ii,j,isp) = (-gp(3,ii,j,isp)+v0*gamp)*gam0
                else
                   gp(3,ii,j,isp) = (+gp(3,ii,j,isp)+v0*gamp)*gam0
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

  end subroutine init__inject


end module init
