module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,                                &
                           np,nsp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye, &
                           q,c,delx,delt,gfac,                      &
                           nup,ndown,mnpr,opsum,nstat,ncomw,nerr)

    use boundary, only : boundary__dfield, boundary__curre
 
    integer, intent(in)    :: np, nsp, nxgs, nxge, nxs, nxe, nys, nye
    integer, intent(in)    :: nup, ndown, opsum, mnpr, ncomw
    integer, intent(in)    :: cumcnt(nxgs:nxge,nys:nye,nsp)
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(5,np,nys:nye,nsp)
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    logical, save              :: lflag=.true.
    integer                    :: i, j, ieq
    real(8)                    :: pi, f1, f2, f3
    real(8)                    :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    real(8), save, allocatable :: df(:,:,:), gkl(:,:,:)

    pi = 4.0*atan(1.0)

    if(lflag)then
       allocate(df(6,nxgs-2:nxge+2,nys-2:nye+2))
       allocate(gkl(3,nxgs-2:nxge+2,nys-2:nye+2))
!$OMP PARALLEL WORKSHARE
       df(1:6,nxgs-2:nxge+2,nys-2:nye+2) = 0.0D0
!$OMP END PARALLEL WORKSHARE
!$OMP PARALLEL WORKSHARE
       gkl(1:3,nxgs-2:nxge+2,nys-2:nye+2) = 0.0D0
!$OMP END PARALLEL WORKSHARE
       lflag=.false.
    endif

    call ele_cur(uj,up,gp, &
                 np,nsp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,q,c,delx,delt)

    call boundary__curre(uj,nxs,nxe,nys,nye, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    !calculation
    f1 = c*delt/delx
    f2 = gfac*f1*f1
    f3 = 4.0*pi*delx/c
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe
       gkl(1,i,j) = +f2*(+uf(1,i,j-1)                          &
                         +uf(1,i-1,j)-4.*uf(1,i,j)+uf(1,i+1,j) &
                         +uf(1,i,j+1)                          &
                         +f3*(-uj(3,i,j-1)+uj(3,i,j)) )        &
                    -f1*(-uf(6,i,j-1)+uf(6,i,j))
       gkl(2,i,j) = +f2*(+uf(2,i,j-1)                          &
                         +uf(2,i-1,j)-4.*uf(2,i,j)+uf(2,i+1,j) &
                         +uf(2,i,j+1)                          &
                         -f3*(-uj(3,i-1,j)+uj(3,i,j)) )        &
                    +f1*(-uf(6,i-1,j)+uf(6,i,j))
       gkl(3,i,j) = +f2*(+uf(3,i,j-1)                          &
                         +uf(3,i-1,j)-4.*uf(3,i,j)+uf(3,i+1,j) &
                         +uf(3,i,j+1)                          &
                         +f3*(-uj(2,i-1,j)+uj(2,i,j)           &
                              +uj(1,i,j-1)-uj(1,i,j)) )        &
                    -f1*(-uf(5,i-1,j)+uf(5,i,j)+uf(4,i,j-1)-uf(4,i,j))
    enddo
    enddo
!$OMP END PARALLEL DO

    !solve  < bx, by & bz >
    call cgm(df,gkl,                    &
             nxgs,nxge,nxs,nxe,nys,nye, &
             c,delx,delt,gfac,          &
             nup,ndown,mnpr,opsum,nstat,ncomw,nerr)

    call boundary__dfield(df,                        &
                          nxgs,nxge,nxs,nxe,nys,nye, &
                          nup,ndown,mnpr,nstat,ncomw,nerr)


    !solve  < ex, ey & ez >
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys,nye
    do i=nxs,nxe
       df(4,i,j) = +f1*(+gfac*(-df(3,i,j)+df(3,i,j+1))   &
                        +     (-uf(3,i,j)+uf(3,i,j+1)) ) &
                   -4.*pi*delt*uj(1,i,j)
       df(5,i,j) = -f1*(+gfac*(-df(3,i,j)+df(3,i+1,j))   &
                        +     (-uf(3,i,j)+uf(3,i+1,j)) ) &
                   -4.*pi*delt*uj(2,i,j)

       df(6,i,j) = +f1*(+gfac*(-df(2,i,j)+df(2,i+1,j)    &
                               +df(1,i,j)-df(1,i,j+1))   &
                        +     (-uf(2,i,j)+uf(2,i+1,j)    &
                               +uf(1,i,j)-uf(1,i,j+1)) ) &
                   -4.*pi*delt*uj(3,i,j)
    enddo
    enddo
!$OMP END PARALLEL DO

    call boundary__dfield(df,                        &
                          nxgs,nxge,nxs,nxe,nys,nye, &
                          nup,ndown,mnpr,nstat,ncomw,nerr)

    !===== Update fields and particles ======
    
!$OMP PARALLEL DO PRIVATE(i,j,ieq)
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       do ieq=1,6
          uf(ieq,i,j) = uf(ieq,i,j)+df(ieq,i,j)
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp, &
                     np,nsp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,q,c,delx,delt)

    integer, intent(in)  :: np, nsp, nxgs,nxge, nxs, nxe, nys, nye
    integer, intent(in)  :: cumcnt(nxgs:nxge,nys:nye,nsp)
    real(8), intent(in)  :: q(nsp), c, delx, delt
    real(8), intent(in)  :: gp(5,np,nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(out) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)

    integer            :: ii, i, j, isp, i2, inc, ip, jp
    real(8), parameter :: fac = 1.D0/3.D0
    real(8)            :: idelx, idelt, dh, gvz, s1_1, s1_2, s1_3, smo_1, smo_2, smo_3
    real(8)            :: s0(-2:2,2), ds(-2:2,2)
    real(8)            :: pjx(-2:2,-2:2), pjy(-2:2,-2:2), pjz(-2:2,-2:2), pjtmp(-2:2,-2:2)

!$OMP PARALLEL WORKSHARE
    uj(1:3,nxs-2:nxe+2,nys-2:nye+2) = 0.D0
!$OMP END PARALLEL WORKSHARE    

    idelt = 1.D0/delt
    idelx = 1.D0/delx

    !--------------Charge Conservation Method -------------!
    !---- Density Decomposition (Esirkepov, CPC, 2001) ----!
!$OMP PARALLEL DO PRIVATE(ii,i,j,i2,isp,inc,ip,jp,dh,s0,ds,pjx,pjy,pjz,pjtmp, &
!$OMP                     gvz,s1_1,s1_2,s1_3,smo_1,smo_2,smo_3) & 
!$OMP REDUCTION(+:uj) 
    do j=nys,nye
    do i=nxs,nxe-1

       pjx(-2:2,-2:2) = 0.D0
       pjy(-2:2,-2:2) = 0.D0
       pjz(-2:2,-2:2) = 0.D0

       isp=1
     
       do ii=cumcnt(i,j,isp)+1,cumcnt(i+1,j,isp)

          !second order shape function
          dh = up(1,ii,j,isp)*idelx-0.5-i
          s0(-2,1) = 0.D0
          s0(-1,1) = 0.5*(0.5-dh)*(0.5-dh)
          s0( 0,1) = 0.75-dh*dh
          s0(+1,1) = 0.5*(0.5+dh)*(0.5+dh)
          s0(+2,1) = 0.D0

          dh = up(2,ii,j,isp)*idelx-0.5-j
          s0(-2,2) = 0.D0
          s0(-1,2) = 0.5*(0.5-dh)*(0.5-dh)
          s0( 0,2) = 0.75-dh*dh
          s0(+1,2) = 0.5*(0.5+dh)*(0.5+dh)
          s0(+2,2) = 0.D0

          i2 = int(gp(1,ii,j,isp)*idelx)
          dh = gp(1,ii,j,isp)*idelx-0.5-i2
          inc = i2-i
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          ds(-2,1) = s1_1*smo_1
          ds(-1,1) = s1_1*smo_2+s1_2*smo_1
          ds( 0,1) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          ds(+1,1) = s1_3*smo_2+s1_2*smo_3
          ds(+2,1) = s1_3*smo_3

          i2 = int(gp(2,ii,j,isp)*idelx)
          dh = gp(2,ii,j,isp)*idelx-0.5-i2
          inc = i2-j
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          ds(-2,2) = s1_1*smo_1
          ds(-1,2) = s1_1*smo_2+s1_2*smo_1
          ds( 0,2) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          ds(+1,2) = s1_3*smo_2+s1_2*smo_3
          ds(+2,2) = s1_3*smo_3

          ds(-2:2,1:2) = ds(-2:2,1:2)-s0(-2:2,1:2)

          gvz = gp(5,ii,j,isp)/dsqrt(1.+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                         +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                         +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c) )

          pjtmp(-2:2,-2:2) = 0.D0
          do jp=-2,2
             do ip=-2,1
                pjtmp(ip+1,jp) = pjtmp(ip,jp) &
                                -q(isp)*delx*idelt*ds(ip,1)*(s0(jp,2)+0.5*ds(jp,2))
             enddo
          enddo
          pjx = pjx+pjtmp

          pjtmp(-2:2,-2:2) = 0.D0
          do jp=-2,1
             do ip=-2,2
                pjtmp(ip,jp+1) = pjtmp(ip,jp) &
                                -q(isp)*delx*idelt*ds(jp,2)*(s0(ip,1)+0.5*ds(ip,1))
             enddo
          enddo
          pjy = pjy+pjtmp

          do jp=-2,2
             do ip=-2,2
                pjz(ip,jp) = pjz(ip,jp)                                           &
                            +q(isp)*gvz*(+s0(ip,1)*s0(jp,2)+0.5*ds(ip,1)*s0(jp,2) &
                                         +0.5*s0(ip,1)*ds(jp,2)+fac*ds(ip,1)*ds(jp,2))
             enddo
          enddo

       enddo

       isp = 2

       do ii=cumcnt(i,j,isp)+1,cumcnt(i+1,j,isp)

          !second order shape function
          dh = up(1,ii,j,isp)*idelx-0.5-i
          s0(-2,1) = 0.D0
          s0(-1,1) = 0.5*(0.5-dh)*(0.5-dh)
          s0( 0,1) = 0.75-dh*dh
          s0(+1,1) = 0.5*(0.5+dh)*(0.5+dh)
          s0(+2,1) = 0.D0

          dh = up(2,ii,j,isp)*idelx-0.5-j
          s0(-2,2) = 0.D0
          s0(-1,2) = 0.5*(0.5-dh)*(0.5-dh)
          s0( 0,2) = 0.75-dh*dh
          s0(+1,2) = 0.5*(0.5+dh)*(0.5+dh)
          s0(+2,2) = 0.D0

          i2 = int(gp(1,ii,j,isp)*idelx)
          dh = gp(1,ii,j,isp)*idelx-0.5-i2
          inc = i2-i
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          ds(-2,1) = s1_1*smo_1
          ds(-1,1) = s1_1*smo_2+s1_2*smo_1
          ds( 0,1) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          ds(+1,1) = s1_3*smo_2+s1_2*smo_3
          ds(+2,1) = s1_3*smo_3

          i2 = int(gp(2,ii,j,isp)*idelx)
          dh = gp(2,ii,j,isp)*idelx-0.5-i2
          inc = i2-j
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          ds(-2,2) = s1_1*smo_1
          ds(-1,2) = s1_1*smo_2+s1_2*smo_1
          ds( 0,2) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          ds(+1,2) = s1_3*smo_2+s1_2*smo_3
          ds(+2,2) = s1_3*smo_3

          ds(-2:2,1:2) = ds(-2:2,1:2)-s0(-2:2,1:2)

          gvz = gp(5,ii,j,isp)/dsqrt(1.+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                          +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                          +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c) )

          pjtmp(-2:2,-2:2) = 0.D0
          do jp=-2,2
             do ip=-2,1
                pjtmp(ip+1,jp) = pjtmp(ip,jp) &
                                -q(isp)*delx*idelt*ds(ip,1)*(s0(jp,2)+0.5*ds(jp,2))
             enddo
          enddo
          pjx = pjx+pjtmp

          pjtmp(-2:2,-2:2) = 0.D0
          do jp=-2,1
             do ip=-2,2
                pjtmp(ip,jp+1) = pjtmp(ip,jp) &
                                -q(isp)*delx*idelt*ds(jp,2)*(s0(ip,1)+0.5*ds(ip,1))
             enddo
          enddo
          pjy = pjy+pjtmp

          do jp=-2,2
             do ip=-2,2
                pjz(ip,jp) = pjz(ip,jp)                                           &
                            +q(isp)*gvz*(+s0(ip,1)*s0(jp,2)+0.5*ds(ip,1)*s0(jp,2) &
                                         +0.5*s0(ip,1)*ds(jp,2)+fac*ds(ip,1)*ds(jp,2))
             enddo
          enddo

       enddo

       do jp=-2,2
          do ip=-2,2
             uj(1,i+ip,j+jp) = uj(1,i+ip,j+jp)+pjx(ip,jp)
             uj(2,i+ip,j+jp) = uj(2,i+ip,j+jp)+pjy(ip,jp)
             uj(3,i+ip,j+jp) = uj(3,i+ip,j+jp)+pjz(ip,jp)
          enddo
       enddo

    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine ele_cur


  subroutine cgm(db,gkl,                    &
                 nxgs,nxge,nxs,nxe,nys,nye, &
                 c,delx,delt,gfac,          &
                 nup,ndown,mnpr,opsum,nstat,ncomw,nerr)

    use boundary, only : boundary__phi

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !  #  this routine will be stoped after itaration number = ite_max
    !-----------------------------------------------------------------------

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye
    integer, intent(in)    :: nup, ndown, mnpr, opsum, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: c, delx, delt, gfac
    real(8), intent(in)    :: gkl(3,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: db(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer, parameter :: ite_max = 100 ! maximum number of interation
    integer            :: i, ii, j, l, ite
    real(8), parameter :: err = 1d-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: phi(nxs-1:nxe+1,nys-1:nye+1), b(nxs:nxe,nys:nye)
    real(8)            :: r(nxs:nxe,nys:nye), p(nxs-1:nxe+1,nys-1:nye+1)
    real(8)            :: ap(nxs:nxe,nys:nye)
    real(8)            :: bff_snd(nxe-nxs+1), bff_rcv(nxe-nxs+1)

    do l=1,1

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sum)
       do j=nys,nye
       do i=nxs,nxe-1
          phi(i,j) = db(l,i,j)
          b(i,j) = f2*gkl(l,i,j)
          sum = sum+b(i,j)*b(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !------ boundary condition of phi ------
       call boundary__phi(phi,               &
                          nxs,nxe,nys,nye,l, &
                          nup,ndown,mnpr,nstat,ncomw,nerr)
       !------ end of  ------

       f1 = 4.0+(delx/(c*delt*gfac))**2
       sumr = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sumr)
       do j=nys,nye
       do i=nxs,nxe-1
          r(i,j) = b(i,j)+phi(i,j-1)                    &
                         +phi(i-1,j)-f1*phi(i,j)+phi(i+1,j) &
                         +phi(i,j+1)
          p(i,j) = r(i,j)
          sumr = sumr+r(i,j)*r(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !------boundary condition of p------
             call boundary__phi(p,                 &
                                nxs,nxe,nys,nye,l, &
                                nup,ndown,mnpr,nstat,ncomw,nerr)
             !------ end of --------       
             sumr = 0.0
             sum2 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sumr,sum2)
             do j=nys,nye
             do i=nxs,nxe-1
                ap(i,j) = -p(i,j-1)                    &
                          -p(i-1,j)+f1*p(i,j)-p(i+1,j) &
                          -p(i,j+1)
                sumr = sumr+r(i,j)*r(i,j)
                sum2 = sum2+p(i,j)*ap(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g

!$OMP PARALLEL DO PRIVATE(i,j)
             do j=nys,nye
             do i=nxs,nxe-1
                phi(i,j) = phi(i,j)+av* p(i,j)
                r(i,j) = r(i,j)-av*ap(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO
             
             sum_g = dsqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif

             sum1 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sum1)
             do j=nys,nye
             do i=nxs,nxe-1
                sum1 = sum1+r(i,j)*r(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO

             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
!$OMP PARALLEL DO PRIVATE(i,j)
             do j=nys,nye
             do i=nxs,nxe-1
                p(i,j) = r(i,j)+bv*p(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO
             
          enddo
       endif

!$OMP PARALLEL WORKSHARE
       db(l,nxs:nxe-1,nys:nye) = phi(nxs:nxe-1,nys:nye)
!$OMP END PARALLEL WORKSHARE

    end do

    do l=2,3

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sum)
       do j=nys,nye
       do i=nxs,nxe
          phi(i,j) = db(l,i,j)
          b(i,j) = f2*gkl(l,i,j)
          sum = sum+b(i,j)*b(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !------ boundary condition of phi ------
       call boundary__phi(phi,               &
                          nxs,nxe,nys,nye,l, &
                          nup,ndown,mnpr,nstat,ncomw,nerr)
       !------ end of

       f1 = 4.0+(delx/(c*delt*gfac))**2
       sumr = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sumr)
       do j=nys,nye
       do i=nxs,nxe
          r(i,j) = b(i,j)+phi(i,j-1)                    &
                         +phi(i-1,j)-f1*phi(i,j)+phi(i+1,j) &
                         +phi(i,j+1)
          p(i,j) = r(i,j)
          sumr = sumr+r(i,j)*r(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !------boundary condition of p------
             call boundary__phi(p,                 &
                                nxs,nxe,nys,nye,l, &
                                nup,ndown,mnpr,nstat,ncomw,nerr)
             !------ end of

             sumr = 0.0
             sum2 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sumr,sum2)      
             do j=nys,nye
             do i=nxs,nxe
                ap(i,j) = -p(i,j-1)                    &
                          -p(i-1,j)+f1*p(i,j)-p(i+1,j) &
                          -p(i,j+1)
                sumr = sumr+r(i,j)*r(i,j)
                sum2 = sum2+p(i,j)*ap(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g

!$OMP PARALLEL DO PRIVATE(i,j)
             do j=nys,nye
             do i=nxs,nxe
                phi(i,j) = phi(i,j)+av* p(i,j)
                r(i,j) = r(i,j)-av*ap(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO
             
             sum_g = dsqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:sum1)
             do j=nys,nye
             do i=nxs,nxe
                sum1 = sum1+r(i,j)*r(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO
             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g

!$OMP PARALLEL DO PRIVATE(i,j)
             do j=nys,nye
             do i=nxs,nxe
                p(i,j) = r(i,j)+bv*p(i,j)
             enddo
             enddo
!$OMP END PARALLEL DO
             
          enddo
       endif

!$OMP PARALLEL WORKSHARE
       db(l,nxs:nxe,nys:nye) = phi(nxs:nxe,nys:nye)
!$OMP END PARALLEL WORKSHARE

    end do
    
  end subroutine cgm


end module field
