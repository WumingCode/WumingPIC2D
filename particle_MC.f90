module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,up,uf,                     &
                            np,nsp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye, &
                            c,q,r,delt,delx)

    integer, intent(in)  :: np, nxgs, nxge, nxs, nxe, nys, nye, nsp
    integer, intent(in)  :: cumcnt(nxgs:nxge,nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(in)  :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)

    integer :: i, j, ii, isp
    real(8) :: idelx, dh, sh(-1:1,2)
    real(8) :: fac1, fac2, txxx, fac1r, fac2r, gam, igam
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8) :: tmp(1:6,nxs-1:nxe+1,nys-1:nye+1)

    idelx = 1./delx

    !fields at (i+1/2, j+1/2)
!$OMP PARALLEL DO PRIVATE(i,j) 
    do j=nys-1,nye+1
    do i=nxs-1,nxe+1
       tmp(1,i,j) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(2,i,j) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(3,i,j) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(4,i,j) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(5,i,j) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(6,i,j) = uf(6,i,j)
    enddo
    enddo
!OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,i,j,isp,sh,dh,gam,igam,fac1,fac2,txxx,fac1r,fac2r, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,uvm1,uvm2,uvm3,uvm4,uvm5,uvm6) 
    do j=nys,nye
    do i=nxs,nxe-1

       isp=1

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

       do ii=cumcnt(i,j,isp)+1,cumcnt(i+1,j,isp)

          !second order shape function
          dh = up(1,ii,j,isp)*idelx-0.5-i
          sh(-1,1) = 0.5*(0.5-dh)*(0.5-dh)
          sh( 0,1) = 0.75-dh*dh
          sh(+1,1) = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(2,ii,j,isp)*idelx-0.5-j
          sh(-1,2) = 0.5*(0.5-dh)*(0.5-dh)
          sh( 0,2) = 0.75-dh*dh
          sh(+1,2) = 0.5*(0.5+dh)*(0.5+dh)

          bpx = +(+tmp(1,i-1,j-1)*sh(-1,1)+tmp(1,i,j-1)*sh(0,1)+tmp(1,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(1,i-1,j  )*sh(-1,1)+tmp(1,i,j  )*sh(0,1)+tmp(1,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(1,i-1,j+1)*sh(-1,1)+tmp(1,i,j+1)*sh(0,1)+tmp(1,i+1,j+1)*sh(+1,1))*sh(+1,2)

          bpy = +(+tmp(2,i-1,j-1)*sh(-1,1)+tmp(2,i,j-1)*sh(0,1)+tmp(2,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(2,i-1,j  )*sh(-1,1)+tmp(2,i,j  )*sh(0,1)+tmp(2,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(2,i-1,j+1)*sh(-1,1)+tmp(2,i,j+1)*sh(0,1)+tmp(2,i+1,j+1)*sh(+1,1))*sh(+1,2)

          bpz = +(+tmp(3,i-1,j-1)*sh(-1,1)+tmp(3,i,j-1)*sh(0,1)+tmp(3,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(3,i-1,j  )*sh(-1,1)+tmp(3,i,j  )*sh(0,1)+tmp(3,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(3,i-1,j+1)*sh(-1,1)+tmp(3,i,j+1)*sh(0,1)+tmp(3,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epx = +(+tmp(4,i-1,j-1)*sh(-1,1)+tmp(4,i,j-1)*sh(0,1)+tmp(4,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(4,i-1,j  )*sh(-1,1)+tmp(4,i,j  )*sh(0,1)+tmp(4,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(4,i-1,j+1)*sh(-1,1)+tmp(4,i,j+1)*sh(0,1)+tmp(4,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epy = +(+tmp(5,i-1,j-1)*sh(-1,1)+tmp(5,i,j-1)*sh(0,1)+tmp(5,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(5,i-1,j  )*sh(-1,1)+tmp(5,i,j  )*sh(0,1)+tmp(5,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(5,i-1,j+1)*sh(-1,1)+tmp(5,i,j+1)*sh(0,1)+tmp(5,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epz = +(+tmp(6,i-1,j-1)*sh(-1,1)+tmp(6,i,j-1)*sh(0,1)+tmp(6,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(6,i-1,j  )*sh(-1,1)+tmp(6,i,j  )*sh(0,1)+tmp(6,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(6,i-1,j+1)*sh(-1,1)+tmp(6,i,j+1)*sh(0,1)+tmp(6,i+1,j+1)*sh(+1,1))*sh(+1,2)

          !accel.
          uvm1 = up(3,ii,j,isp)+fac1*epx
          uvm2 = up(4,ii,j,isp)+fac1*epy
          uvm3 = up(5,ii,j,isp)+fac1*epz

          !rotate
          gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
          igam = 1./gam
          fac1r = fac1*igam
          fac2r = fac2/(gam+txxx*(bpx*bpx+bpy*bpy+bpz*bpz)*igam)

          uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
          uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
          uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

          uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
          uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
          uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

          !accel.
          gp(3,ii,j,isp) = uvm1+fac1*epx
          gp(4,ii,j,isp) = uvm2+fac1*epy
          gp(5,ii,j,isp) = uvm3+fac1*epz

          !move
          gam = 1./dsqrt(1.0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                               +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                               +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))

          gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
          gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam

       enddo

       isp=2

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

       do ii=cumcnt(i,j,isp)+1,cumcnt(i+1,j,isp)

          !second order shape function
          dh = up(1,ii,j,isp)*idelx-0.5-i
          sh(-1,1) = 0.5*(0.5-dh)*(0.5-dh)
          sh( 0,1) = 0.75-dh*dh
          sh(+1,1) = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(2,ii,j,isp)*idelx-0.5-j
          sh(-1,2) = 0.5*(0.5-dh)*(0.5-dh)
          sh( 0,2) = 0.75-dh*dh
          sh(+1,2) = 0.5*(0.5+dh)*(0.5+dh)

          bpx = +(+tmp(1,i-1,j-1)*sh(-1,1)+tmp(1,i,j-1)*sh(0,1)+tmp(1,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(1,i-1,j  )*sh(-1,1)+tmp(1,i,j  )*sh(0,1)+tmp(1,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(1,i-1,j+1)*sh(-1,1)+tmp(1,i,j+1)*sh(0,1)+tmp(1,i+1,j+1)*sh(+1,1))*sh(+1,2)

          bpy = +(+tmp(2,i-1,j-1)*sh(-1,1)+tmp(2,i,j-1)*sh(0,1)+tmp(2,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(2,i-1,j  )*sh(-1,1)+tmp(2,i,j  )*sh(0,1)+tmp(2,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(2,i-1,j+1)*sh(-1,1)+tmp(2,i,j+1)*sh(0,1)+tmp(2,i+1,j+1)*sh(+1,1))*sh(+1,2)

          bpz = +(+tmp(3,i-1,j-1)*sh(-1,1)+tmp(3,i,j-1)*sh(0,1)+tmp(3,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(3,i-1,j  )*sh(-1,1)+tmp(3,i,j  )*sh(0,1)+tmp(3,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(3,i-1,j+1)*sh(-1,1)+tmp(3,i,j+1)*sh(0,1)+tmp(3,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epx = +(+tmp(4,i-1,j-1)*sh(-1,1)+tmp(4,i,j-1)*sh(0,1)+tmp(4,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(4,i-1,j  )*sh(-1,1)+tmp(4,i,j  )*sh(0,1)+tmp(4,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(4,i-1,j+1)*sh(-1,1)+tmp(4,i,j+1)*sh(0,1)+tmp(4,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epy = +(+tmp(5,i-1,j-1)*sh(-1,1)+tmp(5,i,j-1)*sh(0,1)+tmp(5,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(5,i-1,j  )*sh(-1,1)+tmp(5,i,j  )*sh(0,1)+tmp(5,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(5,i-1,j+1)*sh(-1,1)+tmp(5,i,j+1)*sh(0,1)+tmp(5,i+1,j+1)*sh(+1,1))*sh(+1,2)

          epz = +(+tmp(6,i-1,j-1)*sh(-1,1)+tmp(6,i,j-1)*sh(0,1)+tmp(6,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                +(+tmp(6,i-1,j  )*sh(-1,1)+tmp(6,i,j  )*sh(0,1)+tmp(6,i+1,j  )*sh(+1,1))*sh( 0,2) &
                +(+tmp(6,i-1,j+1)*sh(-1,1)+tmp(6,i,j+1)*sh(0,1)+tmp(6,i+1,j+1)*sh(+1,1))*sh(+1,2)

          !accel.
          uvm1 = up(3,ii,j,isp)+fac1*epx
          uvm2 = up(4,ii,j,isp)+fac1*epy
          uvm3 = up(5,ii,j,isp)+fac1*epz

          !rotate
          gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
          igam = 1./gam
          fac1r = fac1*igam
          fac2r = fac2/(gam+txxx*(bpx*bpx+bpy*bpy+bpz*bpz)*igam)

          uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
          uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
          uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

          uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
          uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
          uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

          !accel.
          gp(3,ii,j,isp) = uvm1+fac1*epx
          gp(4,ii,j,isp) = uvm2+fac1*epy
          gp(5,ii,j,isp) = uvm3+fac1*epz

          !move
          gam = 1./dsqrt(1.0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                              +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                              +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))

          gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
          gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam

       enddo

    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine particle__solv


end module particle
