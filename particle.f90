module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,up,uf,                     &
                            np,nsp,np2,nxgs,nxge,nys,nye, &
                            delt,c,q,r)

    integer, intent(in)  :: np, nxgs, nxge, nys, nye, nsp
    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxgs-1:nxge+1,nys-1:nye+1)
    real(8), intent(in)  :: delt, c, q(nsp), r(nsp)
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)
    integer :: i, j, ii, isp, ih, jh
    real(8) :: dx, dxm, dy, dym
    real(8) :: fac1, fac1r, fac2, fac2r, gam, igam, txxx, bt2
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

!$OMP PARALLEL DO PRIVATE(ii,i,j,ih,jh,dx,dxm,dy,dym,bt2,gam,igam,fac1r,fac2r, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
       do j=nys,nye
          do ii=1,np2(j,isp)
             !interpolate fields to particles

             !Bx at (i+1/2, j)
             i  = int(up(1,ii,j,isp))
             ih = int(up(1,ii,j,isp)-0.5)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             bpx = +(dxm*uf(1,ih,j  )+dx*uf(1,ih+1,j  ))*dym &
                   +(dxm*uf(1,ih,j+1)+dx*uf(1,ih+1,j+1))*dy

             !By at (i, j+1/2)
             jh = int(up(2,ii,j,isp)-0.5)
             dx = up(1,ii,j,isp)-i
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             bpy = +(dxm*uf(2,i,jh  )+dx*uf(2,i+1,jh  ))*dym &
                   +(dxm*uf(2,i,jh+1)+dx*uf(2,i+1,jh+1))*dy

             !Bz at (i, j)
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             bpz = +(dxm*uf(3,i,j  )+dx*uf(3,i+1,j  ))*dym &
                   +(dxm*uf(3,i,j+1)+dx*uf(3,i+1,j+1))*dy


             !Ex at (i, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             epx = +(dxm*uf(4,i,jh  )+dx*uf(4,i+1,jh  ))*dym &
                   +(dxm*uf(4,i,jh+1)+dx*uf(4,i+1,jh+1))*dy

             !Ey at (i+1/2, j)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             epy = +(dxm*uf(5,ih,j  )+dx*uf(5,ih+1,j  ))*dym &
                   +(dxm*uf(5,ih,j+1)+dx*uf(5,ih+1,j+1))*dy

             !Ez at (i+1/2, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             epz = +(dxm*uf(6,ih,jh  )+dx*uf(6,ih+1,jh  ))*dym &
                   +(dxm*uf(6,ih,jh+1)+dx*uf(6,ih+1,jh+1))*dy

             !accel.
             uvm1 = up(3,ii,j,isp)+fac1*epx
             uvm2 = up(4,ii,j,isp)+fac1*epy
             uvm3 = up(5,ii,j,isp)+fac1*epz

             !rotate
             gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
             igam = 1./gam
             bt2 = bpx*bpx+bpy*bpy+bpz*bpz
             fac1r = fac1*igam
             fac2r = fac2/(gam+txxx*bt2*igam)

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
!$OMP END PARALLEL DO

    enddo

  end subroutine particle__solv


end module particle
