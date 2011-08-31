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
    real(8) :: fac1, fac2, fac2r, fac3, fac3r, gam, igam, txxx, bt2
    real(8) :: pf(6)
    real(8) :: uvm(6)

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       fac2 = fac1/c
       txxx = fac1*fac1
       fac3 = q(isp)*delt/r(isp)

!$OMP PARALLEL DO PRIVATE(ii,i,j,ih,jh,dx,dxm,dy,dym,bt2,gam,igam,fac2r,fac3r,pf,uvm)
       do j=nys,nye
          do ii=1,np2(j,isp)
             !interpolate fields to particles

             ih = floor(up(1,ii,j,isp)-0.5)

             !Bx at (i+1/2, j)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(1) = +(dxm*uf(1,ih,j  )+dx*uf(1,ih+1,j  ))*dym &
                     +(dxm*uf(1,ih,j+1)+dx*uf(1,ih+1,j+1))*dy

             i  = floor(up(1,ii,j,isp))
             jh = floor(up(2,ii,j,isp)-0.5)
             !By at (i, j+1/2)
             dx = up(1,ii,j,isp)-i
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(2) = +(dxm*uf(2,i,jh  )+dx*uf(2,i+1,jh  ))*dym &
                     +(dxm*uf(2,i,jh+1)+dx*uf(2,i+1,jh+1))*dy

             !Bz at (i, j)
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(3) = +(dxm*uf(3,i,j  )+dx*uf(3,i+1,j  ))*dym &
                     +(dxm*uf(3,i,j+1)+dx*uf(3,i+1,j+1))*dy

             !Ex at (i, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(4) = +(dxm*uf(4,i,jh  )+dx*uf(4,i+1,jh  ))*dym &
                     +(dxm*uf(4,i,jh+1)+dx*uf(4,i+1,jh+1))*dy

             !Ey at (i+1/2, j)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(5) = +(dxm*uf(5,ih,j  )+dx*uf(5,ih+1,j  ))*dym &
                     +(dxm*uf(5,ih,j+1)+dx*uf(5,ih+1,j+1))*dy

             !Ez at (i+1/2, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(6) = +(dxm*uf(6,ih,jh  )+dx*uf(6,ih+1,jh  ))*dym &
                     +(dxm*uf(6,ih,jh+1)+dx*uf(6,ih+1,jh+1))*dy

             bt2 = pf(1)*pf(1)+pf(2)*pf(2)+pf(3)*pf(3)

             !accel.
             uvm(1) = up(3,ii,j,isp)+fac1*pf(4)
             uvm(2) = up(4,ii,j,isp)+fac1*pf(5)
             uvm(3) = up(5,ii,j,isp)+fac1*pf(6)

             !rotate
             gam = dsqrt(1.0+(+uvm(1)*uvm(1) &
                              +uvm(2)*uvm(2) &
                              +uvm(3)*uvm(3))/(c*c))
             igam = 1./gam
             fac2r = fac2*igam
             fac3r = fac3/(gam+txxx*bt2*igam)

             uvm(4) = uvm(1)+fac2r*(+uvm(2)*pf(3)-uvm(3)*pf(2))
             uvm(5) = uvm(2)+fac2r*(+uvm(3)*pf(1)-uvm(1)*pf(3))
             uvm(6) = uvm(3)+fac2r*(+uvm(1)*pf(2)-uvm(2)*pf(1))

             uvm(1) = uvm(1)+fac3r*(+uvm(5)*pf(3)-uvm(6)*pf(2))
             uvm(2) = uvm(2)+fac3r*(+uvm(6)*pf(1)-uvm(4)*pf(3))
             uvm(3) = uvm(3)+fac3r*(+uvm(4)*pf(2)-uvm(5)*pf(1))

             !accel.
             gp(3,ii,j,isp) = uvm(1)+fac1*pf(4)
             gp(4,ii,j,isp) = uvm(2)+fac1*pf(5)
             gp(5,ii,j,isp) = uvm(3)+fac1*pf(6)

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
