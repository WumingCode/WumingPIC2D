module sort

  implicit none

  private

!  public :: sort__insert, sort__insert2, sort__bucket
  public :: sort__bucket

contains


!  subroutine sort__insert(up,np,nsp,np2,nys,nye)

!    integer, intent(in)    :: np, nsp, nys, nye
!    integer, intent(in)    :: np2(nys:nye,nsp)
!    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
!    integer                :: j, ii, iii, isp
!    real(8)                :: insrt(5)

!    do isp=1,nsp

!      !SIMPLE INSERT SORT ALGORITHM FOR ORDERD PARTICLES IN X
!!$OMP PARALLEL DO PRIVATE(iii,ii,j,insrt)
!       do j=nys,nye
!          do ii=2,np2(j,isp)
!             insrt(1:5) = up(1:5,ii,j,isp)
!             if(int(up(1,ii-1,j,isp)) > int(insrt(1)))then
!                iii = ii
!                do while(int(up(1,iii-1,j,isp)) > int(insrt(1)) .and. iii > 1)
!                   up(1:5,iii,j,isp) = up(1:5,iii-1,j,isp)
!                   iii = iii-1
!                enddo
!                up(1:5,iii,j,isp) = insrt(1:5)
!             endif
!          enddo
!       enddo
!!$OMP END PARALLEL DO

!!       if(isp==1)then
!!          do ii=1,np2(4,isp)
!!            write(*,*)ii,int(up(1,ii,4,isp)),up(1,ii,4,isp),np2(4,isp)
!!          enddo
!!       endif

!    enddo
!  end subroutine sort__insert


!  subroutine sort__insert2(up,np,nsp,np2,nys,nye)

!    integer, intent(in)    :: np, nsp, nys, nye
!    integer, intent(in)    :: np2(nys:nye,nsp)
!    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
!    integer                :: j, ii, iii, isp, left, right, mid
!    real(8)                :: insrt(5)

!    do isp=1,nsp

!      !SIMPLE INSERT SORT ALGORITHM FOR ORDERD PARTICLES IN X
!!$OMP PARALLEL DO PRIVATE(iii,ii,j,left,right,mid,insrt)
!       do j=nys,nye
!          do ii=2,np2(j,isp)

!             insrt(1:5) = up(1:5,ii,j,isp)
!             left = 1
!             right = ii

!             do while(left < right)

!                mid = (left+right)/2

!                if(int(up(1,mid,j,isp)) < int(insrt(1)))then 
!                   left = mid+1
!                else
!                   right = mid
!                endif

!             enddo
!  
!             iii = ii
!             do while(iii > left)
!                up(1:5,iii,j,isp) = up(1:5,iii-1,j,isp)
!                iii = iii-1
!             enddo
!             up(1:5,left,j,isp) = insrt(1:5)

!          enddo
!       enddo
!!$OMP END PARALLEL DO
!       
!!       if(isp==1)then
!!          do ii=1,np2(4,isp)
!!            write(*,*)ii,int(up(1,ii,4,isp)),up(1,ii,4,isp),np2(4,isp)
!!          enddo
!!       endif

!    enddo

!  end subroutine sort__insert2


  subroutine sort__bucket(up,np,nsp,np2,nxs,nxe,nys,nye)

    integer, intent(in)    :: np, nsp, nxs, nxe, nys, nye
    integer, intent(in)    :: np2(nys:nye,nsp)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
    logical, save              :: lflag=.true.
    integer                    :: i, j, ii, isp, np2max
    integer                    :: cnt(nxs:nxe-1), sum_cnt(nxs:nxe-1) 
    real(8), save, allocatable :: tmp(:,:,:)

    if(lflag)then
       allocate(tmp(1:5,np,nys:nye))
       lflag=.false.
    endif

    do isp=1,nsp

!      !BUCKET SORT FOR PARTICLES IN X
!$OMP PARALLEL DO PRIVATE(ii,i,j,cnt,sum_cnt)
       do j=nys,nye
          
          cnt(nxs:nxe-1) = 0
          tmp(1:5,1:np2(j,isp),j) = up(1:5,1:np2(j,isp),j,isp)

          do ii=1,np2(j,isp)
             i = int(tmp(1,ii,j))
             cnt(i) = cnt(i)+1
          enddo

          sum_cnt(nxs) = 0
          do i=nxs+1,nxe-1
             sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
          enddo

          do ii=1,np2(j,isp)
             i = int(tmp(1,ii,j))
             up(1:5,sum_cnt(i)+1,j,isp) = tmp(1:5,ii,j)
             sum_cnt(i) = sum_cnt(i)+1
          enddo

       enddo
!$OMP END PARALLEL DO
       
!       if(isp==1)then
!          do ii=1,np2(4,isp)
!            write(*,*)ii,int(up(1,ii,4,isp)),up(1,ii,4,isp),np2(4,isp)
!          enddo
!       endif

    enddo

  end subroutine sort__bucket


end module sort
