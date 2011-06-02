module boundary

  implicit none

  private

  public  :: boundary__field
  public  :: boundary__particle
  public  :: boundary__curre


contains


  subroutine boundary__particle(up,                                        &
                                np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                                nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)
    integer, intent(in)        :: np, nsp, nxgs, nxge, nygs, nyge, nys, nye, bc
    integer, intent(in)        :: nup, ndown, mnpi, mnpr, ncomw
    integer, intent(inout)     :: nerr, nstat(:)
    integer, intent(inout)     :: np2(nys:nye,nsp)
    real(8), intent(inout)     :: up(5,np,nys:nye,nsp)
    logical, save              :: lflag=.true.
    integer                    :: j, ii, iii, isp, ipos, jpos
    integer                    :: cnt(nys-1:nye+1), cnt2(nys:nye), cnt_tmp
    integer, save, allocatable :: flag(:,:)
    real(8), save, allocatable :: bff_ptcl(:,:)

    if(lflag)then
       allocate(flag(np,nys:nye))
       allocate(bff_ptcl(np,nys-1:nye+1))
       lflag=.false.
    endif

    do isp=1,nsp

       cnt(nys-1:nye+1) = 0
       cnt2(nys:nye) = 0

       do j=nys,nye
          do ii=1,np2(j,isp)

             ipos = floor(up(1,ii,j,isp))
             jpos = floor(up(2,ii,j,isp))

             if(bc==0)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)+(nxge-nxgs+1)
                endif
                if(ipos >= nxge+1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)-(nxge-nxgs+1)
                endif
             else if(bc==-1)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = 2.0*nxgs-up(1,ii,j,isp)
                   up(3,ii,j,isp) = -up(3,ii,j,isp)
                endif
                if(ipos >= nxge)then
                   up(1,ii,j,isp) = 2.0*nxge-up(1,ii,j,isp)
                   up(3,ii,j,isp) = -up(3,ii,j,isp)
                endif
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif

             if(jpos /= j)then
                if(jpos <= nygs-1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)+(nyge-nygs+1)
                endif
                if(jpos >= nyge+1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)-(nyge-nygs+1)
                endif
                bff_ptcl(1+5*cnt(jpos),jpos) = up(1,ii,j,isp)
                bff_ptcl(2+5*cnt(jpos),jpos) = up(2,ii,j,isp)
                bff_ptcl(3+5*cnt(jpos),jpos) = up(3,ii,j,isp)
                bff_ptcl(4+5*cnt(jpos),jpos) = up(4,ii,j,isp)
                bff_ptcl(5+5*cnt(jpos),jpos) = up(5,ii,j,isp)
                cnt(jpos) = cnt(jpos)+1
                cnt2(j) = cnt2(j)+1
                flag(cnt2(j),j) = ii
             endif
          enddo
       enddo

       !transfer to rank-1
       call MPI_SENDRECV(cnt(nys-1),1,mnpi,ndown,100, &
                         cnt_tmp   ,1,mnpi,nup  ,100, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1           ,nys-1),5*cnt(nys-1),mnpr,ndown,101, &
                         bff_ptcl(5*cnt(nye)+1,nye  ),5*cnt_tmp   ,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)
       cnt(nye) = cnt(nye)+cnt_tmp

       !transfer to rank+1
       call MPI_SENDRECV(cnt(nye+1),1,mnpi,nup  ,200, &
                         cnt_tmp   ,1,mnpi,ndown,200, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1           ,nye+1),5*cnt(nye+1),mnpr,nup  ,201, &
                         bff_ptcl(5*cnt(nys)+1,nys  ),5*cnt_tmp   ,mnpr,ndown,201, &
                         ncomw,nstat,nerr)
       cnt(nys) = cnt(nys)+cnt_tmp

       do j=nys,nye
          iii=0
          cnt_tmp = cnt2(j)
          loop1 :do ii=1,cnt2(j)
             if(cnt(j) == 0)then
                if(np2(j,isp) < flag(ii,j)) exit loop1
                do while(np2(j,isp) == flag(cnt_tmp,j))
                   np2(j,isp) = np2(j,isp)-1
                   if(np2(j,isp) < flag(ii,j)) exit loop1
                   cnt_tmp = cnt_tmp-1
                enddo
                up(1:6,flag(ii,j),j,isp) = up(1:6,np2(j,isp),j,isp)
                np2(j,isp) = np2(j,isp)-1
             else
                up(1,flag(ii,j),j,isp) = bff_ptcl(1+5*iii,j)
                up(2,flag(ii,j),j,isp) = bff_ptcl(2+5*iii,j)
                up(3,flag(ii,j),j,isp) = bff_ptcl(3+5*iii,j)
                up(4,flag(ii,j),j,isp) = bff_ptcl(4+5*iii,j)
                up(5,flag(ii,j),j,isp) = bff_ptcl(5+5*iii,j)
                iii = iii+1
                cnt(j) = cnt(j)-1
             endif
          enddo loop1
          
          if(cnt(j) > 0)then
             do ii=1,cnt(j)
                up(1,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+1+5*(ii-1),j)
                up(2,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+2+5*(ii-1),j)
                up(3,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+3+5*(ii-1),j)
                up(4,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+4+5*(ii-1),j)
                up(5,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+5+5*(ii-1),j)
             enddo
          endif
       enddo

       do j=nys,nye
          np2(j,isp) = np2(j,isp)+cnt(j)
          if(np2(j,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(j,isp),j,isp
             stop
          endif
       enddo

    enddo

  end subroutine boundary__particle


  subroutine boundary__field(uf,                 &
                             nxs,nxe,nys,nye,bc, &
                             nup,ndown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nys, nye, bc
    integer, intent(in)    :: nup, ndown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uf(6,nxs-1:nxe+1,nys-1:nye+1)
    integer                :: i, ii
    real(8)                :: bff_snd(6*(nxe-nxs+1)), bff_rcv(6*(nxe-nxs+1))

    do i=nxs,nxe
       ii = 6*(i-nxs)
       bff_snd(ii+1) = uf(1,i,nys)
       bff_snd(ii+2) = uf(2,i,nys)
       bff_snd(ii+3) = uf(3,i,nys)
       bff_snd(ii+4) = uf(4,i,nys)
       bff_snd(ii+5) = uf(5,i,nys)
       bff_snd(ii+6) = uf(6,i,nys)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+1),mnpr,ndown,110, &
                      bff_rcv(1),6*(nxe-nxs+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = 6*(i-nxs)
       uf(1,i,nye+1) = bff_rcv(ii+1)   
       uf(2,i,nye+1) = bff_rcv(ii+2)
       uf(3,i,nye+1) = bff_rcv(ii+3)
       uf(4,i,nye+1) = bff_rcv(ii+4)   
       uf(5,i,nye+1) = bff_rcv(ii+5)
       uf(6,i,nye+1) = bff_rcv(ii+6)
    enddo

    do i=nxs,nxe
       ii = 6*(i-nxs)
       bff_snd(ii+1) = uf(1,i,nye)
       bff_snd(ii+2) = uf(2,i,nye)
       bff_snd(ii+3) = uf(3,i,nye)
       bff_snd(ii+4) = uf(4,i,nye)
       bff_snd(ii+5) = uf(5,i,nye)
       bff_snd(ii+6) = uf(6,i,nye)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+1),mnpr,nup  ,100, &
                      bff_rcv(1),6*(nxe-nxs+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = 6*(i-nxs)
       uf(1,i,nys-1) = bff_rcv(ii+1)   
       uf(2,i,nys-1) = bff_rcv(ii+2)
       uf(3,i,nys-1) = bff_rcv(ii+3)
       uf(4,i,nys-1) = bff_rcv(ii+4)   
       uf(5,i,nys-1) = bff_rcv(ii+5)
       uf(6,i,nys-1) = bff_rcv(ii+6)
    enddo

    if(bc == 0)then
       uf(1:6,nxs-1,nys-1:nye+1) = uf(1:6,nxe,nys-1:nye+1)
       uf(1:6,nxe+1,nys-1:nye+1) = uf(1:6,nxs,nys-1:nye+1)
    else if(bc == -1)then
       uf(1  ,nxs-1,nys-1:nye+1) = -uf(1  ,nxs  ,nys-1:nye+1)
       uf(2:4,nxs-1,nys-1:nye+1) = +uf(2:4,nxs+1,nys-1:nye+1)
       uf(5:6,nxs-1,nys-1:nye+1) = -uf(5:6,nxs  ,nys-1:nye+1)

       uf(1  ,nxe  ,nys-1:nye+1) = -uf(1  ,nxe-1,nys-1:nye+1)
       uf(2:4,nxe+1,nys-1:nye+1) = +uf(2:4,nxe-1,nys-1:nye+1)
       uf(5:6,nxe  ,nys-1:nye+1) = -uf(5:6,nxe-1,nys-1:nye+1)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__field


  subroutine boundary__curre(uj,nxs,nxe,nys,nye,bc, &
                             nup,ndown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nys, nye, bc
    integer, intent(in)    :: nup, ndown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    integer                :: i, ii
    real(8)                :: bff_rcv(6*(nxe-nxs+4+1)), bff_snd(6*(nxe-nxs+4+1))

    !send to rank-1
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys-2)
       bff_snd(ii+2) = uj(2,i,nys-2)
       bff_snd(ii+3) = uj(3,i,nys-2)
       bff_snd(ii+4) = uj(1,i,nys-1)
       bff_snd(ii+5) = uj(2,i,nys-1)
       bff_snd(ii+6) = uj(3,i,nys-1)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,ndown,110, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nye-1) = uj(1,i,nye-1)+bff_rcv(ii+1)
       uj(2,i,nye-1) = uj(2,i,nye-1)+bff_rcv(ii+2)
       uj(3,i,nye-1) = uj(3,i,nye-1)+bff_rcv(ii+3)
       uj(1,i,nye  ) = uj(1,i,nye  )+bff_rcv(ii+4)
       uj(2,i,nye  ) = uj(2,i,nye  )+bff_rcv(ii+5)
       uj(3,i,nye  ) = uj(3,i,nye  )+bff_rcv(ii+6)
    enddo

    !send to rank+1
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye+1)
       bff_snd(ii+2) = uj(2,i,nye+1)
       bff_snd(ii+3) = uj(3,i,nye+1)
       bff_snd(ii+4) = uj(1,i,nye+2)
       bff_snd(ii+5) = uj(2,i,nye+2)
       bff_snd(ii+6) = uj(3,i,nye+2)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,nup  ,100, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nys  ) = uj(1,i,nys  )+bff_rcv(ii+1)
       uj(2,i,nys  ) = uj(2,i,nys  )+bff_rcv(ii+2)
       uj(3,i,nys  ) = uj(3,i,nys  )+bff_rcv(ii+3)
       uj(1,i,nys+1) = uj(1,i,nys+1)+bff_rcv(ii+4)
       uj(2,i,nys+1) = uj(2,i,nys+1)+bff_rcv(ii+5)
       uj(3,i,nys+1) = uj(3,i,nys+1)+bff_rcv(ii+6)
    enddo

    !boundary condition in x
    if(bc == 0)then
       uj(1:3,nxs  ,nys:nye) = uj(1:3,nxs  ,nys:nye)+uj(1:3,nxe+1,nys:nye)
       uj(1:3,nxs+1,nys:nye) = uj(1:3,nxs+1,nys:nye)+uj(1:3,nxe+2,nys:nye)

       uj(1:3,nxe-1,nys:nye) = uj(1:3,nxe-1,nys:nye)+uj(1:3,nxs-2,nys:nye)
       uj(1:3,nxe  ,nys:nye) = uj(1:3,nxe  ,nys:nye)+uj(1:3,nxs-1,nys:nye)
    else if(bc == -1)then
       uj(2:3,nxs  ,nys:nye) = uj(2:3,nxs  ,nys:nye)-uj(2:3,nxs-1,nys:nye)
       uj(2:3,nxs+1,nys:nye) = uj(2:3,nxs+1,nys:nye)-uj(2:3,nxs-2,nys:nye)

       uj(2:3,nxe-2,nys:nye) = uj(2:3,nxe-2,nys:nye)-uj(2:3,nxe+1,nys:nye)
       uj(2:3,nxe-1,nys:nye) = uj(2:3,nxe-1,nys:nye)-uj(2:3,nxe  ,nys:nye)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__curre


end module boundary
