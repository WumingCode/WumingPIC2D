program main

  use const
  use mpi_set
  use init
  use boundary
  use fio
  use particle
  use field
  use sort, only : sort__bucket
  use mom_calc
  use h5io
  use h5util

  implicit none

  integer :: it=0
  real(8) :: etime, etime0, omp_get_wtime

!**********************************************************************c
!
!    two-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS    (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90   (by Y. Matsumoto, STEL)  2008/10/21
!    MPI parallelization (by Y. Matsumoto, STEL)  2009/4/1
!    2-D code            (by Y. Matsumoto, STEL)  2009/6/5
!    2-D code w SIMD     (by Y. Matsumoto, Chiba-U)  2013/4/1
!
!**********************************************************************c

  etime0 = omp_get_wtime()
  call init__set_param
  call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)

  loop: do it=1,itmax-it0

     if(nrank == nroot) etime = omp_get_wtime()
     call MPI_BCAST(etime,1,mnpr,nroot,ncomw,nerr)
     if(etime-etime0 >= etlim) then
        call h5io__output(up,uf,np2,nxs,nxe,it-1+it0,.true.)
        if(nrank == nroot) write(*,*) '*** elapse time over ***',it-1+it0,etime-etime0
        exit loop
     endif

     call particle__solv(gp,up,uf,cumcnt,nxs,nxe)
     call boundary__particle_injection(gp,np2,nxs,nxe, &
                                       mod(it+it0-1,intvl2),mod(it+it0-1,intvl3))
     call field__fdtd_i(uf,up,gp,cumcnt,nxs,nxe)
     call boundary__particle_y(gp,np2)
     call sort__bucket(up,gp,cumcnt,np2,nxs,nxe)

     if(mod(it+it0,intvl2) == 0) call init__inject(it+it0)
     if(mod(it+it0,intvl3) == 0) call init__relocate(it+it0)

     if(mod(it+it0,intvl1) == 0) call h5io__output(up,uf,np2,nxs,nxe,it+it0,.false.)
     if(ndim == 6 .and. mod(it+it0,intvl5) == 0) call h5io__orb(up,uf,np2,it+it0)

     if(mod(it+it0,intvl4) == 0)then
        call mom_calc__accl(gp,up,uf,cumcnt,nxs,nxe)
        call mom_calc__nvt(den,vel,temp,gp,np2)
        call boundary__mom(den,vel,temp)
        call h5io__mom(den,vel,temp,uf,it+it0)
     endif

  enddo loop

  call h5util_finalize()
  call MPI_FINALIZE(nerr)

end program main
