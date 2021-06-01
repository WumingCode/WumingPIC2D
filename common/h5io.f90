module h5io
  use hdf5
  use mpi
  implicit none

  private

  public :: h5io__init
  public :: h5io__output
  public :: h5io__input
  public :: h5io__param
  public :: h5io__mom
  public :: h5io__orb
  public :: h5util_finalize

  integer, parameter :: MAXDIM = 32
  integer :: hdferr

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)
  character(len=256), save :: dir

  ! initialization of hdf5 library
  interface h5util_init
     module procedure init
  end interface h5util_init

  ! finalization of hdf5 library
  interface h5util_finalize
     module procedure finalize
  end interface h5util_finalize

  ! create hdf5 file
  interface h5util_create_file
     module procedure create_file
  end interface h5util_create_file

  ! open existing file
  interface h5util_open_file
     module procedure open_file
  end interface h5util_open_file

  ! close file
  interface h5util_close_file
     module procedure close_file
  end interface h5util_close_file

  ! put attribute
  interface h5util_put_attribute
     module procedure &
          & put_attribute_char,   &
          & put_attribute_int,    &
          & put_attribute_r4,     &
          & put_attribute_r8,     &
          & put_attribute_arr_int,&
          & put_attribute_arr_r4, &
          & put_attribute_arr_r8
  end interface h5util_put_attribute

  ! create global dataset
  interface h5util_create_dataset
     module procedure create_dataset
  end interface h5util_create_dataset

  ! create global dataset with explicitly specified chunk
  interface h5util_create_dataset_chunk
     module procedure create_dataset_chunk
  end interface h5util_create_dataset_chunk

  ! write global dataset
  interface h5util_write_dataset
     module procedure &
          & write_dataset_int,&
          & write_dataset_r4, &
          & write_dataset_r8
  end interface h5util_write_dataset

  ! extend last (i.e., time) dimension
  interface h5util_extend_dimension
     module procedure extend_dimension
  end interface h5util_extend_dimension

  ! get array dimension
  interface h5util_get_dimension
     module procedure get_dimension
  end interface h5util_get_dimension

  ! read attribute
  interface h5util_read_attribute
     module procedure &
          & read_attribute_char,   &
          & read_attribute_int,    &
          & read_attribute_r4,     &
          & read_attribute_r8,     &
          & read_attribute_arr_int,&
          & read_attribute_arr_r4, &
          & read_attribute_arr_r8
  end interface h5util_read_attribute

  ! read dataset
  interface h5util_read_dataset
     module procedure &
          & read_dataset_int,&
          & read_dataset_r4, &
          & read_dataset_r8
  end interface h5util_read_dataset

contains

  subroutine h5io__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
                       nproc_in,nrank_in,                                                   &
                       delx_in,delt_in,c_in,q_in,r_in,dir_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    call h5util_init()

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in
    nproc = nproc_in
    nrank = nrank_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    dir   = dir_in

    is_init = .true.

  end subroutine h5io__init


  subroutine h5io__output(up,uf,np2,nxs,nxe,it0,lflag)

    logical, intent(in) :: lflag
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it0
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    character(len=256) :: filename

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(4), gdim(4), cdim(4), loff(4), goff(4)
    integer :: nd, nxg, nyg, nyl

    if(.not.is_init)then
       write(6,*)'Initialize first by calling h5io__init()'
       stop
    endif

    !filename
    if(lflag)then
       write(filename,'(a,i7.7,a)')trim(dir),9999999,'.h5'
    else
       write(filename,'(a,i7.7,a)')trim(dir),it0,'.h5'
    endif
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it0)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxs", nxs)
    call h5util_put_attribute(file_id, "nxe", nxe)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !particle
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 4
    ldim = (/ndim, np, nyl, nsp/)
    gdim = (/ndim, np, nyg, nsp/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)

    call h5util_create_dataset_chunk(file_id, "up", H5T_NATIVE_DOUBLE, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "up", nd, ldim, ldim, &
         & loff, goff, reshape(up, (/size(up)/)))

    nd   = 2
    ldim = (/nyl, nsp, 0, 0/)
    gdim = (/nyg, nsp, 0, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/nyl*nrank, 0, 0, 0/)

    call h5util_create_dataset_chunk(file_id, "np2", H5T_NATIVE_INTEGER, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "np2", nd, ldim, ldim, &
         & loff, goff, reshape(np2, (/size(np2)/)))

    !field
    nxg  = nxge-nxgs+5
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/6, nxg, nyl, 0/)
    gdim = (/6, nxg, nyg, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)

    call h5util_create_dataset_chunk(file_id, "uf", H5T_NATIVE_DOUBLE, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, &
         & reshape(uf(1:6,nxgs-2:nxge+2,nys:nye), (/6*nxg*nyl/)))


    call h5util_close_file(file_id)

  end subroutine h5io__output

  subroutine h5io__input(up,uf,np2,indim,nxs,nxe,it0,file)
    use mpi_set

    character(len=*), intent(in) :: file
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it0, indim
    real(8), intent(out) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, insp, inproc

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(4), loff(4), goff(4)
    integer :: i, ii, nd, nxg, nyl
    real(8) :: bff_snd(12*(nxge-nxgs+5)),bff_rcv(12*(nxge-nxgs+5))
    real(8),allocatable :: tmpp(:), tmpf(:)
    integer,allocatable :: tmpn(:)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling h5io__init()'
       stop
    endif

    !open file
    call h5util_open_file(trim(dir)//trim(file), file_id)

    !time & parameters
    call h5util_read_attribute(file_id, "it", it0)
    call h5util_read_attribute(file_id, "ndim", indim)
    call h5util_read_attribute(file_id, "np", inp)
    call h5util_read_attribute(file_id, "nxs", nxs)
    call h5util_read_attribute(file_id, "nxe", nxe)
    call h5util_read_attribute(file_id, "nxgs", inxgs)
    call h5util_read_attribute(file_id, "nxge", inxge)
    call h5util_read_attribute(file_id, "nygs", inygs)
    call h5util_read_attribute(file_id, "nyge", inyge)
    call h5util_read_attribute(file_id, "nsp", insp)
    call h5util_read_attribute(file_id, "nproc", inproc)
    call h5util_read_attribute(file_id, "delx", delx)
    call h5util_read_attribute(file_id, "delt", delt)
    call h5util_read_attribute(file_id, "c", c)
    call h5util_read_attribute(file_id, "m", r)
    call h5util_read_attribute(file_id, "q", q)

    if((inxgs /= nxgs) .or. (inxge /= nxge) .or.(inygs /= nygs) .or. (inyge /= nyge) &
        .or. (inp /= np) .or. (insp /= nsp) .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    !particle data
    nyl  = nye-nys+1

    nd   = 4
    ldim = (/indim, np, nyl, nsp/)
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)

    allocate(tmpp(indim*np*nyl*nsp))
    call h5util_read_dataset(file_id, "up", nd, ldim, ldim, loff, goff, tmpp)
    up(1:indim,1:np,nys:nye,1:nsp) = reshape(tmpp,ldim)

    nd   = 2
    ldim = (/nyl, nsp, 0, 0/)
    loff = (/0, 0, 0, 0/)
    goff = (/nyl*nrank, 0, 0, 0/)

    allocate(tmpn(nyl*nsp))
    call h5util_read_dataset(file_id, "np2", nd, ldim, ldim, loff, goff, tmpn)
    np2 = reshape(tmpn,ldim(1:2))

    !field data
    nxg  = nxge-nxgs+5
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/6, nxg, nyl, 0/)
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)

    allocate(tmpf(6*nxg*nyl))
    call h5util_read_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, tmpf)
    uf(1:6,nxgs-2:nxge+2,nys:nye) = reshape(tmpf,ldim(1:3))

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxgs-2,nxge+2
       ii = 12*(i-nxgs+2)
       bff_snd(ii+1)  = uf(1,i,nys)
       bff_snd(ii+2)  = uf(2,i,nys)
       bff_snd(ii+3)  = uf(3,i,nys)
       bff_snd(ii+4)  = uf(4,i,nys)
       bff_snd(ii+5)  = uf(5,i,nys)
       bff_snd(ii+6)  = uf(6,i,nys)
       bff_snd(ii+7)  = uf(1,i,nys+1)
       bff_snd(ii+8)  = uf(2,i,nys+1)
       bff_snd(ii+9)  = uf(3,i,nys+1)
       bff_snd(ii+10) = uf(4,i,nys+1)
       bff_snd(ii+11) = uf(5,i,nys+1)
       bff_snd(ii+12) = uf(6,i,nys+1)
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1),12*(nxge-nxgs+5),mnpr,ndown,100, &
                      bff_rcv(1),12*(nxge-nxgs+5),mnpr,nup  ,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i=nxgs-2,nxge+2
       ii = 12*(i-nxgs+2)
       uf(1,i,nye+1) = bff_rcv(ii+1)
       uf(2,i,nye+1) = bff_rcv(ii+2)
       uf(3,i,nye+1) = bff_rcv(ii+3)
       uf(4,i,nye+1) = bff_rcv(ii+4)
       uf(5,i,nye+1) = bff_rcv(ii+5)
       uf(6,i,nye+1) = bff_rcv(ii+6)
       uf(1,i,nye+2) = bff_rcv(ii+7)
       uf(2,i,nye+2) = bff_rcv(ii+8)
       uf(3,i,nye+2) = bff_rcv(ii+9)
       uf(4,i,nye+2) = bff_rcv(ii+10)
       uf(5,i,nye+2) = bff_rcv(ii+11)
       uf(6,i,nye+2) = bff_rcv(ii+12)
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii)
    do i=nxgs-2,nxge+2
       ii = 12*(i-nxgs+2)
       bff_snd(ii+1)  = uf(1,i,nye-1)
       bff_snd(ii+2)  = uf(2,i,nye-1)
       bff_snd(ii+3)  = uf(3,i,nye-1)
       bff_snd(ii+4)  = uf(4,i,nye-1)
       bff_snd(ii+5)  = uf(5,i,nye-1)
       bff_snd(ii+6)  = uf(6,i,nye-1)
       bff_snd(ii+7)  = uf(1,i,nye)
       bff_snd(ii+8)  = uf(2,i,nye)
       bff_snd(ii+9)  = uf(3,i,nye)
       bff_snd(ii+10) = uf(4,i,nye)
       bff_snd(ii+11) = uf(5,i,nye)
       bff_snd(ii+12) = uf(6,i,nye)
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1),12*(nxge-nxgs+5),mnpr,nup  ,110, &
                      bff_rcv(1),12*(nxge-nxgs+5),mnpr,ndown,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxgs-2,nxge+2
       ii = 12*(i-nxgs+2)
       uf(1,i,nys-2) = bff_rcv(ii+1)
       uf(2,i,nys-2) = bff_rcv(ii+2)
       uf(3,i,nys-2) = bff_rcv(ii+3)
       uf(4,i,nys-2) = bff_rcv(ii+4)
       uf(5,i,nys-2) = bff_rcv(ii+5)
       uf(6,i,nys-2) = bff_rcv(ii+6)
       uf(1,i,nys-1) = bff_rcv(ii+7)
       uf(2,i,nys-1) = bff_rcv(ii+8)
       uf(3,i,nys-1) = bff_rcv(ii+9)
       uf(4,i,nys-1) = bff_rcv(ii+10)
       uf(5,i,nys-1) = bff_rcv(ii+11)
       uf(6,i,nys-1) = bff_rcv(ii+12)
    enddo
!$OMP END PARALLEL DO

    call h5util_close_file(file_id)

  end subroutine h5io__input


  subroutine h5io__param(n0,np2,temp,rtemp,fpe,fge,ls,file,nroot)

    integer, intent(in)          :: n0, nroot
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: temp, rtemp, fpe, fge, ls
    character(len=*), intent(in) :: file
    integer :: isp
    real(8) :: pi, vti, vte, va

    integer(hid_t) :: file_id
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling h5io__init()'
       stop
    endif

    pi = 4.0D0*datan(1.0D0)

    vti  = sqrt(2.*temp/r(1))
    vte  = sqrt(2.*temp*rtemp/r(2))
    va   = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    filename = trim(dir)//trim(file)//trim(".h5")
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    call h5util_put_attribute(file_id, "nx", nxge-nxgs+1)
    call h5util_put_attribute(file_id, "ny", nyge-nygs+1)
    call h5util_put_attribute(file_id, "ls", ls)
    call h5util_put_attribute(file_id, "np2", np2(nys,:))
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)
    call h5util_put_attribute(file_id, "fpe", fpe)
    call h5util_put_attribute(file_id, "fge", fge)
    call h5util_put_attribute(file_id, "fpi", fpe*sqrt(r(2)/r(1)))
    call h5util_put_attribute(file_id, "fgi", fge*r(2)/r(1))
    call h5util_put_attribute(file_id, "va", va)
    call h5util_put_attribute(file_id, "vte", vte)
    call h5util_put_attribute(file_id, "vti", vti)
    call h5util_put_attribute(file_id, "beta", (vti/va)**2)
    call h5util_put_attribute(file_id, "rtemp", rtemp)
    call h5util_put_attribute(file_id, "rgi", vti/(fge*r(2)/r(1)))
    call h5util_put_attribute(file_id, "n0", n0)

    call h5util_close_file(file_id)

  end subroutine h5io__param


  subroutine h5io__mom(den,vel,temp,uf,it0)

    integer, intent(in)    :: it0
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,3,nsp)
    integer :: i, j, isp
    real(8) :: tmp(1:6,nxgs:nxge,nys:nye)
    character(len=256) :: filename

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(3), gdim(3), cdim(3), loff(3), goff(3)
    integer :: nd, nxg, nyg, nyl


    if(.not.is_init)then
       write(6,*)'Initialize first by calling h5io__init()'
       stop
    endif

    write(filename,'(a,i7.7,a)')trim(dir),it0,'_mom.h5'
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it0)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !density
    nxg  = nxge-nxgs+1
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/nxg, nyl, nsp/)
    gdim = (/nxg, nyg, nsp/)
    cdim = ldim
    loff = (/0, 0, 0/)
    goff = (/0, nyl*nrank, 0/)
    call h5util_create_dataset_chunk(file_id, "den", H5T_NATIVE_REAL, &
        & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "den", nd, ldim, ldim, loff, goff, &
        & reshape(sngl(den(nxgs:nxge,nys:nye,1:nsp)), (/nxg*nyl*nsp/)))

    !velocity
    !$OMP PARALLEL DO PRIVATE(i,j,isp)
    do isp = 1, nsp
      do j = nys,nye
        do i = nxgs,nxge
          vel(i,j,1,isp) = vel(i,j,1,isp)/den(i,j,isp)
          vel(i,j,2,isp) = vel(i,j,2,isp)/den(i,j,isp)
          vel(i,j,3,isp) = vel(i,j,3,isp)/den(i,j,isp)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call h5util_create_dataset_chunk(file_id, "vx", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "vx", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(vel(nxgs:nxge,nys:nye,1,1:nsp)), (/nxg*nyl*nsp/)))
    call h5util_create_dataset_chunk(file_id, "vy", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "vy", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(vel(nxgs:nxge,nys:nye,2,1:nsp)), (/nxg*nyl*nsp/)))
    call h5util_create_dataset_chunk(file_id, "vz", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "vz", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(vel(nxgs:nxge,nys:nye,3,1:nsp)), (/nxg*nyl*nsp/)))

    !temperature
    !$OMP PARALLEL DO PRIVATE(i,j,isp)
    do isp = 1, nsp
      do j = nys,nye
        do i = nxgs,nxge
          temp(i,j,1,isp) = temp(i,j,1,isp)/den(i,j,isp)
          temp(i,j,2,isp) = temp(i,j,2,isp)/den(i,j,isp)
          temp(i,j,3,isp) = temp(i,j,3,isp)/den(i,j,isp)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call h5util_create_dataset_chunk(file_id, "txx", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "txx", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(temp(nxgs:nxge,nys:nye,1,1:nsp)), (/nxg*nyl*nsp/)))
    call h5util_create_dataset_chunk(file_id, "tyy", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "tyy", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(temp(nxgs:nxge,nys:nye,2,1:nsp)), (/nxg*nyl*nsp/)))
    call h5util_create_dataset_chunk(file_id, "tzz", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "tzz", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(temp(nxgs:nxge,nys:nye,3,1:nsp)), (/nxg*nyl*nsp/)))

    !field
    do j=nys,nye
      do i=nxgs,nxge
        tmp(1,i,j) = 0.5*(uf(1,i,j)+uf(1,i,j+1))
        tmp(2,i,j) = 0.5*(uf(2,i,j)+uf(2,i+1,j))
        tmp(3,i,j) = 0.25*( uf(3,i,j)+uf(3,i+1,j)      &
                           +uf(3,i+1,j)+uf(3,i+1,j+1))
        tmp(4,i,j) = 0.5*(uf(4,i,j)+uf(4,i+1,j))
        tmp(5,i,j) = 0.5*(uf(5,i,j)+uf(5,i,j+1))
        tmp(6,i,j) = uf(6,i,j)
      enddo
    enddo

    nd   = 2
    ldim = (/nxg, nyl, 0/)
    gdim = (/nxg, nyg, 0/)
    cdim = ldim
    loff = (/0, 0, 0/)
    goff = (/0, nyl*nrank, 0/)

    call h5util_create_dataset_chunk(file_id, "bx", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "bx", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(1,nxgs:nxge,nys:nye)), (/nxg*nyl/)))
    call h5util_create_dataset_chunk(file_id, "by", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "by", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(2,nxgs:nxge,nys:nye)), (/nxg*nyl/)))
    call h5util_create_dataset_chunk(file_id, "bz", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "bz", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(3,nxgs:nxge,nys:nye)), (/nxg*nyl/)))
    call h5util_create_dataset_chunk(file_id, "ex", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "ex", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(4,nxgs:nxge,nys:nye)), (/nxg*nyl/)))
    call h5util_create_dataset_chunk(file_id, "ey", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "ey", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(5,nxgs:nxge,nys:nye)), (/nxg*nyl/)))
    call h5util_create_dataset_chunk(file_id, "ez", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "ez", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(6,nxgs:nxge,nys:nye)), (/nxg*nyl/)))

    call h5util_close_file(file_id)

  end subroutine h5io__mom

  subroutine h5io__orb(up,uf,np2,it0)
    use mpi_set

    integer, intent(in) :: it0
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: isp, ii, i, j
    character(len=256) :: filename

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(3), gdim(3), cdim(3), loff(3), goff(3), maxd(3)
    integer :: nd, nxg, nyg, nyl, nrec
    integer :: cnt(0:nproc-1),sum_cnt(0:nproc),cnt_tmp
    real(8) :: tmp(ndim,5*np)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling h5io__init()'
       stop
    endif

    !filename
    write(filename,'(a,i7.7,a)')trim(dir),it0,'_orb.h5'
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it0)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !field
    nxg  = nxge-nxgs+5
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/6, nxg, nyl/)
    gdim = (/6, nxg, nyg/)
    cdim = ldim
    loff = (/0, 0, 0/)
    goff = (/0, 0, nyl*nrank/)

    call h5util_create_dataset_chunk(file_id, "uf", H5T_NATIVE_DOUBLE, &
        & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, &
        & reshape(uf(1:6,nxgs-2:nxge+2,nys:nye), (/6*nxg*nyl/)))

    !particle
    !ion
    isp=1
    cnt_tmp=0
    do j=nys,nye
      do ii=1,np2(j,isp)
         if(up(ndim,ii,j,isp) > 0.D0)then
           tmp(1:ndim,1+cnt_tmp) = up(1:ndim,ii,j,isp)
           cnt_tmp = cnt_tmp+1
         endif
      enddo
    end do

    call MPI_ALLGATHER(cnt_tmp,1,mnpi,cnt,1,mnpi,ncomw,nerr)

    sum_cnt(0) = 0
    do i=1,nproc
      sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
    enddo

    nd   = 2
    ldim = (/ndim-1, cnt(nrank), 0/)
    gdim = (/ndim-1, sum_cnt(nproc), 0/)
    loff = (/0, 0, 0/)
    goff = (/0, sum_cnt(nrank), 0/)
    call h5util_create_dataset(file_id, "upi", H5T_NATIVE_DOUBLE, &
         & nd, gdim)
    call h5util_write_dataset(file_id, "upi", nd, ldim, ldim, loff, goff, &
         & reshape(tmp(1:ndim-1,1:cnt(nrank)), (/(ndim-1)*cnt(nrank)/)))

    nd   = 1
    ldim = (/cnt(nrank), 0, 0/)
    gdim = (/sum_cnt(nproc), 0, 0/)
    loff = (/0, 0, 0/)
    goff = (/sum_cnt(nrank), 0, 0/)
    call h5util_create_dataset(file_id, "idi", H5T_NATIVE_DOUBLE, &
        & nd, gdim)
    call h5util_write_dataset(file_id, "idi", nd, ldim, ldim, loff, goff, &
        & reshape(tmp(ndim,1:cnt(nrank)), (/cnt(nrank)/)))

    !electron
    isp=2
    cnt_tmp=0
    do j=nys,nye
      do ii=1,np2(j,isp)
         if(up(ndim,ii,j,isp) > 0.D0)then
           tmp(1:ndim,1+cnt_tmp) = up(1:ndim,ii,j,isp)
           cnt_tmp = cnt_tmp+1
         endif
      enddo
    end do

    call MPI_ALLGATHER(cnt_tmp,1,mnpi,cnt,1,mnpi,ncomw,nerr)

    sum_cnt(0) = 0
    do i=1,nproc
      sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
    enddo

    nd   = 2
    ldim = (/ndim-1, cnt(nrank), 0/)
    gdim = (/ndim-1, sum_cnt(nproc), 0/)
    loff = (/0, 0, 0/)
    goff = (/0, sum_cnt(nrank), 0/)
    call h5util_create_dataset(file_id, "upe", H5T_NATIVE_DOUBLE, &
         & nd, gdim)
    call h5util_write_dataset(file_id, "upe", nd, ldim, ldim, loff, goff, &
         & reshape(tmp(1:ndim-1,1:cnt(nrank)), (/(ndim-1)*cnt(nrank)/)))

    nd   = 1
    ldim = (/cnt(nrank), 0, 0/)
    gdim = (/sum_cnt(nproc), 0, 0/)
    loff = (/0, 0, 0/)
    goff = (/sum_cnt(nrank), 0, 0/)
    call h5util_create_dataset(file_id, "ide", H5T_NATIVE_DOUBLE, &
        & nd, gdim)
    call h5util_write_dataset(file_id, "ide", nd, ldim, ldim, loff, goff, &
        & reshape(tmp(ndim,1:cnt(nrank)), (/cnt(nrank)/)))


    call h5util_close_file(file_id)

  end subroutine h5io__orb


  !
  ! Parallel HDF5 ultility routines
  !
  ! written by Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
  !

  !
  ! initialize library
  !
  subroutine init()
    implicit none

    call h5open_f(hdferr)

  end subroutine init

  !
  ! finalize library
  !
  subroutine finalize()
    implicit none

    call h5close_f(hdferr)

  end subroutine finalize

  !
  ! create file (replace if already exists)
  !
  subroutine create_file(filename)
    implicit none
    character(len=*), intent(in) :: filename

    integer(hid_t) :: prop, file

    call h5pcreate_f(H5P_FILE_ACCESS_F, prop, hdferr)
    call h5pset_fapl_mpio_f(prop, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr, &
         & H5P_DEFAULT_F, prop)
    call h5pclose_f(prop, hdferr)
    call h5fclose_f(file, hdferr)

  end subroutine create_file

  !
  ! open file
  !
  subroutine open_file(filename, file)
    implicit none
    character(len=*), intent(in) :: filename
    integer(hid_t), intent(out)  :: file

    integer(hid_t) :: prop

    call h5pcreate_f(H5P_FILE_ACCESS_F, prop, hdferr)
    call h5pset_fapl_mpio_f(prop, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file, hdferr, prop)
    call h5pclose_f(prop, hdferr)

  end subroutine open_file

  !
  ! close file
  !
  subroutine close_file(file)
    implicit none
    integer(hid_t), intent(in)  :: file

    call h5fclose_f(file, hdferr)

  end subroutine close_file

  !
  ! put scalar character attribute
  !
  subroutine put_attribute_char(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)
    integer(hsize_t) :: length

    length = len(data)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type, hdferr)
    call h5tset_size_f(type, length, hdferr)

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5tclose_f(type, hdferr)
    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_char

  !
  ! put scalar integer attribute
  !
  subroutine put_attribute_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_int

  !
  ! put scalar real(4) attribute
  !
  subroutine put_attribute_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_REAL

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_r4

  !
  ! put scalar real(8) attribute
  !
  subroutine put_attribute_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(in)          :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_DOUBLE

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_r8

  !
  ! put array integer attribute
  !
  subroutine put_attribute_arr_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_int

  !
  ! put array integer attribute
  !
  subroutine put_attribute_arr_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_REAL
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_r4

  !
  ! put array integer attribute
  !
  subroutine put_attribute_arr_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(in)          :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_DOUBLE
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_r8

  !
  ! create dataset
  !
  subroutine create_dataset(dest, name, type, rank, dims, maxd)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: type
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in), optional :: maxd(rank)

    integer(hsize_t) :: max_dims(rank)

    if( present(maxd) ) then
       max_dims = maxd
    else
       max_dims = dims
    end if

    ! by default, assume chunk == dims
    call create_dataset_chunk(dest, name, type, rank, dims, dims, max_dims)

  end subroutine create_dataset

  !
  ! create dataset with explicitly specified chunk
  !
  subroutine create_dataset_chunk(dest, name, type, rank, dims, cdim, maxd)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: type
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: cdim(rank)
    integer(hsize_t), intent(in), optional :: maxd(rank)

    integer(hid_t) :: space, plist, dtype, dataset

    if( present(maxd) ) then
       call h5screate_simple_f(rank, dims, space, hdferr, maxd)
    else
       call h5screate_simple_f(rank, dims, space, hdferr)
    end if

    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, hdferr)
    call h5tcopy_f(type, dtype, hdferr)

    call h5pset_chunk_f(plist, rank, cdim, hdferr)
    call h5tset_order_f(dtype, H5T_ORDER_LE_F, hdferr)

    call h5dcreate_f(dest, name, dtype, space, dataset, hdferr, &
         & plist, H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5sclose_f(space, hdferr)
    call h5pclose_f(plist, hdferr)
    call h5tclose_f(dtype, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine create_dataset_chunk

  !
  ! write integer global array
  !
  subroutine write_dataset_int(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    integer, intent(in)          :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_int

  !
  ! write real(4) global array
  !
  subroutine write_dataset_r4(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(4), intent(in)          :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_REAL, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_r4

  !
  ! write real(8) global array
  !
  subroutine write_dataset_r8(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(8), intent(in)          :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_r8

  !
  ! extend last dimension by one
  !
  subroutine extend_dimension(dest, name)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name

    integer :: rank
    integer(hid_t) :: dataset, space
    integer(hsize_t) :: dims(MAXDIM), dimm(MAXDIM)

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, space, hdferr)
    call h5sget_simple_extent_ndims_f(space, rank, hdferr)
    call h5sget_simple_extent_dims_f(space, dims, dimm, hdferr)

    dims(rank) = dims(rank) + 1
    call h5dset_extent_f(dataset, dims, hdferr)

    call h5sclose_f(space, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine extend_dimension

  !
  ! get array dimension
  !
  subroutine get_dimension(dest, name, ndim, dims)
    implicit none
    integer(hid_t), intent(in)    :: dest
    character(len=*), intent(in)  :: name
    integer, intent(out)          :: ndim
    integer(hsize_t), intent(out) :: dims(:)

    integer(hid_t) :: dataset, space
    integer(hsize_t) :: dimm(MAXDIM)

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, space, hdferr)
    call h5sget_simple_extent_ndims_f(space, ndim, hdferr)

    if( ndim <= size(dims) ) then
       ! get dimensions
       call h5sget_simple_extent_dims_f(space, dims, dimm, hdferr)
    else
       ! show error message but try to continue
       write(*,*) 'Error in subroutine get_dimension'
    end if

    call h5sclose_f(space, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine get_dimension

  !
  ! read scalar character attribute
  !
  subroutine read_attribute_char(dest, name, data)
    implicit none
    integer(hid_t), intent(in)    :: dest
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_CHARACTER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_char

  !
  ! read scalar integer attribute
  !
  subroutine read_attribute_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(out)         :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_INTEGER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_int

  !
  ! read scalar real(4) attribute
  !
  subroutine read_attribute_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(out)         :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_REAL, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_r4

  !
  ! read scalar real(8) attribute
  !
  subroutine read_attribute_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(out)         :: data

    integer(hid_t) :: attr, type
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_DOUBLE, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_r8

  !
  ! read array integer attribute
  !
  subroutine read_attribute_arr_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(out)         :: data(:)

    integer(hid_t) :: attr, type
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_INTEGER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_arr_int

  !
  ! read array real(4) attribute
  !
  subroutine read_attribute_arr_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(out)          :: data(:)

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_REAL, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_arr_r4

  !
  ! read array real(8) attribute
  !
  subroutine read_attribute_arr_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(out)         :: data(:)

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_DOUBLE, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine read_attribute_arr_r8

  !
  ! read integer global array
  !
  subroutine read_dataset_int(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    integer, intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_INTEGER, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_int

  !
  ! read real(4) global array
  !
  subroutine read_dataset_r4(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(4), intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_REAL, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_r4

  !
  ! read real(8) global array
  !
  subroutine read_dataset_r8(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(8), intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_DOUBLE, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_r8

end module h5io
