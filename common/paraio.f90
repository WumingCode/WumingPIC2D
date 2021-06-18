module paraio
  use iso_fortran_env, only: int64
  use mpi
  use jsonio
  use mpiio
  implicit none
  private

  public :: paraio__init
  public :: paraio__finalize
  public :: paraio__output
  public :: paraio__input
  public :: paraio__param
  public :: paraio__mom
  public :: paraio__ptcl
  public :: paraio__orb


  integer, parameter :: MOK = MPI_OFFSET_KIND ! assume to be the same as int64

  character(len=256), save :: dir
  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)

  integer :: mpierr
  real(8), allocatable :: mpibuf(:)

contains

  !
  ! initialize module
  !
  subroutine paraio__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
       & nproc_in,nrank_in,delx_in,delt_in,c_in,q_in,r_in,dir_in)
    implicit none
    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    integer :: psize, fsize

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

    ! allocate MPI buffer
    psize = ndim*np*(nye-nys+5)*nsp
    fsize = 6*(nxge-nxgs+5)*(nye-nys+5)
    allocate(mpibuf(max(psize, fsize)))

    is_init = .true.

  end subroutine paraio__init

  !
  ! finalize module
  !
  subroutine paraio__finalize()
    implicit none

    deallocate(mpibuf)

  end subroutine paraio__finalize

  !
  ! output data for re-calculation
  !
  subroutine paraio__output(up,uf,np2,nxs,nxe,it,lflag)
    implicit none
    logical, intent(in) :: lflag
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(4), gshape(4), offset(4)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    ! filename
    if ( lflag ) then
       write(filename,'(a,i7.7)') trim(dir), 9999999
    else
       write(filename,'(a,i7.7)') trim(dir), it
    endif

    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'w')

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! attribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call jsonio_put_attribute(json, p, nxs, 'nxs', disp, '')
    call mpiio_write_atomic(fh, disp, nxs)

    call jsonio_put_attribute(json, p, nxe, 'nxe', disp, '')
    call mpiio_write_atomic(fh, disp, nxe)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !

    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! particle
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd     = 4
    lshape = (/ndim, np, nyl, nsp/)
    gshape = (/ndim, np, nyg, nsp/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'particles'
    mpibuf(1:lsize) = reshape(up, (/lsize/))
    call jsonio_put_metadata(json, p, 'up', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    nd     = 2
    lshape = (/nyl, nsp, 0, 0/)
    gshape = (/nyg, nsp, 0, 0/)
    offset = (/nyl*nrank, 0, 0, 0/)
    lsize  = product(gshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'number of active particles'
    mpibuf(1:lsize) = reshape(up, (/lsize/))
    call jsonio_put_metadata(json, p, 'np2', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    ! field: including ghost cells
    nxg = size(uf, 2) ! nxge - nxgs + 5
    nyl = size(uf, 3) ! nye  - nys  + 5
    nyg = nyl * nproc

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic fields including ghost cells'
    mpibuf(1:lsize) = reshape(uf, (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__output

  !
  ! input data for re-calculation
  !
  subroutine paraio__input(up,uf,np2,indim,nxs,nxe,it,filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it, indim
    real(8), intent(out) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, inys, inye, insp, inproc, ibc

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp, dsize, lsize, gsize
    integer :: fh, nxg, nyg, nyl
    integer :: nd, lshape(4), gshape(4), offset(4)

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'r')

    ! json file
    call json%initialize()
    call file%initialize()
    call file%load(jsonfile)
    call file%get(root)

    !
    ! atttribute
    !
    call json%get(root, 'attribute', p)

    call jsonio_get_metadata(json, p, 'it', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, it)

    call jsonio_get_metadata(json, p, 'nxs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxs)

    call jsonio_get_metadata(json, p, 'nxe', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxe)

    call jsonio_get_metadata(json, p, 'ndim', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, indim)

    call jsonio_get_metadata(json, p, 'np', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inp)

    call jsonio_get_metadata(json, p, 'nxgs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inxgs)

    call jsonio_get_metadata(json, p, 'nxge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inxge)

    call jsonio_get_metadata(json, p, 'nygs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inygs)

    call jsonio_get_metadata(json, p, 'nyge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inyge)

    call jsonio_get_metadata(json, p, 'nsp', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, insp)

    call jsonio_get_metadata(json, p, 'nproc', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inproc)

    call jsonio_get_metadata(json, p, 'delx', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, delx)

    call jsonio_get_metadata(json, p, 'delt', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, delt)

    call jsonio_get_metadata(json, p, 'c', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, c)

    call jsonio_get_metadata(json, p, 'r', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, r)

    call jsonio_get_metadata(json, p, 'q', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, q)

    if( .not. ( &
         & inxgs == nxgs .and. &
         & inxge == nxge .and. &
         & inygs == nygs .and. &
         & inyge == nyge .and. &
         & inp == np .and. &
         & insp == nsp .and. &
         & inproc == nproc) ) then
       write(0,*) '*** Error: invalid input parameters ***'
       stop
    endif

    !
    ! dataset
    !
    call json%get(root, 'dataset', p)

    ! * particle
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'up', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 4 .and. &
         &   gshape(1) == ndim .and. &
         &   gshape(2) == np   .and. &
         &   gshape(3) == nyg  .and. &
         &   gshape(4) == nsp ) ) then
       write(0, *) 'Fatail error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lshape = (/ndim, np, nyl, nsp/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:4))
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)
    up = reshape(mpibuf(1:lsize), lshape(1:4))

    ! * particle number

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'np2', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 2 .and. &
         &   gshape(1) == nyg  .and. &
         &   gshape(2) == nsp ) ) then
       write(0, *) 'Fatail error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lshape = (/nyl, nsp, 0, 0/)
    offset = (/nyl*nrank, 0, 0, 0/)
    lsize  = product(lshape(1:nd))
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)
    np2 = reshape(mpibuf(1:lsize), lshape(1:2))

    ! * field
    nxg = size(uf, 2) ! nxge - nxgs + 5
    nyl = size(uf, 3) ! nye  - nys  + 5
    nyg = nyl * nproc

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'uf', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 3 .and. &
         &   gshape(1) == 6  .and. &
         &   gshape(2) == nxg  .and. &
         &   gshape(3) == nyg ) ) then
       write(0, *) 'Fatail error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lshape = (/6, nxg, nyl, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)
    uf = reshape(mpibuf(1:lsize), lshape(1:3))

    !
    ! finalize
    !

    ! close json
    call json%destroy()
    call file%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__input

  !
  ! output parameters
  !
  subroutine paraio__param(n0,np2,temp,rtemp,fpe,fge,ls,filename,nroot)
    implicit none
    integer, intent(in)          :: n0, nroot
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: temp, rtemp, fpe, fge, ls
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp
    integer :: endian, fh, nx, ny
    real(8) :: pi, vti, vte, vai, vae, fpi, fgi

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    pi   = 4*atan(1.0D0)
    vti  = sqrt(2*temp/r(1))
    vte  = sqrt(2*temp*rtemp/r(2))
    vai  = abs(fge*r(1)*c/q(1))/sqrt(4*pi*r(1)*n0) * r(2)/r(1)
    vae  = abs(fge*r(2)*c/q(2))/sqrt(4*pi*r(2)*n0)
    fpi  = fpe * sqrt(r(2)/r(1))
    fgi  = fge * r(2)/r(1)

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'w')

    !
    ! attribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    nx = nxge - nxgs + 1
    call jsonio_put_attribute(json, p, nx, 'nx', disp, '')
    call mpiio_write_atomic(fh, disp, nx)

    ny = nyge - nygs + 1
    call jsonio_put_attribute(json, p, ny, 'ny', disp, '')
    call mpiio_write_atomic(fh, disp, ny)

    call jsonio_put_attribute(json, p, ls, 'ls', disp, '')
    call mpiio_write_atomic(fh, disp, ls)

    call jsonio_put_attribute(json, p, np2(nys,:), 'np2', disp, '')
    call mpiio_write_atomic(fh, disp, np2(nys,:))

    call jsonio_put_attribute(json, p, np, 'np', disp, '')
    call mpiio_write_atomic(fh, disp, np)

    call jsonio_put_attribute(json, p, delx, 'delx', disp, '')
    call mpiio_write_atomic(fh, disp, delx)

    call jsonio_put_attribute(json, p, delt, 'delt', disp, '')
    call mpiio_write_atomic(fh, disp, delt)

    call jsonio_put_attribute(json, p, c, 'c', disp, '')
    call mpiio_write_atomic(fh, disp, c)

    call jsonio_put_attribute(json, p, r, 'r', disp, '')
    call mpiio_write_atomic(fh, disp, r)

    call jsonio_put_attribute(json, p, q, 'q', disp, '')
    call mpiio_write_atomic(fh, disp, q)

    call jsonio_put_attribute(json, p, fpe, 'fpe', disp, '')
    call mpiio_write_atomic(fh, disp, fpe)

    call jsonio_put_attribute(json, p, fge, 'fge', disp, '')
    call mpiio_write_atomic(fh, disp, fge)

    call jsonio_put_attribute(json, p, fpi, 'fpi', disp, '')
    call mpiio_write_atomic(fh, disp, fpi)

    call jsonio_put_attribute(json, p, fgi, 'fgi', disp, '')
    call mpiio_write_atomic(fh, disp, fgi)

    call jsonio_put_attribute(json, p, vai, 'vai', disp, '')
    call mpiio_write_atomic(fh, disp, vai)

    call jsonio_put_attribute(json, p, vae, 'vae', disp, '')
    call mpiio_write_atomic(fh, disp, vae)

    call jsonio_put_attribute(json, p, vte, 'vte', disp, '')
    call mpiio_write_atomic(fh, disp, vte)

    call jsonio_put_attribute(json, p, vti, 'vti', disp, '')
    call mpiio_write_atomic(fh, disp, vti)

    call jsonio_put_attribute(json, p, vte, 'vte', disp, '')
    call mpiio_write_atomic(fh, disp, vte)

    call jsonio_put_attribute(json, p, n0, 'n0', disp, '')
    call mpiio_write_atomic(fh, disp, n0)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__param

  !
  ! output field and moment quantities
  !
  subroutine paraio__mom(den,vel,temp,uf,it)
    implicit none
    integer, intent(in)    :: it
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,3,nsp)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize
    integer :: i, j, k
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(4), gshape(4), offset(4)
    real(8) :: tmp(1:6,nxgs:nxge,nys:nye)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    write(filename,'(a, i7.7, a)') trim(dir), it, '_mom'
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'w')

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! atttribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !
    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! density
    nxg  = nxge - nxgs + 1
    nyg  = nyge - nygs + 1
    nyl  = nye  - nys  + 1

    ! density
    nd    = 3
    lshape = (/nxg, nyl, nsp, 0/)
    gshape = (/nxg, nyg, nsp, 0/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'density'
    mpibuf(1:lsize) = reshape(den(nxgs:nxge,nys:nye,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'den', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    ! velocity
    vel(:,:,1,:) = vel(:,:,1,:) / den(:,:,:)
    vel(:,:,2,:) = vel(:,:,2,:) / den(:,:,:)
    vel(:,:,3,:) = vel(:,:,3,:) / den(:,:,:)

    nd     = 4
    lshape = (/nxg, nyl, 3, nsp/)
    gshape = (/nxg, nyg, 3, nsp/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'velocity'
    mpibuf(1:lsize) = reshape(vel(nxgs:nxge,nys:nye,1:3,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'vel', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    ! temperature
    temp(:,:,1,:) = temp(:,:,1,:) / den(:,:,:)
    temp(:,:,2,:) = temp(:,:,2,:) / den(:,:,:)
    temp(:,:,3,:) = temp(:,:,3,:) / den(:,:,:)

    nd    = 4
    lshape = (/nxg, nyl, 3, nsp/)
    gshape = (/nxg, nyg, 3, nsp/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'temperature'
    mpibuf(1:lsize) = reshape(temp(nxgs:nxge,nys:nye,1:3,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'temp', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    ! electromagnetic field on cell centers
    do j=nys,nye
       do i=nxgs,nxge
          tmp(1,i,j) = (uf(1,i,j)+uf(1,i,j+1)) / 2
          tmp(2,i,j) = (uf(2,i,j)+uf(2,i+1,j)) / 2
          tmp(3,i,j) = (uf(3,i,j)+uf(3,i+1,j)+uf(3,i+1,j)+uf(3,i+1,j+1)) / 2
          tmp(4,i,j) = (uf(4,i,j)+uf(4,i+1,j)) / 2
          tmp(5,i,j) = (uf(5,i,j)+uf(5,i,j+1)) / 2
          tmp(6,i,j) = uf(6,i,j)
       enddo
    enddo

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf(1:lsize) = reshape(tmp(1:6,nxgs:nxge,nys:nye), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__mom

  !
  ! output all the active particles
  !
  subroutine paraio__ptcl(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 0, '_ptcl')

  end subroutine paraio__ptcl

  !
  ! output tracer particles
  !
  subroutine paraio__orb(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 1, '_orb')

  end subroutine paraio__orb

  !
  ! output particles
  !
  subroutine write_particle(up,uf,np2,it,mode,suffix)
    implicit none
    integer, intent(in)          :: it
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer, intent(in)          :: mode
    character(len=*), intent(in) :: suffix

    character(len=256) :: filename, jsonfile, datafile, desc, name
    integer(int64) :: disp, dsize, lsize, gsize, goffset
    integer(int64) :: cumsum(nproc+1,nsp), ip1, ip2
    integer :: i, j, k, isp, irank
    integer :: fh, endian, nxg, nyg, nyl, npl, npg, npo
    integer :: nd, lshape(4), gshape(4), offset(4)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    ! filename
    write(filename,'(a, i7.7, a)') trim(dir), it, trim(suffix)
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'w')

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! atttribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !
    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! field
    nxg = nxge - nxgs + 1
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf(1:lsize) = reshape(uf(1:6,nxgs:nxge,nys:nye), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf)

    ! particle
    call get_particle_count(up, np2, mpibuf, cumsum, mode)

    irank = nrank + 1
    ip2   = 0
    do isp = 1, nsp
       write(desc, '("particle species #", i2.2)') isp
       write(name, '("up", i2.2)') isp

       npl     = cumsum(irank+1,isp) - cumsum(irank,isp)
       npg     = cumsum(nproc+1,isp)
       npo     = cumsum(irank,isp)
       ip1     = ip2 + 1
       ip2     = ip1 + ndim*npl - 1

       nd      = 2
       gshape  = (/ndim, npg, 0, 0/)
       lsize   = ndim*npl
       gsize   = ndim*npg
       goffset = ndim*npo
       dsize   = gsize * 8
       call jsonio_put_metadata(json, p, trim(name), 'f8', disp, &
            & dsize, nd, gshape, desc)
       call mpiio_write_collective(fh, disp, gsize, lsize, goffset, &
            & mpibuf(ip1:ip2))
    end do

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine write_particle

  !
  ! get particles distribution and pack to buffer
  !
  subroutine get_particle_count(up, np2, buf, cumsum, mode)
    implicit none
    integer, intent(in)       :: np2(nys:nye,nsp)
    real(8), intent(in)       :: up(ndim,np,nys:nye,nsp)
    real(8), intent(inout)    :: buf(:)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)
    integer, intent(in)       :: mode

    integer :: i, j, ip, jp, isp
    integer(8) :: lcount(nsp), gcount(nsp, nproc)
    integer(8) :: pid

    ! count number of particles and pack into buffer
    if ( mode == 0 ) then
       ! * mode 0: all the active particles
       ip = 1
       do isp = 1, nsp
          lcount(isp) = 0
          do j = nys, nye
             do i = 1, np2(j,isp)
                lcount(isp) = lcount(isp) + 1
                ! packing
                do jp = 1, ndim
                   buf(ip) = up(jp,i,j,isp)
                   ip = ip + 1
                end do
             enddo
          end do
       end do

    else if ( mode == 1 ) then
       ! * mode 1: tracer particles with positive IDs
       ip = 1
       do isp = 1, nsp
          lcount(isp) = 0
          do j = nys, nye
             do i = 1, np2(j,isp)
                ! get particle ID as 64bit integer
                pid = transfer(up(ndim,i,j,isp), 1_8)

                ! count positive
                if( pid > 0 ) then
                   lcount(isp) = lcount(isp) + 1
                   ! packing
                   do jp = 1, ndim
                      buf(ip) = up(jp,i,j,isp)
                      ip = ip + 1
                   end do
                endif
             enddo
          end do
       end do

    else
       ! error
       write(0,*) 'Error: invalid mode specified for get_particle_count'
       call MPI_Finalize(mpierr)
       stop
    end if

    call MPI_Allgather(lcount, nsp, MPI_INTEGER8, gcount, nsp, MPI_INTEGER8, &
         & MPI_COMM_WORLD, mpierr)

    ! calculate cumulative sum
    do isp = 1, nsp
       cumsum(1,isp) = 0
       do i = 1, nproc
          cumsum(i+1,isp) = cumsum(i,isp) + gcount(isp,i)
       end do
    end do

  end subroutine get_particle_count

  !
  ! put common metadata
  !
  subroutine put_metadata(json, p, file, disp)
    implicit none
    type(json_core), intent(inout)        :: json
    type(json_value), pointer, intent(in) :: p
    integer, intent(in)                   :: file
    integer(int64), intent(inout)         :: disp

    call jsonio_put_attribute(json, p, ndim, 'ndim', disp, '')
    call mpiio_write_atomic(file, disp, ndim)

    call jsonio_put_attribute(json, p, np, 'np', disp, '')
    call mpiio_write_atomic(file, disp, np)

    call jsonio_put_attribute(json, p, nxgs, 'nxgs', disp, '')
    call mpiio_write_atomic(file, disp, nxgs)

    call jsonio_put_attribute(json, p, nxge, 'nxge', disp, '')
    call mpiio_write_atomic(file, disp, nxge)

    call jsonio_put_attribute(json, p, nygs, 'nygs', disp, '')
    call mpiio_write_atomic(file, disp, nygs)

    call jsonio_put_attribute(json, p, nyge, 'nyge', disp, '')
    call mpiio_write_atomic(file, disp, nyge)

    call jsonio_put_attribute(json, p, nsp, 'nsp', disp, '')
    call mpiio_write_atomic(file, disp, nsp)

    call jsonio_put_attribute(json, p, nproc, 'nproc', disp, '')
    call mpiio_write_atomic(file, disp, nproc)

    call jsonio_put_attribute(json, p, delx, 'delx', disp, '')
    call mpiio_write_atomic(file, disp, delx)

    call jsonio_put_attribute(json, p, delt, 'delt', disp, '')
    call mpiio_write_atomic(file, disp, delt)

    call jsonio_put_attribute(json, p, c, 'c', disp, '')
    call mpiio_write_atomic(file, disp, c)

    call jsonio_put_attribute(json, p, r, 'r', disp, '')
    call mpiio_write_atomic(file, disp, r)

    call jsonio_put_attribute(json, p, q, 'q', disp, '')
    call mpiio_write_atomic(file, disp, q)

  end subroutine put_metadata

end module paraio
