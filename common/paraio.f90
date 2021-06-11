module paraio
  use iso_fortran_env, only: int64
  use mpi
  use jsonio
  use mpiio
  implicit none
  private

  public :: paraio__init
  public :: paraio__output
  public :: paraio__input
  public :: paraio__param
  public :: paraio__mom
  public :: paraio__orb


  integer, parameter :: MOK = MPI_OFFSET_KIND ! assume to be the same as int64

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)
  character(len=256), save :: dir

contains

  subroutine paraio__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
       & nproc_in,nrank_in,delx_in,delt_in,c_in,q_in,r_in,dir_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

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

  end subroutine paraio__init


  subroutine paraio__output(up,uf,np2,nxs,nxe,it0,lflag)

    logical, intent(in) :: lflag
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it0
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize, bufsize
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(ndim), gshape(ndim), offset(ndim)
    real(8), allocatable :: buf(:)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(6,*)'Initialize first by calling paraio__init()'
       stop
    endif

    ! filename
    if ( lflag ) then
       write(filename,'(a,i7.7)') trim(dir), 9999999
    else
       write(filename,'(a,i7.7)') trim(dir), it0
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

    call jsonio_put_attribute(json, p, it0, 'it0', disp, '')
    call mpiio_write_atomic(fh, disp, it0)

    call jsonio_put_attribute(json, p, ndim, 'ndim', disp, '')
    call mpiio_write_atomic(fh, disp, ndim)

    call jsonio_put_attribute(json, p, np, 'np', disp, '')
    call mpiio_write_atomic(fh, disp, np)

    call jsonio_put_attribute(json, p, nxs, 'nxs', disp, '')
    call mpiio_write_atomic(fh, disp, nxs)

    call jsonio_put_attribute(json, p, nxe, 'nxe', disp, '')
    call mpiio_write_atomic(fh, disp, nxe)

    call jsonio_put_attribute(json, p, nxgs, 'nxgs', disp, '')
    call mpiio_write_atomic(fh, disp, nxgs)

    call jsonio_put_attribute(json, p, nxge, 'nxge', disp, '')
    call mpiio_write_atomic(fh, disp, nxge)

    call jsonio_put_attribute(json, p, nygs, 'nygs', disp, '')
    call mpiio_write_atomic(fh, disp, nygs)

    call jsonio_put_attribute(json, p, nyge, 'nyge', disp, '')
    call mpiio_write_atomic(fh, disp, nyge)

    call jsonio_put_attribute(json, p, nsp, 'nsp', disp, '')
    call mpiio_write_atomic(fh, disp, nsp)

    call jsonio_put_attribute(json, p, nproc, 'nproc', disp, '')
    call mpiio_write_atomic(fh, disp, nproc)

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

    !
    ! dataset
    !
    bufsize = max(size(up), size(uf))
    allocate(buf(bufsize))

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
    buf(1:lsize) = reshape(up, (/lsize/))
    call jsonio_put_metadata(json, p, 'up', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    nd     = 2
    lshape = (/nyl, nsp, 0, 0/)
    gshape = (/nyg, nsp, 0, 0/)
    offset = (/nyl*nrank, 0, 0/)
    lsize  = product(gshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'number of active particles'
    buf(1:lsize) = reshape(up, (/lsize/))
    call jsonio_put_metadata(json, p, 'np2', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    ! field
    nxg = nxge - nxgs + 5
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'number of active particles'
    buf(1:lsize) = reshape(uf, (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    deallocate(buf)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()
    nullify(p)
    nullify(root)

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__output


  subroutine paraio__input(up,uf,np2,indim,nxs,nxe,it0,filename)

    character(len=*), intent(in) :: filename
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it0, indim
    real(8), intent(out) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, inys, inye, insp, inproc, ibc

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp, dsize, lsize, gsize, bufsize
    integer :: fh, nxg, nyg, nyl
    integer :: nd, lshape(ndim), gshape(ndim), offset(ndim)
    real(8), allocatable :: buf(:)

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(6,*)'Initialize first by calling paraio__init()'
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

    call jsonio_get_metadata(json, p, 'it0', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, it0)

    call jsonio_get_metadata(json, p, 'ndim', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, ndim)

    call jsonio_get_metadata(json, p, 'np', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, np)

    call jsonio_get_metadata(json, p, 'nxs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxs)

    call jsonio_get_metadata(json, p, 'nxe', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxe)

    call jsonio_get_metadata(json, p, 'nxgs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxge)

    call jsonio_get_metadata(json, p, 'nygs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nygs)

    call jsonio_get_metadata(json, p, 'nyge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nyge)

    call jsonio_get_metadata(json, p, 'nsp', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nsp)

    call jsonio_get_metadata(json, p, 'nproc', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nproc)

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

    if((inxgs /= nxgs) .or. (inxge /= nxge) .or.(inygs /= nygs) .or. (inyge /= nyge) &
        .or. (inp /= np) .or. (insp /= nsp) .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    !
    ! dataset
    !
    call json%get(root, 'dataset', p)

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
    call jsonio_get_metadata(json, p, 'up', disp, dsize, nd, gshape)
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, buf)
    up     = reshape(buf(1:lsize), (/ndim, np, nyl, nsp/))

    nd     = 2
    lshape = (/nyl, nsp, 0, 0/)
    gshape = (/nyg, nsp, 0, 0/)
    offset = (/nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    call jsonio_get_metadata(json, p, 'np2', disp, dsize, nd, gshape)
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, buf)
    np2    = reshape(buf(1:lsize), (/nyl, nsp/))

    ! field
    nxg = nxge - nxgs + 5
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    call jsonio_get_metadata(json, p, 'uf', disp, dsize, nd, gshape)
    call mpiio_read_collective(fh, disp, nd, gshape, lshape, offset, buf)
    uf     = reshape(buf(1:lsize), (/6, nxg, nyl/))

    deallocate(buf)

    !
    ! finalize
    !

    ! close json
    call json%destroy()
    call file%destroy()
    nullify(p)
    nullify(root)

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__input


  subroutine paraio__param(n0,np2,temp,rtemp,fpe,fge,ls,filename,nroot)

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
       write(6,*) 'Initialize first by calling paraio__init()'
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
    nullify(p)
    nullify(root)

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__param


  subroutine paraio__mom(den,vel,temp,uf,it0)

    integer, intent(in)    :: it0
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,3,nsp)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize, bufsize
    integer :: i, j, k
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(ndim), gshape(ndim), offset(ndim)
    real(8) :: tmp(nxgs:nxge,nys:nye,1:6)
    real(8), allocatable :: buf(:)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(6,*)'Initialize first by calling paraio__init()'
       stop
    endif

    write(filename,'(a, i7.7, a)') trim(dir), it0, '_mom'
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
    ! dataset
    !
    call json%create_object(p, 'dataset')
    call json%add(root, p)

    bufsize = max(size(uf), size(den), size(vel), size(temp))
    allocate(buf(bufsize))

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
    buf(1:lsize) = reshape(den(nxgs:nxge,nys:nye,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'den', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    ! velocity
    vel(:,:,1,:) = vel(:,:,1,:) / den(:,:,:)
    vel(:,:,2,:) = vel(:,:,2,:) / den(:,:,:)
    vel(:,:,3,:) = vel(:,:,3,:) / den(:,:,:)

    nd    = 4
    lshape = (/nxg, nyl, 3, nsp/)
    gshape = (/nxg, nyg, 3, nsp/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'velocity'
    buf(1:lsize) = reshape(vel(nxgs:nxge,nys:nye,1:3,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'vel', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

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
    buf(1:lsize) = reshape(temp(nxgs:nxge,nys:nye,1:3,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'temp', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    ! electromagnetic field
    do j=nys,nye
       do i=nxgs,nxge
          tmp(i,j,1) = (uf(1,i,j)+uf(1,i,j+1)) / 2
          tmp(i,j,2) = (uf(2,i,j)+uf(2,i+1,j)) / 2
          tmp(i,j,3) = (uf(3,i,j)+uf(3,i+1,j)+uf(3,i+1,j)+uf(3,i+1,j+1)) / 2
          tmp(i,j,4) = (uf(4,i,j)+uf(4,i+1,j)) / 2
          tmp(i,j,5) = (uf(5,i,j)+uf(5,i,j+1)) / 2
          tmp(i,j,6) = uf(6,i,j)
       enddo
    enddo

    nd     = 3
    lshape = (/nxg, nyl, 6, 0/)
    gshape = (/nxg, nyg, 6, 0/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    buf(1:lsize) = reshape(tmp(nxgs:nxge,nys:nye,1:6), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, buf)

    deallocate(buf)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()
    nullify(p)
    nullify(root)

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__mom


  subroutine paraio__orb(up,uf,np2,it0)

    integer, intent(in) :: it0
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: isp, ii, j
    character(len=256) :: filename

    if ( .not. is_init ) then
       write(6,*)'Initialize first by calling paraio__init()'
       stop
    endif

    ! TODO; implement me!

  end subroutine paraio__orb


end module paraio
