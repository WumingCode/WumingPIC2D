module mpi_set

  implicit none
  private
  public :: mpi_set__init, MPI_WTIME

  include 'mpif.h'

  integer, public, parameter   :: mnpi  = MPI_INTEGER
  integer, public, parameter   :: mnpr  = MPI_DOUBLE_PRECISION
  integer, public, parameter   :: mnpc  = MPI_CHARACTER
  integer, public, parameter   :: opsum = MPI_SUM
  integer, public              :: nerr, ncomw, nsize, nrank
  integer, public              :: nstat(MPI_STATUS_SIZE)
  integer, public              :: nxs, nxe, nxs1, nxe1, nys, nye, nys1, nye1, nup, ndown


contains


  subroutine mpi_set__init(nxgs,nxge,nygs,nyge,nproc)

    integer, intent(in) :: nxgs, nxge, nygs, nyge, nproc
    integer             :: iwork1, iwork2

    !*********** Initialization for MPI  ***************!
    call MPI_INIT(nerr)
    ncomw = MPI_COMM_WORLD
    call MPI_COMM_SIZE(ncomw, nsize, nerr)
    call MPI_COMM_RANK(ncomw, nrank, nerr)

    if(nsize /= nproc) then
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
       stop 'error in proc no.'
    endif

    !start and end of loop counter
    nxs = nxgs
    nxe = nxge
    iwork1 = (nyge-nygs+1)/nsize
    iwork2 = mod(nyge-nygs+1,nsize)
    nys = nrank*iwork1+nygs+min(nrank,iwork2)
    nye = nys+iwork1-1
    if(iwork2 > nrank) nye = nye+1

    !For MPI_SENDRECV
    nxs1 = nxs-1
    nxe1 = nxe+1
    nys1 = nys-1
    nye1 = nye+1
    nup   = nrank+1
    ndown = nrank-1
    if(nrank == nsize-1) nup = 0    ! periodic boundary condition
    if(nrank == 0) ndown = nsize-1  ! periodic boundary condition

  end subroutine mpi_set__init


end module mpi_set

