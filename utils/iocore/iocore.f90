!
! MPI-IO routines
!
! written by Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
!
module iocore
  use mpi
  implicit none
  private

  ! public routiens
  public :: iocore_get_endian_flag
  public :: iocore_open_file
  public :: iocore_close_file
  public :: iocore_write_atomic
  public :: iocore_write_collective
  public :: iocore_read_atomic
  public :: iocore_read_collective


  integer, parameter :: MOK = MPI_OFFSET_KIND
  integer :: mpierr
  integer :: mpistat(MPI_STATUS_SIZE)

  ! open file
  interface iocore_open_file
     module procedure open_file
  end interface iocore_open_file

  ! close file
  interface iocore_close_file
     module procedure close_file
  end interface iocore_close_file

  ! atomic write
  interface iocore_write_atomic
     module procedure &
          & write_atomic_scalar_char,&
          & write_atomic_scalar_i4, &
          & write_atomic_scalar_i8, &
          & write_atomic_scalar_r4, &
          & write_atomic_scalar_r8, &
          & write_atomic_array_i4, &
          & write_atomic_array_i8, &
          & write_atomic_array_r4, &
          & write_atomic_array_r8, &
          & write_atomic_array_type
  end interface iocore_write_atomic

  ! collective write
  interface iocore_write_collective
     module procedure &
          & write_collective_i4, &
          & write_collective_i8, &
          & write_collective_r4, &
          & write_collective_r8, &
          & write_collective_type
  end interface iocore_write_collective

  ! atomic wread
  interface iocore_read_atomic
     module procedure &
          & read_atomic_scalar_char,&
          & read_atomic_scalar_i4, &
          & read_atomic_scalar_i8, &
          & read_atomic_scalar_r4, &
          & read_atomic_scalar_r8, &
          & read_atomic_array_i4, &
          & read_atomic_array_i8, &
          & read_atomic_array_r4, &
          & read_atomic_array_r8, &
          & read_atomic_array_type
  end interface iocore_read_atomic

  ! collective read
  interface iocore_read_collective
     module procedure &
          & read_collective_i4, &
          & read_collective_i8, &
          & read_collective_r4, &
          & read_collective_r8, &
          & read_collective_type
  end interface iocore_read_collective

contains

  !
  ! get 4 byte integer representing endian (byte order)
  !
  ! On a little endian system, the return value will be 1,
  ! while it will be 16777216 on a big endian system.
  !
  function iocore_get_endian_flag() result(flag)
    integer(4) :: flag

    integer(1) :: byte(4) = (/1_1, 0_1, 0_1, 0_1/)

    flag = transfer(byte, flag)

  end function iocore_get_endian_flag

  !
  ! open file
  !
  subroutine open_file(filename, file, disp, mode)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: file
    integer(MOK), intent(out)    :: disp
    character(len=*), intent(in) :: mode

    integer :: comm = MPI_COMM_WORLD

    select case(mode)
       case('r')
          !
          ! open for read with the pointer at the beginning
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_RDONLY, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_seek(file, disp, MPI_SEEK_SET, mpierr)
       case('w')
          !
          ! open for write with the pointer at the beginning
          ! (truncated if the file already exists)
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_CREATE + MPI_MODE_WRONLY, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_set_size(file, disp, mpierr)
          call MPI_File_seek(file, disp, MPI_SEEK_SET, mpierr)
       case('a')
          !
          ! open for write with the pointer at the end
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_CREATE + MPI_MODE_WRONLY, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_seek(file, disp, MPI_SEEK_END, mpierr)
          call MPI_File_get_position(file, disp, mpierr)
       case('r+')
          !
          ! open for read/write with the pointer at the beginning
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_RDWR, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_seek(file, disp, MPI_SEEK_SET, mpierr)
       case('w+')
          !
          ! open for read/write with the pointer at the beginning
          ! (truncated if the file already exists)
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_CREATE + MPI_MODE_RDWR, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_set_size(file, disp, mpierr)
          call MPI_File_seek(file, disp, MPI_SEEK_SET, mpierr)
       case('a+')
          !
          ! open for read/write with the pointer at the beginning
          ! (truncated if the file already exists)
          !
          call MPI_File_open(comm, filename, &
               & MPI_MODE_CREATE + MPI_MODE_RDWR, &
               & MPI_INFO_NULL, file, mpierr)

          if( mpierr /= 0 ) then
             write(0, &
                  & '("open file ", a, " with mode ", a ," failed !")') &
                  & filename, mode
             call MPI_Finalize(mpierr)
             stop
          end if

          ! initialize
          disp = 0
          call MPI_File_seek(file, disp, MPI_SEEK_END, mpierr)
          call MPI_File_get_position(file, disp, mpierr)
       case default
          write(0, &
               & '("unknown mode ", a, " specifield at open_file for ", a)') &
               & mode, filename
          call MPI_Finalize(mpierr)
          stop
    end select

  end subroutine open_file

  !
  ! close file
  !
  subroutine close_file(file)
    implicit none
    integer, intent(inout) :: file

    call MPI_File_close(file, mpierr)

  end subroutine close_file

  !
  ! generic routine for atomic writing
  !
  subroutine write_atomic_array_type(file, disp, data, byte)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: data(:)
    integer, intent(in)         :: byte

    integer :: rank

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)

    if( rank == 0 ) then
       call MPI_File_write_at(file, disp, data, byte, MPI_BYTE, mpistat, mpierr)
    end if

    disp = disp + byte

  end subroutine write_atomic_array_type

  !
  ! routine for atomic writing of character scalar
  !
  subroutine write_atomic_scalar_char(file, disp, data)
    implicit none
    integer, intent(in)          :: file
    integer(MOK), intent(inout)  :: disp
    character(len=*), intent(in) :: data

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 1*len_trim(data))

  end subroutine write_atomic_scalar_char

  !
  ! routine for atomic writing of integer(4) scalar
  !
  subroutine write_atomic_scalar_i4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(4), intent(in)      :: data

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 4)

  end subroutine write_atomic_scalar_i4

  !
  ! routine for atomic writing of integer(8) scalar
  !
  subroutine write_atomic_scalar_i8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(8), intent(in)      :: data

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 8)

  end subroutine write_atomic_scalar_i8

  !
  ! routine for atomic writing of real(4) scalar
  !
  subroutine write_atomic_scalar_r4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(4), intent(in)         :: data

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 4)

  end subroutine write_atomic_scalar_r4

  !
  ! routine for atomic writing of real(8) scalar
  !
  subroutine write_atomic_scalar_r8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(8), intent(in)         :: data

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 8)

  end subroutine write_atomic_scalar_r8

  !
  ! routine for atomic writing of integer(4) array
  !
  subroutine write_atomic_array_i4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(4), intent(in)      :: data(:)

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 4*size(data))

  end subroutine write_atomic_array_i4

  !
  ! routine for atomic writing of integer(8) array
  !
  subroutine write_atomic_array_i8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(8), intent(in)      :: data(:)

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 8*size(data))

  end subroutine write_atomic_array_i8

  !
  ! routine for atomic writing of real(4) array
  !
  subroutine write_atomic_array_r4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(4), intent(in)         :: data(:)

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 4*size(data))

  end subroutine write_atomic_array_r4

  !
  ! routine for atomic writing of real(8) array
  !
  subroutine write_atomic_array_r8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(8), intent(in)         :: data(:)

    integer :: dummy_type(1)

    call write_atomic_array_type(file, disp, &
         & transfer(data, dummy_type), 8*size(data))

  end subroutine write_atomic_array_r8

  !
  ! generic routine for collective writing of array
  !
  subroutine write_collective_type(file, disp, rank, gshape, lshape, offset, data, &
       & mpitype)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    integer, intent(in)         :: data(:)
    integer, intent(in)         :: mpitype

    integer :: filetype

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)

    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_write_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)

    call MPI_Type_free(filetype, mpierr)

  end subroutine write_collective_type

  !
  ! routine for collective writing of integer(4) array
  !
  subroutine write_collective_i4(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    integer(4), intent(in)      :: data(:)

    integer :: dummy_type(1)

    call write_collective_type(file, disp, rank, gshape, lshape, offset, &
         & transfer(data, dummy_type), MPI_INTEGER4)
    disp = disp + 4*product(gshape)

  end subroutine write_collective_i4

  !
  ! routine for collective writing of integer(8) array
  !
  subroutine write_collective_i8(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    integer(8), intent(in)      :: data(:)

    integer :: dummy_type(1)

    call write_collective_type(file, disp, rank, gshape, lshape, offset, &
         & transfer(data, dummy_type), MPI_INTEGER8)
    disp = disp + 8*product(gshape)

  end subroutine write_collective_i8

  !
  ! routine for collective writing of real(4) array
  !
  subroutine write_collective_r4(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    real(4), intent(in)         :: data(:)

    integer :: dummy_type(1)

    call write_collective_type(file, disp, rank, gshape, lshape, offset, &
         & transfer(data, dummy_type), MPI_REAL4)
    disp = disp + 4*product(gshape)

  end subroutine write_collective_r4

  !
  ! routine for collective writing of real(8) array
  !
  subroutine write_collective_r8(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    real(8), intent(in)         :: data(:)

    integer :: dummy_type(1)

    call write_collective_type(file, disp, rank, gshape, lshape, offset, &
         & transfer(data, dummy_type), MPI_REAL8)
    disp = disp + 8*product(gshape)

  end subroutine write_collective_r8

  !
  ! generic routine for atomic rading of array
  !
  subroutine read_atomic_array_type(file, disp, data, byte)
    implicit none
    integer, intent(in)          :: file
    integer(MOK), intent(inout)  :: disp
    integer, pointer, intent(in) :: data
    integer, intent(in)          :: byte

    call MPI_File_read_at(file, disp, data, byte, MPI_BYTE, mpistat, mpierr)
    disp = disp + byte

  end subroutine read_atomic_array_type

  !
  ! routine for atomic rading of character scalar
  !
  subroutine read_atomic_scalar_char(file, disp, data, length)
    implicit none
    integer, intent(in)             :: file
    integer(MOK), intent(inout)     :: disp
    character(len=*), intent(inout) :: data
    integer, optional, intent(in)   :: length

    integer :: byte

    if ( .not. present(length) ) then
       byte = len(data)
    else
       byte = length
    end if

    data = ''
    call MPI_File_read_at(file, disp, data, byte, MPI_BYTE, mpistat, mpierr)
    disp = disp + byte

  end subroutine read_atomic_scalar_char

  !
  ! routine for atomic rading of integer(4) scalar
  !
  subroutine read_atomic_scalar_i4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(4), intent(inout)   :: data

    call MPI_File_read_at(file, disp, data, 4, MPI_BYTE, mpistat, mpierr)
    disp = disp + 4

  end subroutine read_atomic_scalar_i4

  !
  ! routine for atomic rading of integer(8) scalar
  !
  subroutine read_atomic_scalar_i8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer(8), intent(inout)   :: data

    call MPI_File_read_at(file, disp, data, 8, MPI_BYTE, mpistat, mpierr)
    disp = disp + 8

  end subroutine read_atomic_scalar_i8

  !
  ! routine for atomic rading of real(4) scalar
  !
  subroutine read_atomic_scalar_r4(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(4), intent(inout)      :: data

    call MPI_File_read_at(file, disp, data, 4, MPI_BYTE, mpistat, mpierr)
    disp = disp + 4

  end subroutine read_atomic_scalar_r4

  !
  ! routine for atomic rading of real(8) scalar
  !
  subroutine read_atomic_scalar_r8(file, disp, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    real(8), intent(inout)      :: data

    call MPI_File_read_at(file, disp, data, 8, MPI_BYTE, mpistat, mpierr)
    disp = disp + 8

  end subroutine read_atomic_scalar_r8

  !
  ! routine for atomic rading of integer(4) array
  !
  subroutine read_atomic_array_i4(file, disp, data)
    implicit none
    integer, intent(in)               :: file
    integer(MOK), intent(inout)       :: disp
    integer(4), target, intent(inout) :: data(:)

    integer, pointer :: dummy(:)

    dummy => data
    call read_atomic_array_type(file, disp, dummy(1), 4*size(data))

  end subroutine read_atomic_array_i4

  !
  ! routine for atomic rading of integer(8) array
  !
  subroutine read_atomic_array_i8(file, disp, data)
    implicit none
    integer, intent(in)            :: file
    integer(MOK), intent(inout)    :: disp
    integer(8), target, intent(in) :: data(:)

    integer, pointer :: dummy(:)
    integer, target  :: dummy_data(1)

    dummy_data = transfer(data, dummy_data)
    dummy => dummy_data
    call read_atomic_array_type(file, disp, dummy(1), 8*size(data))

  end subroutine read_atomic_array_i8

  !
  ! routine for atomic rading of real(4) array
  !
  subroutine read_atomic_array_r4(file, disp, data)
    implicit none
    integer, intent(in)            :: file
    integer(MOK), intent(inout)    :: disp
    real(4), target, intent(in)    :: data(:)

    integer, pointer :: dummy(:)
    integer, target  :: dummy_data(1)

    dummy_data = transfer(data, dummy_data)
    dummy => dummy_data
    call read_atomic_array_type(file, disp, dummy(1), 4*size(data))

  end subroutine read_atomic_array_r4

  !
  ! routine for atomic rading of real(8) array
  !
  subroutine read_atomic_array_r8(file, disp, data)
    implicit none
    integer, intent(in)            :: file
    integer(MOK), intent(inout)    :: disp
    real(8), target, intent(in)    :: data(:)

    integer, pointer :: dummy(:)
    integer, target  :: dummy_data(1)

    dummy_data = transfer(data, dummy_data)
    dummy => dummy_data
    call read_atomic_array_type(file, disp, dummy(1), 8*size(data))

  end subroutine read_atomic_array_r8

  !
  ! generic routine for collective reading of array
  !
  ! <<< Note >>>
  ! This routine does not work for some reasons if called from other routines.
  ! Consequently, similar routines with other data types are hard coded.
  !
  subroutine read_collective_type(file, disp, rank, gshape, lshape, offset, data, &
       & mpitype)
    implicit none
    integer, intent(in)          :: file
    integer(MOK), intent(inout)  :: disp
    integer, intent(in)          :: rank
    integer, intent(in)          :: gshape(rank)
    integer, intent(in)          :: lshape(rank)
    integer, intent(in)          :: offset(rank)
    integer, intent(in)          :: data(:)
    integer, intent(in)          :: mpitype

    integer :: filetype

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)

    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_read_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)

    call MPI_Type_free(filetype, mpierr)

  end subroutine read_collective_type

  !
  ! collective reading of integer(4) array
  !
  subroutine read_collective_i4(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)               :: file
    integer(MOK), intent(inout)       :: disp
    integer, intent(in)               :: rank
    integer, intent(in)               :: gshape(rank)
    integer, intent(in)               :: lshape(rank)
    integer, intent(in)               :: offset(rank)
    integer(4), intent(in)            :: data(:)

    integer :: filetype
    integer :: mpitype = MPI_INTEGER4

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)
    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_read_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)
    call MPI_Type_free(filetype, mpierr)

    disp = disp + 4*product(gshape)

  end subroutine read_collective_i4

  !
  ! collective reading of integer(8) array
  !
  subroutine read_collective_i8(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)               :: file
    integer(MOK), intent(inout)       :: disp
    integer, intent(in)               :: rank
    integer, intent(in)               :: gshape(rank)
    integer, intent(in)               :: lshape(rank)
    integer, intent(in)               :: offset(rank)
    integer(8), intent(in)            :: data(:)

    integer :: filetype
    integer :: mpitype = MPI_INTEGER8

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)
    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_read_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)
    call MPI_Type_free(filetype, mpierr)

    disp = disp + 8*product(gshape)

  end subroutine read_collective_i8

  !
  ! collective reading of real(4) array
  !
  subroutine read_collective_r4(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    real(4), intent(in)         :: data(:)

    integer :: filetype
    integer :: mpitype = MPI_REAL4

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)
    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_read_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)
    call MPI_Type_free(filetype, mpierr)

    disp = disp + 4*product(gshape)

  end subroutine read_collective_r4

  !
  ! collective reading of real(8) array
  !
  subroutine read_collective_r8(file, disp, rank, gshape, lshape, offset, data)
    implicit none
    integer, intent(in)         :: file
    integer(MOK), intent(inout) :: disp
    integer, intent(in)         :: rank
    integer, intent(in)         :: gshape(rank)
    integer, intent(in)         :: lshape(rank)
    integer, intent(in)         :: offset(rank)
    real(8), intent(in)         :: data(:)

    integer :: filetype
    integer :: mpitype = MPI_REAL8

    call MPI_Type_create_subarray(rank, gshape, lshape, offset, &
         & MPI_ORDER_FORTRAN, mpitype, filetype, mpierr)
    call MPI_Type_commit(filetype, mpierr)
    call MPI_File_set_view(file, disp, mpitype, filetype, "native", &
         & MPI_INFO_NULL, mpierr)
    call MPI_File_read_all(file, data, product(lshape), mpitype, &
         & mpistat, mpierr)
    call MPI_Type_free(filetype, mpierr)

    disp = disp + 8*product(gshape)

  end subroutine read_collective_r8

end module iocore
