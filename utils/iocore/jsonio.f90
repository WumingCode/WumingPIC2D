!
! JSON IO routines
!
! written by Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
!
module jsonio
  use json_module, IK => json_IK, RK => json_RK
  use iso_fortran_env
  implicit none
  private

  public :: jsonio_check_error
  public :: jsonio_put_metadata
  public :: jsonio_put_attribute
  public :: jsonio_get_metadata
  public :: jsonio_get_attribute

  !
  ! check error
  !
  interface jsonio_check_error
     module procedure &
          & check_error_core, &
          & check_error_file
  end interface jsonio_check_error

  !
  ! put metadata
  !
  interface jsonio_put_metadata
     module procedure &
          & put_metadata
  end interface jsonio_put_metadata

  !
  ! put attribute; get metadata and data
  !
  interface jsonio_put_attribute
     module procedure &
          & put_scalar_i4, &
          & put_array_i4,  &
          & put_scalar_i8, &
          & put_array_i8,  &
          & put_scalar_r4, &
          & put_array_r4,  &
          & put_scalar_r8, &
          & put_array_r8
  end interface jsonio_put_attribute

  !
  ! get metadata
  !
  interface jsonio_get_metadata
     module procedure &
          & get_metadata
  end interface jsonio_get_metadata

  !
  ! get attribute; get metadata and data
  !
  interface jsonio_get_attribute
     module procedure &
          & get_scalar_i4, &
          & get_array_i4,  &
          & get_scalar_r4, &
          & get_array_r4,  &
          & get_scalar_r8, &
          & get_array_r8
  end interface jsonio_get_attribute

contains

  ! put metadata
  subroutine put_metadata(json, dst, name, datatype, offset, &
       & dsize, dshape, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: datatype
    integer(IK), intent(in)        :: offset
    integer(IK), intent(in)        :: dsize
    integer(IK), intent(in)        :: dshape(:)
    character(len=*), intent(in)   :: desc

    type(json_value), pointer :: obj
    integer(IK) :: ndim

    ndim = size(dshape)

    call json%create_object(obj, name)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'datatype', datatype)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'offset', offset)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'size', dsize)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'ndim', ndim)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'shape', dshape)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'description', trim(desc))
    call jsonio_check_error(json, 'put_metadata')

    ! add to dst
    call json%add(dst, obj)
    call jsonio_check_error(json, 'put_metadata')

    nullify(obj)

  end subroutine put_metadata

  ! put scalar of integer(4)
  subroutine put_scalar_i4(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(4), intent(in)         :: data
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    integer(IK) :: x, s, d(1)

    x = data
    s = 4
    d = (/1/)
    call jsonio_put_metadata(json, dst, name, 'i4', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_i4

  ! put array of integer(4)
  subroutine put_array_i4(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(4), intent(in)         :: data(:)
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    integer(IK) :: x(size(data)), s, d(rank(data))

    x = data
    s = size(data) * 4
    d = shape(data)
    call jsonio_put_metadata(json, dst, name, 'i4', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_i4

  ! put scalar of integer(8)
  subroutine put_scalar_i8(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(8), intent(in)         :: data
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    integer(IK) :: x, s, d(1)

    x = data
    s = 8
    d = (/1/)
    call jsonio_put_metadata(json, dst, name, 'i8', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_i8

  ! put array of integer(8)
  subroutine put_array_i8(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(8), intent(in)         :: data(:)
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    integer(IK) :: x(size(data)), s, d(rank(data))

    x = data
    s = size(data) * 8
    d = shape(data)
    call jsonio_put_metadata(json, dst, name, 'i8', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_i8

  ! put scalar of real(4)
  subroutine put_scalar_r4(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(4), intent(in)            :: data
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    real(RK) :: x
    integer(IK) :: s, d(1)

    x = data
    s = 4
    d = (/1/)
    call jsonio_put_metadata(json, dst, name, 'f4', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_r4

  ! put array of real(4)
  subroutine put_array_r4(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(4), intent(in)            :: data(:)
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    real(RK) :: x(size(data))
    integer(IK) :: s, d(rank(data))

    x = data
    s = size(data) * 4
    d = shape(data)
    call jsonio_put_metadata(json, dst, name, 'f4', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_r4

  ! put scalar of real(8)
  subroutine put_scalar_r8(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(8), intent(in)            :: data
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    real(RK) :: x
    integer(IK) :: s, d(1)

    x = data
    s = 8
    d = (/1/)
    call jsonio_put_metadata(json, dst, name, 'f8', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_r8

  ! put array of real(8)
  subroutine put_array_r8(json, dst, data, name, offset, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(8), intent(in)            :: data(:)
    character(len=*), intent(in)   :: name
    integer(IK), intent(in)        :: offset
    character(len=*), intent(in)   :: desc

    real(RK) :: x(size(data))
    integer(IK) :: s, d(rank(data))

    x = data
    s = size(data) * 8
    d = shape(data)
    call jsonio_put_metadata(json, dst, name, 'f8', offset, s, d, desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_r8

  !
  ! get metadata
  !
  subroutine get_metadata(json, src, name, offset, dsize, ndim, dshape)
    implicit none
    type(json_core), intent(inout)    :: json
    type(json_value), pointer         :: src
    character(len=*), intent(in)      :: name
    integer(IK), intent(out)          :: offset
    integer(IK), intent(out)          :: dsize
    integer(IK), intent(out)          :: ndim
    integer(IK), intent(out)          :: dshape(:)

    character(len=128) :: dataname
    integer(IK) :: i

    call json%get(src, name // '.offset', offset)
    call jsonio_check_error(json, 'get_metadata')

    call json%get(src, name // '.size', dsize)
    call jsonio_check_error(json, 'get_metadata')

    call json%get(src, name // '.ndim', ndim)
    call jsonio_check_error(json, 'get_metadata')

    do i = 1, ndim
       write(dataname,'(A,".shape(",i2,")")') name, i
       call json%get(src, trim(dataname), dshape(i))
       call jsonio_check_error(json, 'get_metadata')
    end do

  end subroutine get_metadata

  ! get scalar of integer(4)
  subroutine get_scalar_i4(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    integer(4), intent(out)        :: data

    integer(IK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_i4

  ! get array of integer(4)
  subroutine get_array_i4(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    integer(4), intent(out)        :: data(:)

    character(len=128) :: dataname
    integer(IK) :: i, x, dsize, ndim, dshape(1)

    call jsonio_get_metadata(json, src, name, offset, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_i4

  ! get scalar of integer(8)
  subroutine get_scalar_i8(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    integer(8), intent(out)        :: data

    integer(IK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_i8

  ! get array of integer(8)
  subroutine get_array_i8(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    integer(8), intent(out)        :: data(:)

    character(len=128) :: dataname
    integer(IK) :: i, x, dsize, ndim, dshape(1)

    call jsonio_get_metadata(json, src, name, offset, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_i8

  ! get scalar of real(4)
  subroutine get_scalar_r4(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    real(4), intent(out)           :: data

    real(RK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_r4

  ! get array of real(4)
  subroutine get_array_r4(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    real(4), intent(out)           :: data(:)

    character(len=128) :: dataname
    integer(IK) :: i, ndim, dsize, dshape(1)
    real(RK) :: x

    call jsonio_get_metadata(json, src, name, offset, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_r4

  ! get scalar of real(8)
  subroutine get_scalar_r8(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    real(8), intent(out)           :: data

    real(RK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_r8

  ! get array of real(8)
  subroutine get_array_r8(json, src, name, offset, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(IK), intent(out)       :: offset
    real(8), intent(out)           :: data(:)

    character(len=128) :: dataname
    integer(IK) :: i, ndim, dsize, dshape(1)
    real(RK) :: x

    call jsonio_get_metadata(json, src, name, offset, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_r8

  ! error check code
  subroutine check_error_core(json, message)
    implicit none
    type(json_core), intent(inout) :: json
    character(len=*), intent(in)   :: message

    if( json%failed() ) then
       write(0, *) 'Some error detected during json operation: ', message
       call json%print_error_message(0)
    end if

  end subroutine check_error_core

  ! error check code
  subroutine check_error_file(json, message)
    implicit none
    type(json_file), intent(inout) :: json
    character(len=*), intent(in)   :: message

    if( json%failed() ) then
       write(0, *) 'Some error detected during json operation: ', message
       call json%print_error_message(0)
    end if

  end subroutine check_error_file

end module jsonio
