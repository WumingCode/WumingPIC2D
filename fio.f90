module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__energy


contains


  subroutine fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,bc,nproc,nrank, &
                         c,q,r,delt,delx,it,it0,dir)

    integer, intent(in) :: np, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, nsp, bc
    integer, intent(in) :: np2(nys:nye,nsp)
    integer, intent(in) :: nproc, nrank
    integer, intent(in) :: it, it0
    real(8), intent(in) :: up(5,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxs-1:nxe+1,nys-1:nye+1)
    real(8), intent(in) :: c, q(nsp), r(nsp), delt, delx
    character(len=*), intent(in) :: dir
    integer :: it2
    character(len=256) :: filename

    it2=it+it0

    !filename
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it2,'_rank=',nrank,'.dat'
    open(100+nrank,file=filename,form='unformatted')

    !time & parameters
    write(100+nrank)it2,np,nxgs,nxge,nygs,nyge,nsp,nproc,bc,delt,delx,c
    write(100+nrank)np2
    write(100+nrank)q
    write(100+nrank)r

    !field data
    write(100+nrank)uf

    !particle data
    write(100+nrank)up
       
    close(100+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,c,q,r,delt,delx,it0,                                     &
                        np,np2,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,bc,nproc,nrank, &
                        dir,file)
    integer, intent(in)  :: np, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, nsp, bc, nproc, nrank
    character(len=*), intent(in) :: dir, file
    integer, intent(out) :: np2(nys:nye,nsp), it0
    real(8), intent(out) :: up(5,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxs-1:nxe+1,nys-1:nye+1)
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer :: inp, inxgs, inxge, inygs, inyge, insp, ibc, inproc

    !filename
    open(101+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(101+nrank)it0,inp,inxgs,inxge,inygs,inyge,insp,inproc,ibc,delt,delx,c
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge) &
       .or. (inp /= np) .or. (insp /= nsp) .or. (ibc /= bc) .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif
    read(101+nrank)np2
    read(101+nrank)q
    read(101+nrank)r

    !field data
    read(101+nrank)uf

    !particle data
    read(101+nrank)up

    close(101+nrank)

  end subroutine fio__input


  subroutine fio__param(np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye, &
                        c,q,r,n0,temp,rtemp,fpe,fge,            &
                        ldb,delt,delx,dir,file,                 &
                        nroot,nrank)

    integer, intent(in)          :: np, nsp 
    integer, intent(in)          :: nxgs, nxge, nygs, nyge, nys, nye
    integer, intent(in)          :: nroot, nrank
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: c, q(nsp), r(nsp), n0, temp, rtemp, fpe, fge, ldb, delt, delx
    character(len=*), intent(in) :: dir, file
    integer :: isp
    real(8) :: pi, vti, vte, va

    pi = 4.0*atan(1.0)

    vti = sqrt(2.*temp/r(1))
    vte = sqrt(2.*temp*rtemp/r(2))
    va  = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    if(nrank == nroot)then

       !filename
       open(9,file=trim(dir)//trim(file),status='unknown')

       write(9,610) nxge-nxgs+1,'x',nyge-nygs+1, ldb
       write(9,620) (np2(nys,isp),isp=1,nsp),np
       write(9,630) delx,delt,c
       write(9,640) (r(isp),isp=1,nsp)
       write(9,650) (q(isp),isp=1,nsp)
       write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
       write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
       write(9,*)
610    format(' grid size, debye lngth ============> ',i4,a,i4,f8.4)
620    format(' particle number in cell============> ',i8,i8,'/',i8)
630    format(' dx, dt, c =========================> ',f8.4,3x,f8.4,3x,f8.4)
640    format(' Mi, Me  ===========================> ',2(1p,e10.2,1x))
650    format(' Qi, Qe  ===========================> ',2(1p,e10.2,1x))
660    format(' Fpe, Fge, Fpi Fgi =================> ',4(1p,e10.2,1x))
670    format(' Va, Vi, Ve, beta, Te/Ti, rgi     ==> ',6(1p,e10.2,1x))
       close(9)

    endif
    
  end subroutine fio__param


  subroutine fio__energy(up,uf,np,nsp,np2,nxs,nxe,nys,nye,bc, &
                         c,r,delt,it,it0,dir,file,            &
                         nroot,nrank,mnpr,opsum,ncomw,nerr)

    integer, intent(in)          :: nxs, nxe, nys, nye, bc
    integer, intent(in)          :: nroot, nrank, mnpr, opsum, ncomw
    integer, intent(in)          :: it, it0, np, nsp, np2(nys:nye,nsp)
    integer, intent(inout)       :: nerr
    real(8), intent(in)          :: c, r(nsp), delt
    real(8), intent(in)          :: up(5,np,nys:nye,nsp)
    real(8), intent(in)          :: uf(6,nxs-1:nxe+1,nys-1:nye+1)
    character(len=*), intent(in) :: dir, file
    integer :: i, j, ii, isp
    integer, save :: iflag
    real(8) :: pi
    real(8) :: vene(nsp), vene_g(nsp)
    real(8) :: efield, bfield, gam, total, u2
    real(8) :: efield_g, bfield_g

    pi = 4.0*atan(1.0)

    !filename
    if(iflag /= 1)then
       if(nrank == nroot)then
          if(it == 0) open(12,file=trim(dir)//trim(file),status='unknown')
       endif
       iflag = 1
    endif

    !energy
    vene(1:nsp) = 0.0
    efield = 0.0
    bfield = 0.0

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             u2 =  up(3,ii,j,isp)*up(3,ii,j,isp) &
                  +up(4,ii,j,isp)*up(4,ii,j,isp) &
                  +up(5,ii,j,isp)*up(5,ii,j,isp)
             gam = dsqrt(1.0+u2/(c*c))
             vene(isp) = vene(isp)+r(isp)*u2/(gam+1.)
          enddo
       enddo
    enddo

    do isp=1,nsp
       call MPI_REDUCE(vene(isp),vene_g(isp),1,mnpr,opsum,nroot,ncomw,nerr)
    enddo

    do j=nys,nye
    do i=nxs,nxe+bc
       bfield = bfield+uf(1,i,j)*uf(1,i,j)+uf(2,i,j)*uf(2,i,j)+uf(3,i,j)*uf(3,i,j)
       efield = efield+uf(4,i,j)*uf(4,i,j)+uf(5,i,j)*uf(5,i,j)+uf(6,i,j)*uf(6,i,j)
    enddo
    enddo
    if(bc == -1)then
       i=nxe
       do j=nys,nye
          bfield = bfield+uf(2,i,j)*uf(2,i,j)+uf(3,i,j)*uf(3,i,j)
          efield = efield+uf(4,i,j)*uf(4,i,j)
       enddo
    endif

    efield = efield/(8.0*pi)
    bfield = bfield/(8.0*pi)
    call MPI_REDUCE(efield,efield_g,1,mnpr,opsum,nroot,ncomw,nerr)
    call MPI_REDUCE(bfield,bfield_g,1,mnpr,opsum,nroot,ncomw,nerr)

    if(nrank == nroot)then
       total=vene_g(1)+vene_g(2)+efield_g+bfield_g
       
       write(12,610) (it+it0)*delt,vene_g(1),vene_g(2),efield_g,bfield_g,total
610    format(f10.2,5(e12.4))
    endif

  end subroutine fio__energy


end module fio
