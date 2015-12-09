module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param


contains


  subroutine fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,nproc,nrank, &
                         c,q,r,delt,delx,it,it0,dir,lflag)

    logical, intent(in) :: lflag
    integer, intent(in) :: np, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, nsp
    integer, intent(in) :: np2(nys:nye,nsp)
    integer, intent(in) :: nproc, nrank
    integer, intent(in) :: it, it0
    real(8), intent(in) :: up(5,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(in) :: c, q(nsp), r(nsp), delt, delx
    character(len=*), intent(in) :: dir
    integer :: it2
    character(len=256) :: filename

    it2=it+it0

    !filename
    if(lflag)then
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),9999999,'_rank=',nrank,'.dat'
    else
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it2,'_rank=',nrank,'.dat'
    endif
    open(200+nrank,file=filename,form='unformatted')

    !time & parameters
    write(200+nrank)it2,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,nproc,-1,delt,delx,c
    write(200+nrank)np2
    write(200+nrank)q
    write(200+nrank)r

    !field data
    write(200+nrank)uf

    !particle data
    write(200+nrank)up

    close(200+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,nxs,nxe,c,q,r,delt,delx,it0,          &
                        np,nxgs,nxge,nygs,nyge,nys,nye,nsp,nproc,nrank, &
                        dir,file)
    integer, intent(in)  :: np, nxgs, nxge, nygs, nyge, nys, nye, nsp, nproc, nrank
    character(len=*), intent(in) :: dir, file
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it0
    real(8), intent(out) :: up(5,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer :: inp, inxgs, inxge, inygs, inyge, inys, inye, insp, inproc, ibc

    !filename
    open(201+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(201+nrank)it0,inp,inxgs,inxge,inygs,inyge,nxs,nxe,inys,inye,insp,inproc,ibc,delt,delx,c
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge) &
        .or. (inys /= nys) .or. (inye /= nye) .or. (inp /= np) .or. (insp /= nsp) &
        .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    read(201+nrank)np2
    read(201+nrank)q
    read(201+nrank)r

    !field data
    read(201+nrank)uf

    !particle data
    read(201+nrank)up

    close(201+nrank)

  end subroutine fio__input


  subroutine fio__param(np,n0,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye, &
                        c,q,r,temp,rtemp,fpe,fge,            &
                        ldb,delt,delx,dir,file,                 &
                        nroot,nrank)

    integer, intent(in)          :: np, n0, nsp
    integer, intent(in)          :: nxgs, nxge, nygs, nyge, nys, nye
    integer, intent(in)          :: nroot, nrank
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: c, q(nsp), r(nsp), temp, rtemp, fpe, fge, ldb, delt, delx
    character(len=*), intent(in) :: dir, file
    integer :: isp
    real(8) :: pi, vti, vte, va

    pi = 4.0D0*datan(1.0D0)

    vti = sqrt(2.*temp/r(1))
    vte = sqrt(2.*temp*rtemp/r(2))
    va  = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    if(nrank == nroot)then

       !filename
       open(9,file=trim(dir)//trim(file),status='unknown')

       write(9,610) nxge-nxgs+1,' x ',nyge-nygs+1, ldb
       write(9,620) (np2(nys,isp),isp=1,nsp),np
       write(9,630) delx,delt,c
       write(9,640) (r(isp),isp=1,nsp)
       write(9,650) (q(isp),isp=1,nsp)
       write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
       write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
       write(9,*)
610    format(' grid size, debye lngth ============> ',i6,a,i6,f8.4)
620    format(' particle number in cell============> ',i8,i8,'/',i8)
630    format(' dx, dt, c =========================> ',f8.4,3x,f8.4,3x,f8.4)
640    format(' Mi, Me  ===========================> ',2(1p,e10.2,1x))
650    format(' Qi, Qe  ===========================> ',2(1p,e10.2,1x))
660    format(' Fpe, Fge, Fpi Fgi =================> ',4(1p,e10.2,1x))
670    format(' Va, Vi, Ve, beta, Te/Ti, rgi     ==> ',6(1p,e10.2,1x))
       close(9)

    endif

  end subroutine fio__param


end module fio
