module set_gate_list
  implicit none
contains
  !
  subroutine read_2qgate_list(list_p,tq,nol,nog,ifile)
    integer, allocatable, intent(out) :: list_p(:,:)
    integer, intent(out) :: tq, nol, nog
    character(*), intent(in) :: ifile
    integer :: p1,p2,ln,i
    open(10,file=ifile,action='read')
    nog = 0
    do 
       read(10,*,end=9) p1, p2, ln
       nog = nog + 1
    end do
9   rewind(10)
    allocate(list_p(3,nog))
    do i = 1, nog
       read(10,*) p1, p2, ln
       list_p(1:3,i) = (/p1, p2, ln/)
    end do 
    close(10)
    nol = maxval(list_p(3,:))
    tq = maxval(list_p(1:2,:))
  end subroutine read_2qgate_list
  !
  subroutine read_so4_theta_list(ths,noth,nol,ifile)
    integer, intent(in) :: nol
    character(*), intent(in) :: ifile
    integer, intent(out) :: noth
    real(8), allocatable, intent(out) :: ths(:,:)
    real(8) :: th(6*nol)
    integer :: i
    open(10,file=ifile,action='read')
    noth = 0
    do 
       read(10,*,end=9) th
       noth = noth + 1
    end do
9   rewind(10)
    allocate(ths(6*nol,noth))
    do i = 1, noth
       read(10,*) th
       ths(:,i) = th
    end do 
    close(10)
  end subroutine read_so4_theta_list
  !
end module set_gate_list
