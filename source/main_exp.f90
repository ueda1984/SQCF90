program main
  implicit none
  call sub_main('./input/HC_1.dat','ene_HC_1_03-06.dat',3,6)
  call sub_main('./input/HC_2.dat','ene_HC_2_06-07.dat',6,7)
  call sub_main('./input/HC_3.dat','ene_HC_3_07-08.dat',7,8)
  call sub_main('./input/K_1.dat','ene_K_1_07-08.dat',7,8)
  call sub_main('./input/K_2.dat','ene_K_2_10-12.dat',10,12)
  call sub_main('./input/K_3.dat','ene_K_3_05-06.dat',5,6)
  call sub_main('./input/K_4.dat','ene_K_4_02-07.dat',2,7)
  call sub_main('./input/K_4.dat','ene_K_4_06-07.dat',6,7)
  call sub_main('./input/K_4.dat','ene_K_4_02-06.dat',2,6)
  call sub_main('./input/S_1.dat','ene_S_1_09-10.dat',9,10)
  call sub_main('./input/S_1.dat','ene_S_1_15-16.dat',15,16)
  call sub_main('./input/S_2.dat','ene_S_2_10-14.dat',10,14)
  call sub_main('./input/S_3.dat','ene_S_3_11-15.dat',11,15)
  call sub_main('./input/S_4.dat','ene_S_4_06-07.dat',6,7)
  call sub_main('./input/S_4.dat','ene_S_4_10-11.dat',10,11)
  call sub_main('./input/S_5.dat','ene_S_5_07-11.dat',7,11)
  call sub_main('./input/S_6.dat','ene_S_6_06-10.dat',6,10)
contains 
  subroutine sub_main(ifile,ofile,p1,p2)
    use omp_lib
    use prod_gate
    use set_gate_list
    implicit none
    character(*), intent(in) :: ifile, ofile
    integer, intent(in) :: p1, p2
    integer :: nol, nog, tq, noth
    character(2) :: cnol
    real(8) :: e(3), st, en, mtime
    real(8), allocatable :: ths(:,:)
    integer, allocatable :: list_p(:,:)
    complex(8), allocatable :: v1(:), v2(:)
    complex(8) :: zdotc
    integer :: lv, i
    call set_mat
    mtime = 1.0d9
    call read_2qgate_list(list_p,tq,nol,nog,ifile)
    write(cnol,'(i2.2)') nol
    call read_so4_theta_list(ths,noth,nol,'./input/thetas_'//cnol//'.dat')
    lv = 2**tq
    allocate(v1(lv), v2(lv))
    open(20,file=ofile)
    write(20,*) "nol, nog, noth, tq, lv", nol, nog, noth, tq, lv
    do i = 1, noth
       st = omp_get_wtime()
       call phi_1_i(v1,lv,list_p,nog,ths(:,i),nol)
       call gate2_vec_o(v2,v1,lv,pxx,p1,p2)
       e(1)=dble(zdotc(lv,v1,1,v2,1))
       call gate2_vec_o(v2,v1,lv,pyy,p1,p2)
       e(2)=dble(zdotc(lv,v1,1,v2,1))
       call gate2_vec_o(v2,v1,lv,pzz,p1,p2)
       e(3)=dble(zdotc(lv,v1,1,v2,1))
       en = omp_get_wtime()
       if(en-st < mtime) mtime = en-st
       write(20,'(3es23.15)') e
    end do
99  write(20,*) "Min. of elapsed time in second is:", mtime
    close(10)
    close(20)
    deallocate(ths,list_p,v1,v2)
  end subroutine sub_main
  !
end program main
