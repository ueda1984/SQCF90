module prod_gate
  use def_gate
  implicit none
contains
  !
  subroutine gate2_vec_o(vo,v,lv,g,p1,p2)
    use omp_lib
    integer, intent(in) :: lv, p1, p2
    complex(8), intent(in) :: v(0:lv-1)
    complex(8), intent(out) :: vo(0:lv-1)
    complex(8), intent(in) :: g(0:1,0:1,0:1,0:1)
    complex(8) :: gt(0:1,0:1,0:1,0:1)
    integer :: p ,i, j, id, jd, p1m, p2m, u0, up1, up2, up12
    if(p1 < p2)then
       forall(i=0:1,j=0:1,id=0:1,jd=0:1) gt(id,jd,i,j) = g(i,j,id,jd)
    else
       forall(i=0:1,j=0:1,id=0:1,jd=0:1) gt(jd,id,j,i) = g(i,j,id,jd)
    end if
    p1m = p1 - 1 
    p2m = p2 - 1 
    u0  = 1
    up1 = 2**p1m
    up2 = 2**p2m
    up12= up1 + up2
    !$omp parallel do private(i,j,id,jd)
    do p = 0, lv-1
       i = ibits(p,p1m,1)
       j = ibits(p,p2m,1)
       id= ieor(i,u0)
       jd= ieor(j,u0)
       vo(p) = gt(i,j,i,j) * v(p) + &
            gt(id,j,i,j) * v(ieor(p,up1)) + &
            gt(i,jd,i,j) * v(ieor(p,up2)) + &
            gt(id,jd,i,j) * v(ieor(p,up12))
    end do
  end subroutine gate2_vec_o
  !
  subroutine gate2_vec_i(v,lv,g,p1,p2)
    integer, intent(in) :: lv, p1, p2
    complex(8), intent(inout) :: v(0:lv-1)
    complex(8):: vo(0:lv-1)
    complex(8), intent(in) :: g(2,2,2,2)
    call gate2_vec_o(vo,v,lv,g,p1,p2)
    call zcopy(lv,vo,1,v,1)
  end subroutine gate2_vec_i
  !
  subroutine phi_1_i(v1,lv,list_p,nog,ths,nol)
    integer, intent(in) :: lv, nog, nol
    integer, intent(in) :: list_p(3,nog)
    real(8), intent(in) :: ths(6,nol)
    complex(8), intent(out) :: v1(lv)
    complex(8) :: v2(lv)
    complex(8) :: g(4,4)
    integer :: i, ip
    !
    v1(1) = zu
    !$omp parallel do
    do i = 2, lv
       v1(i) = zz
    end do
    !
    do i = 1, nog-1, 2
       g = f_so4(ths(:,list_p(3,i)))
       call gate2_vec_o(v2,v1,lv,g,list_p(1,i),list_p(2,i))
       ip = i + 1
       g = f_so4(ths(:,list_p(3,ip)))
       call gate2_vec_o(v1,v2,lv,g,list_p(1,ip),list_p(2,ip))
    end do
    if(ip==nog)then
       return
    else
       g = f_so4(ths(:,list_p(3,nog)))
       call gate2_vec_i(v1,lv,g,list_p(1,nog),list_p(2,nog))
    end if
    !
  end subroutine phi_1_i
  !
end module prod_gate
