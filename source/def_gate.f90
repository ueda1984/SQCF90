module def_gate
  implicit none
  complex(8), parameter :: zu=(1.0d0,0.0d0)
  complex(8), parameter :: zi=(0.0d0,1.0d0)
  complex(8), parameter :: zz=(0.0d0,0.0d0)
  complex(8) :: mm(4,4), tcmm(4,4)
  complex(8) :: px(2,2), py(2,2), pz(2,2), p0(2,2)
  complex(8) :: pxx(4,4), pyy(4,4), pzz(4,4)
  real(8) :: pi
  integer :: lo(2)=(/2,2/)
  integer :: lo2(2)=(/4,4/)
  integer :: lo3(2)=(/8,8/)
contains
  !
  subroutine set_mat
    pi = acos(-1.0d0)
    !
    px(1,1:2)=(/ zz,  zu /)
    px(2,1:2)=(/ zu,  zz /)
    !
    py(1,1:2)=(/ zz, -zi /)
    py(2,1:2)=(/ zi,  zz /)
    !
    pz(1,1:2)=(/ zu,  zz /)
    pz(2,1:2)=(/ zz, -zu /)
    !
    p0(1,1:2)=(/zu, zz/)
    p0(2,1:2)=(/zz, zu/)
    !
    mm(1,:)=(/zu, zi, zz, zz/)
    mm(2,:)=(/zz, zz, zi, zu/)
    mm(3,:)=(/zz, zz, zi,-zu/)
    mm(4,:)=(/zu,-zi, zz, zz/)
    mm = mm/sqrt(2.0d0)
    !
    tcmm = transpose(conjg(mm))
    !
    pxx = f_tens_prod(px,lo,px,lo)
    pyy = f_tens_prod(py,lo,py,lo)
    pzz = f_tens_prod(pz,lo,pz,lo)
    !
  end subroutine set_mat
  !
  function f_rz(th) result(mat)
    real(8), intent(in) :: th
    complex(8) :: mat(2,2)
    mat = p0*cos(th/2.0d0) + zi * pz *sin(th/2.0d0)
  end function f_rz
  !
  function f_ry(th) result(mat)
    real(8), intent(in) :: th
    complex(8) :: mat(2,2)
    mat = p0*cos(th/2.0d0) + zi * py *sin(th/2.0d0)
  end function f_ry
  !
  function f_su2(th) result(mat)
    real(8) :: th(3)
    complex(8) :: mat(2,2)
    mat = matmul(matmul(f_rz(th(3)),f_ry(th(2))),f_rz(th(1)))
  end function f_su2
  !
  function f_su2_sub(th) result(mat)
    real(8) :: th(3)
    complex(8) :: mat(2,2)
    mat(1,:) = (/exp(zi/2*(th(1)+th(3)))*cos(th(2)/2), exp(zi/2*(-th(1)+th(3)))*sin(th(2)/2)/)
    mat(2,:) = (/-exp(zi/2*(th(1)-th(3)))*sin(th(2)/2), exp(-zi/2*(th(1)+th(3)))*cos(th(2)/2)/)
  end function f_su2_sub
  !
  function f_tens_prod(mat1,l1,mat2,l2) result(mat)
    integer :: l1(2), l2(2)
    complex(8) :: mat1(l1(1),l1(2)), mat2(l2(1),l2(1))
    complex(8) :: mat(l1(1)*l2(1),l1(2)*l2(2))
    integer :: i1,j1,j2,j
    !$omp parallel do private(j1,j2,i1)
    do j = 1, l1(2)*l2(2)
      j1 = (j-1) / l2(2) + 1
      j2 = j - l2(2)*(j1-1)
      do i1 = 1, l1(1)
        mat(l1(2)*(i1-1)+1:l1(2)*i1,j) = mat1(i1,j1)*mat2(:,j2)
      end do
    end do
  end function f_tens_prod  
  !
  function f_so4(th) result(mat)
    real(8) :: th(6)
    complex(8) :: mat(4,4),su2(2,2),su2_2(2,2)
    su2 = f_su2(th(1:3))
    su2_2 = f_su2(th(4:6))
    mat = matmul(tcmm,matmul(f_tens_prod(su2,lo,su2_2,lo),mm))
  end function f_so4  
  !
end module def_gate
