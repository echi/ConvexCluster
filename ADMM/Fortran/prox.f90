module proximal
  use constants
  use hilbert
  
contains
  
  subroutine proj_L2(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    lv = norm(x,2)
    px = x
    if (lv > tau) then
       px = (tau/lv)*x
    end if
  end subroutine proj_L2

  subroutine proj_Linf(x,n,px,tau)
    implicit none
    integer :: i,n
    real(kind=dble_prec) :: ax(n), px(n), tau, x(n)

    ax = abs(x)
    do i=1,n
       px(i) = min(ax(i),tau)
    end do
    px = sign(px,x)
  end subroutine proj_Linf
  
  subroutine prox_L1(x,n,px,tau)
    implicit none
    integer :: i,n
    real(kind=dble_prec) :: px(n), tau, x(n)
      
    px = abs(x) - tau
    do i=1,n
       px(i) = max(px(i),zero)
    end do
    px = sign(px,x)
    
  end subroutine prox_L1
  
  subroutine prox_L2(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    px = zero
 !   lv = dsqrt(dot_product(x,x))                                                                                                                                                                     
    lv = norm(x,2)
    
    if (lv.eq.zero) then
       px = x
    else
       px = max(zero,one-(tau/lv))*x
    end if
  end subroutine prox_L2

  subroutine prox_Linf(x,n,px,tau)
!
! This function performs the proximal mapping of tau * L-infinity norm.
! It is computed via Moreau decomposition and Projection onto
! an L1-ball of radius tau.
!
!   px = x - project_L1(x,tau)
!
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    px = (one/tau)*x
    call proj_L1(x,n,px,tau)
    px = x - px
  end subroutine prox_Linf

  subroutine project_to_simplex(x,n,z)
    use SORT
    implicit none
    integer :: j, n, rho
    real(kind=dble_prec) :: x(n), z
    real(kind=dble_prec) :: cumsum, mu(n), theta
    mu = x
    call quick_sort(mu)
    cumsum = mu(1)
    do j = 2,n
       cumsum = cumsum + mu(j)
       if (dble(j)*mu(j) - cumsum + z .LE. zero) then
          rho = j-1
          exit
       end if
       rho = j
    end do
    theta = (sum(mu(1:rho)) - z) / dble(rho)
    do j = 1,n
       x(j) = max(x(j) - theta, zero)
    end do
  end subroutine project_to_simplex

  subroutine proj_L1(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: px(n), tau, x(n)
    px = abs(x)
    call project_to_simplex(px,n,tau)
    px = sign(px,x)
  end subroutine proj_L1

  subroutine prox(vector_in,n,vector_out,tau,type)
!
! This function is a wrapper for performing various proximal and projection mappings.
!
! Menu of types:
!  type = 1 : proximal : L1
!  type = 2 : proximal : L2
!  type = 3 : proximal : L-infinity
!  type = 4 : project  : L-infinity
!  type = 5 : project  : L2
!  type = 6 : project  : L1
!  type = 7 : project  : L12 *** Not implemented
!  type = 8 : proximal : L12 *** Not implememted
!
! Notes:
!  March 13, 2013: Added proximal mapping of L-infinity norm
!  March 13, 2013: Swapped order 1->4, 2->5, 3->6, 4->1, 5->2, 6->3
!
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau

  select case(type)
    case(1)
       call prox_L1(vector_in,n,vector_out,tau)
    case(2)
       call prox_L2(vector_in,n,vector_out,tau)
    case(3)
       call prox_Linf(vector_in,n,vector_out,tau)
    case(4)
       call proj_Linf(vector_in,n,vector_out,tau)
    case(5)
       call proj_L2(vector_in,n,vector_out,tau)
    case(6)
       call proj_L1(vector_in,n,vector_out,tau)
    end select

  end subroutine prox

end module proximal
