module hilbert
!
!    This module provides basic operations with vectors: inner products and common norms.
!
  use constants
  implicit none
  include "mkl_blas.fi"
  
contains

  function norm(vector,L)
!
!     This function returns the l_1 norm (L=1), the l_2 norm (L=2), 
!     or the l_infinity norm (L=3) of a vector.
!
    implicit none
    integer :: L, n
    real(kind=dble_prec) :: norm
    real(kind=dble_prec) :: vector(:)

    n = size(vector)

    select case(L)
    case(1)
       norm = sum(abs(vector))
    case(2)
       norm = dnrm2(n,vector,1)
!       dsqrt(ddot(n,vector,1,vector,1))
    case(3)
       norm = maxval(abs(vector))
    end select
  end function norm

  function inner_product(a,b)
    implicit none
    integer :: n
    real(kind=dble_prec) :: inner_product
    real(kind=dble_prec) :: a(:), b(:)

    n = size(a)
    inner_product = ddot(n,a,1,b,1)

  end function inner_product

end module hilbert
