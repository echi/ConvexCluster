subroutine test_prox(vector_in,n,vector_out,tau,type)
  use constants
  use proximal
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau
  call prox(vector_in,n,vector_out,tau,type)
end subroutine test_prox
