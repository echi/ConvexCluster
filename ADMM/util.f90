module constants
  implicit none
  integer, parameter :: dble_prec = kind(0.0d0)
  real(kind=dble_prec), parameter :: zero  = 0.0_dble_prec
  real(kind=dble_prec), parameter :: one   = 1.0_dble_prec
  real(kind=dble_prec), parameter :: two   = 2.0_dble_prec
  real(kind=dble_prec), parameter :: three = 3.0_dble_prec
  real(kind=dble_prec), parameter :: four  = 4.0_dble_prec
  real(kind=dble_prec), parameter :: five  = 5.0_dble_prec
  real(kind=dble_prec), parameter :: six   = 6.0_dble_prec
  real(kind=dble_prec), parameter :: seven = 7.0_dble_prec
  real(kind=dble_prec), parameter :: eight = 8.0_dble_prec
  real(kind=dble_prec), parameter :: nine  = 9.0_dble_prec
  real(kind=dble_prec), parameter :: half  = 0.5_dble_prec
end module constants

subroutine vec2tri(ix,k,n,p)
  use constants
  implicit none
  integer :: n, p
  integer :: ix(n,2), k(n)
  ix(:,1) = ceiling(half*(two*p-one - sqrt((two*p-one)**2 - eight*k)))
  ix(:,2) = k - p*(ix(:,1)-1) + ix(:,1)*(ix(:,1)-1)/2 + ix(:,1)
end subroutine vec2tri

subroutine tri2vecA(i,j,p,k,n)
  implicit none
  integer :: j,p,n
  integer :: i(n), k(n)
  k = p*(i-1) - i*(i-1)/2 + j - i
end subroutine tri2vecA

subroutine tri2vecB(i,j,p,k,n)
  implicit none
  integer :: i,p,n
  integer :: j(n), k(n)
  k = p*(i-1) - i*(i-1)/2 + j - i
end subroutine tri2vecB

subroutine prox_L2(x,n,px,tau)
  use constants
  implicit none
  include "mkl_blas.fi"
  integer :: n
  real(kind=dble_prec) :: lv, px(n), tau, x(n), ddot_result
  px = zero
  lv = sqrt(dot_product(x,x))                                                                                                        
  if (lv.eq.zero) then
     px = x
  else
     px = max(zero,one-tau/lv)*x
  end if
end subroutine prox_L2

subroutine loss_primal(X,U,gamma,ix,p,q,nK,w,output)
  use constants
  implicit none
  integer :: k,p,q,nK
  integer :: ix(nK,2)
  real(kind=dble_prec) :: DU(q), D(q,p), X(q,p), U(q,p), gamma, w(nK), output, penalty
  penalty = zero
  do k=1,nK
     DU = U(:,ix(k,1)) - U(:,ix(k,2))
     penalty = penalty + w(k)*sqrt(dot_product(DU,DU))
  end do
  output = zero
  D = X - U
  do k=1,p
     output = output + dot_product(D(:,k),D(:,k))
  end do
  output = 0.5*output + gamma*penalty
end subroutine loss_primal

subroutine loss_dual(X,Lambda,p,q,nK,output)
  use constants
  implicit none
  integer :: ii,jj,k,p,q,nK
  integer :: ix(p), seq(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK)
  real(kind=dble_prec) :: L1(q), L2(q), diff(q)
  real(kind=dble_prec) :: first_term, second_term, output
  seq = (/(ii, ii=1,p, 1)/)

! Loop unrolling: ii=1
  ii = 1
  call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
  L2 = sum(Lambda(:,ix(1:(p-ii))),2)
  first_term = dot_product(L2,L2) 
! Loop over ii=2:p-1
  do ii=2,p-1
     call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
     L1 = sum(Lambda(:,ix(1:(ii-1))),2)
     call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
     L2 = sum(Lambda(:,ix(1:(p-ii))),2)
     diff = L1 - L2
     first_term = first_term + dot_product(diff,diff)
  end do
! Loop unrolling: ii = p
  ii = p
  call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
  L1 = sum(Lambda(:,ix(1:(ii-1))),2)
  first_term = first_term + dot_product(L1,L1)
  k = 1
  second_term = zero
  do ii=1,p-1
     do jj=ii+1,p
        second_term = second_term + dot_product(X(:,ii)-X(:,jj),Lambda(:,k))
        k = k + 1
     end do
  end do
  output = -half*first_term - second_term
end subroutine loss_dual

subroutine update_U(X,Lambda,U,p,q,nK)
  use constants
  implicit none
  integer :: ii,p,q,nK
  integer :: ix(p), seq(p)
  real(kind=dble_prec) :: X(q,p), U(q,p), Lambda(q,nK)
  real(kind=dble_prec) :: L1(q), L2(q)
  seq = (/(ii, ii=1,p, 1)/)

! Loop unrolling: ii=1
  ii = 1
  call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
  L2 = sum(Lambda(:,ix(1:(p-ii))),2)
  U(:,ii) = X(:,ii) + L2
! Loop over ii=2:p-1
  do ii=2,p-1
     call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
     L1 = sum(Lambda(:,ix(1:(ii-1))),2)
     call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
     L2 = sum(Lambda(:,ix(1:(p-ii))),2)
     U(:,ii) = X(:,ii) - L1 + L2
  end do
! Loop unrolling: ii = p
  ii = p
  call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
  L1 = sum(Lambda(:,ix(1:(ii-1))),2)
  U(:,ii) = X(:,ii) - L1
  
end subroutine update_U

subroutine update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK)
  use constants
  implicit none
  integer :: kk,nK,p,q
  integer :: i,j,ix(nK,2)
  real(kind=dble_prec) :: gamma, nu
  real(kind=dble_prec) :: U(q,p), Lambda(q,nK), V(q,nK), w(nK)
  real(kind=dble_prec) :: z(q)
  
  do kk=1,nK
     i = ix(kk,1)
     j = ix(kk,2)
     z = U(:,i) - U(:,j) - (one/nu)*Lambda(:,kk)
     call prox_L2(z,q,V(:,kk),w(kk)*gamma/nu)
  end do
  
end subroutine update_V

subroutine update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
  use constants
  implicit none
  integer :: nK, p, q
  integer :: ix(nK,2)
  real(kind=dble_prec) :: nu, Lambda(q,nK), U(q,p), V(q,nK)
  real(kind=dble_prec) :: Z(q,nK)
  Z = U(:,ix(:,1)) - U(:,ix(:,2)) - V
  Lambda = Lambda - nu*Z
end subroutine update_Lambda
