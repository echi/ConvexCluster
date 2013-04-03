subroutine get_xbar(X,xbar,p,q)
  use constants
  implicit none
  integer :: p,q
  real(kind=dble_prec) :: X(q,p), xbar(q)
  xbar = (one/p)*sum(X,2)
end subroutine get_xbar

!subroutine residual_primal(U,V,ix,p,q,nK,residual)
!
! This subroutine computes the L-infinity norm of the primal residual.
!
!  use constants
!  implicit none
!  integer :: ix(nK,2),p,q,nK
!  real(kind=dble_prec) :: U(q,p), V(q,nK), residual
!  
!  residual = maxval(abs(U(:,ix(:,1)) - U(:,ix(:,2)) - V))
!
!end subroutine residual_primal

subroutine residual_primal(U,V,ix,p,q,nK,residual)
!
! This subroutine computes the L2 norm of the primal residual.
!
  use constants
  implicit none
  integer :: ix(nK,2),p,q,nK
  real(kind=dble_prec) :: U(q,p), V(q,nK), residual
 
  residual = sqrt(sum( (U(:,ix(:,1)) - U(:,ix(:,2)) - V)**2))

end subroutine residual_primal

!subroutine residual_dual(V,V_old,M1,M2,s1,s2,mix1,mix2,p,q,nK,nu,residual)
!
! This subroutine computes the L-infinity norm of the dual residual.
!
!  use constants
!  implicit none
!  integer :: ii,mix1,mix2,nK,p,q
!  integer :: s1(p), s2(p), M1(mix1,p), M2(mix2,p)
!  real(kind=dble_prec) :: nu, residual, V(q,nK), V_old(q,nK), dual_residual(q)
!  real(kind=dble_prec) :: L1(q), L2(q)
!
!  residual = zero
!  do ii=1,p
!     dual_residual = zero
!     if (s1(ii) > 0) then
!        dual_residual = sum(V(:,M1(1:s1(ii),ii))-V_old(:,M1(1:s1(ii),ii)),2)
!     end if
!     if (s2(ii) > 0) then
!        dual_residual = dual_residual - sum(V(:,M2(1:s2(ii),ii)) - V_old(:,M2(1:s2(ii),ii)),2)
!     end if
!     residual = max(residual, maxval(abs(nu*dual_residual)))
!  end do
!za
!end subroutine residual_dual

!subroutine update_U(X,Lambda,U,V,xbar,M1,M2,s1,s2,mix1,mix2,p,q,nK,nu)
!  use constants
!  implicit none
!  integer :: ii,mix1,mix2,nK,p,q
!  integer :: M1(mix1,p), M2(mix2,p), s1(p), s2(p) 
!  real(kind=dble_prec) :: Lambda(q,nK), nu, U(q,p), V(q,nK), X(q,p), xbar(q), omega
!
!  omega = one/(one+p*nu)
!
!  U = X
!  do ii=1,p
!     if (s1(ii) > 0) then
!        U(:,ii) = U(:,ii) + sum(Lambda(:,M1(1:s1(ii),ii)),2) + nu*sum(V(:,M1(1:s1(ii),ii)),2)
!     end if
!     if (s2(ii) > 0) then
!        U(:,ii) = U(:,ii) - sum(Lambda(:,M2(1:s2(ii),ii)),2) - nu*sum(V(:,M2(1:s2(ii),ii)),2)
!     end if
!     U(:,ii) = omega*U(:,ii) + (one-omega)*xbar
!  end do
!
!end subroutine update_U

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

subroutine residual_dual(V,V_old,p,q,nK,nu,residual)
!
! This subroutine computes the L2 norm of the dual residual
!
  use hilbert
  use constants
  implicit none
  integer :: ii, p, q, nK
  integer :: ix(p), seq(p)
  real(kind=dble_prec) :: nu, residual, V(q,nK), V_old(q,nK)
  real(kind=dble_prec) :: L1(q), L2(q)

  seq = (/(ii, ii=1,p, 1)/)
! Loop unrolling: ii=1
  ii = 1
  call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
  L2 = sum(V(:,ix(1:(p-ii))),2)
  L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
  residual = norm(L2,2)**2
! Loop over ii=2:p-1
  do ii=2,p-1
     call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
     L1 = sum(V(:,ix(1:(ii-1))),2)
     L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
     call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
     L2 = sum(V(:,ix(1:(p-ii))),2)
     L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
     residual = residual + norm(L1-L2,2)**2
  end do
! Loop unrolling: ii=p
  ii = p
  call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
  L1 = sum(V(:,ix(1:(ii-1))),2)
  L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
  residual = nu*(residual + norm(L1,2)**2)
end subroutine residual_dual

!subroutine residual_dual(V,V_old,p,q,nK,nu,residual)
!
! This subroutine computes the L-infinity norm of the dual residual
!
!  use constants
!  implicit none
!  integer :: ii, p, q, nK
!  integer :: ix(p), seq(p)
!  real(kind=dble_prec) :: nu, residual, V(q,nK), V_old(q,nK)
!  real(kind=dble_prec) :: L1(q), L2(q)
!
!  seq = (/(ii, ii=1,p, 1)/)
! Loop unrolling: ii=1
!  ii = 1
!  call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
!  L2 = sum(V(:,ix(1:(p-ii))),2)
!  L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
!  residual = maxval(abs(L2))
! Loop over ii=2:p-1
!  do ii=2,p-1
!     call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
!     L1 = sum(V(:,ix(1:(ii-1))),2)
!     L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
!     call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
!     L2 = sum(V(:,ix(1:(p-ii))),2)
!     L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
!     residual = max(residual, maxval(abs(L1-L2)))
!  end do
! Loop unrolling: ii=p
!  ii = p
!  call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
!  L1 = sum(V(:,ix(1:(ii-1))),2)
!  L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
!  residual = nu*max(residual, maxval(abs(L1)))
!end subroutine residual_dual

subroutine update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
  use constants
  implicit none
  integer :: ii,p,q,nK
  integer :: ix(p), seq(p)
  real(kind=dble_prec) :: X(q,p), U(q,p), V(q,nK), Lambda(q,nK), xbar(q)
  real(kind=dble_prec) :: L1(q), L2(q), nu, omega
  seq = (/(ii, ii=1,p, 1)/)

  omega = one/(one + p*nu)
! Loop unrolling: ii=1
  ii = 1
  call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
  L2 = sum(Lambda(:,ix(1:(p-ii))),2)
  L2 = L2 + nu*sum(V(:,ix(1:(p-ii))),2)
  U(:,ii) = omega*(X(:,ii) + L2) + (one-omega)*xbar
! Loop over ii=2:p-1
  do ii=2,p-1
     call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
     L1 = sum(Lambda(:,ix(1:(ii-1))),2)
     L1 = L1 + nu*sum(V(:,ix(1:(ii-1))),2)
     call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
     L2 = sum(Lambda(:,ix(1:(p-ii))),2)
     L2 = L2 + nu*sum(V(:,ix(1:(p-ii))),2)
     U(:,ii) = omega*(X(:,ii) - L1 + L2) + (one-omega)*xbar
  end do
! Loop unrolling: ii = p
  ii = p
  call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
  L1 = sum(Lambda(:,ix(1:(ii-1))),2)
  L1 = L1 + nu*sum(V(:,ix(1:(ii-1))),2)
  U(:,ii) = omega*(X(:,ii) - L1) + (one-omega)*xbar

end subroutine update_U

subroutine update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
  use constants
  use proximal
  implicit none
  integer :: kk,nK,p,q
  integer :: i,j,ix(nK,2),type
  real(kind=dble_prec) :: gamma, nu
  real(kind=dble_prec) :: U(q,p), Lambda(q,nK), V(q,nK), w(nK)
  real(kind=dble_prec) :: z(q)
  
  do kk=1,nK
     i = ix(kk,1)
     j = ix(kk,2)
     z = U(:,i) - U(:,j) - (one/nu)*Lambda(:,kk)
     call prox(z,q,V(:,kk),w(kk)*gamma/nu,type)
!     call prox_L2(z,q,V(:,kk),w(kk)*gamma/nu)
  end do
  
end subroutine update_V

subroutine update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
  use constants
  implicit none
  integer :: nK, p, q
  integer :: ix(nK,2)
  real(kind=dble_prec) :: nu, Lambda(q,nK), U(q,p), V(q,nK)
  Lambda = Lambda - nu*(U(:,ix(:,1)) - U(:,ix(:,2)) - V)
end subroutine update_Lambda

!subroutine convex_cluster(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
!  use constants
!  implicit none
!  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
!  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
!  real(kind=dble_prec) :: gamma, Lambda(q,nK), nu, U(q,p), V(q,nK), V_old(q,nK), w(nK), X(q,p)
!  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
!  real(kind=dble_prec) :: xbar(q), rp, rd
!
!  print *, "Inside!"
!  print *, "iter", iter, "max_iter", max_iter
!  call get_xbar(X,xbar,p,q)
!  print *, xbar
!  do iter=1,max_iter
!     print *, iter
!     V_old = V
!     call update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
!     call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
!     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
!     call residual_primal(U,V,ix,p,q,nK,rp)
!     primal(iter) = rp
!     call residual_dual(V,V_old,p,q,nK,nu,rd)
!     dual(iter) = rd
!     print *, "primal residual", rp, "dual residual", rd     
!     if (max(rp,rd) < tol) exit
!  end do
!  print *, "iter", iter, "max_iter", max_iter
!  print *, "goodbye!"
!end subroutine convex_cluster

subroutine convex_clusterB_admm(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  use hilbert
  use constants
  implicit none
  integer :: iter, max_iter, nK, p, q, type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: gamma, Lambda(q,nK), nu, U(q,p), V(q,nK), V_old(q,nK), w(nK), X(q,p)
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: xbar(q), rp, rd, upsilon, LambdaNorm

  call get_xbar(X,xbar,p,q)

  upsilon = zero
  do iter=1,p
     upsilon = upsilon + norm(X(:,iter) - xbar,2)**2
  end do
  upsilon = sqrt(upsilon)

  do iter=1,max_iter
     V_old = V
     call update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
     call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
     LambdaNorm = sqrt(sum(Lambda**2))
     call residual_primal(U,V,ix,p,q,nK,rp)
     primal(iter) = rp
     call residual_dual(V,V_old,p,q,nK,nu,rd)
     dual(iter) = rd
!     if (max(rp,rd) < tol) exit
     if (rp*LambdaNorm + upsilon*rd < tol) exit
  end do
  iter = min(iter,max_iter)
end subroutine convex_clusterB_admm

subroutine convex_cluster_admm_acc(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  use hilbert
  use constants
  implicit none
  integer :: iter,q,p,nK,max_iter,type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), V_old(q,nK), w(nK), gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: xbar(q), rp, rd
  real(kind=dble_prec) :: alpha, alpha_old, S_old(q,nK), Lambda_old(q,nK), Ld(q,nK), r_last
  real(kind=dble_prec) :: upsilon, LambdaNorm
  
  call get_xbar(X,xbar,p,q)
  upsilon = zero
  do iter=1,p
     upsilon = upsilon + norm(X(:,iter) - xbar,2)**2
  end do
  upsilon = sqrt(upsilon)

  V_old = V
  Lambda_old = Lambda
  alpha_old = one
  r_last = huge(one)
  do iter=1,max_iter
     call update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
     call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
     LambdaNorm = sqrt(sum(Lambda**2))
     call residual_primal(U,V,ix,p,q,nK,rp)
     primal(iter) = rp
     call residual_dual(V,V_old,p,q,nK,nu,rd)
     dual(iter) = rd
!     if (max(rp,rd) < tol) exit
     if (rp*LambdaNorm + upsilon*rd < tol) exit
     if (max(rp,rd) > r_last) then
        alpha_old = one
        r_last = huge(one)
     else
        alpha = half*(one + sqrt(1 + four*(alpha_old**2)))
        V = V + ((alpha_old-one)/alpha)*(V-V_old)
        Lambda = Lambda + ((alpha_old-one)/alpha)*(Lambda-Lambda_old)
        r_last = max(rp,rd)
        alpha_old = alpha
     end if
     V_old = V
     Lambda_old = Lambda
  end do
end subroutine convex_cluster_admm_acc
