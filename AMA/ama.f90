subroutine loss_primal(X,U,gamma,ix,p,q,nK,w,output,type)
  use constants
  use hilbert
  implicit none
  integer :: k,p,q,nK,type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: DU(q), X(q,p), U(q,p), gamma, w(nK), output, penalty

  penalty = zero
  do k=1,nK
     DU = U(:,ix(k,1)) - U(:,ix(k,2))
     penalty = penalty + w(k)*norm(DU,type)
   end do
  output = half*sum((X - U)**2) + gamma*penalty
end subroutine loss_primal

subroutine loss_primal_delta(X,Delta,gamma,ix,p,q,nK,w,output,type)
  use constants
  use hilbert
  implicit none
  integer :: i,j,k,p,q,nK,type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: X(q,p), Delta(q,p), gamma, w(nK)
  real(kind=dble_prec) :: g(q), output, penalty

  penalty = zero
  do k=1,nK
     i = ix(k,1)
     j = ix(k,2)
     g = X(:,i) + Delta(:,i) - X(:,j) - Delta(:,j)
     penalty = penalty + w(k)*norm(g,type)
  end do
  output = half*sum(Delta**2) + gamma*penalty
end subroutine loss_primal_delta

subroutine loss_dual_delta(X,Lambda,Delta,ix,p,q,nK,output)
  use constants
  use hilbert
  implicit none
  integer :: ii,nK,p,q
  integer :: ix(nK,2)
  real(kind=dble_prec) :: Delta(q,p),Lambda(q,nK),X(q,p)
  real(kind=dble_prec) :: output

  output = -half*sum(Delta**2)
  do ii=1,nK
     output = output - inner_product(X(:,ix(ii,1))-X(:,ix(ii,2)),Lambda(:,ii))
  end do

end subroutine loss_dual_delta

subroutine loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,output)
  use constants
  use hilbert
  implicit none
  integer :: ii, mix1, mix2, nK, p, q
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: output, Lambda(q,nK), X(q,p)
  real(kind=dble_prec) :: first_term, second_term, L1(q), L2(q)
 
  first_term = zero
  do ii=1,p
     L1 = zero
     L2 = zero
     if (s1(ii) > 0) then
        L1 = sum(Lambda(:,M1(1:s1(ii),ii)),2) 
     end if
     if (s2(ii) > 0) then
        L2 = sum(Lambda(:,M2(1:s2(ii),ii)),2)
     end if
     first_term = first_term + norm(L1-L2,2)**2
  end do
  second_term = zero
  do ii=1,nK
     second_term = second_term + inner_product(X(:,ix(ii,1))-X(:,ix(ii,2)),Lambda(:,ii))
  end do
  output = -half*first_term - second_term

end subroutine loss_dual

subroutine update_Delta(Lambda,Delta,M1,M2,s1,s2,mix1,mix2,p,q,nK)
  use constants
  implicit none
  integer :: ii,mix1,mix2,nK,p,q
  integer :: s1(p), s2(p), M1(mix1,p), M2(mix2,p)
  real(kind=dble_prec) :: d(q), Delta(q,p), Lambda(q,nK)

  do ii=1,p
     d = zero
     if (s1(ii) > 0) then
        d = d + sum(Lambda(:,M1(1:s1(ii),ii)),2)
     end if
     if (s2(ii) > 0) then
        d = d - sum(Lambda(:,M2(1:s2(ii),ii)),2)
     end if
     Delta(:,ii) = d
  end do
end subroutine update_Delta

subroutine update_U(X,Lambda,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
  use constants
  implicit none
  integer :: ii,mix1,mix2,nK,p,q
  integer :: s1(p), s2(p), M1(mix1,p), M2(mix2,p)
  real(kind=dble_prec) :: Lambda(q,nK), U(q,p), X(q,p)
  
  do ii=1,p
     U(:,ii) = X(:,ii)
     if (s1(ii) > 0) then
        U(:,ii) = U(:,ii) + sum(Lambda(:,M1(1:s1(ii),ii)),2)
     end if
     if (s2(ii) > 0) then
        U(:,ii) = U(:,ii) - sum(Lambda(:,M2(1:s2(ii),ii)),2)
     end if
  end do

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
!
! Uncomment next line for L2-norm
!
     call prox_L2(z,q,V(:,kk),w(kk)*gamma/nu)
!     call prox_L1(z,q,V(:,kk),w(kk)*gamma/nu)  
     call prox(z,q,V(:,kk),w(kk)*gamma/nu,type)
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

subroutine update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
  use constants
  use proximal
  implicit none
  integer :: ii, nK, p, q, type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: gamma, Lambda(q,nK), nu, U(q,p), w(nK), Z(q,nK), buffer(q)
  Z = nu*(U(:,ix(:,1)) - U(:,ix(:,2))) - Lambda
  do ii=1,nK
     call prox(Z(:,ii),q,buffer,gamma*w(ii),type)
     Lambda(:,ii) = -buffer
  end do
end subroutine update_LambdaDP

subroutine update_Lambda_Delta(Lambda,G,nu,gamma,ix,q,p,nK,w,type)
!
! G equals X + Delta
!
  use constants
  use proximal
  implicit none
  integer :: ii,nK,p,q,type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: G(q,nK),Lambda(q,nK),w(nK)
  real(kind=dble_prec) :: backward(q),forward(q),gamma,nu
  do ii=1,nK
     forward = nu*(G(:,ix(ii,1)) - G(:,ix(ii,2))) - Lambda(:,ii)
!     call proj_L2(forward,q,backward,gamma*w(ii))
     call prox(forward,q,backward,gamma*w(ii),type)
     Lambda(:,ii) = -backward
  end do
end subroutine update_Lambda_Delta

subroutine convex_cluster(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), w(nK), gamma, nu, eta
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: fp, fd

  do iter=1,max_iter
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
!     call update_V(U,Lambda_old,V,w,gamma,nu,ix,q,p,nK)
!     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
     call update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
     if (fp-fd < tol) exit
  end do
!
! I add an offset to the type since V is updated using proximal mapping of
! the conjugate function.
!
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type+3)
end subroutine convex_cluster

subroutine convex_cluster_fista_fixed(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
!
! This subroutine solves the convex clustering problem with a fixed step size using FISTA extrapolation steps.
!
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), w(nK), gamma, nu, eta
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), S(q,nK)
  real(kind=dble_prec) :: fp, fd

  do iter=1,2
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
  end do
  do iter=3,max_iter
     S = Lambda + (dble(iter-2)/dble(iter+1))*(Lambda - Lambda_old)     
     call update_U(X,S,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call update_LambdaDP(S,U,nu,gamma,ix,q,p,nK,w,type)
     Lambda_old = Lambda
     Lambda = S
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
!     if (fp-fd < tol*(one+half*(fp+fd))) exit
     if (fp-fd < tol) exit
  end do
!
! I add an offset to the type since V is updated using proximal mapping of
! the conjugate function.
!
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type+3)
end subroutine convex_cluster_fista_fixed

subroutine convex_cluster_backtrack(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: Lambda(q,nK), U(q,p), V(q,nK), X(q,p), w(nK)
  real(kind=dble_prec) :: eta, gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: f, f_last, fp, fd, qloss

  do iter=1,max_iter
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call loss_dual(X,Lambda_old,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f_last)
     do
        Lambda = Lambda_old
        call update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
        Ld = Lambda - Lambda_old
        qloss = sum(Ld*(U(:,ix(:,1)) - U(:,ix(:,2)))) + (half/nu)*sum(Ld**2) - f_last
        call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f)
        if (qloss+f.ge.zero) then
           exit
        end if
        nu = nu/eta
     end do
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
!     if (fp-fd < tol) exit
     if (fp-fd < tol*(one + half*(fp+fd))) exit
  end do
  iter = min(iter,max_iter)
!
! I add an offset to the type since V is updated using proximal mapping of
! the conjugate function.
!
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type+3)
end subroutine convex_cluster_backtrack

subroutine check_dual_feasibility(Lambda, q, nK, w, gamma, infeasible_flag)
  use constants
  use hilbert
  implicit none
  integer :: ii, nK, q
  real(kind=dble_prec) :: Lambda(q,nK), w(nK), gamma
  logical :: infeasible_flag
  
  infeasible_flag = .FALSE.
  do ii=1,nK
     if (norm(Lambda(:,ii),2) > w(ii)*gamma) then
        infeasible_flag = .TRUE.
        exit
     end if
  end do
end subroutine check_dual_feasibility

subroutine convex_cluster_fista_backtrack(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: Lambda(q,nK), U(q,p), V(q,nK), X(q,p), w(nK)
  real(kind=dble_prec) :: eta, gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: f, f_last, fp, fd, qloss
  real(kind=dble_prec) :: alpha, alpha_old, S_old(q,nK), S(q,nK)

  alpha_old = one
  S_old = Lambda
  Lambda_old = S_old
  do iter=1,max_iter
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call loss_dual(X,Lambda_old,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f_last)
     do
        Lambda = Lambda_old
        call update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
        Ld = Lambda - Lambda_old
        qloss = sum(Ld*(U(:,ix(:,1)) - U(:,ix(:,2)))) + (half/nu)*sum(Ld**2) - f_last
        call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f)
        if (qloss+f.ge.zero) then
           exit
        end if
        nu = nu/eta
     end do
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
     if (fp-fd < tol*(one + half*(fp+fd))) exit
     S = Lambda
     alpha = half*(one + sqrt(one + four*(alpha_old**2)))
     Lambda = S + ((alpha_old-one)/alpha)*(S-S_old)
     S_old = S
     alpha_old = alpha
     Lambda_old = Lambda
  end do
  iter = min(iter,max_iter)
!
! I add an offset to the type since V is updated using proximal mapping of
! the conjugate function.
!
  call update_V(U,S_old,V,w,gamma,nu,ix,q,p,nK,type+3)
end subroutine convex_cluster_fista_backtrack

subroutine convex_cluster_fista_backtrackB(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), w(nK), gamma, nu, eta
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: S(q,nK), fp, fd, f, qloss
  real(kind=dble_prec) :: f_last

!  print *, nu

  S = Lambda
  Lambda_old = S
  do iter=1,max_iter
     S = Lambda + (dble(iter-2)/dble(iter+1))*(Lambda-Lambda_old)
     Lambda_old = Lambda
     call update_U(X,S,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call loss_dual(X,S,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f_last)
     do
        Lambda = S
        call update_LambdaDP(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
        Ld = Lambda - S
        qloss = -f_last + sum(Ld*(U(:,ix(:,1)) - U(:,ix(:,2)))) + (half/nu)*sum(Ld**2)
        call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f)
        if (qloss+f.ge.zero) exit
        nu = nu/eta
     end do
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
!     print *, iter, nu, fd, fp, fp - fd
     if (fp-fd < tol) exit
  end do
  iter = min(iter,max_iter)
!
! I add an offset to the type since V is updated using proximal mapping of
! the conjugate function.
!
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type+3)
end subroutine convex_cluster_fista_backtrackB
