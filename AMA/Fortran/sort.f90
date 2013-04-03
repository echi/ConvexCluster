module sort
  use constants

contains
  recursive subroutine quick_sort(list)
    !
    !     This subroutine sorts the real array LIST into nonincreasing
    !     order by a combination of quick sort and exchange sort.
    !
    !     Taken from Ken's code.
    !
    implicit none
    integer :: i,j,n
    real(kind=dble_prec) :: split,uniform
    real(kind=dble_prec) :: list(:)
    !
    n = size(list)
    !
    !     For small lists do an exchange sort.
    !
    if (n<=6) then
       do i = 1,n
          do j = i+1,n
             if (list(i)<list(j)) then
                call swap(list(i),list(j))
             end if
          end do
       end do
       return
    end if
    !
    !     Otherwise carry out a quick sort.           
    !
    call random_number(uniform)
    i = int(n*uniform)+1
    split = list(i)
    list(i) = list(1)
    list(1) = split
    i = 1
    do j = 2,n
       if (list(j)>=split) then
          i = i+1
          call swap(list(i),list(j))
       end if
    end do
    list(1) = list(i)
    list(i) = split
    if (i>2) call quick_sort(list(1:i-1))
    if (i+1<n) call quick_sort(list(i+1:n))
  end subroutine quick_sort
  
  subroutine swap(a,b)
    !
    !     This subroutine swaps two reals.
    !
    implicit none
    real(kind=dble_prec) :: a,b,temp
    !
    temp = a
    a = b
    b = temp
  end subroutine swap
end module sort
