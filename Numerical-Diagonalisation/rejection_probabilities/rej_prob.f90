subroutine rej_prob(rprob,N,dx,w)

    implicit none

    integer, intent(in) :: N
    real*8, intent(in) :: dx
    real*8, intent(out) :: rprob

    double precision, intent(in) :: w(N)

    integer :: i
    real*8 :: rr
    real*8 :: diffU
    real*8 :: heaviside_fun, U_fun

    ! Note: it computes R(x=0), not R(x) !!

    rprob = 0.d0
    do i=1,N
        rr = dx*(i)
        diffU = U_fun(rr) - U_fun(0.d0)
        rprob = rprob + 2.d0*w(i)*(1-Exp(-diffU))*heaviside_fun(diffu)*dx
    enddo

end subroutine rej_prob


