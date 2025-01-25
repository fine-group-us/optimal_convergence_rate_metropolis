program test5_dichotomy
    use omp_lib
    implicit none

    ! DECLARATION SECTION   

    integer :: i, j, k
    integer :: N_alpha, N, N_dic, Ns
    integer :: m

    real*8 :: epsilon, sum
    real*8 :: alpha
    real*8 :: a_min, a_max, alpha_min, alpha_max
    real*8 :: aa, ab, ac
    real*8 :: x, dx, L, da
    real*8 :: rprobc
    real*8 :: lambdac, lambdacdw
    real*8 :: fc, knum
    real*8 :: heaviside_fun, hermite_fun

    parameter(N=2000, N_alpha=20, N_dic=200, epsilon=0.0001d0, L=10.d0)
    parameter(a_min=2.d0, a_max=3.5d0, alpha_min=0.0d0, alpha_max=1.22d0)

    double precision :: wc(N), lambda_min(N_dic)
    double precision :: phi(N)

    ! INITIALISATION SECTION

    open(1,file="HRbranch.dat",status="unknown")

    ! CODE

    dx = L/dble(N)
    do i=1,N_alpha+1

        alpha = alpha_min + dble(i-1)*(alpha_max - alpha_min)/dble(N_alpha)
        aa = a_min
        ab = a_max
        
        do k=1,N_dic
            ac = (aa + ab)*0.5d0
            call eigen_num_dia_al(lambdac,phi,alpha,ac,L,N)
            da = ac/dble(N)
            sum = 0.d0
            do j=1,N
                x = da*(j)
                wc(j) = x**(alpha)
                sum = sum + 2.d0*wc(j)*da
            enddo
             
            do j=1,N
                wc(j) = wc(j)/sum
            enddo
            call rej_prob(rprobc,N,dx,wc)
            
            ! Change fc depending on the branch

             fc = lambdac - rprobc
             
            !lambda_min(k) = max(lambdac,rprobc)

            if (fc.ge.0.d0) then
                aa = ac
            else
                ab = ac
            endif
            if ((abs(fc).le.epsilon)) then
                !print*, k
                goto 10
             endif             
        enddo

10      continue

        !max(rprobc,lambdac)
        !max(lambdacdw,lambdac)
        !max(rprobc,lambdacdw)
        
        write(1,199) ac, alpha, max(rprobc,lambdac)
        print*, ac, alpha, lambdac

199     format(10(e14.7,2x))

    enddo






end program test5_dichotomy
