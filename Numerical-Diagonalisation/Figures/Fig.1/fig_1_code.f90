! This program generates the data depicted in Fig. 1 from the paper. 
! It calculates the eigenvalues of a system parameterized by the variable 'a' 
! from the (alpha,a)-family of jump-distributionsfor a range of values and stores
! them in a file. The results can then be used to produce a figure showing the behavior
! of the eigenvalues as a function of the parameter 'a'.

! Additional functions that it employs:
!   - eigen_num_dia_alphaa_multiple.f90

program fig_1_code
    use omp_lib            ! Include OpenMP library for parallelization
    implicit none          ! Ensure all variables are explicitly declared

    ! Local variables and parameters
    integer :: i, N, Na, Neigen, j        ! Loop counters, number of discretization points, and number of eigenvalues
    real*8 :: alpha, a, L, lambda1, lambda2, lambda3, dl, w, sum, x, R_fun_alpha, da
    real*8 :: amin, amax                     ! Minimum and maximum values of parameter 'a'

    ! Parameters
    parameter(N = 2000, alpha = 3.0d0, L = 10.d0, amin = 0.9d0, amax = 3.5d0, Na = 120)
    parameter(Neigen = 9)                   ! Number of eigenvalues to be calculated

    double precision :: lambdavec(N)         ! Array to store eigenvalues

    ! Open output file to store eigenvalues
    open(1, file = "eigenvalues-alpha_3.dat", status = "unknown")

    ! Step size for parameter 'a'
    da = dble((amax - amin) / Na)

    ! Loop over values of parameter 'a' to calculate eigenvalues
    do i = 0, Na
        a = amin + i * da                     ! Increment 'a' within the specified range

        ! Call a subroutine to calculate the eigenvalues for the current 'a'
        call eigen_num_dia_alphaa_multiple(lambdavec, alpha, a, L, N)

        ! Output the last two eigenvalues for the current 'a' (just for debugging)
        print*, lambdavec(N-1), lambdavec(N-2)

        ! Write the current value of 'a' and the first 'Neigen' eigenvalues to the file
        write(1, 199) a, (lambdavec(N-j), j = 1, Neigen)
    enddo

    ! Format for writing data to the file
    199 format(10(e14.7, 2x))  

end program fig_1_code

