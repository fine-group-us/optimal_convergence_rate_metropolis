! This program generates the data depicted in Fig. 2 from the paper.
! It calculates the eigenvector corresponding to a given set of parameters
! for the \( (a, b) \)-family of jump-distributions. The program saves the
! eigenvector data to a file, which is later used to generate a plot showing
! the spatial profile of the eigenvector.

! Additional functions that it employs:
!   - eigen_num_dia_ab.f90

program fig2_ab
    use omp_lib            ! Include OpenMP library for parallelization
    implicit none          ! Ensure all variables are explicitly declared

    ! Local variables and parameters
    integer :: i, N        ! Loop counter and grid size
    real*8 :: b, a, L, lambda, dl, w_ab, sum, x, R_fun_ab
    double precision :: phi(N)     ! Array to store the eigenvector

    ! Parameters
    parameter(N = 2000, a = 2.d0, b = 1.7d0, L = 10.d0)  ! Set grid size, parameter values, and system size

    ! Open output file to store eigenvector data
    open(1, file = "eigenvector_ab.dat", status = "unknown")

    ! Calculate grid spacing
    dl = 2.d0 * L / N

    ! Call a subroutine to compute the eigenvector corresponding to the given parameters
    call eigen_num_dia_ab(lambda, phi, b, a, L, N)

    ! Output the eigenvalue for debugging purposes
    print*, lambda

    ! Write the spatial profile of the eigenvector to the output file
    do i = 1, N
        write(1, 199) dl * i - L, phi(i)    ! Write position and corresponding eigenvector value
    enddo

    ! Format for writing data to the file
    199 format(10(e14.7, 2x))  

end program fig2_ab
