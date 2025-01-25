! This function computes the Inverse Participation Ratio (IPR) for a given eigenfunction (phi).
! The IPR is a measure of the spatial localization of the eigenfunction, with larger values
! indicating greater localization. It is particularly useful in studying phenomena such as
! Anderson localization and disorder-induced transitions.

function IPR(phi, N) result(ipr_val)
    use omp_lib            ! Include OpenMP library for parallelization
    implicit none          ! Ensure all variables are explicitly declared

    ! Inputs
    integer, intent(in) :: N             ! Number of discretization points
    double precision, intent(in) :: phi(N)  ! Input eigenfunction in real space

    ! Local variables
    integer :: i                         ! Loop counter
    real*8 :: sum1, sum2                 ! Summation variables for numerator and denominator
    real*8 :: sum1_par, sum2_par         ! Partial sums for parallel computation
    real*8 :: ipr_val                    ! Final IPR value

    ! Initialize summations
    sum1 = 0.d0
    sum2 = 0.d0

    !$omp parallel private(sum1_par, sum2_par) shared(sum1, sum2)
    sum1_par = 0.d0       ! Local partial sum for |phi(i)|^4
    sum2_par = 0.d0       ! Local partial sum for |phi(i)|^2
    !$omp do
    do i = 1, N
        sum1_par = sum1_par + abs(phi(i))**4.d0    ! Compute |phi(i)|^4
        sum2_par = sum2_par + abs(phi(i))**2.d0    ! Compute |phi(i)|^2
    enddo
    !$omp enddo

    !$omp critical
    sum1 = sum1 + sum1_par    ! Safely update the shared numerator sum
    sum2 = sum2 + sum2_par    ! Safely update the shared denominator sum
    !$omp end critical
    !$omp end parallel

    ! Compute the Inverse Participation Ratio (IPR)
    ipr_val = sum1 / sum2**2.d0

end function IPR
