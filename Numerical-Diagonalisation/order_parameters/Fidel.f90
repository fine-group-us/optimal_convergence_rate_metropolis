! This function computes the quantum fidelity between a given input eigenfunction (phi) 
! and a set of Charge Density Wave (CDW) eigenfunctions. The fidelity measures the overlap 
! between the input eigenfunction and the normalized CDW eigenfunctions in momentum space, 
! returning the maximum fidelity value.

function Fidel(phi,N,L,a) result(fid)
    use omp_lib            ! Include OpenMP library for parallelization
    implicit none          ! Ensure all variables are explicitly declared

    ! Inputs
    integer, intent(in) :: N       ! Number of discretization points
    real*8, intent(in) :: L        ! Half-width of the domain
    real*8, intent(in) :: a        ! Lattice spacing
    double precision, intent(in) :: phi(N)  ! Input eigenfunction in real space

    ! Local variables
    integer :: i, j, Nk            ! Loop counters and number of momentum points
    real*8 :: fid                 ! Final fidelity result
    real*8 :: k, k_low, k_up      ! Momentum values and limits
    real*8 :: pi, dk, dl          ! Constants: pi, momentum step, and real-space step
    real*8 :: x                   ! Real-space position
    real*8 :: fid_aux, norm       ! Auxiliary variables for fidelity and normalization
    real*8 :: norm_tot, norm_aux  ! Totals and auxiliary values for normalization
    real*8 :: fid_tot             ! Total fidelity in the parallel region

    parameter(Nk = 1000)          ! Number of momentum points in the summation

    double precision :: fid_vec(Nk)  ! Array to store fidelity values for each k

    dl = 2.d0 * L / N             ! Real-space discretization step
    pi = 4.d0 * datan(1.d0)       ! Value of π (using arctan identity)
    k_low = 0.d0                  ! Minimum momentum value
    k_up = 4.d0 * pi / a          ! Maximum momentum value
    dk = (k_up - k_low) / Nk      ! Momentum discretization step

    ! Loop over momentum values
    do i = 1, Nk
        k = dk * i + k_low        ! Compute the current momentum value

        ! Calculate normalization factor for the CDW eigenfunction
        !$omp parallel private(x, norm_aux) shared(norm_tot)
        norm_aux = 0.d0           ! Initialize local normalization
        norm_tot = 0.d0           ! Initialize shared normalization
        !$omp do
        do j = 1, N
            x = dl * j - L        ! Compute the real-space position
            norm_aux = norm_aux + dl * exp(-x**2 / 2.d0) * cos(k * x)**2
        enddo
        !$omp enddo
        !$omp critical
        norm_tot = norm_tot + norm_aux
        !$omp end critical
        !$omp end parallel
        norm = sqrt(norm_tot)     ! Compute the normalization

        ! Calculate fidelity for the input eigenfunction with the current k
        !$omp parallel private(fid_aux, x, norm_aux) shared(fid_tot)
        fid_aux = 0.d0            ! Initialize local fidelity
        fid_tot = 0.d0            ! Initialize shared fidelity
        !$omp do
        do j = 1, N
            x = dl * j - L        ! Compute the real-space position
            fid_aux = fid_aux + dl * cos(k * x) * exp(-x**2 / 4.d0) * phi(j) / norm
        enddo
        !$omp enddo
        !$omp critical
        fid_tot = fid_tot + fid_aux
        !$omp end critical
        !$omp end parallel
        fid_vec(i) = abs(fid_tot) ! Store the absolute fidelity for this k
    enddo

    fid = maxval(fid_vec)         ! Find the maximum fidelity across all k values

end function Fidel

