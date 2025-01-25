! This program generates the data for the heat maps depicted in Figs. 3 and 4 of the paper.
! It computes the eigenvalue (`lambda`), Inverse Participation Ratio (IPR), and quantum fidelity
! for a range of values of \( a \) and \( \alpha \), and stores the results in three separate output files.
! The data from these computations is used to create heat maps of these quantities as functions of \( a \) and \( \alpha \).

! Additional functions that it employs:
!   - eigen_num_dia_alphaa.f90
!   - IPR.f90
!   - Fidel.f90

program phase_diagram
    use omp_lib           ! Include OpenMP library for parallelization
    implicit none         ! Ensure all variables are explicitly declared

    ! Local variables and parameters
    integer :: i, j, Na, Nal, N   ! Loop counters and grid dimensions
    real*8 :: a_min, a_max, alpha_min, alpha_max, da, dal, a, alpha  ! Range and step size for \( a \) and \( \alpha \)
    real*8 :: L, ipr_num, fid_num, lambda  ! Variables for the system's properties (IPR, fidelity, eigenvalue)
    real*8 :: IPR, Fidel  ! Variables for the functions to compute IPR and fidelity

    ! Parameters
    parameter(N = 2000, L = 10.d0)                   ! Grid size and system size
    parameter(a_min = 1.8d0, a_max = 2.3d0, alpha_min = 0.8d0, alpha_max = 1.8d0)  ! Ranges for \( a \) and \( \alpha \)
    parameter(Na = 100, Nal = 100)                   ! Number of \( a \)-values and \( \alpha \)-values

    ! Array to store the eigenfunction
    double precision :: phi(N)

    ! Open output files to store eigenvalue, IPR, and fidelity data
    open(1, file = "lambda_map.dat", status = "unknown")
    open(2, file = "ipr_map.dat", status = "unknown")
    open(3, file = "fidel_map.dat", status = "unknown")

    ! Calculate the step size for \( a \) and \( \alpha \) values
    da = (a_max - a_min) / Na
    dal = (alpha_max - alpha_min) / Nal

    ! Loop over \( a \)-values and \( \alpha \)-values to compute and store the results
    do i = 1, Na
        a = da * i + a_min  ! Compute the current value of \( a \)
        do j = 1, Nal
            alpha = dal * j + alpha_min  ! Compute the current value of \( \alpha \)
            
            ! Call a subroutine to compute the eigenfunction corresponding to the parameters
            call eigen_num_dia_alphaa(lambda, phi, alpha, a, L, N)
            
            ! Compute the Inverse Participation Ratio (IPR) and quantum fidelity
            ipr_num = IPR(phi, N)
            fid_num = Fidel(phi, N, L, a)
            
            ! Write the results to the output files
            write(1, 199) a, alpha, lambda    ! Store eigenvalue data
            write(2, 199) a, alpha, ipr_num  ! Store IPR data
            write(3, 199) a, alpha, fid_num  ! Store fidelity data
            
            ! Print results for debugging or direct display
            print*, a, alpha
        enddo
    enddo

    ! Format for writing data to the files
    199 format(10(e14.7, 2x))

end program phase_diagram
