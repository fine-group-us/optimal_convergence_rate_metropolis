! This program generates the data depicted in Fig. 5 from the paper.
! It calculates the quantum fidelity and the Inverse Participation Ratio (IPR)
! for a given set of parameters and stores the results in two separate files.
! These quantities are computed for different values of the parameter \( a \),
! and the results are written to two output files for further analysis and plotting.

! Additional functions that it employs:
!   - eigen_num_dia_alphaa.f90
!   - IPR.f90
!   - Fidel.f90

program fig5_Fid_IPR
    use omp_lib            ! Include OpenMP library for parallelization
    implicit none          ! Ensure all variables are explicitly declared

    ! Local variables and parameters
    integer :: i, N, Na    ! Loop counters and number of \( a \)-values
    real*8 :: alpha, a, L, lambda, dl, w_alphaa, sum, x, R_fun_alphaa, La, da
    real*8 :: ipr_num, IPR, fid_num, Fidel  ! Variables for IPR and quantum fidelity

    ! Parameters
    parameter(N = 2000, La = 2.5d0, Na = 100, alpha = 3.d0, L = 15.d0)  ! Set grid size, system size, and parameters

    ! Array to store eigenfunction
    double precision :: phi(N)

    ! Open output files to store IPR and fidelity data
    open(1, file = "IPR.dat", status = "unknown")
    open(2, file = "Fid.dat", status = "unknown")

    ! Calculate the step size for \( a \) values
    da = La / Na

    ! Loop over \( a \)-values to compute IPR and fidelity
    do i = 1, Na
        a = da * i + 1.d0   ! Compute the current value of \( a \)
        
        ! Call a subroutine to compute the eigenfunction corresponding to the parameters
        call eigen_num_dia_alphaa(lambda, phi, alpha, a, L, N)
        
        ! Compute the Inverse Participation Ratio (IPR) and quantum fidelity
        ipr_num = IPR(phi, N)  
        fid_num = Fidel(phi, N, L, a)
        
        ! Write the results to the output files
        write(1, 199) a, ipr_num
        write(2, 199) a, fid_num
        
        ! Print results for debugging or display
        print*, a, ipr_num, fid_num
    enddo

    ! Format for writing data to the file
    199 format(10(e14.7, 2x))  

end program fig5_Fid_IPR

