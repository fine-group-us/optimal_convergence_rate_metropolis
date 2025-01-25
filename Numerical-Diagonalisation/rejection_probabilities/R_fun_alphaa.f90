! This function computes the rejection probability for the (alpha, a)-family of jump-distributions.
! The rejection probability quantifies the likelihood of a proposed jump being rejected due to
! energy constraints in the presence of a harmonic potential. It incorporates the effects of the 
! jump distribution and the energy landscape defined by U_fun.

! Additional functions that it employs:
!   - heaviside_fun.f90
!   - U_fun.f90
!   - w_alphaa.f90

function R_fun_alphaa(x, alpha, a) result(rej)
  
  use omp_lib            ! Include OpenMP library for parallelization
  implicit none          ! Ensure all variables are explicitly declared

  ! Inputs
  real*8, intent(in) :: x       ! Position of the particle
  real*8, intent(in) :: alpha   ! Power-law exponent for the jump distribution
  real*8, intent(in) :: a       ! Cutoff for the jump distribution

  ! Local variables
  integer :: i, N                       ! Loop counter and number of discretization points
  real*8 :: rej_part, rej_tot, rej      ! Partial and total rejection probabilities
  real*8 :: U_fun, heaviside_fun        ! Functions used for energy and constraints
  real*8 :: L, dl, da, y, sum, diffU    ! Parameters and intermediate variables
  real*8 :: w_alphaa                    ! Jump-distribution
  parameter(N = 1000)                   ! Number of discretization points

  double precision :: w_vec(N)          ! Array to store normalized jump weights

  ! Step size for discretizing the jump distribution and potential
  da = 2.0d0 * a / N
  L = 30.0d0        ! System size or boundary
  dl = 2.0d0 * L / N

  ! Compute the unnormalized weights for the (alpha, a)-family of jump distributions
  sum = 0.d0
  do i = 1, N
     y = da * i - a
     w_vec(i) = abs(y)**alpha
     sum = sum + w_vec(i) * da
  enddo

  ! Normalize the weights
  do i = 1, N
     w_vec(i) = w_vec(i) / sum
  enddo

  ! Parallel computation of the rejection probability
  !$omp parallel private(y, diffU, rej_part) shared(rej_tot)
  rej_part = 0.d0
  rej_tot = 0.d0
  !$omp do
  do i = 1, N
     y = da * i - a
     diffU = U_fun(y + x) - U_fun(x)   ! Energy difference for the proposed jump
     rej_part = rej_part + w_vec(i) * da * (1.d0 - exp(-diffU) * &
                  heaviside_fun(L - (y + x)) * heaviside_fun((y + x) + L)) * &
                  heaviside_fun(diffU)      ! Include constraints and rejection criteria
  enddo
  !$omp enddo
  !$omp critical
  rej_tot = rej_tot + rej_part          ! Safely update the shared rejection total
  !$omp end critical
  !$omp end parallel

  rej = rej_tot                         ! Assign the total rejection probability to the result

end function R_fun_alphaa

