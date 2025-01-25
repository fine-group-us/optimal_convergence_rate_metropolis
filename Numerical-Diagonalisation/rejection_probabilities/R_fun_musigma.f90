! This function computes the rejection probability for the (mu, sigma)-family of jump distributions.
! The (mu, sigma)-family is based on a Gaussian distribution with mean mu and standard deviation sigma.
! The rejection probability reflects whether a proposed jump from position x to y + x is accepted, 
! based on the energy difference and Metropolis criterion.

! Additional functions that it employs:
!   - heaviside_fun.f90
!   - U_fun.f90
!   - w_musigma.f90

function R_fun_musigma(x, mu, sigma) result(rej)
  use omp_lib            ! Include OpenMP library for parallelization
  implicit none          ! Ensure all variables are explicitly declared

  ! Inputs
  real*8, intent(in) :: x   ! Current position of the particle
  real*8, intent(in) :: mu  ! Mean of the Gaussian jump distribution
  real*8, intent(in) :: sigma   ! Standard deviation of the Gaussian jump distribution

  ! Local variables
  integer :: i, N                       ! Loop counter and number of discretization points
  real*8 :: rej_part, rej_tot, rej      ! Partial and total rejection probabilities
  real*8 :: U_fun, heaviside_fun        ! Functions for energy and boundary constraints
  real*8 :: L, dl, da, y, sum           ! Parameters and intermediate variables
  real*8 :: diffU, w_musigma            ! Gaussian jump distribution and energy difference
  parameter(N = 1000)                   ! Number of discretization points

  double precision :: w_vec(N)          ! Array to store normalized jump weights

  L = 30.0d0        ! System size or boundary
  dl = 2.0d0 * L / N  ! Step size for the discretization of the distribution

  ! Parallel computation of the rejection probability
  !$omp parallel private(y, diffU, rej_part) shared(rej_tot)
  rej_part = 0.d0
  rej_tot = 0.d0
  !$omp do
  do i = 1, N
     y = dl * i - L                          ! Discretize the jump positions
     diffU = U_fun(y + x) - U_fun(x)         ! Compute the energy difference for the jump
     
     ! Compute rejection probability for the jump from x to y + x, considering
     ! the Gaussian jump distribution, energy difference, and system boundaries
     rej_part = rej_part + w_musigma(y, mu, sigma, L, N) * dl * (1.d0 - exp(-diffU) * &
                  heaviside_fun(L - (y + x)) * heaviside_fun((y + x) + L)) * &
                  heaviside_fun(diffU)
  enddo
  !$omp enddo
  !$omp critical
  rej_tot = rej_tot + rej_part             ! Safely update the shared rejection total
  !$omp end critical
  !$omp end parallel

  rej = rej_tot                            ! Assign the total rejection probability to the result

end function R_fun_musigma
