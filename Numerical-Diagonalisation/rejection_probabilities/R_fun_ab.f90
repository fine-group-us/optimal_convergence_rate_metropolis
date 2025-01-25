! This function computes the rejection probability for the (a, b)-family of jump-distributions.
! The (a, b)-family restricts jumps to the interval [b, a] and defines a uniform distribution 
! over this range, with additional constraints imposed by the energy difference (Metropolis criterion)
! and system boundaries.

! Additional functions that it employs:
!   - heaviside_fun.f90
!   - U_fun.f90
!   - w_ab.f90

function R_fun_ab(x, b, a) result(rej)
  use omp_lib            ! Include OpenMP library for parallelization
  implicit none          ! Ensure all variables are explicitly declared

  ! Inputs
  real*8, intent(in) :: x  ! Position of the particle
  real*8, intent(in) :: b  ! Inner cutoff of the jump distribution
  real*8, intent(in) :: a  ! Outer cutoff of the jump distribution

  ! Local variables
  integer :: i, N                         ! Loop counter and number of discretization points
  real*8 :: rej_part, rej_tot, rej        ! Partial and total rejection probabilities
  real*8 :: U_fun, heaviside_fun          ! Functions used for energy and constraints
  real*8 :: L, dl, da, y, sum             ! Parameters and intermediate variables
  real*8 :: diffU1, diffU2, y1, y2        ! Variables for split jump cases
  real*8 :: w_ab                           ! Jump distribution function
  parameter(N = 1000)                     ! Number of discretization points

  double precision :: w_vec(N)            ! Array to store normalized jump weights

  ! Step size for the (a, b)-jump distribution
  da = (a - b) / N
  L = 30.0d0          ! System size or boundary
  dl = 2.0d0 * L / N

  ! Compute the unnormalized weights for the (a, b)-family of jump distributions
  sum = 0.d0
  do i = 1, N
     y = dl * i - L
     w_vec(i) = heaviside_fun(abs(y) - b) * heaviside_fun(a - abs(y))
     sum = sum + w_vec(i) * dl
  enddo

  ! Normalize the weights
  do i = 1, N
     w_vec(i) = w_vec(i) / sum
  enddo

  ! Parallel computation of the rejection probability
  !$omp parallel private(y1, y2, diffU1, diffU2, rej_part) shared(rej_tot)
  rej_part = 0.d0
  rej_tot = 0.d0
  !$omp do
  do i = 1, N + 1
     y1 = da * (i - 1) - a                  ! First case: negative jump
     diffU1 = U_fun(y1 + x) - U_fun(x)      ! Energy difference for the negative jump
     y2 = da * (i - 1) + b                  ! Second case: positive jump
     diffU2 = U_fun(y2 + x) - U_fun(x)      ! Energy difference for the positive jump

     ! Compute rejection for negative and positive jumps using the normalized jump distribution
     rej_part = rej_part + w_ab(y1, b, a, L, N) * da * (1.d0 - exp(-diffU1) * &
                  heaviside_fun(L - (y1 + x)) * heaviside_fun((y1 + x) + L)) * &
                  heaviside_fun(diffU1)

     rej_part = rej_part + w_ab(y2, b, a, L, N) * da * (1.d0 - exp(-diffU2) * &
                  heaviside_fun(L - (y2 + x)) * heaviside_fun((y2 + x) + L)) * &
                  heaviside_fun(diffU2)
  enddo
  !$omp enddo
  !$omp critical
  rej_tot = rej_tot + rej_part             ! Safely update the shared rejection total
  !$omp end critical
  !$omp end parallel

  rej = rej_tot                            ! Assign the total rejection probability to the result

end function R_fun_ab
