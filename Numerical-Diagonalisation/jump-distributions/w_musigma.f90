! This function defines the jump-distribution function w_musigma(x) for the (\mu, \sigma)-family.
! The function describes the probability distribution of jumps characterized by a Gaussian-like
! profile centered at \(\mu\) and with width parameter \(\sigma\). The result is normalized over
! the specified domain.

! Additional functions that it employs:
!    - heaviside_fun.f90

function w_musigma(x,mu,sigma,L,N) result(w_dis)
  implicit none                ! Ensure all variables are explicitly declared

  ! Inputs
  integer, intent(in) :: N     ! Number of discretization points
  real*8, intent(in) :: x      ! The distance or relative position
  real*8, intent(in) :: mu     ! Center of the Gaussian-like distribution
  real*8, intent(in) :: sigma  ! Width parameter (standard deviation) of the distribution
  real*8, intent(in) :: L      ! Half-width of the domain

  ! Local variables
  integer :: i                 ! Loop counter
  real*8 :: w_dis              ! The resulting weight for the given x
  real*8 :: sum                ! Normalization factor for the weight
  real*8 :: y, dl              ! Current position and discretization step
  real*8 :: heaviside_fun      ! External heaviside function (not used directly here)

  dl = 2.d0 * L / N            ! Compute the discretization step size

  ! Calculate the normalization sum over the domain
  sum = 0.d0                   ! Initialize the sum
  do i = 1, N
    y = dl * i - L             ! Compute the position in the domain
    sum = sum + dl * exp(-0.5d0 * (abs(y) - mu)**2 / sigma**2)
  enddo

  ! Compute the weight function w_musigma(x)
  ! The Gaussian-like distribution is evaluated at x and normalized by the sum
  w_dis = exp(-0.5d0 * (abs(x) - mu)**2 / sigma**2) / sum

end function w_musigma

