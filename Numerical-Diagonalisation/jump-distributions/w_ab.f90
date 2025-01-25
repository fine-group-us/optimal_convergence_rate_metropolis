! This function defines the jump-distribution function wab(x) for the (a,b)-family.
! The function describes the weight distribution of jumps occurring within a bounded
! interval [b, a], where b < a. The heaviside function ensures that the distribution
! is non-zero only within this interval, while the result is normalized by a computed
! sum to ensure the total weight is consistent.

! Additional functions that it employs:
!    - heaviside_fun.f90

function w_ab(x,b,a,L,N) result(w_dis)
  implicit none                ! Ensure all variables are explicitly declared

  ! Inputs
  integer, intent(in) :: N     ! Number of discretization points
  real*8, intent(in) :: x      ! The distance or relative position
  real*8, intent(in) :: b      ! Lower bound of the jump interval
  real*8, intent(in) :: a      ! Upper bound of the jump interval
  real*8, intent(in) :: L      ! Half-width of the domain

  ! Local variables
  integer :: i                 ! Loop counter
  real*8 :: w_dis              ! The resulting weight for the given x
  real*8 :: sum                ! Normalization factor for the weight
  real*8 :: y, dl              ! Current position and discretization step
  real*8 :: heaviside_fun      ! External heaviside function

  dl = 2.d0 * L / N            ! Compute the discretization step size

  ! Calculate the normalization sum over the domain
  sum = 0.d0                   ! Initialize the sum
  do i = 1, N
    y = dl * i - L             ! Compute the position in the domain
    sum = sum + heaviside_fun(abs(y) - b) * heaviside_fun(a - abs(y)) * dl
  enddo

  ! Compute the weight function wab(x)
  ! Only consider values within the interval [b, a] using heaviside functions
  w_dis = heaviside_fun(abs(x) - b) * heaviside_fun(a - abs(x)) / sum

end function w_ab

