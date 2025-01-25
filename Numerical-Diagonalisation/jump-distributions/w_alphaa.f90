! This function defines the jump-distribution function w_alphaa(x) for the (alpha, a)-family.
! The function describes the probability distribution of jumps characterized by a power-law 
! dependence with exponent \\alpha\), constrained to an interval \([0, a]\). The heaviside 
! function ensures that jumps are restricted to the specified range, and the result is normalized.

! Additional functions that it employs:
!    - heaviside_fun.f90

function w_alphaa(x,alpha,a,L,N) result(w_dis)
  implicit none                ! Ensure all variables are explicitly declared

  ! Inputs
  integer, intent(in) :: N     ! Number of discretization points
  real*8, intent(in) :: x      ! The distance or relative position
  real*8, intent(in) :: alpha  ! Power-law exponent for the jump distribution
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
    sum = sum + abs(y)**(alpha) * heaviside_fun(a - abs(y)) * dl
  enddo

  ! Compute the weight function w_alphaa(x)
  ! Only consider values within the interval [0, a] using heaviside functions
  w_dis = abs(x)**(alpha) * heaviside_fun(a - abs(x)) / sum

end function w_alphaa

