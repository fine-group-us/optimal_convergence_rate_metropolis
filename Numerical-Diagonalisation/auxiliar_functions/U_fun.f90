! This function defines the harmonic potential U(x) = 0.5 * x^2.

function U_fun(x) result(u)
    implicit none                ! Ensure all variables are explicitly declared

    ! Input
    real*8, intent(in) :: x      ! Independent variable (position)

    ! Output
    real*8 :: u                  ! Harmonic potential value at x

    ! Compute the harmonic potential
    u = 0.5d0 * x**(2.d0)

end function U_fun
