! This function defines the Heaviside step function H(x), 
! which is a piecewise function commonly used in mathematical and physical applications 
! to model thresholds or activation functions. The function returns:
!   - H(x) = 1 if x > 0
!   - H(x) = 1 if x = 0 (conventionally included)
!   - H(x) = 0 if x < 0

function heaviside_fun(x) result(h)
    implicit none                ! Ensure all variables are explicitly declared

    ! Input
    real*8, intent(in) :: x      ! Independent variable (input to the Heaviside function)

    ! Output
    real*8 :: h                  ! Heaviside function value at x

    ! Local variable
    real*8 :: dif                ! Difference variable (set to x for clarity)

    dif = x                      ! Assign x to dif

    ! Define the Heaviside function behavior
    if (dif.gt.0.d0) then
        h = 1.d0                 ! H(x) = 1 for x > 0
    elseif (dif.eq.0.d0) then
        h = 1.d0                 ! H(x) = 1 for x = 0 (conventional choice)
    else
        h = 0.d0                 ! H(x) = 0 for x < 0
    endif

end function heaviside_fun


