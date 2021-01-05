!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module contains two functions used to evaluate a volume integral. 
!! The integral method is Boole's quadrature which utilizes 5 equally
!! spaced points to approximate an integral across a given interval.
!! The corresponding integrals for each interval are summed to produce
!! the full integral.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function receives incoming arrays and breaks them into intervals
!! to be sent to the booles_rule function, which then formulates
!! and evaluates the expansion approximation given for each chunk. The 
!! function then sums over the results as well as provides a check that 
!! the incoming array is of correct size.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(:), delta_x
    integer :: fx_size, i
    fx_size = size(fx)

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. Implement a test 
    ! to make sure that the number of evaluated points in the fx array
    ! is the correct one

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4


    if (modulo(fx_size - 1, 4) /= 0) then
         print *, 'fx array size plus 1 has to be divisible by 4'
         print *, fx_size
         stop
     endif

    ! We define a smaller function that returns Boole's 
    ! five point rule and pass slices (1:5), (5:9), (9:13), 
    ! ... of fx to such function to then add all the results. 

    s = 0._dp

    do i = 1, ((fx_size-1)/4)
        s = s + booles_rule(fx((4*i - 3):(4*i + 1)), delta_x)
    enddo
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function receives discrete slices sent from booles_quadrature and
!! evaluates each location by utilizing the 5 equally spaced point formula.
!! Thus, it is evaluating f(x) for each x to then be sent back to 
!! booles_quadrature and summed over all points to give the full integral.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x
    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Making an additional test to make sure that the array
    ! received has 5 and only 5 points 

    if (fx_size /= 5) then
        print*, 'Array fx must have 5 points. Currently has', fx_size
        stop
    endif
    
     fx0 = fx(1)
     fx1 = fx(2)
     fx2 = fx(3)
     fx3 = fx(4)
     fx4 = fx(5)
    
     s = (delta_x*2/45._dp)*(7*fx0 + 32*fx1 +12*fx2 + 32*fx3 + 7*fx4)

end function booles_rule
end module quadrature