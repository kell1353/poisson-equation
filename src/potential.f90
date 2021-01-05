!-----------------------------------------------------------------------
!! Module: potential
!-----------------------------------------------------------------------
!! By: Noah, Austin, & Heng
!!
!! This module is responsible solving Poisson's equation for a charge
!! distribution held within a 2D conducting box. Ultimately, this 
!! module finds the electrostatic potential as a function of x and y.
!! To accomplish this, we construct the source term, rho_xy, and then
!! perform a 2D integration over x and y to find rho_mn. From rho_mn,
!! we find the fourier coefficients c_mn. We then sum c_mn multiplied
!! by the basis functions over m,n to construct the potential, phi(x,y).
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_coefficients
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! rho_xy
!! rho_mn_int
!! phi_potential
!!
!-----------------------------------------------------------------------
module potential
use types
use OMP_LIB
use quadrature, only : booles_quadrature
implicit none

private
public :: calculate_coefficients, phi_potential


contains

!-----------------------------------------------------------------------
!! Subroutine: rho_xy
!-----------------------------------------------------------------------
!! By: Noah, Austin, & Heng
!!
!! This function defines the source charge density. The exponentials
!! describe the rate of spatial decay as an observer moves away from 
!! the center. This function is separated form the main calculation 
!! to clean up the rho_mn_integrate function. 
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! x            real       x position to sample the charge distribution
!! y            real       y position to sample the charge distribution
!! center(:)    real       array containing the location of the center of the charge distribution
!! width(:)     real       array containing the length/width of the charge distribution    
!! rho_zero     real       charge density magnitude at (x_0, y_0)
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real(dp)        value of the function at the desired x and y values
!!
!-----------------------------------------------------------------------
real(dp) function rho_xy(x, y, center, width) result(k)
    implicit none
    real(dp), intent(in) :: center(1:2), width(1:2)
    real(dp), intent(in) :: x, y
    real(dp) :: r_x, r_y, x_zero, y_zero

    ! Charge distrubution width values 
    r_x = width(1)
    r_y = width(2)
    ! Charge distribution center values 
    x_zero = center(1)
    y_zero = center(2)

    k = (1/(pi*r_x*r_y))*exp(-((x-x_zero)/r_x)**2)*exp(-((y-y_zero)/r_y)**2)

end function rho_xy

!-----------------------------------------------------------------------
!! Subroutine: rho_mn_int
!-----------------------------------------------------------------------
!! By: Noah, Austin & Heng
!!
!! This subroutine computes the double integral over the box volume to 
!! produce rho_mn, the charge distribution coefficients. From this, we
!! may calculate the Fourier coefficients for each index, m and n. 
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! n_max        integer         The total number of sampling points
!! m            integer         the current index over the x dimension cosine function
!! n            integer         The current index over the y dimension cosine function
!! length(:)    real(dp)        array containing the length/width of the box 
!! center(:)    real(dp)        array containing the location of the center of the charge distribution
!! width(:)     real(dp)        array containing the length/width of the charge distribution   
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! rho          real(dp)        Result of integration from Booles Quadrature    
!-----------------------------------------------------------------------
real(dp) function rho_mn_int(length, center, width, n_max, n_samples, m, n) result(rho)
    implicit none
    real(dp), intent(in) :: length(1:2), center(1:2), width(1:2)
    integer, intent(in) :: n_max, n_samples, m, n
    real(dp), allocatable :: f_xy(:), f_x(:)
    real(dp) :: l_x, l_y, x, y, delta_x, delta_y
    integer :: i_y, i_x

    allocate(f_xy(1:n_samples))
    allocate(f_x(1:n_samples))

    ! Length and width of the box
    l_x = length(1)
    l_y = length(2)

    ! Define step sizes in the x and y directions 
    delta_x = l_x/(real(n_samples-1, kind=dp))
    delta_y = l_y/(real(n_samples-1, kind=dp))

    ! i_x index identifies index in the x-direction
    ! Outer loop
    do i_x = 1, n_samples
        x = (i_x - 1)*delta_x

        ! i_y index identifies index in the y-direction
        ! Inner loop
        do i_y = 1, n_samples 
            y = (i_y - 1)*delta_y

            ! Defines a 1-D array which contains the product of exponentials 
            ! and cosine functions for a specific x-value. The different values
            ! of the array correspond to the different y-values indexed by
            ! i_y.
            f_xy(i_y) = rho_xy(x, y, center, width)*cos((m*pi*x)/(l_x))*cos((n*pi*y)/(l_y))

       
        enddo

        ! Now, use booles_quadrature to itegrate over all y-values for
        ! a given i_x. Essentially, we are collapsing (integrating) 
        ! the 1-D array of y-values to form an element of the new array 
        ! where each element is indexed by i_x. Thus, we have a new
        ! array which contains the value of the products at each x
        ! indexed by i_x
        
        f_x(i_x) = booles_quadrature(f_xy, delta_y)

    enddo

    ! Perform the final integration over x.
    rho = booles_quadrature(f_x, delta_x)

end function rho_mn_int

!-----------------------------------------------------------------------
!! Subroutine: calculate_coefficients
!-----------------------------------------------------------------------
!! By: Noah, Austin, & Heng
!! 
!! This subroutine calculates the fourier coefficients needed for 
!! ultiamtely calculating the electrostatic potential, phi(x,y).
!!----------------------------------------------------------------------
!! Input:
!!
!! length(:)    real        contains the length/width of the conducting box 
!! center(:)    real        array containing the location of the center of the charge distribution
!! width(:)     real        array containing the length/width of the charge distribution 
!! n_max        integer     total number of sampling points
!! rho_zero     real        charge magnitude at x_0,y_0
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! num_threads      integer         The number of threads used in the parallelization 
!! c_fourier(:,:)   real(dp)        2D array containing the Fourier coefficients for each pair of m,n
!!
!-----------------------------------------------------------------------
subroutine calculate_coefficients(length, rho_zero, center, width, n_max, n_samples, c_fourier, num_threads)
    implicit none
    real(dp), intent(in) :: length(1:2), center(1:2), width(1:2), rho_zero
    integer, intent(in) :: n_max, n_samples
    integer, intent(out) :: num_threads
    real(dp), allocatable, intent(out) :: c_fourier(:,:)
    real(dp) :: l_x, l_y, r_x, r_y, norm, rho_mn
    integer :: m, n, OMP_GET_NUM_THREADS

    allocate(c_fourier(1:n_max, 1:n_max))

    ! Define length and width of box, and length and width of charge source
    l_x = length(1)
    l_y = length(2)
    r_x = width(1)
    r_y = width(2)

    ! Define normalization constant
    norm = 2/sqrt(l_x*l_y)

    ! Parallelization implementation
    !$omp parallel default(none) private(m,n,rho_mn) & 
    !$omp & shared(length,rho_zero,norm,center,width,l_x,l_y,r_x,r_y,c_fourier,n_max, n_samples, num_threads)

    ! Get number of threads for time tests
    num_threads = OMP_GET_NUM_THREADS()

    !$omp do schedule(dynamic)

    ! Calculate the fourier coefficients by constructing a 2D array
    ! the contains the values of the coefficients for each m,n.
    do m = 1, n_max
        do n = 1, n_max
            rho_mn = norm*rho_zero*rho_mn_int(length, center, width, n_max, n_samples, m, n)
            c_fourier(m, n) = 4*pi*rho_mn/((m*pi/l_x)**2 + (n*pi/l_y)**2)
        enddo
    enddo

    !$omp end do
    !$omp end parallel

end subroutine calculate_coefficients

!-----------------------------------------------------------------------
!! Function: phi_potential
!-----------------------------------------------------------------------
!! By: Noah, Austin, & Heng
!!
!! Computes the electrostatic potential at a given x,y by summing over
!! the indices m,n. The function receives the fourier coefficients 
!! array that contains the values for a given m,n.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_fourier(:,:)   real        2D array containing the Fourier coefficients (m,n)
!! length(:)        real        array containing the length/width of conducting box
!! x                real        x position of sampling point 
!! y                real        y position of sampling point 
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! r                real        value of potential for a given x,y
!-----------------------------------------------------------------------
real(dp) function phi_potential(c_fourier, length, x, y) result(r)
    implicit none
    real(dp), intent(in) :: c_fourier(:,:), length(1:2), x, y
    real(dp) :: l_x, l_y, norm
    integer :: m, n

    ! Box length and width
    l_x = length(1)
    l_y = length(2)

    ! Start sum at zero
    r = 0._dp

    ! Define normalization constant
    norm= 2/sqrt(l_x*l_y)

    ! Calculate phi(x,y) potential 
    do m=1,size(c_fourier, 1)
        do n=1,size(c_fourier, 2)
            r = r + c_fourier(m,n)*norm*cos(m*pi*x/l_x)*cos(n*pi*y/l_y)
        enddo
    enddo

end function phi_potential


end module potential