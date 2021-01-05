!!-----------------------------------------------------------------------
!! Program: poisson_potential
!!-----------------------------------------------------------------------
!! By: Noah, Austin, & Heng
!!
!! This program solves Poisson's equation within a conducting box and
!! obtains the electrostatic potential as a function of position x,y
!! within the box. In addition, this program implements a parallelization
!! via OpenMP in order to decrease the computation time. A script is 
!! provided to compare program run times as a function of thread number.
!! 
!!-----------------------------------------------------------------------
program poisson_potential

use types
use read_write, only : read_input, write_potential, print_wall_clock_time, write_wall_clock_time
use potential, only : calculate_coefficients

implicit none

real(dp) :: length(1:2)
real(dp) :: rho_zero, center(1:2), width(1:2)
integer :: n_max, n_samples, num_threads
character(len = 1024) :: file_name
character(len = 1024) :: time_file
logical :: run_write_time, series, parallel

integer :: count_1, count_2, count_rate, count_max

real(dp), allocatable :: c_fourier(:,:)

call system_clock(count_1, count_rate, count_max)

call read_input(length, rho_zero, center, width, n_max, n_samples, file_name)
call calculate_coefficients(length, rho_zero, center, width, n_max, n_samples, c_fourier, num_threads)

call write_potential(c_fourier, length, file_name)

call system_clock(count_2, count_rate, count_max)

call print_wall_clock_time(count_1, count_2, count_rate, num_threads)

! Make true if one wishes to run time tests and write results to a file.
! Otherwise, make false to go with one single run
run_write_time = .false.

! If run_write_time is true, decide if going with series or parallel time test:
! If series is true, comment out OMP_GET_NUM_THREADS() from potential module line 270
! and remove -fopenmp flag from makefile
series = .false.

! If parallel is true, make series false 
! Include -fopenmp flag in makefile
! Include OMP_GET_NUM_THREADS() in potential module line 270
parallel = .true.

if (run_write_time) then
    if (series) then
        time_file = 'serial_time.dat'
    else if (parallel) then
        time_file = 'thread_v_time.dat'
    endif
    call write_wall_clock_time(count_1, count_2, count_rate, time_file, num_threads)
endif

end program poisson_potential