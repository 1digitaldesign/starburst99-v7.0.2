!> Submodule for error handling functionality in galaxy_module
!>
!> This submodule implements all error handling procedures
!> defined in the parent galaxy_module.
submodule (galaxy_module) galaxy_module_error
   use, intrinsic :: iso_fortran_env, only: real32, error_unit
   implicit none

contains

   !> Handle errors with consistent formatting and optional termination
   module procedure error_handler
      character(len=:), allocatable :: prefix
      character(len=19) :: timestamp
      logical :: should_stop
      
      ! Determine if error is fatal
      should_stop = .false.
      if (present(fatal)) should_stop = fatal
      
      ! Create error prefix based on severity
      if (should_stop) then
         prefix = "FATAL ERROR"
      else
         prefix = "WARNING"
      end if
      
      ! Get current timestamp
      call date_and_time(time=timestamp)
      
      ! Write error message to stderr with timestamp
      write(error_unit, '(A,A,"[",A,"]: ",A,A,": ",A)') &
            "Galaxy: ", prefix, timestamp, " in ", routine, message
      
      ! Store the latest error message in the module variable
      error_message = trim(prefix) // " in " // trim(routine) // ": " // trim(message)
      
      ! Terminate program if error is fatal
      if (should_stop) then
         error stop
      end if
   end procedure error_handler
   
   !> Implement the linear interpolation function
   module procedure flin
      if (abs(x2 - x1) < epsilon(1.0_real32)) then
         y = y1  ! Avoid division by zero
      else
         y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
      end if
   end procedure flin

end submodule galaxy_module_error