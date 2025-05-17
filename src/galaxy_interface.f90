!==============================================================================
! GALAXY_INTERFACE - Interface module for external functions
!==============================================================================
!> Module providing interfaces for functions in the main module
!> This allows external programs to call these functions without redefining
!> the interfaces.
!==============================================================================
module galaxy_interface
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   
   public :: error_handler
   
   ! Interface for error_handler from the main module
   interface
      subroutine error_handler(routine, message, fatal)
         character(len=*), intent(in) :: routine
         character(len=*), intent(in) :: message
         logical, intent(in), optional :: fatal
      end subroutine error_handler
   end interface
   
end module galaxy_interface