!==============================================================================
! DATA_PROFILES - Module for spectrum data profile handling
!==============================================================================
!> Module containing data structures and routines for handling spectral data
!> profiles used in the galaxy/starburst population synthesis code.
!>
!> Extracted from starburst_main.f90 to improve code organization
!==============================================================================
module data_profiles
   use, intrinsic :: iso_fortran_env, only: real32
   implicit none
   private

   ! Public types and procedures
   public :: initialize_data_profiles

   ! Interface blocks for module procedures
   interface
      module subroutine initialize_data_profiles()
      end subroutine initialize_data_profiles
   end interface

contains
   ! Implementation will be in submodule

end module data_profiles