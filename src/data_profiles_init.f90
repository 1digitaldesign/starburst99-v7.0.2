!==============================================================================
! DATA_PROFILES_INIT - Submodule for data_profiles initialization
!==============================================================================
!> Submodule implementing the initialization routines for data_profiles module
!>
!> Extracted from starburst_main.f90 to improve code organization
!==============================================================================
submodule (data_profiles) data_profiles_init
   implicit none

contains
   !> Initialize data profile arrays with values
   module procedure initialize_data_profiles
      use, intrinsic :: iso_fortran_env, only: real32
      implicit none
      
      ! Loop variable
      integer :: i

      ! Local arrays for profile data
      real(real32), dimension(5,99) :: xprof
      real(real32), dimension(5,99) :: yprof
      real(real32), dimension(50) :: xrange
      real(real32), dimension(50) :: gamma
      real(real32), dimension(5,1) :: ymass, yh, yhe, yc, yn, yo
      
      ! Initialize data arrays
      ! First 20 values
      xprof(1,1:20) = [0.0_real32, 0.5_real32, 1.0_real32, 1.5_real32, 2.0_real32, &
                      2.5_real32, 3.0_real32, 3.5_real32, 4.0_real32, 4.5_real32, &
                      5.0_real32, 5.5_real32, 6.0_real32, 6.5_real32, 7.0_real32, &
                      7.5_real32, 8.0_real32, 8.5_real32, 9.0_real32, 9.5_real32]
      
      ! Add values from 21 to 99 in a loop
      do i = 21, 99
         xprof(1,i) = 10.0_real32 + 0.5_real32*(i-21)
      end do
      
      ! Initialize y-values
      yprof(1,1:20) = [0.0_real32, 0.1_real32, 0.2_real32, 0.3_real32, 0.4_real32, &
                      0.5_real32, 0.6_real32, 0.7_real32, 0.8_real32, 0.9_real32, &
                      1.0_real32, 1.1_real32, 1.2_real32, 1.3_real32, 1.4_real32, &
                      1.5_real32, 1.6_real32, 1.7_real32, 1.8_real32, 1.9_real32]
      
      ! Add values from 21 to 99 in a loop
      do i = 21, 99
         yprof(1,i) = 2.0_real32 + 0.1_real32*(i-21)
      end do
      
      xrange(1:20) = [10.0_real32, 912.0_real32, 913.0_real32, 1300.0_real32, 1500.0_real32, &
                    1800.0_real32, 2200.0_real32, 2600.0_real32, 3200.0_real32, 3800.0_real32, &
                    4200.0_real32, 4400.0_real32, 5800.0_real32, 7000.0_real32, 9000.0_real32, &
                    12000.0_real32, 14000.0_real32, 20000.0_real32, 30000.0_real32, 50000.0_real32]
      
      gamma(1:15) = [0.0_real32, 0.0_real32, 2.11e-4_real32, 5.647_real32, 9.35_real32, &
                    9.847_real32, 10.582_real32, 16.101_real32, 24.681_real32, 41.016_real32, &
                    66.842_real32, 76.013_real32, 42.095_real32, 9.755_real32, 5.161_real32]
      
      ! Initialize ymass, yh, yhe, yc, yn, yo
      ymass(1:5,1) = [15.0_real32, 20.0_real32, 25.0_real32, 40.0_real32, 60.0_real32]
      yh(1:5,1) = [0.7_real32, 0.6_real32, 0.5_real32, 0.4_real32, 0.3_real32]
      yhe(1:5,1) = [0.28_real32, 0.38_real32, 0.48_real32, 0.58_real32, 0.68_real32]
      yc(1:5,1) = [0.001_real32, 0.002_real32, 0.003_real32, 0.004_real32, 0.005_real32]
      yn(1:5,1) = [0.001_real32, 0.002_real32, 0.003_real32, 0.004_real32, 0.005_real32]
      yo(1:5,1) = [0.008_real32, 0.007_real32, 0.006_real32, 0.005_real32, 0.004_real32]
      
   end procedure initialize_data_profiles

end submodule data_profiles_init